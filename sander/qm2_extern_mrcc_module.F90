#include "../include/dprec.fh"
module qm2_extern_mrcc_module
! ----------------------------------------------------------------
! Interface for MRCC based QM/MM MD 
!
! Currently supports:
! -pure QM
! -QM/MM with mechanical and electronic embedding. 
! -QM/MM with cutoff for QM-MM electrostatics
! -QM/MM with exact QM/QM embedding techniques:
!  -Huzinaga embedding
!  -Embedding wavefunction into local correlation methods
!  -combined embedding of the latter schemes.
!   (for the description and performance of these methods, see:
!    B. Hégely, P. Nagy, Gy. Ferenczy, and M. Kállay
!    J. Chem. Phys. 2016, 145, 064107
! -Path Integral Molecular Dynamics (PIMD) with QM/MM (not tested)
! -Replica Exchange Molecular Dynamics (REMD) with QM/MM (not tested)
!
! Implementation by
! Bence Hégely (hegelyb@mail.bme.hu)
!
! Based on Matthew Clark / Andreas Goetz's work 
! qm2_extern_gau_module.F90
!
! Date: August 2016
!
! Integration into Amber release version and tests by
! Andreas Goetz, February 2017
!
! ----------------------------------------------------------------
!
  use qm2_extern_util_module,only:debug_enter_function,&
                                  debug_exit_function
  use qmmm_module,           only:qmmm_nml,qmmm_struct 
   
  implicit none
  !
  private
  public :: get_mrcc_forces
  !
  character(len=*),parameter :: module_name="qm2_extern_mrcc_module"
   
  type mrcc_nml_type
     character(len=20) :: calc
     character(len=20) :: dft
     character(len=20) :: basis
     character(len=50) :: mem
     character(len=20) :: embed
     character(len=20) :: corembed
     integer,pointer   :: embedatoms(:) =>null()
     integer,pointer   :: corembedatoms(:) =>null()
     integer           :: nmo_embed
     integer           :: nmo_corembed
     integer           :: nprintlog
     integer           :: ntpr
     integer           :: debug
     integer           :: verbosity
     logical           :: use_template
     logical           :: do_dipole
   
     ! Deprecated
     integer :: charge
     integer :: spinmult
  end type mrcc_nml_type

contains

  ! --------------------------------------------
  ! Get QM energy and forces from MRCC
  ! --------------------------------------------
  subroutine get_mrcc_forces(do_grad,&
                             nstep,&
                             ntpr_default,&
                             id,&
                             nqmatoms,&
                             qmcoords,&
                             qmtypes,&
                             nclatoms,&
                             clcoords,&
                             escf,&
                             dxyzqm,&
                             dxyzcl,&
                             charge,&
                             spinmult)
  !
    use qm2_extern_util_module, only: print_results,&
                                      check_installation,&
                                      write_dipole
    use constants,              only: CODATA08_AU_TO_KCAL,&
                                      CODATA08_A_TO_BOHRS,&
                                      ZERO
    use file_io_dat
    !     
    use qmmm_module
    !     
    implicit none
    !     
    logical, intent(in)         :: do_grad      ! Return gradient/not
    integer, intent(in)         :: nstep        ! MD step number
    integer, intent(in)         :: ntpr_default ! frequency of printing
    character(len=3),intent(in) :: id           ! ID number for PIMD or 
    ! REMD
    integer, intent(in)         :: nqmatoms     ! Number of QM atoms
    ! QM atoms coordinates
    _REAL_,  intent(in)         :: qmcoords(3,nqmatoms) 
    integer, intent(in)         :: qmtypes(nqmatoms)    ! QM atom types 
    ! (nuclear charge in au)
    ! Number of MM point charges
    integer, intent(in)         :: nclatoms             ! Number of MM atoms
    ! MM point charge coordinates and magnitudes (in Angs and AU)
    _REAL_,  intent(in)         :: clcoords(4,nclatoms)
    _REAL_, intent(out)         :: escf                 ! SCF energy 
    _REAL_, intent(out)         :: dxyzqm(3,nqmatoms)   ! QM force
    _REAL_, intent(out)         :: dxyzcl(3,nclatoms)   ! MM force
    _REAL_                      :: dipxyz(3), dipole    ! Dipole moment 
                                                        ! direction and 
                                                        ! magnitude
    integer, intent(in)         :: charge, spinmult     ! Charge and spin 
                                                        ! multiplicity
    
    type(mrcc_nml_type), save   :: mrcc_nml
    logical, save               :: first_call = .true.
    logical                     :: exist
    integer                     :: i
    integer                     :: printed =-1 ! Used to tell if we 
                                               ! have printed this step yet 
                                               ! since the same step may 
                                               ! be called multiple times.
    ! buffer variable for system call
    character(len=1024)         :: call_buffer 
    character(len=80), save     :: program  = 'dmrcc' 
    ! File names and extensions
    character(len=*),parameter  :: basename = 'mrcc_job'
    character(len=*),parameter  :: inpext   = '.inp'
    character(len=*),parameter  :: logext   = '.log'
    character(len=*),parameter  :: datext   = '.dat'! for MRCC output 
                                                    ! data file
    character(len=*),parameter  :: dipext   = '.dip'
    character(len=*),parameter  :: tplext   = '.tpl'
    ! HLSCF model calculation output
    character(len=*),parameter  :: modext   = '.mod' 
    character(len=14)           :: inpfile,datfile,logfile,dipfile,&
                                   tplfile
    character(len=256)          :: inttochr
    integer, save               :: nprintlogstep
    ! Need to prepend subdirectory if doing REMD, PIMD or multi-region 
    ! QM/MM. 
    ! This is triggered if 'id' is defined (not empty). 
    character(len=25)           :: subdir 
    
    ! Embedding   
    integer,save                :: max_embed,max_corembed
    character(len=1024),save    :: iembed,icorembed
    ! for system call
    integer :: system
    integer :: stat
  
    ! assemble input - / output data filenames
    inpfile=basename//trim(id)//inpext
    logfile=basename//trim(id)//logext
    dipfile=basename//trim(id)//dipext
    datfile=basename//trim(id)//datext
    tplfile=basename//tplext
  
    ! Initalize subdir or sander_bomb will be called if single point
    ! calculation is called.
    subdir=''
  
    ! Setup on first call
    if(first_call) then
       call get_namelist( ntpr_default,mrcc_nml,max_embed,max_corembed)
       nprintlogstep=0
       iembed=' '
       icorembed=' '
       call embed_list(max_embed,max_corembed,iembed,icorembed,&
                       mrcc_nml)
       write(6,'(/,a,/)') '  >>> Running calculations with MRCC <<<'
       call print_namelist(mrcc_nml,max_embed,max_corembed) 
       ! Check for MRCC's driver program dmrcc's full path; 
       ! store as 'program'
       call check_installation(trim(program),id,.true., &
                               mrcc_nml%debug,found=exist)
       ! Remove old MINP, inpfile, logfile, dipfile and datfile at the 
       ! beginning of a run so only the latest run is stored.
       stat=system('rm -f '//inpfile//' '//dipfile//' '&
                   //logfile//' MINP mrcc_job.dat')
       if(stat/=0) then
          call sander_bomb('get_mrcc_forces (qm2_extern_mrcc_module)', & 
               'Error with system call (removing files)', &
               'Will quit now.')
       endif
       write(6,'(80a)') ('-', i=1,80)
       write(6,'(a)') '   4.  RESULTS'
       write(6,'(80a)') ('-', i=1,80)
    endif
  
    if(do_grad.eqv..false.) then
       if(first_call.eqv..false.) then
          if(mrcc_nml%debug>0) then
             write (6,'(/,a,/)') '>>> Skipping second MRCC call &
                                 & for single point calculation <<<'
          endif
          goto 1000
       endif
    endif
  
    call write_inpfile(trim(inpfile),&
                       trim(tplfile),&
                       nqmatoms,&
                       qmcoords,&
                       qmtypes,&
                       nclatoms,&
                       clcoords,&
                       mrcc_nml,&
                       do_grad,&
                       charge,&
                       spinmult,&
                       max_embed,&
                       max_corembed,&
                       iembed,&
                       icorembed)
  
    if(mrcc_nml%debug>0) then
       write(6,'(a)') ' Input file is written successfully'
       write(6,'(a)') ' Calling MRCC driver program: dmrcc...'
    endif
  
    ! Run dmrcc.
    ! Separate runs into different directories if we are doing PIMD or REMD.
    ! Note that dmrcc only running with with the input filename MINP, hence
    ! we have to rename the inpfile to MINP.
    subdir=''
    call_buffer='cp '//inpfile//' MINP;'
    if(trim(id)/='') then 
       subdir='./'//trim(id)//'/'
       call_buffer=' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//&
                   '; cp ../'//inpfile//' .; cp '//inpfile//' MINP;'
    endif
    call_buffer=trim(call_buffer)//' '//trim(program)//' MINP >'//&
                logfile
    stat=system(call_buffer) 
    if(stat/=0) then
       call sander_bomb('get_mrcc_forces (qm2_extern_mrcc_module)', & 
                        'Error with system call (executing dmrcc)', &
                        'Will quit now.')
    endif
  
    if(mrcc_nml%debug>0) then    
       write(6,'(a)') ' MRCC execution success. &
                        & Processing MRCC results...'
    endif
  
    ! Call read_results - retrieve data from MRCC .dat file.
    ! Search in subdir for logfile if doing PIMD.
    ! Otherwise, search current directory.
    1000 continue
    call read_results(datfile,&
                      nqmatoms,&
                      escf,&
                      dxyzqm,&
                      nclatoms,&
                      dxyzcl,&
                      dipxyz,&
                      dipole,&
                      do_grad,&
                      mrcc_nml%debug,&
                      mrcc_nml%do_dipole)
  
    ! Call write_dipole. 
    ! write direction (dipxyz) and magnitude (dipole)
    ! of the QM region's dipole moment to dipfile.
    if(mrcc_nml%ntpr>0.and.mod(nstep,mrcc_nml%ntpr)==0) then
       if(printed/=nstep.and.mrcc_nml%do_dipole) then
          call write_dipole(trim(dipfile),dipxyz,dipole,mrcc_nml%debug)
          printed=nstep
       endif
    endif
  
    ! Save copy of the last input and log files.
    stat=system('cp '//trim(subdir)//inpfile//' '//trim(subdir)//&
                'old.'//inpfile)
    stat=stat+system('cp '//trim(subdir)//logfile//' '&
                      //trim(subdir)//'old.'//logfile)
    if(mrcc_nml%nprintlog.ne.0) then
       nprintlogstep=nprintlogstep+1
       if(nprintlogstep.eq.mrcc_nml%nprintlog) then
          write(inttochr,'(i256)') nstep
          stat=system('cp '//trim(subdir)//logfile//' '&
                     //trim(subdir)//'step.'//trim(adjustl(inttochr))&
                     //'.'//logfile)
          nprintlogstep=0
       endif
    endif
  
    if(stat/=0) then
       call sander_bomb('get_mrcc_forces (qm2_extern_mrcc_module)', & 
                        'Error with system call (moving / &
                         & removing files)','Will quit now.')
    endif
  
    ! F = E*q to get gradients
    ! Note dxyzcl is currently holding the electric field strength.
    do i=1,nclatoms
       dxyzcl(:,i)=-dxyzcl(:,i)*clcoords(4,i)
    enddo
  
    ! Convert gradient from au to kcal/(mol*A).
    if(do_grad) then
       dxyzqm(:,:)=dxyzqm(:,:)*CODATA08_AU_TO_KCAL*CODATA08_A_TO_BOHRS
       if(nclatoms>0) then
          dxyzcl(:,:)=dxyzcl(:,:)*CODATA08_AU_TO_KCAL*CODATA08_A_TO_BOHRS
       endif
    else
       dxyzqm=ZERO
       if(nclatoms>0) dxyzcl=ZERO
    endif
  
    ! Convert final QM energy from au to kcal/(mol*A).
    escf=escf*CODATA08_AU_TO_KCAL
  
    call print_results('qm2_extern_mrcc_module',escf,nqmatoms,dxyzqm,&
                       mrcc_nml%debug,nclatoms,dxyzcl)
    if(first_call) first_call=.false.

  end subroutine get_mrcc_forces

  ! ---------------------------------------------
  ! Read MRCC namelist values from file mdin,
  ! use default values if none are present.
  ! ---------------------------------------------
  subroutine get_namelist(ntpr_default, mrcc_nml,&
                              max_embed,max_corembed)

    use qmmm_module,           only:qmmm_struct

    implicit none
    integer, intent(in)              :: ntpr_default
    type(mrcc_nml_type), intent(out) :: mrcc_nml

    character(len=20)                :: calc,dft,basis
    character(len=20)                :: embed,corembed
    integer                          :: embedatoms(10000)
    integer                          :: corembedatoms(10000)
    integer                          :: nmo_embed,nmo_corembed
    character(len=50)                :: mem
    integer                          :: debug,verbosity,nprintlog
    integer                          :: ntpr,use_template
    integer                          :: charge,spinmult ! deprecated
    integer                          :: do_dipole

    namelist /mrcc/ calc,dft,basis,mem,ntpr,&
                    debug,use_template,charge,spinmult,do_dipole,&
                    verbosity,embed,corembed,embedatoms,&
                    corembedatoms,nmo_embed,nmo_corembed,nprintlog

    integer                          :: ierr,i,j,k
    integer, intent(inout)           :: max_embed,max_corembed

    ! Set default values for mrcc namelist values.
    calc         = 'SCF'
    basis        = '6-31G*'
    mem          = '256MB'
    ntpr         = ntpr_default
    nprintlog    = 0
    debug        = 0
    use_template = 0
    do_dipole    = 0
    dft          = 'off'
    verbosity    = 2
    embed        = 'off'
    corembed     = 'off'

    ! These are now deprecated and should be specified in the .
    ! &qmmmm namelist
    charge   = -351
    spinmult = -351
    ! These are now deprecated and should be specified in the .
    ! &mrcc namelist
    embedatoms   = 0
    corembedatoms= 0
    max_embed    = 0
    max_corembed = 0
    nmo_embed    = 0
    nmo_corembed = 0

    ! Read namelist
    rewind 5
    read(5,nml=mrcc,iostat=ierr)

    if(ierr>0) then
       call sander_bomb('get_namelist (qm2_extern_mrcc_module)', &
                        'mrcc namelist read error', &
                        'Please check your input.')
    elseif (ierr<0) then
       write(6,'(a,/,a)')'mrcc namelist read encountered end of file.',&
                         'Please check your input if the calculation &
                          & encounters a problem.'
    endif

    if(charge/=-351.or.spinmult/=-351 ) then
       call sander_bomb('get_namelist (qm2_extern_mrcc_module)', &
                        'The charge and spin keywords are deprecated.',&
                        'Please specify charge (qmcharge) and spin &
                         & multiplicity (spin) in the qmmm namelist.')
    endif

    if(embed.ne.'off') allocate(mrcc_nml%embedatoms(10000))
    if(corembed.ne.'off') allocate(mrcc_nml%corembedatoms(10000))
    ! Count the number of embed and corembedatoms
    if(embed.ne.'off'.or.corembed.ne.'off') then
       do i=1,qmmm_struct%nquant
          if(embedatoms(i).ne.0) max_embed=max_embed+1
          if(corembedatoms(i).ne.0) max_corembed=max_corembed+1
       enddo
    endif
    if(embed.ne.'off'.and.max_embed.gt.0) then 
       allocate(mrcc_nml%embedatoms(max_embed))
    endif
    if(corembed.ne.'off'.and.max_corembed.gt.0) then
       allocate(mrcc_nml%corembedatoms(max_corembed))
    endif
    ! Assign namelist values to mrcc_nml data type
    mrcc_nml%calc         = calc
    mrcc_nml%dft          = dft
    mrcc_nml%basis        = basis
    mrcc_nml%mem          = mem
    mrcc_nml%ntpr         = ntpr
    mrcc_nml%nprintlog    = nprintlog
    mrcc_nml%verbosity    = verbosity
    mrcc_nml%debug        = debug
    mrcc_nml%embed        = embed
    mrcc_nml%corembed     = corembed
    if(max_embed.ne.0) then
       call sortatoms(max_embed,embedatoms)
       mrcc_nml%embedatoms(1:max_embed)=embedatoms(1:max_embed)
       mrcc_nml%nmo_embed = nmo_embed
    endif
    if(max_corembed.ne.0) then
       call sortatoms(max_corembed,corembedatoms)
       mrcc_nml%corembedatoms(1:max_corembed)=&
            corembedatoms(1:max_corembed)
       mrcc_nml%nmo_corembed = nmo_corembed
    endif
    ! Check if the user specified embedatoms and corembedatoms
    ! correctly in case of 3 layer calculation.
    if(embed.ne.'off'.and.max_embed.eq.0) then
       call sander_bomb('get_namelist (qm2_extern_mrcc_module)', &
                        'Illegal embedatoms specification!','Please &
                         & check your input. embedatoms are not specified.')
    endif
    if(embed.ne.'off'.and.max_embed.ne.0) then
       k=0
       do i=1,qmmm_struct%nquant
          do j=1,max_embed
             if(qmmm_struct%iqmatoms(i).eq.embedatoms(j)) k=k+1
          enddo
       enddo
       if(k.ne.max_embed) then
          call sander_bomb('get_namelist (qm2_extern_mrcc_module)', &
                           'Illegal embedatoms specification!', &
                           'Please check your input. embedatoms have &
                            & to be a subset of iqmatoms.')
       endif
    endif
    if(corembed.ne.'off'.and.max_corembed.eq.0) then
       call sander_bomb('get_namelist (qm2_extern_mrcc_module)', &
                        'Illegal corembedatoms specification!', &
                        'Please check your input. corembedatoms &
                         & are not specified.')
    endif
    if(corembed.ne.'off'.and.max_corembed.ne.0) then
       k=0
       do i=1,qmmm_struct%nquant
          do j=1,max_corembed
             if(qmmm_struct%iqmatoms(i).eq.corembedatoms(j)) k=k+1
          enddo
       enddo
       if(k.ne.max_corembed) then
          call sander_bomb('get_namelist (qm2_extern_mrcc_module)', &
                           'Illegal corembedatoms specification!', &
                           'Please check your input. corembedatoms have &
                            & to be a subset of iqmatoms.')
       endif
    endif
    ! Check: in case of 4 layer calculation, the corembed array has to be
    !        the subset of embed array.
    if(embed.ne.'off'.and.corembed.ne.'off'.and.&
         max_embed.ne.0.and.max_corembed.ne.0) then
       k=0
       do i=1,max_embed
          do j=1,max_corembed
             if(embedatoms(i).eq.corembedatoms(j)) k=k+1
          enddo
       enddo
       if(k.ne.max_corembed) then
          call sander_bomb('get_namelist (qm2_extern_mrcc_module)', &
                           'Illegal corembedatoms specification!', &
                           'Please check your input. corembedatoms &
                            & have to be a subset of embedatoms in &
                            & case of 4 layer calculation.')
       endif
    endif

    if(do_dipole==0 ) then
       mrcc_nml%do_dipole=.false.
    else if( do_dipole==1) then
       mrcc_nml%do_dipole=.true.
    else
       call sander_bomb('get_namelist (qm2_extern_mrcc_module)', &
                        'mrcc dipole value is not allowed!', &
                        'Please check your input. do_dipole can &
                         & only be 0 or 1.')
    endif

    if(use_template==0) then
       mrcc_nml%use_template=.false.
    elseif(use_template==1) then
       mrcc_nml%use_template=.true.
    else
       call sander_bomb('get_namelist (qm2_extern_mrcc_module)', &
                        'mrcc use_template value is not allowed', &
                        'Please check your input. use_template &
                         & can only be 0 or 1.')
    endif

  end subroutine get_namelist

  ! --------------------------------
  ! Print MRCC namelist settings
  ! --------------------------------
  subroutine print_namelist(mrcc_nml,max_embed,max_corembed)
    
    implicit none
    type(mrcc_nml_type),intent(in) :: mrcc_nml
    integer, intent(in) :: max_embed,max_corembed
    integer             :: i
    
    if(mrcc_nml%use_template.eqv..true.) then
       write(6,'(a)')&
            '  Warning! The following keywords maybe overwritten &
            & by the template file!'
    endif
    write(6,'(a)')       '| &mrcc'
    write(6,'(2a)')      '|   method       = ', mrcc_nml%calc
    write(6,'(2a)')      '|   dft          = ', mrcc_nml%dft
    write(6,'(2a)')      '|   basis        = ', mrcc_nml%basis
    write(6,'(2a)')      '|   mem          = ', trim(mrcc_nml%mem)
    write(6,'(a,i0)')    '|   ntpr         = ', mrcc_nml%ntpr
    write(6,'(a,i0)')    '|   nprintlog    = ', mrcc_nml%nprintlog
    write(6,'(a,i2)')    '|   verbosity    = ', mrcc_nml%verbosity
    write(6,'(a,i2)')    '|   debug        = ', mrcc_nml%debug
    write(6,'(a,l)')     '|   do_dipole    = ', mrcc_nml%do_dipole
    write(6,'(a,l)')     '|   use_template = ', mrcc_nml%use_template
    write(6,'(2a)')      '|   embed        = ', mrcc_nml%embed
    if(mrcc_nml%embed.ne.'off') then
       write(6,'(a,i6)')    '|   nmo_embed    = ', mrcc_nml%nmo_embed
    endif
    write(6,'(2a)')      '|   corembed     = ', mrcc_nml%corembed
    if(mrcc_nml%corembed.ne.'off') then
       write(6,'(a,i6)')    '|   nmo_corembed = ', mrcc_nml%nmo_corembed
    endif
    write(6,'(a)')       '| /'
    write(6,'(a)')
    if(max_embed.gt.0.and.max_corembed.gt.0) then
       write(6,'(a)')      '|   Performing 4 layer calculation.'
       write(6,'(a)')
       write(6,'(2a)')     '|   Method of the 2nd layer: ',&
            mrcc_nml%corembed
       write(6,'(a)')      '|   Atoms of the 2nd layer (embedatoms):'
       write(6,'(10i5)') (mrcc_nml%embedatoms(i),i=1,max_embed)
       write(6,'(a)')
       write(6,'(2a)')     '|   Method of the 1st layer:',mrcc_nml%calc
       write(6,'(a)')      '|   Atoms of the 1st layer (corembedatoms):'
       write(6,'(10i5)') (mrcc_nml%corembedatoms(i),i=1,max_corembed)
       write(6,'(a)')
    endif
    if(max_embed.eq.0.and.max_corembed.gt.0) then
       write(6,'(a)')      '|   Performing 3 layer calculation.'
       write(6,'(a)')
       write(6,'(2a)')     '|   Method of the 1st layer: ',mrcc_nml%calc
       write(6,'(a)')      '|   Atoms of the 1st layer (corembedatoms):'
       write(6,'(10i5)') (mrcc_nml%corembedatoms(i),i=1,max_corembed)
       write(6,'(a)')
    endif
    if(max_embed.gt.0.and.max_corembed.eq.0) then
       write(6,'(a)')      '|   Performing 3 layer calculation.'
       write(6,'(a)')
       write(6,'(2a)')     '|   Method of the 1st layer: ',mrcc_nml%calc
       write(6,'(a)')      '|   Atoms of the 1st layer (embedatoms):'
       write(6,'(10i5)') (mrcc_nml%embedatoms(i),i=1,max_embed)
       write(6,'(a)')
    endif

  end subroutine print_namelist

  ! -----------------------------
  ! Write input file for MRCC
  ! -----------------------------
  subroutine write_inpfile(inpfile,&
                           tplfile,&
                           nqmatoms,&
                           qmcoords,&
                           qmtypes,&
                           nclatoms,&
                           clcoords,&
                           mrcc_nml,&
                           do_grad,&
                           charge,&
                           spinmult,&
                           max_embed,&
                           max_corembed,&
                           iembed,&
                           icorembed)

    use ElementOrbitalIndex, only : elementSymbol

    implicit none

    character(len=*), intent(in)   :: inpfile,tplfile
    integer, intent(in)            :: nqmatoms
    _REAL_,  intent(in)            :: qmcoords(:,:)
    integer, intent(in)            :: qmtypes(:)
    integer, intent(in)            :: nclatoms
    _REAL_,  intent(in)            :: clcoords(:,:)
    type(mrcc_nml_type), intent(in):: mrcc_nml
    logical, intent(in)            :: do_grad
    integer, intent(in)            :: charge, spinmult

    integer, parameter             :: iunit=351,tplunit=352
    integer                        :: i,ierr
    logical, save                  :: first_call=.true.

    character(len=8)               :: inttochr

    ! embed and corembed
    integer, intent(in) :: max_embed,max_corembed
    character(len=1024) :: iembed,icorembed

    call debug_enter_function('write_inpfile',module_name,&
                             mrcc_nml%debug)

    ! If we are using template, copy tplfile to inpfile.
    ! Note that MRCC only using the first read keywords.
    if(mrcc_nml%use_template) then
       call system(&
            'echo "! User defined keywords from MRCC template file" >'&
            //inpfile//'; cat '//tplfile//' >> '//inpfile)
       open(iunit, file=inpfile, iostat=ierr, position='append')
    else
       open(iunit, file=inpfile, iostat=ierr, status='replace')
    endif
    ! Write namelist parameters to inpfile
    write(iunit,'(a)')'! User defined keywords from AMBER', &
                         ' control file'
    write(iunit,'(a,a)')'mem=',mrcc_nml%mem
    write(iunit,'(a,a)')'calc=',mrcc_nml%calc
    write(iunit,'(a,a)')'dft=',mrcc_nml%dft
    write(iunit,'(a,a)')'basis=',mrcc_nml%basis
    write(inttochr,'(i5)') mrcc_nml%verbosity
    write(iunit,'(2a)')'verbosity=',trim(adjustl(inttochr))
    write(inttochr,'(i5)') charge
    write(iunit,'(2a)')'charge=',trim(adjustl(inttochr))
    write(inttochr,'(i5)') spinmult
    write(iunit,'(2a)')'mult=',trim(adjustl(inttochr))
    write(iunit,'(a)')'! Default keywords for AMBER-MRCC run'
    write(iunit,'(a)')'qmmm=amber'
    ! Gradient or single point
    if(do_grad) then
       write(iunit,'(a)')'dens=2'
       if(first_call.eqv..false.) then
          write(iunit,'(a)') 'scfiguess=restart'
       endif
    else
       write(iunit,'(a,i1)')'dens=0'
    endif
    if(max_embed.gt.0) then
       write(iunit,'(a)')'embed=huzinaga'
       write(iunit,'(a)') trim(adjustl(iembed))
       write(iunit,'(a)') trim(adjustl(mrcc_nml%embed))
       write(inttochr,'(i5)') mrcc_nml%nmo_embed
       write(iunit,'(a)') trim(adjustl(inttochr))
    endif
    if(max_corembed.gt.0) then
       write(iunit,'(a)')'corembed=on'
       write(iunit,'(a)') trim(adjustl(icorembed))
       write(iunit,'(a)') trim(adjustl(mrcc_nml%corembed))
       write(inttochr,'(i5)') mrcc_nml%nmo_corembed
       write(iunit,'(a)') trim(adjustl(inttochr))
    endif
    ! blank line for MINP format
    write(iunit,'(a)')
    ! Write comment line than make a blank line for MINP format
    if(mrcc_nml%use_template) then
       write(iunit,'(a,/)') 'MRCC run using SANDER external interface &
            & with template '//tplfile
    else
       write(iunit,'(a,/)') 'MRCC run using SANDER external interface.'
    endif
    ! Write the geom keyword
    write(iunit,'(a)') 'geom=xyz'
    ! Link Atom approach
    write(iunit,'(i5,/)') nqmatoms
    do i=1,nqmatoms
       write(iunit,'(a2,1x,3f25.16)') &
            elementSymbol(qmtypes(i)), &
            qmcoords(1:3,i)
    enddo
    write(iunit,'()')
    ! When electrostatic embedding QM/MM is in use 
    ! write MM coordinates with point charges.
    if(nclatoms.gt.0) then
       write(iunit,'(a)') 'pointcharges'
       write(iunit,*) nclatoms
       ! Write point charges, but write nuclear charge for qm ghost atom.
       do i=1,nclatoms
          write(iunit,'(4f25.16)') clcoords(:,i) 
       enddo
       write(iunit,'()')
    else
       write(iunit,'(a)') 'pointcharges'
       write(iunit,'(a1)') '0'
    endif
    write(iunit,'()')

    close(iunit,iostat=ierr)

    if(ierr/=0) then
       call sander_bomb('write_inpfile (qm2_extern_mrcc_module)', &
                        'Error closing MRCC input file after writing',&
                        'Will quit now.')
    endif
    first_call=.false.

    call debug_exit_function('write_inpfile',module_name,&
                              mrcc_nml%debug )

  end subroutine write_inpfile

  subroutine read_results(datfile,&
                          nqmatoms,&
                          escf,&
                          dxyzqm,&
                          nclatoms,&
                          dxyzcl,&
                          dipxyz,&
                          dipole,&
                          do_grad,&
                          debug,&
                          do_dipole)

    implicit none
    
    character(len=*), intent(in)  :: datfile
    integer, intent(in)           :: nqmatoms, nclatoms
    _REAL_,  intent(out)          :: escf,dxyzqm(3,nqmatoms),&
                                     dxyzcl(3,nclatoms) ! dxyzcl returns 
                                                        ! containing the 
                                                        ! electric field 
                                                        ! at x,y,z
    _REAL_,  intent(out)          :: dipxyz(3),dipole
    logical, intent(in)           :: do_grad
    integer, intent(in)           :: debug
    
    _REAL_                        :: pc_self_energy=0 ! Temporary variable to hold self energy of point charges
    integer                       :: ios, i
    integer, parameter            :: iunit = 351
    character(len=8)              :: read_buffer
    logical, intent(in)           :: do_dipole

    call debug_enter_function('read_results',module_name,debug)
    ! Energy and gradient; file will be empty if not do_grad
    open(iunit,file=datfile,iostat=ios)
    if(ios/=0 ) then
       call sander_bomb('read_results (qm2_extern_mrcc_module)', &
            'Error opening MRCC data file '//datfile//' for reading.',&
            ' Will quit now.')
    endif
    do
       read(iunit,'(a)',iostat=ios) read_buffer
       ! End of file; data not found
       if(ios<0) then
          exit
       endif
       ! Read final QM energy.
       if(read_buffer(2:5)==' Fin') then
          read(iunit,*) escf
       endif
       ! Read the self energy of the point charges.
       if(read_buffer(2:5)==' Sel') then
          read(iunit,*) pc_self_energy
       endif
       ! Read gradients of the QM region.
       if(do_grad) then
          if(read_buffer(2:5)=='Forc') then
             do i=1,nqmatoms
                read(iunit,*) dxyzqm(1,i),dxyzqm(2,i),dxyzqm(3,i)
             enddo
          endif
       endif
       ! QM/MM electrostatic embedding
       ! Read electric field of the MM region.
       if(nclatoms>0) then
          if(do_grad) then
             if(read_buffer(2:5)=='Elec') then
                do i=1,nclatoms
                   read(iunit,*) dxyzcl(1,i),dxyzcl(2,i),dxyzcl(3,i)
                enddo
             endif
          endif
       endif
       ! Read magnitude of dipole moment vector.
       if(do_dipole) then
          if(read_buffer(2:5) =='Magn') then
             read(iunit,*) dipole
          end if
          ! Read direction of dipole moment vector
          if(read_buffer(2:5) == 'Dire' ) then
             read(iunit,*) dipxyz(1), dipxyz(2), dipxyz(3)
          endif
       endif
    enddo
    close(iunit)
    ! Return the corrected QM energy
    escf=escf-pc_self_energy

    call debug_exit_function('read_results',module_name,debug)

  end subroutine read_results

  subroutine embed_list(max_embed,&
                        max_corembed,&
                        iembed,&
                        icorembed,&
                        mrcc_nml)
    use qmmm_module,           only:qmmm_struct
    implicit none
    integer, intent(in)             :: max_embed,max_corembed
    character(len=1024),intent(out) :: iembed,icorembed
    character(len=1024)             :: embed_line,corembed_line
    type(mrcc_nml_type),intent(in)  :: mrcc_nml
    integer                         :: int_embed(max_embed)
    integer                         :: int_corembed(max_corembed)
    character(len=1)                :: embed_vline(1024)
    character(len=1)                :: corembed_vline(1024)
    character(len=16)               :: buffer_line
    character(len=1)                :: buffer_vline(16)
    equivalence(embed_line,embed_vline)
    equivalence(corembed_line,corembed_vline)
    equivalence(buffer_line,buffer_vline)
    !
    integer                      :: i,j,k

    ! Embed: -determine the number of embedatoms
    !        -convert natom indexes to qmlist indexes
    !        -convert integer array to character array 
    !         (mrcc requirement)
    if(mrcc_nml%embed.ne.'off'.and.max_embed.ne.0) then
       k=0
       do i=1,qmmm_struct%nquant
          do j=1,max_embed
             if(qmmm_struct%iqmatoms(i).eq.mrcc_nml%embedatoms(j)) then
                k=k+1
                int_embed(k)=i
             endif
          enddo
       enddo
       ! convert integers to character array
       j=0
       embed_line=''
       do i=1,max_embed
          buffer_line=''
          write(buffer_line,*) int_embed(i)
          buffer_line=trim(adjustl(buffer_line))
          k=1
          j=j+1
          do while(buffer_vline(k).ne.' ')
             embed_vline(j)=buffer_vline(k)
             k=k+1 
             j=j+1
          enddo
          if(i.ne.max_embed) embed_vline(j)=','
       enddo
       iembed=''
       iembed=embed_line
    endif
    ! Corembed: same operations for corembedatoms
    if(mrcc_nml%corembed.ne.'off'.and.max_corembed.ne.0) then
       k=0
       do i=1,qmmm_struct%nquant
          do j=1,max_corembed
             if(qmmm_struct%iqmatoms(i).eq.mrcc_nml%corembedatoms(j)) then
                k=k+1
                int_corembed(k)=i
             endif
          enddo
       enddo
       j=0
       corembed_line=''
       do i=1,max_corembed
          buffer_line=''
          write(buffer_line,*) int_corembed(i)
          buffer_line=trim(adjustl(buffer_line))
          k=1
          j=j+1
          do while(buffer_vline(k).ne.' ')
             corembed_vline(j)=buffer_vline(k)
             k=k+1
             j=j+1
          enddo
          if(i.ne.max_corembed) corembed_vline(j)=','
       enddo
       icorembed=''
       icorembed=corembed_line
    endif

  end subroutine embed_list

  subroutine sortatoms(nat,arr1)
    implicit none
    integer,intent(inout) :: nat
    integer,intent(inout) :: arr1(nat)
    integer :: locmin,n,tmp
    do n=1,nat-1
       locmin=n-1+minloc(arr1(n:nat),1)
       tmp=arr1(n)
       arr1(n)=arr1(locmin)
       arr1(locmin)=tmp
    enddo

  end subroutine sortatoms

end module qm2_extern_mrcc_module
