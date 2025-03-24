#include "../include/dprec.fh"
module qm2_extern_quick_module
! ----------------------------------------------------------------
! Interface for QUICK based QM MD
!
! Currently supports:
! pure QM
!
! Initial implementation by
! --
! under supervision of
! Andreas Goetz (SDSC)
!
! Date: February 2011
!       April 2015
!       September 2019
! Extensions by Andreas Goetz (SDSC)
!
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_quick_forces

  character(len=*), parameter :: module_name = "qm2_extern_quick_module"

  type quick_nml_type
     character(len=40) :: method
     character(len=20) :: basis
     character(len=100) :: executable
     character(len=100) :: do_parallel
     integer :: scf_cyc
     _REAL_  :: denserms
     integer :: ntpr
     integer :: debug
     logical :: dipole
     logical :: use_template

     ! Deprecated
     integer :: charge
     integer :: spinmult
  end type quick_nml_type

contains

  ! --------------------------------------------
  ! Get QM energy and forces from QUICK
  ! --------------------------------------------
  subroutine get_quick_forces( do_grad, nstep, ntpr_default, id, &
       nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
       escf, dxyzqm, dxyzcl, charge, spinmult )

    use qm2_extern_util_module, only: print_results, check_installation, write_dipole
    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO
    use file_io_dat

    implicit none

    logical, intent(in) :: do_grad              ! Return gradient/not
    integer, intent(in) :: nstep                ! MD step number
    integer, intent(in) :: ntpr_default         ! frequency of printing
    character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
    integer, intent(in) :: nqmatoms             ! Number of QM atoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
    integer, intent(in) :: qmtypes(nqmatoms)    ! QM atom types (nuclear charge in au)
    integer, intent(in) :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
    _REAL_, intent(out) :: escf                 ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
    _REAL_              :: dipxyz(3), dipole    ! Dipole moment
    integer, intent(in) :: charge, spinmult     ! Charge and spin multiplicity

    type(quick_nml_type), save   :: quick_nml
    logical, save                :: first_call = .true.
    logical                      :: exist
    integer                      :: i
    integer                      :: printed =-1 ! Used to tell if we have printed this step yet
                                                ! since the same step may be called multiple times
    character(len=1024)          :: call_buffer
    character(len=80), save      :: program     = 'quick'
    character(len=*), parameter  :: basename = 'QUICK_job'
    character(len=*), parameter  :: inpext = '.in'
    character(len=*), parameter  :: logext = '.log'
    character(len=*), parameter  :: outext = '.out'
    character(len=*), parameter  :: dipext = '.dip'
    character(len=*), parameter  :: tplext = '.tpl'
    character(len=14)            :: inpfile, rstfile, frstfile, outfile, dipfile, tplfile
    ! Need to prepend subdirectory if doing REMD, PIMD or multi-region QM/MM.
    !   This is triggered if 'id' is defined (not empty).
    character(len=25)            :: subdir

    ! for system call
    integer :: system
    integer :: stat

    ! assemble input - / output data filenames
    inpfile = basename//trim(id)//inpext
    outfile = basename//trim(id)//outext
    dipfile = basename//trim(id)//dipext
    tplfile = basename//tplext

    ! Setup on first call
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') '  >>> Running QM calculation with QUICK <<<'
      call get_namelist( ntpr_default, quick_nml )
      call print_namelist( quick_nml )

      ! Check for version of QUICK to use; store in program
      if ( trim(quick_nml%executable) /= '' ) then
        ! Will quit if this is not found
        call check_installation( trim(quick_nml%executable), id, .true., quick_nml%debug )
        program = quick_nml%executable
      else
        ! Try in order quick, quick.cuda producing an error if none are found.
        call check_installation( 'quick', id, .false., quick_nml%debug, found=exist )
        if ( exist ) then
          program = 'quick'
        else
          call check_installation( 'quick.cuda', id, .true., quick_nml%debug, found=exist )
          program = 'quick.cuda'
        end if
     end if

      write (6,'(80a)') ('-', i=1,80)
      write (6,'(a)') '   4.  RESULTS'
      write (6,'(80a)') ('-', i=1,80)
      ! Remove old inpfile, outfile and dipfile at the
      ! beginning of a run so only the latest run is stored.
      stat = system('rm -f '//inpfile//' '//outfile//' '//dipfile)

      if ( stat /= 0 ) then
        call sander_bomb('get_quick_forces (qm2_extern_quick_module)', &
          'Error with system call (removing files)', &
          'Will quit now.')
      end if
    end if

    call write_inpfile( trim(inpfile), trim(tplfile), &
         nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
         quick_nml, do_grad, charge, spinmult)

    if ( quick_nml%debug > 0 ) then
      write(6,'(a)') ' Input file written successfully; calling QUICK...'
    end if


!----
    ! Run QUICK
    ! Separate runs into different directories if we are doing PIMD or REMD
    subdir=''
    call_buffer=''
    if (trim(id)/='') then
      subdir='./'//trim(id)//'/'
      call_buffer=' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//'; mv ../'//inpfile//' .;'
    end if
    call_buffer = trim(call_buffer)//' '//trim(quick_nml%do_parallel)//' '//trim(program)//' '//inpfile
    stat = system(call_buffer)
    if ( stat /= 0 ) then
      call sander_bomb('get_quick_forces (qm2_extern_quick_module)', &
        'Error with system call (executing QUICK)', &
        'Will quit now.')
    end if

    if ( quick_nml%debug > 0 ) then
      write(6,'(a)') ' QUICK execution success; Processing QUICK results...'
    end if


    ! Call read_results - retrieve data from QUICK .out file
    ! Will output data to escf and dxyqm for pure QM or MM runs
    ! For QM/MM runs will also return dxyzcl containing electric field strength
!----


    ! Search in subdir for outfile if doing PIMD
    ! Otherwise, search current directory
    call read_results( trim(subdir)//trim(outfile), &
         nqmatoms, escf, dxyzqm, nclatoms, dxyzcl, dipxyz, dipole, &
         do_grad, quick_nml%debug )

    ! Call write_dipole with dipfile, dipxyz, and magnitude of dipxyz;
    ! will write output to dipfile
    if ( quick_nml%ntpr > 0 .and. mod(nstep, quick_nml%ntpr) == 0 ) then
      if ( printed /= nstep .and. quick_nml%dipole ) then
        call write_dipole( trim(dipfile), dipxyz, dipole, quick_nml%debug )
        printed = nstep
      end if
    end if

    ! Save copy of last input and out files
    stat = system('mv '//trim(subdir)//inpfile//' '//trim(subdir)//'old.'//inpfile)
    stat = stat + system('mv '//trim(subdir)//outfile//' '//trim(subdir)//'old.'//outfile)
    if ( stat /= 0 ) then
      call sander_bomb('get_quick_forces (qm2_extern_quick_module)', &
        'Error with system call (moving / removing files)', &
        'Will quit now.')
    end if

    ! Convert gradient from au to kcal/(mol*A)
    if ( do_grad ) then
      dxyzqm(:,:) = dxyzqm(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
      if ( nclatoms > 0 ) then
        dxyzcl(:,:) = dxyzcl(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
      end if
    else
      dxyzqm = ZERO
      if ( nclatoms > 0 ) dxyzcl = ZERO
    end if

    escf = escf * CODATA08_AU_TO_KCAL

    call print_results( 'qm2_extern_quick_module', escf, nqmatoms, dxyzqm, &
      quick_nml%debug, nclatoms, dxyzcl )

  end subroutine get_quick_forces

  ! ---------------------------------------------
  ! Read QUICK namelist values from file mdin,
  ! use default values if none are present.
  ! ---------------------------------------------


  subroutine get_namelist(ntpr_default, quick_nml)

    implicit none
    integer, intent(in) :: ntpr_default
    type(quick_nml_type), intent(out) :: quick_nml

    character(len=40) :: method, basis, executable, do_parallel
    integer :: debug
    integer :: scf_cyc, ntpr, dipole, use_template
    _REAL_  :: denserms
    integer :: charge, spinmult ! deprecated

    namelist /quick/ method, basis, executable, do_parallel, scf_cyc, denserms, &
      ntpr, debug, dipole, use_template, charge, spinmult

    integer :: ierr

    ! Set default values for quick namelist values
    method       = 'BLYP'
    basis        = '6-31G'
    executable   = 'quick'
    do_parallel  = ''
    scf_cyc      = 200
    denserms     = 1.0E-6
    ntpr         = ntpr_default
    debug        = 0
    dipole       = 0
    use_template = 0

    ! These are now deprecated and should be specified in the &qmmmm namelist
    charge   = -351
    spinmult = -351

    ! Read namelist
    rewind 5
    read(5,nml=quick,iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_namelist (qm2_extern_quick_module)', &
            '&quick namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a,/,a)') '&quick namelist read encountered end of file', &
            'Please check your input if the calculation encounters a problem'
    end if

    if ( charge /= -351 .or. spinmult /= -351 ) then
      call sander_bomb('get_namelist (qm2_extern_quick_module)', &
        'The charge and spin keywords are deprecated', &
        'Please specify charge (qmcharge) and spin multiplicity (spin) in the &qmmm namelist.')
    end if

    ! Assign namelist values to quick_nml data type
    quick_nml%method       = method
    quick_nml%basis        = basis
    quick_nml%executable   = executable
    quick_nml%do_parallel  = do_parallel
    quick_nml%scf_cyc      = scf_cyc
    quick_nml%denserms     = denserms
    quick_nml%ntpr         = ntpr
    quick_nml%debug        = debug
    if ( dipole == 0 ) then
       quick_nml%dipole = .false.
    else if ( dipole == 1 ) then
       quick_nml%dipole = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_quick_module)', &
            '&quick dipole value not allowed', &
            'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
       quick_nml%use_template = .false.
    else if ( use_template == 1 ) then
       quick_nml%use_template = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_quick_module)', &
            '&quick use_template value not allowed', &
            'Please check your input. use_template can only be 0 or 1.')
    end if

  end subroutine get_namelist

  ! --------------------------------
  ! Print QUICK namelist settings
  ! --------------------------------
  subroutine print_namelist(quick_nml)

    implicit none
    type(quick_nml_type), intent(in) :: quick_nml

    write(6, '(a)')       '| &quick'
    write(6, '(2a)')      '|   method       = ', quick_nml%method
    write(6, '(2a)')      '|   basis        = ', quick_nml%basis
    write(6, '(2a)')      '|   executable   = ', quick_nml%executable
    write(6, '(2a)')      '|   do_parallel  = ', quick_nml%do_parallel
    write(6, '(a,i0)')    '|   scf_cyc      = ', quick_nml%scf_cyc
    write(6, '(a,f15.10)')'|   denserms     = ', quick_nml%denserms
    write(6, '(a,i0)')    '|   ntpr         = ', quick_nml%ntpr
    write(6, '(a,i2)')    '|   debug        = ', quick_nml%debug
    write(6, '(a,l)')     '|   dipole       = ', quick_nml%dipole
    write(6, '(a,l)')     '|   use_template = ', quick_nml%use_template
    write(6,'(a)')        '| /'

  end subroutine print_namelist

  ! -----------------------------
  ! Write input file for QUICK
  ! -----------------------------

  subroutine write_inpfile( inpfile, tplfile, &
       nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
       quick_nml, do_grad, charge, spinmult )

    use ElementOrbitalIndex, only : elementSymbol
    use qm2_extern_util_module, only : is_empty_string
    use UtilitiesModule, only: Upcase

    implicit none

    character(len=*), intent(in)   :: inpfile, tplfile
    integer, intent(in)            :: nqmatoms
    _REAL_,  intent(in)            :: qmcoords(:,:)
    integer, intent(in)            :: qmtypes(:)
    integer, intent(in)            :: nclatoms
    _REAL_,  intent(in)            :: clcoords(:,:)
    type(quick_nml_type), intent(in) :: quick_nml
    logical, intent(in)            :: do_grad
    integer, intent(in)            :: charge, spinmult

    intrinsic :: is_iostat_end

    integer, parameter :: iunit = 351, tplunit = 352
    integer            :: i, ierr
    integer            :: tplerr, ios
    character(len=256) :: route
    character(len=256) :: read_buffer
    logical, save      :: first_call = .true.
    logical            :: empty_line, done

    call debug_enter_function( 'write_inpfile', module_name, quick_nml%debug )

    open(iunit, file=inpfile, iostat=ierr)
    if ( ierr /= 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_quick_module)', &
        'Error opening QUICK inpfile '//inpfile//' for writing', &
        'Will quit now.')
    end if

    ! Assemble route
    ! If using template, write route information in template file to route card
    route = ''
    if ( quick_nml%use_template) then

      open(tplunit, file=tplfile, iostat=tplerr)
      if ( tplerr /= 0 ) then
        call sander_bomb('write_inpfile (qm2_extern_quick_module)', &
          'Error opening QUICK template file '//tplfile//' for reading', &
          'Will quit now.')
      end if
      read (tplunit, '(a)', iostat = ios) read_buffer
      if (ios < 0) then
        call sander_bomb('write_inpfile (qm2_extern_quick_module)', &
          'Error reading QUICK template file '//tplfile, &
          'Will quit now.')
      end if
      close(tplunit)
      write(route,'(a)') trim(route)//' '//trim(read_buffer)//' '

    else

      ! Method/Basis
      route = trim(quick_nml%method)//' BASIS='//trim(quick_nml%basis)
      ! SCF convergence setting
      write(route,'(a,i0,a,f15.10,a,i0,a,i0)') trim(route)//' SCF=', quick_nml%scf_cyc,' DENSERMS=', quick_nml%denserms, &
                                               ' CHARGE=', charge,' MULT=', spinmult
    end if
    ! Gradient or single point
    if ( do_grad ) then
      route = trim(route)//' GRADIENT DIPOLE'
    end if
    ! If doing electrostatic embedding QM/MM,
    ! read external point charges, print to .out
    if ( nclatoms > 0 ) then
      route = trim(route)//' EXTCHARGES'
    end if

    write(iunit,'(a)') trim(route)
    write(iunit,'(a)')

    ! Write QM atoms and coordinates
    do i = 1, nqmatoms
      write(iunit,'(a2,1x,3f25.16)') elementSymbol(qmtypes(i)), qmcoords(1:3,i)
    end do
    write(iunit,'()')

    ! When electrostatic embedding QM/MM is in use
    ! write MM coordinates with point charges
    if ( nclatoms > 0 ) then
      do i = 1, nclatoms
        write(iunit,'(4f21.16)') clcoords(:,i)
      end do
      write(iunit,'()')
    end if

!-------------should be removed
    ! Write generic basis set / ecp information from template file
    ! if available (check for EOF or 2 empty lines)
    if ( quick_nml%use_template) then

      open(tplunit, file=tplfile, iostat=tplerr)
      if ( tplerr /= 0 ) then
        call sander_bomb('write_inpfile (qm2_extern_quick_module)', &
          'Error opening QUICK template file '//tplfile//' for reading', &
          'Will quit now.')
      end if

      ! discard first two lines (route plus empty line)
      done = .false.
      i = 0
      do
         read (tplunit, '(a)', iostat = ios) read_buffer
         i = i + 1
         done = is_iostat_end(ios) .or. (i == 3)
         if (done) exit
      end do

      empty_line = is_empty_string(read_buffer)

      ! we read the third line and it was not empty
      if ( (i == 3) .and. (.not. empty_line) ) then

         write(iunit, '(a)') trim(read_buffer)

         ! now read/write until we hit eof or 2 empty lines
         done = .false.
         do
            read (tplunit, '(a)', iostat = ios) read_buffer
            done = (empty_line .and. is_empty_string(read_buffer))
            done = done .or. is_iostat_end(ios)
            if (done) exit
            empty_line = is_empty_string(read_buffer)
            write(iunit, '(a)') trim(read_buffer)

         end do

         if (.not. empty_line) then
            ! last written line was not empty, add empty line
            write(iunit,'()')
         end if

      end if
!------------------------------


      close(tplunit)

   end if

   close(iunit, iostat=ierr)
   if ( ierr /= 0 ) then
     call sander_bomb('write_inpfile (qm2_extern_quick_module)', &
       'Error closing QUICK runfile after writing', &
       'Will quit now.')
   end if

   first_call = .false.

   call debug_exit_function( 'write_inpfile', module_name, quick_nml%debug )

  end subroutine write_inpfile

  subroutine read_results( outfile, nqmatoms, escf, dxyzqm, &
       nclatoms, dxyzcl, dipxyz, dipole, do_grad, debug )

    implicit none

    character(len=*), intent(in)  :: outfile
    integer, intent(in)           :: nqmatoms, nclatoms
    _REAL_, intent(out)           :: escf, dxyzqm(3,nqmatoms), &
                                     dxyzcl(3,nclatoms) ! dxyzcl returns containing the electric field at x,y,z
    _REAL_, intent(out)           :: dipxyz(3), dipole
    REAL, DIMENSION(3*nqmatoms)   :: grad
    logical, intent(in)           :: do_grad
    integer, intent(in)           :: debug

    _REAL_ :: self_energy = 0 ! Temporary variable to hold self energy of point charges
    integer :: ios, ixyz, iqm, i, j
    integer, parameter :: iunit = 351
    character(len=256) :: read_buffer


    call debug_enter_function( 'read_results', module_name, debug )

    open(iunit, file=outfile, status='old', iostat=ios)
    if ( ios /= 0 ) then
      call sander_bomb('read_results (qm2_extern_quick_module)', &
        'Error opening QUICK formatted checkpoint file '//outfile//' (expected in same dir as input file).', &
        'Will quit now')
    end if

    do ! START INFINITE LOOP OVER OUTPUT FILE LINES

       read (iunit, '(a)', iostat = ios) read_buffer
       ! End of file; data not found
       if (ios < 0) then
          call sander_bomb('read_results (qm2_extern_quick_module)', &
               'Error reading Quick output from file '//outfile, &
               '(Current total energy or gradient not found.)')
       end if

       ! looking and store SCF energy to escf
       if (read_buffer(1:20) == ' TOTAL ENERGY       ') then
          read (read_buffer(index(read_buffer,'=')+1:),*) escf
       end if

       ! look for QM gradients and read into dxyzqm
       if ( do_grad ) then

          if (read_buffer(1:21) == ' ANALYTICAL GRADIENT:' ) then
             do i = 1, 3
                read (iunit, '(a)', iostat = ios) read_buffer
             end do
             do iqm = 1, nqmatoms
                do ixyz = 1, 3
                   read (iunit, '(a)', iostat = ios) read_buffer
                   if (ios < 0) then
                      call sander_bomb('read_results (qm2_extern_quick_module)', &
                           'Error reading Quick output from file '//outfile, &
                           '(Problem reading QM gradient.)')
                   end if
                   read (read_buffer(24:),*) dxyzqm(ixyz,iqm)
                end do
             end do
          end if

       end if

       ! looking for QM dipoles
       if ((read_buffer(1:18) == '    DIPOLE (DEBYE)' )) then
          read (iunit, '(a)', iostat = ios) read_buffer
          read (iunit, '(a)', iostat = ios) read_buffer
          if (ios < 0) then
             call sander_bomb('read_results (qm2_extern_quick_module)', &
                  'Error reading Quick output from file '//outfile, &
                  '(Current dipole not found.)')
          end if
          read(read_buffer,*) dipxyz, dipole
          ! after dipole moment, we are done
          exit
       end if

       ! QM/MM electrostatic embedding
       if ( nclatoms > 0 .and. do_grad ) then

          ! Read MM gradients
          if (read_buffer(1:23) == ' POINT CHARGE GRADIENT:') then
             do i = 1, 3
                read (iunit, '(a)', iostat = ios) read_buffer
             end do
             ! Read into dxyzcl
             do i = 1, nclatoms
                do j = 1, 3
                   read (iunit, '(a)', iostat = ios) read_buffer
                   read(read_buffer(25:),*) dxyzcl(j,i)
                end do
             end do
          end if

      end if

    end do ! Reach end of file

    close(iunit)

    call debug_exit_function( 'read_results', module_name, debug )

  end subroutine read_results

end module qm2_extern_quick_module
