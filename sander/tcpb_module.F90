#include "../include/dprec.fh"
! QM/MM interface to TeraChem Protocol Buffers (TCPB) library
!
! Written by Vinicius Wilian D. Cruzeiro
! Date: October 2021

module tcpb_module

  use UtilitiesModule, only: Upcase
  use qmmm_module, only: qmmm_nml, qmmm_mpi

  implicit none

#ifdef MPI
  include 'mpif.h'
#endif

  private

  public :: get_tcpb_qmmm_forces, tcpb_finalize

#ifdef TCPB
  ! TCPB namelist
  type tcpb_nml_type
    character(len=256) :: tcfile, host 
    integer :: port
    character(len=20) :: basis
    character(len=20) :: method
    character(len=20) :: precision
    character(len=20) :: dftd
    character(len=20) :: cis
    _REAL_ :: threall
    _REAL_ :: convthre
    integer :: maxit
    integer :: dftgrid
    integer :: cisnumstates
    integer :: cistarget
    logical :: recycleinitguess ! Controls if we will try to use the wavefunction from previous step
    logical :: debug ! print debug output
  end type tcpb_nml_type

  type(tcpb_nml_type) :: tcpb_nml
#endif

contains

  ! Calculate the TCPB qmmm forces and energy
  subroutine get_tcpb_qmmm_forces(nstep, nqmatoms, qmcoords, qmtypes, &
       & nclatoms, clcoords, escf, dxyzqm, dxyzcl, id)

    use ElementOrbitalIndex, only : elementSymbol

    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS
    ! use constants, only: AU_TO_KCAL, A_TO_BOHRS

    implicit none

    integer, intent(in) :: nqmatoms
    _REAL_, intent(in)  :: qmcoords(3,nqmatoms)
    integer, intent(in) :: qmtypes(nqmatoms)
    integer, intent(in) :: nclatoms
    _REAL_, intent(in)  :: clcoords(4,nclatoms)     ! MM atom coordinates and charges
    _REAL_, intent(out) :: escf
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)
    character(len=3), intent(in) :: id
    integer, intent(in) :: nstep

    ! local variables
    _REAL_ :: tcpb_qmcoords(3*nqmatoms)
    character(len=5) :: qmatelsymb(nqmatoms)
    _REAL_ :: tcpb_mmcoords(3*nclatoms), tcpb_mmcharges(nclatoms)
    _REAL_ :: tcpb_dxyzqm(3*nqmatoms), tcpb_dxyzcl(3*nclatoms)
    integer :: i, ierr, globaltreatment
    logical, save :: first_call = .true.

    ! Initialize output values as zero
    escf = 0.0
    dxyzqm(:,:) = 0.0
    dxyzcl(:,:) = 0.0

#ifdef TCPB
    if (first_call) then
        first_call = .false.
        call tcpb_input_setting(nqmatoms, qmtypes, nclatoms, id)
    end if

    ! Only master thread works
    if (qmmm_mpi%commqmmm_master) then

      if (tcpb_nml%debug) then
        write(6,*) '>>> entered subroutine get_tcpb_qmmm_forces'
        write(6,*)
        write(6,*) ' TCPB input coordinates (Atomic number, and X, Y, and Z in A):'
        do i = 1, nqmatoms
            write(6,'(i3,3(f16.8,2x))') qmtypes(i), qmcoords(:,i)
        end do

        if (nclatoms > 0) then
          write(6,*)
          write(6,*) ' TCPB external point coordinates and charges (X, Y, and Z in A, and charge in a.u.):'
          do i = 1, nclatoms
            write(6,'(4(f16.8,2x))') clcoords(:,i)
          end do
        end if
      end if

      ! Prepare input/output variables for tc_compute_energy_gradient
      do i = 1, nqmatoms
        tcpb_qmcoords(3*(i-1)+1) = qmcoords(1,i) * CODATA08_A_TO_BOHRS
        tcpb_qmcoords(3*(i-1)+2) = qmcoords(2,i) * CODATA08_A_TO_BOHRS
        tcpb_qmcoords(3*(i-1)+3) = qmcoords(3,i) * CODATA08_A_TO_BOHRS
        tcpb_dxyzqm(3*(i-1)+1) = 0.0
        tcpb_dxyzqm(3*(i-1)+2) = 0.0
        tcpb_dxyzqm(3*(i-1)+3) = 0.0
      end do
      do i = 1, nclatoms
        tcpb_mmcoords(3*(i-1)+1) = clcoords(1,i) * CODATA08_A_TO_BOHRS
        tcpb_mmcoords(3*(i-1)+2) = clcoords(2,i) * CODATA08_A_TO_BOHRS
        tcpb_mmcoords(3*(i-1)+3) = clcoords(3,i) * CODATA08_A_TO_BOHRS
        tcpb_mmcharges(i) = clcoords(4,i)
        tcpb_dxyzcl(3*(i-1)+1) = 0.0
        tcpb_dxyzcl(3*(i-1)+2) = 0.0
        tcpb_dxyzcl(3*(i-1)+3) = 0.0
      end do
      do i = 1, nqmatoms
        qmatelsymb(i) = trim(elementSymbol(qmtypes(i)))//CHAR(0)
      end do
      ! globaltreatment handles how the initial guess for the wavefunction will be constructed
      ! globaltreatment = 1 means that a new condition will be done from scratch
      if (nstep<1 .or. .not. tcpb_nml%recycleinitguess) then
        globaltreatment = 1
      else
        globaltreatment = 0
      end if
      ! Compute energy and gradients with TCPB
      call tc_compute_energy_gradient(qmatelsymb,tcpb_qmcoords,nqmatoms,escf,tcpb_dxyzqm,tcpb_mmcoords,&
                                      tcpb_mmcharges,nclatoms,tcpb_dxyzcl,globaltreatment,ierr)
      if (ierr == 0 .and. tcpb_nml%debug) then
        write(6,*) ">>> Computed energy and gradients with success with tc_compute_energy_gradient."
        write(6,*)
      else if (ierr == 1) then
        write(6,*) "ERROR: Mismatch in the variables passed to the tc_compute_energy_gradient function!"
        call mexit(6, 1)
      else if (ierr /= 0) then
        write(6,*) "ERROR: Problem to compute energy and gradient with tc_compute_energy_gradient!"
        call mexit(6, 1)
      end if
      ! Convert output variables
      do i = 1, nqmatoms
        dxyzqm(1,i) = tcpb_dxyzqm(3*(i-1)+1)
        dxyzqm(2,i) = tcpb_dxyzqm(3*(i-1)+2)
        dxyzqm(3,i) = tcpb_dxyzqm(3*(i-1)+3)
      end do
      do i = 1, nclatoms
        dxyzcl(1,i) = tcpb_dxyzcl(3*(i-1)+1)
        dxyzcl(2,i) = tcpb_dxyzcl(3*(i-1)+2)
        dxyzcl(3,i) = tcpb_dxyzcl(3*(i-1)+3)
      end do

      ! Convert the Energy unit from Hartree to (kcal/mol)
      escf = escf * CODATA08_AU_TO_KCAL

      ! Convert the gradient unit from Hartree/Bohr to (kcal/mol)/A
      dxyzqm(:,:) = dxyzqm(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS

      if (nclatoms > 0) then
        ! Convert the gradient unit from Hartree/Bohr to (kcal/mol)/A
        dxyzcl(:,:) = dxyzcl(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
      end if

      if (tcpb_nml%debug) then
        write(6,*)
        write(6,*) ' TCPB energy (in kcal/mol):', escf
        write(6,*) ' TCPB gradients in the QM region (in (kcal/mol)/A):'
        do i = 1, nqmatoms
            write(6,*) dxyzqm(:,i)
        end do
        write(6,*) ' TCPB gradients in the MM region (in (kcal/mol)/A):'
        do i = 1, nclatoms
            write(6,*) dxyzcl(:,i)
        end do
        write(6,*)
        write(6,*) '<<< leaving subroutine get_tcpb_qmmm_forces'
      end if

    end if ! qmmm_mpi%commqmmm_master

#else
     call sander_bomb('get_tcpb_qmmm_forces','TCPB is not enabled', &
           'Amber needs to be installed with the -terachem', &
           ' flag (if using the configure script) or BUILD_TCPB=TRUE (if using CMake).')
#endif

  end subroutine get_tcpb_qmmm_forces

! ----------------
! Private routines
! ----------------

#ifdef TCPB
  ! Set up input data for the TCPB library
  subroutine tcpb_input_setting(nqmatoms, qmtypes, nclatoms, id)

    use ElementOrbitalIndex, only : elementSymbol

    implicit none

    integer, intent(in) :: nqmatoms, nclatoms
    integer, intent(in) :: qmtypes(nqmatoms)
    character(len=3), intent(in) :: id
    character(len=5) :: qmatelsymb(nqmatoms)
    integer :: ierr, i
    logical :: exstfl
    character(len=256) :: tcfile
    integer, parameter :: iurun = 10

    ! Only master reads the mdin input file
    if (qmmm_mpi%commqmmm_master) then
      ! read input from mdin file
      call read_tcpb_nml(tcpb_nml)

      ! Check if input flags are ok
      if (tcpb_nml%tcfile .eq. '' .and. tcpb_nml%method .eq. '') then
        write(6, '(a,a)') 'ERROR: Please specify the tcfile flag or method and basis',&
                          ' flags for TCPB in the tc namelist!'
        call mexit(6,1)
      end if
      if (tcpb_nml%method .ne. '' .and. tcpb_nml%basis .eq. '') then
        write(6, '(a,a)') 'ERROR: You specified a method, but not a basis set for TCPB. Please set it using the',&
                          ' basis flag in the tc namelist!'
        call mexit(6,1)
      end if
      if (tcpb_nml%host .eq. '') then
        write(6, '(a)') 'ERROR: Please specify the host for TCPB in the tc namelist!'
        call mexit(6,1)
      end if
      if (tcpb_nml%port < 0) then
        write(6, '(a)') 'ERROR: Please specify the port for TCPB in the tc namelist!'
        call mexit(6,1)
      end if

      ! Constructing the TeraChem input file, if needed
      if (tcpb_nml%method .ne. '') then
        if (tcpb_nml%tcfile .eq. '') then
          tcfile = "terachem"//trim(id)//".inp"
        else
          tcfile = tcpb_nml%tcfile
        end if
        write(6, '(2a)') '| Writing TeraChem input file to ', trim(tcfile)
        open(iurun, file=tcfile, iostat=ierr)
        if ( ierr > 0 ) then
          call sander_bomb('tcpb_input_setting', &
              'Error opening TeraChem input file '//trim(tcfile)//' for writing', &
              'Will quit now')
        end if
        write(iurun, '(a,/,a)')'# Run using SANDER TCPB interface for TeraChem','#'
        write(iurun, '(2a)')       'method       ', trim(tcpb_nml%method)
        write(iurun, '(2a)')       'basis        ', trim(tcpb_nml%basis)
        write(iurun, '(2a)')       'precision    ', trim(tcpb_nml%precision)
        write(iurun, '(a,E22.16)') 'threall      ', tcpb_nml%threall
        write(iurun, '(a,E22.16)') 'convthre     ', tcpb_nml%convthre
        write(iurun, '(2a)')       'dftd         ', trim(tcpb_nml%dftd)
        write(iurun, '(a,i4)')     'maxit        ', tcpb_nml%maxit
        write(iurun, '(a,i3)')     'dftgrid      ', tcpb_nml%dftgrid
        write(iurun, '(2a)')       'cis          ', trim(tcpb_nml%cis)
        if ( trim(tcpb_nml%cis) == 'yes' ) then
          write(iurun, '(a,i3)')     'cisnumstates ', tcpb_nml%cisnumstates
          write(iurun, '(a,i3)')     'cistarget    ', tcpb_nml%cistarget
        end if
        write(iurun, '(a,i3)')     'charge       ', qmmm_nml%qmcharge
        write(iurun, '(a,i3)')     'spinmult     ', qmmm_nml%spin
        flush(iurun)
        close(iurun)
      else
        tcfile = tcpb_nml%tcfile
        inquire(file=tcfile, exist=exstfl)
        if (.not. exstfl) then
          call sander_bomb('tcpb_input_setting', &
            'TeraChem input file '//trim(tcfile)//' does not exist!', &
            'Will quit now.')
        end if
        write(6, '(2a)') '| Reading TeraChem input file from ', trim(tcfile)
      end if

      ! Attempts to connect to the TeraChem server
      call tc_connect(trim(tcpb_nml%host)//CHAR(0), tcpb_nml%port, ierr)
      if (ierr == 0 .and. tcpb_nml%debug) then
        write(6,*) ">>> Successfully connected to TeraChem server using host ", &
        trim(tcpb_nml%host), " and port ", tcpb_nml%port, "."
        write(6,*)
      else if (ierr == 2) then
        write (6,*) "ERROR: Connection to TeraChem server succeed in the tc_connect function, but the server is not available!"
        call mexit(6,1)
      else if (ierr /= 0) then
        write (6,*) "ERROR: Connection to TeraChem server failed in the tc_connect function!"
        call mexit(6,1)
      end if

      ! Setup TeraChem
      do i = 1, nqmatoms
        qmatelsymb(i) = trim(elementSymbol(qmtypes(i)))//CHAR(0)
      end do
      call tc_setup(trim(tcfile)//CHAR(0),qmatelsymb,nqmatoms,ierr)
      if (ierr == 0 .and. tcpb_nml%debug) then
        write(6,*) ">>> TeraChem setup completed with success."
        write(6,*)
      else if (ierr == 1) then
        write(6,*) "ERROR: No options read from TeraChem input file in the tc_setup function!"
        write(6,*)
      else if (ierr /= 0) then
        write(6,*) "ERROR: Failed to setup TeraChem in the tc_setup."
        write(6,*)
      end if

      if (tcpb_nml%debug) then
        call print_tcpb_nml(tcpb_nml)
        write(6,*) ''
        write(6,*) '<<< leaving subroutine tcpb_input_setting'
      end if
    ! End master
    end if
  end subroutine tcpb_input_setting

  ! Read TCPB namelist from mdin file
  subroutine read_tcpb_nml(tcpb_nml)

    implicit none
    type(tcpb_nml_type) :: tcpb_nml
    ! local variables
    integer, parameter :: iu_mdin = 5  ! assume mdin file connected to unit 5
    integer :: ierr
    logical :: is_open

    ! namelist variables
    character(len=256) :: tcfile, host
    integer :: port
    character(len=20) :: basis
    character(len=20) :: method
    character(len=20) :: precision
    character(len=20) :: dftd
    character(len=20) :: cis
    _REAL_ :: threall
    _REAL_ :: convthre
    integer :: maxit
    integer :: dftgrid
    integer :: cisnumstates
    integer :: cistarget
    integer :: recycleinitguess
    integer :: debug ! print debug output

    namelist /tc/ tcfile, host, port, method, basis, precision, dftd, cis, threall,&
                  convthre, maxit, dftgrid, cisnumstates, cistarget, recycleinitguess, debug

    ! Default namelist variable values
    tcfile          = ''
    host            = ''
    port            = -1
    basis           = ''
    method          = ''
    dftd            = 'no'
    precision       = 'mixed'
    cis             = 'no'
    threall         = 1.0d-11
    convthre        = 3.0d-05
    maxit           = 100
    dftgrid         = 1
    cisnumstates    = 1
    cistarget       = 0
    recycleinitguess= 1 ! Controls if we will try to use the wavefunction from previous step
    debug           = 0

    ! Read namelist
    inquire(unit=iu_mdin, opened=is_open)
    if ( .not. is_open) then
      call sander_bomb('read_tcpb_nml', &
        'mdin file not connected to unit 5', &
        'Stopping now.')
    end if
    rewind(unit=iu_mdin)
    read(unit=iu_mdin, nml=tc, iostat=ierr)
    if ( ierr /= 0 ) then
        call sander_bomb('read_tcpb_nml', &
          '&tc namelist read error', &
          'Please check your input.')
    end if

    ! Setting variables from what was read in the namelist
    tcpb_nml%tcfile       = tcfile
    tcpb_nml%host         = host
    tcpb_nml%port         = port
    tcpb_nml%method       = method
    tcpb_nml%basis        = basis
    tcpb_nml%precision    = precision
    tcpb_nml%dftd         = dftd
    tcpb_nml%cis          = cis
    tcpb_nml%threall      = threall
    tcpb_nml%convthre     = convthre
    tcpb_nml%maxit        = maxit
    tcpb_nml%dftgrid      = dftgrid
    tcpb_nml%cisnumstates = cisnumstates
    tcpb_nml%cistarget    = cistarget
    if ( recycleinitguess == 0) then
       tcpb_nml%recycleinitguess = .false.
    else if ( recycleinitguess > 0) then
       tcpb_nml%recycleinitguess = .true.
    else
       call sander_bomb('read_tcpb_nml', &
          '&tc recycleinitguess read error', &
          'Please check your input.')
    end if
    if ( debug == 0) then
       tcpb_nml%debug = .false.
    else if ( debug > 0) then
       tcpb_nml%debug = .true.
    else
       call sander_bomb('read_tcpb_nml', &
          '&tc debug read error', &
          'Please check your input.')
    end if

    if (tcpb_nml%debug .and. qmmm_mpi%commqmmm_master) then
      write(6,*) '<<< leaving subroutine read_tcpb_nml'
    end if

  end subroutine read_tcpb_nml

  ! Print the TCPB namelist
  subroutine print_tcpb_nml(self)

    implicit none
    type(tcpb_nml_type), intent(in) :: self

    write(6,'(/,a)')      '     ======== TCPB settings ======== '
    write(6,'(a,a)')      ' tcfile                     : ', trim(self%tcfile)
    write(6,'(a,a)')      ' host                       : ', trim(self%host)
    write(6,'(a,i0)')     ' port                       : ', self%port
    write(6,'(a,a)')      ' method                     : ', trim(self%method)
    write(6,'(a,a)')      ' basis                      : ', trim(self%basis)
    write(6,'(a,i0)')     ' charge (from qmmm namelist): ', qmmm_nml%qmcharge
    write(6,'(a,i0)')     ' mult   (from qmmm namelist): ', qmmm_nml%spin
    write(6,'(a,a)')      ' precision                  : ', trim(self%precision)
    write(6,'(a,a)')      ' dftd                       : ', trim(self%dftd)
    write(6,'(a,es10.2)') ' threall                    : ', self%threall
    write(6,'(a,es10.2)') ' convthre                   : ', self%convthre
    write(6,'(a,i0)')     ' maxit                      : ', self%maxit
    write(6,'(a,i0)')     ' dftgrid                    : ', self%dftgrid
    write(6,'(a,a)')      ' cis                        : ', trim(self%cis)
    write(6,'(a,i0)')     ' cisnumstates               : ', self%cisnumstates
    write(6,'(a,i0)')     ' cistarget                  : ', self%cistarget
    write(6,'(a,l)')      ' recycleinitguess           : ', self%recycleinitguess
    write(6,'(a,l)')      ' debug                      : ', self%debug

  end subroutine print_tcpb_nml
#endif

  subroutine tcpb_finalize()

    IMPLICIT NONE

    integer :: ierr

#ifdef TCPB
    ! Finalizes variables on the TeraChem side
    call tc_finalize()
#endif

  end subroutine

end module tcpb_module
