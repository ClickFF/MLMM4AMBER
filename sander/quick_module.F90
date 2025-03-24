#include "../include/dprec.fh"
! QM/MM interface to Quick library
!
! Written by Vinicius Cruzeiro, Andreas W. Goetz
! Date: July 2020

module quick_module

  use UtilitiesModule, only: Upcase
  use qmmm_module, only: qmmm_nml, qmmm_mpi
#ifdef QUICK
  use quick_api_module, only : setQuickJob, getQuickEnergyGradients, deleteQuickJob
#ifdef MPI
  use quick_api_module, only : setQuickMPI
#endif
#endif

  implicit none

#ifdef MPI
  include 'mpif.h'
#endif

  private

  public :: get_quick_qmmm_forces, quick_finalize

#ifdef QUICK
  ! QUICK namelist
  type quick_nml_type
    character(len=256) :: keywords, outfprefix
    character(80) :: method, basis
    integer :: scf_cyc
    _REAL_  :: denserms, intcutoff, xccutoff, gradcutoff, basiscutoff
    logical :: debug ! print debug output
    logical :: reuse_dmx
  end type quick_nml_type

  type(quick_nml_type) :: quick_nml
#endif

  !
  ! Amber-consistent conversion factors
  ! Parts of the Ewald are performed by sander, cew, and quick;
  ! they need to use consistent conversions in order for real
  ! and reciprocal space interactions to properly cancel.
  !
  ! NOTE: These constants are probably in cewmod.f90 (module cewmod)
  !       and we should use them from there
  !       but we cannot because cewmod uses quick_module (sigh)
  _REAL_,parameter :: AMBERELE = 18.2223d0
  _REAL_,parameter :: BOHRS_TO_A = 0.529177249d0
  _REAL_,parameter :: AU_PER_AMBER_CHARGE = 1.d0 / AMBERELE
  _REAL_,parameter :: AU_PER_AMBER_ANGSTROM = 1.d0 / BOHRS_TO_A
  _REAL_,parameter :: AU_PER_AMBER_KCAL_PER_MOL = AU_PER_AMBER_CHARGE *  AU_PER_AMBER_CHARGE / AU_PER_AMBER_ANGSTROM
  _REAL_,parameter :: AU_PER_AMBER_FORCE = AU_PER_AMBER_KCAL_PER_MOL / AU_PER_AMBER_ANGSTROM

contains

  ! Calculate the QUICK qmmm force and energy
  subroutine get_quick_qmmm_forces(nqmatoms, qmcoords, qmtypes, &
       & nclatoms, clcoords, escf, dxyzqm, dxyzcl, use_cew)

    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS

    implicit none

    integer, intent(in) :: nqmatoms
    _REAL_, intent(in)  :: qmcoords(3,nqmatoms)
    integer, intent(in) :: qmtypes(nqmatoms)
    integer, intent(in) :: nclatoms
    _REAL_, intent(in)  :: clcoords(4,nclatoms)     ! MM atom coordinates
    _REAL_, intent(out) :: escf
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)
    logical, intent(in) :: use_cew

    ! local variables
    integer :: i, ierr
    logical, save :: first_call = .true.

    ! Define local conversion factors
    ! Use CODATA08 values without CEw
    ! Use Amber internal conversion factors if using CEw to be consistent across Amber functions
    _REAL_ :: AU_TO_KCAL_PER_MOL = CODATA08_AU_TO_KCAL
    _REAL_ :: AU_TO_KCAL_PER_MOL_PER_ANGSTROM = CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS

    ! Assign CEw consistent conversion factors
    if (use_cew) then
       AU_TO_KCAL_PER_MOL = 1.0D0 / AU_PER_AMBER_KCAL_PER_MOL
       AU_TO_KCAL_PER_MOL_PER_ANGSTROM = 1.0D0 / AU_PER_AMBER_FORCE
    end if

    ! Initialize output values as zero
    escf = 0.0d0
    dxyzqm(:,:) = 0.0d0
    dxyzcl(:,:) = 0.0d0

#ifdef QUICK
    if (quick_nml%debug .and. qmmm_mpi%commqmmm_master) then
      write(6,*) '>>> entered subroutine get_quick_qmmm_forces'
    end if

    if (first_call) then
        first_call = .false.
        call quick_input_setting(nqmatoms, qmtypes, nclatoms)
    end if

    if (quick_nml%debug .and. qmmm_mpi%commqmmm_master) then
       write(6,*)
       write(6,*) ' QUICK input coordinates (Atomic number, and X, Y, and Z in A):'
       do i = 1, nqmatoms
          write(6,'(i3,3(f16.8,2x))') qmtypes(i), qmcoords(:,i)
       end do
    end if

    if (nclatoms > 0 .and. quick_nml%debug .and. qmmm_mpi%commqmmm_master) then
      write(6,*)
      write(6,*) ' QUICK external point coordinates and charges (X, Y, and Z in A, and charge in a.u.):'
      do i = 1, nclatoms
         write(6,'(4(f16.8,2x))') clcoords(:,i)
      end do
    end if
    
    ! Compute energy and gradients with Quick
    call getQuickEnergyGradients(qmcoords, nclatoms, clcoords, escf, dxyzqm, dxyzcl, ierr)
    if ( ierr /= 0 ) then
        call sander_bomb('Error getting energy and gradient with Quick in getQuickEnergyGradients at quick_module.F90')
     end if

#ifdef MPI
    ! Quick returns the same energies and forces to all threads,
    ! but these are going to be added up by the QM/MM module later on.
    ! For this reason we set the energy and force to zero if we are not
    ! the master task
    if ( .not. qmmm_mpi%commqmmm_master ) then
       escf = 0.0d0
       dxyzqm(:,:) = 0.0d0
       if (nclatoms > 0) dxyzcl(:,:) = 0.0d0
    end if
#endif

    ! Convert the Energy unit from Hartree to (kcal/mol)
    escf = escf * AU_TO_KCAL_PER_MOL

    ! Convert the gradient unit from Hartree/Bohr to (kcal/mol)/A
    dxyzqm(:,:) = dxyzqm(:,:) * AU_TO_KCAL_PER_MOL_PER_ANGSTROM

    if (nclatoms > 0) then
       ! Convert the gradient unit from Hartree/Bohr to (kcal/mol)/A
       dxyzcl(:,:) = dxyzcl(:,:) * AU_TO_KCAL_PER_MOL_PER_ANGSTROM
    end if

    if (quick_nml%debug .and. qmmm_mpi%commqmmm_master) then
       write(6,*)
       write(6,*) ' QUICK energy (in kcal/mol):', escf
       write(6,*) ' QUICK gradients in the QM region (in (kcal/mol)/A):'
       do i = 1, nqmatoms
          write(6,*) dxyzqm(:,i)
       end do
       write(6,*) ' QUICK gradients in the MM region (in (kcal/mol)/A):'
       do i = 1, nclatoms
          write(6,*) dxyzcl(:,i)
       end do
       write(6,*)
       write(6,*) '<<< leaving subroutine get_quick_qmmm_forces'
    end if

#else
     call sander_bomb('get_quick_qmmm_forces','QUICK is not enabled', &
           'Check your installation or reconfigure with the -quick option.')
#endif

  end subroutine get_quick_qmmm_forces

! ----------------
! Private routines
! ----------------

#ifdef QUICK
  ! Set up input data for the QUICK library
  subroutine quick_input_setting(nqmatoms, qmtypes, nclatoms)

    implicit none

    integer, intent(in) :: nqmatoms, nclatoms
    integer, intent(in) :: qmtypes(nqmatoms)
    character(len=256) :: keys
    integer :: ierr
    character(len=4) :: dftkey

#ifdef MPI
    ! Initialize Quick MPI
    call setQuickMPI(qmmm_mpi%mytaskid,qmmm_mpi%numthreads, ierr)
    if ( ierr /= 0 ) then
        call sander_bomb('Error setting up Quick MPI in setQuickMPI at quick_module.F90')
    end if
#endif

    ! Only master reads the mdin input file
    if (qmmm_mpi%commqmmm_master) then
      ! read input from mdin file
      call read_quick_nml(quick_nml)

      ! Check if input flags are ok
      if (quick_nml%keywords .eq. '' .and. quick_nml%method .eq. '') then
        write(6, '(a,a)') 'ERROR: Please specify the keywords flag or method and basis',&
                          ' flags for QUICK in the quick namelist!'
        call mexit(6, 1)
      end if
      if (quick_nml%keywords .eq. '' .and. quick_nml%method .ne. '' .and. quick_nml%basis .eq. '') then
        write(6, '(a,a)') 'ERROR: Please specify the basis set for QUICK using the',&
                          ' basis flag in the quick namelist!'
        call mexit(6, 1)
      end if

      ! Constructing the keywords flag, if needed
      if (quick_nml%keywords .eq. '') then
         if (Upcase(quick_nml%method) .eq. 'HF') then
            dftkey = ' HF '
         else
            dftkey = 'DFT '
         end if
         write(keys, '(5a,i0,4(a,f0.10),a,2(i0,a))') &
              dftkey, trim(quick_nml%method), &
              ' BASIS=', trim(quick_nml%basis), &
              ' SCF=', quick_nml%scf_cyc, &
              ' DENSERMS=', quick_nml%denserms, &
              ' CUTOFF=', quick_nml%intcutoff, &
              ' GRADCUTOFF=', quick_nml%gradcutoff, &
              ' BASISCUTOFF=', quick_nml%basiscutoff, &
              ' CHARGE=', qmmm_nml%qmcharge, &
              ' MULT=', qmmm_nml%spin, &
              ' GRADIENT'
         if (nclatoms > 0) then
            write(keys, '(a,a)') trim(keys), ' EXTCHARGES'
         end if
      else
         write(keys, '(a)') trim(quick_nml%keywords)
      end if
    ! End master
    end if

#ifdef MPI
    ! Broadcast from master to the other threads
    call mpi_bcast(quick_nml%outfprefix, 256, MPI_CHARACTER, 0, qmmm_mpi%commqmmm, ierr)
    call mpi_bcast(keys, 256, MPI_CHARACTER, 0, qmmm_mpi%commqmmm, ierr)
    call mpi_bcast(quick_nml%reuse_dmx, 1, MPI_LOGICAL, 0, qmmm_mpi%commqmmm, ierr)
#endif

    ! Initialize QUICK
    call setQuickJob(quick_nml%outfprefix, keys, nqmatoms, qmtypes, quick_nml%reuse_dmx, ierr)
    if ( ierr /= 0 ) then
        call sander_bomb('Error setting up Quick in setQuickJob at quick_module.F90')
    end if

    if (quick_nml%debug .and. qmmm_mpi%commqmmm_master) then
      call print_quick_nml(quick_nml)
      write(6,*) ''
      write(6,*) '<<< leaving subroutine quick_input_setting'
    end if
  end subroutine quick_input_setting

  ! Read QUICK namelist from mdin file
  subroutine read_quick_nml(quick_nml)

    implicit none
    type(quick_nml_type) :: quick_nml
    ! local variables
    integer, parameter :: iu_mdin = 5  ! assume mdin file connected to unit 5
    integer :: ierr
    logical :: is_open

    ! namelist variables
    character(len=256) :: keywords, outfprefix
    character(80) :: method, basis
    integer :: scf_cyc, debug, reuse_dmx
    _REAL_  :: denserms, intcutoff, xccutoff, gradcutoff, basiscutoff

    namelist /quick/ outfprefix, method, basis, keywords, scf_cyc, denserms, &
                     intcutoff, xccutoff, gradcutoff, basiscutoff, debug, reuse_dmx

    ! Default namelist variable values
    outfprefix = 'quick'
    method = 'BLYP'
    basis = '6-31G'
    keywords = ''
    scf_cyc = 200
    denserms = 1.0E-6
    intcutoff = 1.0E-8
    xccutoff = 1.0E-8
    gradcutoff = 1.0E-7
    basiscutoff = 1.0E-6
    debug = 0
    reuse_dmx=1

    ! Read namelist
    inquire(unit=iu_mdin, opened=is_open)
    if ( .not. is_open) then
      call sander_bomb('read_quick_nml', &
        'mdin file not connected to unit 5', &
        'Stopping now.')
    end if
    rewind(unit=iu_mdin)
    read(unit=iu_mdin, nml=quick, iostat=ierr)
    if ( ierr /= 0 ) then
        call sander_bomb('read_quick_nml', &
          '&quick namelist read error', &
          'Please check your input.')
    end if

    ! Setting variables from what was read in the namelist
    quick_nml%outfprefix = outfprefix
    quick_nml%method     = method
    quick_nml%basis      = basis
    quick_nml%keywords   = keywords
    quick_nml%scf_cyc    = scf_cyc
    quick_nml%denserms   = denserms
    quick_nml%intcutoff  = intcutoff
    quick_nml%xccutoff   = xccutoff
    quick_nml%basiscutoff  = basiscutoff
    quick_nml%gradcutoff   = gradcutoff
    quick_nml%reuse_dmx  = reuse_dmx
    
    if ( debug == 0) then
       quick_nml%debug = .false.
    else if ( debug > 0) then
       quick_nml%debug = .true.
    else
       call sander_bomb('read_quick_nml', &
          '&quick debug read error', &
          'Please check your input.')
    end if

    if (quick_nml%debug .and. qmmm_mpi%commqmmm_master) then
      write(6,*) '<<< leaving subroutine read_quick_nml'
    end if

  end subroutine read_quick_nml

  ! Print the QUICK namelist
  subroutine print_quick_nml(self)

    implicit none
    type(quick_nml_type), intent(in) :: self

    write(6,'(/,a)') '     ======== QUICK settings ======== '
    write(6,'(a,a)')  ' outfprefix                 : ', trim(self%outfprefix)
    write(6,'(a,a)')  ' method                     : ', trim(self%method)
    write(6,'(a,a)')  ' basis                      : ', trim(self%basis)
    write(6,'(a,i0)') ' charge (from qmmm namelist): ', qmmm_nml%qmcharge
    write(6,'(a,i0)') ' mult   (from qmmm namelist): ', qmmm_nml%spin
    write(6,'(a,i0)') ' scf_cyc                    : ', self%scf_cyc
    write(6,'(a,f15.10)') ' denserms               : ', self%denserms
    write(6,'(a,f15.10)') ' intcutoff              : ', self%intcutoff
    write(6,'(a,f15.10)') ' xccutoff               : ', self%xccutoff
    write(6,'(a,f15.10)') ' basiscutoff            : ', self%basiscutoff
    write(6,'(a,f15.10)') ' gradcutoff             : ', self%gradcutoff
    write(6,'(a,l)')  ' reuse_dmx                  : ', self%reuse_dmx
    write(6,'(a,a)')  ' keywords                   : ', trim(self%keywords)
    write(6,'(a,l)')  ' debug                      : ', self%debug

  end subroutine print_quick_nml
#endif

  subroutine quick_finalize()

    IMPLICIT NONE

    integer :: ierr

    ! Finalize Quick
#ifdef QUICK
    call deleteQuickJob(ierr)
    if ( ierr /= 0 ) then
        call sander_bomb('Error ending Quick in deleteQuickJob at quick_module.F90')
    end if
#endif

  end subroutine

end module quick_module
