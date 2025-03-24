!*******************************************************************************
!
! Module: external_module
!
! Description: This module houses the capabilities of calling external libraries
!              to compute energies and forces, instead of using the force field ones.
!
!              Adapted here by Vinicius Wilian D. Cruzeiro and Delaram Ghoreishi
!
!*******************************************************************************

#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

module external_module

  use UtilitiesModule, only: Upcase
#ifdef QUICK
  use quick_api_module, only : setQuickJob, getQuickEnergyGradients, deleteQuickJob
#endif

  private

  character(256), public :: extprog, json, keywords, outfprefix
  character(80), public :: method, basis
  integer, public :: charge, mult

  ! The active namelist:

  private       :: extpot

  namelist /extpot/      extprog, json, outfprefix, method, basis, charge, mult, keywords

  public external_init, gb_external, pme_external, external_cleanup

  contains

  subroutine external_init(ix,ih,xx)

    use memory_module,   only : natom, m04, nres, m02, i02, i100, lmass
    use qmmm_module,     only : get_atomic_number

    IMPLICIT NONE

    integer :: i, ifind, cnt, ierr
    _REAL_ :: coord(3*natom)
    character(len=5), allocatable :: monomers(:)
    integer, allocatable :: nats(:)
    character(len=5), allocatable :: at_name(:)
    integer :: nmon, val
    integer :: ix(*)
    character(len=4) :: ih(*)
    _REAL_ :: xx(*)
    integer :: atomicnumber(natom)
    logical :: isvalid, errFlag
#ifdef QUICK
    character(256) :: keys
#endif

!    write (6,*) " -> Printing atom names"
!    do i=1,natom
!      write (6,*) "  * ", ih(m04+i-1)
!    end do
!    write (6,*) " -> Printing residue names"
!    do i=1,nres
!      write (6,*) "  * ", ih(m02+i-1)
!    end do
!    write (6,*) " -> Printing residue pointers"
!    do i=1,nres
!      write (6,*) "  * ", ix(i02+i-1)
!    end do

    extprog = ''
    json = ''
    outfprefix = 'quick'
    method = ''
    basis = ''
    charge = 0
    mult = 1
    keywords = ''

    ! Read input in namelist format:

    rewind(5)                  ! Insurance against maintenance mods.

    call nmlsrc('extpot', 5, ifind)

    if (ifind .ne. 0) then        ! Namelist found. Read it:
      read(5, nml = extpot)
    else                          ! No namelist was present,
      write(6, '(a)') 'ERROR: Could not find the namelist extpot in the mdin file'
      call mexit(6, 1)
    end if

    if (upcase(extprog) .eq. 'KMMD') then
        call kmmd_init(trim(json)//CHAR(0))
        return
    end if
    
    
    do i=1, natom
      if(ix(i100) .eq. 1) then
        atomicnumber(i) = ix(i100+i)
      else
        call get_atomic_number(ih(m04+i-1), xx(lmass+i-1), atomicnumber(i), errFlag)
      end if
    end do

    ! For MBX
    if (Upcase(extprog) .eq. 'MBX') then
#ifdef MBX
      nmon = nres

      allocate(nats(nmon),monomers(nmon),at_name(natom))

      cnt = 1
      do i=1,nmon
         isvalid = .True.
         if (ih(m02+i-1) == "WAT") then
            val = 3
            monomers(i) = "h2o"
            at_name(cnt+0) = "O"
            if (atomicnumber(cnt+0) .ne. 8) isvalid = .False.
            at_name(cnt+1) = "H"
            if (atomicnumber(cnt+1) .ne. 1) isvalid = .False.
            at_name(cnt+2) = "H"
            if (atomicnumber(cnt+2) .ne. 1) isvalid = .False.
            cnt = cnt + val
         else if (ih(m02+i-1) == "N2O") then
            val = 7
            monomers(i) = "n2o5"
            at_name(cnt+0) = "O"
            if (atomicnumber(cnt+0) .ne. 8) isvalid = .False.
            at_name(cnt+1) = "N"
            if (atomicnumber(cnt+1) .ne. 7) isvalid = .False.
            at_name(cnt+2) = "N"
            if (atomicnumber(cnt+2) .ne. 7) isvalid = .False.
            at_name(cnt+3) = "O"
            if (atomicnumber(cnt+3) .ne. 8) isvalid = .False.
            at_name(cnt+4) = "O"
            if (atomicnumber(cnt+4) .ne. 8) isvalid = .False.
            at_name(cnt+5) = "O"
            if (atomicnumber(cnt+5) .ne. 8) isvalid = .False.
            at_name(cnt+6) = "O"
            if (atomicnumber(cnt+6) .ne. 8) isvalid = .False.
            cnt = cnt + val
         else
            write(6, '(a,a,a)') 'ERROR: The residue ',ih(m02+i-1),' is not recognized by MBX!'
            call mexit(6, 1)
         end if
         if (val == ix(i02+i)-ix(i02+i-1)) then
            nats(i) = val
         else
            write(6, '(a,a,a)') 'ERROR: The number of atoms in residue ',ih(m02+i-1),' does not match the expected by MBX!'
            call mexit(6, 1)
         end if
         if (.not. isvalid) then
            write(6, '(a,a,a)') 'ERROR: The order or type of the atoms in residue ',ih(m02+i-1),&
                                    ' does not match the expected by MBX!'
            call mexit(6, 1)
         end if
      end do

      do i =1, natom
        at_name(i)=trim(at_name(i))//CHAR(0)
      end do
      do i=1,nmon
        monomers(i) = trim(monomers(i))//CHAR(0)
      end do

      if (json /= '') then
        call initialize_system(coord, nats, at_name, monomers, nmon, trim(json)//CHAR(0))
      else
        call initialize_system(coord, nats, at_name, monomers, nmon)
      end if
#endif
    ! For Quick
    else if (Upcase(extprog) .eq. 'QUICK') then
#ifdef QUICK
      ! Check if input flags are ok
      if (keywords .eq. '' .and. method .eq. '') then
        write(6, '(a,a)') 'ERROR: Please specify the keywords flag or method and basis',&
                              ' flags for QUICK in the extpot namelist!'
        call mexit(6, 1)
      end if
      if (keywords .eq. '' .and. method .ne. '' .and. basis .eq. '') then
        write(6, '(a,a)') 'ERROR: Please specify the basis set for QUICK using the',&
                              ' basis flag in the extpot namelist!'
        call mexit(6, 1)
      end if

      ! Constructing the keywords flag, if needed
      if (keywords .eq. '') then
        if (Upcase(method) .eq. 'HF') then
          write(keys, '(4a,2(i0,a))') trim(method), ' BASIS=', trim(basis),&
                                      ' CHARGE=', charge, ' MULT=', mult,' GRADIENT'
        else
          write(keys, '(5a,2(i0,a))') 'DFT ', trim(method), ' BASIS=', trim(basis),&
                                      ' CHARGE=', charge, ' MULT=', mult,' GRADIENT'
        end if
      else
        write(keys, '(a)') trim(keywords)
      end if
      call setQuickJob(outfprefix, keys, natom, atomicnumber, .true., ierr)
      if ( ierr /= 0 ) then
        write(6, '(a)') 'ERROR: setting up Quick in setQuickJob at external.F90'
        call mexit(6, 1)
      end if
#endif
    else
      write(6, '(a,a,a)') 'ERROR: External program ',trim(extprog),&
             ' is not valid! Please set a valid value in the extprog flag'
      call mexit(6, 1)
    end if

  end subroutine

  subroutine gb_external(crd, frc, pot_ene)

    use memory_module,   only : natom
    use constants, only :  EV_TO_KCAL, AU_TO_EV, BOHRS_TO_A

    IMPLICIT NONE

    _REAL_   ::  crd(*)
    _REAL_   ::  frc(*)
    _REAL_   ::  pot_ene
    integer  ::  i, i3
    _REAL_   ::  coord(3*natom), grads(3*natom)
#ifdef QUICK
    integer  ::  ierr
    _REAL_, dimension(3,natom) :: crd_quick, frc_quick
    _REAL_, dimension(:,:), pointer :: xc_coord       => null()
    _REAL_, dimension(:,:), pointer :: ptchgGrad      => null()
#endif

    ! For MBX
    if (Upcase(extprog) .eq. 'MBX') then
#ifdef MBX
      do i = 1, natom
        i3=3*(i-1)
        coord(3*(i-1)+1) = crd(i3+1)
        coord(3*(i-1)+2) = crd(i3+2)
        coord(3*(i-1)+3) = crd(i3+3)
      end do

      call get_energy_g(coord, natom, pot_ene, grads)

      do i = 1, natom
        i3=3*(i-1)
        frc(i3+1) = frc(i3+1) - grads(3*(i-1)+1)
        frc(i3+2) = frc(i3+2) - grads(3*(i-1)+2)
        frc(i3+3) = frc(i3+3) - grads(3*(i-1)+3)
      end do
#endif
    ! For Quick
    else if (Upcase(extprog) .eq. 'QUICK') then
#ifdef QUICK
      do i = 1, natom
        crd_quick(1,i) = crd(3*(i-1)+1)
        crd_quick(2,i) = crd(3*(i-1)+2)
        crd_quick(3,i) = crd(3*(i-1)+3)
        frc_quick(1,i) = 0.0
        frc_quick(2,i) = 0.0
        frc_quick(3,i) = 0.0
      end do
      call getQuickEnergyGradients(crd_quick, 0, xc_coord, pot_ene, frc_quick, ptchgGrad, ierr)
      if ( ierr /= 0 ) then
        write(6, '(a)') 'ERROR: getting energy and gradient with Quick in getQuickEnergyGradients at external.F90'
        call mexit(6, 1)
      end if
      do i = 1, natom
        frc(3*(i-1)+1) = - frc_quick(1,i) * AU_TO_EV * EV_TO_KCAL / BOHRS_TO_A
        frc(3*(i-1)+2) = - frc_quick(2,i) * AU_TO_EV * EV_TO_KCAL / BOHRS_TO_A
        frc(3*(i-1)+3) = - frc_quick(3,i) * AU_TO_EV * EV_TO_KCAL / BOHRS_TO_A
      end do
      pot_ene = pot_ene * AU_TO_EV * EV_TO_KCAL
#endif
    else if (upcase(extprog) .eq. 'KMMD') then
      do i = 1, natom
        coord(3*(i-1)+1) = crd(3*(i-1)+1)
        coord(3*(i-1)+2) = crd(3*(i-1)+2)
        coord(3*(i-1)+3) = crd(3*(i-1)+3)
        grads(3*(i-1)+1) = 0.0
        grads(3*(i-1)+2) = 0.0
        grads(3*(i-1)+3) = 0.0
      end do

      call kmmd_frccalc(coord, grads, pot_ene)

      do i = 1, natom
        frc(3*(i-1)+1) = frc(3*(i-1)+1) + grads(3*(i-1)+1)
        frc(3*(i-1)+2) = frc(3*(i-1)+2) + grads(3*(i-1)+2)
        frc(3*(i-1)+3) = frc(3*(i-1)+3) + grads(3*(i-1)+3)
      end do
    end if

  end subroutine

  subroutine pme_external(crd, frc, pot_ene)

    use memory_module, only : natom
    use nblist, only: ucell

    IMPLICIT NONE

    double precision           ::  crd(*)
    double precision           ::  frc(*)
    double precision           ::  pot_ene
    integer                    ::  i, i3
    double precision           ::  coord(3*natom), grads(3*natom), box(9)

    ! For MBX
    if (Upcase(extprog) .eq. 'MBX') then
#ifdef MBX
      do i = 1, natom
        i3=3*(i-1)
        coord(3*(i-1)+1) = crd(i3+1)
        coord(3*(i-1)+2) = crd(i3+2)
        coord(3*(i-1)+3) = crd(i3+3)
      end do

      box(1) = ucell(1,1)
      box(2) = ucell(2,1)
      box(3) = ucell(3,1)
      box(4) = ucell(1,2)
      box(5) = ucell(2,2)
      box(6) = ucell(3,2)
      box(7) = ucell(1,3)
      box(8) = ucell(2,3)
      box(9) = ucell(3,3)

      call get_energy_pbc_g(coord, natom, box, pot_ene, grads)

      do i = 1, natom
        i3=3*(i-1)
        frc(i3+1) = frc(i3+1) - grads(3*(i-1)+1)
        frc(i3+2) = frc(i3+2) - grads(3*(i-1)+2)
        frc(i3+3) = frc(i3+3) - grads(3*(i-1)+3)
      end do
#endif
    end if
    
    if (upcase(extprog) .eq. 'KMMD') then
      do i = 1, natom
        i3=3*(i-1)
        coord(3*(i-1)+1) = crd(i3+i)
        coord(3*(i-1)+2) = crd(i3+i+1)
        coord(3*(i-1)+3) = crd(i3+i+2)
        grads(3*(i-1)+1) = 0.0
        grads(3*(i-1)+2) = 0.0
        grads(3*(i-1)+3) = 0.0
      end do

      call kmmd_frccalc(coord, grads, pot_ene)

      do i = 1, natom
        i3=3*(i-1)
        frc(i3+1) = frc(i3+1) + grads(3*(i-1)+1)
        frc(i3+2) = frc(i3+2) + grads(3*(i-1)+2)
        frc(i3+3) = frc(i3+3) + grads(3*(i-1)+3)
      end do
    end if

  end subroutine

  subroutine external_cleanup()

    IMPLICIT NONE

    integer :: ierr

    ! For Quick
    if (Upcase(extprog) .eq. 'QUICK') then
#ifdef QUICK
       call deleteQuickJob(ierr)
       if ( ierr /= 0 ) then
         write(6, '(a)') 'ERROR: ending Quick in deleteQuickJob at external.F90'
         call mexit(6, 1)
       end if
#endif
    end if

  end subroutine

end module external_module
