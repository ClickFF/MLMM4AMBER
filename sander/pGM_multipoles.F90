#include "../include/dprec.fh"
#include "../include/assert.fh"

!---------------------------------------------------------
module pol_gauss_multipoles

  implicit none
  private

  integer,save :: do_flag
  _REAL_,parameter :: coulomb_const_kcal_per_mole = 332.05382d0 !Tinker

  integer, save :: num_multipoles
  integer, save :: start_multipoles,end_multipoles
  _REAL_, save, allocatable :: local_multipole(:,:)
  _REAL_, save, allocatable :: global_multipole(:,:)
  integer,parameter :: MAXMP=35

  !Gaussian radii
  _REAL_, parameter :: betamax = sqrt(2.0d0)/0.7147597d0
  _REAL_, dimension(:), allocatable  :: r_gauss, b_polar

  !Covalent local frame lists
  integer, dimension(:), allocatable  :: covalent_ptr, covalent_atm
  _REAL_, dimension(:), allocatable  :: covalent_monop, covalent_dip

  public global_multipole,r_gauss,betamax,covalent_ptr,covalent_atm,covalent_dip, &
       pGM_MPOLE_readparm,&
       pGM_MPOLE_deallocate,pGM_MPOLE_global_to_local,pGM_MPOLE_local_to_global, &
       coulomb_const_kcal_per_mole, &
       pGM_MPOLE_bcast

  contains
!---------------------------------------------------------
subroutine pGM_MPOLE_get_start_end_lists()

  integer :: siztask
  call POL_GAUSS_get_startlist_endlist(num_multipoles,  &
                             start_multipoles,end_multipoles,siztask)
! find local mpoles boundaries
!write(6,*)'start_multipoles,end_multipoles = ', start_multipoles,end_multipoles
end subroutine pGM_MPOLE_get_start_end_lists
!---------------------------------------------------------
subroutine pGM_MPOLE_bcast()
  implicit none
#ifdef MPI
  integer ier

  include 'mpif.h'
# include "extra.h"
# include "parallel.h"
  call mpi_bcast(do_flag,1,MPI_INTEGER,0,commsander,ier)
  call mpi_bcast(num_multipoles,1,MPI_INTEGER,0,commsander,ier)

  if(.not.master) then
    allocate(global_multipole(10,num_multipoles),stat=ier)
    REQUIRE(ier==0)
  end if
     
  call mpi_bcast(local_multipole,10*num_multipoles,MPI_DOUBLE_PRECISION,0,commsander,ier)
 
  call pGM_MPOLE_get_start_end_lists()
#endif 
end subroutine pGM_MPOLE_bcast
!---------------------------------------------------------
function pGM_MPOLE_readparm(nf,num_atoms)
  use pol_gauss_mdin, only : pol_gauss_verbose
  integer :: pGM_MPOLE_readparm
  integer, intent(in) :: nf,num_atoms
# include "do_flag.h"  

  integer :: n,ier,dim1,ncovalent

  pGM_MPOLE_readparm = 0

  num_multipoles = num_atoms

  ! first. read in covalent atom/local frame pointer
  allocate(covalent_ptr(0:num_atoms), stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call POL_GAUSS_read_list_data('POL_GAUSS_COVALENT_POINTERS_',nf, &
                 dim1,num_atoms,covalent_ptr(1))
  ncovalent = sum (covalent_ptr(1:num_atoms))
  if (pol_gauss_verbose == 2 ) then
    write(116,'(///)')
    write(116,*) 'Finding covalent pointers', num_atoms
    write(116,'(10i8)') covalent_ptr(1:num_atoms)
    write(116,'(///)')
  end if
  covalent_ptr(0) = 0
  do n = 1, num_atoms
    covalent_ptr(n) = covalent_ptr(n) + covalent_ptr(n-1)
  end do

  allocate(covalent_atm(1:ncovalent),stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call POL_GAUSS_read_list_data('POL_GAUSS_COVALENT_ATOMS_',nf, &
                 dim1,ncovalent,covalent_atm)
  if (pol_gauss_verbose == 2 ) then
    write(116,'(///)')
    write(116,*) 'Finding covalent atoms', ncovalent
    write(116,'(10i8)') covalent_atm(1:ncovalent)
    write(116,'(///)')
  end if

  ! third. read in throughbond moments (in electron Angstrom)
  allocate(covalent_dip(1:ncovalent),stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call POL_GAUSS_read_real_list_data('POL_GAUSS_COVALENT_DIPOLES_',nf, &
                 dim1,ncovalent,covalent_dip)
  if (pol_gauss_verbose == 2 ) then
    write(116, '(///)')
    write(116,*) 'Finding throughbond moments',ncovalent 
    write(116,'(5E16.8)') covalent_dip(1:ncovalent)
    write(116, '(///)')
  end if

  ! fourth. read in gaussian monopole
  allocate(covalent_monop(num_atoms),stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call POL_GAUSS_read_real_list_data('POL_GAUSS_MONOPOLES_',nf, &
                 dim1,num_atoms,covalent_monop)
  if (pol_gauss_verbose == 2 ) then
    write(116,'(///)')
    write(116,*) 'Finding gaussian monopoles', num_atoms
    write(116,'(5E16.8)') covalent_monop(1:num_atoms)
    write(116,'(///)')
  end if

  ! allocate global multipoles
  allocate(global_multipole(10,num_atoms),stat=ier)
  REQUIRE(ier==0)

  ! fifth. read in gaussian multipole radius
  allocate(b_polar(num_atoms),stat=ier)
  REQUIRE(ier==0)
  allocate(r_gauss(num_atoms),stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call POL_GAUSS_read_real_list_data('POL_GAUSS_RADII_',nf, &
                 dim1,num_atoms,b_polar)
  if (pol_gauss_verbose == 2 ) then
    write(116,'(///)')
    write(116,*) 'Finding gaussian radii', num_atoms
    write(116,'(5E16.8)') b_polar(1:num_atoms)
    write(116,'(///)')
  end if
  ! b_polar is read in as radii
  ! preprocess radius and beta for efficient combination
  r_gauss(1:num_atoms) = 2.0d0*b_polar(1:num_atoms)**2 ! this is to be consistent with pyresp
  !r_gauss(1:num_atoms) = 1.0d0*b_polar(1:num_atoms)**2 ! this is the original condition used for MD testing.
  !r_polar(1:num_atoms) = 1.0d0/b_polar(1:num_atoms)

  call pGM_MPOLE_get_start_end_lists()

  pGM_MPOLE_readparm = 1
  do_flag = ibset(do_flag,VALID_BIT)
end function pGM_MPOLE_readparm
!---------------------------------------------------------
subroutine pGM_MPOLE_rescale_multipoles()

  integer j,n
  _REAL_ bohr

  ! change units from Bohr to Angstroms
  !bohr = 0.5291772083d0 ! R. Luo: if input uses Bohr as the length unit
  bohr = 1.0d0/4.803239377d0 ! R. Luo: if input dipoles are read in as debyes
  !bohr = 1.0d0 ! R. Luo: if input dipoles are read in as electron*Angstrom

  !the following is changed to scale input global moments only
  do n = 1,num_multipoles
    do j = 2,4
      global_multipole(j,n) = bohr*global_multipole(j,n)
    end do
  end do
end subroutine pGM_MPOLE_rescale_multipoles
!---------------------------------------------------------
subroutine pGM_MPOLE_deallocate()
  if ( allocated(local_multipole) )deallocate(local_multipole)
end subroutine pGM_MPOLE_deallocate
!---------------------------------------------------------
subroutine pGM_MPOLE_global_to_local(crd)
  use pol_gauss_mdin, only : pol_gauss_verbose
  use constants, only : ZERO, ONE
  use stack
  _REAL_,intent(in) :: crd(3,*)
# include "do_flag.h"

  integer :: iatm, jbond, jatm, start_bonds, end_bonds
  integer, parameter :: mp = 3, np = 4
  integer i, j, info, ier
  integer m, n
  _REAL_, parameter :: TOL = 1d-6
  _REAL_ a(mp,np), u(mp,np), w(np), v(np,np), vt(np,np), b(mp), t(np)
  _REAL_ wmax, thresh
  _REAL_ work(512)
  _REAL_, allocatable :: vector(:,:), norm(:)

  ! if multipoles are read from prmtop file in the global frame
  ! the input multipoles can be converted into covalent dipoles here
  start_bonds = 1
  do iatm = start_multipoles, end_multipoles
    start_bonds = start_bonds + covalent_ptr(iatm-1)
    end_bonds = start_bonds + covalent_ptr(iatm) - 1

    ! compute through-bond basis vectors
    allocate(vector(1:3,start_bonds:end_bonds), stat=ier)
    REQUIRE(ier==0)
    allocate(norm(start_bonds:end_bonds), stat=ier)
    REQUIRE(ier==0)
    do jbond = start_bonds, end_bonds
      jatm = covalent_atm(jbond)
      vector(1:3,jbond) = crd(1:3,jatm) - crd(1:3,iatm)
      norm(jbond) = ONE/sqrt(sum(vector(1:3,jbond)**2))
    end do
    do jbond = start_bonds, end_bonds
      vector(1:3,jbond) = vector(1:3,jbond)*norm(jbond)
    end do

    ! set up linear system for least square fitting.
    ! Here the goal is to obtain the coefficients T for X, the basis vectors
    ! so that X * T = b. Thus linear system matrix A is X as written when
    ! solving for T.
    a = ZERO; u = ZERO; w = ZERO; vt = ZERO; b = ZERO; t = ZERO
    m = 3
    n = covalent_ptr(iatm)
    if (n > np) then
      write(6,*) ' pGM_MPOLE_global_to_local(): Too many bonds for SVD'
      call mexit(6,1)
    end if
    do jbond = start_bonds, end_bonds
      a(1:3, jbond-start_bonds+1) = vector(1:3, jbond)
    end do
    b(1:3) = global_multipole(2:4, iatm)
    deallocate(vector, stat=ier)
    REQUIRE(ier==0)
    deallocate(norm, stat=ier)
    REQUIRE(ier==0)

    ! svd call for U, W, V as in A = U * 1/W * V^T 
    ! DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call dgesvd('A','A',m,n,a,mp,w,u,mp,vt,np,work,512,info)
    do i = 1, np
      do j = 1, np
        v(i,j) = vt(j,i)
      end do
    end do
    wmax = ZERO
    do j = 1, n
      if ( w(j) > wmax ) wmax = w(j)
    end do
    thresh = TOL * wmax
    do j = 1, n
      if ( w(j) < thresh ) w(j) = ZERO
    end do

    ! back substitution call for A^-1 * b
    call svbksb(u,w,v,m,n,mp,np,b,t)

    ! store through-bond dipole moment in T for later
    if ( pol_gauss_verbose == 2 ) write(6,*) ' Atom info', iatm, info
    do jbond = start_bonds, end_bonds
      covalent_dip(jbond) = t(jbond-start_bonds+1)
      if ( pol_gauss_verbose == 2 ) write(6,*) &
        ' through-bond moments', covalent_dip(jbond)
    end do
  end do ! end of iatm = start_multipoles, end_multipoles

  if ( pol_gauss_verbose == 2 ) then
    write(6,*) 'global multipoles usd to compute through-bond moments (D)'
    write(6,'(f15.6)') global_multipole(1:1,1:num_multipoles)
    write(6,'(3f15.6)') global_multipole(2:4,1:num_multipoles)*4.803239377d0
  end if
end subroutine pGM_MPOLE_global_to_local
!---------------------------------------------------------
subroutine pGM_MPOLE_local_to_global(crd)
  use pol_gauss_mdin, only : pol_gauss_verbose
  use constants, only : ZERO, ONE
  use stack
  _REAL_,intent(in) :: crd(3,*)
#include "do_flag.h"

  integer ier
  integer :: iatm, jbond, jatm, start_bonds, end_bonds
  _REAL_, allocatable :: vector(:,:), norm(:)

  if ( .not. btest(do_flag,VALID_BIT) )return

  global_multipole(1,start_multipoles:end_multipoles) = &
  covalent_monop(start_multipoles:end_multipoles)
  global_multipole(2:4,start_multipoles:end_multipoles) = 0.0d0

  do iatm = start_multipoles, end_multipoles
    start_bonds = covalent_ptr(iatm-1) + 1
    end_bonds = covalent_ptr(iatm)

    ! compute through-bond unit vectors
    allocate(vector(1:3,start_bonds:end_bonds), stat=ier)
    REQUIRE(ier==0)
    allocate(norm(start_bonds:end_bonds), stat=ier)
    REQUIRE(ier==0)
    do jbond = start_bonds, end_bonds
      jatm = covalent_atm(jbond)
      vector(1:3,jbond) = crd(1:3,jatm) - crd(1:3,iatm)
      norm(jbond) = ONE/sqrt(sum(vector(1:3,jbond)**2))
    enddo
    do jbond = start_bonds, end_bonds
      vector(1:3,jbond) = vector(1:3,jbond)*norm(jbond)
    enddo

    ! recover global dipoles
    do jbond = start_bonds, end_bonds
      global_multipole(2:4,iatm) = global_multipole(2:4,iatm) + covalent_dip(jbond)*vector(1:3,jbond)
    enddo
    deallocate(vector, stat=ier)
    REQUIRE(ier==0)
    deallocate(norm, stat=ier)
    REQUIRE(ier==0)
  enddo
  
  if ( pol_gauss_verbose == 2 ) then
    write(6,*) 'converted global multipoles from through-bond dipoles (D)'
    write(6,'(f15.6)') global_multipole(1:1,1:num_multipoles)
    write(6,'(3f15.6)') global_multipole(2:4,1:num_multipoles)*4.803239377d0
  end if
end subroutine pGM_MPOLE_local_to_global
!---------------------------------------------------------
end module pol_gauss_multipoles
