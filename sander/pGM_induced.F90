#include "../include/dprec.fh"
#include "../include/assert.fh"

!------------------------------------------------------------
module pol_gauss_induced
  implicit none

  private

  integer, save :: do_flag
  logical, save :: initial

  logical,save, allocatable :: is_polarizable(:)
  _REAL_, save, allocatable :: polarizability(:)
  _REAL_, save, allocatable :: ind_dip(:,:)
  _REAL_, save, allocatable :: dip_field(:,:)
  _REAL_, save, allocatable :: cart_dipole_field(:,:)
  _REAL_, save, allocatable :: record1(:,:)
  _REAL_, save, allocatable :: record2(:,:)
  _REAL_, save, allocatable :: record3(:,:)
  _REAL_, save, allocatable :: record4(:,:)
  _REAL_, save, allocatable :: oldrec(:,:)
  _REAL_, save, allocatable :: ap(:)
  _REAL_, save, allocatable :: ip(:)
  _REAL_, save, allocatable :: xv(:)
  _REAL_, save, allocatable :: bv(:)
  _REAL_, save, allocatable :: zv(:)
  _REAL_, save, allocatable :: zu(:)
  _REAL_, save, allocatable :: rv(:)
  _REAL_, save, allocatable :: ru(:)
  _REAL_, save, allocatable :: pv(:)
  _REAL_, save, allocatable :: pu(:)

  private :: scf_solv_init, pGM_cg, pGM_pcg, pGM_sor, pGM_pkcg, local_precondition

  public pGM_INDUCED_eval,pGM_INDUCED_set_user_bit,pGM_INDUCED_readparm, &
         polarizability,is_polarizable,ind_dip
# ifdef MPI
  public pGM_INDUCED_bcast
# endif

contains
!------------------------------------------------------------------------
subroutine pGM_INDUCED_set_user_bit(do_this)
  integer,intent(in) :: do_this
#include "do_flag.h"

  if ( do_this == 1 )then
    do_flag = ibset(do_flag,USER_BIT)
    call pGM_INDUCED_set_induced_flags()
  else
    do_flag = ibclr(do_flag,USER_BIT)
  endif
end subroutine pGM_INDUCED_set_user_bit
!------------------------------------------------------------------------
#ifdef MPI
subroutine pGM_INDUCED_bcast(num_atoms)
   implicit none
   integer num_atoms,ier,n,j

   include 'mpif.h'
#  include "extra.h"
#  include "parallel.h"
   call mpi_bcast(do_flag,1,MPI_INTEGER,0,commsander,ier)

   if(.not.master) then
      allocate(polarizability(num_atoms),stat=ier)
      REQUIRE(ier==0)
      allocate(is_polarizable(num_atoms),stat=ier)
      REQUIRE(ier==0)
      allocate(ind_dip(3,num_atoms),stat=ier)
      REQUIRE(ier==0)
      allocate(dip_field(3,num_atoms),stat=ier)
      REQUIRE(ier==0)
      allocate(cart_dipole_field(3,num_atoms),stat=ier)
      REQUIRE(ier==0)
   end if  

   if(.not.master) then
      do n = 1,num_atoms
         is_polarizable(n) = .false.
         polarizability(n) = 0.d0
         do j = 1,3
            ind_dip(j,n) = 0.d0
         end do
      end do
   end if

   call mpi_bcast(polarizability,num_atoms,MPI_DOUBLE_PRECISION,0,commsander,ier)
  
end subroutine pGM_INDUCED_bcast
#endif
!------------------------------------------------------------------------
function pGM_INDUCED_readparm(nf,num_atoms)
  use pol_gauss_mdin, only : dipole_scf_init_step,pol_gauss_verbose
  implicit none

  integer :: pGM_INDUCED_readparm
  integer,intent(in) :: nf,num_atoms

#include "do_flag.h"  
  integer :: j,n,ier,dim1,numpolar
  _REAL_ :: polmin

  pGM_INDUCED_readparm = 0
  polmin = 1.d-6

  ! initialize things
  allocate(polarizability(num_atoms),stat=ier)
  REQUIRE(ier==0)
  allocate(is_polarizable(num_atoms),stat=ier)
  REQUIRE(ier==0)
  allocate(ind_dip(3,num_atoms),stat=ier)
  REQUIRE(ier==0)
  allocate(dip_field(3,num_atoms),stat=ier)
  REQUIRE(ier==0)
  allocate(cart_dipole_field(3,num_atoms),stat=ier)
  REQUIRE(ier==0)
  allocate(record1(3*num_atoms,dipole_scf_init_step+1),stat=ier)
  REQUIRE(ier==0)
  allocate(record2(3*num_atoms,dipole_scf_init_step+1),stat=ier)
  REQUIRE(ier==0)
  allocate(record3(3*num_atoms,dipole_scf_init_step+1),stat=ier)
  REQUIRE(ier==0)
  allocate(record4(3*num_atoms,dipole_scf_init_step+1),stat=ier)
  REQUIRE(ier==0)
  allocate(oldrec(3*num_atoms,2),stat=ier)
  REQUIRE(ier==0)

  is_polarizable = .false.
  polarizability = 0.d0
  ind_dip = 0.0d0
  record1 = 0.0d0
  record2 = 0.0d0
  record3 = 0.0d0
  record4 = 0.0d0

  dim1 = 1
  numpolar = num_atoms
  call POL_GAUSS_read_real_list_data('POL_GAUSS_POLARIZABILITY_',nf,dim1,numpolar,polarizability)
  if (pol_gauss_verbose == 2) then
    write(116, '(///)')
    write(116,*) 'Finding polarizability',numpolar
    write(116,'(5E16.8)') polarizability(1:numpolar)
    write(116, '(///)')
  end if
 
  do n = 1,numpolar
    if ( polarizability(n) > polmin )then
      is_polarizable(n) = .true.
    endif
  enddo
  pGM_INDUCED_readparm = 1
  initial = .true.

  call scf_solv_init(num_atoms)

  do_flag = ibset(do_flag,VALID_BIT)
end function pGM_INDUCED_readparm
!------------------------------------------------------------------------
subroutine pGM_INDUCED_set_induced_flags()
#  include "ew_cntrl.h"
#  include "ew_mpole.h"
  induced = 1
  indmeth = 1
end subroutine pGM_INDUCED_set_induced_flags
!------------------------------------------------------------------------
subroutine pGM_INDUCED_eval(numatoms,crd,x,ipairs,diprms,dipiter)
  use pol_gauss_mdin, only : dipole_scf_init,dipole_scf_init_order,use_average,&
                             dipole_scf_init_step,scf_solv_opt,pol_gauss_verbose
  use constants, only : TEN_TO_MINUS6, ZERO 

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: crd(3,*)
  _REAL_,intent(in) :: x(*)
  integer, intent(in) :: ipairs(*)
  _REAL_,intent(out) :: diprms, dipiter

  integer :: n

#include "do_flag.h"  

  if ( do_flag /= PROCEED ) return

  if (pol_gauss_verbose == 2) then
    if (use_average == 1) then
      write (6, '(a)') ' pGM_Induced_eval: Use average error to check convergence' 
    else
      write (6, '(a)') ' pGM_Induced_eval: Use maximum error to check convergence' 
    end if
    write (6, '(a)') ' pGM_Induced_eval: Setting up rhs of the linear system' 
  end if
  ! set up the linear system AM * xv = bv, in short as A * x = b
  ! the off-diagonal elements of A (both full and short) are the
  ! dip-dip tensor elements set up in the call, the diagoanl elements
  ! of A is 1/pol, the right hand side bv() is E_perm, which is
  ! obtained after adding the negative sign to the returned gradient zv
  call pGM_NonBond_perm_fields(numatoms, is_polarizable, crd, x, ipairs, &
           zv, cart_dipole_field)
  bv = -zv ! this is E_perm

  ! return zero induced dipoles if there is no E_perm field
  if ( sum(abs(bv)) .lt. TEN_TO_MINUS6 ) then
     ind_dip = ZERO 
     return
  end if

  ! initialize induced dipoles
  select case ( dipole_scf_init )
    case (1)
      if (pol_gauss_verbose == 2) write(6,'(/,a)') ' pGM Induced_eval: initial dipole set as pol*E_perm'
      xv = ap*bv
    case (2)
      if (pol_gauss_verbose == 2) write(6,'(/,a)') ' pGM Induced_eval: initial dipole set from last step'
      if ( initial ) then
        xv = ap*bv
        initial = .false.
      end if
    case (3)
      if (pol_gauss_verbose == 2) write(6,'(/,a)') ' pGM Induced_eval: initial dipole extropolated from multiple steps'
      if ( dipole_scf_init_order > 3 ) then
        write(6,'(/,a)') ' pGM_Induced_eval: dipole_scf_init_order too high'
        call mexit(6,1)
      end if
      if ( dipole_scf_init_step > 8 ) then
        write(6,'(/,a)') ' pGM_Induced_eval: dipole_scf_init_step too high'
        call mexit(6,1)
      end if
      call interp_ind_dipole(numatoms,dipole_scf_init_order,dipole_scf_init_step,ap,bv,xv)
    case default
      ! invalid dipole_scf_init
      ASSERT( .false. )
  end select

  ! solve the linear system
  select case ( scf_solv_opt )
    case ( 1 )
      if (pol_gauss_verbose /= 0) write(6,'(/,a)') ' pGM CG selected for induced dipole solution '
      call pGM_cg(numatoms, x, diprms, dipiter, bv, xv)
    case ( 2 )
      if (pol_gauss_verbose /= 0) write(6,'(/,a)') ' pGM PCG selected for induced dipole solution '
      call pGM_pcg(numatoms, x, diprms, dipiter, bv, xv)
    case ( 3 )
      if (pol_gauss_verbose /= 0) write(6,'(/,a)') ' pGM PCG with Peek selected for induced dipole solution '
      call pGM_pkcg(numatoms, x, diprms, dipiter, bv, xv)
    case ( 4 )
      if (pol_gauss_verbose /= 0) write(6,'(/,a)') ' pGM SOR selected for induced dipole solution '
      call pGM_sor(numatoms, x, diprms, dipiter, bv, xv)
    case ( 5 )
      if (pol_gauss_verbose /= 0) write(6,'(/,a)') ' pGM Ehanced SOR selected for induced dipole solution '
      call pGM_esor(numatoms, x, diprms, dipiter, bv, xv)
    case default
      ! invalid scf_solv_opt
      ASSERT( .false. )
  end select

  ! save some data for initial guess of induced dipoles
  if ( dipole_scf_init == 3 ) then
    call update_record(numatoms,dipole_scf_init_step,record2)
    record2(:,1) = xv(:)
    if ( dipole_scf_init_order >= 2 .and. &
         sum(abs(record2(:,dipole_scf_init_step+1))) > 0d0 ) then
      call update_record(numatoms,dipole_scf_init_step,record3)
      record3(:,1) = xv(:) - oldrec(:,1)
    end if
    if ( dipole_scf_init_order == 3 .and. &
         sum(abs(record3(:,dipole_scf_init_step+1))) > 0d0 ) then
      call update_record(numatoms,dipole_scf_init_step,record4)
      record4(:,1) = xv(:) - oldrec(:,2)
    end if
  end if

  ! return final induced dipoles for energy/force processing
  do n = 1, numatoms
    ind_dip(1,n) = xv(1+3*(n-1))
    ind_dip(2,n) = xv(2+3*(n-1))
    ind_dip(3,n) = xv(3+3*(n-1))
  end do

end subroutine pGM_INDUCED_eval
!----------------------------------------------------------
subroutine scf_solv_init(numatoms)
  use pol_gauss_mdin, only : pol_gauss_verbose,scf_cg_niter
  implicit none

  integer, intent(in) :: numatoms

#if defined(MPI) && defined(POL_GAUSS_WORKS_WITH_IT)
  include 'mpif.h'
# include "extra.h"
# include "parallel.h"
#endif /* MPI */

  integer :: ierr, n

#ifdef MPI
!   call mpi_bcast(1,MPI_INTEGER,0,commsander,ierr)
#endif /* MPI */

  if (pol_gauss_verbose /= 0) write (6, '(/a)') ' pGM_INDUCED_readparm(): Allocating working arrays'

  allocate(ap(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)
  allocate(ip(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)
  allocate(xv(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)
  allocate(bv(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)
  allocate(rv(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)
  allocate(ru(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)
  allocate(pv(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)
  allocate(pu(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)
  allocate(zv(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)
  allocate(zu(3*numatoms), stat = ierr)
  REQUIRE(ierr==0)

  scf_cg_niter = max(2, scf_cg_niter)

  ! set up polarizability arrays for vectorized process
  do n = 1, numatoms
    ap(1+3*(n-1)) = polarizability(n)
    ap(2+3*(n-1)) = polarizability(n)
    ap(3+3*(n-1)) = polarizability(n)
  end do
  do n = 1, numatoms
    ip(1+3*(n-1)) = 1.0d0/polarizability(n)
    ip(2+3*(n-1)) = 1.0d0/polarizability(n)
    ip(3+3*(n-1)) = 1.0d0/polarizability(n)
  end do

end subroutine scf_solv_init
!----------------------------------------------------------
subroutine interp_ind_dipole(numatoms,dipole_scf_init_order,&
                             dipole_scf_init_step,pol,e_perm,ind_dip)

  integer :: numatoms, dipole_scf_init_order, dipole_scf_init_step
  _REAL_ :: pol(1:3*numatoms)
  _REAL_ :: e_perm(1:3*numatoms)
  _REAL_ :: ind_dip(1:3*numatoms)

  integer :: n, i, j
  _REAL_ :: c(8)

  ! record1 is used to save pol*E_perm from previous steps to
  ! Prepare for interp by obtaining interp coefficients c(:)
  call update_record(numatoms,dipole_scf_init_step,record1)
  record1(:,1) = pol*e_perm
  if (sum(abs(record1(:,dipole_scf_init_step+1))) > 0d0) then
    call least_square_interp(numatoms,dipole_scf_init_step,record1,c)
  end if

  ! record2 is used to save total ind_dip()
  ! First interp the total ind_dip from previous steps
  if ( dipole_scf_init_order >= 1 ) then
    if (sum(abs(record2(:,dipole_scf_init_step))) > 0d0) then
    ! we have enough data, so proceed with interpolation
      ind_dip = 0d0
      do j = 1, dipole_scf_init_step
        ind_dip = ind_dip + c(j)*record2(:,j)
      end do
      oldrec(:,1) = ind_dip ! the predicted ind_dip is saved in oldrec(,1)
    else
    ! don't have enough data, just use pol*E_perm
      if ( dipole_scf_init_order == 1 ) ind_dip = record1(:,1)
    end if
  end if

  ! record3 is used to save the difference of ind_dip(), or record2()
  ! Second interp the difference from previous steps and
  ! add to total ind_dip
  if ( dipole_scf_init_order >= 2 ) then
    if (sum(abs(record3(:,dipole_scf_init_step))) > 0d0) then
    ! we have enough data, so proceed with interpolation
      do j = 1, dipole_scf_init_step
        ind_dip = ind_dip + c(j)*record3(:,j)
      end do
      ! the predicted ind_dip is saved in oldrec(,2)
      oldrec(:,2) = ind_dip
    else
    ! don't have enough data, just use pol*E_perm
      if ( dipole_scf_init_order == 2 ) ind_dip = record1(:,1)
    end if
  end if

  ! record4 is used to save the difference of the difference of ind_dip
  ! Third interp the difference of the difference from previous steps and
  ! add to total ind_dip
  if ( dipole_scf_init_order == 3 ) then
    if (sum(abs(record4(:,dipole_scf_init_step))) > 0d0) then
    ! we have enough data, so proceed with interpolation
      do j = 1, dipole_scf_init_step
        ind_dip = ind_dip + c(j)*record4(:,j)
      end do
    else
    ! don't have enough data, just use pol*E_perm
      ind_dip = record1(:,1)
    end if
  end if

  contains
!----------------------------------------------------------
subroutine least_square_interp(numatoms,dipole_scf_init_step,record,c)

  integer :: numatoms, dipole_scf_init_step
  _REAL_ :: record(3*numatoms,dipole_scf_init_step+1),c(8)

  _REAL_ :: matrix(dipole_scf_init_step,dipole_scf_init_step),b(dipole_scf_init_step)
  integer :: i,j,k

  matrix = 0d0
  b = 0d0
  c = 0d0
  do i = 1, dipole_scf_init_step
    do j = i, dipole_scf_init_step
      do k = 1, 3*numatoms
        matrix(i,j) = matrix(i,j) + record(k,i+1)*record(k,j+1)
        if (j == dipole_scf_init_step) b(i) = b(i) + record(k,i+1)*record(k,1)
      end do
      if (i < j) matrix(j,i) = matrix(i,j)
    end do
  end do

  call solve_least_square(dipole_scf_init_step,matrix,b,c)

end subroutine least_square_interp
!----------------------------------------------------------
subroutine solve_least_square(dipole_scf_init_step,matrix,b,c)

  integer :: dipole_scf_init_step
  _REAL_ :: matrix(dipole_scf_init_step,dipole_scf_init_step),b(dipole_scf_init_step),c(8)

  integer :: i,j,k
  _REAL_ :: cg_r(dipole_scf_init_step), cg_p(dipole_scf_init_step), cg_matrix_p(dipole_scf_init_step)
  _REAL_ :: alpha, beta

  c = 0d0
  cg_r = b
  cg_p = cg_r

  do i = 1, dipole_scf_init_step
    cg_matrix_p = 0d0
    do j = 1, dipole_scf_init_step
      do k = 1, dipole_scf_init_step
        cg_matrix_p(j) = cg_matrix_p(j) + matrix(j,k)*cg_p(k)
      end do
    end do
    alpha = dot_product(cg_r,cg_r)/dot_product(cg_matrix_p,cg_p)
    c = c + alpha*cg_p
    beta = dot_product(cg_r,cg_r)
    cg_r = cg_r - alpha*cg_matrix_p
    beta = dot_product(cg_r,cg_r)/beta
    cg_p = cg_r + beta*cg_p
  end do

end subroutine solve_least_square
!----------------------------------------------------------
end subroutine interp_ind_dipole
!----------------------------------------------------------
subroutine update_record(numatoms,dipole_scf_init_step,record)
  implicit none

  integer :: numatoms,dipole_scf_init_step
  _REAL_ record(3*numatoms,dipole_scf_init_step+1)
  integer :: i

  do i = dipole_scf_init_step, 1, -1
    record(:,i+1) = record(:,i)
  end do

end subroutine update_record
!----------------------------------------------------------
! To solve A*x = b with the SOR solver, we intend to utilize the
! iterative approach: x(i+1) = x(i) + wsor*(b - A*x(i))/D
! When wsor = 1, this is just the Gauss-Seidel iteration.
! The code is documented for educational purposes only and
! for comparison with other faster solvers.
!
! When A can be easily written as L+D+U, there are additional
! tricks possible, but we can only reasonably use A = T + D,
! D=1/pol, so we have x(i+1) = x(i) + wsor*r(i))*pol,
! with r(i) = b - A*x(i)
! 
! In PB we are using the following, which is to use the latest
! x's to update the subsequent x's as
! x(i+1) = x(i) + wsor*[b - (T + 1/pol)*x(i)]*pol
! x(i+1) = x(i) + wsor*[b - T*x(i)]*pol - wsor*1/pol*x(i)*pol
! x(i+1) = x(i) + wsor*[b - T*x(i)]*pol - wsor*x(i)
! x(i+1) = [1-wsor]*x(i) + wsor*[b - T*x(i)]*pol
! In the above loop, the latest x's are used to speed up convergence.
! since T*x must be realized in one call in pGM, we can't use
! the same trick.
!----------------------------------------------------------
subroutine pGM_sor(numatoms, x, r_norm, iter, bv, xv)
  use pol_gauss_mdin, only : scf_sor_niter,scf_sor_coefficient,&
                             dipole_scf_tol,pol_gauss_verbose
  implicit none

  integer, intent(in) :: numatoms
  _REAL_,  intent(in) :: x(*)
  _REAL_,  intent(out) :: r_norm, iter
  _REAL_,  intent(in) :: bv(1:3*numatoms)
  _REAL_,  intent(inout) :: xv(1:3*numatoms)

  logical :: unconvg

  unconvg = .true.
  iter = 0.0d0

  ! main loop
  do while ( unconvg )

    ! Note that A's off-diagonal elements are those of the dipole-dipole tensor and
    ! its diagonal elements are 1/pol (1/polarizability)
    call pGM_NonBond_dip_dip_fields(numatoms, x, xv, zv)
    zv = zv + ip*xv
    rv = bv - zv

    iter = iter + 1.0d0

    ! check convergence
    r_norm = maxval(abs(rv/bv))
    if ( pol_gauss_verbose == 2 ) write (6, '(a,i3,a,f14.10)') &
            ' pGM SOR: iter = ', int(iter), ', residue norm: ', r_norm
    if ( int(iter) >= scf_sor_niter .or. r_norm <= dipole_scf_tol ) then
      unconvg = .false.
      if ( int(iter) >= scf_sor_niter ) then
        write(6,'(a)') ' pGM SOR: FATAL maxitn exceeded!'
        call mexit(6,1)
      end if
    end if
    xv = xv + scf_sor_coefficient * rv * ap

  end do ! do while (unconvg)

  if ( pol_gauss_verbose /= 0 )then
    write(6,'(/a,i5,1x,f14.9/)') ' pGM SOR: out of loop, iter, norm = ',int(iter),r_norm
  endif

end subroutine pGM_sor
!----------------------------------------------------------
! To solve A*x = b with the ESOR solver, we intend to utilize the
! relation of A = T + D, D=1/pol, and T = T_short + T_long
! so we have T_short*x + T_long*x +D*x = b
! (T_short+D)*x = b - T_long*x
! This shows that we can use the local_preconditoner matrix
! to do the iteration given a revised righ hand side.
!
! Note in this method, A*x can still be computed as T*x has been computed.
! So r = b - A*x is ready available. This allows us to use it to
! check convergence just like other solvers. 
!----------------------------------------------------------
subroutine pGM_esor(numatoms, x, r_norm, iter, bv, xv)
  use pol_gauss_mdin, only : scf_sor_niter,scf_sor_coefficient,&
                             dipole_scf_tol,scf_local_niter,pol_gauss_verbose
  implicit none

  integer, intent(in) :: numatoms
  _REAL_,  intent(in) :: x(*)
  _REAL_,  intent(out) :: r_norm, iter
  _REAL_,  intent(in) :: bv(1:3*numatoms)
  _REAL_,  intent(inout) :: xv(1:3*numatoms)

  logical :: unconvg

  unconvg = .true.
  iter = 0.0d0

  ! main loop
  do while ( unconvg )

    ! compute T_long*xv for local_precondition()
    ! compute residue rv
    ! Note that A's off-diagonal elements are those of the dipole-dipole tensor and
    ! its diagonal elements are 1/pol (1/polarizability)
    call pGM_NonBond_dip_dip_fields(numatoms, x, xv, zv)
    call pGM_NonBond_dip_dip_fields_short(numatoms, x, xv, zu)
    pv = zv - zu ! this is T_long*xv
    zv = zv + ip*xv ! this is A*xv
    rv = bv - zv
    r_norm = maxval(abs(rv/bv))
    ! check convergence
    if ( pol_gauss_verbose == 2 ) write (6, '(a,i3,a,f14.10)') &
            ' pGM ESOR: iter = ', int(iter), ', residue norm: ', r_norm
    if ( int(iter) >= scf_sor_niter .or. r_norm <= dipole_scf_tol ) then
      unconvg = .false.
      if ( int(iter) >= scf_sor_niter ) then
        write(6,'(a)') ' pGM ESOR: FATAL maxitn exceeded!'
        call mexit(6,1)
      end if
    end if

    iter = iter + 1.0d0

    ! Set up rhs as E_perm-T_long*xv
    pv = bv - pv
    ! The preconditioner is solving (T_short+1/pol)*xv = E_perm-T_long*xv
    call local_precondition(numatoms, x, pv, xv, scf_local_niter)

  end do ! do while (unconvg)

  if ( pol_gauss_verbose /= 0 )then
    write(6,'(/a,i5,1x,f14.9/)') ' pGM ESOR: out of loop, iter, norm = ',int(iter),r_norm
  end if

end subroutine pGM_esor
!----------------------------------------------------------
! The standard CG for symmetric positive-definite matrix
! will be implemented here as pGM_cg().
!
! Working arrays are all declared as 1-d for vectoerization.
!
! The notation will also follow that of standard practice
! See inline comments on how related vectors are obtained.
!
! The original two dipole systems are no longer necessary so
! only the d system will be retained.
!----------------------------------------------------------
subroutine pGM_cg(numatoms, x, r_norm, iter, bv, xv)
  use pol_gauss_mdin, only : use_average,scf_cg_niter,dipole_scf_tol,pol_gauss_verbose
  implicit none

  integer, intent(in) :: numatoms
  _REAL_,  intent(in) :: x(*)
  _REAL_,  intent(out) :: r_norm, iter
  _REAL_,  intent(in) :: bv(1:3*numatoms)
  _REAL_,  intent(inout) :: xv(1:3*numatoms)

  logical :: unconvg
  _REAL_ :: alpha, beta, rdotr1, rdotr2, pdotz
  _REAL_ :: bnorm, a_norm

  bnorm = sum(abs(bv))

  ! r(0) = b(0) - A*x(0)
  ! Note that A's off-diagonal elements are those of the dipole-dipole tensor and
  ! its diagonal elements are 1/pol (1/polarizability)
  call pGM_NonBond_dip_dip_fields(numatoms, x, xv, zv)
  zv = zv + ip*xv
  rv = bv - zv

  ! p(0) = r(0)
  pv = rv

  ! iteration 0:
  ! compute <r(0),r(0)>
  iter = 0.0d0
  rdotr1 = dot_product(rv, rv)

  ! main loop
  unconvg = .true.
  do while (unconvg)

    ! compute Ap(i) = A * p(i)
    ! compute alpha(i) = <r(i),r(i)>/<p(i),Ap(i)>
    call pGM_NonBond_dip_dip_fields(numatoms, x, pv, zv)
    zv = zv + ip*pv
    pdotz = dot_product(pv, zv)
    alpha = rdotr1/pdotz

    ! iteration i:
    iter = iter + 1.0d0

    ! update x(i) = x(i-1) + alpha(i) p(i-1)
    !        r(i) = r(i-1) - alpha(i) Ap(i-1)
    xv = xv + alpha*pv
    rv = rv - alpha*zv
    if ( use_average == 1 ) then
      r_norm = sum(abs(rv))/bnorm
    else
      r_norm = maxval(abs(rv/bv))
      a_norm = sum(abs(rv))/bnorm
    end if

    ! check convergence
    if ( pol_gauss_verbose == 2 ) write (6, '(a,i3,a,f14.10)') &
            ' pGM CG: iter = ', int(iter), ', residue norm: ', r_norm
    if ( iter >= scf_cg_niter .or. r_norm <= dipole_scf_tol ) then
      unconvg = .false.
      if ( iter >= scf_cg_niter ) then
        write(6,'(a)') ' pGM CG: FATAL maxitn exceeded!'
        call mexit(6,1)
      end if
      exit
    end if

    ! compute beta(i) = <r(i),r(i)>/<r(i-1),r(i-1)>
    rdotr2 = dot_product(rv, rv)
    beta = rdotr2/rdotr1
    rdotr1 = rdotr2

    ! update p(i) = r(i) + beta(i) p(i-1)
    pv = rv + beta*pv

  end do ! do while (unconvg)

  if ( pol_gauss_verbose /= 0 ) then
    !write(6,'(/a,i5,1x,f14.9/)') ' pGM CG: out of loop, iter, rel norm = ', int(iter), r_norm
    write(6,'(/a,i5,1x,2e20.9/)') ' pGM CG: out of loop, iter, rel norm = ', int(iter), r_norm, a_norm
  end if

end subroutine pGM_cg
!----------------------------------------------------------
! The preconditioned CG for symmetric positive-definite matrix
! will be implemented here as pGM_pcg(). The preconditioner
! is simply the short-range dipole-dipole tensor elements and
! the diagonal elements of the A matrix.
!
! Working arrays are all declared as 1-d for vectoerization.
!
! The notation will also follow that of standard practice
! See inline comments on how related vectors are obtained.
!
! The original two dipole systems are no longer necessary so
! only the d system will be retained.
!----------------------------------------------------------
subroutine pGM_pcg(numatoms, x, r_norm, iter, bv, xv)
  use pol_gauss_mdin, only : scf_cg_niter,scf_local_niter,dipole_scf_tol,pol_gauss_verbose
  implicit none

  integer, intent(in) :: numatoms
  _REAL_,  intent(in) :: x(*)
  _REAL_,  intent(out) :: r_norm, iter
  _REAL_,  intent(in) :: bv(1:3*numatoms)
  _REAL_,  intent(inout) :: xv(1:3*numatoms)

  logical :: unconvg
  _REAL_ :: alpha, beta, rdotr1, rdotr2, pdotz
  _REAL_ :: bnorm, a_norm

  bnorm = sum(abs(bv))
  if ( pol_gauss_verbose /= 0 ) write(6,'(a, f14.9)') 'initial dip norm', sum(abs(xv))

  ! r(0) = b(0) - A*x(0)
  ! Note that A's off-diagonal elements are those of the dipole-dipole tensor and
  ! its diagonal elements are 1/polarizability
  call pGM_NonBond_dip_dip_fields(numatoms, x, xv, zv)
  zv = zv + ip*xv
  rv = bv - zv

  ! apply preconditoner to rv
  ! this means solving a new linear system of M * z = r, i.e. r is the rhs
  ! z(0) = M^-1 * r(0)
  ! z(0) is initialized as pol*r(0)
  !zv = 0.0d0
  zv = ap*rv
  call local_precondition(numatoms, x, rv, zv, scf_local_niter)

  ! p(0) = z(0)
  pv = zv

  !  iteration 0:
  !  compute <r(0),z(0)>
  iter = 0.0d0
  rdotr1 = dot_product(rv, zv)

  ! main loop
  unconvg = .true.
  do while (unconvg)

    ! compute Ap(i) = A * p(i)
    ! compute alpha(i) = <r(i),z(i)>/<p(i),Ap(i)>
    call pGM_NonBond_dip_dip_fields(numatoms, x, pv, zv)
    zv = zv + ip*pv
    pdotz = dot_product(pv, zv)
    alpha = rdotr1/pdotz

    ! iteration i:
    iter = iter + 1.0d0

    ! update x(i) = x(i-1) + alpha(i) p(i-1)
    !        r(i) = r(i-1) - alpha(i) Ap(i-1)
    xv = xv + alpha*pv
    rv = rv - alpha*zv
    r_norm = maxval(abs(rv/bv))
    a_norm = sum(abs(rv))/bnorm

    ! check convergence
    if ( pol_gauss_verbose == 2 ) write (6, '(a,i3,a,e14.9)') &
            ' pGM PCG: iter = ', int(iter), ', residue norm: ', r_norm
    if ( iter >= scf_cg_niter .or. r_norm <= dipole_scf_tol ) then
      unconvg = .false.
      if ( iter >= scf_cg_niter ) then
        write(6,'(a)') ' pGM PCG: FATAL maxitn exceeded!'
        call mexit(6,1)
      end if
      exit
    end if

    ! apply preconditoner to rv
    ! z(i) = M^-1 * r(i)
    ! z(i) is initialized as pol*r(i)
    zv = ap*rv
    call local_precondition(numatoms, x, rv, zv, scf_local_niter)

    ! compute beta(i) = <r(i),z(i)>/<r(i-1),z(i-1)>
    rdotr2 = dot_product(rv, zv)
    beta = rdotr2/rdotr1
    rdotr1 = rdotr2

    ! update p(i) = r(i) + beta(i) p(i-1)
    pv = zv + beta*pv

  end do ! do while (unconvg)

  if ( pol_gauss_verbose /= 0 )then
    !write(6,'(/a,i5,1x,f14.9/)') ' pGM PCG: out of loop, iter, rel norm = ', int(iter), r_norm
    write(6,'(/a,i5,1x,2e20.9/)') ' pGM PCG: out of loop, iter, rel norm = ', int(iter), r_norm, a_norm
  endif

end subroutine PGM_pcg
!----------------------------------------------------------
! The preconditioned&peeked CG for symmetric positive-definite matrix
! will be implemented here as pGM_pkcg(). The preconditioner
! is simply the short-range dipole-dipole tensor elements and
! the diagonal elements of the A matrix. Here a peek step is
! implemented as one extra SOR step to save one full CG iteration
! worth of work for majority of the cases.
!
! Working arrays are all declared as 1-d for vectoerization.
!
! The notation will also follow that of standard practice
! See inline comments on how related vectors are obtained.
!
! The original two dipole systems are no longer necessary so
! only the d system will be retained.
!----------------------------------------------------------
subroutine pGM_pkcg(numatoms, x, r_norm, iter, bv, xv)
  use pol_gauss_mdin, only : scf_cg_niter,scf_local_niter,scf_sor_coefficient,&
                             dipole_scf_tol,pol_gauss_verbose
  implicit none

  integer, intent(in) :: numatoms
  _REAL_,  intent(in) :: x(*)
  _REAL_,  intent(out) :: r_norm, iter
  _REAL_,  intent(in) :: bv(1:3*numatoms)
  _REAL_,  intent(inout) :: xv(1:3*numatoms)

  logical :: unconvg
  _REAL_ :: alpha, beta, rdotr1, rdotr2, pdotz
  _REAL_ :: bnorm, a_norm

  bnorm = sum(abs(bv))
  if ( pol_gauss_verbose /= 0 ) write(6,'(a, f14.9)') 'initial dip norm', sum(abs(xv))

  ! r(0) = b(0) - A*x(0)
  ! Note that A's off-diagonal elements are those of the dipole-dipole tensor and
  ! its diagonal elements are 1/polarizability
  call pGM_NonBond_dip_dip_fields(numatoms, x, xv, zv)
  zv = zv + ip*xv
  rv = bv - zv

  ! apply preconditoner to rv
  ! this means solving a new linear system of M * z = r, i.e. r is the rhs
  ! z(0) = M^-1 * r(0)
  ! z(0) is initialized as pol*r(0)
  zv = ap*rv
  call local_precondition(numatoms, x, rv, zv, scf_local_niter)

  ! p(0) = z(0)
  pv = zv

  !  iteration 0:
  !  compute <r(0),z(0)>
  iter = 0.0d0
  rdotr1 = dot_product(rv, zv)

  ! main loop
  unconvg = .true.
  do while (unconvg)

    ! compute Ap(i) = A * p(i)
    ! compute alpha(i) = <r(i),z(i)>/<p(i),Ap(i)>
    call pGM_NonBond_dip_dip_fields(numatoms, x, pv, zv)
    zv = zv + ip*pv
    pdotz = dot_product(pv, zv)
    alpha = rdotr1/pdotz

    ! iteration i:
    iter = iter + 1.0d0

    ! update x(i) = x(i-1) + alpha(i) p(i-1)
    !        r(i) = r(i-1) - alpha(i) Ap(i-1)
    xv = xv + alpha*pv
    rv = rv - alpha*zv
    r_norm = maxval(abs(rv/bv))
    a_norm = sum(abs(rv))/bnorm

    ! check convergence
    if ( pol_gauss_verbose == 2 ) write (6, '(a,i3,a,f14.10)') &
            ' pGM PKCG: iter = ', int(iter), ', residue norm: ', r_norm
    if ( iter >= scf_cg_niter .or. r_norm <= dipole_scf_tol ) then
      unconvg = .false.
      if ( iter >= scf_cg_niter ) then
        write(6,'(a)') ' pGM PKCG: FATAL maxitn exceeded!'
        call mexit(6,1)
      end if
      exit
    end if

    ! apply preconditoner to rv
    ! z(i) = M^-1 * r(i)
    ! z(i) is initialized as pol*r(i)
    zv = ap*rv
    call local_precondition(numatoms, x, rv, zv, scf_local_niter)

    ! compute beta(i) = <r(i),z(i)>/<r(i-1),z(i-1)>
    rdotr2 = dot_product(rv, zv)
    beta = rdotr2/rdotr1
    rdotr1 = rdotr2

    ! update p(i) = r(i) + beta(i) p(i-1)
    pv = zv + beta*pv

  end do ! do while (unconvg)

  if ( pol_gauss_verbose > 0 ) then
    !write(6,'(/a,i5,1x,f14.9)') ' pGM PKCG: out of loop, iter, rel norm = ', int(iter), r_norm
    write(6,'(/a,i5,1x,2e20.9/)') ' pGM PCG: out of loop, iter, rel norm = ', int(iter), r_norm, a_norm
  end if

  ! peek using a SOR step to save some time
  ! a) Haixin's ESOR step
  ! The following call sets zv = -E^pol_short
  ! Note also xv*ip = E^tot = E^perm + E^pol_short + E^pol_long
  ! so xv*ip + zv = E^perm + E^pol_long, which is the rhs of
  ! the linear system in ESOR.
  ! The difference here is that we are trying to solve
  ! A*xv^i = bv + rv^i, which is different from the original
  ! target equation A*xv^i = bv in ESOR, so the rhs here is
  ! E^perm + E^pol_long + rv = zv + xv*ip + rv
  !call pGM_NonBond_dip_dip_fields_short(numatoms, x, xv, zv)
  !zv = zv + xv*ip + rv
  !call local_precondition(numatoms, x, zv, xv, 3)

  ! b) Standard SOR step seems to be good enough
  xv = xv + scf_sor_coefficient*rv*ap

  ! let's check to see how good it is ...
  if ( pol_gauss_verbose /= 0 ) then
    call pGM_NonBond_dip_dip_fields(numatoms, x, xv, zv)
    zv = zv + ip*xv
    rv = bv - zv
    r_norm = maxval(abs(rv/bv))
    a_norm = sum(abs(rv))/bnorm
    write(6,'(/a,2e20.9/)') ' pGM PKCG: after peek, rel norm = ', r_norm, a_norm
  end if

end subroutine pGM_pkcg
!----------------------------------------------------------
! Preconditioning using short-range dip-dip tensor elements only
! This is simply another CG solver with limited iteration, so
! convergence check is not necessary.
! The job type has been removed by simply passing in the appopriate
! rhs vector bv and unknown vector xv() as in AM * xv = bv, or
! xv = AM^-1 * bv
!----------------------------------------------------------
subroutine local_precondition(numatoms, x, bv, xv, niter)
  use pol_gauss_mdin, only : pol_gauss_verbose, scf_sor_coefficient

  integer :: numatoms, niter
  _REAL_, intent(in) :: x(*)
  _REAL_ :: bv(1:3*numatoms), xv(1:3*numatoms)

  _REAL_ :: alpha, beta, rdotr1, rdotr2, pdotz, r_norm
  integer :: i

  ! r(0) = b(0) - A*x(0)
  ! Note that A's off-diagonal elements are those of the dipole-dipole tensor
  ! and its diagonal elements are 1/pol (1/polarizability)
  call pGM_NonBond_dip_dip_fields_short(numatoms, x, xv, zu)
  zu = zu + ip*xv
  ru = bv - zu

  ! do a quick diagonal scaling as preconditioning
  zu = ap*ru

  ! p(0) = z(0)
  pu = zu

  !  iteration 0:
  !  compute <r(0),z(0)>
  i = 0
  rdotr1 = dot_product(ru, zu)

  ! main loop
  do i = 1, niter

    ! compute Ap(i) = A * p(i)
    ! compute alpha(i) = <r(i),r(i)>/<p(i),Ap(i)>
    call pGM_NonBond_dip_dip_fields_short(numatoms, x, pu, zu)
    zu = zu + ip*pu
    pdotz = dot_product(pu, zu)
    alpha = rdotr1/pdotz

    ! update x(i) = x(i-1) + alpha(i) p(i-1)
    !        r(i) = r(i-1) - alpha(i) Ap(i-1)
    xv = xv + alpha*pu
    ru = ru - alpha*zu

    ! do not check convergence, just for info only
    if ( pol_gauss_verbose == 2 ) then
      r_norm = maxval(abs(ru/bv))
      write (6, '(a,i3,a,f14.10)') ' pGM precond: iter = ', i, ', residue norm: ', r_norm
    end if

    ! do a quick diagonal scaling as preconditioning
    zu = ap*ru

    ! compute beta(i) = <r(i),z(i)>/<r(i-1),z(i-1)>
    rdotr2 = dot_product(ru, zu)
    beta = rdotr2/rdotr1
    rdotr1 = rdotr2

    ! update p(i) = z(i) + beta(i) p(i-1)
    pu = zu + beta*pu

  end do ! do while (unconvg)

  ! Also do a peek step here with the the undamped SOR step
  xv = xv + ru*ap

end subroutine local_precondition
!----------------------------------------------------------
end module pol_gauss_induced
