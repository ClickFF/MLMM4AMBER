! <compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

!
! input:
!
! cv%i = (a1, a2, a3, 0, b1, b2, b3, b4, 0, ..., c1, c2, c3, 0)
!
!     (a[1-3] - 1st group, b[1-4] - 2nd group, etc; an atom may
!      enter a few groups simultaneously; last zero is optional;
!      empty groups [e.g., 2+ zeros in a row] are not allowed)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, b1x, ...)
!
!        (reference coordinates without '0' sentinel(s))
!
! value = sqrt[(M_1*rmsd1^2 + ... + M_N*rmsdN^2)/(M_1 + ... + M_N)]
!         M_1 - mass of group 1, ..., M_N - mass of group N
!

module nfe_cv_MULTI_RMSD

!=============================================================================

implicit none

private

!=============================================================================

!
! V-table
!

public :: colvar_value
public :: colvar_force

public :: colvar_bootstrap
public :: colvar_cleanup

public :: print_details

!=============================================================================

type, private :: group_t

   integer :: i0, i1, r0

   NFE_REAL, pointer :: mass(:) => null()

   NFE_REAL, pointer :: cm_crd(:) => null()
   NFE_REAL, pointer :: ref_crd(:) => null()

   NFE_REAL :: total_mass
   NFE_REAL :: ref_nrm
   NFE_REAL :: quaternion(4)
   NFE_REAL :: rmsd2

end type group_t

private :: group_bootstrap
private :: group_finalize
private :: group_print
private :: group_evaluate

type, private :: priv_t
   NFE_REAL :: value
   NFE_REAL :: total_mass
   type(group_t), pointer :: groups(:) => null()
#ifdef MPI
   integer :: first_cpu
#endif /* MPI */
#  include "nfe-cv-priv.type"
end type priv_t

#include "nfe-cv-priv.decl"

!=============================================================================

contains

!=============================================================================

#include "nfe-cv-priv.impl"

!=============================================================================

subroutine group_evaluate(cv, grp, x)

   use nfe_rmsd
   use nfe_utils
   use nfe_constants
   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in)    :: cv
   type(group_t),  intent(inout) :: grp

   NFE_REAL, intent(in) :: x(*)

   NFE_REAL :: cm(3), cm_nrm, lambda
   integer :: i, n
! gfortran 4.2 and lower does not support the volatile keyword, but this is
! necessary to work around an intel compiler optimization bug. So turn off
! the volatile keyword for qualifying compilers
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
   integer, volatile :: a, a3 ! work around ifort optimization bug
#else
   integer :: a, a3
#endif /* __GFORTRAN__ */

   ! compute the center of mass (of the moving atoms)
   cm = ZERO
   n = 1
   do a = grp%i0, grp%i1
      a3 = 3*cv%i(a)
      cm = cm + grp%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/grp%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = grp%i0, grp%i1
      a3 = 3*(n - 1)
      do i = 1, 3
         grp%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + grp%mass(n)*grp%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call rmsd_q(grp%i1 - grp%i0 + 1, grp%mass, &
         grp%cm_crd, grp%ref_crd, lambda, grp%quaternion)

   grp%rmsd2 = max(ZERO, ((grp%ref_nrm + cm_nrm) - 2*lambda))

end subroutine group_evaluate

!=============================================================================

function colvar_value(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_constants
   use nfe_colvar_type

   implicit none

   NFE_REAL :: value

   type(colvar_t), intent(inout) :: cv

   NFE_REAL, intent(in) :: x(*)

#  include "nfe-mpi.h"

   type(priv_t), pointer :: priv
   integer :: g

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)

   priv => get_priv(cv)
   nfe_assert(associated(priv%groups))
   nfe_assert(size(priv%groups).gt.0)

   priv%value = ZERO

   do g = 1, size(priv%groups)
#ifdef MPI
      if (mod(g + priv%first_cpu - 1, sandersize).eq.sanderrank) then
#endif /* MPI */
      call group_evaluate(cv, priv%groups(g), x)
      priv%value = priv%value + priv%groups(g)%rmsd2
#ifdef MPI
      endif
#endif /* MPI */
   end do

#ifdef MPI
   call mpi_reduce(priv%value, value, 1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, 0, commsander, error)
   nfe_assert(error.eq.0)
   priv%value = sqrt(value/priv%total_mass)
#else
   priv%value = sqrt(priv%value/priv%total_mass)
#endif /* MPI */

   value = priv%value

end function colvar_value

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine colvar_force(cv, fcv, f)

   use nfe_rmsd
   use nfe_utils
   use nfe_constants
   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   NFE_REAL, intent(in) :: fcv
   NFE_REAL, intent(inout) :: f(*)

#  include "nfe-mpi.h"

   integer :: n, g, a, a3, f3
   NFE_REAL :: U(3,3), tmp
   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)
   nfe_assert(associated(cv%i))

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%groups))
   nfe_assert(size(priv%groups).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, commsander, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

   do g = 1, size(priv%groups)

#ifdef MPI
      if (mod(g + priv%first_cpu - 1, sandersize).ne.sanderrank) &
         cycle
#endif /* MPI */

      call rmsd_q2u(priv%groups(g)%quaternion, U)
      tmp = fcv/max(priv%value, NFE_TO_REAL(0.000001))/priv%total_mass

      n = 1
      do a = priv%groups(g)%i0, priv%groups(g)%i1
         a3 = 3*n
         f3 = 3*cv%i(a)

         f(f3 - 2:f3) = f(f3 - 2:f3) &
            + tmp*priv%groups(g)%mass(n)*(priv%groups(g)%cm_crd(a3 - 2:a3) &
               - matmul(U, priv%groups(g)%ref_crd(a3 - 2:a3)))
         n = n + 1
      end do
   end do

#ifdef MPI
   priv%first_cpu = mod(priv%first_cpu + 1, sandersize)
#endif /* MPI */

end subroutine colvar_force

!=============================================================================

subroutine group_bootstrap(grp, cv, cvno, amass, i0, i1, r0)

   use nfe_utils
   use nfe_constants
   use nfe_colvar_type
   use nfe_sander_proxy

   implicit none

   type(group_t),  intent(inout) :: grp

   type(colvar_t), intent(in) :: cv
   integer,        intent(in) :: cvno
   NFE_REAL,      intent(in) :: amass(*)

   integer, intent(in) :: i0, i1, r0

#  include "nfe-mpi.h"

   integer :: a, b, n_atoms, error
   NFE_REAL :: cm(3)

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)
   nfe_assert(associated(cv%i).and.associated(cv%r))
   nfe_assert(i0.gt.0.and.i0.le.size(cv%i))
   nfe_assert(i1.gt.0.and.i1.le.size(cv%i))
   nfe_assert(r0.gt.0.and.r0.lt.size(cv%r))

   ! basic checks
   if (i0 + 2 > i1) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (i0)//',a,'//pfmt(i1)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (MULTI_RMSD) : too few integers in group (', i0, ':', i1, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = i0, i1
      nfe_assert(a .le. size(cv%i))
      if (cv%i(a) < 1 .or. cv%i(a) > sander_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(sander_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (MULTI_RMSD) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', sander_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = i0, i1
      do b = a + 1, i1
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (MULTI_RMSD) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! cv%i(a) == cv%i(b)
      end do
   end do

   ! allocate/setup

   n_atoms = i1 - i0 + 1

   grp%i0 = i0
   grp%i1 = i1
   grp%r0 = r0

   allocate(grp%mass(n_atoms), grp%cm_crd(3*n_atoms), &
            grp%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   grp%total_mass = ZERO

   do a = 1, n_atoms
      grp%mass(a) = amass(cv%i(a + grp%i0 - 1))
      grp%total_mass = grp%total_mass + grp%mass(a)
      cm = cm + grp%mass(a)*cv%r(r0 + 3*(a - 1):r0 + 3*a - 1)
   end do

   cm = cm/grp%total_mass

   ! translate reference coordinates to CM frame
   grp%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         grp%ref_crd(3*(a - 1) + b) = cv%r(r0 + 3*(a - 1) + b - 1) - cm(b)
         grp%ref_nrm = grp%ref_nrm + grp%mass(a)*grp%ref_crd(3*(a - 1) + b)**2
      end do
   end do

end subroutine group_bootstrap

!=============================================================================

subroutine group_finalize(grp)

   NFE_USE_AFAILED

   implicit none

   type(group_t), intent(inout) :: grp

   nfe_assert(associated(grp%mass))
   nfe_assert(associated(grp%cm_crd))
   nfe_assert(associated(grp%ref_crd))

   deallocate(grp%mass, grp%cm_crd, grp%ref_crd)

end subroutine group_finalize

!=============================================================================

subroutine group_print(cv, grp, lun)

   use nfe_utils
   use nfe_colvar_type
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(in) :: cv
   type(group_t),  intent(in) :: grp
   integer,        intent(in) :: lun

   integer :: a, c
   character(4) :: aname

   nfe_assert(associated(cv%i).and.associated(cv%r))
   nfe_assert(grp%i0.gt.0.and.grp%i0.le.size(cv%i))
   nfe_assert(grp%i1.gt.0.and.grp%i1.le.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = grp%i0, grp%i1

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= sander_natoms())
      aname = sander_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == grp%i1) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 0
   do a = grp%i0, grp%i1
      nfe_assert(grp%r0 + c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(grp%r0 + c), ', ', cv%r(grp%r0 + c + 1), &
         ', ', cv%r(grp%r0 + c + 2)
      c = c + 3
   end do

end subroutine group_print

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   use nfe_utils
   use nfe_constants
   use nfe_colvar_type
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   NFE_REAL,      intent(in)    :: amass(*)

   integer :: i, i0, n_atoms, n_groups, error

   type(priv_t), pointer :: priv

#  include "nfe-mpi.h"

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (MULTI_RMSD) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (MULTI_RMSD) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   ! count the groups (number of zeros in the cv%i array)
   n_atoms = 0
   n_groups = 0
   i0 = 1

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         n_groups = n_groups + 1
         if (i .eq. i0) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (MULTI_RMSD) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! i .eq. i0
         i0 = i + 1
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (MULTI_RMSD) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (size(cv%i) .gt. 0 .and. cv%i(size(cv%i)) .gt. 0) &
      n_groups = n_groups + 1

   nfe_assert(n_groups .gt. 0)

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (MULTI_RMSD) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r) /= 3*n_atoms

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! allocate/setup groups
   allocate(priv%groups(n_groups), stat = error)
   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   i0 = 1 ! first atom
   n_atoms = 1
   n_groups = 0

   priv%total_mass = ZERO

   do i = 1, size(cv%i)
      if (cv%i(i) .eq. 0) then
         n_groups = n_groups + 1
         nfe_assert(n_groups .le. size(priv%groups))
         call group_bootstrap(priv%groups(n_groups), &
                              cv, cvno, amass, i0, i - 1, n_atoms)
         n_atoms = n_atoms + 3*(i - i0)
         i0 = i + 1
         priv%total_mass = priv%total_mass &
            + priv%groups(n_groups)%total_mass
      end if
   end do

   if (size(cv%i) .gt. 0 .and. cv%i(size(cv%i)) .gt. 0) then
      n_groups = n_groups + 1
      nfe_assert(n_groups .le. size(priv%groups))
      call group_bootstrap(priv%groups(n_groups), &
                           cv, cvno, amass, i0, size(cv%i), n_atoms)
      priv%total_mass = priv%total_mass &
         + priv%groups(n_groups)%total_mass
   end if

#ifdef MPI
   priv%first_cpu = 0
#endif /* MPI */

end subroutine colvar_bootstrap

!=============================================================================

subroutine colvar_cleanup(cv)

   NFE_USE_AFAILED

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(inout) :: cv

   integer :: g
   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)
   nfe_assert(associated(priv%groups))

   do g = 1, size(priv%groups)
      call group_finalize(priv%groups(g))
   end do

   deallocate(priv%groups)
   call del_priv(cv)

end subroutine colvar_cleanup

!=============================================================================

subroutine print_details(cv, lun)

   use nfe_utils
   use nfe_colvar_type
   use nfe_colvar_utils
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: g
   type(priv_t), pointer :: priv

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_MULTI_RMSD)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(associated(priv%groups))

   do g = 1, size(priv%groups)
      write (unit = lun, fmt = '(a,a,'//pfmt(g)//',a)') &
         NFE_INFO, '<> group <> #', g, ':'
      call group_print(cv, priv%groups(g), lun)
   end do

end subroutine print_details

!=============================================================================

end module nfe_cv_MULTI_RMSD
