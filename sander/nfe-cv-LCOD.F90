! <compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

!
! LCOD (Linear Combination Of Distances)
!
! cv%i = (a11, a12, a21, a22, ..., aN1, aN2)
!
!     indexes of the participating atoms
!
! cv%r = (r1, r2, ..., rN)
!
!     non-zero weights
!
! value = r1*d1 + r2*d2 + ... + rN*dN
!
!     d1 -- distance between atoms #a11 and #a12
!     d2 -- distance between atoms #a21 and #a22
!                       . . .
!     dN -- distance between atoms #aN1 and #aN2
!

module nfe_cv_LCOD

!=============================================================================

implicit none

private

!=============================================================================

public :: colvar_value
public :: colvar_force

public :: colvar_bootstrap
public :: print_details

!=============================================================================

contains

!=============================================================================

function colvar_value(cv, x) result(value)

#  ifndef NFE_DISABLE_ASSERT
   use nfe_utils
   use nfe_sander_proxy
#  endif /* NFE_DISABLE_ASSERT */

   use nfe_constants
   use nfe_colvar_type
   use nfe_colvar_math

   implicit none

   NFE_REAL :: value

   type(colvar_t), intent(in) :: cv

   NFE_REAL, intent(in) :: x(*)

#  include "nfe-mpi.h"

   integer :: a1, a2, n

   nfe_assert(cv%type.eq.COLVAR_LCOD)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   nfe_assert(size(cv%i).gt.0)
   nfe_assert(size(cv%r).gt.0)

   nfe_assert(mod(size(cv%i), 2).eq.0)
   nfe_assert((size(cv%i)/2).eq.size(cv%r))

   value = ZERO

   NFE_MASTER_ONLY_BEGIN

      do n = 1, size(cv%r)

         nfe_assert(cv%i(2*n - 1).gt.0.and.cv%i(2*n - 1).le.sander_natoms())
         a1 = 3*cv%i(2*n - 1) - 2

         nfe_assert(cv%i(2*n).gt.0.and.cv%i(2*n).le.sander_natoms())
         a2 = 3*cv%i(2*n) - 2

         nfe_assert(a1.ne.a2)

         value = value + cv%r(n)*distance(x(a1:a1 + 2), x(a2:a2 + 2))

      end do

   NFE_MASTER_ONLY_END

end function colvar_value

!=============================================================================

subroutine colvar_force(cv, x, fcv, f)

#  ifndef NFE_DISABLE_ASSERT
   use nfe_utils
   use nfe_sander_proxy
#  endif /* NFE_DISABLE_ASSERT */

   use nfe_colvar_type
   use nfe_colvar_math

   implicit none

   type(colvar_t), intent(in) :: cv

   NFE_REAL, intent(in) :: x(*), fcv

   NFE_REAL, intent(inout) :: f(*)

#  include "nfe-mpi.h"

   NFE_REAL :: d1(3), d2(3)

   integer :: a1, a2, n

   nfe_assert(cv%type.eq.COLVAR_LCOD)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   nfe_assert(size(cv%i).gt.0)
   nfe_assert(size(cv%r).gt.0)

   nfe_assert(mod(size(cv%i), 2).eq.0)
   nfe_assert((size(cv%i)/2).eq.size(cv%r))

   NFE_MASTER_ONLY_BEGIN

      do n = 1, size(cv%r)

         nfe_assert(cv%i(2*n - 1).gt.0.and.cv%i(2*n - 1).le.sander_natoms())
         a1 = 3*cv%i(2*n - 1) - 2

         nfe_assert(cv%i(2*n).gt.0.and.cv%i(2*n).le.sander_natoms())
         a2 = 3*cv%i(2*n) - 2

         nfe_assert(a1.ne.a2)

         call distance_d(x(a1:a1 + 2), x(a2:a2 + 2), d1, d2)

         f(a1:a1 + 2) = f(a1:a1 + 2) + cv%r(n)*fcv*d1
         f(a2:a2 + 2) = f(a2:a2 + 2) + cv%r(n)*fcv*d2

      end do

   NFE_MASTER_ONLY_END

end subroutine colvar_force

!=============================================================================

subroutine colvar_bootstrap(cv, cvno)

   use nfe_utils
   use nfe_constants
   use nfe_colvar_type
   use nfe_colvar_utils
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno

#  include "nfe-mpi.h"

   integer :: n, j

   nfe_assert(cv%type.eq.COLVAR_LCOD)

   if (.not.associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%i)

   if (.not.associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : no reals found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%r)

   if (mod(size(cv%i), 2).ne.0) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : odd number of integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! mod(size(cv%i), 2).ne.0

   if (size(cv%i).lt.2) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%i).lt.2

   if (size(cv%r).ne.size(cv%i)/2) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : size(cv%r).ne.size(cv%i)/2'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r).ne.size(cv%i)/2

   do n = 1, size(cv%r)
      if (abs(cv%r(n)).lt.1.0D-8) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(n)//',a/)') &
               NFE_ERROR, 'CV #', cvno, ' (LCOD) : real number #', n, &
               ' is too small'
         NFE_MASTER_ONLY_END
         call terminate()
      end if ! abs(cv%r(n)).lt.1.0D-8

      do j = 0, 1
         if (cv%i(2*n - j).lt.1.or.cv%i(2*n - j).gt.sander_natoms()) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(2*n - j)//',a,'//pfmt &
                  (cv%i(2*n - j))//',a,'//pfmt(sander_natoms())//',a/)') &
                  NFE_ERROR, 'CV #', cvno, ' (LCOD) : integer #', &
                  (2*n - j), ' (', cv%i(2*n - j), ') is out of range [1, ', &
                  sander_natoms(), ']'
            NFE_MASTER_ONLY_END
            call terminate()
         end if
      end do

      if (cv%i(2*n - 1).eq.cv%i(2*n)) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(2*n - 1)//',a,'//pfmt &
               (2*n)//',a,'//pfmt(cv%i(2*n))//',a/)') NFE_ERROR, 'CV #', &
               cvno, ' (LCOD) : integers #', (2*n - 1), ' and #', (2*n), &
               ' are equal (', cv%i(2*n), ')'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

end subroutine colvar_bootstrap

!=============================================================================

subroutine print_details(cv, lun)

   use nfe_utils
   use nfe_colvar_type
   use nfe_colvar_utils
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: n, a1, a2
   character(4) :: aname1, aname2

   nfe_assert(is_master())
   nfe_assert(cv%type.eq.COLVAR_LCOD)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   nfe_assert(size(cv%i).gt.0)
   nfe_assert(size(cv%r).gt.0)

   nfe_assert(mod(size(cv%i), 2).eq.0)
   nfe_assert((size(cv%i)/2).eq.size(cv%r))

   do n = 1, size(cv%r)

      a1 = cv%i(2*n - 1)
      a2 = cv%i(2*n)

      aname1 = sander_atom_name(a1)
      aname2 = sander_atom_name(a2)

      write (unit = lun, fmt = '(a,4x,f8.3,a,'//pfmt &
         (a1)//',a,a,a,'//pfmt(a2)//',a,a,a)') NFE_INFO, cv%r(n), ' * (', &
         a1, ' [', trim(aname1), '] <=> ', a2, ' [', trim(aname2), '])'

   end do

end subroutine print_details

!=============================================================================

end module nfe_cv_LCOD
