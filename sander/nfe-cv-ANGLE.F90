! <compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

!
! ANGLE "subclass"
!

module nfe_cv_ANGLE

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

#ifndef NFE_DISABLE_ASSERT
   use nfe_utils
   use nfe_sander_proxy
#endif /* NFE_DISABLE_ASSERT */

   use nfe_constants
   use nfe_colvar_type
   use nfe_colvar_math

   implicit none

   NFE_REAL :: value

   type(colvar_t), intent(in) :: cv

   NFE_REAL, intent(in) :: x(*)

#include "nfe-mpi.h"

   integer :: a1, a2, a3

   nfe_assert(cv%type.eq.COLVAR_ANGLE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).eq.3)

   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= sander_natoms())
   a1 = 3*cv%i(1) - 2

   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= sander_natoms())
   a2 = 3*cv%i(2) - 2

   nfe_assert(cv%i(3) > 0 .and. cv%i(3) <= sander_natoms())
   a3 = 3*cv%i(3) - 2

   nfe_assert(a1 /= a2 .and. a1 /= a3 .and. a2 /= a3)

#ifdef MPI
   if (sanderrank.eq.0) then
#endif /* MPI */
      value = angle(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2))
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function colvar_value

!=============================================================================

subroutine colvar_force(cv, x, fcv, f)

#ifndef NFE_DISABLE_ASSERT
   use nfe_utils
   use nfe_sander_proxy
#endif /* NFE_DISABLE_ASSERT */

   use nfe_colvar_type
   use nfe_colvar_math

   implicit none

   type(colvar_t), intent(in) :: cv

   NFE_REAL, intent(in) :: x(*), fcv

   NFE_REAL, intent(inout) :: f(*)

#  include "nfe-mpi.h"

   NFE_REAL :: d1(3), d2(3), d3(3)

   integer :: a1, a2, a3

   nfe_assert(cv%type == COLVAR_ANGLE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i) == 3)

   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= sander_natoms())
   a1 = 3*cv%i(1) - 2

   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= sander_natoms())
   a2 = 3*cv%i(2) - 2

   nfe_assert(cv%i(3) > 0 .and. cv%i(3) <= sander_natoms())
   a3 = 3*cv%i(3) - 2

   nfe_assert(a1 /= a2 .and. a1 /= a3 .and. a2 /= a3)

   NFE_MASTER_ONLY_BEGIN

   call angle_d(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2), d1, d2, d3)

   f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1
   f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2
   f(a3:a3 + 2) = f(a3:a3 + 2) + fcv*d3

   NFE_MASTER_ONLY_END

end subroutine colvar_force

!=============================================================================

subroutine colvar_bootstrap(cv, cvno)

   NFE_USE_AFAILED

   use nfe_colvar_type
   use nfe_colvar_utils

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno

   nfe_assert(cv%type == COLVAR_ANGLE)
   call check_i(cv%i, cvno, 'ANGLE', 3)

end subroutine colvar_bootstrap

!=============================================================================

subroutine print_details(cv, lun)

#ifndef NFE_DISABLE_ASSERT
   use nfe_utils
   use nfe_sander_proxy
#endif /* NFE_DISABLE_ASSERT */

   use nfe_colvar_type
   use nfe_colvar_utils

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_ANGLE)
   nfe_assert(associated(cv%i))

   call print_i(cv%i, lun)

end subroutine print_details

!=============================================================================

end module nfe_cv_ANGLE
