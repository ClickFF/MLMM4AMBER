! <compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

! fpan3-at-ncsu-dot-edu, 03/30/2016
!
! DF_COM_DISTANCE -- difference of distances between centers of mass of two atom groups
!
! input: i = (a1, ..., aN, 0, b1, ..., bM, 0, c1, ..., ck, d1, ..., dL)
! last zero is optional; a?/b?/c?/d? -- indices of the participating atoms

module nfe_cv_DF_COM_DISTANCE

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

   NFE_USE_AFAILED

   use nfe_constants, only : ZERO
   use nfe_colvar_type
   use nfe_colvar_math
   use nfe_colvar_utils

   implicit none

   NFE_REAL :: value

   type(colvar_t), intent(in) :: cv
   NFE_REAL, intent(in) :: x(*)

#  include "nfe-mpi.h"

   integer :: n

   NFE_REAL :: cm1(3), cm2(3), cm3(3), cm4(4)

   nfe_assert(cv%type == COLVAR_DF_COM_DISTANCE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.6)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

#ifdef MPI
   if (sanderrank.eq.0) then
#endif /* MPI */
   n = 1
   call group_com(cv, x, n, cm1)
   nfe_assert(n.le.size(cv%i))
   call group_com(cv, x, n, cm2)
   nfe_assert(n.le.size(cv%i))
   call group_com(cv, x, n, cm3)
   nfe_assert(n.le.size(cv%i))   
   call group_com(cv, x, n, cm4)   

   value = distance(cm1, cm2) - distance(cm3, cm4)
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function colvar_value

!=============================================================================

subroutine colvar_force(cv, x, fcv, f)

   NFE_USE_AFAILED

   use nfe_constants, only : ERR_UNIT
   use nfe_colvar_type
   use nfe_colvar_math
   use nfe_colvar_utils
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(in) :: cv
   NFE_REAL, intent(in) :: x(*), fcv
   NFE_REAL, intent(inout) :: f(*)

#  include "nfe-mpi.h"

   NFE_REAL :: d1(3), d2(3), d3(3), d4(3)
   NFE_REAL :: cm1(3), cm2(3), cm3(3), cm4(3)

   integer :: n

   nfe_assert(cv%type == COLVAR_DF_COM_DISTANCE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.6)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

   NFE_MASTER_ONLY_BEGIN

   n = 1
   call group_com(cv, x, n, cm1)
   call group_com(cv, x, n, cm2)
   call group_com(cv, x, n, cm3)
   call group_com(cv, x, n, cm4)   

   if ((distance(cm1, cm2).lt.1.0d-8).or.(distance(cm3, cm4).lt.1.0d-8)) then
      write (unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
         'DF_COM_DISTANCE : the distance got too small'
      call terminate()
   end if ! distance too small

   call distance_d(cm1, cm2, d1, d2)
   call distance_d(cm3, cm4, d3, d4)
   d3 = -d3
   d4 = -d4   

   d1 = fcv*d1
   d2 = fcv*d2
   d3 = fcv*d3
   d4 = fcv*d4   

   n = 1
   call group_com_d(cv, f, d1, n)
   call group_com_d(cv, f, d2, n)
   call group_com_d(cv, f, d3, n)
   call group_com_d(cv, f, d4, n)   

   NFE_MASTER_ONLY_END

end subroutine colvar_force

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   use nfe_utils
   use nfe_constants
   use nfe_colvar_type
   use nfe_colvar_utils
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   NFE_REAL,      intent(in)    :: amass(*)

#  include "nfe-mpi.h"

   integer :: error, ngroups

   nfe_assert(cv%type == COLVAR_DF_COM_DISTANCE)

   call com_check_i(cv%i, cvno, 'DF_COM_DISTANCE', ngroups)
   if (ngroups.ne.4) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (DF_COM_DISTANCE) : number of atom groups is not four'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! ngroups.ne.4

   if (associated(cv%r)) &
      deallocate(cv%r)

   allocate(cv%r(size(cv%i)), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   cv%r = ZERO
   do error = 1, size(cv%i)
      if (cv%i(error).gt.0) &
         cv%r(error) = amass(cv%i(error))
   end do

   call com_init_weights(cv%r)

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
   nfe_assert(cv%type == COLVAR_DF_COM_DISTANCE)
   nfe_assert(associated(cv%i))

   call com_print_i(cv%i, lun)

end subroutine print_details

!=============================================================================

end module nfe_cv_DF_COM_DISTANCE
