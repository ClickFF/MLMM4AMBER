! <compile=optimized>
#include "nfe-utils.h"
#include "nfe-config.h"

!
! cv%i = (i1, ..., iN) -- list of participating atoms
!

module nfe_cv_R_OF_GYRATION

!=============================================================================

implicit none

private

!=============================================================================

public :: colvar_value
public :: colvar_force

public :: colvar_bootstrap
public :: print_details

public :: colvar_cleanup

type, private :: priv_t
   NFE_REAL :: Rg, cm(3) ! value & center of mass
   NFE_REAL, pointer :: weights(:) ! m_i/total_mass
#  include "nfe-cv-priv.type"
end type priv_t

#include "nfe-cv-priv.decl"

!=============================================================================

contains

!=============================================================================

#include "nfe-cv-priv.impl"

!=============================================================================

function colvar_value(cv, x) result(value)

   use nfe_utils
   use nfe_constants
   use nfe_colvar_type

   implicit none

   NFE_REAL :: value

   type(colvar_t), intent(inout) :: cv

   NFE_REAL, intent(in) :: x(*)

   integer :: natoms, a, a3, i
   NFE_REAL :: tmp

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_R_OF_GYRATION)
   nfe_assert(associated(cv%i))

   natoms = size(cv%i)
   nfe_assert(natoms > 2)

   priv => get_priv(cv)

   priv%Rg = ZERO
   priv%cm(1:3) = ZERO

   do a = 1, natoms
      a3 = 3*(cv%i(a) - 1)
      tmp = ZERO
      do i = 1, 3
         tmp = tmp + x(a3 + i)**2
         priv%cm(i) = priv%cm(i) + priv%weights(a)*x(a3 + i)
      end do
      priv%Rg = priv%Rg + priv%weights(a)*tmp
   end do

   priv%Rg = sqrt(priv%Rg - priv%cm(1)**2 - priv%cm(2)**2 - priv%cm(3)**2)
   value = priv%Rg

end function colvar_value

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine colvar_force(cv, x, fcv, f)

   use nfe_utils
   use nfe_constants
   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   NFE_REAL, intent(in) :: x(*), fcv
   NFE_REAL, intent(inout) :: f(*)

#ifdef MPI
#  include "nfe-mpi.h"
   integer :: a_first, a_last
#endif

   integer :: natoms, a, a3
   NFE_REAL :: tmp

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_R_OF_GYRATION)
   nfe_assert(associated(cv%i))

   natoms = size(cv%i)
   nfe_assert(natoms > 2)

   priv => get_priv(cv)
   nfe_assert(priv%Rg > ZERO)

   tmp = fcv/priv%Rg

#ifdef MPI
   a = natoms/sandersize
   if (a.gt.0) then
      if (sanderrank.ne.(sandersize - 1)) then
         a_first = 1 + sanderrank*a
         a_last = (sanderrank + 1)*a
      else
         a_first = 1 + sanderrank*a
         a_last = natoms
      end if
   else
      if (sanderrank.eq.0) then
         a_first = 1
         a_last = natoms
      else
         a_first = 1
         a_last = 0
      end if
   end if
   do a = a_first, a_last
#else
   do a = 1, natoms
#endif /* MPI */
      a3 = 3*cv%i(a)

      f(a3 - 2:a3) = f(a3 - 2:a3) &
         + tmp*priv%weights(a)*(x(a3 - 2:a3) - priv%cm(1:3))
   end do

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

   integer   :: natoms, a, error
   NFE_REAL :: total_mass

   type(priv_t), pointer :: priv

#  include "nfe-mpi.h"

   nfe_assert(cv%type == COLVAR_R_OF_GYRATION)

   natoms = size(cv%i)

   call check_i(cv%i, cvno, 'R_OF_GYRATION')
   if (.not. natoms > 2) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (R_OF_GYRATION) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   total_mass = ZERO

   do a = 1, natoms
      total_mass = total_mass + amass(cv%i(a))
   end do

   nfe_assert(total_mass > ZERO)

   priv => new_priv(cv)

   allocate(priv%weights(natoms), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   do a = 1, natoms
      priv%weights(a) = amass(cv%i(a))/total_mass
   end do

end subroutine colvar_bootstrap

!=============================================================================

subroutine print_details(cv, lun)

   NFE_USE_AFAILED

   use nfe_colvar_type
   use nfe_colvar_utils

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(cv%type == COLVAR_R_OF_GYRATION)
   nfe_assert(associated(cv%i))

   call print_i(cv%i, lun)

end subroutine print_details

!=============================================================================

subroutine colvar_cleanup(cv)

   NFE_USE_AFAILED

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type.eq.COLVAR_R_OF_GYRATION)

   priv => get_priv(cv)
   nfe_assert(associated(priv))

   deallocate(priv%weights)
   call del_priv(cv)

end subroutine colvar_cleanup

!=============================================================================

end module nfe_cv_R_OF_GYRATION
