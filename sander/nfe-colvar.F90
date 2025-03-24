! <compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

module nfe_colvar

implicit none

private

!=============================================================================

public :: colvar_value
public :: colvar_force

public :: colvar_difference
public :: colvar_interpolate

public :: colvar_is_periodic

public :: colvar_has_min
public :: colvar_min

public :: colvar_has_max
public :: colvar_max

public :: colvar_print
public :: colvar_cleanup

!public :: colvar_mdread  old subroutine
public :: colvar_bootstrap
public :: colvar_nlread

public :: colvar_is_quaternion
public :: colvar_has_axis
public :: colvar_has_refcrd

!=============================================================================

contains

!=============================================================================

! the value is needed only on master
function colvar_value(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_colvar_type

   use nfe_cv_ANGLE,              only : v_ANGLE            => colvar_value
   use nfe_cv_TORSION,            only : v_TORSION          => colvar_value
   use nfe_cv_DISTANCE,           only : v_DISTANCE         => colvar_value
   use nfe_cv_MULTI_RMSD,         only : v_MULTI_RMSD       => colvar_value
   use nfe_cv_R_OF_GYRATION,      only : v_R_OF_GYRATION    => colvar_value
   use nfe_cv_HANDEDNESS,         only : v_HANDEDNESS       => colvar_value
   use nfe_cv_N_OF_BONDS,         only : v_N_OF_BONDS       => colvar_value
   use nfe_cv_N_OF_STRUCTURES,    only : v_N_OF_STRUCTURES  => colvar_value
   use nfe_cv_LCOD,               only : v_LCOD             => colvar_value
   use nfe_cv_COS_OF_DIHEDRAL,    only : v_COS_OF_DIHEDRAL  => colvar_value
   use nfe_cv_COM_ANGLE,          only : v_COM_ANGLE        => colvar_value
   use nfe_cv_COM_TORSION,        only : v_COM_TORSION      => colvar_value
   use nfe_cv_COM_DISTANCE,       only : v_COM_DISTANCE     => colvar_value
   use nfe_cv_PCA,                only : v_PCA              => colvar_value
   use nfe_cv_SIN_OF_DIHEDRAL,    only : v_SIN_OF_DIHEDRAL  => colvar_value
   use nfe_cv_PAIR_DIHEDRAL,      only : v_PAIR_DIHEDRAL    => colvar_value
   use nfe_cv_PATTERN_DIHEDRAL,   only : v_PATTERN_DIHEDRAL => colvar_value
   use nfe_cv_DF_COM_DISTANCE,    only : v_DF_COM_DISTANCE  => colvar_value
   use nfe_cv_ORIENTATION_ANGLE,  only : v_ORIENTATION_ANGLE=> colvar_value
   use nfe_cv_ORIENTATION_PROJ,   only : v_ORIENTATION_PROJ => colvar_value
   use nfe_cv_SPINANGLE,          only : v_SPINANGLE        => colvar_value
   use nfe_cv_TILT,               only : v_TILT             => colvar_value
   use nfe_cv_QUATERNION1,        only : v_QUATERNION1      => colvar_value
   use nfe_cv_QUATERNION2,        only : v_QUATERNION2      => colvar_value
   use nfe_cv_QUATERNION3,        only : v_QUATERNION3      => colvar_value
   use nfe_cv_QUATERNION0,        only : v_QUATERNION0      => colvar_value

   implicit none

   NFE_REAL :: value

   type(colvar_t) :: cv ! mutable
   NFE_REAL, intent(in) :: x(*)

   select case(cv%type)
      case(COLVAR_ANGLE)
         value = v_ANGLE(cv, x)
      case(COLVAR_TORSION)
         value = v_TORSION(cv, x)
      case(COLVAR_DISTANCE)
         value = v_DISTANCE(cv, x)
      case(COLVAR_MULTI_RMSD)
         value = v_MULTI_RMSD(cv, x)
      case(COLVAR_R_OF_GYRATION)
         value = v_R_OF_GYRATION(cv, x)
      case(COLVAR_HANDEDNESS)
         value = v_HANDEDNESS(cv, x)
      case(COLVAR_N_OF_BONDS)
         value = v_N_OF_BONDS(cv, x)
      case(COLVAR_N_OF_STRUCTURES)
         value = v_N_OF_STRUCTURES(cv, x)
      case(COLVAR_LCOD)
         value = v_LCOD(cv, x)
      case(COLVAR_COS_OF_DIHEDRAL)
         value = v_COS_OF_DIHEDRAL(cv, x)
      case(COLVAR_COM_ANGLE)
         value = v_COM_ANGLE(cv, x)
      case(COLVAR_COM_TORSION)
         value = v_COM_TORSION(cv, x)
      case(COLVAR_COM_DISTANCE)
         value = v_COM_DISTANCE(cv, x)   
      case(COLVAR_PCA)
         value = v_PCA(cv,x)
      case(COLVAR_SIN_OF_DIHEDRAL)
         value = v_SIN_OF_DIHEDRAL(cv,x) 
      case(COLVAR_PAIR_DIHEDRAL)
         value = v_PAIR_DIHEDRAL(cv,x) 
      case(COLVAR_PATTERN_DIHEDRAL)
         value = v_PATTERN_DIHEDRAL(cv,x)
      case(COLVAR_DF_COM_DISTANCE)
         value = v_DF_COM_DISTANCE(cv, x)
      case(COLVAR_ORIENTATION_ANGLE)
         value = v_ORIENTATION_ANGLE(cv,x)
      case(COLVAR_ORIENTATION_PROJ)
         value = v_ORIENTATION_PROJ(cv,x)
      case(COLVAR_SPINANGLE)
         value = v_SPINANGLE(cv,x)
      case(COLVAR_TILT)
         value = v_TILT(cv,x)
      case(COLVAR_QUATERNION0)
         value = v_QUATERNION0(cv,x)
      case(COLVAR_QUATERNION1)
         value = v_QUATERNION1(cv,x)
      case(COLVAR_QUATERNION2)
         value = v_QUATERNION2(cv,x)
      case(COLVAR_QUATERNION3)
         value = v_QUATERNION3(cv,x)
      case default
         nfe_assert_not_reached()
         value = NFE_TO_REAL(0)
   end select

end function colvar_value

!=============================================================================

subroutine colvar_force(cv, x, fcv, f)

   NFE_USE_AFAILED

   use nfe_colvar_type

   use nfe_cv_ANGLE,              only : f_ANGLE            => colvar_force
   use nfe_cv_TORSION,            only : f_TORSION          => colvar_force
   use nfe_cv_DISTANCE,           only : f_DISTANCE         => colvar_force
   use nfe_cv_MULTI_RMSD,         only : f_MULTI_RMSD       => colvar_force
   use nfe_cv_R_OF_GYRATION,      only : f_R_OF_GYRATION    => colvar_force
   use nfe_cv_HANDEDNESS,         only : f_HANDEDNESS       => colvar_force
   use nfe_cv_N_OF_BONDS,         only : f_N_OF_BONDS       => colvar_force
   use nfe_cv_N_OF_STRUCTURES,    only : f_N_OF_STRUCTURES  => colvar_force
   use nfe_cv_LCOD,               only : f_LCOD             => colvar_force
   use nfe_cv_COS_OF_DIHEDRAL,    only : f_COS_OF_DIHEDRAL  => colvar_force
   use nfe_cv_COM_ANGLE,          only : f_COM_ANGLE        => colvar_force
   use nfe_cv_COM_TORSION,        only : f_COM_TORSION      => colvar_force
   use nfe_cv_COM_DISTANCE,       only : f_COM_DISTANCE     => colvar_force
   use nfe_cv_PCA,                only : f_PCA              => colvar_force
   use nfe_cv_SIN_OF_DIHEDRAL,    only : f_SIN_OF_DIHEDRAL  => colvar_force
   use nfe_cv_PAIR_DIHEDRAL,      only : f_PAIR_DIHEDRAL    => colvar_force
   use nfe_cv_PATTERN_DIHEDRAL,   only : f_PATTERN_DIHEDRAL => colvar_force
   use nfe_cv_DF_COM_DISTANCE,    only : f_DF_COM_DISTANCE  => colvar_force
   use nfe_cv_ORIENTATION_ANGLE, only : f_ORIENTATION_ANGLE=> colvar_force
   use nfe_cv_ORIENTATION_PROJ,  only : f_ORIENTATION_PROJ => colvar_force
   use nfe_cv_SPINANGLE,         only : f_SPINANGLE        => colvar_force
   use nfe_cv_TILT,              only : f_TILT             => colvar_force
   use nfe_cv_QUATERNION0,       only : f_QUATERNION0      => colvar_force
   use nfe_cv_QUATERNION1,       only : f_QUATERNION1      => colvar_force
   use nfe_cv_QUATERNION2,       only : f_QUATERNION2      => colvar_force
   use nfe_cv_QUATERNION3,       only : f_QUATERNION3      => colvar_force


   implicit none

   type(colvar_t) :: cv ! mutable

   NFE_REAL, intent(in) :: x(*), fcv
   NFE_REAL, intent(inout) :: f(*)

   select case(cv%type)
      case(COLVAR_ANGLE)
         call f_ANGLE(cv, x, fcv, f)
      case(COLVAR_TORSION)
         call f_TORSION(cv, x, fcv, f)
      case(COLVAR_DISTANCE)
         call f_DISTANCE(cv, x, fcv, f)
      case(COLVAR_MULTI_RMSD)
         call f_MULTI_RMSD(cv, fcv, f)
      case(COLVAR_R_OF_GYRATION)
         call f_R_OF_GYRATION(cv, x, fcv, f)
      case(COLVAR_HANDEDNESS)
         call f_HANDEDNESS(cv, x, fcv, f)
      case(COLVAR_N_OF_BONDS)
         call f_N_OF_BONDS(cv, x, fcv, f)
      case(COLVAR_N_OF_STRUCTURES)
         call f_N_OF_STRUCTURES(cv, fcv, f)
      case(COLVAR_LCOD)
         call f_LCOD(cv, x, fcv, f)
      case(COLVAR_COS_OF_DIHEDRAL)
         call f_COS_OF_DIHEDRAL(cv, x, fcv, f)
      case(COLVAR_COM_ANGLE)
         call f_COM_ANGLE(cv, x, fcv, f)
      case(COLVAR_COM_TORSION)
         call f_COM_TORSION(cv, x, fcv, f)
      case(COLVAR_COM_DISTANCE)
         call f_COM_DISTANCE(cv, x, fcv, f)     
      case(COLVAR_PCA)
         call f_PCA(cv, fcv, f)
      case(COLVAR_SIN_OF_DIHEDRAL)
         call f_SIN_OF_DIHEDRAL(cv, x, fcv, f)
      case(COLVAR_PAIR_DIHEDRAL)
         call f_PAIR_DIHEDRAL(cv, x, fcv, f)
      case(COLVAR_PATTERN_DIHEDRAL)
         call f_PATTERN_DIHEDRAL(cv, x, fcv, f)      
      case(COLVAR_DF_COM_DISTANCE)
         call f_DF_COM_DISTANCE(cv, x, fcv, f)            
         case(COLVAR_ORIENTATION_ANGLE)
         call f_ORIENTATION_ANGLE(cv, fcv, f)
      case(COLVAR_ORIENTATION_PROJ)
         call f_ORIENTATION_PROJ(cv, fcv, f)
      case(COLVAR_SPINANGLE)
         call f_SPINANGLE(cv, fcv, f)
      case(COLVAR_TILT)
         call f_TILT(cv, fcv, f)
      case(COLVAR_QUATERNION0)
         call f_QUATERNION0(cv, fcv, f)
      case(COLVAR_QUATERNION1)
         call f_QUATERNION1(cv, fcv, f)
      case(COLVAR_QUATERNION2)
         call f_QUATERNION2(cv, fcv, f)
      case(COLVAR_QUATERNION3)
         call f_QUATERNION3(cv, fcv, f) 
      case default
         nfe_assert_not_reached()
   end select

end subroutine colvar_force

!=============================================================================

function colvar_difference(cv, v1, v2) result(diff)

   NFE_USE_AFAILED

   use nfe_colvar_type
   use nfe_constants, only : ZERO
   use constants, only : PI

   implicit none

   NFE_REAL :: diff
   type(colvar_t), intent(in) :: cv
   NFE_REAL, intent(in) :: v1, v2

   NFE_REAL :: t1, t2

   t1 = fix_value(v1)
   t2 = fix_value(v2)

   diff = t1 - t2
  
   if (cv%type.eq.COLVAR_TORSION.or.cv%type.eq.COLVAR_COM_TORSION) then
      nfe_assert(- PI.le.t1.and.t1.le.PI)
      nfe_assert(- PI.le.t2.and.t2.le.PI)
      if (diff.gt.PI) then
         diff = diff - PI - PI
      else if (diff.lt.-PI) then
         diff = diff + PI + PI
      end if
   end if

   if (cv%type.eq.COLVAR_SPINANGLE) then
      !nfe_assert(-180.0 < t1 < 180.0)
      !nfe_assert(-180.0 < t2 < 180.0)
      !print*, 't1', t1
      !print*, 't2', t2 
      if (diff > 180.0) then
         diff = diff - 360.0
      else if (diff < -180.0) then
         diff = diff + 360.0
      end if
   end if


contains

function fix_value(v) result(t)

   implicit none

   NFE_REAL :: t
   NFE_REAL, intent(in) :: v

   if (cv%type.eq.COLVAR_ANGLE.or.cv%type.eq.COLVAR_COM_ANGLE) then
      if (ZERO.le.v.and.v.le.PI) then
         t = v
      else
         t = acos(cos(v))
      end if
   else if (cv%type.eq.COLVAR_TORSION.or.cv%type.eq.COLVAR_COM_TORSION) then
      if (- PI.le.v.and.v.le.PI) then
        t = v
      else
         t = atan2(sin(v), cos(v))
      endif
   else
      t = v
   !print*, 'vvv', t  
   end if

end function fix_value

end function colvar_difference

!=============================================================================

function colvar_interpolate(cv, a1, v1, a2, v2) result(interp)

   NFE_USE_AFAILED

   use nfe_colvar_type

   implicit none

   NFE_REAL :: interp
   type(colvar_t), intent(in) :: cv
   NFE_REAL, intent(in) :: a1, v1, a2, v2

   NFE_REAL :: ts, tc

   if (cv%type.eq.COLVAR_TORSION.or.cv%type.eq.COLVAR_COM_TORSION) then
      ts = a1*sin(v1)+a2*sin(v2)
      tc = a1*cos(v1)+a2*cos(v2)
      interp = atan2(ts,tc)
   else
      interp = a1*v1+a2*v2
   end if

end function colvar_interpolate

!=============================================================================


logical function colvar_is_periodic(cv)

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   if (cv%type.eq.COLVAR_TORSION.or.cv%type.eq.COLVAR_COM_TORSION) then
      colvar_is_periodic = .true.
   else
      colvar_is_periodic = .false.
   end if
 
end function colvar_is_periodic

!=============================================================================


logical function colvar_has_axis(cv)

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   if (cv%type.eq.COLVAR_TILT.or.cv%type.eq.COLVAR_SPINANGLE) then
      colvar_has_axis = .true.
   else
      colvar_has_axis = .false.
   end if

end function colvar_has_axis

!=============================================================================

logical function colvar_has_refcrd(cv)

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   if (cv%type.eq.COLVAR_ORIENTATION_ANGLE.or.cv%type.eq.COLVAR_ORIENTATION_PROJ &
       .or.cv%type.eq.COLVAR_TILT.or.cv%type.eq.COLVAR_SPINANGLE &
       .or.cv%type.eq.COLVAR_QUATERNION0.or.cv%type.eq.COLVAR_QUATERNION1 &
       .or.cv%type.eq.COLVAR_QUATERNION2.or.cv%type.eq.COLVAR_QUATERNION3) then
      colvar_has_refcrd = .true.
   else
      colvar_has_refcrd = .false.
   end if

end function colvar_has_refcrd

!=============================================================================

logical function colvar_is_quaternion(cv)

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   if (cv%type.eq.COLVAR_QUATERNION0.or.cv%type.eq.COLVAR_QUATERNION1 &
      .or.cv%type.eq.COLVAR_QUATERNION2.or.cv%type.eq.COLVAR_QUATERNION3) then
       colvar_is_quaternion = .true.
   else
      colvar_is_quaternion = .false.
   end if

end function colvar_is_quaternion

!=============================================================================

logical function colvar_has_min(cv)

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   select case(cv%type)
      case(COLVAR_ORIENTATION_ANGLE)
         colvar_has_min = .true.
      case(COLVAR_ORIENTATION_PROJ)
         colvar_has_min = .true.
      case(COLVAR_SPINANGLE)
         colvar_has_min = .true.
      case(COLVAR_TILT)
         colvar_has_min = .true.
      case(COLVAR_ANGLE:COLVAR_R_OF_GYRATION)
         colvar_has_min = .true.
      case(COLVAR_N_OF_BONDS:COLVAR_N_OF_STRUCTURES)
         colvar_has_min = .true.
      case(COLVAR_COM_ANGLE:COLVAR_COM_DISTANCE)
         colvar_has_min = .true.
      case default
         colvar_has_min = .false.
   end select

end function colvar_has_min

!=============================================================================

NFE_REAL function colvar_min(cv)

   NFE_USE_AFAILED

   use constants, only : PI
   use nfe_constants, only : zero
   use nfe_constants, only : ONE
   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   nfe_assert(colvar_has_min(cv))

   select case(cv%type)
      case(COLVAR_ORIENTATION_ANGLE)
         colvar_min = ZERO
      case(COLVAR_ORIENTATION_PROJ)
         colvar_min = -ONE
      case(COLVAR_SPINANGLE)
         colvar_min = -PI
      case(COLVAR_TILT)
         colvar_min = -ONE
      case(COLVAR_ANGLE)
         colvar_min = ZERO
      case(COLVAR_TORSION)
         colvar_min = -PI
      case(COLVAR_DISTANCE:COLVAR_R_OF_GYRATION)
         colvar_min = ZERO
      case(COLVAR_N_OF_BONDS:COLVAR_N_OF_STRUCTURES)
         colvar_min = ZERO
      case(COLVAR_COM_ANGLE)
         colvar_min = ZERO
      case(COLVAR_COM_TORSION)
         colvar_min = -PI
      case(COLVAR_COM_DISTANCE)
         colvar_min = ZERO
      case default
         nfe_assert_not_reached()
         colvar_min = ZERO
   end select

end function colvar_min

!=============================================================================

logical function colvar_has_max(cv)

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   select case(cv%type)
      case(COLVAR_ORIENTATION_ANGLE)
         colvar_has_max = .true.
      case(COLVAR_ORIENTATION_PROJ)
         colvar_has_max = .true.
      case(COLVAR_SPINANGLE)
         colvar_has_max = .true.
      case(COLVAR_TILT)
         colvar_has_max = .true.
      case(COLVAR_ANGLE:COLVAR_TORSION)
         colvar_has_max = .true.
      case(COLVAR_COM_ANGLE:COLVAR_COM_TORSION)
         colvar_has_max = .true. 
      case default
         colvar_has_max = .false.
   end select

end function colvar_has_max

!=============================================================================

NFE_REAL function colvar_max(cv)

   NFE_USE_AFAILED

   use constants, only : PI
   use nfe_constants, only : ZERO
   use nfe_constants, only : ONE
   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   nfe_assert(colvar_has_max(cv))

   select case(cv%type)
      case(COLVAR_ANGLE:COLVAR_TORSION)
         colvar_max = PI
      case(COLVAR_ORIENTATION_ANGLE)
         colvar_max = PI
      case(COLVAR_ORIENTATION_PROJ)
         colvar_max = ONE
      case(COLVAR_SPINANGLE)
         colvar_max = PI
      case(COLVAR_TILT)
         colvar_max = ONE
      case(COLVAR_COM_ANGLE:COLVAR_COM_TORSION)
         colvar_max = PI
      case default
         nfe_assert_not_reached()
         colvar_max = ZERO
   end select

end function colvar_max

!=============================================================================

subroutine colvar_print(cv, lun)

   use nfe_utils
   use nfe_constants
   use nfe_colvar_type
   use nfe_sander_proxy

   use nfe_cv_ANGLE,              only : p_ANGLE            => print_details
   use nfe_cv_TORSION,            only : p_TORSION          => print_details
   use nfe_cv_DISTANCE,           only : p_DISTANCE         => print_details
   use nfe_cv_MULTI_RMSD,         only : p_MULTI_RMSD       => print_details
   use nfe_cv_R_OF_GYRATION,      only : p_R_OF_GYRATION    => print_details
   use nfe_cv_HANDEDNESS,         only : p_HANDEDNESS       => print_details
   use nfe_cv_N_OF_BONDS,         only : p_N_OF_BONDS       => print_details
   use nfe_cv_N_OF_STRUCTURES,    only : p_N_OF_STRUCTURES  => print_details
   use nfe_cv_LCOD,               only : p_LCOD             => print_details
   use nfe_cv_COS_OF_DIHEDRAL,    only : p_COS_OF_DIHEDRAL  => print_details
   use nfe_cv_COM_ANGLE,          only : p_COM_ANGLE        => print_details
   use nfe_cv_COM_TORSION,        only : p_COM_TORSION      => print_details
   use nfe_cv_COM_DISTANCE,       only : p_COM_DISTANCE     => print_details
   use nfe_cv_PCA,                only : p_PCA              => print_details
   use nfe_cv_SIN_OF_DIHEDRAL,    only : p_SIN_OF_DIHEDRAL  => print_details
   use nfe_cv_PAIR_DIHEDRAL,      only : p_PAIR_DIHEDRAL    => print_details
   use nfe_cv_PATTERN_DIHEDRAL,   only : p_PATTERN_DIHEDRAL => print_details
   use nfe_cv_DF_COM_DISTANCE,    only : p_DF_COM_DISTANCE  => print_details
   use nfe_cv_ORIENTATION_ANGLE, only : p_ORIENTATION_ANGLE=> print_details
   use nfe_cv_ORIENTATION_PROJ,  only : p_ORIENTATION_PROJ => print_details
   use nfe_cv_SPINANGLE,         only : p_SPINANGLE        => print_details
   use nfe_cv_TILT,              only : p_TILT             => print_details
   use nfe_cv_QUATERNION0,       only : p_QUATERNION0      => print_details
   use nfe_cv_QUATERNION1,       only : p_QUATERNION1      => print_details
   use nfe_cv_QUATERNION2,       only : p_QUATERNION2      => print_details
   use nfe_cv_QUATERNION3,       only : p_QUATERNION3      => print_details

   
   implicit none

   type(colvar_t), intent(in) :: cv
   integer,        intent(in) :: lun

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, '  type = '''

   select case(cv%type)
      case(COLVAR_ANGLE)
         write (unit = lun, fmt = '(a)') 'ANGLE'''
         call p_ANGLE(cv, lun)
      case(COLVAR_TORSION)
         write (unit = lun, fmt = '(a)') 'TORSION'''
         call p_TORSION(cv, lun)
      case(COLVAR_DISTANCE)
         write (unit = lun, fmt = '(a)') 'DISTANCE'''
         call p_DISTANCE(cv, lun)
      case(COLVAR_MULTI_RMSD)
         write (unit = lun, fmt = '(a)') 'MULTI_RMSD'''
         call p_MULTI_RMSD(cv, lun)
      case(COLVAR_R_OF_GYRATION)
         write (unit = lun, fmt = '(a)') 'R_OF_GYRATION'''
         call p_R_OF_GYRATION(cv, lun)
      case(COLVAR_HANDEDNESS)
         write (unit = lun, fmt = '(a)') 'HANDEDNESS'''
         call p_HANDEDNESS(cv, lun)
      case(COLVAR_N_OF_BONDS)
         write (unit = lun, fmt = '(a)') 'N_OF_BONDS'''
         call p_N_OF_BONDS(cv, lun)
      case(COLVAR_N_OF_STRUCTURES)
         write (unit = lun, fmt = '(a)') 'N_OF_STRUCTURES'''
         call p_N_OF_STRUCTURES(cv, lun)
      case(COLVAR_LCOD)
         write (unit = lun, fmt = '(a)') &
            'LCOD'' (Linear Combination Of Distances)'
         call p_LCOD(cv, lun)
      case(COLVAR_COS_OF_DIHEDRAL)
         write (unit = lun, fmt = '(a)') 'COS_OF_DIHEDRAL'''
         call p_COS_OF_DIHEDRAL(cv, lun)
      case(COLVAR_COM_ANGLE)
         write (unit = lun, fmt = '(a)') 'COM_ANGLE'''
         call p_COM_ANGLE(cv, lun)
      case(COLVAR_COM_TORSION)
         write (unit = lun, fmt = '(a)') 'COM_TORSION'''
         call p_COM_TORSION(cv, lun)
      case(COLVAR_COM_DISTANCE)
         write (unit = lun, fmt = '(a)') 'COM_DISTANCE'''
         call p_COM_DISTANCE(cv, lun)
      case(COLVAR_PCA) 
         write (unit = lun, fmt = '(a)') 'PRINCIPAL COMPONENT'''
         call p_PCA(cv, lun)
      case(COLVAR_SIN_OF_DIHEDRAL)
         write (unit = lun, fmt = '(a)') 'SIN_OF_DIHEDRAL'''
         call p_SIN_OF_DIHEDRAL(cv, lun)
      case(COLVAR_PAIR_DIHEDRAL)
         write (unit = lun, fmt = '(a)') 'PAIR_DIHEDRAL'''
         call p_PAIR_DIHEDRAL(cv, lun)
      case(COLVAR_PATTERN_DIHEDRAL)
         write (unit = lun, fmt = '(a)') 'PATTERN_DIHEDRAL'''
         call p_PATTERN_DIHEDRAL(cv, lun)
      case(COLVAR_DF_COM_DISTANCE)
         write (unit = lun, fmt = '(a)') 'DF_COM_DISTANCE'''
         call p_DF_COM_DISTANCE(cv, lun)
      case(COLVAR_ORIENTATION_ANGLE)
         write (unit = lun, fmt = '(a)') 'ORIENTATION_ANGLE'''
         call p_ORIENTATION_ANGLE(cv, lun)
      case(COLVAR_ORIENTATION_PROJ)
         write (unit = lun, fmt = '(a)') 'ORIENTATION_PROJ'''
         call p_ORIENTATION_PROJ(cv, lun)
      case(COLVAR_SPINANGLE)
         write (unit = lun, fmt = '(a)') 'SPINANGLE'''
         call p_SPINANGLE(cv, lun)
      case(COLVAR_TILT)
         write (unit = lun, fmt = '(a)') 'TILT'''
         call p_TILT(cv, lun)
      case(COLVAR_QUATERNION0)
         write (unit = lun, fmt = '(a)') 'QUATERNION0'''
         call p_QUATERNION0(cv, lun)
      case(COLVAR_QUATERNION1)
         write (unit = lun, fmt = '(a)') 'QUATERNION1'''
         call p_QUATERNION1(cv, lun)
      case(COLVAR_QUATERNION2)
         write (unit = lun, fmt = '(a)') 'QUATERNION2'''
         call p_QUATERNION2(cv, lun)
      case(COLVAR_QUATERNION3)
         write (unit = lun, fmt = '(a)') 'QUATERNION3'''
         call p_QUATERNION3(cv, lun)
      case default
         nfe_assert_not_reached()
         continue
   end select

end subroutine colvar_print

!=============================================================================

subroutine colvar_cleanup(cv)

   use nfe_colvar_type

   use nfe_cv_MULTI_RMSD,         only : c_MULTI_RMSD       => colvar_cleanup
   use nfe_cv_R_OF_GYRATION,      only : c_R_OF_GYRATION    => colvar_cleanup
   use nfe_cv_N_OF_STRUCTURES,    only : c_N_OF_STRUCTURES  => colvar_cleanup
   use nfe_cv_PCA,                only : c_PCA              => colvar_cleanup
   use nfe_cv_ORIENTATION_ANGLE, only : c_ORIENTATION_ANGLE=> colvar_cleanup
   use nfe_cv_ORIENTATION_PROJ,  only : c_ORIENTATION_PROJ => colvar_cleanup
   use nfe_cv_SPINANGLE,         only : c_SPINANGLE        => colvar_cleanup
   use nfe_cv_TILT,              only : c_TILT             => colvar_cleanup
   use nfe_cv_QUATERNION0,       only : c_QUATERNION0      => colvar_cleanup
   use nfe_cv_QUATERNION1,       only : c_QUATERNION1      => colvar_cleanup
   use nfe_cv_QUATERNION2,       only : c_QUATERNION2      => colvar_cleanup
   use nfe_cv_QUATERNION3,       only : c_QUATERNION3      => colvar_cleanup

   implicit none

   type(colvar_t), intent(inout) :: cv

   select case(cv%type)
      case(COLVAR_MULTI_RMSD)
         call c_MULTI_RMSD(cv)
      case(COLVAR_R_OF_GYRATION)
         call c_R_OF_GYRATION(cv)
      case(COLVAR_N_OF_STRUCTURES)
         call c_N_OF_STRUCTURES(cv)
      case (COLVAR_PCA) 
         call c_PCA(cv)
      case (COLVAR_ORIENTATION_ANGLE)
         call c_ORIENTATION_ANGLE(cv)
      case (COLVAR_ORIENTATION_PROJ)
         call c_ORIENTATION_PROJ(cv)
      case (COLVAR_SPINANGLE)
         call c_SPINANGLE(cv)
      case (COLVAR_TILT)
         call c_TILT(cv)
      case (COLVAR_QUATERNION0)
         call c_QUATERNION0(cv)
      case (COLVAR_QUATERNION1)
         call c_QUATERNION1(cv)
      case (COLVAR_QUATERNION2)
         call c_QUATERNION2(cv)
      case (COLVAR_QUATERNION3)
         call c_QUATERNION3(cv) 
      case default
         continue
   end select

   if (associated(cv%i)) &
      deallocate(cv%i)

   if (associated(cv%r)) &
      deallocate(cv%r)

   if (associated(cv%avgcrd)) &
      deallocate(cv%avgcrd)

   if (associated(cv%evec)) &
      deallocate(cv%evec)
 
   if (associated(cv%state_ref)) &
      deallocate(cv%state_ref)
      
   if (associated(cv%state_pca)) &
      deallocate(cv%state_pca) 
   
   if (associated(cv%ipca_to_i)) &
      deallocate(cv%ipca_to_i) 
   
   if (associated(cv%q_index)) &
      deallocate(cv%q_index) 
  
   if (associated(cv%axis)) &
      deallocate(cv%axis)

   cv%type = -1

end subroutine colvar_cleanup

!=============================================================================

!subroutine colvar_mdread(cv, node, cvno)

!   use nfe_utils
!   use nfe_value
!   use nfe_cftree
!   use nfe_constants
!   use nfe_colvar_type
!   use nfe_sander_proxy
!   use nfe_read_pca

!   implicit none

!   type(colvar_t), intent(inout) :: cv
!   type(node_t),   intent(in)    :: node
!   integer,        intent(in)    :: cvno

!   integer :: n, error
!   logical :: found2
!   integer :: nsolut
!   
!   type(value_node_t), pointer :: alist, aiter

!   character(len = STRING_LENGTH) :: type

!   ! declare three strings: ref_file avg_file evec_file 
!   character(len = STRING_LENGTH) :: ref_file
!   character(len = STRING_LENGTH) :: avg_file
!   character(len = STRING_LENGTH) :: evec_file
!   character(len = STRING_LENGTH) :: index_file 

!! DEBUG
!!  integer :: first, last 
!!  integer :: found , crdsize 

!   nfe_assert(is_master())

!   nfe_assert(.not. associated(cv%i))
!   nfe_assert(.not. associated(cv%r))
!   nfe_assert(.not. associated(cv%avgcrd)) 
!   nfe_assert(.not. associated(cv%evec))
!   nfe_assert(.not. associated(cv%state_ref))
!   nfe_assert(.not. associated(cv%state_pca)) 

!   nfe_assert(node_title(node) == 'variable')

!   !
!   ! type
!   !

!   if (.not.node_lookup_string(node, 'type', type)) then
!      write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//'/)') &
!            NFE_ERROR, 'type is not specified for CV #', cvno
!      call terminate()
!   end if

!   if (type == 'ANGLE') then
!      cv%type = COLVAR_ANGLE
!   else if (type == 'TORSION') then
!      cv%type = COLVAR_TORSION
!   else if (type == 'DISTANCE') then
!      cv%type = COLVAR_DISTANCE
!   else if (type == 'MULTI_RMSD') then
!      cv%type = COLVAR_MULTI_RMSD
!   else if (type == 'R_OF_GYRATION') then
!      cv%type = COLVAR_R_OF_GYRATION
!   else if (type == 'HANDEDNESS') then
!      cv%type = COLVAR_HANDEDNESS
!   else if (type == 'N_OF_BONDS') then
!      cv%type = COLVAR_N_OF_BONDS
!   else if (type == 'N_OF_STRUCTURES') then
!      cv%type = COLVAR_N_OF_STRUCTURES
!   else if (type == 'LCOD') then
!      cv%type = COLVAR_LCOD
!   else if (type == 'COS_OF_DIHEDRAL') then
!      cv%type = COLVAR_COS_OF_DIHEDRAL
!   else if (type == 'COM_ANGLE') then
!      cv%type = COLVAR_COM_ANGLE
!   else if (type == 'COM_TORSION') then
!      cv%type = COLVAR_COM_TORSION
!   else if (type == 'COM_DISTANCE') then
!      cv%type = COLVAR_COM_DISTANCE
!   else if (type == 'PCA') then
!      cv%type = COLVAR_PCA
!      nsolut = sander_nsolut()
!   else if (type == 'SIN_OF_DIHEDRAL') then
!      cv%type = COLVAR_SIN_OF_DIHEDRAL
!   else if (type == 'PAIR_DIHEDRAL') then
!      cv%type = COLVAR_PAIR_DIHEDRAL
!   else if (type == 'PATTERN_DIHEDRAL') then
!      cv%type = COLVAR_PATTERN_DIHEDRAL
!   else if (type == 'DF_COM_DISTANCE') then
!      cv%type = COLVAR_DF_COM_DISTANCE            
!   else
!      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,'//pfmt(cvno)//',a/)') &
!            NFE_ERROR, 'CV type ''', trim(type), &
!            ''' is not supported so far (CV #', cvno, ')'
!      call terminate()
!   end if

!   !
!   ! cv%i
!   !

!   if (node_lookup_list(node, 'i', alist)) then

!      n = 0
!      aiter => alist

!      do while (associated(aiter))
!         n = n + 1
!         if (.not. value_is_integer(aiter%value)) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
!               NFE_ERROR, 'CV #', cvno, ' : unexpected &
!               &(not an integer) element of ''i'' list'
!            call terminate()
!         end if
!         aiter => aiter%next
!      end do

!      if (n > 0) then
!         allocate(cv%i(n), stat = error)
!         if (error /= 0) &
!            NFE_OUT_OF_MEMORY

!         n = 0
!         aiter => alist
!         do while (associated(aiter))
!            n = n + 1
!            cv%i(n) = value_get_integer(aiter%value)
!            aiter => aiter%next
!         end do
!      end if

!   end if ! node_lookup_list(vnode, 'i', alist))

!   !
!   ! cv%r
!   !

!  if (type /= 'PCA') then 
!   if (node_lookup_list(node, 'r', alist)) then

!      n = 0
!      aiter => alist

!      do while (associated(aiter))
!         n = n + 1
!         if (.not. value_is_real(aiter%value)) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
!               NFE_ERROR, 'CV #', cvno, ' : unexpected &
!               &(not a real number) element of ''r'' list'
!            call terminate()
!         end if
!         aiter => aiter%next
!      end do

!      if (n > 0) then
!         allocate(cv%r(n), stat = error)
!         if (error /= 0) &
!            NFE_OUT_OF_MEMORY

!         n = 0
!         aiter => alist
!         do while (associated(aiter))
!            n = n + 1
!            cv%r(n) = value_get_real(aiter%value)
!            aiter => aiter%next
!         end do
!      end if

!   end if ! node_lookup_list(vnode, 'r', alist))
!  endif  
!  
!   ! read refcrd from "refcrd" 
!   ! store in cv% r 
!   if (type == 'PCA') then
!!       first = 1
!!       last  = sander_nsolut()
!!       crdsize = (last-first+1)*3
!       
!!       write(*,*) "sander_nsolut=", sander_nsolut();
!       
!!       nfe_assert(cv%i(1)*3 == crdsize)
!         
!       found2 = node_lookup_string(node, 'refcrd', ref_file)
!       if (.not. found2) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
!               NFE_ERROR, 'could not find ''reference coordinates'' for CV #', n
!            call terminate()
!       else
!       ! write (unit = OUT_UNIT, fmt = *) 'crdsize is', crdsize
!       ! allocate and read reference crd file into CV      
!           if(cv%i(1) > 0) then 
!              allocate(cv%r(cv%i(1)*3), stat = error)
!              if (error /= 0) &
!                 NFE_OUT_OF_MEMORY
!               call read_refcrd(cv, ref_file)
!               write (unit = OUT_UNIT, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
!               'ref_file = ', trim(ref_file), ' ('
!               write (unit = OUT_UNIT, fmt = '(a)') 'loaded)'
!              end if
!       endif

!      ! cv % evec 
!      found2 = node_lookup_string(node, 'evec', evec_file)
!      if (.not. found2 ) then
!           write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
!              NFE_ERROR, 'could not find ''reference coordinates'' for CV #', n
!            call terminate()
!      else
!      ! allocate and read evec crd into CV 
!!            if(crdsize > 0) then
!!              allocate(cv%evec(crdsize), stat = error)
!            if(cv%i(3)>0) then
!              allocate( cv%evec(cv%i(3)*3), stat = error)
!              ! write (unit = OUT_UNIT, fmt = *) 'crdsize is', crdsize
!              if (error /= 0) &
!                 NFE_OUT_OF_MEMORY
!!              call read_evec(cv, evec_file, first, last)
!              call read_evec(cv, evec_file)
!              write (unit = OUT_UNIT, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
!              'evec_file = ', trim(evec_file), ' ('
!              write (unit = OUT_UNIT, fmt = '(a)') 'loaded)'
!            endif
!      endif

!     
!      ! cv % avgcrd
!      found2 = node_lookup_string(node, 'avgcrd', avg_file)
!      if (.not. found2) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
!               NFE_ERROR, 'could not find ''average coordinates'' for CV #', n
!            call terminate()
!      else

!      ! allocate and read average crd into CV 
!!          if(crdsize > 0) then
!!              allocate(cv%avgcrd(crdsize), stat = error)
!           if(cv%i(3)>0) then
!              allocate( cv%avgcrd(cv%i(3)*3), stat = error)
!              if (error /= 0) &
!                 NFE_OUT_OF_MEMORY
!!             call read_avgcrd(cv, avg_file, first, last)
!              
!              call read_avgcrd(cv, avg_file)
!              write (unit = OUT_UNIT, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
!              'avg_file = ', trim(avg_file), ' ('
!              write (unit = OUT_UNIT, fmt = '(a)') 'loaded)'
!           endif
!       endif
!   
!   
!      ! cv % index
!      found2 = node_lookup_string(node, 'index', index_file)
!      if (.not. found2) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
!               NFE_ERROR, 'could not find ''index for ref and pca part'' for CV #', n
!            call terminate()
!      else
!      ! allocate and read average crd into CV 
!          if(cv%i(1) > 0 .and. cv%i(3) > 0) then
!              allocate(cv%state_ref(cv%i(1)), cv%state_pca(cv%i(1)), cv%ipca_to_i(cv%i(3)), stat = error)
!              if (error /= 0) &
!                 NFE_OUT_OF_MEMORY
!              call read_index(cv, in, only :: ERR_UNITdex_file)
!              write (unit = OUT_UNIT, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
!              'index_file = ', trim(index_file), ' ('
!              write (unit = OUT_UNIT, fmt = '(a)') 'loaded)'
!           endif
!       endif

!   
!  endif 
!     
! 
!end subroutine colvar_mdread

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   use nfe_utils
   use nfe_colvar_type
   use nfe_constants, only : ERR_UNIT
   use nfe_sander_proxy

   use nfe_cv_ANGLE,              only : b_ANGLE            => colvar_bootstrap
   use nfe_cv_TORSION,            only : b_TORSION          => colvar_bootstrap
   use nfe_cv_DISTANCE,           only : b_DISTANCE         => colvar_bootstrap
   use nfe_cv_MULTI_RMSD,         only : b_MULTI_RMSD       => colvar_bootstrap
   use nfe_cv_R_OF_GYRATION,      only : b_R_OF_GYRATION    => colvar_bootstrap
   use nfe_cv_HANDEDNESS,         only : b_HANDEDNESS       => colvar_bootstrap
   use nfe_cv_N_OF_BONDS,         only : b_N_OF_BONDS       => colvar_bootstrap
   use nfe_cv_N_OF_STRUCTURES,    only : b_N_OF_STRUCTURES  => colvar_bootstrap
   use nfe_cv_LCOD,               only : b_LCOD             => colvar_bootstrap
   use nfe_cv_COS_OF_DIHEDRAL,    only : b_COS_OF_DIHEDRAL  => colvar_bootstrap
   use nfe_cv_COM_ANGLE,          only : b_COM_ANGLE        => colvar_bootstrap
   use nfe_cv_COM_TORSION,        only : b_COM_TORSION      => colvar_bootstrap
   use nfe_cv_COM_DISTANCE,       only : b_COM_DISTANCE     => colvar_bootstrap
   use nfe_cv_PCA,                only : b_PCA              => colvar_bootstrap
   use nfe_cv_SIN_OF_DIHEDRAL,    only : b_SIN_OF_DIHEDRAL  => colvar_bootstrap
   use nfe_cv_PAIR_DIHEDRAL,      only : b_PAIR_DIHEDRAL    => colvar_bootstrap
   use nfe_cv_PATTERN_DIHEDRAL,   only : b_PATTERN_DIHEDRAL => colvar_bootstrap
   use nfe_cv_DF_COM_DISTANCE,    only : b_DF_COM_DISTANCE  => colvar_bootstrap
   use nfe_cv_ORIENTATION_ANGLE, only : b_ORIENTATION_ANGLE=> colvar_bootstrap
   use nfe_cv_ORIENTATION_PROJ,  only : b_ORIENTATION_PROJ => colvar_bootstrap
   use nfe_cv_SPINANGLE,         only : b_SPINANGLE        => colvar_bootstrap
   use nfe_cv_TILT,              only : b_TILT             => colvar_bootstrap
   use nfe_cv_QUATERNION0,       only : b_QUATERNION0      => colvar_bootstrap
   use nfe_cv_QUATERNION1,       only : b_QUATERNION1      => colvar_bootstrap
   use nfe_cv_QUATERNION2,       only : b_QUATERNION2      => colvar_bootstrap
   use nfe_cv_QUATERNION3,       only : b_QUATERNION3      => colvar_bootstrap



   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   NFE_REAL,      intent(in)    :: amass(*)

#ifdef MPI
   integer :: bcastdata(10), ierr
   ! original: integer :: bcastdata(3), ierr
#include "nfe-mpi.h"

   !
   ! bcast type/i/r first
   !

#ifndef NFE_DISABLE_ASSERT
   if (sanderrank == 0) then
      nfe_assert(cv%type > 0)
   else
      nfe_assert(.not. associated(cv%i))
      nfe_assert(.not. associated(cv%r))
      if (cv%type.eq.COLVAR_QUATERNION0.or.cv%type.eq.COLVAR_QUATERNION1 &
          .or.cv%type.eq.COLVAR_QUATERNION2.or.cv%type.eq.COLVAR_QUATERNION3) then
          nfe_assert(.not. associated(cv%q_index))
      endif 
      if (cv%type.eq.COLVAR_TILT.or.cv%type.eq.COLVAR_SPINANGLE) then
         nfe_assert(.not. associated(cv%axis))
      endif 
      if (cv%type == COLVAR_PCA) then 
        nfe_assert(.not. associated(cv%avgcrd))
        nfe_assert(.not. associated(cv%evec))
        nfe_assert(.not. associated(cv%state_ref))
        nfe_assert(.not. associated(cv%state_pca))
        nfe_assert(.not. associated(cv%ipca_to_i)) 
      endif
   end if
#endif /* NFE_DISABLE_ASSERT */

   if (sanderrank == 0) then
      bcastdata(1) = cv%type

      bcastdata(2) = 0
      if (associated(cv%i)) &
         bcastdata(2) = size(cv%i)

      bcastdata(3) = 0
      if (associated(cv%r)) &
         bcastdata(3) = size(cv%r)
     
      bcastdata(4) = 0 
      bcastdata(5) = 0
      bcastdata(6) = 0
      bcastdata(7) = 0
      bcastdata(8) = 0
      bcastdata(9) = 0               !q_index 
      bcastdata(10) = 0              !axis

     if (cv%type == COLVAR_PCA) then 
       
        if (associated(cv%avgcrd)) &
          bcastdata(4) = size(cv%avgcrd) 
        
        if (associated(cv%evec)) &
          bcastdata(5) = size(cv%evec)  
          
        if (associated(cv%state_ref)) &
          bcastdata(6) = size(cv%state_ref)
          
        if (associated(cv%state_pca)) &
          bcastdata(7) = size(cv%state_pca)
        
        if (associated(cv%ipca_to_i)) &
          bcastdata(8) = size(cv%ipca_to_i)  
          
     endif
     if (cv%type.eq.COLVAR_QUATERNION0.or.cv%type.eq.COLVAR_QUATERNION1 &
          .or.cv%type.eq.COLVAR_QUATERNION2.or.cv%type.eq.COLVAR_QUATERNION3) then
        if (associated(cv%q_index)) &
           bcastdata(9) = cv%q_index
     endif
     if (cv%type.eq.COLVAR_TILT.or.cv%type.eq.COLVAR_SPINANGLE) then
       if (associated(cv%axis)) &
         bcastdata(10) = size(cv%axis)       
     endif
   end if ! sanderrank == 0

   call mpi_bcast(bcastdata, size(bcastdata), MPI_INTEGER, 0, commsander, ierr)
   nfe_assert(ierr == 0)

   if (sanderrank /= 0) &
      cv%type = bcastdata(1)

   !
   ! cv%i
   !
   if (bcastdata(2) > 0) then
      if (.not. associated(cv%i)) then
         nfe_assert(sanderrank > 0)
         allocate(cv%i(bcastdata(2)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%i)

      call mpi_bcast(cv%i, bcastdata(2), MPI_INTEGER, 0, commsander, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%i)
   end if ! bcastdata(2) > 0

   !
   ! cv%r
   !
   if (bcastdata(3) > 0) then
      if (.not. associated(cv%r)) then
         nfe_assert(sanderrank > 0)
         allocate(cv%r(bcastdata(3)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%r)

      call mpi_bcast(cv%r, bcastdata(3), MPI_DOUBLE_PRECISION, &
                     0, commsander, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%r)
   end if ! bcastdata(3) > 0

   ! 
   ! cv % avgcrd    
   ! 
   if (bcastdata(4) > 0) then
      if (.not. associated(cv%avgcrd)) then
         nfe_assert(sanderrank > 0)
         allocate(cv%avgcrd(bcastdata(4)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%avgcrd)

      call mpi_bcast(cv%avgcrd, bcastdata(4), MPI_DOUBLE_PRECISION, &
                     0, commsander, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%avgcrd)
   end if ! bcastdata(4) > 0
  
   ! 
   ! cv % evec 
   ! 
   if (bcastdata(5) > 0) then
      if (.not. associated(cv%evec)) then
         nfe_assert(sanderrank > 0)
         allocate(cv%evec(bcastdata(5)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%evec)

      call mpi_bcast(cv%evec, bcastdata(5), MPI_DOUBLE_PRECISION, &
                     0, commsander, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%evec)
   end if ! bcastdata(5) > 0


   ! 
   ! cv % state_ref 
   ! 
   if (bcastdata(6) > 0) then
      if (.not. associated(cv%state_ref)) then
         nfe_assert(sanderrank > 0)
         allocate(cv%state_ref(bcastdata(6)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%state_ref)

      call mpi_bcast(cv%state_ref, bcastdata(6), MPI_INTEGER, &
                     0, commsander, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%state_ref)
   end if ! bcastdata(6) > 0


   ! 
   ! cv % state_pca
   ! 
   if (bcastdata(7) > 0) then
      if (.not. associated(cv%state_pca)) then
         nfe_assert(sanderrank > 0)
         allocate(cv%state_pca(bcastdata(7)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%state_pca)

      call mpi_bcast(cv%state_pca, bcastdata(7), MPI_INTEGER, &
                     0, commsander, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%state_pca)
   end if ! bcastdata(7) > 0
   
   ! 
   ! cv % ipca_to_i()
   ! 
   if (bcastdata(8) > 0) then
      if (.not. associated(cv%ipca_to_i)) then
         nfe_assert(sanderrank > 0)
         allocate(cv%ipca_to_i(bcastdata(8)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%state_pca)

      call mpi_bcast(cv%ipca_to_i, bcastdata(8), MPI_INTEGER, &
                     0, commsander, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%ipca_to_i)
   end if ! bcastdata(8) > 0
   
   !
   ! cv%q_index
   !
   if (bcastdata(9) > 0) then
     if (.not. associated(cv%q_index)) then
         nfe_assert(sanderrank > 0)
         allocate(cv%q_index, stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%q_index)
      !cv%q_index = bcastdata(9) ?
      call mpi_bcast(cv%q_index, 1, MPI_INTEGER, &
                     0, commsander, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%q_index)
   end if ! bcastdata(9) > 0

   !
   ! cv%axis
   !
   if (bcastdata(10) > 0) then
      if (.not. associated(cv%axis)) then
         nfe_assert(sanderrank > 0)
         allocate(cv%axis(bcastdata(10)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%axis)
      call mpi_bcast(cv%axis, bcastdata(10), MPI_DOUBLE_PRECISION, &
                     0, commsander, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%axis)
   end if ! bcastdata(10) > 0
#endif /* MPI */

   !
   ! dispatch according to the type
   !

   select case(cv%type)
      case(COLVAR_ANGLE)
         call b_ANGLE(cv, cvno)
      case(COLVAR_TORSION)
         call b_TORSION(cv, cvno)
      case(COLVAR_DISTANCE)
         call b_DISTANCE(cv, cvno)
      case(COLVAR_MULTI_RMSD)
         call b_MULTI_RMSD(cv, cvno, amass)
      case(COLVAR_R_OF_GYRATION)
         call b_R_OF_GYRATION(cv, cvno, amass)
      case(COLVAR_HANDEDNESS)
         call b_HANDEDNESS(cv, cvno)
      case(COLVAR_N_OF_BONDS)
         call b_N_OF_BONDS(cv, cvno)
      case(COLVAR_N_OF_STRUCTURES)
         call b_N_OF_STRUCTURES(cv, cvno, amass)
      case(COLVAR_LCOD)
         call b_LCOD(cv, cvno)
      case(COLVAR_COS_OF_DIHEDRAL)
         call b_COS_OF_DIHEDRAL(cv, cvno)
      case(COLVAR_COM_ANGLE)
         call b_COM_ANGLE(cv, cvno, amass)
      case(COLVAR_COM_TORSION)
         call b_COM_TORSION(cv, cvno, amass)
      case(COLVAR_COM_DISTANCE)
         call b_COM_DISTANCE(cv, cvno, amass)
      case(COLVAR_PCA)
         call b_PCA(cv, cvno, amass)
      case(COLVAR_SIN_OF_DIHEDRAL)
         call b_SIN_OF_DIHEDRAL(cv, cvno)
      case(COLVAR_PAIR_DIHEDRAL)
         call b_PAIR_DIHEDRAL(cv, cvno)
      case(COLVAR_PATTERN_DIHEDRAL)
         call b_PATTERN_DIHEDRAL(cv, cvno)
      case(COLVAR_DF_COM_DISTANCE)
         call b_DF_COM_DISTANCE(cv, cvno, amass)
      case(COLVAR_ORIENTATION_ANGLE)
         call b_ORIENTATION_ANGLE(cv, cvno, amass)
      case(COLVAR_ORIENTATION_PROJ)
         call b_ORIENTATION_PROJ(cv, cvno, amass)
      case(COLVAR_SPINANGLE)
         call b_SPINANGLE(cv, cvno, amass)
      case(COLVAR_TILT)
         call b_TILT(cv, cvno, amass)
      case(COLVAR_QUATERNION0)
         call b_QUATERNION0(cv, cvno, amass)
      case(COLVAR_QUATERNION1)
         call b_QUATERNION1(cv, cvno, amass)
      case(COLVAR_QUATERNION2)
         call b_QUATERNION2(cv, cvno, amass)
      case(COLVAR_QUATERNION3)
         call b_QUATERNION3(cv, cvno, amass)
      case default
         nfe_assert_not_reached()
         continue
   end select

end subroutine colvar_bootstrap

!=============================================================================

subroutine colvar_nlread(cv_unit,cv)

  use nfe_colvar_type
  use nfe_sander_proxy
  use nfe_constants
  use nfe_utils
  use sander_lib, only : upper

  implicit none

  integer, intent(in)           :: cv_unit
  type(colvar_t), intent(inout) :: cv

  integer                    :: ifind, nsolut, error, i 
  
  call nmlsrc('colvar', cv_unit, ifind)

  if (ifind .eq. 0) then
     write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
                 'cannot find collective variables info...'
     call terminate()
  end if
  
  ! initiate 
  cv_min = 0.0
  cv_max = 5.0
  resolution = 0.0
  cv_type = ''
  cv_i(:) = 0
  cv_r(:) = 0.0
  cv_ni = 0
  cv_nr = 0
  path(:) = 0.0
  harm(:) = 0.0
  path_mode = ' '
  harm_mode = ' '
  anchor_position(:) = 0.0
  anchor_strength(:) = 0.0
  q_index = 1
  axis = [0.0, 0.0, 1.0]
  refcrd_file = 'inpcrd'

  read(cv_unit,nml=colvar,err=666)
  
  if (cv_ni.le.0) then
     write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
           'number of elements in cv_i array must be positive!'
     call terminate()
  else
     allocate(cv%i(cv_ni), stat = error)
     if (error /= 0) &
     NFE_OUT_OF_MEMORY
     
     i = 1
     do while (i.le.cv_ni)
       cv%i(i) = cv_i(i)
       i = i + 1
     end do
  end if
  
  if (cv_nr.lt.0) then
     write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
           'number of elements in cv_r array cannot be negative!'
     call terminate()
  end if
  
  if (cv_nr.gt.0) then
     allocate(cv%r(cv_nr), stat = error)
     if (error /= 0) &
     NFE_OUT_OF_MEMORY
     
     i = 1
     do while (i.le.cv_nr)
       cv%r(i) = cv_r(i)
       i = i + 1
     end do
  end if
  
  call upper(cv_type)
  if (cv_type == 'ANGLE') then
      cv%type = COLVAR_ANGLE
  else if (cv_type == 'TORSION') then
      cv%type = COLVAR_TORSION
  else if (cv_type == 'DISTANCE') then
      cv%type = COLVAR_DISTANCE
  else if (cv_type == 'MULTI_RMSD') then
      cv%type = COLVAR_MULTI_RMSD
  else if (cv_type == 'R_OF_GYRATION') then
      cv%type = COLVAR_R_OF_GYRATION
  else if (cv_type == 'HANDEDNESS') then
      cv%type = COLVAR_HANDEDNESS
  else if (cv_type == 'N_OF_BONDS') then
      cv%type = COLVAR_N_OF_BONDS
  else if (cv_type == 'N_OF_STRUCTURES') then
      cv%type = COLVAR_N_OF_STRUCTURES
  else if (cv_type == 'LCOD') then
      cv%type = COLVAR_LCOD
  else if (cv_type == 'COS_OF_DIHEDRAL') then
      cv%type = COLVAR_COS_OF_DIHEDRAL
  else if (cv_type == 'COM_ANGLE') then
      cv%type = COLVAR_COM_ANGLE
  else if (cv_type == 'COM_TORSION') then
      cv%type = COLVAR_COM_TORSION
  else if (cv_type == 'COM_DISTANCE') then
      cv%type = COLVAR_COM_DISTANCE
  else if (cv_type == 'PCA') then
      cv%type = COLVAR_PCA
      nsolut = sander_nsolut()
  else if (cv_type == 'SIN_OF_DIHEDRAL') then
      cv%type = COLVAR_SIN_OF_DIHEDRAL
  else if (cv_type == 'PAIR_DIHEDRAL') then
      cv%type = COLVAR_PAIR_DIHEDRAL
  else if (cv_type == 'PATTERN_DIHEDRAL') then
      cv%type = COLVAR_PATTERN_DIHEDRAL
  else if (cv_type == 'DF_COM_DISTANCE') then
      cv%type = COLVAR_DF_COM_DISTANCE
  else if (cv_type == 'ORIENTATION_ANGLE') then
      cv%type = COLVAR_ORIENTATION_ANGLE
  else if (cv_type == 'ORIENTATION_PROJ') then
      cv%type = COLVAR_ORIENTATION_PROJ
  else if (cv_type == 'SPINANGLE') then
      cv%type = COLVAR_SPINANGLE
  else if (cv_type == 'TILT') then
      cv%type = COLVAR_TILT
  else if (cv_type == 'QUATERNION0') then
      cv%type = COLVAR_QUATERNION0
  else if (cv_type == 'QUATERNION1') then
      cv%type = COLVAR_QUATERNION1
  else if (cv_type == 'QUATERNION2') then
      cv%type = COLVAR_QUATERNION2
  else if (cv_type == 'QUATERNION3') then
      cv%type = COLVAR_QUATERNION3
  else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
            NFE_ERROR, 'CV type ''', trim(cv_type), &
            ''' is not supported so far '
      call terminate()
  end if
  
  return
  
  
  
  
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &colvar namelist!'
    call terminate()

end subroutine colvar_nlread

end module nfe_colvar
