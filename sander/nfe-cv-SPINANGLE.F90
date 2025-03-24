! <compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"
!
! input:
!
! cv%i = (a1, a2, ..., aN)
!
!     (list of participating atoms)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, ...)
!
!        (reference coordinates)
!
! axis = (e1, e2, e3)
!      (axis which defined by user)
! value = spinAngle[-180,180]: phi = 2.0 * atan (alpha)
! where alpha= q.e/q0
! in which q=(q1,q2,q3) is orientation vector and e=axis
!

module nfe_cv_SPINANGLE

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
public :: read_refcrd
public :: print_details

!=============================================================================
   NFE_REAL :: get_vec(3)
   NFE_REAL :: alpha, beta
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

type, private :: priv_t
   NFE_REAL :: value
   NFE_REAL :: total_mass
   integer :: n_atoms
   NFE_REAL, pointer :: mass(:) => null()
   NFE_REAL, pointer :: cm_crd(:) => null()
   NFE_REAL, pointer :: ref_crd(:) => null()
   NFE_REAL :: ref_nrm
   NFE_REAL :: quaternion(4,4)
   NFE_REAL :: lambda(4)

#  include "nfe-cv-priv.type"
end type priv_t

#include "nfe-cv-priv.decl"

!=============================================================================

contains

!=============================================================================

#include "nfe-cv-priv.impl"

!=============================================================================

function colvar_value(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_rmsd
   use nfe_utils
   use nfe_constants
   use nfe_colvar_type
   use nfe_colvar_math, only : norm3

   implicit none

   NFE_REAL :: value

   type(colvar_t), intent(inout) :: cv

   NFE_REAL, intent(in) :: x(*)
   NFE_REAL :: e(3)

#  include "nfe-mpi.h"

   type(priv_t), pointer :: priv

#ifdef MPI
   !integer :: error
#endif /* MPI */

   NFE_REAL :: cm(3), cm_nrm
   integer :: i, n
! gfortran 4.2 and lower does not support the volatile keyword, but this is
! necessary to work around an intel compiler optimization bug. So turn off
! the volatile keyword for qualifying compilers
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
   integer, volatile :: a, a3 !work around ifort optimization bug      !VOLATILE attribute and statement specifies
                                                                        !that the value of an object is unpredictable.
#else
   integer :: a, a3
#endif /* __GFORTRAN__ */

   nfe_assert(cv%type == COLVAR_SPINANGLE)

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

   priv%value = ZERO

   cm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*cv%i(a)
      cm = cm + priv%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/priv%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + priv%mass(n)*priv%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call orientation_q(priv%n_atoms, priv%mass, &
         priv%cm_crd, priv%ref_crd, priv%lambda, priv%quaternion)
   get_vec = [priv%quaternion(2,1), priv%quaternion(3,1), priv%quaternion(4,1)]
   if (norm3(cv%axis).ne.ZERO) then
       e = cv%axis/norm3(cv%axis)
   else
       e = [0.0, 0.0, 1.0]                                           !defualt value for the axis: [0,0, 0,0, 1.0]
   end if
   alpha = (dot_product(get_vec, e)) 
   priv%value = 2.0 * atan2 (alpha, priv%quaternion(1,1)) * 180.0/PI
   do while (priv%value > 180.0)                                     !check the value is between [-180,180]
         priv%value = 360.0 - priv%value
   enddo
   do while (priv%value < - 180.0)
         priv%value = 360.0 + priv%value
   enddo

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
   use nfe_colvar_math

   implicit none

   type(colvar_t), intent(in) :: cv

   NFE_REAL, intent(in) :: fcv
   NFE_REAL, intent(inout) :: f(*)
   NFE_REAL :: e(3)

#  include "nfe-mpi.h"

   integer :: n, a, a3, f3
   NFE_REAL :: dSx(4,4), dSy(4,4), dSz(4,4), QL(3)
   NFE_REAL :: dSxq1(4), dSyq1(4), dSzq1(4)
   NFE_REAL :: a1x, a1y, a1z, tmp
   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_SPINANGLE)
   nfe_assert(associated(cv%i))

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, commsander, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

   get_vec = [priv%quaternion(2,1), priv%quaternion(3,1), priv%quaternion(4,1)]
   if (norm3(cv%axis).ne.ZERO) then
       e = cv%axis/norm3(cv%axis)
   else
       e = [0.0, 0.0, 1.0]
   end if
   alpha = (dot_product(get_vec, e)) / priv%quaternion(1,1)
   beta = dot_product(get_vec, e)
   if ( priv%quaternion(1,1).ne.ZERO ) then                                     !factor out q0 to avoid being zero
      tmp = (2.0/priv%quaternion(1,1))*(1.0/(1.0 + alpha*alpha)) * 180.0/PI     ! +/- q11 !XX
    
                                                                                !computing QL
      QL(1) = (fcv*tmp*(alpha*priv%quaternion(1,2) - &
                 (e(1) * priv%quaternion(2,2) - &
                  e(2) * priv%quaternion(3,2) - &
                  e(3) * priv%quaternion(4,2))))/ &
                 max(priv%lambda(1)-priv%lambda(2),NFE_TO_REAL(0.00000001))

      QL(2) = (fcv*tmp*(alpha*priv%quaternion(1,3) - &
                 (e(1) * priv%quaternion(2,3) - &
                  e(2) * priv%quaternion(3,3) - &
                  e(3) * priv%quaternion(4,3))))/ &
                 max(priv%lambda(1)-priv%lambda(3),NFE_TO_REAL(0.00000001))

      QL(3) = (fcv*tmp*(alpha*priv%quaternion(1,4) - &
                 (e(1) * priv%quaternion(2,4) - &
                  e(2) * priv%quaternion(3,4) - &
                  e(3) * priv%quaternion(4,4))))/ &
                 max(priv%lambda(1)-priv%lambda(4),NFE_TO_REAL(0.00000001))
   else
      tmp = 2.0/dot_product(get_vec, e)                                         !1/1+x^2~1/x^2 and ignoring the second part!
      QL(1) = fcv*tmp*priv%quaternion(1,2) / &
              max(priv%lambda(1)-priv%lambda(2),NFE_TO_REAL(0.00000001))
      QL(2) = fcv*tmp*priv%quaternion(1,3) / &
              max(priv%lambda(1)-priv%lambda(3),NFE_TO_REAL(0.00000001))
      QL(3) = fcv*tmp*priv%quaternion(1,4) / &
              max(priv%lambda(1)-priv%lambda(4),NFE_TO_REAL(0.00000001))        !XX to do! what if dot_product(get_vec, e) = 0
   end if

      n = 1
      do a = 1, priv%n_atoms
         a3 = 3*n
         f3 = 3*cv%i(a)

         a1x = priv%ref_crd(a3 - 2)
         a1y = priv%ref_crd(a3 - 1)
         a1z = priv%ref_crd(a3)

         dSx(1,1) = a1x
         dSy(1,1) = a1y
         dSz(1,1) = a1z

         dSx(2,1) = 0.0
         dSy(2,1) = -a1z
         dSz(2,1) = a1y

         dSx(1,2) = 0.0
         dSy(1,2) = -a1z
         dSz(1,2) = a1y

         dSx(3,1) = a1z
         dSy(3,1) = 0.0
         dSz(3,1) = -a1x

         dSx(1,3) = a1z
         dSy(1,3) = 0.0
         dSz(1,3) = -a1x

         dSx(4,1) = -a1y
         dSy(4,1) = a1x
         dSz(4,1) = 0.0

         dSx(1,4) = -a1y
         dSy(1,4) = a1x
         dSz(1,4) = 0.0

         dSx(2,2) = a1x
         dSy(2,2) = -a1y
         dSz(2,2) = -a1z

         dSx(3,2) = a1y
         dSy(3,2) = a1x
         dSz(3,2) = 0.0

         dSx(2,3) = a1y
         dSy(2,3) = a1x
         dSz(2,3) = 0.0

         dSx(4,2) = a1z
         dSy(4,2) = 0.0
         dSz(4,2) = a1x

         dSx(2,4) = a1z
         dSy(2,4) = 0.0
         dSz(2,4) = a1x

         dSx(3,3) = -a1x
         dSy(3,3) = a1y
         dSz(3,3) = -a1z

         dSx(4,3) = 0.0
         dSy(4,3) = a1z
         dSz(4,3) = a1y

         dSx(3,4) = 0.0
         dSy(3,4) = a1z
         dSz(3,4) = a1y

         dSx(4,4) = -a1x
         dSy(4,4) = -a1y
         dSz(4,4) = a1z

         dSxq1 = matmul(dSx,priv%quaternion(:,1))
         dSyq1 = matmul(dSy,priv%quaternion(:,1))
         dSzq1 = matmul(dSz,priv%quaternion(:,1))

          f(f3 - 2) = f(f3 - 2) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSxq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSxq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSxq1,priv%quaternion(:,4)) )

          f(f3 - 1) = f(f3 - 1) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSyq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSyq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSyq1,priv%quaternion(:,4)) )

          f(f3) = f(f3) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSzq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSzq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSzq1,priv%quaternion(:,4)) )
         n = n + 1
      end do

end subroutine colvar_force

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   use nfe_utils
   use nfe_constants
   use nfe_colvar_type
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(inout) :: cv
   NFE_REAL, pointer :: coor(:) => null() 
   integer,        intent(in)    :: cvno
   NFE_REAL,      intent(in)    :: amass(*)

   integer :: i, n_atoms, error
   integer :: a, b, a1
   NFE_REAL :: cm(3)

   type(priv_t), pointer :: priv

#  include "nfe-mpi.h"

   nfe_assert(cv%type == COLVAR_SPINANGLE)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (cv_nr.gt.0) then
    if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (ORIENTATION) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
    end if
   end if 

   ! count the number of atoms
   n_atoms = 0

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (cv_nr.eq.0) then
      if (associated(coor)) &
          deallocate(coor)
          allocate(coor(3*sander_natoms()), stat = error)     
          if (error.ne.0) &
           NFE_OUT_OF_MEMORY
           call read_refcrd(coor, refcrd_file)
           allocate(cv%r(3*size(cv%i)), stat = error)
             if (error.ne.0) &
              NFE_OUT_OF_MEMORY
              do a = 1, n_atoms
               a1=cv%i(a)
               cv%r(1 + 3*(a - 1): 3*a) = coor(1 + 3*(a1 - 1): 3*a1)
              end do
          deallocate(coor)
    end if

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r) /= 3*n_atoms

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! basic checks
   if (n_atoms < 3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (1)//',a,'//pfmt(n_atoms)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (ORIENTATION) : too few integers in reference (', 1, ':', n_atoms, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = 1, n_atoms
      if (cv%i(a) < 1 .or. cv%i(a) > sander_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(sander_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', sander_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = 1, n_atoms
      do b = a + 1, n_atoms
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (ORIENTATION) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! cv%i(a) == cv%i(b)
      end do
   end do

   ! allocate/setup

   priv%n_atoms = n_atoms

   allocate(priv%mass(n_atoms), priv%cm_crd(3*n_atoms), &
            priv%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   priv%total_mass = ZERO

   do a = 1, n_atoms
      priv%mass(a) = amass(cv%i(a))
      priv%total_mass = priv%total_mass + priv%mass(a)
      cm = cm + priv%mass(a)*cv%r(1 + 3*(a - 1): 3*a)
   end do

   cm = cm/priv%total_mass

   ! translate reference coordinates to CM frame
   priv%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b) - cm(b)
         priv%ref_nrm = priv%ref_nrm + priv%mass(a)*priv%ref_crd(3*(a - 1) + b)**2
      end do
   end do


end subroutine colvar_bootstrap

!=============================================================================


subroutine read_refcrd(coor, refcrd)

  use nfe_colvar_type
  use nfe_constants
  use nfe_sander_proxy

  implicit none

  NFE_REAL, pointer :: coor(:) !=> null()

  character(len = *), intent(in)  :: refcrd

  integer ::i, j, error
  character :: dummy

  open(REF_UNIT1, FILE = refcrd, iostat = error, status = 'old')
     if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,a/)') &
         NFE_ERROR, 'Failed to open reference coordinates file :''', trim(refcrd), ''''
      call terminate()
     end if
  ! skip the first line of the CRD files 
  read(REF_UNIT1, '(A1)') dummy
  read(REF_UNIT1, *) j                             ! Read total number of atoms
  read(REF_UNIT1, '(6F12.7)') (coor(i), i=1, 3*j)  ! Read coordinates 

  close(REF_UNIT1)

end subroutine read_refcrd

!=============================================================================
subroutine colvar_cleanup(cv)

   NFE_USE_AFAILED

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_SPINANGLE)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)

   nfe_assert(associated(priv%mass))
   nfe_assert(associated(priv%cm_crd))
   nfe_assert(associated(priv%ref_crd))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd)

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

   type(priv_t), pointer :: priv
   integer :: a, c
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_SPINANGLE)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(priv%n_atoms.gt.0.and.priv%n_atoms.eq.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = 1, priv%n_atoms

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= sander_natoms())
      aname = sander_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == priv%n_atoms) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   if (cv_nr.eq.0) &
      write (unit = lun, fmt = '(a,a,a)') &
         NFE_INFO, 'reference coordinates is loaded from : ', trim(refcrd_file)

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 1
   do a = 1, priv%n_atoms
      nfe_assert(c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(c), ', ', cv%r(c + 1), &
         ', ', cv%r(c + 2)
      c = c + 3
   end do

end subroutine print_details

!=============================================================================

end module nfe_cv_SPINANGLE
