! <compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

!
! written by mmoradi-at-ncsu-dot-edu (12/2008)
!
! SIN_OF_DIHEDRAL (sum of sines of dihedral angles)
!
! cv%i = (a11, a12, a13, a14,
!         a21, a22, a23, a24,
!         ...,
!         aN1, aN2, aN3, aN4)
!
!     indexes of the participating atoms
!
! value = sin(d1) + sin(d2) + ... + sin(dN)
!
!   where d1 is dihedral formed by the atoms #a11, #a13, #a13, #a14;
!         d2 is formed by #a21, #a22, #a23, #a24, and so on
!

module nfe_cv_SIN_OF_DIHEDRAL

implicit none

private

!=============================================================================

private :: cos3, cos4, sin4, sin4_d

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

   implicit none

   NFE_REAL :: value

   type(colvar_t), intent(in) :: cv

   NFE_REAL, intent(in) :: x(*)

#  include "nfe-mpi.h"

   integer :: n, a1, a2, a3, a4

   nfe_assert(cv%type == COLVAR_SIN_OF_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   value = ZERO

   do n = 1, size(cv%i)/4

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.sander_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.sander_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.sander_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.sander_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      NFE_MASTER_ONLY_BEGIN
      value = value &
         + sin4(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2), x(a4:a4 + 2))
      NFE_MASTER_ONLY_END

   end do

end function colvar_value

!=============================================================================

subroutine colvar_force(cv, x, fcv, f)

#ifndef NFE_DISABLE_ASSERT
   use nfe_utils
   use nfe_sander_proxy
#endif /* NFE_DISABLE_ASSERT */

   use nfe_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   NFE_REAL, intent(in) :: x(*), fcv

   NFE_REAL, intent(inout) :: f(*)

#  include "nfe-mpi.h"

   NFE_REAL :: d1(3), d2(3), d3(3), d4(3)

   integer :: n, a1, a2, a3, a4

   nfe_assert(cv%type == COLVAR_SIN_OF_DIHEDRAL)
   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   do n = 1, size(cv%i)/4

#     ifdef MPI
      if (mod(n, sandersize).ne.sanderrank) &
         cycle
#     endif /* MPI */

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.sander_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.sander_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.sander_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.sander_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      call sin4_d(x(a1:a1 + 2), x(a2:a2 + 2), &
         x(a3:a3 + 2), x(a4:a4 + 2), d1, d2, d3, d4)

      f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1
      f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2
      f(a3:a3 + 2) = f(a3:a3 + 2) + fcv*d3
      f(a4:a4 + 2) = f(a4:a4 + 2) + fcv*d4

   end do

end subroutine colvar_force

!=============================================================================

subroutine colvar_bootstrap(cv, cvno)

   use nfe_utils
   use nfe_constants
   use nfe_colvar_type
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno

#  include "nfe-mpi.h"

   integer :: n, j, l

   nfe_assert(cv%type.eq.COLVAR_SIN_OF_DIHEDRAL)

   if (.not.associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (SIN_OF_DIHEDRAL) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%i)

   if (mod(size(cv%i), 4).ne.0) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (SIN_OF_DIHEDRAL) : number of integers is not a multiple of 4'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! mod(size(cv%i), 4).ne.0

   if (size(cv%i).lt.4) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (SIN_OF_DIHEDRAL) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%i).lt.4

   do n = 1, size(cv%i)/4
      do j = 0, 3
         if (cv%i(4*n - j).lt.1.or.cv%i(4*n - j).gt.sander_natoms()) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                  (cv%i(4*n - j))//',a,'//pfmt(sander_natoms())//',a/)') &
                  NFE_ERROR, 'CV #', cvno, ' (SIN_OF_DIHEDRAL) : integer #', &
                  (4*n - j), ' (', cv%i(4*n - j), ') is out of range [1, ', &
                  sander_natoms(), ']'
            NFE_MASTER_ONLY_END
            call terminate()
         end if

         do l = 0, 3
            if (l.ne.j.and.cv%i(4*n - l).eq.cv%i(4*n - j)) then
               NFE_MASTER_ONLY_BEGIN
                  write (unit = ERR_UNIT, &
                     fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
                     (4*n - l)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                     (cv%i(4*n - j))//',a/)') NFE_ERROR, 'CV #', &
                     cvno, ' (SIN_OF_DIHEDRAL) : integers #', (4*n - l), &
                     ' and #', (4*n - j), ' are equal (', cv%i(4*n - j), ')'
               NFE_MASTER_ONLY_END
               call terminate()
            end if
         end do
      end do
   end do

end subroutine colvar_bootstrap

!=============================================================================

subroutine print_details(cv, lun)

   use nfe_utils
   use nfe_colvar_type
   use nfe_sander_proxy

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: n, a1, a2, a3, a4
   character(4) :: aname1, aname2, aname3, aname4

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_SIN_OF_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   do n = 1, size(cv%i)/4

      a1 = cv%i(4*n - 3)
      a2 = cv%i(4*n - 2)
      a3 = cv%i(4*n - 1)
      a4 = cv%i(4*n - 0)

      aname1 = sander_atom_name(a1)
      aname2 = sander_atom_name(a2)
      aname3 = sander_atom_name(a3)
      aname4 = sander_atom_name(a4)

      write (unit = lun, fmt = '(a,8x,'//pfmt(a1)//',a,a,a,'//pfmt &
            (a2)//',a,a,a,'//pfmt(a3)//',a,a,a,'//pfmt(a4)//',a,a,a)') &
         NFE_INFO, a1, ' [', trim(aname1), '] ==> ', &
                    a2, ' [', trim(aname2), '] ==> ', &
                    a3, ' [', trim(aname3), '] ==> ', &
                    a4, ' [', trim(aname4), ']'
   end do

end subroutine print_details

!=============================================================================

!
! for A == B == C == D, a torsion along B == C is arccos([ABxBC]*[CDx(-BC)])
!        and its sign is given by sign(BC*[[ABxBC]x[CDx(-BC)]])
!

pure NFE_REAL function sin4(r1, r2, r3, r4)

   use nfe_colvar_math, only : torsion

   implicit none

   NFE_REAL, intent(in) :: r1(3), r2(3), r3(3), r4(3)

   sin4 = sin(torsion(r1, r2, r3, r4))

end function sin4

!=============================================================================

subroutine sin4_d(r1, r2, r3, r4, d1, d2, d3, d4)

   use nfe_colvar_math, only : torsion_d

   implicit none

   NFE_REAL, intent(in) :: r1(3), r2(3), r3(3), r4(3)
   NFE_REAL, intent(out) :: d1(3), d2(3), d3(3), d4(3)

   NFE_REAL :: cos_tor

   call torsion_d(r1, r2, r3, r4, d1, d2, d3, d4)
   
   cos_tor = cos4(r1, r2, r3, r4)

   d1 = d1 * cos_tor
   d2 = d2 * cos_tor
   d3 = d3 * cos_tor
   d4 = d4 * cos_tor

end subroutine sin4_d

!=============================================================================

pure NFE_REAL function cos4(r1, r2, r3, r4)

   implicit none

   NFE_REAL, intent(in) :: r1(3), r2(3), r3(3), r4(3)

   cos4 = cos3(r2 - r1, r3 - r2, r4 - r3)

end function cos4

!=============================================================================

pure NFE_REAL function cos3(v1, v2, v3)

   use nfe_colvar_math, only : cross3, cosine

   implicit none

   NFE_REAL, intent(in) :: v1(3), v2(3), v3(3)

   NFE_REAL :: n1(3), n2(3)

   n1 = cross3(v1, v2)
   n2 = cross3(v2, v3)

   cos3 = cosine(n1, n2)

end function cos3

end module nfe_cv_SIN_OF_DIHEDRAL
