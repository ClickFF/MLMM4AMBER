! <compile=optimized>
!----------------------------------------------------------------------
! Copyright (C) 2004, 2005 Chaok Seok, Evangelos Coutsias and Ken Dill
!      UCSF, Univeristy of New Mexico, Seoul National University
! Witten by Chaok Seok and Evangelos Coutsias 2004.
! modified by Mahmoud Moradi and Volodymyr Babin at NFE, 2007
! The original version is at http://www.dillgroup.ucsf.edu/rmsd

! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!

! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
! MA  02110-1301  USA
!-----------------------------------------------------------------------

!
! make up by Mahmoud Moradi and Volodymyr Babin at NFE, Cox 308
!

#include "nfe-utils.h"
#include "nfe-config.h"

module nfe_rmsd

!=============================================================================

implicit none

private

!=============================================================================

#ifdef NFE_ENABLE_RMSD_CANNED
public :: rmsd_canned
#endif /* NFE_ENABLE_RMSD_CANNED */

public :: rmsd_q
public :: orientation_q
public :: rmsd_q1
public :: rmsd_q2u
public :: rmsd_q3u
!=============================================================================

contains

!=============================================================================

!
! a "driver" in FORTRAN slang (for testing only)
!

#ifdef NFE_RMSD_CANNED
NFE_REAL function rmsd_canned(m, x1, x2, g)

   use nfe_utils
   use nfe_constants

   implicit none

   NFE_REAL, intent(in)    :: m(:)
   NFE_REAL, intent(inout) :: x1(:), x2(:) ! destroyed upon return

   NFE_REAL, intent(out), optional :: g(:)

   integer   :: n, a, a3, i
   NFE_REAL :: mass, c1(3), c2(3), q(4), lambda, x1n, x2n, tmp, U(3,3)

   n = size(m)

   nfe_assert(n > 1)

   nfe_assert(3*n == size(x1))
   nfe_assert(3*n == size(x2))

   !
   ! find centers of mass
   !

   c1 = ZERO
   c2 = ZERO

   mass = ZERO

   do a = 1, n
      a3 = 3*(a - 1)
      mass = mass + m(a)
      do i = 1, 3
         c1(i) = c1(i) + m(a)*x1(a3 + i)
         c2(i) = c2(i) + m(a)*x2(a3 + i)
      end do
   end do

   nfe_assert(mass > ZERO)

   c1 = c1/mass
   c2 = c2/mass

   ! center x1, x2 && find "norms"

   x1n = ZERO
   x2n = ZERO

   do a = 1, n
      a3 = 3*(a - 1)
      do i = 1, 3
         x1(a3 + i) = x1(a3 + i) - c1(i)
         x1n = x1n + m(a)*x1(a3 + i)**2
         x2(a3 + i) = x2(a3 + i) - c2(i)
         x2n = x2n + m(a)*x2(a3 + i)**2
      end do
   end do

   call rmsd_q(n, m, x1, x2, lambda, q)

   rmsd_canned = sqrt(max(ZERO, ((x1n + x2n) - 2*lambda))/mass)

   ! g is w.r.t x1

   if (present(g)) then
      nfe_assert(3*n == size(g))
      call rmsd_q2u(q, U)

      tmp = ONE/mass/max(rmsd_canned, NFE_TO_REAL(0.000001))
      do a = 1, n
         a3 = 3*a
         g(a3 - 2:a3) = m(a)*tmp*(x1(a3 - 2:a3) - matmul(U, x2(a3 - 2:a3)))
      end do
   end if ! present g

end function rmsd_canned
#endif /* NFE_RMSD_CANNED */

!=============================================================================

!
! finds optimal rotation; size(w) == n; size(x?) == 3*n
!

subroutine rmsd_q(n, w, x1, x2, lambda, q)

   use nfe_constants

   implicit none

   integer,   intent(in) :: n
   NFE_REAL, intent(in) :: w(*), x1(*), x2(*)

   NFE_REAL, intent(out) :: lambda, q(4)

   integer   :: a, a3, i, j
   NFE_REAL :: R(4,4), S(4,4)

   ! calculate the R matrix

   R = ZERO

   do a = 0, n - 1
      a3 = 3*a
      do i = 1, 3
         do j = 1, 3
            R(i,j) = R(i,j) + w(a + 1)*x1(a3 + i)*x2(a3 + j)
         end do
      end do
   end do   


   ! S matrix

   S(1,1) = R(1,1) + R(2,2) + R(3,3)
   S(2,1) = R(2,3) - R(3,2)
   S(3,1) = R(3,1) - R(1,3)
   S(4,1) = R(1,2) - R(2,1)

   S(1,2) = S(2,1)
   S(2,2) = R(1,1) - R(2,2) - R(3,3)
   S(3,2) = R(1,2) + R(2,1)
   S(4,2) = R(1,3) + R(3,1)

   S(1,3) = S(3,1)
   S(2,3) = S(3,2)
   S(3,3) =-R(1,1) + R(2,2) - R(3,3)
   S(4,3) = R(2,3) + R(3,2)

   S(1,4) = S(4,1)
   S(2,4) = S(4,2)
   S(3,4) = S(4,3)
   S(4,4) =-R(1,1) - R(2,2) + R(3,3) 

   call dstmev(S, lambda, q)

end subroutine rmsd_q

subroutine orientation_q(n, w, x1, x2, lambda, q)   ! no need for w -- we have already considered 
   use nfe_constants

   implicit none

   integer,   intent(in) :: n
   NFE_REAL, intent(in) :: w(*), x1(*), x2(*)
   NFE_REAL, intent(out) :: lambda(4), q(4,4)
   NFE_REAL :: R(4,4), S(4,4)
   integer   :: a, a3, i, j

   ! calculate the R matrix (the correlation matrix)
   R = ZERO
   do a = 0, n - 1
      a3 = 3*a
      do i = 1, 3
         do j = 1, 3
            !R(i,j) = R(i,j) + w(a+1)*x1(a3 + i)*x2(a3 + j)
             R(i,j) = R(i,j) + x1(a3 + i)*x2(a3 + j)
         end do
      end do
   end do
   
   ! S matrix (overlap matrix)

   S(1,1) = R(1,1) + R(2,2) + R(3,3)
   S(2,1) = R(2,3) - R(3,2)
   S(3,1) = R(3,1) - R(1,3)
   S(4,1) = R(1,2) - R(2,1)

   S(1,2) = S(2,1)
   S(2,2) = R(1,1) - R(2,2) - R(3,3)
   S(3,2) = R(1,2) + R(2,1)
   S(4,2) = R(1,3) + R(3,1)

   S(1,3) = S(3,1)
   S(2,3) = S(3,2)
   S(3,3) =-R(1,1) + R(2,2) - R(3,3)
   S(4,3) = R(2,3) + R(3,2)

   S(1,4) = S(4,1)
   S(2,4) = S(4,2)
   S(3,4) = S(4,3)
   S(4,4) =-R(1,1) - R(2,2) + R(3,3) 
   call dstmev4(S, lambda, q)


end subroutine orientation_q

!=============================================================================

! add variable control

subroutine rmsd_q1(n, state, w, x1, x2, lambda, q)

   use nfe_constants

   implicit none

   integer,   intent(in) :: n, state(*)
   NFE_REAL, intent(in) :: w(*), x1(*), x2(*)

   NFE_REAL, intent(out) :: lambda, q(4)

   integer   :: a, a3, i, j
   NFE_REAL :: R(4,4), S(4,4)

   ! calculate the R matrix

   R = ZERO

   do a = 0, n - 1
      if(state(a+1)==0) cycle
      a3 = 3*a
      do i = 1, 3
         do j = 1, 3
            R(i,j) = R(i,j) + w(a + 1)*x1(a3 + i)*x2(a3 + j)
         end do
      end do
   end do

   ! S matrix

   S(1,1) = R(1,1) + R(2,2) + R(3,3)
   S(2,1) = R(2,3) - R(3,2)
   S(3,1) = R(3,1) - R(1,3)
   S(4,1) = R(1,2) - R(2,1)

   S(1,2) = S(2,1)
   S(2,2) = R(1,1) - R(2,2) - R(3,3)
   S(3,2) = R(1,2) + R(2,1)
   S(4,2) = R(1,3) + R(3,1)

   S(1,3) = S(3,1)
   S(2,3) = S(3,2)
   S(3,3) =-R(1,1) + R(2,2) - R(3,3)
   S(4,3) = R(2,3) + R(3,2)

   S(1,4) = S(4,1)
   S(2,4) = S(4,2)
   S(3,4) = S(4,3)
   S(4,4) =-R(1,1) - R(2,2) + R(3,3) 

   call dstmev(S, lambda, q)

end subroutine rmsd_q1

!=============================================================================


!
! This subroutine constructs (transposed) rotation matrix U from quaternion q.
! 

subroutine rmsd_q2u(q, U)

   use nfe_constants

   implicit none

   NFE_REAL, intent(in)  :: q(*) ! 4
   NFE_REAL, intent(out) :: U(3,3)

   NFE_REAL :: b0, b1, b2, b3
   NFE_REAL :: q00, q01, q02, q03
   NFE_REAL :: q11, q12, q13, q22, q23, q33

   b0 = q(1) + q(1)
   b1 = q(2) + q(2)
   b2 = q(3) + q(3)
   b3 = q(4) + q(4)

   q00 = b0*q(1) - ONE
   q01 = b0*q(2)
   q02 = b0*q(3)
   q03 = b0*q(4)

   q11 = b1*q(2)
   q12 = b1*q(3)
   q13 = b1*q(4)

   q22 = b2*q(3)
   q23 = b2*q(4)

   q33 = b3*q(4)

   U(1,1) = q00 + q11
   U(2,1) = q12 - q03
   U(3,1) = q13 + q02

   U(1,2) = q12 + q03
   U(2,2) = q00 + q22
   U(3,2) = q23 - q01

   U(1,3) = q13 - q02
   U(2,3) = q23 + q01
   U(3,3) = q00 + q33

end subroutine rmsd_q2u

!=============================================================================
!
! This subroutine constructs (UNtransposed) rotation matrix U from quaternion q.
! 

subroutine rmsd_q3u(q, U)

   use nfe_constants

   implicit none

   NFE_REAL, intent(in)  :: q(*) ! 4
   NFE_REAL, intent(out) :: U(3,3)

   NFE_REAL :: b0, b1, b2, b3
   NFE_REAL :: q00, q01, q02, q03
   NFE_REAL :: q11, q12, q13, q22, q23, q33

   b0 = q(1) + q(1)
   b1 = q(2) + q(2)
   b2 = q(3) + q(3)
   b3 = q(4) + q(4)

   q00 = b0*q(1) - ONE
   q01 = b0*q(2)
   q02 = b0*q(3)
   q03 = b0*q(4)

   q11 = b1*q(2)
   q12 = b1*q(3)
   q13 = b1*q(4)

   q22 = b2*q(3)
   q23 = b2*q(4)

   q33 = b3*q(4)

   U(1,1) = q00 + q11
   U(1,2) = q12 - q03
   U(1,3) = q13 + q02

   U(2,1) = q12 + q03
   U(2,2) = q00 + q22
   U(2,3) = q23 - q01

   U(3,1) = q13 - q02
   U(3,2) = q23 + q01
   U(3,3) = q00 + q33

end subroutine rmsd_q3u

!=============================================================================
!
! a simple subroutine to compute the leading eigenvalue and eigenvector
! of a symmetric, traceless 4x4 matrix A by an inverse power iteration:
! (1) the matrix is converted to tridiagonal form by 3 Givens
! rotations;  V*A*V' = T
! (2) Gershgorin's theorem is used to estimate a lower
! bound for the leading negative eigenvalue:
! lambda_1 > g=min(T11-t12,-t21+T22-t23,-t32+T33-t34,-t43+T44)
!          =
! where tij=abs(Tij)
! (3) Form the positive definite matrix 
!     B = T-gI
! (4) Use svd (algorithm svdcmp from "Numerical Recipes")
!     to compute eigenvalues and eigenvectors for SPD matrix B
! (5) Shift spectrum back and keep leading singular vector
!     and largest eigenvalue.
! (6) Convert eigenvector to original matrix A, through 
!     multiplication by V'.  
!

subroutine dstmev(A, lambda, evec)

   implicit none

   NFE_REAL, intent(in)  :: A(4,4)
   NFE_REAL, intent(out) :: lambda, evec(4)

   NFE_REAL :: T(4,4), V(4,4), SV(4,4), SW(4)
   NFE_REAL :: U(4,4), SVT(4,4), work(20)
   integer :: lwork, info, iwork

   integer :: i!, max_loc(1)

   ! (I) Convert to tridiagonal form, keeping similarity transform
   !            (a product of 3 Givens rotations)
   call givens4(A, T, V)

   ! (II) Estimate lower bound of smallest eigenvalue by Gershgorin's theorem
   lambda = min(T(1,1) - abs(T(1,2)), &
                T(2,2) - abs(T(2,1)) - abs(T(2,3)), &
                T(3,3) - abs(T(3,2)) - abs(T(3,4)), &
                T(4,4) - abs(T(4,3)))

   ! (III) Form positive definite matrix  T = lambda*I - T
   do i = 1, 4
      T(i,i) = T(i,i) - lambda
   end do

   lwork=size(work)
   iwork=size(work)
   info=0

   ! Use lapack routine instead, svdcmp is not numerically stable at
   ! higher optimization levels
   call D_OR_S()gesvd('A', 'A', 4, 4, T, 4, SW, U, 4, SVT, 4, work, lwork, info)
   ! (IV) Compute singular values/vectors of SPD matrix B
   ! call svdcmp(T, SW, SV)

   ! (V) Shift spectrum back
   ! xgesvd always returns largest value in SW(1)
   ! max_loc = maxloc(SW) 
   !lambda = lambda + SW(max_loc(1))
   lambda = lambda + SW(1)

   ! xgesvd returns a transposed SV compared to the svdcmp routine
   SV = transpose(SVT)

   ! (VI) Convert eigenvector to original coordinates: (V is transposed!)
   evec = matmul(V, SV(:, 1))
   !evec = matmul(V, SV(:, max_loc(1)))

end subroutine dstmev

subroutine dstmev4(A, lambda, evec)

   implicit none

   NFE_REAL, intent(in)  :: A(4,4)
   NFE_REAL, intent(out) :: lambda(4), evec(4,4)

   !NFE_REAL :: T(4,4), V(4,4), SVT(4,4), SW(4)
   NFE_REAL :: U(4,4), SV(4,4), work(20), SW(4), norm1, norm2
   integer :: lwork, info, iwork, nrot

   integer :: ip, iq

   ! (I) Convert to tridiagonal form, keeping similarity transform
   !            (a product of 3 Givens rotations)
!   call givens4(A, T, V)

   ! (II) Estimate lower bound of smallest eigenvalue by Gershgorin's theorem    !smallest EV
!   lambda(1) = min(T(1,1) - abs(T(1,2)), &
!                T(2,2) - abs(T(2,1)) - abs(T(2,3)), &
!                T(3,3) - abs(T(3,2)) - abs(T(3,4)), &
!                T(4,4) - abs(T(4,3)))

   ! (III) Form positive definite matrix  T = lambda*I - T
!   do i = 1, 4
!      T(i,i) = T(i,i) - lambda(1)
!   end do
   lwork=size(work)
   iwork=size(work)
   info=0

   ! Use lapack routine instead, svdcmp is not numerically stable at
   ! higher optimization levels

   ! SW are E-vals, U (or SVT) are E-vecs, and T(or A-dont need to reduce to T) our matrix
   !call D_OR_S()gesvj('G', 'U', 'A', 4, 4, A, 4, SW, 4, SVT, 4, work, lwork, info)
   !call D_OR_S()gesvd('A', 'A', 4, 4, A, 4, SW, U, 4, SVT, 4, work, lwork, info)
   !call D_OR_S()gejsv('F', 'F', 'N', 'R', 'N', 'N', 4, 4, A, 4, SW, U, 4, SVT, 4, work, lwork, iwork, info)
   ! Be independent -- define our own jacobi
   call jacobi(A, 4, 4, SW, U, nrot)
   ! sort, transpose, and normalize
   call eigsrt (SW, U, 4, 4)
   !gejsvd
   !write(*,*), 'info= ', info
   ! (IV) Compute singular values/vectors of SPD matrix B
   !call svdcmp(T, SW, SV)

   ! (V) Shift spectrum back
   ! xgesvd always returns largest value in SW(1)
   
   !IF( info.GT.0 ) THEN
   !      WRITE(*,*)'The algorithm computing SVD failed to converge.'
   !      STOP
   !END IF

   !SV = transpose(SVT)
   !SV = A
   lambda(1) =  SW(1)
   lambda(2) =  SW(2)
   lambda(3) =  SW(3)
   lambda(4) =  SW(4)
   SV = transpose (U)
   !SV = transpose(SVT)
   !SV = U
  
   ! Normalization
   do ip = 1, 4
      norm2 = 0.0
      do iq = 1, 4 
         norm2 = norm2 + SV(ip,iq)*SV(ip,iq)
      enddo
      norm1 = sqrt(norm2) 
      do iq = 1, 4
         SV(:,iq) = SV(:,iq)/norm1 
      enddo
   enddo
 
   !norm1 = sqrt( SV(1,1)*SV(1,1) + SV(2,1)*SV(2,1)+ &
   !        SV(3,1)*SV(3,1) + SV(4,1)*SV(4,1))
   !SV(:,1) = SV(:,1)/norm1
   
   !norm2 = sqrt(SV(1,2)*SV(1,2) + SV(2,2)*SV(2,2)+ &
   !        SV(3,2)*SV(3,2) + SV(4,2)*SV(4,2))
   !SV(:,2) = SV(:,2)/norm2
   
   !norm3 = sqrt(SV(1,3)*SV(1,3) + SV(2,3)*SV(2,3)+ &
   !        SV(3,3)*SV(3,3) + SV(4,3)*SV(4,3))
   !SV(:,3) = SV(:,3)/norm3
   
   !norm4 = sqrt(SV(1,4)*SV(1,4) + SV(2,4)*SV(2,4)+ &
   !        SV(3,4)*SV(3,4) + SV(4,4)*SV(4,4))
   !SV(:,4) = SV(:,4)/norm4

   ! (VI) Convert eigenvector to original coordinates: (V is transposed!)
   evec(:,1) = SV(:, 1)
   evec(:,2) = SV(:, 2)
   evec(:,3) = SV(:, 3)
   evec(:,4) = SV(:, 4)
   !if (evec(1,1) >= 0.0) then             !This is in go 'nfe-cv-QUATERNION'! 
   !    evec(:,1) = evec(:,1)
   !else
   !   evec(:,1) = -1.0 * evec(:,1)
   !endif
   
end subroutine dstmev4


!=============================================================================

!
! performs givens rotations to reduce symmetric 4x4 matrix to tridiagonal
!

subroutine givens4(S, T, V)

   use nfe_constants

   implicit none

   NFE_REAL, dimension(4,4), intent(in)  :: S
   NFE_REAL, dimension(4,4), intent(out) :: T,V

   NFE_REAL :: c1, c2, c3, s1, s2, s3, r1, r2, r3, c1c2, s1c2

   T = S
   V = ZERO

   ! zero out entries T(4,1) and T(1,4)
   ! compute cos and sin of rotation angle in the 3-4 plane

   r1 = pythag(T(3,1), T(4,1))

   if (r1 .ne. ZERO) then
      c1 = T(3,1)/r1
      s1 = T(4,1)/r1

      V(3,3) = c1
      V(3,4) = s1
      V(4,3) =-s1
      V(4,4) = c1

      T(3,1) = r1
      T(4,1) = ZERO

      T(3:4,2:4) = matmul(V(3:4,3:4), T(3:4,2:4))
      T(1:2,3:4) = transpose(T(3:4,1:2))
      T(3:4,3:4) = matmul(T(3:4,3:4), transpose(V(3:4,3:4)))
   else
      c1 = ONE
      s1 = ZERO
   end if

   ! zero out entries T(3,1) and T(1,3)
   ! compute cos and sin of rotation angle in the 2-3 plane

   r2 = pythag(T(3,1), T(2,1))

   if (r2 .ne. ZERO) then
      c2 = T(2,1)/r2
      s2 = T(3,1)/r2

      V(2,2) = c2
      V(2,3) = s2
      V(3,2) =-s2
      V(3,3) = c2

      T(2,1) = r2
      T(3,1) = ZERO

      T(2:3,2:4) = matmul(V(2:3,2:3), T(2:3,2:4))
      T(1,2:3)   = T(2:3,1)
      T(4,2:3)   = T(2:3,4)
      T(2:3,2:3) = matmul(T(2:3,2:3), transpose(V(2:3,2:3)))
   else
      c2 = ONE
      s2 = ZERO
   end if

   ! zero out entries T(4,2) and T(2,4)
   ! compute cos and sin of rotation angle in the 3-4 plane

   r3 = pythag(T(4,2), T(3,2))

   if (r3 .ne. ZERO) then
      c3 = T(3,2)/r3
      s3 = T(4,2)/r3

      V(3,3) = c3
      V(3,4) = s3
      V(4,3) =-s3
      V(4,4) = c3

      T(3,2) = r3
      T(4,2) = ZERO

      T(3:4,3:4) = matmul(V(3:4,3:4), T(3:4,3:4))
      T(1:2,3:4) = transpose(T(3:4,1:2))
      T(3:4,3:4) = matmul(T(3:4,3:4), transpose(V(3:4,3:4)))
   else
      c3 = ONE
      s3 = ZERO
   end if

   ! compute net rotation matrix (accumulate similarity for evec. computation)
   ! To save transposing later, This is the transpose!

   V(1,1)   = ONE
   V(1,2:4) = ZERO
   V(2:4,1) = ZERO

   V(2,2) = c2
   V(3,2) = c1*s2
   V(4,2) = s1*s2

   c1c2 = c1*c2
   s1c2 = s1*c2

   V(2,3) =-s2*c3
   V(3,3) = c1c2*c3 - s1*s3
   V(4,3) = s1c2*c3 + c1*s3
   V(2,4) = s2*s3
   V(3,4) =-c1c2*s3 - s1*c3
   V(4,4) =-s1c2*s3 + c1*c3

end subroutine givens4

!============================================================================
subroutine jacobi(a, n, np, d, v, nrot)

   use nfe_constants
   implicit none

   integer :: n, np, nrot
   NFE_REAL :: a(np,np), d(np), v(np,np)
   integer :: i, ip, iq, j
   NFE_REAL :: c, g, h, s, sm, t, tau, theta, tresh, b(n), z(n) 

   do ip = 1, n                !initialze to the identity
      do iq = 1, n
         v(ip, iq) = 0.0
      enddo
      v(ip, ip) = 1.0
   enddo
   do ip = 1, n
      b(ip) = a(ip,ip)       ! initialze b and d to the diagnol of a
      d(ip) = b(ip)
      z(ip) = 0.0
   enddo
   nrot = 0
   do i = 1, 50
      sm = 0.0
      do ip = 1, n-1
         do iq = ip+1, n
            sm = sm + abs(a(ip,iq))
         enddo
      enddo 
      if (sm .eq. ZERO) return
      if (i .lt. 4) then
         tresh = 0.2*sm/n**2
      else
         tresh = 0.0
      endif
      do ip =1, n-1
         do iq = ip+1, n 
            g = 100.0*abs(a(ip,iq))
            if ((i .gt. 4) .and. (abs(d(ip))+ &
            g.eq.abs(d(ip))).and.(abs(d(iq))+ &
            g .eq. abs(d(iq)))) then
               a(ip,iq) = 0.0
            else if (abs(a(ip,iq)) .gt. tresh) then
                 h = d(iq) - d(ip)
                 if (abs(h) + g .eq. abs(h)) then
                    t = a(ip,iq)/h
                 else
                    theta = 0.5*h/a(ip,iq)
                    t = 1.0/(abs(theta) + sqrt(1.0+theta**2))
                    if (theta .lt. ZERO) t = -t
                 endif
                 c = 1.0/sqrt(1.0+t**2)
                 s = t*c
                 tau = s/(1.0+c)
                 h = t*a(ip,iq)
                 z(ip) = z(ip) - h
                 z(iq) = z(iq) + h
                 d(ip) = d(ip) - h
                 d(iq) = d(iq) + h
                 a(ip,iq) = 0.0
                 do j =1, ip-1
                    g = a(j,ip)
                    h = a(j,iq)
                    a(j,ip) = g - s*(h+g*tau)
                    a(j,iq) = h + s*(g-h*tau)
                 enddo
                 do j = ip+1, iq-1
                    g = a(ip,j)
                    h = a(j,iq)
                    a(ip,j) = g-s*(h+g*tau)
                    a(j, iq) = h+s*(g-h*tau)
                 enddo
                 do j = iq+1,n
                    g = a(ip,j)
                    h = a(iq,j)
                    a(ip,j) = g-s*(h+g*tau)
                    a(iq,j) = h+s*(g-h*tau)
                 enddo
                 do j =1, n
                 g = v(j,ip)
                 h = v(j,iq)
                 v(j,ip) = g-s*(h+g*tau)
                 v(j, iq) = h+s*(g-h*tau)
                 enddo
                 nrot = nrot + 1
            endif
         enddo
      enddo
      do ip = 1, n
         b(ip) = b(ip) + z(ip)
         d(ip) = b(ip)
         z(ip) = 0.0
      enddo
   enddo
   return
end subroutine jacobi


!=============================================================================

subroutine eigsrt(d, v, n, np)

   use nfe_constants
   implicit none

   integer :: n, np, i, j, k
   NFE_REAL :: d(np), v(np,np), p

   do i = 1, n-1
      k = i
      p = d(i)
      do j = i+1, n
         if (d(j) .ge. p)then
            k = j
            p = d(j)
         endif
      enddo
      if (k .ne. i)then
          d(k) = d(i)
          d(i) = p
          do j = 1, n
             p = v(j,i)
             v(j,i) = v(j,k)
             v(j,k) = p
          enddo
      endif
   enddo
   return
end subroutine eigsrt

!=============================================================================

subroutine svdcmp(a, w, v)

   use nfe_utils
   use nfe_constants

   implicit none

   integer, parameter :: N = 4

   NFE_REAL, intent(inout) :: a(N,*)
   NFE_REAL, intent(out)   :: v(N,*), w(*)

   integer :: i, its, j, jj, k, l, nm

   NFE_REAL :: anorm, c, f, g, h, s, scale, x, y, z, rv1(2*N)

   g = ZERO
   scale = ZERO
   anorm = ZERO

   nm = 0 ! for g95

   do i = 1, N

      l = i + 1
      rv1(i) = scale*g

      g = ZERO
      s = ZERO
      scale = ZERO

      do k = i, N
         scale = scale + abs(a(k,i))
      end do

      if (scale .ne. ZERO) then
         do k = i, N
            a(k,i) = a(k,i)/scale
            s = s + a(k,i)*a(k,i)
         end do

         f = a(i,i)
         g =-sign(sqrt(s),f)
         h = f*g - s
         a(i,i) = f - g

         do j = l, N 
            s = ZERO
            do k = i, N
               s = s + a(k,i)*a(k,j)
            end do

            f = s/h
            do k = i, N
               a(k,j) = a(k,j) + f*a(k,i)
            end do
         end do

         do k = i, N
            a(k,i) = scale*a(k,i)
         end do
      endif ! scale .ne. ZERO

      w(i) = scale*g
      g = ZERO
      s = ZERO
      scale = ZERO

      if (i .ne. N) then
         do k = l, N
            scale = scale + abs(a(i,k))
         end do
         if (scale .ne. ZERO) then
            do k = l, N
               a(i,k) = a(i,k)/scale
               s = s + a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g =-sign(sqrt(s),f)
            h = f*g - s
            a(i,l) = f - g
            do k = l, N
               rv1(k) = a(i,k)/h
            end do
            do j = l, N
               s = ZERO
               do k = l, N
                  s = s + a(j,k)*a(i,k)
               end do
               do k = l, N
                  a(j,k) = a(j,k) + s*rv1(k)
               end do
            end do
            do k = l, N
               a(i,k) = scale*a(i,k)
            end do
         endif
      endif
      anorm = max(anorm, (abs(w(i)) + abs(rv1(i))))
   end do

   do i = N, 1, -1
      if (i .lt. N) then
         if (g .ne. ZERO) then
            do j = l, N
               v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j = l, N
               s = ZERO
               do k = l, N
                  s = s + a(i,k)*v(k,j)
               end do
               do k = l, N
                  v(k,j) = v(k,j) + s*v(k,i)
               end do
            end do
         endif ! g .ne. ZERO
         do j = l, N
            v(i,j) = ZERO
            v(j,i) = ZERO
         end do
      endif
      v(i,i) = ONE
      g = rv1(i)
      l = i
   end do

   do i = N, 1, -1
      l = i + 1
      g = w(i)
      do j = l, N
         a(i,j) = ZERO
      end do
      if (g .ne. ZERO) then
         g = ONE/g
         do j = l, N
            s = ZERO
            do k = l, N
               s = s + a(k,i)*a(k,j)
            end do
            f = (s/a(i,i))*g
            do k = i, N
               a(k,j) = a(k,j) + f*a(k,i)
            end do
         end do
         do j = i, N
            a(j,i) = a(j,i)*g
         end do
      else
         do j = i, N
            a(j,i) = ZERO
         end do
      endif ! g .ne. ZERO
      a(i,i) = a(i,i) + ONE
   end do

   do k = N, 1, -1
      do its = 1, 30
         do l = k, 1, -1
            nm = l - 1
            if ((abs(rv1(l)) + anorm) .eq. anorm) &
               goto 2
            if ((abs(w(nm)) + anorm) .eq. anorm) &
               goto 1
         end do
1        c = ZERO
         s = ONE
         do i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((abs(f) + anorm) .eq. anorm) &
               goto 2
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            h = ONE/h
            c = (g*h)
            s =-(f*h)
            do j = 1, N     
               y = a(j,nm)
               z = a(j,i)
               a(j,nm) = (y*c) + (z*s)
               a(j,i) = -(y*s) + (z*c)
            end do
         end do
2        z = w(k)
         if (l .eq. k) then
            if (z .lt. ZERO) then
               w(k) = -z
               do j = 1, N
                  v(j,k) = -v(j,k)
               end do
            end if
            goto 3
         end if
         if (its .eq. 30) &
            call fatal('nfe-rmsd: no convergence in svdcmp()')
         x = w(l)
         nm = k - 1
         y = w(nm)
         g = rv1(nm)
         h = rv1(k)
         f = ((y - z)*(y + z) + (g - h)*(g + h))/(2*h*y)
         g = pythag(f, ONE)
         f = ((x - z)*(x + z) + h*((y/(f + sign(g, f))) - h))/x
         c = ONE
         s = ONE
         do j = l, nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g    
            z = pythag(f, h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c) + (g*s)
            g =-(x*s) + (g*c)
            h = y*s
            y = y*c
            do jj = 1, N
               x = v(jj,j)
               z = v(jj,i)
               v(jj,j) = (x*c) + (z*s)
               v(jj,i) =-(x*s) + (z*c)
            end do
            z = pythag(f, h)
            w(j) = z
            if (z .ne. ZERO) then
               z = ONE/z
               c = f*z
               s = h*z
            end if
            f = (c*g) + (s*y)
            x =-(s*g) + (c*y)
            do jj = 1, N
               y = a(jj,j)
               z = a(jj,i)
               a(jj,j) = (y*c) + (z*s)
               a(jj,i) =-(y*s) + (z*c)
            end do
         end do
         rv1(l) = ZERO
         rv1(k) = f
         w(k) = x
      end do
3     continue
   end do

end subroutine svdcmp

!============================================================================

!
! computes sqrt(a**2 + b**2) carefully
!

pure NFE_REAL function pythag(a, b)

   use nfe_constants

   implicit none

   NFE_REAL, intent(in) :: a, b

   NFE_REAL :: absa, absb

   absa = abs(a)
   absb = abs(b)

   if (absa .gt. absb) then
      pythag = absa*sqrt(ONE + (absb/absa)**2)
   else
      if (absb .eq. ZERO) then
         pythag = ZERO
      else
         pythag = absb*sqrt(ONE + (absa/absb)**2)
     endif
   endif

end function pythag
!============================================================================

!
! simple E-value sorting 
!

!return p,q in ascending order
Subroutine Order(p,q)

   NFE_REAL :: p,q,temp

    if (p<q) then
       temp=p
       p=q
       q=temp
    end if
  return

end

!sorting of array A
Subroutine sort(A, n)

   NFE_REAL :: A(1:n)
   integer :: i,j,n 

    do i=1, n
      do j=n, i+1, -1
        call Order(A(j-1), A(j))
      end do
    end do
  return
end


!============================================================================

end module nfe_rmsd
