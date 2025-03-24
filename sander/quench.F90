!<compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
!------------------------------------------------------------------------------
! quench: subroutine to compute quenched MD in Nudged Elastic Band simulations.
!
! Arguments:
!   f:      forces on all particles
!   v:      velocities of all particles
!------------------------------------------------------------------------------
subroutine quench(f, v) 
   implicit none
   
#include "../include/md.h" 
!need access to vv - temp verlet scaling
#include "../include/memory.h" 
!need access to natom

  _REAL_ f(*), v(*), dotproduct, force
  ! f is the forces and v is the velocity
 
  integer index
  dotproduct = 0.d0
  force = 0.d0

  do index = 1, 3*natom
    force = force + f(index)**2
    dotproduct = dotproduct + v(index)*f(index)
  enddo
 
  if (force .ne. 0.0d0) then 
    force = 1.0d0 / sqrt(force)
    dotproduct = dotproduct*force
  end if
   
  if (dotproduct > 0.0d0) then
    v(1:3*natom) = dotproduct * f(1:3*natom) * force
  else 
    v(1:3*natom) = vfac * dotproduct * f(1:3*natom) * force
  end if
   
end subroutine quench
