! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! RATTLEV: correct velocities with bond length constraints
subroutine rattlev( nrr,nbonh,nbona,nbper,ib,jb,igrp,winv,conp, &
      skip,x,v,niter,belly,ifstwt,noshake,qspatial )

   implicit none
   logical skip(*),ready,first
   logical belly
   
   !     ----- ALL THE BONDS INVOLVING HYDROGEN ATOMS ARE LOADED FIRST
   !           IN ARRAYS IB AND JB FOLLOWED BY THOSE INVOLVING NON-
   !           HYDROGEN ATOMS AND THE PERTURBED ATOMS -----
   
#ifdef MPI
#  include "parallel.h"
#endif
#include "../include/md.h"
   
   integer nrr,nbonh,nbona,nbper ! nrr will be number of atoms;
                                 ! nbonh number of bonds to H
                                 ! nbona other bonds
                                 ! nbper should (now) always be zero
   _REAL_ winv(*)
   integer ib(*),jb(*),igrp(*)
   _REAL_ conp(*),x(*),v(*)
   _REAL_ xij(3),vij(3)
   integer ifstwt(*),noshake(*)
   logical :: qspatial

   _REAL_ winvi,winvj,diff,acor,xh,toler,rvdot, &
         tol2,zero
   integer m,i,j,i3,j3,niter,nbt,nit,k,ll
   integer nbt_new,bnd_first,bnd_last
   data first /.true./
   data zero /0.0d0/
   save first,zero,nbt_new
   integer llind
#ifdef MPI
   integer istart, iend
#endif

   if(qspatial)then
      if(ntc == 1)then
         niter=1
         return
      end if
   else
      if(ntc == 1)then
         niter=1
         return
      else if(ntc == 2)then
         nbt = nbonh
      else
         nbt = nbonh+nbona+nbper
      end if
      
      !     first time through, figure out last bond that is not a 3-point water
      !     bond:
      
      if (first) then
         nbt_new = 0
         do i=1,nbt
            if(ifstwt(i) == 0 ) nbt_new = i
         end do
         first= .false.
      end if
      if( nbt_new == 0 ) then
         niter=1
         return
      end if
   end if

   niter = 0
   tol2 = tol / dt
   
   nit = 0
   do k = 1,nrr
      skip(k) = .true.
      skip(nrr+k) = .false.
   end do
   ready = .false.
   
   !===========================================================================
   !        BIG ITERATION LOOP
   
   do while( .not. ready )
      if( nit > 3000 ) then
         write(6,311)
         write(6,322)
         if (imin == 1) write(6,*) &
               ' *** Especially for minimization, try ntc=1 (no shake)'
         call mexit(6,1)
      end if

      ready = .true.
#ifdef MPI
      
      !     ---- divide atoms up among the processors, always splitting on
      !     residue boundaries:

      ! DAN ROE: Modify this so that REMD with 1 cpu / replica is ok      
      if (ntc /= 2 .and. numtasks.gt.1) then
         write(6,*) 'this parallel version only works for ntc < 3'
         write(6,*) '      only fast water and bonds with H.'
         call mexit(6,1)
      end if
      istart = iparpt(mytaskid) + 1
      iend = iparpt(mytaskid+1)
#endif
      
      !     ----- LOOP OVER ALL THE BONDS that are not 3-point waters: -----
      
      bnd_first=1
      bnd_last=nbt_new

      bonds: do llind = bnd_first,bnd_last
         !  Skip bonds constrained in fast 3-point routine:
         if (ifstwt(llind) == 1) cycle bonds
         !  Skip bonds identified by the noshake user input:
         if (noshake(llind) == 1) cycle bonds
         ll=llind
         i3 = ib(ll)
         i  = i3/3+1
         j3 = jb(ll)
         j  = j3/3+1
         if (skip(nrr+i) .and. skip(nrr+j)) cycle bonds ! (already converged)

#ifdef MPI
         !     --- skip constrained bonds not destined for this processor:
         
         if( .not. mpi_orig ) then
            if (i < istart .or. i > iend) cycle bonds
            if (j < istart .or. j > iend) then
               write(6,*) 'partition error in shake on processor ', mytaskid
               write(6,*)'this processor has atoms ',istart,' through ',iend
               write(6,*) 'atom ', i, 'is within this range'
               write(6,*) 'atom ', j, 'is not within this range !'
               call mexit(6,1)
            end if
         end if
#endif
         
         toler = conp(ll)  ! square of bond length
         
         !     ----- IF BELLY OPTION IS ON THEN THE RESETTING OF THE FROZEN ATOM
         !              IS TO BE PREVENTED -----
         
         winvi = winv(i)
         winvj = winv(j)
         if(belly) then
            if(igrp(i) <= 0) winvi = zero
            if(igrp(j) <= 0) winvj = zero
         end if
         
         ! -- dot(x_{ij}, v_{ij}) = dcdt/2

         rvdot = zero
         do m = 1, 3
           xij(m) = x(i3+m)-x(j3+m)
           vij(m) = v(i3+m)-v(j3+m)
           rvdot = rvdot + xij(m)*vij(m)
         end do
         diff = -rvdot

         ! -- check if constraint for velocity is satisfied
         acor = diff/(toler*(winv(i)+winv(j)))
         if (abs(acor) < tol2) cycle bonds
         do m = 1,3
            xh = xij(m)*acor
            v(i3+m) = v(i3+m)+xh*winvi
            v(j3+m) = v(j3+m)-xh*winvj
         end do
         skip(i) = .false.
         skip(j) = .false.
         ready = .false.
      end do bonds
      
      nit = nit+1
      do k = 1,nrr
         skip(nrr+k) = skip(k)
         skip(k) = .true.
      end do
      
   end do  !  ( .not. ready )

   !          END OF BIG ITERATION LOOP
   !==================================================================

   niter = niter+nit
   return

   311 format(/5x,'Velocity resetting (RATTLEV) was not accomplished', &
         /5x,'within 3000 iterations')
   321 format(/5x,'Velocity resetting (RATTLEV) cannot be accomplished,', &
         /5x,'deviation is too large', &
         /5x,'NITER, NIT, LL, I and J are :',5i7)
   322 format(/5x,'Note: This is usually a symptom of some deeper', &
         /5x,'problem with the energetics of the system.' )
end subroutine
