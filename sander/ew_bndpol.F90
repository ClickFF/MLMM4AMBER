! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

module bndpol

   integer, parameter :: MAXBND = 8

   integer, allocatable :: nbnd(:)
   integer, allocatable :: ibnd(:,:)
   _REAL_, allocatable :: pbnd(:,:)

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ read and set up working arrays for bond polarizabilities
subroutine ew_bndpol( natom,nbonh,nbona,ibh,jbh,iba,jba,igraph,isymbl )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Ray Luo, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! Passed variables

   integer natom, nbonh, nbona
   integer ibh(*), jbh(*), iba(*), jba(*)
   character (len=4) :: igraph(*), isymbl(*)

   integer iatm, jatm, idum, ntmp
   _REAL_, allocatable :: bona_pol(:), bonh_pol(:)
   
   ! temperary read statements for bonh_pol and bona_pol arrays
   ! note these two arrays should be formatted in the same orders of bonh and
   ! bona arrays in prmtop 

   allocate(bonh_pol(1:nbonh))
   allocate(bona_pol(1:nbona))

   read(100, *) ntmp
   if ( ntmp /= nbonh ) then
      write(6,*) 'FATAL in EW_BNDPOL: mismatch in nbonh'
      call mexit(6,1)
   end if
   read(100, *) bonh_pol(1:ntmp)

   read(100, *) ntmp
   if ( ntmp /= nbona ) then
      write(6,*) 'FATAL in EW_BNDPOL: mismatch in nbona'
      call mexit(6,1)
   end if
   read(100, *) bona_pol(1:ntmp)

   ! these are the arrays to be used later

   allocate(nbnd(natom))
   allocate(ibnd(MAXBND,natom))
   allocate(pbnd(MAXBND,natom))

   ! accumulate no. of bonds attached to each atom

   nbnd(1:natom) = 0

   ! h-containing bonds

   do idum = 1, nbonh
      iatm = ibh(idum)/3 + 1
      jatm = jbh(idum)/3 + 1

      ! find a new bond for each bonded atom

      nbnd(iatm) = nbnd(iatm) + 1
      nbnd(jatm) = nbnd(jatm) + 1
     
      ! save the atom no. for this bond 

      ibnd(nbnd(iatm),iatm) = jatm
      ibnd(nbnd(jatm),jatm) = iatm

      ! save the bond pol for this bond

      pbnd(nbnd(iatm),iatm) = bonh_pol(idum)
      pbnd(nbnd(jatm),jatm) = bonh_pol(idum)

   end do

   ! other bonds

   do idum = 1, nbona
      iatm = iba(idum)/3 + 1
      jatm = jba(idum)/3 + 1

      ! find a new bond for each bonded atom

      nbnd(iatm) = nbnd(iatm) + 1
      nbnd(jatm) = nbnd(jatm) + 1
      if ( nbnd(iatm) > MAXBND .or. nbnd(jatm) > MAXBND ) then
         write(6,*) 'FATAL in EW_BNDPOL: MAXBND breached at bond I/J', iatm, jatm
         call mexit(6,1)
      end if
     
      ! save the atom no. for this bond 

      ibnd(nbnd(iatm),iatm) = jatm
      ibnd(nbnd(jatm),jatm) = iatm

      ! save the bond pol for this bond

      pbnd(nbnd(iatm),iatm) = bona_pol(idum)
      pbnd(nbnd(jatm),jatm) = bona_pol(idum)
   end do

   ! debugging setup arrays

   do iatm = 1, natom
      write(200, '(a12,i6,3a12)') 'BND POL:', iatm, isymbl(iatm), igraph(iatm), '::::'
      do idum = 1, nbnd(iatm)
         write(200, '(a12,i6,2a12,f15.6)') '  bonded to', ibnd(idum,iatm), &
         isymbl(ibnd(idum,iatm)), igraph(ibnd(idum,iatm)), pbnd(idum,iatm)
      end do
   end do

end subroutine ew_bndpol 

end module bndpol
