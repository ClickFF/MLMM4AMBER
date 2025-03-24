! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "nfe-config.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write restart data to file 'restrt' and possibly a numbered restrt
subroutine mdwrit(nstep,nr,ntxo,ntb,x,v,tt,temp0, solvph, solve)

   use file_io_dat
#if !defined(DISABLE_NFE)
   use nfe_sander_hooks, only : nfe_on_mdwrit => on_mdwrit
   use nfe_sander_proxy, only : infe
#endif /* DISABLE_NFE */
   use binrestart, only: write_nc_restart
#ifdef MPI
   use remd, only : rem, remd_dimension, remd_types, group_num, &
                    replica_indexes, stagid, remd_repidx, remd_crdidx
   use sgld, only : trxsgld
#endif

   implicit none
   integer nstep,nr,ntxo,ntb
   _REAL_ x(*),v(*),tt,temp0, solvph, solve, temp0_or_stagid
   character(len=89) restrt2
   character(len=12) num
   integer istart,iend
   logical first
   save first
   data first/.true./
   
   !     -- open/write/close the restrt file:

#if !defined(DISABLE_NFE)
   if ( infe == 1 ) then
     call nfe_on_mdwrit()
   end if
#endif /* DISABLE_NFE */

   ! For a NetCDF restart a RXSGLD stage ID can be labeled as such.
   ! However, since the non NetCDF file formats are fixed, for RXSGLD store
   ! the stage ID in the temp0 slot.  There is potential for confusion,
   ! but temperatures and stage IDs are likely to be numerically distinct.
   temp0_or_stagid = temp0
#ifdef MPI
   if (rem == 1) then
      temp0_or_stagid = temp0
   else if (rem == 4) then
      temp0_or_stagid = solvph
   else if (rem == 5) then
      temp0_or_stagid = solve
   end if
   if (trxsgld) then
      temp0_or_stagid = stagid
   end if
#endif

   if ( ntxo == 2) then
   ! Netcdf restart
#     ifdef MPI
      call write_nc_restart(restrt, title,owrite,nr,ntb,first, x,v, tt, .true. &
               , temp0, rem, remd_dimension, remd_types, group_num &
               , replica_indexes, stagid, remd_repidx, remd_crdidx &
               , solvph, solve)
#     else
      call write_nc_restart(restrt, title,owrite,nr,ntb,first, x,v, tt, .true.)
#     endif
      first = .false.
   else
   ! Standard formatted restart   
      if( first ) then
         call amopen(16,restrt,owrite,'F','W')
         first = .false.
      else
         call amopen(16,restrt,'O','F','W')
      end if
      call mdwri2(16,nr,ntb,x,v,tt,temp0_or_stagid)
      close(16)
   endif
   
   if (ntwr >= 0) return

   !     -- write unique restart files: restrt_<nstep> 
   
   do iend=1,80
      if (restrt(iend:iend) <= ' ') goto 1
   end do
   1 continue
   iend = iend - 1
   write(num,'(i12)') nstep
   do istart=1,12
      if (num(istart:istart) /= ' ') goto 2
   end do
   2 continue
   write(restrt2, '(a,a,a)') restrt(1:iend), '_', num(istart:12)
   write(6,'(a,a)') ' writing ', restrt2
   if (ntxo == 2) then
#     ifdef MPI
      call write_nc_restart(restrt2,title,owrite,nr,ntb,.true.,x,v, tt, .true. &
               , temp0, rem, remd_dimension, remd_types, group_num &
               , replica_indexes, stagid, remd_repidx, remd_crdidx &
               , solvph, solve)
#     else
      call write_nc_restart(restrt2,title,owrite,nr,ntb,.true.,x,v, tt, .true.)
#     endif
      return
   else
      call amopen(17,restrt2,owrite,'F','W')
   end if
   call mdwri2(17,nr,ntb,x,v,tt,temp0_or_stagid)
   close(17)
   return

   entry mdwrit_reset
   first = .true.

end subroutine mdwrit

!------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write dynamics-restart data to an arbitrary filehandle
subroutine mdwri2(nf,nr,ntb,x,v,tt,temp0)
   use nblist, only: a,b,c,alpha,beta,gamma
   use file_io_dat
#ifdef MPI
   use remd, only : rem
#endif
   implicit none
   integer nf,nr,ntb
#ifndef MPI
   integer, parameter :: rem = 0
#endif
   
   !     ----- ROUTINE TO WRITE FINAL COORDINATES AND VELOCITIES -----
   
   _REAL_ x(*),v(*),tt,temp0
   integer nr3,i
   
   nr3 = 3*nr
   !     ----- FORMATTED WRITING -----
      
   write(nf,9008) title

   if(rem > 0) then
      if( nr < 100000 ) then
        write(nf,9018) nr,tt,temp0
      elseif ( nr < 1000000 ) then
        write(nf,9019) nr,tt,temp0 ! sander 7/8/9/10 large system format...
      elseif ( nr < 10000000 ) then
        write(nf,9020) nr,tt,temp0 ! Sander 11 - 1 mil+ format
      else
        write(nf,9021) nr,tt,temp0 ! assume amber 11 VERY large system format. 10 mil+
      end if
   else
      if( nr < 100000 ) then
        write(nf,9018) nr,tt
      elseif ( nr < 1000000 ) then ! sander 7/8/9/10 large system format...
        write(nf,9019) nr,tt
      elseif ( nr < 10000000 ) then ! Sander 11 - 1 mil+ format
        write(nf,9020) nr,tt
      else
        write(nf,9021) nr,tt ! assume amber 11 VERY large system format. 10 mil+
       end if
   end if
   write(nf,9028) (x(i),i=1,nr3)
   write(nf,9028) (v(i),i=1,nr3)

   if ( ntb /= 0 ) write(nf,9028) a,b,c,alpha,beta,gamma
   
   9008 format(a80)
   9018 format(i5,2e15.7)
   9019 format(i6,2e15.7)
   9020 format(i7,2e15.7)
   9021 format(i8,2e15.7)
   9028 format(6f12.7)
   return
end subroutine mdwri2 
