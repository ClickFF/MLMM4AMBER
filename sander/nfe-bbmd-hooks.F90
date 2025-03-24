#include "nfe-utils.h"
#include "nfe-config.h"

!---------------------!
! Bizarrely Biased MD !
!---------------------!

!------------------------------------------------------------------
! Updated in Sep. 18, since the compatible H-REMD is implemented by
! F. Pan (with -rem 3), this module seems redundant, but we still
! keep it for records.
!------------------------------------------------------------------

module nfe_bbmd_hooks

#ifdef NFE_ENABLE_BBMD

use mt19937, only : mt19937_t
use nfe_bbmd_ctxt, only : bbmd_ctxt_t
use nfe_constants, only : SL => STRING_LENGTH, LOG_UNIT => BBMD_LOG_UNIT

implicit none

public :: on_sander_init
public :: on_sander_exit

public :: on_force
public :: on_mdstep
public :: on_mdwrit

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

character(*), private, parameter :: SECTION = 'bbmd'

integer, private, save :: active = 0 ! integer for MPI

type(bbmd_ctxt_t), private, allocatable, save :: contexts(:)

integer, private, save :: mdstep
integer, private, save :: exchno

integer, private, save :: current
integer, private, save :: partner
integer, private, save :: partner_masterrank

NFE_REAL, private, save :: U_mm, U_mo, U_om, U_oo

type(mt19937_t), private, save :: mersenne_twister
integer, private, allocatable, save :: permutation(:)

! symmetric matrices with 0s on the diagonal
! [row major : (2,1), (3,1), (3,2), (4,1), ... ]
integer, private, allocatable, save :: local_exchg_attempts(:)
integer, private, allocatable, save :: local_exchg_accepted(:)

! on masterrank.eq.0 (layout as above)
integer, private, allocatable, save :: exchg_attempts(:)
integer, private, allocatable, save :: exchg_accepted(:)

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

contains

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

subroutine on_sander_init(mdin_unit, amass)

   use mt19937, only : &
      init_by_seed, mt19937_load, mt19937_bcast

   use nfe_utils
   use nfe_constants
   use nfe_sander_proxy
   use file_io_dat

   use nfe_bbmd_ctxt

   implicit none

   integer, intent(in) :: mdin_unit

   NFE_REAL, intent(in) :: amass(*)

#  include "nfe-mpi.h"

   integer :: n, error, mt19937_loaded, ifind
   integer, allocatable :: active_all(:)

#  ifndef NFE_DISABLE_ASSERT
   logical, save :: called = .false.
#  endif /* NFE_DISABLE_ASSERT */

   logical :: found

   nfe_assert(active.eq.0)
   nfe_assert(multisander_rem().eq.0)
   nfe_assert(multisander_numgroup().gt.1)

   nfe_assert(.not.called)
   nfe_assert(.not.allocated(contexts))

#  ifndef NFE_DISABLE_ASSERT
   called = .true.
#  endif /* NFE_DISABLE_ASSERT */
   
   if (sanderrank.eq.0) then
         rewind(mdin_unit)
         call nmlsrc('bbmd', mdin_unit, ifind)
      if (ifind.ne.0) then
         active = 1
      else
         active = 0
      end if
   end if ! sanderrank.eq.0

   if (sanderrank.eq.0) then
      nfe_assert(commmaster.ne.MPI_COMM_NULL)

      allocate(active_all(mastersize), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY

      call mpi_gather(active, 1, MPI_INTEGER, active_all, &
         1, MPI_INTEGER, 0, commmaster, error)
      nfe_assert(error.eq.0)

      if (masterrank.eq.0) then
         do n = 1, mastersize
            if (active.ne.active_all(n)) &
               call fatal('either all or none of MDIN files &
                          &should contain the '''//SECTION//''' section')
         end do
      end if ! masterrank.eq.0

      deallocate(active_all)
   end if ! sanderrank.eq.0

   call mpi_bcast(active, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   if (active.eq.0) &
      return

   ! <<>> BBMD is requested in all replicas <<>>

   call mpi_bcast(mastersize, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(masterrank, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   allocate(contexts(mastersize), &
      permutation(mastersize + mod(mastersize, 2)), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   call ctxt_init(contexts(masterrank + 1), amass)

   do n = 1, mastersize
      call ctxt_bcast(contexts(n), n - 1, amass)
   end do

   call mpi_barrier(commworld, error) ! unneeded actually
   nfe_assert(error.eq.0)

   ! exchange frequency, logfile, etc

   if (sanderrank.eq.0) then
      exchange_freq = 500
      mt19937_file = 'nfe-mt19937.nc'
      mt19937_seed = 5489
      exchange_log_file = 'nfe-bbmd-exchange-log'
      exchange_log_freq = 100
      
      rewind(mdin_unit)
      read(mdin_unit,nml=bbmd,err=666)

      ! try load the twister first
      mt19937_loaded = 0

#     ifndef NFE_NO_NETCDF
      if (masterrank.eq.0) then
         inquire (file = mt19937_file, exist = found)
         if (found) then
            call mt19937_load(mersenne_twister, mt19937_file)
            mt19937_loaded = 1
         end if ! found
      end if ! masterrank.eq.0
#     else
      write (unit = ERR_UNIT, fmt = '(a,a)') NFE_WARNING, &
         'mt19937 : netCDF is not available (try ''-bintraj'' configure option)'
#     endif /* NFE_NO_NETCDF */

      if (mt19937_loaded.eq.0) then
         call init_by_seed(mersenne_twister, mt19937_seed)
      end if ! .not.mt19937_loaded

      call mt19937_bcast(mersenne_twister, commmaster, 0)

      call mpi_bcast(mt19937_file, len(mt19937_file), &
                     MPI_CHARACTER, 0, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(mt19937_seed, 1, MPI_INTEGER, 0, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(mt19937_loaded, 1, MPI_INTEGER, 0, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(exchange_freq, 1, MPI_INTEGER, 0, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(exchange_log_file, len(exchange_log_file), &
                     MPI_CHARACTER, 0, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(exchange_log_freq, 1, MPI_INTEGER, 0, commmaster, error)
      nfe_assert(error.eq.0)

      allocate(local_exchg_attempts(mastersize*(mastersize - 1)/2), &
         local_exchg_accepted(mastersize*(mastersize - 1)/2), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY

      local_exchg_attempts(:) = 0
      local_exchg_accepted(:) = 0

      if (masterrank.eq.0) then
         allocate(exchg_attempts(mastersize*(mastersize - 1)/2), &
            exchg_accepted(mastersize*(mastersize - 1)/2), stat = error)
         if (error.ne.0) &
            NFE_OUT_OF_MEMORY

         open (unit = LOG_UNIT, file = exchange_log_file, iostat = error, &
               form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')
         if (error.ne.0) then
            write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
               NFE_ERROR, 'failed to open ''', &
               trim(exchange_log_file), ''' for writing'
            call terminate()
         end if
      end if ! masterrank.eq.0
   end if ! sanderrank.eq.0

   mdstep = 0
   exchno = 0

   current = masterrank + 1
   partner = -1

   call mpi_bcast(exchange_freq, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   if (sanderrank.gt.0) &
      return

   write (unit = OUT_UNIT, fmt = '(/a,a)') NFE_INFO, &
      '/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/ B. B. M. D. /=/=/=/=/=/=/=/=/=/=/=/=/=/'
   write (unit = OUT_UNIT, &
      fmt = '(a,/a,a,'//pfmt(current)//',a,'//pfmt(mastersize)//')') &
      NFE_INFO, NFE_INFO, &
      'this is replica ', current, ' of ', mastersize
   write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(exchange_freq)//')') &
      NFE_INFO, 'number of MD steps between exchange attempts = ', exchange_freq
#  ifdef NFE_NO_NETCDF
   write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(mt19937_seed)//',a)') &
      NFE_INFO, 'mersenne twister seed = ', mt19937_seed, ' (no netCDF)'
#  else
   write (unit = OUT_UNIT, fmt = '(a,a,a,a)') NFE_INFO, &
      'mersenne twister file = ''', trim(mt19937_file), ''''
   if (mt19937_loaded.eq.1) then
      write (unit = OUT_UNIT, fmt = '(a,a)') NFE_INFO, &
         'mersenne twister state loaded'
   else
      write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(mt19937_seed)//')') &
         NFE_INFO, 'mersenne twister seed = ', mt19937_seed
   end if
#  endif /* NFE_NO_NETCDF */
   write (unit = OUT_UNIT, fmt = '(a,a,a,a)') NFE_INFO, &
      'exchange log file = ''', trim(exchange_log_file), ''''
   write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(exchange_log_freq)//',/a)') &
      NFE_INFO, 'exchange log update frequency = ', &
      exchange_log_freq, NFE_INFO

   call ctxt_print(contexts(masterrank + 1), OUT_UNIT)

   write (unit = OUT_UNIT, fmt = '(a,/a,a/)') NFE_INFO, NFE_INFO, &
      '/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/'
      
   return
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &bbmd namelist!'
    call terminate()      

end subroutine on_sander_init

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

subroutine on_sander_exit()

   use nfe_utils
   use nfe_bbmd_ctxt

   implicit none

#  include "nfe-mpi.h"

   integer :: n

   if (active.eq.0) then
      nfe_assert(.not.allocated(contexts))
      nfe_assert(.not.allocated(permutation))
      return
   end if ! active.eq.0

   nfe_assert(mastersize.gt.1)

   do n = 1, mastersize
      call ctxt_fini(contexts(n))
   end do

   deallocate(contexts, permutation)
   if (sanderrank.eq.0) then
      deallocate(local_exchg_attempts, local_exchg_accepted)

      if (masterrank.eq.0) then
         deallocate(exchg_attempts, exchg_accepted)
         close (LOG_UNIT)
      end if ! masterrank.eq.0
   end if ! sanderrank.eq.0

   active = 0

end subroutine on_sander_exit

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

! Modified by M Moradi
! Driven ABMD
subroutine on_force(x, f, work, udr, pot)
! Moradi end

   NFE_USE_AFAILED

   use mt19937
   use nfe_bbmd_ctxt

   implicit none

   NFE_REAL, intent(in) :: x(*)

   NFE_REAL, intent(inout) :: f(*)
   
! Modified by M Moradi
! Driven ABMD
   NFE_REAL, intent(in) :: work
   NFE_REAL, intent(in) :: udr
! Moradi end

   NFE_REAL, intent(inout) :: pot
   
#  include "nfe-mpi.h"

   integer :: n, k, error

   if (active.eq.0) &
      return

   nfe_assert(mastersize.gt.1)

   nfe_assert(current.gt.0)
   nfe_assert(current.le.mastersize)
   nfe_assert(allocated(contexts))
   nfe_assert(allocated(permutation))

   call ctxt_on_force(contexts(current), x, f, mdstep, work, udr, pot)

   nfe_assert(partner.lt.0)

   ! prepare for an exchange attempt if needed
   if (mod(mdstep + 1, exchange_freq).ne.0) &
      return

   ! generate a random permutation
   if (sanderrank.eq.0) then
      do n = 1, mastersize + mod(mastersize, 2)
         permutation(n) = n
      end do
      do n = 1, mastersize + mod(mastersize, 2) - 1
         k = 1 + n + mod(random_int31(mersenne_twister), &
            mastersize + mod(mastersize, 2) - n)
         error = permutation(k)
         permutation(k) = permutation(n)
         permutation(n) = error
      end do
      partner_masterrank = -1
      do n = 1, mastersize + mod(mastersize + 1, 2), 2
         if (permutation(n).le.mastersize.and. &
            permutation(n + 1).le.mastersize) then
            if (permutation(n).eq.(masterrank + 1)) then
               partner_masterrank = permutation(n + 1) - 1
               exit
            end if
            if (permutation(n + 1).eq.(masterrank + 1)) then
               partner_masterrank  = permutation(n) - 1
               exit
            end if
         end if
      end do
   end if ! sanderrank.eq.0

   call mpi_bcast(partner_masterrank, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   nfe_assert(masterrank.ne.partner_masterrank)

   if (partner_masterrank.lt.0) & ! for mod(mastersize, 2).ne.0
      return

   ! get the partner's context number
   if (sanderrank.eq.0) then
      call mpi_sendrecv(current, 1, MPI_INTEGER, partner_masterrank, mdstep, &
                        partner, 1, MPI_INTEGER, partner_masterrank, mdstep, &
                        commmaster, MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
      nfe_assert(partner.ge.1.and.partner.le.mastersize)
   end if ! sanderrank.eq.0

   call mpi_bcast(partner, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   ! U_?? are needed only in the replica with smaller masterrank
   if (masterrank.lt.partner_masterrank) then
      call ctxt_Um(contexts(current), partner_masterrank, x, U_mm, U_mo)
      call ctxt_Uo(contexts(partner), partner_masterrank, x, U_om, U_oo)
   else
      call ctxt_Um(contexts(partner), partner_masterrank, x, U_mm, U_mo)
      call ctxt_Uo(contexts(current), partner_masterrank, x, U_om, U_oo)
   end if ! masterrank.lt.partner_masterrank

end subroutine on_force

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

subroutine on_mdstep(E_p, v, ekmh)

   use mt19937
   use nfe_utils
   use nfe_constants
   use nfe_bbmd_ctxt
   use nfe_sander_proxy

   implicit none

   NFE_REAL, intent(in) :: E_p
   NFE_REAL, intent(inout) :: v(*)
   NFE_REAL, intent(inout) :: ekmh

   NFE_REAL, parameter :: TINY = 0.00001d0

#  include "nfe-mpi.h"

   NFE_REAL :: v_scale
   NFE_REAL :: E_px, partner_E_px, temp0, partner_temp0
   NFE_REAL :: beta_m, beta_o, delta, random_number

   integer :: error, exchange, n

   if (active.eq.0) &
      return

   nfe_assert(mastersize.gt.1)

   nfe_assert(current.gt.0)
   nfe_assert(current.le.mastersize)
   nfe_assert(allocated(contexts))

   mdstep = mdstep + 1

   ! added by F Pan, since we already include the potential contribution
   ! from bbmd in E_p, we need to deduct it
   E_px = E_p - nfe_pot_ene%bbmd

   if (mod(mdstep, exchange_freq).ne.0) then
      if (mdstep.eq.sander_nstlim()) &
         goto 2 ! update_log()
      return
   end if ! mod(mdstep, exchange_freq).ne.0

   random_number = random_res53(mersenne_twister)

   exchno = exchno + 1

   if (partner_masterrank.lt.0) & ! this replica does not participate
      goto 2 ! update_log()

   nfe_assert(masterrank.ne.partner_masterrank)
   nfe_assert(partner.ge.1.and.partner.le.mastersize)
   nfe_assert(partner.ne.current)

   temp0 = sander_temp0()

   if (sanderrank.ne.0) &
      goto 1 ! bcast(exchange)

   ! get partner's E_p
   if (masterrank.lt.partner_masterrank) then
      call mpi_recv(partner_E_px, 1, MPI_DOUBLE_PRECISION, &
                    partner_masterrank, mdstep, commmaster, &
                    MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
   else
      call mpi_send(E_px, 1, MPI_DOUBLE_PRECISION, &
                    partner_masterrank, mdstep, commmaster, error)
      nfe_assert(error.eq.0)
   end if

   call mpi_sendrecv &
      (temp0, 1, MPI_DOUBLE_PRECISION, partner_masterrank, mdstep, &
       partner_temp0, 1, MPI_DOUBLE_PRECISION, partner_masterrank, mdstep, &
       commmaster, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)

   if (masterrank.lt.partner_masterrank) then

      nfe_assert(temp0.gt.ZERO)
      nfe_assert(partner_temp0.gt.ZERO)

      beta_m = 503.01D0/temp0
      beta_o = 503.01D0/partner_temp0

      delta = (beta_o - beta_m)*(E_px - partner_E_px) &
         + beta_m*(U_mo - U_mm) - beta_o*(U_oo - U_om)

      if (delta.lt.ZERO) then
         exchange = 1
      else if (exp(-delta).gt.random_number) then
         exchange = 1
      else
         exchange = 0
      end if

      call mpi_send(exchange, 1, MPI_INTEGER, partner_masterrank, &
                    mdstep, commmaster, error)
      nfe_assert(error.eq.0)
   else
      call mpi_recv(exchange, 1, MPI_INTEGER, partner_masterrank, &
                    mdstep, commmaster, MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
   end if ! masterrank.lt.partner_masterrank

   if (current.lt.partner) then
      error = (partner - 1)*(partner - 2)/2 + current
      local_exchg_attempts(error) = local_exchg_attempts(error) + 1
      if (exchange.eq.1) &
         local_exchg_accepted(error) = local_exchg_accepted(error) + 1
   end if ! current.lt.partner

1  call mpi_bcast(exchange, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   if (exchange.eq.1) then

      if (sanderrank.eq.0) then
         write (unit = OUT_UNIT, &
         fmt ='(a,a,'//pfmt(partner)//',a,'//pfmt(sander_mdtime(), 3)//')') &
            NFE_INFO, 'BBMD : exchanged coordinates/velocities with ', &
            partner, ' at t = ', sander_mdtime()
         call ctxt_close_units(contexts(current))
         if (masterrank.lt.partner_masterrank) then
            call ctxt_recv(contexts(partner), partner_masterrank)
            call ctxt_send(contexts(current), partner_masterrank)
         else
            call ctxt_send(contexts(current), partner_masterrank)
            call ctxt_recv(contexts(partner), partner_masterrank)
         end if ! masterrank.lt.partner_masterrank
         call ctxt_open_units(contexts(partner))
      end if ! sanderrank.eq.0

      current = partner

      call mpi_bcast(partner_temp0, 1, MPI_DOUBLE_PRECISION, &
                     0, commsander, error)
      nfe_assert(error.eq.0)

      if (abs(temp0 - partner_temp0).gt.TINY) then
         v_scale = sqrt(partner_temp0/temp0)
         do n = 1, 3*sander_natoms()
            v(n) = v_scale*v(n)
         end do
         ekmh = ekmh*v_scale*v_scale
      end if

      call set_sander_temp0(partner_temp0)
   end if ! exchange.eq.1

2  if (sanderrank.eq.0) then
      if (mod(exchno, exchange_log_freq).eq.0.or.mdstep.eq.sander_nstlim()) &
         call update_log()
   end if ! sanderrank.eq.0

#  ifndef NFE_DISABLE_ASSERT
   partner = -1
#  endif /* NFE_DISABLE_ASSERT */

contains

subroutine update_log()

   implicit none

   integer :: n1, n2, idx
   character(32) :: log_fmt

   nfe_assert(sanderrank.eq.0)

   if (masterrank.eq.0) then
      call mpi_reduce(local_exchg_attempts, exchg_attempts, &
         mastersize*(mastersize - 1)/2, MPI_INTEGER, MPI_SUM, &
         0, commmaster, error)
   else
      call mpi_reduce(local_exchg_attempts, 0, &
         mastersize*(mastersize - 1)/2, MPI_INTEGER, MPI_SUM, &
         0, commmaster, error)
   end if
   nfe_assert(error.eq.0)

   if (masterrank.eq.0) then
      call mpi_reduce(local_exchg_accepted, exchg_accepted, &
         mastersize*(mastersize - 1)/2, MPI_INTEGER, MPI_SUM, &
         0, commmaster, error)
   else
      call mpi_reduce(local_exchg_accepted, 0, &
         mastersize*(mastersize - 1)/2, MPI_INTEGER, MPI_SUM, &
         0, commmaster, error)
   end if
   nfe_assert(error.eq.0)

   if (masterrank.ne.0) &
      return

   nfe_assert(allocated(exchg_attempts))
   nfe_assert(allocated(exchg_accepted))

   if (exchno.gt.exchange_log_freq) &
      write (unit = LOG_UNIT, fmt = '(80(''=''))')

   write (unit = LOG_UNIT, &
      fmt = '(/a,'//pfmt((mastersize - mod(mastersize, 2)/2)/2)//',a,'//pfmt &
      (exchno)//',a/)') ' <> exchange statistics over ', &
      ((mastersize - mod(mastersize, 2)/2)/2),' x ', exchno, ' attempts <>'

   idx = 2 + int(floor(log10(ONE + NFE_TO_REAL(maxval(exchg_attempts)))))
   write (unit = log_fmt, fmt = '(a,'//pfmt(idx)//',a)') '(i', idx, ')'

   write (unit = LOG_UNIT, fmt = '(a/)') ' * number of trials:'
   write (unit = LOG_UNIT, fmt = '(6x)', advance = 'NO')

   do n1 = 1, mastersize - 1
      write (unit = LOG_UNIT, fmt = log_fmt, advance = 'NO') n1
   end do

   do n1 = 2, mastersize
      write (unit = LOG_UNIT, fmt = '(/i3,1x,'':'',1x)', advance = 'NO') n1
      do n2 = 1, n1 - 1
         idx = (n1 - 1)*(n1 - 2)/2 + n2
         write (unit = LOG_UNIT, fmt = log_fmt, advance = 'NO') &
            exchg_attempts(idx)
      end do
   end do

   write (unit = LOG_UNIT, fmt = '(//a/)') ' * number of exchanges:'
   write (unit = LOG_UNIT, fmt = '(6x)', advance = 'NO')

   do n1 = 1, mastersize - 1
      write (unit = LOG_UNIT, fmt = log_fmt, advance = 'NO') n1
   end do

   do n1 = 2, mastersize
      write (unit = LOG_UNIT, fmt = '(/i3,1x,'':'',1x)', advance = 'NO') n1
      do n2 = 1, n1 - 1
         idx = (n1 - 1)*(n1 - 2)/2 + n2
         write (unit = LOG_UNIT, fmt = log_fmt, advance = 'NO') &
            exchg_accepted(idx)
      end do
   end do

   write (unit = LOG_UNIT, fmt = '(//a/)') ' * acceptance rates:'
   write (unit = LOG_UNIT, fmt = '(6x)', advance = 'NO')

   do n1 = 1, mastersize - 1
      write (unit = LOG_UNIT, fmt = '(i6)', advance = 'NO') n1
   end do

   do n1 = 2, mastersize
      write (unit = LOG_UNIT, fmt = '(/i3,1x,'':'',1x)', advance = 'NO') n1
      do n2 = 1, n1 - 1
         idx = (n1 - 1)*(n1 - 2)/2 + n2
         if (exchg_attempts(idx).gt.0) then
            write (unit = LOG_UNIT, fmt = '(f6.1)', advance = 'NO') &
              ((NFE_TO_REAL(100)*exchg_accepted(idx))/exchg_attempts(idx))
         else
            write (unit = LOG_UNIT, fmt = '(a)', advance = 'NO') '  XX.X'
         end if ! exchg_attempts(idx).gt.0
      end do
   end do

   write (unit = LOG_UNIT, fmt = '(/a)') ''
   call flush_UNIT(LOG_UNIT)

end subroutine update_log

end subroutine on_mdstep

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

subroutine on_mdwrit()

   NFE_USE_AFAILED
   use nfe_bbmd_ctxt
#  ifndef NFE_NO_NETCDF
   use mt19937, only : mt19937_save
#  endif /* NFE_NO_NETCDF */

   implicit none

#  include "nfe-mpi.h"

   if (active.eq.0) &
      return

   nfe_assert(mastersize.gt.1)

   nfe_assert(current.gt.0)
   nfe_assert(current.le.mastersize)
   nfe_assert(allocated(contexts))

   call ctxt_on_mdwrit(contexts(current))

#  ifndef NFE_NO_NETCDF
   if (sanderrank.eq.0.and.masterrank.eq.0) &
      call mt19937_save(mersenne_twister, mt19937_file)
#  endif /* NFE_NO_NETCDF */

end subroutine on_mdwrit

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

#endif /* NFE_ENABLE_BBMD */

end module nfe_bbmd_hooks

