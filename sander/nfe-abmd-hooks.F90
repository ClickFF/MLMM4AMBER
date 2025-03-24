!<compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

module nfe_abmd_hooks

use nfe_constants, only : SL => STRING_LENGTH, ABMD_MONITOR_UNIT, ABMD_CV_UNIT

use nfe_umbrella, only : umbrella_t, &
   MAX_NUMBER_OF_COLVARS => UMBRELLA_MAX_NEXTENTS

use nfe_colvar_type, only : colvar_t

use random

implicit none

private

!-----------------------------------------------------------------------------
!                          T H E    H O O K S
!-----------------------------------------------------------------------------

! old REMD subroutines, not used now --------
#ifdef MPI
public :: on_delta
public :: on_exchange
#endif /* MPI */
! -------------------------------------------
public :: on_multisander_exit

!-----------------------------------------------------------------------------

! sander.f
public :: on_sander_init
public :: on_sander_exit

! force.f
public :: on_force

! Modified by M Moradi
! Selection algorithm
! runmd.f
#ifdef NFE_ENABLE_BBMD
public :: on_mdstep
#endif /* NFE_ENABLE_BBMD */
! Moradi end

! mdwrit.f
public :: on_mdwrit

!-----------------------------------------------------------------------------
!                             * P R I V A T E *
!-----------------------------------------------------------------------------

integer, private, parameter :: MONITOR_UNIT = ABMD_MONITOR_UNIT
integer, private, parameter :: CV_UNIT = ABMD_CV_UNIT

character(*), private, parameter :: &
   DEFAULT_MONITOR_FILE = 'nfe-abmd-monitor', &
   DEFAULT_UMBRELLA_FILE = 'nfe-abmd-umbrella', &
   DEFAULT_WT_UMBRELLA_FILE = 'nfe-abmd-wt-umbrella', &
   DEFAULT_SNAPSHOTS_BASENAME = 'nfe-abmd-umbrella-snapshot', &
   DEFAULT_CV_FILE = 'nfe-abmd-cv'


integer, private, parameter :: MODE_NONE = 0

integer, private, parameter :: MODE_ANALYSIS = 123
integer, private, parameter :: MODE_UMBRELLA = 234
integer, private, parameter :: MODE_FLOODING = 345

integer, private, save :: imode = MODE_NONE
character(len = SL), private,save :: mode = 'NONE'

!-----------------------------------------------------------------------------
integer, private, save :: ncolvars = 0

type(colvar_t),   private, save :: colvars(MAX_NUMBER_OF_COLVARS)
type(umbrella_t), private, save :: umbrella ! master only
type(umbrella_t), private, save :: wt_umbrella
!-----------------------------------------------------------------------------
#ifndef MPI
NFE_REAL, private, save :: instantaneous(MAX_NUMBER_OF_COLVARS) ! master only
#else
NFE_REAL, private, save :: instantaneous(MAX_NUMBER_OF_COLVARS + 1)
#endif /* MPI */
!-----------------------------------------------------------------------------
character(len = SL), private, save :: monitor_file ! master only
character(len = SL), private, save :: monitor_fmt ! master only
integer,             private, save :: monitor_freq = 50
character(len = SL), private, save :: cv_file = DEFAULT_CV_FILE
!-----------------------------------------------------------------------------
character(len = SL), private, save :: snapshots_basename ! master only
integer,             private, save :: snapshots_freq = -1 ! master only
!-----------------------------------------------------------------------------
character(len = SL), private, save :: umbrella_file ! master only
character(len = SL), private, save :: wt_umbrella_file ! master only
!-----------------------------------------------------------------------------
NFE_REAL, private, save :: timescale = 1 ! master only
!-----------------------------------------------------------------------------
integer, private, save :: mdstep ! = runmd.f::nstep + 1
!-----------------------------------------------------------------------------
integer, private, save :: nhills ! not zeroed on exchange
!-----------------------------------------------------------------------------
! Modified by M Moradi
! Selection algorithm
integer,            private, save :: selection_freq = 0
NFE_REAL,           private, save :: selection_constant = 1.0
NFE_REAL,           private, save :: selection_epsilon = 0.0
type(rand_gen_state),        save :: nfe_abmd_gen
! Well-tempered ABMD
NFE_REAL,           private, save :: pseudo = 0.0
NFE_REAL,           private, save :: wt_temperature = 0.0
! Driven ABMD
integer,            private, save :: drivenw
integer,            private, save :: drivenu
NFE_REAL,           private, save :: driven_cutoff = 0.0
character(len = SL), private, save :: driven_weight = 'NONE'
! Moradi end
!-----------------------------------------------------------------------------
#ifdef MPI

!
! multiwalk/rem
!

NFE_REAL, private, allocatable, save :: all_hills(:)

! Modified by M Moradi
! Selection algorithm
#ifdef NFE_ENABLE_BBMD
NFE_REAL, private, allocatable, save :: all_x(:)
NFE_REAL, private, allocatable, save :: all_v(:)
NFE_REAL, private, allocatable, save :: all_w(:)
NFE_REAL, private, allocatable, save :: all_box(:)
integer, private, allocatable, save :: all_n(:)
NFE_REAL, private, save :: my_w
integer, private, save :: my_n, n_tot
character(len = SL), private, save :: MDOUT_FILE
#endif /* NFE_ENABLE_BBMD */
! Moradi end

!
! rem-specific, for old subroutiens, not used now
!

! rem_* arrays live on sander masters; indexed by masterrank

character(len = SL), private, allocatable, save :: rem_monitor_files(:)
integer, private, allocatable, save :: rem_monitor_freqs(:)

character(len = SL), private, allocatable, save :: rem_umbrella_files(:)
type(umbrella_t), private, allocatable, save :: rem_umbrellas(:)

character(len = SL), private, allocatable, save :: rem_snapshots_basenames(:)
integer, private, allocatable, save :: rem_snapshots_freqs(:)

private :: rem_postinit
private :: rem_cleanup

integer, private, save :: LOG_UNIT
integer, private, save :: exchange_partner

#else
#  define LOG_UNIT OUT_UNIT
#endif /* MPI */

NFE_REAL, parameter, private :: TINY = 0.000001d0 ! NFE_TO_REAL(0.00001)

namelist / abmd /    mode, monitor_file, monitor_freq, timescale,&
                     umbrella_file, snapshots_basename, snapshots_freq,&
                     selection_freq, selection_constant, selection_epsilon,&
                     wt_temperature, wt_umbrella_file, cv_file, &
                     driven_weight, driven_cutoff

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

#ifdef MPI
subroutine on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)

   NFE_USE_AFAILED

   use nfe_umbrella
   use nfe_sander_proxy

   implicit none

   integer, intent(in) :: o_masterrank
   logical, intent(in) :: need_U_xx

   NFE_REAL, intent(inout) :: U_mm, U_mo, U_om, U_oo

#  include "nfe-mpi.h"

   NFE_REAL :: o_instantaneous(MAX_NUMBER_OF_COLVARS), U_o(2), U_m(2)

   integer :: error

   if (imode.eq.MODE_NONE.or.imode.eq.MODE_ANALYSIS) &
      return

   nfe_assert(multisander_rem().eq.1)
   nfe_assert(sanderrank.eq.0) ! master
   nfe_assert(commmaster.ne.mpi_comm_null)

   ! exchange instantaneous(:) with the partner
   call mpi_sendrecv &
      (instantaneous, ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       o_instantaneous, ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       commmaster, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)

   ! evaluate 'my' values
   U_m(1) = umbrella_eval_v(umbrella, instantaneous)   ! U_mm = U_m(x_m)
   U_m(2) = umbrella_eval_v(umbrella, o_instantaneous) ! U_mo = U_m(x_o)

   ! get partner's U_m? (i.e., U_o? in this replica)
   call mpi_sendrecv &
      (U_m, 2, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       U_o, 2, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
       commmaster, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)

   if (need_U_xx) then
      U_mm = U_mm + U_m(1)
      U_mo = U_mo + U_m(2)
      U_om = U_om + U_o(2)
      U_oo = U_oo + U_o(1)
   end if

end subroutine on_delta

!-----------------------------------------------------------------------------

subroutine on_exchange(o_masterrank)

#ifndef NFE_DISABLE_ASSERT
   use nfe_utils
   use nfe_sander_proxy
#endif /* NFE_DISABLE_ASSERT */

   implicit none

   integer, intent(in) :: o_masterrank

#  include "nfe-mpi.h"

   if (imode.eq.MODE_NONE) &
      return

   nfe_assert(multisander_rem().eq.1)
   nfe_assert(sanderrank.eq.0) ! master
   nfe_assert(commmaster.ne.mpi_comm_null)

   exchange_partner = o_masterrank

end subroutine on_exchange
#endif /* MPI */

!-----------------------------------------------------------------------------

subroutine on_multisander_exit()

   NFE_USE_AFAILED

   use nfe_colvar
   use nfe_umbrella
   use nfe_sander_proxy

   implicit none

#  include "nfe-mpi.h"

   integer :: n

   if (imode.ne.MODE_NONE) then
      nfe_assert(ncolvars.gt.0)
      do n = 1, ncolvars
         call colvar_cleanup(colvars(n))
      end do
   end if

   NFE_MASTER_ONLY_BEGIN

#ifdef MPI
   if (allocated(all_hills)) &
      deallocate(all_hills)
! Modified by M Moradi
! Selection algorithm
#ifdef NFE_ENABLE_BBMD
   if (multisander_numgroup().gt.1) then
   if (allocated(all_w)) &
      deallocate(all_w)
   if (allocated(all_x)) &
      deallocate(all_x)
   if (allocated(all_v)) &
      deallocate(all_v)
   if (allocated(all_box)) &
      deallocate(all_box)
   end if
#endif /* NFE_ENABLE_BBMD */
! Moradi end
#endif /* MPI */

   if (imode.eq.MODE_FLOODING.or.imode.eq.MODE_UMBRELLA) &
      call umbrella_fini(umbrella)
   NFE_MASTER_ONLY_END

   imode = MODE_NONE

end subroutine on_multisander_exit

!-----------------------------------------------------------------------------

!
! 'on_sander_init()' is called at the point when
! MPI is initialized and MDIN/PRMTOP/INPCRD are loaded
!

subroutine on_sander_init(mdin_unit, amass)

   use nfe_utils
   use nfe_colvar
   use nfe_colvar_type
   use nfe_constants
   use nfe_umbrella
   use nfe_sander_proxy

   implicit none

   integer, intent(in) :: mdin_unit

   NFE_REAL, intent(in) :: amass(*)

   logical :: umbrella_file_exists, do_transfer
   type(umbrella_t) :: umbrella_from_file

   integer :: n, error, i

   integer :: cv_extents(UMBRELLA_MAX_NEXTENTS)
   logical :: cv_periodicity(UMBRELLA_MAX_NEXTENTS)

   NFE_REAL :: cv_origin(UMBRELLA_MAX_NEXTENTS)
   NFE_REAL :: cv_spacing(UMBRELLA_MAX_NEXTENTS)

   NFE_REAL :: tmp

   integer             :: ifind
   character(len = 80) :: buf   

#ifdef MPI
#  include "nfe-mpi.h"
#  include "../include/md.h"
#endif /* MPI */
   integer, save :: counter = 0

   do_transfer = .false.

   if(counter == 0) mdstep = 0 ! used in on_force() below

#ifdef MPI
   if (counter.gt.0) then
      nfe_assert(multisander_rem().ne.0)

      NFE_MASTER_ONLY_BEGIN
      if (imode.ne.MODE_NONE) then
         ! re-open MONITOR_UNIT after exchange (closed in on_sander_exit())
         open (unit = MONITOR_UNIT, file = monitor_file, &
               iostat = error, form = 'FORMATTED', action = 'WRITE', &
               position = 'APPEND', status = 'OLD')

         if (error.ne.0) then
            write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NFE_ERROR, &
               'could not open ''', trim(monitor_file), ''' file for writing'
            call terminate()
         end if
      end if ! imode.ne.MODE_NONE
      NFE_MASTER_ONLY_END

      return  ! return for REMD
   end if ! counter.gt.0

   counter = counter + 1

   NFE_MASTER_ONLY_BEGIN

   if (multisander_numgroup().gt.1) &
      n = masterrank + 1

   if (multisander_numgroup().gt.1) then
      write (unit = monitor_file, fmt = '(a,a,i3.3)') &
         DEFAULT_MONITOR_FILE, '-', n
      write (unit = umbrella_file, fmt = '(a,a,i3.3,a)') &
         DEFAULT_UMBRELLA_FILE, '-', n, '.nc'
      write (unit = wt_umbrella_file, fmt = '(a,a,i3.3,a)') &
         DEFAULT_WT_UMBRELLA_FILE, '-', n, '.nc'         
      write (unit = snapshots_basename, fmt = '(a,a,i3.3)') &
         DEFAULT_SNAPSHOTS_BASENAME, '-', n
   else
      monitor_file = DEFAULT_MONITOR_FILE
      umbrella_file = DEFAULT_UMBRELLA_FILE//'.nc'
      wt_umbrella_file = DEFAULT_WT_UMBRELLA_FILE//'.nc'      
      snapshots_basename = DEFAULT_SNAPSHOTS_BASENAME
   end if
#else
   monitor_file = DEFAULT_MONITOR_FILE
   umbrella_file = DEFAULT_UMBRELLA_FILE
   wt_umbrella_file = DEFAULT_WT_UMBRELLA_FILE
   snapshots_basename = DEFAULT_SNAPSHOTS_BASENAME
#endif /* MPI */

   ! parse MDIN (on master)
   nfe_assert(ncolvars.eq.0)

   rewind(mdin_unit)
   call nmlsrc('abmd', mdin_unit, ifind)
   ! no abmd section
   if (ifind.eq.0) then
      imode = MODE_NONE
      goto 1
   end if

   rewind(mdin_unit)
   read(mdin_unit,nml=abmd,err=666)
   
   ! discover the run-mode

   if (mode == 'NONE') then
      imode = MODE_NONE
      goto 1
   else if (mode == 'ANALYSIS') then
      imode = MODE_ANALYSIS
   else if (mode == 'UMBRELLA') then
      imode = MODE_UMBRELLA
   else if (mode == 'FLOODING') then
      imode = MODE_FLOODING
   else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'unknown mode ''', trim(mode), ''''
      call terminate()
   end if

#ifdef NFE_NO_NETCDF
   umbrella_file_exists = .false.
   write (unit = ERR_UNIT, fmt = '(a,a)') NFE_WARNING, &
      'netCDF is not available (try ''-bintraj'' configure option)'
#else
   inquire (file = umbrella_file, exist = umbrella_file_exists)
#endif /* NFE_NO_NETCDF */

   if (.not.umbrella_file_exists.and.imode.eq.MODE_UMBRELLA) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NFE_ERROR, '''', &
         trim(umbrella_file), ''' does not exist (required for UMBRELLA mode)'
      call terminate()
   end if

   if (imode.eq.MODE_ANALYSIS) &
      umbrella_file_exists = .false.

#ifndef NFE_NO_NETCDF
   if (umbrella_file_exists) &
      call umbrella_load(umbrella_from_file, umbrella_file)
#endif /* NFE_NO_NETCDF */

   ! collective variables
   nfe_assert(ncolvars.eq.0)

   call amopen(CV_UNIT, cv_file, 'O', 'F', 'R')
   
   do
     call nmlsrc('colvar', CV_UNIT, ifind)
     if (ifind.eq.0) exit
     read(CV_UNIT,'(a80)') buf
     ncolvars = ncolvars + 1
   end do

   if (ncolvars.eq.0) &
      call fatal('no variable(s) in the CV file')

   if (ncolvars.gt.MAX_NUMBER_OF_COLVARS) &
      call fatal('too many variables in the CV file')

   if (umbrella_file_exists) then
      if(umbrella_nextents(umbrella_from_file).ne.ncolvars) &
         call fatal('number of variables in the CV file does not &
                 &match with the number of extents found in the umbrella_file')
   end if ! umbrella_file_exists

   n = 1

   do while (n.le.ncolvars)

         call colvar_nlread(CV_UNIT,colvars(n))

         if (imode.eq.MODE_FLOODING) then
            cv_spacing(n) = resolution
            cv_spacing(n) = cv_spacing(n)/4
            if ((resolution.eq.ZERO).and..not.umbrella_file_exists) then
               write (unit = ERR_UNIT, fmt = '(/a,a,i1/)') NFE_ERROR, &
                  'could not determine ''resolution'' for CV #', n
               call terminate()
            end if

            if (resolution.eq.ZERO) &
               cv_spacing(n) = umbrella_spacing(umbrella_from_file, n)

            cv_periodicity(n) = colvar_is_periodic(colvars(n))

            if (cv_periodicity(n)) then

               nfe_assert(colvar_has_min(colvars(n)))
               nfe_assert(colvar_has_max(colvars(n)))

               cv_origin(n) = colvar_min(colvars(n))

               nfe_assert(cv_spacing(n).gt.ZERO)
               nfe_assert(colvar_max(colvars(n)).gt.cv_origin(n))

               cv_extents(n) = &
                  int((colvar_max(colvars(n)) - cv_origin(n))/cv_spacing(n))

               if (cv_extents(n).lt.UMBRELLA_MIN_EXTENT) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1,a/)') NFE_ERROR, &
                     'CV #', n, ' : ''resolution'' is too big'
                  call terminate()
               end if

               cv_spacing(n) = &
                  (colvar_max(colvars(n)) - cv_origin(n))/cv_extents(n)

            else ! .not.periodic

               cv_origin(n) = cv_min
               tmp = cv_max

               if (cv_origin(n).ge.tmp) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1/)') NFE_ERROR, &
                     'min.ge.max for CV #', n
                  call terminate()
               end if

               nfe_assert(cv_spacing(n).gt.ZERO)
               cv_extents(n) = 1 &
                  + int((tmp - cv_origin(n))/cv_spacing(n))

               if (cv_extents(n).lt.UMBRELLA_MIN_EXTENT) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1,a/)') NFE_ERROR, &
                     'CV #', n, ' : the ''resolution'' is too big'
                  call terminate()
               end if

               cv_spacing(n) = (tmp - cv_origin(n))/(cv_extents(n) - 1)

            end if ! cv_periodicity(n)
         end if ! imode.eq.MODE_FLOODING
         if (colvar_has_refcrd(colvars(n))) then
            refcrd_file = trim(refcrd_file)
            refcrd_len = len_trim(refcrd_file)
         end if
         if (colvar_is_quaternion(colvars(n))) then 
          allocate(colvars(n)%q_index, stat = error)
           if (error.ne.0) &
            NFE_OUT_OF_MEMORY
             colvars(n)%q_index = q_index
         end if 
         if (colvar_has_axis(colvars(n))) then
           allocate(colvars(n)%axis(3), stat = error)
             if (error.ne.0) &
                NFE_OUT_OF_MEMORY
           i = 1
           do while (i.le.3)
             colvars(n)%axis(i) = axis(i)
             i = i + 1
           end do
         end if

         n = n + 1
   end do
   close(unit=CV_UNIT)

   monitor_freq = min(monitor_freq, sander_nstlim())
   monitor_freq = max(1, monitor_freq)

! Modified by M Moradi
! Selection algorithm
#ifdef NFE_ENABLE_BBMD  
   selection_epsilon = min(ONE,selection_epsilon)
#endif /* NFE_ENABLE_BBMD */
! Well-tempered ABMD
   if (wt_temperature .gt. ZERO) then
      pseudo = ONE/wt_temperature
      call umbrella_init(wt_umbrella, ncolvars, cv_extents, &
                         cv_origin, cv_spacing, cv_periodicity)
   end if
!  Driven ABMD
   if (driven_weight == 'NONE') then
      drivenw = 0
      drivenu = 0
   else if (driven_weight == 'CONSTANT') then
      drivenw = 1
      drivenu = 1
   else if (driven_weight == 'PULLING') then
      drivenw = 1
      drivenu = 0
   else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'unknown driven weight scheme ''', trim(driven_weight), ''''
      call terminate()
   end if
! Moradi end

1  continue ! done with parsing

   NFE_MASTER_ONLY_END

#ifdef MPI
   call mpi_bcast(imode, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

   if (imode.eq.MODE_NONE) &
      return

#ifdef MPI
   nfe_assert(.not.is_master().or.ncolvars.gt.0)

   call mpi_bcast(ncolvars, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)
   
   call mpi_bcast(refcrd_len, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(refcrd_file, refcrd_len, MPI_CHARACTER, 0, commsander, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(monitor_freq, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

   nfe_assert(ncolvars.gt.0)
   nfe_assert(ncolvars.le.MAX_NUMBER_OF_COLVARS)

   if (multisander_numwatkeep().gt.0) &
      call fatal('numwatkeep.gt.0 is not supported')

   if (sander_imin().ne.0) &
      call fatal('imin.ne.0 is not supported')
   do n = 1, ncolvars
      call colvar_bootstrap(colvars(n), n, amass)
   end do
#ifdef MPI
   if( multisander_numgroup().gt.1 ) &
      call rem_checks()
#endif /* MPI */

#ifdef NFE_ENABLE_BBMD
   if (multisander_numgroup().gt.1) then
   call mpi_bcast(selection_freq, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)
   MDOUT_FILE = sander_mdout_name()
   call mpi_bcast(MDOUT_FILE, len(MDOUT_FILE), MPI_CHARACTER, 0, commsander, error)
   nfe_assert(error.eq.0)
   end if
#endif /* NFE_ENABLE_BBMD */

   NFE_MASTER_ONLY_BEGIN

   if (imode.eq.MODE_UMBRELLA) then
      nfe_assert(umbrella_file_exists)
      do_transfer = .false.
      call umbrella_swap(umbrella, umbrella_from_file)
! Modified by M Moradi
! Selection algorithm
#ifdef NFE_ENABLE_BBMD
   if (multisander_numgroup().gt.1) then 
      if (selection_freq.gt.0) then
          my_w = ZERO
          allocate(all_x(3*sander_natoms()*multisander_numgroup()), &
             stat = error)
          allocate(all_v(3*sander_natoms()*multisander_numgroup()), &
             stat = error)
          if (sander_ntb().ne.0) &
             allocate(all_box(3*multisander_numgroup()), &
             stat = error)
          allocate(all_w(multisander_numgroup()), &
             stat = error)
          allocate(all_n(multisander_numgroup()), &
             stat = error)
          call amrset_gen(nfe_abmd_gen, ig+masterrank)
      end if ! selection_freq.gt.0
   end if
#endif /* NFE_ENABLE_BBMD */
! Moradi end      
   else if (imode.eq.MODE_FLOODING) then
      if (umbrella_file_exists) then
         do_transfer = .false.
         do n = 1, ncolvars
            do_transfer = do_transfer &
               .or.(cv_extents(n).ne.umbrella_extent(umbrella_from_file, n))
            do_transfer = do_transfer &
               .or.(cv_periodicity(n).neqv.&
                  umbrella_periodicity(umbrella_from_file, n))
            do_transfer = do_transfer &
               .or.(abs(cv_origin(n) - umbrella_origin(umbrella_from_file, n)) &
                  .gt.TINY)
            do_transfer = do_transfer &
               .or.(abs(cv_spacing(n) - umbrella_spacing(umbrella_from_file, &
                  n)).gt.TINY)
            if (do_transfer) &
               exit
         end do
         if (do_transfer) then
            call umbrella_init(umbrella, ncolvars, cv_extents, &
                               cv_origin, cv_spacing, cv_periodicity)
            call umbrella_transfer(umbrella, umbrella_from_file)
            call umbrella_fini(umbrella_from_file)
         else
            call umbrella_swap(umbrella, umbrella_from_file)
         end if ! do_transfer
      else
         call umbrella_init(umbrella, ncolvars, cv_extents, &
                            cv_origin, cv_spacing, cv_periodicity)
      end if ! umbrella_file_exits
#ifdef MPI
      if (multisander_numgroup().gt.1) then
         allocate(all_hills((ncolvars + 1)*multisander_numgroup()), &
            stat = error)
! Modified by M Moradi
! Selection algorithm
#ifdef NFE_ENABLE_BBMD
         if (selection_freq.gt.0) then
            my_w = ZERO
            allocate(all_x(3*sander_natoms()*multisander_numgroup()), &
               stat = error)
            allocate(all_v(3*sander_natoms()*multisander_numgroup()), &
               stat = error)
            if (sander_ntb().ne.0) &
               allocate(all_box(3*multisander_numgroup()), &
                  stat = error)
            allocate(all_w(multisander_numgroup()), &
               stat = error)
            allocate(all_n(multisander_numgroup()), &
               stat = error)
            call amrset_gen(nfe_abmd_gen, ig+masterrank)
         end if ! selection_freq.gt.0
#endif /* NFE_ENABLE_BBMD */
! Moradi end            
         if (error.ne.0) &
            NFE_OUT_OF_MEMORY
         if (multisander_rem().eq.0) then
            if (masterrank.gt.0) &
               call umbrella_fini(umbrella)

            call umbrella_bcast(umbrella, commmaster, 0)

            call mpi_bcast(timescale, 1, &
               MPI_DOUBLE_PRECISION, 0, commmaster, error)
            nfe_assert(error.eq.0)
         end if ! .not.REM
      end if ! ng.gt.1
#endif /* MPI */
   end if

   ! prepare monitor_fmt & open MONITOR_UNIT

   open (unit = MONITOR_UNIT, file = monitor_file, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

   if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'failed to open ''', trim(monitor_file), ''' for writing'
      call terminate()
   end if

   write (unit = MONITOR_UNIT, fmt = '(a,/a)', advance = 'NO') &
      '#', '# MD time (ps), '
   do n = 1, ncolvars - 1
      write (unit = MONITOR_UNIT, fmt = '(a,i1,a)', advance = 'NO') &
         'CV #', n, ', '
   end do
   if (imode == MODE_FLOODING) then
      write (unit = MONITOR_UNIT, fmt = '(a,i1,a,/a)') &
         'CV #', ncolvars, ', E_{bias} (kcal/mol)', '#'
      write (unit = monitor_fmt, fmt = '(a,i1,a)') &
         '(f12.4,', ncolvars, '(1x,f16.10),1x,f16.10)'
   else
      write (unit = MONITOR_UNIT, fmt = '(a,i1,/a)') &
         'CV #', ncolvars, '#'
      write (unit = monitor_fmt, fmt = '(a,i1,a)') &
         '(f12.4,', ncolvars, '(1x,f16.10))'
   end if

   call flush_UNIT(MONITOR_UNIT)

   ! print summary & return

#ifdef MPI
   LOG_UNIT = OUT_UNIT ! write to MDOUT
#endif

   write (unit = LOG_UNIT, fmt = '(/a,a)') NFE_INFO, &
      '() () () () () () () () () ()   A. B. M. D.  () () () () () () () () ()'
   write (unit = LOG_UNIT, fmt = '(a,/a,a)', advance = 'NO') &
      NFE_INFO, NFE_INFO, 'mode = '

   select case(imode)
      case(MODE_ANALYSIS)
         write (unit = LOG_UNIT, fmt = '(a)') 'ANALYSIS'
      case(MODE_UMBRELLA)
         write (unit = LOG_UNIT, fmt = '(a)') 'UMBRELLA'
      case(MODE_FLOODING)
         write (unit = LOG_UNIT, fmt = '(a)') 'FLOODING'
      case default
         nfe_assert_not_reached()
         continue
   end select

   write (unit = LOG_UNIT, fmt = '(a)') NFE_INFO

   do n = 1, ncolvars
      write (unit = LOG_UNIT, fmt = '(a,a,i1)') NFE_INFO, 'CV #', n
      call colvar_print(colvars(n), LOG_UNIT)
      write (unit = LOG_UNIT, fmt = '(a)') NFE_INFO
      if (colvar_is_quaternion(colvars(n))) then
       write (unit = OUT_UNIT, fmt = '(a,a,I3)') NFE_INFO, &
          ' q_index = ', colvars(n)%q_index
      end if
      if (colvar_has_axis(colvars(n))) then
        write (unit = OUT_UNIT, fmt = '(a,a,f8.4,a,f8.4,a,f8.4,a)') NFE_INFO, &
        ' axis = [',colvars(n)%axis(1),', ', colvars(n)%axis(2),',',colvars(n)%axis(3),']'
      end if
   end do

   write (unit = LOG_UNIT, fmt = '(a,a,a)') NFE_INFO, &
      'monitor_file = ', trim(monitor_file)
   write (unit = LOG_UNIT, &
   fmt = '(a,a,'//pfmt(monitor_freq)//',a,'//pfmt(&
   monitor_freq*sander_timestep(), 4)//',a)') NFE_INFO, &
      'monitor_freq = ', monitor_freq, ' (', &
      monitor_freq*sander_timestep(), ' ps)'

   if (imode.eq.MODE_ANALYSIS) &
      goto 3

#  ifdef MPI
   if (multisander_numgroup().gt.1.and.multisander_rem().eq.0) then
      if (masterrank.gt.0) &
         write (unit = LOG_UNIT, fmt = '(a,/a,a,/a)') NFE_INFO, NFE_INFO, &
            'ng.gt.1.and.rem.eq.0 => using umbrella from replica #1', NFE_INFO
   end if
#  endif /* MPI */

   write (unit = LOG_UNIT, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
      'umbrella_file = ', trim(umbrella_file), ' ('

   if (umbrella_file_exists) then
      write (unit = LOG_UNIT, fmt = '(a)') 'loaded)'
   else
      write (unit = LOG_UNIT, fmt = '(a)') 'not found)'
   end if

   write (unit = LOG_UNIT, fmt = '(a)') NFE_INFO
   write (unit = LOG_UNIT, fmt = '(a,a)', advance = 'NO') NFE_INFO, &
      'umbrella discretization '

   if (umbrella_file_exists) then
      if (do_transfer) then
         write (unit = LOG_UNIT, fmt = '(a)') '(modified) :'
      else
         write (unit = LOG_UNIT, fmt = '(a)') '(unchanged) :'
      end if
   else
      write (unit = LOG_UNIT, fmt = '(a)') '(new) :'
   end if

   do n = 1, ncolvars
      write (unit = LOG_UNIT, fmt = '(a,a,i1)', advance = 'NO') &
         NFE_INFO, 'CV #', n
      if (umbrella_periodicity(umbrella, n)) then
         write (unit = LOG_UNIT, fmt = '(a)', advance = 'NO') ' periodic, '
         tmp = umbrella_origin(umbrella, n) &
            + umbrella_spacing(umbrella, n)*umbrella_extent(umbrella, n)
      else
         write (unit = LOG_UNIT, fmt = '(a)', advance = 'NO') ' not periodic, '
         tmp = umbrella_origin(umbrella, n) &
            + umbrella_spacing(umbrella, n)*(umbrella_extent(umbrella, n) - 1)
      end if

      write (unit = LOG_UNIT, &
         fmt = '('//pfmt(umbrella_extent(umbrella, &
         n))//',a,'//pfmt(umbrella_origin(umbrella, n), &
         6)//',a,'//pfmt(tmp, 6)//')') &
         umbrella_extent(umbrella, n), ' points, min/max = ', &
         umbrella_origin(umbrella, n), '/', tmp
   end do

   if (imode.eq.MODE_UMBRELLA) &
      goto 3

   write (unit = LOG_UNIT, fmt = '(a/,a,a,'//pfmt(timescale, 3)//',a)') &
      NFE_INFO, NFE_INFO, 'flooding timescale = ', timescale, ' ps'

   if (snapshots_freq.gt.0) then
      write (unit = LOG_UNIT, fmt = '(a,a,a)') NFE_INFO, &
         'snapshots_basename = ', trim(snapshots_basename)
      write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(snapshots_freq)//',a,'//pfmt(snapshots_freq*sander_timestep(), 4)//',a)')&
          NFE_INFO, 'snapshots_freq = ', snapshots_freq, ' (', &
         snapshots_freq*sander_timestep(), ' ps)'
   end if

   nhills = 0

! Modified by M Moradi
! Well-tempered ABMD
   if (pseudo.gt.ZERO) then
      write (unit = LOG_UNIT, &
      fmt = '(a/,a,a)')&
          NFE_INFO, NFE_INFO,'well-tempered ABMD:'
      write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(ONE/pseudo,6)//')')&
          NFE_INFO, 'pseudo-temperature = ', ONE/pseudo
      write (unit = LOG_UNIT, fmt = '(a,a,a)') NFE_INFO, &
         'wt_umbrella_file = ', trim(wt_umbrella_file)          
   end if
! Driven ABMD
   if (drivenw.gt.ZERO) then
      write (unit = LOG_UNIT, &
      fmt = '(a/,a,a)')&
          NFE_INFO, NFE_INFO,'driven ABMD:'
      if (drivenu.gt.ZERO) then
         write (unit = LOG_UNIT, &
         fmt = '(a,a)')&
             NFE_INFO,'CONSTANT weighting scheme (use delta work)'
      else
         write (unit = LOG_UNIT, &
         fmt = '(a,a)')&
             NFE_INFO,'PULLING weighting scheme (use work only)'
      end if
      write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(driven_cutoff,6)//')')&
          NFE_INFO, 'driven (delta)work cutoff = ', driven_cutoff
   end if
3 continue
! Selection algorithm
#ifdef NFE_ENABLE_BBMD
   if (multisander_numgroup().gt.1) then
     if (selection_freq.gt.0) then
      write (unit = LOG_UNIT, &
      fmt = '(a/,a,a/,a)')&
          NFE_INFO, NFE_INFO,'selection algorithm parameters:', NFE_INFO
      write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(selection_freq)//',a,'//pfmt(selection_freq*sander_timestep(), 4)//',a)')&
          NFE_INFO, 'selection_freq = ', selection_freq, ' (', &
         selection_freq*sander_timestep(), ' ps)'
      write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(selection_constant,6)//')')&
          NFE_INFO, 'selection scoring constant = ', selection_constant
      write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(selection_epsilon,6)//')')&
          NFE_INFO, 'selection criterion epsilon = ', selection_epsilon
     end if
   end if
#endif /* NFE_ENABLE_BBMD */
! Moradi end
   write (unit = OUT_UNIT, fmt = '(a)') NFE_INFO
   write (unit = LOG_UNIT, fmt = '(a,a/)') NFE_INFO, &
      '() () () () () () () () () () () () () () () () () () () () () () () ()'
   call flush_UNIT(LOG_UNIT)

#ifndef MPI
   return
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &abmd namelist!'
    call terminate()
#endif

#ifdef MPI

   NFE_MASTER_ONLY_END

   return
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &abmd namelist!'
    call terminate()
contains

!.............................................................................

subroutine rem_checks()

   implicit none

   integer, allocatable :: int_recv(:)

   if (commmaster.eq.mpi_comm_null) &
      return

   nfe_assert(mastersize.gt.0)
   nfe_assert(masterrank.lt.mastersize)

   ! For H-REMD and Multi-REMD, imode and ncolvar can vary among replicas   
   if (multisander_rem().eq.3 .or. multisander_rem().eq.-1) &
      return

   allocate(int_recv(mastersize), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   call mpi_allgather(imode, 1, MPI_INTEGER, int_recv, 1, MPI_INTEGER, &
                      commmaster, error)
   nfe_assert(error.eq.0)

   do n = 1, mastersize
      if (int_recv(n).ne.imode) &
         call fatal('''mode'' has different values in different replicas')
   end do

   call mpi_allgather(ncolvars, 1, MPI_INTEGER, int_recv, 1, MPI_INTEGER, &
                      commmaster, error)
   nfe_assert(error.eq.0)

   do n = 1, mastersize
      if (int_recv(n).ne.ncolvars) &
         call fatal('number of collective variables is &
                    &different in different replicas')
   end do

   ! FIXME: more here

   deallocate(int_recv)

end subroutine rem_checks
#endif /* MPI */

end subroutine on_sander_init

!-----------------------------------------------------------------------------

subroutine on_sander_exit()

   use nfe_utils, only : close_UNIT

   implicit none

#  include "nfe-mpi.h"

   NFE_MASTER_ONLY_BEGIN
   call close_UNIT(MONITOR_UNIT)
   NFE_MASTER_ONLY_END

end subroutine on_sander_exit

!-----------------------------------------------------------------------------

! Modified by M Moradi
! Driven ABMD
subroutine on_force(x, f, wdriven, udriven, pot)
! Moradi end

   use nfe_utils
   use nfe_colvar
   use nfe_umbrella
   use nfe_constants
   use nfe_sander_proxy
   use nfe_colvar_type

   implicit none

   NFE_REAL, intent(in) :: x(*)

   NFE_REAL, intent(inout) :: f(*)
   
! Modified by M Moradi
! Driven ABMD
   NFE_REAL, intent(in) :: wdriven
   NFE_REAL, intent(in) :: udriven
! Moradi end
   
   NFE_REAL, intent(inout) :: pot

#ifdef MPI
#  include "nfe-mpi.h"
   integer :: error
#endif /* MPI */

   NFE_REAL :: u_value, u_derivative(UMBRELLA_MAX_NEXTENTS), alt
   integer, DIMENSION(4) :: cv_q = (/COLVAR_QUATERNION0, COLVAR_QUATERNION1, &
   COLVAR_QUATERNION2, COLVAR_QUATERNION3/)
   NFE_REAL :: norm4(100), cv_N(4, 100)


   character(len = SL + 16) :: snapshot

   integer :: n, m

   if (imode.eq.MODE_NONE) &
      return

   nfe_assert(ncolvars.gt.0)

   if (imode.eq.MODE_ANALYSIS) then
      if (nfe_real_mdstep.and.mod(mdstep, monitor_freq).eq.0) then
         do n = 1, ncolvars
            instantaneous(n) = colvar_value(colvars(n), x)
         end do

         NFE_MASTER_ONLY_BEGIN
         do n = 1, ncolvars
           if (colvar_is_quaternion(colvars(n))) then
             do m = 1, 4
                if (colvars(n)%type == cv_q(m)) then
                   cv_N(m, colvars(n)%q_index) = instantaneous(n)
                end if
             end do
           end if
         end do
         do n = 1, ncolvars
           if (colvar_is_quaternion(colvars(n))) then
              norm4(colvars(n)%q_index) = sqrt(cv_N(1,colvars(n)%q_index)**2 &
                                             + cv_N(2,colvars(n)%q_index)**2 &
                                             + cv_N(3,colvars(n)%q_index)**2 &
                                             + cv_N(4,colvars(n)%q_index)**2)
           end if
         end do
         do n = 1, ncolvars
           if (colvar_is_quaternion(colvars(n))) then
              instantaneous(n) = instantaneous(n) / norm4(colvars(n)%q_index)
           else
              instantaneous(n) = instantaneous(n)
           end if
         end do

         write (unit = MONITOR_UNIT, fmt = monitor_fmt) &
            sander_mdtime(), instantaneous(1:ncolvars)
         call flush_UNIT(MONITOR_UNIT)
         NFE_MASTER_ONLY_END
      end if
      goto 1
   end if ! imode.eq.MODE_ANALYSIS

   !
   ! either UMBRELLA or FLOODING
   !

   do n = 1, ncolvars
      instantaneous(n) = colvar_value(colvars(n), x)
   end do

   NFE_MASTER_ONLY_BEGIN
   do n = 1, ncolvars
    if (colvar_is_quaternion(colvars(n))) then
     do m = 1, 4
       if (colvars(n)%type == cv_q(m)) then
         cv_N(m, colvars(n)%q_index) = instantaneous(n)
       end if
     end do
    end if
   end do
   do n = 1, ncolvars
    if (colvar_is_quaternion(colvars(n))) then
      norm4(colvars(n)%q_index) = sqrt(cv_N(1,colvars(n)%q_index)**2 &
                                     + cv_N(2,colvars(n)%q_index)**2 &
                                     + cv_N(3,colvars(n)%q_index)**2 &
                                     + cv_N(4,colvars(n)%q_index)**2)
    end if
   end do
   do n = 1, ncolvars
    if (colvar_is_quaternion(colvars(n))) then
      instantaneous(n) = instantaneous(n) / norm4(colvars(n)%q_index)
    else
       instantaneous(n) = instantaneous(n)
    end if
   end do

   call umbrella_eval_vdv(umbrella, instantaneous, u_value, u_derivative)
   pot = umbrella_eval_v(umbrella, instantaneous)
   NFE_MASTER_ONLY_END

#ifdef MPI
   call mpi_bcast(u_derivative, ncolvars, &
      MPI_DOUBLE_PRECISION, 0, commsander, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(pot, 1, &
      MPI_DOUBLE_PRECISION, 0, commsander, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

   ! FIXME: virial
   do n = 1, ncolvars
      call colvar_force(colvars(n), x, -u_derivative(n), f)
   end do

   if (nfe_real_mdstep) then
      NFE_MASTER_ONLY_BEGIN
      if (mod(mdstep, monitor_freq).eq.0) then
         if (imode.eq.MODE_FLOODING) then
            write (unit = MONITOR_UNIT, fmt = monitor_fmt) &
               sander_mdtime(), instantaneous(1:ncolvars), u_value
         else
            write (unit = MONITOR_UNIT, fmt = monitor_fmt) &
               sander_mdtime(), instantaneous(1:ncolvars)
         end if
         call flush_UNIT(MONITOR_UNIT)
      end if

      if (imode.eq.MODE_FLOODING) then
#ifndef NFE_NO_NETCDF
         if (snapshots_freq.gt.0.and.mod(nhills, snapshots_freq).eq.0) then
            write (unit = snapshot, fmt = '(a,a,i10.10,a)') &
               trim(snapshots_basename), '.', nhills, '.nc'
            call umbrella_save(umbrella, snapshot)
            write (unit = OUT_UNIT, fmt = '(/a,a,f16.4,a,/a,a,a,a)') &
               NFE_INFO, 'biasing potential snapshot at t = ', &
               sander_mdtime(), ' ps', NFE_INFO, 'saved as ''', &
               trim(snapshot), ''''
         end if
#endif /* NFE_NO_NETCDF */
         alt = sander_timestep()/timescale
! Modified by M Moradi
! Well-tempered ABMD (based on Barducci, Bussi, and Parrinello, PRL(2008) 100:020603)
         if (pseudo.gt.ZERO) &
            alt = alt * &
            exp(-pseudo*umbrella_eval_v(umbrella,instantaneous)/kB)
! Driven ABMD (based on Moradi and Tajkhorshid, JPCL(2013) 4:1882)
         if (drivenw.ne.0) then
            if (drivenw*wdriven-drivenu*udriven.gt.driven_cutoff) then
               alt = alt * &
                exp(-(drivenw*wdriven-drivenu*udriven)/(kB*sander_temp0()))
            else
               alt = alt * &
                exp(-driven_cutoff/(kB*sander_temp0()))
            end if
         end if ! drivenw.ne.0
! Moradi end         
#ifdef MPI
         if (multisander_numgroup().gt.1) then
            nfe_assert(commmaster.ne.mpi_comm_null)
            ! get all instantaneous/altitudes
            nfe_assert(allocated(all_hills))
            instantaneous(ncolvars + 1) = alt
            call mpi_allgather(instantaneous, ncolvars + 1, &
               MPI_DOUBLE_PRECISION, all_hills, ncolvars + 1, &
               MPI_DOUBLE_PRECISION, commmaster, error)
            nfe_assert(error.eq.0)

            if (multisander_rem().eq.0) then
               do n = 0, multisander_numgroup() - 1
                  call umbrella_hill(umbrella, &
                     all_hills(n*(ncolvars + 1) + 1:), &
                     all_hills((n + 1)*(ncolvars + 1)))
               end do
            else
               call umbrella_hill(umbrella, instantaneous, alt)
            end if
         else
#endif /* MPI */
            call umbrella_hill(umbrella, instantaneous, alt)
#ifdef MPI
         end if ! multisander_numgroup().gt.1
#endif /* MPI */
         nhills = nhills + 1
      end if
      NFE_MASTER_ONLY_END

   end if ! nfe_real_mdstep

1  if (nfe_real_mdstep) &
   mdstep = mdstep + 1

end subroutine on_force

!-----------------------------------------------------------------------------

subroutine on_mdwrit()

   use nfe_utils
   use nfe_umbrella
   use nfe_sander_proxy

   implicit none
   double precision :: t
#ifndef NFE_NO_NETCDF
   if (imode.eq.MODE_FLOODING) then
      nfe_assert(is_master())
      call umbrella_save(umbrella, umbrella_file)
! Modified by F Pan      
      if (pseudo.gt.TINY) then
         t = sander_temp0()
         call umbrella_copy(wt_umbrella,umbrella)
         call umbrella_wt_mod(wt_umbrella,t,1/pseudo)
         call umbrella_save(wt_umbrella, wt_umbrella_file)
      end if
! Pan end            
   end if
#endif /* NFE_NO_NETCDF */

end subroutine on_mdwrit

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

! Added by M Moradi
! Selection algorithm (a bootstrap-like resampling for multiple-walker simulations)
! Based on the method introduced in K. Minoukadeh et al, JCTC 2010, 6, 1008-1017.
! The differences w.r.t the original paper:
! * ABMD used instead of ABF
! * Biasing potential is used instead of histogram
! * Laplacian is used as a generalization of second derivative (for multi-dimensional cases)
! * Having an analytic function for biasing potential, the derivatives are exact
#ifdef NFE_ENABLE_BBMD
subroutine on_mdstep(x, v)

   use nblist
   use nfe_utils
   use nfe_colvar
   use nfe_umbrella
   use nfe_constants
   use nfe_sander_proxy
   use random, only : amrand

   implicit none

   NFE_REAL, intent(inout) :: x(*)
   NFE_REAL, intent(inout) :: v(*)

   character(len=*), parameter :: nullfile='/dev/null'

#  include "nfe-mpi.h"
   integer :: error
#  include "box.h"
#  include "ew_cntrl.h"
   NFE_REAL :: my_rand, w_norm, my_bias, my_laplacian, selection_criterion
   integer :: my_n, n_tot
   integer :: n

   if (((imode.eq.MODE_FLOODING).or.(imode.eq.MODE_UMBRELLA)).and.(selection_freq.gt.0) &
      .and.(multisander_numgroup().gt.1) &
      .and.(sander_init().eq.4)) then

      nfe_assert(ncolvars.gt.0)
      do n = 1, ncolvars
         instantaneous(n) = colvar_value(colvars(n), x)
      end do

      NFE_MASTER_ONLY_BEGIN
        call umbrella_eval_laplacian(umbrella,instantaneous(1:ncolvars),my_bias,my_laplacian)
        if (abs(my_bias) > TINY) &
           my_w = my_w + selection_constant * my_laplacian / my_bias

        if (mdstep.gt.1.and.mod(mdstep, selection_freq).eq.0) then
          if (multisander_numgroup().gt.1) then
            nfe_assert(commmaster.ne.mpi_comm_null)
            ! get all x
            nfe_assert(allocated(all_x))
            call mpi_allgather(x, 3*sander_natoms(), &
               MPI_DOUBLE_PRECISION, all_x, 3*sander_natoms(), &
               MPI_DOUBLE_PRECISION, commmaster,error)
            nfe_assert(error.eq.0)
            ! get all v
            nfe_assert(allocated(all_v))
            call mpi_allgather(v, 3*sander_natoms(), &
               MPI_DOUBLE_PRECISION, all_v, 3*sander_natoms(), &
               MPI_DOUBLE_PRECISION, commmaster,error)
            nfe_assert(error.eq.0)
            if (sander_ntb().ne.0) then
               nfe_assert(allocated(all_box))
               call mpi_allgather(box, 3, &
                  MPI_DOUBLE_PRECISION, all_box, 3, &
                  MPI_DOUBLE_PRECISION, commmaster,error)
               nfe_assert(error.eq.0)
            end if ! sander_ntb().ne.0
            my_w = exp(my_w)
            if (my_w > 1.0d+20) &
               call fatal('Selection weight overflows, try smaller selection constant')
            nfe_assert(allocated(all_w))
            call mpi_allgather(my_w, 1, &
               MPI_DOUBLE_PRECISION, all_w, 1, &
               MPI_DOUBLE_PRECISION, commmaster,error)
            nfe_assert(error.eq.0)
            call mpi_barrier(commmaster,error)
            allocate(all_n(multisander_numgroup()), stat = error)
            if (masterrank.eq.0) then
               call amrand_gen(nfe_abmd_gen, my_rand)
               n_tot = 0
               nfe_assert(allocated(all_n))
               all_n = 0
               w_norm = sum(all_w)
               if (abs(w_norm).gt.0) then
                  selection_criterion = ZERO
                  do n = 1, multisander_numgroup()
                    my_n = int(my_rand+multisander_numgroup()*sum(all_w(1:n))/w_norm)
                    selection_criterion = selection_criterion + (all_w(n)/w_norm) &
                                        * log(multisander_numgroup()*all_w(n)/w_norm)
                    write (unit= LOG_UNIT , fmt =&
                      '(a,a,'//pfmt(n)//',a,'//pfmt(all_w(n),6)//',a,'&
                          //pfmt(w_norm,6)//',a,f5.3,a,'//pfmt(my_n-n_tot)//',a)')&
                          NFE_INFO,'selection score for walker ',n,' is ', &
                          all_w(n),' / ',w_norm,' = ',all_w(n)/w_norm,&
                          ' => ',my_n-n_tot,' walker(s)'
                    if (my_n.gt.n_tot) then
                       all_n(n_tot+1:my_n) = n
                       n_tot = my_n
                    end if
                  end do                  
                  if (selection_epsilon.lt.ONE) then
                     nfe_assert(selection_criterion.ge.-TINY)
                     if (abs(selection_criterion).lt.TINY) &
                         selection_criterion = ZERO
                     if (selection_criterion.ge.selection_epsilon* &
                         log(NFE_TO_REAL(multisander_numgroup()))) then
                        write (unit= LOG_UNIT , fmt =&
                         '(a,a,'//pfmt(selection_criterion,6)//',a,'&
                          //pfmt(selection_epsilon*log(NFE_TO_REAL(multisander_numgroup())),6)//')')&
                          NFE_INFO,'Selection entropy ',selection_criterion,&
                          ' is greater than threshold ', &
                           selection_epsilon*log(NFE_TO_REAL(multisander_numgroup()))
                     else
                        write (unit= LOG_UNIT , fmt =&
                         '(a/,a,a,'//pfmt(selection_criterion,6)//',a,'&
                         //pfmt(selection_epsilon*log(NFE_TO_REAL(multisander_numgroup())),6)//'/,a,a/,a)')&
                         NFE_INFO,NFE_INFO,'selection entropy ',selection_criterion,&
                         ' is below the threshold ', &
                         selection_epsilon*log(NFE_TO_REAL(multisander_numgroup())),&
                         NFE_INFO,'selection process will be stopped.',NFE_INFO
                         selection_freq = 0
                     end if
                  end if
               end if
            end if
            call mpi_bcast(selection_freq, 1, MPI_INTEGER, 0, commmaster,error)
            nfe_assert(error.eq.0)
            call mpi_bcast(all_n, multisander_numgroup(), &
               MPI_INTEGER, 0, commmaster,error)
            nfe_assert(error.eq.0)
            call mpi_barrier(commmaster,error)
            my_w = ZERO

            if (masterrank+1.ne.all_n(masterrank+1)) then
              x(1:3*sander_natoms()) = &
                 all_x((all_n(masterrank+1)-1)*3*sander_natoms()+1:all_n(masterrank+1)*3*sander_natoms())
              v(1:3*sander_natoms()) = &
                 all_v((all_n(masterrank+1)-1)*3*sander_natoms()+1:all_n(masterrank+1)*3*sander_natoms())
              if (sander_ntb().ne.0) &
                 box(1:3) = &
                    all_box((all_n(masterrank+1)-1)*3+1:all_n(masterrank+1)*3)
              write (unit= LOG_UNIT , fmt ='(a,a,'//pfmt(masterrank+1)//',a,'//pfmt(all_n(masterrank+1))//'/)')&
                 NFE_INFO,'Selection resampling : new ',masterrank+1,' comes from ',all_n(masterrank+1)
            end if
            if (allocated(all_n)) deallocate(all_n)
          end if
        call flush(LOG_UNIT)
        end if
      NFE_MASTER_ONLY_END

      call mpi_bcast(selection_freq, 1, MPI_INTEGER, 0, commsander, error)
      nfe_assert(error.eq.0)
      if ((selection_freq.gt.0).and.mdstep.gt.1.and.mod(mdstep, selection_freq).eq.0) then
         call mpi_bcast(x(1:3*sander_natoms()), 3*sander_natoms(), MPI_DOUBLE_PRECISION, 0, commsander, error)
         nfe_assert(error.eq.0)
         call mpi_bcast(v(1:3*sander_natoms()), 3*sander_natoms(), MPI_DOUBLE_PRECISION, 0, commsander, error)
         nfe_assert(error.eq.0)
         if (sander_ntb().ne.0) then
            call mpi_bcast(box(1:3), 3, MPI_DOUBLE_PRECISION, 0, commsander, error)
            nfe_assert(error.eq.0)
            !open(unit=OUT_UNIT, file=nullfile, status='old')
            call fill_ucell(box(1),box(2),box(3),alpha,beta,gamma)
            !close(unit=OUT_UNIT)
            !call amopen(OUT_UNIT, MDOUT_FILE, 'O', 'F', 'A')
         end if
      end if
   end if

end subroutine on_mdstep
#endif /* NFE_ENABLE_BBMD */
! Moradi end

!><><><><><><><><><><><><><><><>< R E M D ><><><><><><><><><><><><><><><><><><

#ifdef MPI

!
! populates rem_* arrays (indexed by masterrank)
!

subroutine rem_postinit()

   use nfe_utils
   use nfe_umbrella
   use nfe_sander_proxy

   implicit none

#  include "nfe-mpi.h"

   integer :: ng, error, i

   nfe_assert(multisander_rem().ne.0)
   nfe_assert(commmaster.ne.mpi_comm_null)

   ng = multisander_numgroup() - 1

   select case(imode)
      case(MODE_ANALYSIS)
         allocate(rem_monitor_files(0:ng), &
                  rem_monitor_freqs(0:ng), stat = error)
      case(MODE_UMBRELLA)
         allocate(rem_monitor_files(0:ng), &
                  rem_monitor_freqs(0:ng), rem_umbrellas(0:ng), stat = error)
      case(MODE_FLOODING)
         allocate(rem_monitor_files(0:ng), &
                  rem_monitor_freqs(0:ng), rem_umbrella_files(0:ng), &
                  rem_snapshots_basenames(0:ng), rem_snapshots_freqs(0:ng), &
                  rem_umbrellas(0:ng), stat = error)
      case default
         nfe_assert_not_reached()
         continue
   end select

   if (error.ne.0) &
       NFE_OUT_OF_MEMORY

   do i = 0, ng

      if (masterrank.eq.i) then
         rem_monitor_files(i) = monitor_file
         rem_monitor_freqs(i) = monitor_freq
      end if ! masterrank.eq.i

      call mpi_bcast(rem_monitor_files(i), SL, &
         MPI_CHARACTER, i, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(rem_monitor_freqs(i), 1, &
         MPI_INTEGER, i, commmaster, error)
      nfe_assert(error.eq.0)

      if (imode.eq.MODE_FLOODING) then

         if (masterrank.eq.i) then
            rem_umbrella_files(i) = umbrella_file
            rem_snapshots_basenames(i) = snapshots_basename
            rem_snapshots_freqs(i) = snapshots_freq
         end if ! masterrank.eq.i

         call mpi_bcast(rem_umbrella_files(i), SL, &
            MPI_CHARACTER, i, commmaster, error)
         nfe_assert(error.eq.0)

         call mpi_bcast(rem_snapshots_basenames(i), SL, &
            MPI_CHARACTER, i, commmaster, error)
         nfe_assert(error.eq.0)

         call mpi_bcast(rem_snapshots_freqs(i), 1, &
            MPI_INTEGER, i, commmaster, error)
         nfe_assert(error.eq.0)

      end if ! FLOODING

      if (imode.eq.MODE_UMBRELLA.or.imode.eq.MODE_FLOODING) then
         if (masterrank.eq.i) &
            call umbrella_swap(rem_umbrellas(i), umbrella)

         call umbrella_bcast(rem_umbrellas(i), commmaster, i)

         if (masterrank.eq.i) &
            call umbrella_swap(rem_umbrellas(i), umbrella)
      end if

   end do

end subroutine rem_postinit

!-----------------------------------------------------------------------------

subroutine rem_cleanup()

   NFE_USE_AFAILED

   implicit none

#  include "nfe-mpi.h"

   nfe_assert(sanderrank.eq.0)

   if (imode.eq.MODE_NONE) then
      continue
      nfe_assert(.not.allocated(rem_umbrellas))
      nfe_assert(.not.allocated(rem_monitor_files))
   else if (imode.eq.MODE_ANALYSIS) then
      nfe_assert(allocated(rem_monitor_files))
      nfe_assert(allocated(rem_monitor_freqs))
      deallocate(rem_monitor_files, rem_monitor_freqs)
   else if (imode.eq.MODE_UMBRELLA) then
      nfe_assert(allocated(rem_monitor_files))
      nfe_assert(allocated(rem_monitor_freqs))
      nfe_assert(allocated(rem_umbrellas))
      call finalize_umbrellas()
      deallocate(rem_monitor_files, rem_monitor_freqs, rem_umbrellas)
   else if (imode.eq.MODE_FLOODING) then
      nfe_assert(allocated(rem_monitor_files))
      nfe_assert(allocated(rem_monitor_freqs))
      nfe_assert(allocated(rem_umbrella_files))
      nfe_assert(allocated(rem_snapshots_basenames))
      nfe_assert(allocated(rem_snapshots_freqs))
      nfe_assert(allocated(rem_umbrellas))
      call finalize_umbrellas()
      deallocate(rem_monitor_files, rem_monitor_freqs, &
                 rem_snapshots_basenames, rem_snapshots_freqs, &
                 rem_umbrella_files, rem_umbrellas)
   else
      continue
      nfe_assert_not_reached()
   end if

contains

subroutine finalize_umbrellas

   use nfe_umbrella, only : umbrella_fini
   use nfe_sander_proxy, only : multisander_numgroup

   implicit none

   integer :: i

   do i = 0, multisander_numgroup() - 1
      if (masterrank.ne.i) &
         call umbrella_fini(rem_umbrellas(i))
   end do

end subroutine finalize_umbrellas

end subroutine rem_cleanup

#endif /* MPI */

end module nfe_abmd_hooks
