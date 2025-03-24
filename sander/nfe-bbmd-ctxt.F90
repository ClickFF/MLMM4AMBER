#include "nfe-utils.h"
#include "nfe-config.h"


module nfe_bbmd_ctxt

#ifdef NFE_ENABLE_BBMD

use nfe_constants, only : SL => STRING_LENGTH, BBMD_MONITOR_UNIT, BBMD_CV_UNIT

use nfe_umbrella, only : umbrella_t, &
   MAX_NUMBER_OF_COLVARS => UMBRELLA_MAX_NEXTENTS

use nfe_colvar_type, only : colvar_t

implicit none

private

integer, private, parameter :: MONITOR_UNIT = BBMD_MONITOR_UNIT
integer, private, parameter :: CV_UNIT = BBMD_CV_UNIT

character(*), private, parameter :: SECTION = 'bbmd'

character(*), private, parameter :: &
   DEFAULT_MONITOR_FILE = 'nfe-bbmd-monitor', &
   DEFAULT_UMBRELLA_FILE = 'nfe-bbmd-umbrella', &
   DEFAULT_WT_UMBRELLA_FILE = 'nfe-bbmd-wt-umbrella', &
   DEFAULT_SNAPSHOTS_BASENAME = 'nfe-bbmd-umbrella-snapshot'

integer, private, parameter :: MODE_NONE = 5432

integer, private, parameter :: MODE_ANALYSIS = 1234
integer, private, parameter :: MODE_UMBRELLA = 2345
integer, private, parameter :: MODE_FLOODING = 3456

character(SL), public, save :: mode = 'NONE'
character(SL), public, save :: monitor_file ! master only
integer,       public, save :: monitor_freq = 50
NFE_REAL,      public, save :: timescale = 1 ! master only
character(SL), public, save :: umbrella_file ! master only
character(SL), public, save :: snapshots_basename ! master only
integer,       public, save :: snapshots_freq = -1 ! master only
NFE_REAL,      public, save :: wt_temperature = 0.0
character(SL), public, save :: wt_umbrella_file ! master only
character(SL), public, save :: cv_file = 'nfe-bbmd-cv'
character(SL), public, save :: driven_weight = 'NONE'
NFE_REAL,      public, save :: driven_cutoff = 0.0

integer, public, save :: exchange_freq
character(SL), public, save :: exchange_log_file
integer, public, save :: exchange_log_freq

integer, public, save :: mt19937_seed
character(SL), public, save :: mt19937_file

!-------------------------------------------------------------------------------

type, public :: bbmd_ctxt_t

   character(SL) :: mdout ! master only
   character(SL) :: restrt ! master only
   character(SL) :: mdvel ! master only
   character(SL) :: mden ! master only
   character(SL) :: mdcrd ! master only
   character(SL) :: mdinfo ! master only

   integer :: ioutfm
   integer :: ntpr
   integer :: ntwr
   integer :: ntwx

   character(SL) :: monitor_file ! master only
   character(SL) :: umbrella_file ! master only
   character(SL) :: wt_umbrella_file ! master only
   character(SL) :: snapshots_basename ! master only

   character(SL) :: monitor_fmt ! master only

   integer :: monitor_freq
   integer :: snapshots_freq ! master only

   NFE_REAL :: timescale ! master only

   integer :: imode
   integer :: ncolvars

   type(colvar_t)   :: colvars(MAX_NUMBER_OF_COLVARS)
   type(umbrella_t) :: umbrella ! master only (sanderrank.eq.0)
   type(umbrella_t) :: wt_umbrella ! master only (sanderrank.eq.0)

   logical :: umbrella_file_existed ! home-master only
   logical :: umbrella_discretization_changed  ! home-master only

! Added by M Moradi
! Well-tempered ABMD
   NFE_REAL :: pseudo
! Driven ABMD
   integer   :: drivenw
   integer   :: drivenu
   NFE_REAL  :: driven_cutoff
! Moradi end

#ifndef NFE_DISABLE_ASERT
   logical :: initialized = .false.
#endif /* NFE_DISABLE_ASSERT */

end type bbmd_ctxt_t

public :: ctxt_init
public :: ctxt_fini

public :: ctxt_print
public :: ctxt_bcast

public :: ctxt_on_force
public :: ctxt_on_mdwrit

public :: ctxt_Um
public :: ctxt_Uo

public :: ctxt_send
public :: ctxt_recv

public :: ctxt_close_units
public :: ctxt_open_units

!-------------------------------------------------------------------------------

NFE_REAL, parameter, private :: TINY = 0.000001d0

public ::  bbmd
namelist / bbmd /    mode, monitor_file, monitor_freq, timescale,&
                     umbrella_file, snapshots_basename, snapshots_freq,&
                     wt_temperature, wt_umbrella_file, cv_file, &
                     driven_weight, driven_cutoff, &
                     exchange_freq, exchange_log_file, exchange_log_freq, &
                     mt19937_seed, mt19937_file

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

subroutine ctxt_init(self, amass)

   use nfe_utils
   use nfe_colvar
   use nfe_colvar_type
   use nfe_constants
   use nfe_umbrella
   use nfe_sander_proxy
   use file_io_dat

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   NFE_REAL, intent(in) :: amass(*)

#  include "nfe-mpi.h"

   logical :: umbrella_file_exists
   type(umbrella_t) :: umbrella_from_file

   logical :: do_transfer
   integer :: n, error, ifind, i
   character(len = 80) :: buf 

   integer :: cv_extents(UMBRELLA_MAX_NEXTENTS)
   logical :: cv_periodicity(UMBRELLA_MAX_NEXTENTS)

   NFE_REAL :: cv_origin(UMBRELLA_MAX_NEXTENTS)
   NFE_REAL :: cv_spacing(UMBRELLA_MAX_NEXTENTS)

   NFE_REAL :: tmp

   nfe_assert(.not.self%initialized)

   if (sanderrank.gt.0) &
      goto 1

   ! store SANDER's filenames

   self%mdout  = mdout
   self%restrt = restrt
   self%mdvel  = mdvel
   self%mden   = mden
   self%mdcrd  = mdcrd
   self%mdinfo = mdinfo

   self%ioutfm = ioutfm
   self%ntpr = ntpr
   self%ntwr = ntwr
   self%ntwx = ntwx

   ! setup defaults

   write (unit = monitor_file, fmt = '(a,a,i3.3)') &
      DEFAULT_MONITOR_FILE, '-', (masterrank + 1)
   write (unit = umbrella_file, fmt = '(a,a,i3.3,a)') &
      DEFAULT_UMBRELLA_FILE, '-', (masterrank + 1), '.nc'
   write (unit = snapshots_basename, fmt = '(a,a,i3.3)') &
      DEFAULT_SNAPSHOTS_BASENAME, '-', (masterrank + 1)
   write (unit = wt_umbrella_file, fmt = '(a,a,i3.3,a)') &
      DEFAULT_WT_UMBRELLA_FILE, '-', (masterrank + 1), '.nc'             

   self%ncolvars = 0
   
   rewind(5)
   read(5,nml=bbmd,err=666)   
   ! discover the run-mode

   if (mode == 'NONE') then
      self%imode = MODE_NONE
      goto 1
   else if (mode == 'ANALYSIS') then
      self%imode = MODE_ANALYSIS
   else if (mode == 'UMBRELLA') then
      self%imode = MODE_UMBRELLA
   else if (mode == 'FLOODING') then
      self%imode = MODE_FLOODING
   else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'unknown mode ''', trim(mode), ''''
      call terminate()
   end if

! Added by M Moradi
! Well-tempered ABMD
   self%pseudo = ZERO
   if (wt_temperature .gt. ZERO) then
      self%pseudo = ONE/wt_temperature
   end if
   self%wt_umbrella_file = wt_umbrella_file
!  Driven ABMD
   if (driven_weight == 'NONE') then
      self%drivenw = 0
      self%drivenu = 0
   else if (driven_weight == 'CONSTANT') then
      self%drivenw = 1
      self%drivenu = 1
   else if (driven_weight == 'PULLING') then
      self%drivenw = 1
      self%drivenu = 0
   else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'unknown driven weight scheme ''', trim(driven_weight), ''''
      call terminate()
   end if

   self%driven_cutoff = driven_cutoff
! Moradi end

   self%umbrella_file = umbrella_file
#ifdef NFE_NO_NETCDF
   umbrella_file_exists = .false.
   write (unit = ERR_UNIT, fmt = '(a,a)') NFE_WARNING, &
      'netCDF is not available (try ''-bintraj'' configure option)'
#else
   inquire (file = self%umbrella_file, exist = umbrella_file_exists)
#endif /* NFE_NO_NETCDF */

   if (.not.umbrella_file_exists.and.self%imode.eq.MODE_UMBRELLA) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NFE_ERROR, '''', &
      trim(self%umbrella_file), ''' does not exist (required for UMBRELLA mode)'
      call terminate()
   end if

   if (self%imode.eq.MODE_ANALYSIS) &
      umbrella_file_exists = .false.

#ifndef NFE_NO_NETCDF
   if (umbrella_file_exists) &
      call umbrella_load(umbrella_from_file, self%umbrella_file)
#endif /* NFE_NO_NETCDF */

   ! collective variables
   nfe_assert(self%ncolvars.eq.0)
   
   call amopen(CV_UNIT, cv_file, 'O', 'F', 'R')
   
   do
     call nmlsrc('colvar', CV_UNIT, ifind)
     if (ifind.eq.0) exit
     read(CV_UNIT,'(a80)') buf
     self%ncolvars = self%ncolvars + 1
   end do   

   if (self%ncolvars.eq.0) &
      call fatal('no variable(s) in the '''//SECTION//''' section')

   if (self%ncolvars.gt.MAX_NUMBER_OF_COLVARS) &
      call fatal('too many variables in the '''//SECTION//''' section')

   if (umbrella_file_exists) then
      if(umbrella_nextents(umbrella_from_file).ne.self%ncolvars) &
         call fatal('number of variables in the '''//SECTION//''' does not &
                 &match with the number of extents found in the umbrella_file')
   end if ! umbrella_file_exists

   n = 1

   do while (n.le.self%ncolvars)

         call colvar_nlread(CV_UNIT, self%colvars(n))
         
         if (self%imode.eq.MODE_FLOODING) then
            cv_spacing(n) = resolution
            cv_spacing(n) = cv_spacing(n)/4
            if ((resolution.eq.ZERO).and..not.umbrella_file_exists) then
               write (unit = ERR_UNIT, fmt = '(/a,a,i1/)') NFE_ERROR, &
                  'could not determine ''resolution'' for CV #', n
               call terminate()
            end if

            if (resolution.eq.ZERO) &
               cv_spacing(n) = umbrella_spacing(umbrella_from_file, n)

            cv_periodicity(n) = colvar_is_periodic(self%colvars(n))

            if (cv_periodicity(n)) then

               nfe_assert(colvar_has_min(self%colvars(n)))
               nfe_assert(colvar_has_max(self%colvars(n)))

               cv_origin(n) = colvar_min(self%colvars(n))

               nfe_assert(cv_spacing(n).gt.ZERO)
               nfe_assert(colvar_max(self%colvars(n)).gt.cv_origin(n))

               cv_extents(n) = &
               int((colvar_max(self%colvars(n)) - cv_origin(n))/cv_spacing(n))

               if (cv_extents(n).lt.UMBRELLA_MIN_EXTENT) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1,a/)') NFE_ERROR, &
                     'CV #', n, ' : ''resolution'' is too big'
                  call terminate()
               end if

               cv_spacing(n) = &
                  (colvar_max(self%colvars(n)) - cv_origin(n))/cv_extents(n)

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
         end if ! mode.eq.MODE_FLOODING
         if (colvar_has_refcrd(self%colvars(n))) then
            refcrd_file = trim(refcrd_file)
            refcrd_len = len_trim(refcrd_file)         
         end if
         if (colvar_is_quaternion(self%colvars(n))) then 
          allocate(self%colvars(n)%q_index, stat = error)
           if (error.ne.0) &
            NFE_OUT_OF_MEMORY
            self%colvars(n)%q_index = q_index
         end if 
         if (colvar_has_axis(self%colvars(n))) then
           allocate(self%colvars(n)%axis(3), stat = error)
             if (error.ne.0) &
                NFE_OUT_OF_MEMORY
           i = 1
           do while (i.le.3)
             self%colvars(n)%axis(i) = axis(i)
             i = i + 1
           end do
         end if

         n = n + 1
   end do

   ! monitor
   self%monitor_file = monitor_file

   self%monitor_freq = monitor_freq

   self%monitor_freq = min(self%monitor_freq, sander_nstlim())
   self%monitor_freq = max(1, self%monitor_freq)

   ! umbrella snapshots
   self%snapshots_basename = snapshots_basename

   self%snapshots_freq = snapshots_freq

   if (self%imode.eq.MODE_FLOODING) then
      self%timescale = timescale
      if (self%timescale.eq.ZERO) &
         call fatal('timescale cannot be zero !')
   end if ! mode.eq.MODE_FLOODING

1  continue ! sanderrank.gt.0 jumps here

   call mpi_bcast(self%imode, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   if (self%imode.eq.MODE_NONE) &
      goto 2

   nfe_assert(.not.is_master().or.self%ncolvars.gt.0)

   call mpi_bcast(self%ncolvars, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(self%monitor_freq, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)

   nfe_assert(self%ncolvars.gt.0)
   nfe_assert(self%ncolvars.le.MAX_NUMBER_OF_COLVARS)
 
   call mpi_bcast(refcrd_len, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(refcrd_file, refcrd_len, MPI_CHARACTER, 0, commsander, error)
   nfe_assert(error.eq.0)

   if (multisander_numwatkeep().gt.0) &
      call fatal('numwatkeep.gt.0 is not supported')

   if (sander_imin().ne.0) &
      call fatal('imin.ne.0 is not supported')

   do n = 1, self%ncolvars
      call colvar_bootstrap(self%colvars(n), n, amass)
   end do

   if (sanderrank.gt.0) &
      goto 2

   do_transfer = .false.

   if (self%imode.eq.MODE_UMBRELLA) then
      nfe_assert(umbrella_file_exists)
      do_transfer = .false.
      call umbrella_swap(self%umbrella, umbrella_from_file)
   else if (self%imode.eq.MODE_FLOODING) then
      if (umbrella_file_exists) then
         do_transfer = .false.
         do n = 1, self%ncolvars
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
            call umbrella_init(self%umbrella, self%ncolvars, cv_extents, &
                               cv_origin, cv_spacing, cv_periodicity)
            call umbrella_transfer(self%umbrella, umbrella_from_file)
            call umbrella_fini(umbrella_from_file)
         else
            call umbrella_swap(self%umbrella, umbrella_from_file)
         end if ! do_transfer
      else
         call umbrella_init(self%umbrella, self%ncolvars, cv_extents, &
                            cv_origin, cv_spacing, cv_periodicity)
      end if ! umbrella_file_exits
      call umbrella_init(self%wt_umbrella, self%ncolvars, cv_extents, &
                            cv_origin, cv_spacing, cv_periodicity)
   end if ! self%mode.eq.MODE_FLOODING

   self%umbrella_file_existed = umbrella_file_exists
   self%umbrella_discretization_changed = do_transfer

   ! prepare monitor_fmt & open MONITOR_UNIT

   open (unit = MONITOR_UNIT, file = self%monitor_file, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

   if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'failed to open ''', trim(self%monitor_file), &
         ''' for writing'
      call terminate()
   end if

   write (unit = MONITOR_UNIT, fmt = '(a,/a)', advance = 'NO') &
      '#', '# MD time (ps), '
   do n = 1, self%ncolvars - 1
      write (unit = MONITOR_UNIT, fmt = '(a,i1,a)', advance = 'NO') &
         'CV #', n, ', '
   end do

   if (self%imode == MODE_FLOODING) then
      write (unit = MONITOR_UNIT, fmt = '(a,i1,a,/a)') &
         'CV #', self%ncolvars, ', E_{bias} (kcal/mol)', '#'
      write (unit = self%monitor_fmt, fmt = '(a,i1,a)') &
         '(f12.4,', self%ncolvars, '(1x,f16.10),1x,f16.10)'
   else
      write (unit = MONITOR_UNIT, fmt = '(a,i1,/a)') &
         'CV #', self%ncolvars, '#'
      write (unit = self%monitor_fmt, fmt = '(a,i1,a)') &
         '(f12.4,', self%ncolvars, '(1x,f16.10))'
   end if

   call flush_UNIT(MONITOR_UNIT)

2  continue ! sanderrank.gt.0 jump here (or mode.eq.MODE_NONE)

#  ifndef NFE_DISABLE_ASSERT
   self%initialized = .true.
#  endif /* NFE_DISABLE_ASSERT */

   return
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &bbmd namelist!'
    call terminate() 

end subroutine ctxt_init

!-------------------------------------------------------------------------------

subroutine ctxt_fini(self)

   NFE_USE_AFAILED

   use nfe_colvar
   use nfe_umbrella
   use nfe_sander_proxy

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

#  include "nfe-mpi.h"

   integer :: n

   nfe_assert(self%initialized)

   if (self%imode.ne.MODE_NONE) then
      nfe_assert(self%ncolvars.gt.0)
      do n = 1, self%ncolvars
         call colvar_cleanup(self%colvars(n))
      end do
   end if

   if (sanderrank.eq.0) then
      if (self%imode.eq.MODE_FLOODING.or.self%imode.eq.MODE_UMBRELLA) &
         call umbrella_fini(self%umbrella)
   end if ! sanderrank.eq.0

   self%imode = MODE_NONE

#  ifndef NFE_DISABLE_ASSERT
   self%initialized = .false.
#  endif /* NFE_DISABLE_ASSERT */

end subroutine ctxt_fini

!-------------------------------------------------------------------------------

subroutine ctxt_print(self, lun)

   use nfe_utils
   use nfe_colvar
   use nfe_umbrella
   use nfe_constants
   use nfe_sander_proxy

   implicit none

   type(bbmd_ctxt_t), intent(in) :: self
   integer, intent(in) :: lun

   integer :: n
   NFE_REAL :: tmp

   nfe_assert(self%initialized)

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'mode = '

   select case(self%imode)
      case(MODE_NONE)
         write (unit = lun, fmt = '(a)') 'NONE'
         goto 1
      case(MODE_ANALYSIS)
         write (unit = lun, fmt = '(a)') 'ANALYSIS'
      case(MODE_UMBRELLA)
         write (unit = lun, fmt = '(a)') 'UMBRELLA'
      case(MODE_FLOODING)
         write (unit = lun, fmt = '(a)') 'FLOODING'
      case default
         nfe_assert_not_reached()
         continue
   end select

   write (unit = lun, fmt = '(a)') NFE_INFO

   do n = 1, self%ncolvars
      write (unit = lun, fmt = '(a,a,i1)') NFE_INFO, 'CV #', n
      call colvar_print(self%colvars(n), lun)
      write (unit = lun, fmt = '(a)') NFE_INFO
      if (colvar_is_quaternion(self%colvars(n))) then
         write (unit = lun, fmt = '(a,a,I3)') NFE_INFO, &
         ' q_index = ', self%colvars(n)%q_index
      end if
      if (colvar_has_axis(self%colvars(n))) then
         write (unit = lun, fmt = '(a,a,f8.4,a,f8.4,a,f8.4,a)') NFE_INFO, &
         ' axis = [', self%colvars(n)%axis(1),', ', self%colvars(n)%axis(2),',', self%colvars(n)%axis(3),']'
      end if
   end do

   write (unit = lun, fmt = '(a,a,a)') NFE_INFO, &
      'monitor_file = ', trim(self%monitor_file)
   write (unit = lun, fmt = '(a,a,'//pfmt &
      (self%monitor_freq)//',a,'//pfmt &
      (self%monitor_freq*sander_timestep(), 4)//',a)') NFE_INFO, &
      'monitor_freq = ', self%monitor_freq, ' (', &
      self%monitor_freq*sander_timestep(), ' ps)'

   if (self%imode.eq.MODE_ANALYSIS) &
      goto 1

   write (unit = lun, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
      'umbrella_file = ', trim(self%umbrella_file), ' ('

   if (self%umbrella_file_existed) then
      write (unit = lun, fmt = '(a)') 'loaded)'
   else
      write (unit = lun, fmt = '(a)') 'not found)'
   end if

   write (unit = lun, fmt = '(a)') NFE_INFO
   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, &
      'umbrella discretization '

   if (self%umbrella_file_existed) then
      if (self%umbrella_discretization_changed) then
         write (unit = lun, fmt = '(a)') '(modified) :'
      else
         write (unit = lun, fmt = '(a)') '(unchanged) :'
      end if
   else
      write (unit = lun, fmt = '(a)') '(new) :'
   end if

   do n = 1, self%ncolvars
      write (unit = lun, fmt = '(a,a,i1)', advance = 'NO') &
         NFE_INFO, 'CV #', n
      if (umbrella_periodicity(self%umbrella, n)) then
         write (unit = lun, fmt = '(a)', advance = 'NO') ' periodic, '
         tmp = umbrella_origin(self%umbrella, n) &
         + umbrella_spacing(self%umbrella, n)*umbrella_extent(self%umbrella, n)
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ' not periodic, '
         tmp = umbrella_origin(self%umbrella, n) &
         + umbrella_spacing(self%umbrella, n)*(umbrella_extent(self%umbrella, n) - 1)
      end if

      write (unit = lun, &
         fmt = '('//pfmt(umbrella_extent(self%umbrella, n))//',a,'//pfmt &
         (umbrella_origin(self%umbrella, n), 6)//',a,'//pfmt(tmp, 6)//')') &
         umbrella_extent(self%umbrella, n), ' points, min/max = ', &
         umbrella_origin(self%umbrella, n), '/', tmp
   end do

   if (self%imode.eq.MODE_UMBRELLA) &
      goto 1

   write (unit = lun, fmt = '(a/,a,a,'//pfmt(self%timescale, 3)//',a)') &
      NFE_INFO, NFE_INFO, 'flooding timescale = ', self%timescale, ' ps'

   if (self%snapshots_freq.gt.0) then
      write (unit = lun, fmt = '(a,a,a)') NFE_INFO, &
         'snapshots_basename = ', trim(self%snapshots_basename)
      write (unit = lun, &
         fmt = '(a,a,'//pfmt(self%snapshots_freq)//',a,'//pfmt &
         (self%snapshots_freq*sander_timestep(), 4)//',a)') NFE_INFO, &
         'snapshots_freq = ', self%snapshots_freq, ' (', &
         self%snapshots_freq*sander_timestep(), ' ps)'
   end if

! Modified by M Moradi
! Well-tempered ABMD
   if (self%pseudo.gt.ZERO) then
      write (unit = lun, &
      fmt = '(a/,a,a)')&
          NFE_INFO, NFE_INFO,'well-tempered ABMD:'
      write (unit = lun, &
      fmt = '(a,a,'//pfmt(ONE/self%pseudo,6)//')')&
          NFE_INFO, 'pseudo-temperature = ', ONE/self%pseudo
      write (unit = lun, fmt = '(a,a,a)') NFE_INFO, &
         'wt_umbrella_file = ', trim(self%wt_umbrella_file)          
   end if
! Driven ABMD
   if (self%drivenw.gt.ZERO) then
      write (unit = lun, &
      fmt = '(a/,a,a)')&
          NFE_INFO, NFE_INFO,'driven ABMD:'
      if (self%drivenu.gt.ZERO) then
         write (unit = lun, &
         fmt = '(a,a)')&
             NFE_INFO,'CONSTANT weighting scheme (use delta work)'
      else
         write (unit = lun, &
         fmt = '(a,a)')&
             NFE_INFO,'PULLING weighting scheme (use work only)'
      end if
      write (unit = lun, &
      fmt = '(a,a,'//pfmt(self%driven_cutoff,6)//')')&
          NFE_INFO, 'driven (delta)work cutoff = ', self%driven_cutoff
   end if
! Moradi end

1  call flush_UNIT(lun)

end subroutine ctxt_print

!-------------------------------------------------------------------------------

subroutine ctxt_bcast(self, masterroot, amass)

   use nfe_utils
   use nfe_umbrella, only : umbrella_bcast

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   integer, intent(in) :: masterroot
   NFE_REAL, intent(in) :: amass(*)

#  include "nfe-mpi.h"

   integer :: n, error

   nfe_assert(masterroot.ge.0.and.masterroot.lt.mastersize)

#  ifndef NFE_DISABLE_ASSERT
   if (masterroot.eq.masterrank) then
      nfe_assert(self%initialized)
   else
      self%initialized = .true.
   end if ! masterroot.eq.masterrank
#  endif /* NFE_DISABLE_ASSERT */

   ! SANDER's files (no matter what the mode is)

   if (sanderrank.eq.0) then
      nfe_assert(commmaster.ne.MPI_COMM_NULL)

      call mpi_bcast(self%mdout, len(self%mdout), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%restrt, len(self%restrt), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%mdvel, len(self%mdvel), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%mden, len(self%mden), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%mdcrd, len(self%mdcrd), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%mdinfo, len(self%mdinfo), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%ioutfm, 1, MPI_INTEGER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%ntpr, 1, MPI_INTEGER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%ntwr, 1, MPI_INTEGER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%ntwx, 1, MPI_INTEGER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%imode, 1, MPI_INTEGER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

   end if ! sanderrank.eq.0

   if (masterrank.ne.masterroot) then
      call mpi_bcast(self%imode, 1, MPI_INTEGER, 0, commsander, error)
      nfe_assert(error.eq.0)
   end if ! masterrank.ne.masterroot

   if (self%imode.eq.MODE_NONE) &
      return

   ! CVs

   nfe_assert(self%ncolvars.gt.0.or.masterrank.ne.masterroot)

   if (sanderrank.eq.0) then
      call mpi_bcast(self%ncolvars, 1, &
                     MPI_INTEGER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)
   end if ! sanderrank.eq.0

   if (masterrank.ne.masterroot) then
      call mpi_bcast(self%ncolvars, 1, MPI_INTEGER, 0, commsander, error)
      nfe_assert(error.eq.0)
   end if ! masterrank.ne.masterroot

   nfe_assert(self%ncolvars.gt.0)

   do n = 1, self%ncolvars
      call bcast_colvar(self%colvars(n), n + 10*masterrank)
   end do

   ! mode = ANALYSIS

   if (sanderrank.eq.0) then
      call mpi_bcast(self%monitor_file, len(self%monitor_file), &
                     MPI_CHARACTER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%monitor_fmt, len(self%monitor_fmt), &
                     MPI_CHARACTER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%monitor_freq, 1, &
                     MPI_INTEGER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)
   end if ! sanderrank.eq.0

   if (masterrank.ne.masterroot) then
      call mpi_bcast(self%monitor_freq, 1, MPI_INTEGER, 0, commsander, error)
      nfe_assert(error.eq.0)
   end if ! masterrank.ne.masterroot

   if (self%imode.eq.MODE_ANALYSIS) &
      return

   ! mode = UMBRELLA | FLOODING (these are on masters only)

   if (sanderrank.eq.0) then
      call mpi_bcast(self%umbrella_file, len(self%umbrella_file), &
                     MPI_CHARACTER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%snapshots_basename, len(self%snapshots_basename), &
                     MPI_CHARACTER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%snapshots_freq, 1, &
                     MPI_INTEGER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%timescale, 1, &
                     MPI_DOUBLE_PRECISION, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      call umbrella_bcast(self%umbrella, commmaster, masterroot)

      if (self%imode.eq.MODE_FLOODING) then
         call mpi_bcast(self%wt_umbrella_file, len(self%wt_umbrella_file), &
                     MPI_CHARACTER, masterroot, commmaster, error)
         nfe_assert(error.eq.0)      
      
         call mpi_bcast(self%pseudo, 1, &
                     MPI_DOUBLE_PRECISION, masterroot, commmaster, error)
         nfe_assert(error.eq.0) 

         call mpi_bcast(self%drivenu, 1, &
                     MPI_INTEGER, masterroot, commmaster, error)
         nfe_assert(error.eq.0)

         call mpi_bcast(self%drivenw, 1, &
                     MPI_INTEGER, masterroot, commmaster, error)
         nfe_assert(error.eq.0)

         call mpi_bcast(self%driven_cutoff, 1, &
                     MPI_DOUBLE_PRECISION, masterroot, commmaster, error)
         nfe_assert(error.eq.0)     

         call umbrella_bcast(self%wt_umbrella, commmaster, masterroot)
       end if

   end if ! sanderrank.eq.0

contains

subroutine bcast_colvar(cv, cvno)

   use nfe_colvar, only : colvar_bootstrap

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer, intent(in) :: cvno

   integer :: bcastdata(3)

   if (sanderrank.eq.0) then

      if (masterrank.eq.masterroot) then

         nfe_assert(cv%type.gt.0)

         bcastdata(1) = cv%type

         bcastdata(2) = 0
         if (associated(cv%i)) &
            bcastdata(2) = size(cv%i)

         bcastdata(3) = 0
         if (associated(cv%r)) &
            bcastdata(3) = size(cv%r)

      end if ! masterrank.eq.masterroot

      call mpi_bcast(bcastdata, 3, MPI_INTEGER, masterroot, commmaster, error)
      nfe_assert(error.eq.0)

      if (masterrank.ne.masterroot) then
         cv%type = bcastdata(1)

         if (bcastdata(2).gt.0) then
            allocate(cv%i(bcastdata(2)), stat = error)
            if (error.ne.0) &
               NFE_OUT_OF_MEMORY
         end if ! bcastdata(2).gt.0

         if (bcastdata(3).gt.0) then
            allocate(cv%r(bcastdata(3)), stat = error)
            if (error.ne.0) &
               NFE_OUT_OF_MEMORY
         end if ! bcastdata(3).gt.0

      end if ! masterrank.ne.masterroot

      if (bcastdata(2).gt.0) then
         call mpi_bcast(cv%i, bcastdata(2), MPI_INTEGER, &
                        masterroot, commmaster, error)
         nfe_assert(error.eq.0)
      end if ! bcastdata(2).gt.0

      if (bcastdata(3).gt.0) then
         call mpi_bcast(cv%r, bcastdata(3), MPI_DOUBLE_PRECISION, &
                        masterroot, commmaster, error)
         nfe_assert(error.eq.0)
      end if ! bcastdata(3).gt.0

   end if ! sanderrank.eq.0

   if (masterrank.ne.masterroot) &
      call colvar_bootstrap(cv, cvno, amass)

end subroutine bcast_colvar

end subroutine ctxt_bcast

!-------------------------------------------------------------------------------

! Modified by M Moradi
! Driven ABMD
subroutine ctxt_on_force(self, x, f, mdstep, wdriven, udriven, pot)
! Moradi end

   use nfe_utils
   use nfe_colvar
   use nfe_colvar_type
   use nfe_umbrella
   use nfe_constants
   use nfe_sander_proxy


   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   NFE_REAL, intent(in) :: x(*)

   NFE_REAL, intent(inout) :: f(*)

   integer, intent(in) :: mdstep
   
! Modified by M Moradi
! Driven ABMD
   NFE_REAL, intent(in) :: wdriven
   NFE_REAL, intent(in) :: udriven
! Moradi end
  
   NFE_REAL, intent(inout) :: pot

#  include "nfe-mpi.h"

   NFE_REAL :: u_derivative(UMBRELLA_MAX_NEXTENTS)
   NFE_REAL :: instantaneous(UMBRELLA_MAX_NEXTENTS)
   NFE_REAL :: alt, u_value

   character(len = SL + 16) :: snapshot

   integer :: n, error, m 
   integer, DIMENSION(4) :: cv_q = (/COLVAR_QUATERNION0, COLVAR_QUATERNION1, &
                                     COLVAR_QUATERNION2, COLVAR_QUATERNION3/)
   NFE_REAL :: norm4(100), cv_N(4, 100)

   nfe_assert(self%initialized)

   if (self%imode.eq.MODE_NONE) &
      return

   nfe_assert(self%ncolvars.gt.0)

   if (self%imode.eq.MODE_ANALYSIS) then
      if (nfe_real_mdstep.and.mod(mdstep, self%monitor_freq).eq.0) then
         do n = 1, self%ncolvars
            instantaneous(n) = colvar_value(self%colvars(n), x)
         end do
         NFE_MASTER_ONLY_BEGIN
         do n = 1, self%ncolvars
          if (colvar_is_quaternion(self%colvars(n))) then
           do m = 1, 4
            if (self%colvars(n)%type == cv_q(m)) then
              cv_N(m, self%colvars(n)%q_index) = instantaneous(n)
            end if
           end do
          end if
         end do
         do n = 1, self%ncolvars
          if (colvar_is_quaternion(self%colvars(n))) then
            norm4(self%colvars(n)%q_index) = sqrt(cv_N(1,self%colvars(n)%q_index)**2 &
                                                + cv_N(2,self%colvars(n)%q_index)**2 &
                                                + cv_N(3,self%colvars(n)%q_index)**2 &
                                                + cv_N(4,self%colvars(n)%q_index)**2)
          end if
         end do
         do n = 1, self%ncolvars
          if (colvar_is_quaternion(self%colvars(n))) then
             instantaneous(n) = instantaneous(n) / norm4(self%colvars(n)%q_index)
          else
             instantaneous(n) = instantaneous(n)
          end if
         end do
         NFE_MASTER_ONLY_END

         if (sanderrank.eq.0) then
            write (unit = MONITOR_UNIT, fmt = self%monitor_fmt) &
               sander_mdtime(), instantaneous(1:self%ncolvars)
            call flush_UNIT(MONITOR_UNIT)
         end if ! sanderrank.eq.0
      end if

      return
   end if ! self%imode.eq.MODE_ANALYSIS

   !
   ! either UMBRELLA or FLOODING
   !

   do n = 1, self%ncolvars
      instantaneous(n) = colvar_value(self%colvars(n), x)
   end do
   NFE_MASTER_ONLY_BEGIN
   do n = 1, self%ncolvars
    if (colvar_is_quaternion(self%colvars(n))) then
      do m = 1, 4
        if (self%colvars(n)%type == cv_q(m)) then
          cv_N(m, self%colvars(n)%q_index) = instantaneous(n)
        end if
      end do
    end if
   end do
   do n = 1, self%ncolvars
    if (colvar_is_quaternion(self%colvars(n))) then
       norm4(self%colvars(n)%q_index) =sqrt(cv_N(1,self%colvars(n)%q_index)**2 &
                                          + cv_N(2,self%colvars(n)%q_index)**2 &
                                          + cv_N(3,self%colvars(n)%q_index)**2 &
                                          + cv_N(4,self%colvars(n)%q_index)**2)
    end if
   end do
   do n = 1, self%ncolvars
    if (colvar_is_quaternion(self%colvars(n))) then
      instantaneous(n) = instantaneous(n) / norm4(self%colvars(n)%q_index)
    else
      instantaneous(n) = instantaneous(n)
    end if
   end do
   NFE_MASTER_ONLY_END

   if (sanderrank.eq.0) then
      call umbrella_eval_vdv(self%umbrella, instantaneous, &
                             u_value, u_derivative)
      pot = umbrella_eval_v(self%umbrella, instantaneous)
   end if

   call mpi_bcast(u_derivative, self%ncolvars, &
      MPI_DOUBLE_PRECISION, 0, commsander, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(pot, 1, &
      MPI_DOUBLE_PRECISION, 0, commsander, error)
   nfe_assert(error.eq.0)

   ! FIXME: virial
   do n = 1, self%ncolvars
      call colvar_force(self%colvars(n), x, -u_derivative(n), f)
   end do

   if (.not.nfe_real_mdstep.or.sanderrank.ne.0) &
      return

   nfe_assert(self%imode.eq.MODE_UMBRELLA.or.self%imode.eq.MODE_FLOODING)

   if (mod(mdstep, self%monitor_freq).eq.0) then
      if (self%imode.eq.MODE_FLOODING) then
         write (unit = MONITOR_UNIT, fmt = self%monitor_fmt) &
            sander_mdtime(), instantaneous(1:self%ncolvars), u_value
      else
         write (unit = MONITOR_UNIT, fmt = self%monitor_fmt) &
            sander_mdtime(), instantaneous(1:self%ncolvars)
      end if
      call flush_UNIT(MONITOR_UNIT)
   end if

   if (self%imode.eq.MODE_FLOODING) then
#  ifndef NFE_NO_NETCDF
      if (self%snapshots_freq.gt.0 &
          .and.mod(mdstep, self%snapshots_freq).eq.0) then
         write (unit = snapshot, fmt = '(a,a,i10.10,a)') &
            trim(self%snapshots_basename), '.', mdstep, '.nc'
         call umbrella_save(self%umbrella, snapshot)
         write (unit = OUT_UNIT, fmt = '(/a,a,'//pfmt &
            (sander_mdtime(), 4)//',a,/a,a,a,a)') &
            NFE_INFO, 'biasing potential snapshot at t = ', &
            sander_mdtime(), ' ps', NFE_INFO, 'saved as ''', &
            trim(snapshot), ''''
      end if
#  endif /* NFE_NO_NETCDF */
! Modified by M Moradi
      alt = sander_timestep()/self%timescale
! Well-tempered ABMD (based on Barducci, Bussi, and Parrinello, PRL(2008) 100:020603)
      if (self%pseudo.gt.ZERO) &
         alt = alt * &
         exp(-self%pseudo*umbrella_eval_v(self%umbrella,instantaneous)/kB)
! Driven ABMD (based on Moradi and Tajkhorshid, JPCL(2013) 4:1882)
      if (self%drivenw.ne.0) then
         if (self%drivenw*wdriven-self%drivenu*udriven.gt.self%driven_cutoff) then
            alt = alt * &
             exp(-(self%drivenw*wdriven-self%drivenu*udriven)/(kB*sander_temp0()))
         else
            alt = alt * &
             exp(-self%driven_cutoff/(kB*sander_temp0()))
         end if
      end if ! self%drivenw.ne.0
! Moradi end
      call umbrella_hill(self%umbrella, instantaneous, alt)
   end if ! self%imode.eq.MODE_FLOODING

end subroutine ctxt_on_force

!-------------------------------------------------------------------------------

subroutine ctxt_on_mdwrit(self)

   use nfe_utils
   use nfe_umbrella
   use nfe_sander_proxy

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   double precision :: t

   nfe_assert(self%initialized)

#ifndef NFE_NO_NETCDF
   if (self%imode.eq.MODE_FLOODING) then
      nfe_assert(is_master())
      call umbrella_save(self%umbrella, self%umbrella_file)
! Modified by F Pan 
      if (self%pseudo.gt.TINY) then
         t = sander_temp0()
         call umbrella_copy(self%wt_umbrella, self%umbrella)
         call umbrella_wt_mod(self%wt_umbrella,t,1/self%pseudo)
         call umbrella_save(self%wt_umbrella, self%wt_umbrella_file)
      end if
! Pan end      
   end if
#endif /* NFE_NO_NETCDF */

end subroutine ctxt_on_mdwrit

!-------------------------------------------------------------------------------

!
! U_m? are valid for masterrank.lt.r_masterrank
!

! assumes that local self is up to date
subroutine ctxt_Um(self, r_masterrank, x, U_mm, U_mo)

   NFE_USE_AFAILED

   use nfe_colvar
   use nfe_umbrella
   use nfe_constants

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   integer, intent(in) :: r_masterrank
   NFE_REAL, intent(in) :: x(*)

   NFE_REAL, intent(out) :: U_mm, U_mo

#  include "nfe-mpi.h"

   NFE_REAL :: local_inst(UMBRELLA_MAX_NEXTENTS)
   NFE_REAL :: remote_inst(UMBRELLA_MAX_NEXTENTS)

   integer :: n, error
   nfe_assert(self%initialized)

   U_mm = ZERO
   U_mo = ZERO

   if (self%imode.eq.MODE_NONE.or.self%imode.eq.MODE_ANALYSIS) &
      return

   nfe_assert(self%ncolvars.gt.0)
   do n = 1, self%ncolvars
      local_inst(n) = colvar_value(self%colvars(n), x)
   end do

   if (sanderrank.gt.0) &
      return

   nfe_assert(self%imode.eq.MODE_UMBRELLA.or.self%imode.eq.MODE_FLOODING)

   if (masterrank.lt.r_masterrank) then
      call mpi_recv(remote_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                    r_masterrank, 0, commmaster, MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
      U_mm = umbrella_eval_v(self%umbrella, local_inst)
      U_mo = umbrella_eval_v(self%umbrella, remote_inst)
   else
      call mpi_send(local_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                    r_masterrank, 0, commmaster, error)
      nfe_assert(error.eq.0)
   end if ! masterrank.lt.r_masterrank

end subroutine ctxt_Um

! assumes that remote self is up to date (if mode.eq.MODE_FLOODING)
subroutine ctxt_Uo(self, r_masterrank, x, U_om, U_oo)

   NFE_USE_AFAILED

   use nfe_colvar
   use nfe_umbrella
   use nfe_constants

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   integer, intent(in) :: r_masterrank
   NFE_REAL, intent(in) :: x(*)

   NFE_REAL, intent(out) :: U_om, U_oo

#  include "nfe-mpi.h"

   NFE_REAL :: local_inst(UMBRELLA_MAX_NEXTENTS)
   NFE_REAL :: remote_inst(UMBRELLA_MAX_NEXTENTS)

   NFE_REAL :: tmp(2)

   integer :: n, error

   nfe_assert(self%initialized)

   U_om = ZERO
   U_oo = ZERO

   if (self%imode.eq.MODE_NONE.or.self%imode.eq.MODE_ANALYSIS) &
      return

   nfe_assert(self%ncolvars.gt.0)
   do n = 1, self%ncolvars
      local_inst(n) = colvar_value(self%colvars(n), x)
   end do

   if (sanderrank.gt.0) &
      return

   if (self%imode.eq.MODE_UMBRELLA) then
      ! remote umbrella is same as local
      if (masterrank.lt.r_masterrank) then
         call mpi_recv(remote_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 0, commmaster, MPI_STATUS_IGNORE, error)
         nfe_assert(error.eq.0)
         U_om = umbrella_eval_v(self%umbrella, local_inst)
         U_oo = umbrella_eval_v(self%umbrella, remote_inst)
      else
         call mpi_send(local_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 0, commmaster, error)
         nfe_assert(error.eq.0)
      end if ! masterrank.lt.r_masterrank
   else
      nfe_assert(self%imode.eq.MODE_FLOODING)
      if (masterrank.gt.r_masterrank) then
         call mpi_recv(remote_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 0, commmaster, MPI_STATUS_IGNORE, error)
         nfe_assert(error.eq.0)
         tmp(1) = umbrella_eval_v(self%umbrella, local_inst)
         tmp(2) = umbrella_eval_v(self%umbrella, remote_inst)
         call mpi_send(tmp, 2, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 1, commmaster, error)
         nfe_assert(error.eq.0)
      else
         call mpi_send(local_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 0, commmaster, error)
         nfe_assert(error.eq.0)
         call mpi_recv(tmp, 2, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 1, commmaster, MPI_STATUS_IGNORE, error)
         nfe_assert(error.eq.0)
         U_om = tmp(2)
         U_oo = tmp(1)
      end if ! masterrank.gt.r_masterrank
   end if ! self%imode.eq.MODE_UMBRELLA

end subroutine ctxt_Uo

!-------------------------------------------------------------------------------

subroutine ctxt_send(self, dst_masterrank)

   NFE_USE_AFAILED

   use nfe_umbrella, only : umbrella_send_coeffs

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   integer, intent(in) :: dst_masterrank

   integer :: error

#  include "nfe-mpi.h"

   nfe_assert(sanderrank.eq.0)
   nfe_assert(self%initialized)

   if (self%imode.eq.MODE_FLOODING) then
      call umbrella_send_coeffs(self%umbrella, dst_masterrank, commmaster)
   else
      ! for synchronization
      call mpi_send(masterrank, 1, MPI_INTEGER, &
                    dst_masterrank, 8, commmaster, error)
      nfe_assert(error.eq.0)
   end if

end subroutine ctxt_send

!-------------------------------------------------------------------------------

subroutine ctxt_recv(self, src_masterrank)

   NFE_USE_AFAILED

   use nfe_umbrella, only : umbrella_recv_coeffs

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   integer, intent(in) :: src_masterrank

#  include "nfe-mpi.h"

   integer :: error, itmp

   nfe_assert(sanderrank.eq.0)
   nfe_assert(self%initialized)

   if (self%imode.eq.MODE_FLOODING) then
      call umbrella_recv_coeffs(self%umbrella, src_masterrank, commmaster)
   else
      ! for synchronization
      call mpi_recv(itmp, 1, MPI_INTEGER, src_masterrank, 8, &
                    commmaster, MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
      nfe_assert(itmp.eq.src_masterrank)
   end if

end subroutine ctxt_recv

!-------------------------------------------------------------------------------

subroutine ctxt_close_units(self)

   NFE_USE_AFAILED

   use nfe_constants
   use nfe_sander_proxy

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

#  include "nfe-mpi.h"

   nfe_assert(self%initialized)
   nfe_assert(sanderrank.eq.0)

   if (self%imode.ne.MODE_NONE) &
      close (MONITOR_UNIT)

   if (self%mdout.ne.'stdout') &
      close (OUT_UNIT)

   call close_dump_files()

end subroutine ctxt_close_units

!-------------------------------------------------------------------------------

subroutine ctxt_open_units(self)

   NFE_USE_AFAILED

   use nfe_constants
   use nfe_sander_proxy
   use file_io_dat

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

#  include "nfe-mpi.h"

   integer :: error

   nfe_assert(self%initialized)
   nfe_assert(sanderrank.eq.0)

   if (self%imode.ne.MODE_NONE) then
      open (unit = MONITOR_UNIT, file = self%monitor_file, iostat = error, &
        form = 'FORMATTED', action = 'WRITE', status = 'OLD', position = 'APPEND')
      if (error.ne.0) then
         write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
            NFE_ERROR, 'failed to open ''', trim(self%monitor_file), &
            ''' for appending'
         call terminate()
      end if
   end if ! self%imode.ne.MODE_NONE

   mdout = self%mdout
   mdinfo = self%mdinfo
   restrt = self%restrt
   mdvel = self%mdvel
   mden = self%mden
   mdcrd = self%mdcrd

   ioutfm = self%ioutfm

   ntpr = self%ntpr
   ntwr = self%ntwr
   ntwx = self%ntwx

   if (self%mdout.ne.'stdout') &
      call amopen(OUT_UNIT, mdout, 'O', 'F', 'A')

   call amopen(MDINFO_UNIT, mdinfo, 'U', 'F', 'W')

   facc = 'A'
   owrite = 'U'
   call open_dump_files()
   facc = 'W'

end subroutine ctxt_open_units

!-------------------------------------------------------------------------------

#endif /* NFE_ENABLE_BBMD */

end module nfe_bbmd_ctxt
