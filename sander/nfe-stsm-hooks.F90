!<compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

module nfe_stsm_hooks

use nfe_constants, only : SL => STRING_LENGTH, STSM_OUTPUT_UNIT, STSM_CV_UNIT
use nfe_colvar_type, only : colvar_t

implicit none

private

public :: on_multisander_exit

public :: on_sander_init
public :: on_sander_exit

public :: on_force

!- - - - - - - - - - - - - - - - P R I V A T E - - - - - - - - - - - - - - - -

integer, private, parameter :: stsm_UNIT = STSM_OUTPUT_UNIT
integer, private, parameter :: CV_UNIT = STSM_CV_UNIT

character(*), private, parameter :: DEFAULT_OUTPUT_FILE = 'nfe-stsm'
character(*), private, parameter :: DEFAULT_CV_FILE = 'nfe-stsm-cv'

integer, private, parameter :: DEFAULT_OUTPUT_FREQ = 50

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer, private, save :: ncolvars = 0 ! .gt.0 means "active"
type(colvar_t), private, allocatable, save :: cv(:)

NFE_REAL, private, pointer, save :: anchor(:) => null() ! master
NFE_REAL, private, pointer, save :: a_position(:) => null() ! master
NFE_REAL, private, pointer, save :: a_strength(:) => null() ! master

NFE_REAL, private, allocatable, save :: cv_inst(:)
NFE_REAL, private, allocatable, save :: cv_copy(:)

NFE_REAL, private, allocatable, save :: f_cv(:)

NFE_REAL, private, allocatable, save :: x_eq(:)

NFE_REAL, private, allocatable, save :: all_cv(:)
NFE_REAL, private, allocatable, save :: all_anchor(:)
integer, private, allocatable, save :: all_image(:)


character(SL), private, save :: output_file
character(SL), private, save :: output_fmt
character(SL), private, save :: output_fmt_centers
character(SL), private, save :: cv_file = DEFAULT_CV_FILE

integer, private, save :: output_freq = DEFAULT_OUTPUT_FREQ
integer, private, save :: mdstep ! = runmd.f::nstep + 1 (not zeroed on exchange)

integer, private, save :: eq_freq
integer, private, save :: re_freq
integer, private, save :: run_freq
integer, private, save :: update_freq
integer, private, save :: num_repeats
integer, private, save :: num_copies
integer, private, save :: num_images

integer, private, save :: equilibration = 0
integer, private, save :: release = 1
integer, private, save :: copies = 1
integer, private, save :: repeats = 1

integer, private, save :: id_image
integer, private, save :: image = 0

logical, private, save :: report_drift
logical, private, save :: report_smooth
logical, private, save :: report_reparam

NFE_REAL, private, save :: smooth

double precision, private, save :: smoothing = 0.0
character(SL), private, save :: report_centers = 'NONE'

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

namelist / stsm /    image, repeats, equilibration, release, &
                     smoothing, report_centers, &
                     output_file, output_freq, cv_file

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

contains

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_multisander_exit()

   use nfe_colvar, only : colvar_cleanup

   implicit none

   integer :: n

#  include "nfe-mpi.h"

   if (ncolvars.gt.0) then
      do n = 1, ncolvars
         call colvar_cleanup(cv(n))
      end do
      deallocate(cv, f_cv, cv_inst)

      NFE_MASTER_ONLY_BEGIN
      nullify(a_position, a_strength)
      deallocate(anchor,all_cv,all_anchor,all_image)
      NFE_MASTER_ONLY_END
   end if

   mdstep = 0
   ncolvars = 0

end subroutine on_multisander_exit

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_sander_init(mdin_unit, amass)

   use nfe_utils
   use nfe_colvar
   use nfe_colvar_type
   use nfe_constants
   use nfe_sander_proxy

   implicit none

   integer, intent(in) :: mdin_unit


   NFE_REAL, intent(in) :: amass(*)

   integer :: n, error, ifind, i
   character(len = 80) :: buf

   integer :: LOG_UNIT

#  ifdef MPI
#  include "nfe-mpi.h"
   nfe_assert(multisander_rem().eq.0)
#  endif /* MPI */

   NFE_MASTER_ONLY_BEGIN

   nfe_assert(ncolvars.eq.0)

#  ifdef MPI
   if (multisander_numgroup().gt.1) then
      n = masterrank + 1
      write (unit = output_file, fmt = '(a,a,i3.3,a)') &
           DEFAULT_OUTPUT_FILE, '.', n, '.dat'
   else
#  endif /* MPI */
      output_file = DEFAULT_OUTPUT_FILE//'.dat'
#  ifdef MPI
   end if
#  endif /* MPI */

   rewind(mdin_unit)
   call nmlsrc('stsm', mdin_unit, ifind)
   ! no stsm section
   if (ifind.eq.0) then
      goto 1
   end if

   ncolvars = 0

   rewind(mdin_unit)
   read(mdin_unit,nml=stsm,err=666)

   call amopen(CV_UNIT, cv_file, 'O', 'F', 'R')

   do
     call nmlsrc('colvar', CV_UNIT, ifind)
     if (ifind.eq.0) exit
     read(CV_UNIT,'(a80)') buf
     ncolvars = ncolvars + 1
   end do

   if (ncolvars.eq.0) &
      call fatal('no variable(s) in the CV file')

   allocate(cv(ncolvars), anchor(2*ncolvars), &
      f_cv(ncolvars), cv_inst(ncolvars), stat = error)
   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   a_position => anchor(0*ncolvars + 1:1*ncolvars)
   a_strength => anchor(1*ncolvars + 1:2*ncolvars)

   n = 1

   do while (n.le.ncolvars)
         nfe_assert(n.le.ncolvars)
         call colvar_nlread(CV_UNIT, cv(n))
         a_strength(n) = anchor_strength(1)
         a_position(n) = anchor_position(1)
         if (colvar_has_refcrd(cv(n))) then
            refcrd_file = trim(refcrd_file)
            refcrd_len = len_trim(refcrd_file)
         end if
         if (colvar_is_quaternion(cv(n))) then
          allocate(cv(n)%q_index, stat = error)
           if (error.ne.0) &
            NFE_OUT_OF_MEMORY
            cv(n)%q_index = q_index
         end if 
         if (colvar_has_axis(cv(n))) then
           allocate(cv(n)%axis(3), stat = error)
             if (error.ne.0) &
                NFE_OUT_OF_MEMORY
           i = 1
           do while (i.le.3)
             cv(n)%axis(i) = axis(i)
             i = i + 1
           end do
         end if

         n = n + 1
   end do

   output_freq = min(output_freq, sander_nstlim())
   output_freq = max(1, output_freq)

   eq_freq = equilibration
   re_freq = release
   num_copies = copies
   num_repeats = repeats

   run_freq = re_freq + eq_freq
   run_freq = min(run_freq, sander_nstlim())
   run_freq = max(1, run_freq)

   update_freq = num_repeats * run_freq
  
   id_image = image
   if (id_image.eq.0) &
#  ifdef MPI
      id_image = mod(masterrank,num_copies)+1
#  else
      id_image = 1
#  endif /* MPI */
    id_image = id_image - 1
   
   smooth = smoothing

   if (report_centers=='NONE') then
      report_drift = .false.
      report_smooth = .false.
      report_reparam = .false.
   else if (report_centers=='ALL') then
      report_drift = .true.
      report_smooth = .true.
      report_reparam = .true.
   else if (report_centers=='DRIFT') then
      report_drift = .true.
      report_smooth = .false.
      report_reparam = .false.
   else if (report_centers=='SMOOTHED') then
      report_drift = .false.
      report_smooth = .true.
      report_reparam = .false.
   else if (report_centers=='REPARAMETRIZED') then
      report_drift = .false.
      report_smooth = .false.
      report_reparam = .true.
   else if (report_centers=='NO_DRIFT') then
      report_drift = .false.
      report_smooth = .true.
      report_reparam = .true.
   else if (report_centers=='NO_SMOOTHED') then
      report_drift = .true.
      report_smooth = .false.
      report_reparam = .true.
   else if (report_centers=='NO_REPARAMETRIZED') then
      report_drift = .true.
      report_smooth = .true.
      report_reparam = .false.
   end if

#  ifdef MPI
      if (multisander_numgroup().gt.1) then
        nfe_assert(commmaster.ne.mpi_comm_null)
        ! get all image ids
        allocate(all_image(multisander_numgroup()), stat = error)
        call mpi_allgather(id_image, 1, &
            MPI_INTEGER, all_image, 1, &
            MPI_INTEGER, commmaster,error)
        nfe_assert(error.eq.0)
        call mpi_barrier(commmaster,error)
        nfe_assert(error.eq.0)
        if (masterrank.eq.0) then
            num_images = maxval(all_image)+1
        end if
        call mpi_bcast(num_images, 1, MPI_INTEGER, 0, commmaster, error)
        nfe_assert(error.eq.0)
        allocate(all_anchor(ncolvars*num_images), stat = error)
        allocate(all_cv(ncolvars*multisander_numgroup()), stat = error)
      end if
#  endif /* MPI */

1  continue

   NFE_MASTER_ONLY_END

#  ifdef MPI
   call mpi_bcast(ncolvars, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)
 
   call mpi_bcast(refcrd_len, 1, MPI_INTEGER, 0, commsander, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(refcrd_file, refcrd_len, MPI_CHARACTER, 0, commsander, error)
   nfe_assert(error.eq.0)
#  endif /* MPI */

   if (ncolvars.eq.0) &
      return

#  ifdef MPI
   if (sanderrank.ne.0) then
      allocate(cv(ncolvars), f_cv(ncolvars), cv_inst(ncolvars), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY
   end if
#  endif /* MPI */

   do n = 1, ncolvars
      call colvar_bootstrap(cv(n), n, amass)
   end do

   mdstep = 0

   NFE_MASTER_ONLY_BEGIN

   open (unit = stsm_UNIT, file = output_file, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

   if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'failed to open ''', trim(output_file), ''' for writing'
      call terminate()
   end if

   write (unit = stsm_UNIT, fmt = '(a,50(''=''))') '# = NFE%STSM Initial Pathway'
   do n = 1, ncolvars
      write (unit = stsm_UNIT, &
         fmt = '(a,'//pfmt(n)//',a,'//pfmt(a_position(n), 6)//',a,'//pfmt(a_strength(n), 6)//',a)') &
         '#   << anchor(', n, ') : position = ', a_position(n), &
         ', strength = ', a_strength(n), ' >>'
   end do
   write (unit = stsm_UNIT, &
      fmt = '(a,77(''-''),/a,'//pfmt(ncolvars)//',a,/a,77(''=''))') '# ', &
      '# MD time (ps), CV(1:', ncolvars, ')', '# '

   call flush_UNIT(stsm_UNIT)

   write (unit = output_fmt, fmt = '(a,'//pfmt(ncolvars)//',a)') &
      '(f12.4,', ncolvars, '(1x,f16.8))'

   write (unit = output_fmt_centers, fmt = '(a,'//pfmt(ncolvars)//',a)') &
      '(i12,', ncolvars, '(1x,f16.8))'

   ! print summary & we'r done

   LOG_UNIT = OUT_UNIT ! write to MDOUT

   write (unit = LOG_UNIT, fmt = '(a,a)') NFE_INFO, &
      '~~ ~~ ~~ ~~ ~~ STRING METHOD WITH SWARMS OF TRAJECTOTIES ~~ ~~ ~~ ~~ ~~'

   write (unit = LOG_UNIT, fmt = '(a,/a,a,a)') NFE_INFO, NFE_INFO, &
      'output_file = ', trim(output_file)
   write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(output_freq)//',a,'//pfmt(output_freq*sander_timestep(), 4)//',a)') &
        NFE_INFO, 'output_freq = ', output_freq, ' (', &
        output_freq*sander_timestep(), ' ps)'

   write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(eq_freq)//',a,'//pfmt(eq_freq*sander_timestep(), 4)//',a)') &
        NFE_INFO, 'equilibration per iteration = ', eq_freq, ' (', &
        eq_freq*sander_timestep(), ' ps)'

   write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(re_freq)//',a,'//pfmt((re_freq)&
        *sander_timestep(), 4)//',a)') &
        NFE_INFO, 'release per iteration = ', re_freq, ' (', &
        re_freq*sander_timestep(), ' ps)'

   write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(num_repeats)//')') &
        NFE_INFO, 'number of repeats per copy = ', num_repeats

   write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(multisander_numgroup()/num_images)//')') &
        NFE_INFO, 'number of copies per image = ', multisander_numgroup()/num_images

   write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(id_image+1)//',a,'//pfmt(num_images)//')') &
        NFE_INFO, 'image id = ', id_image+1,'/',num_images

   write (unit = LOG_UNIT, &
      fmt = '(a,a,'//pfmt(smooth, 4)//')') &
        NFE_INFO, 'smoothing strength = ', smooth

   if (report_drift) &
   write (unit = LOG_UNIT, &
      fmt = '(a,a)') &
        NFE_INFO, 'drifted centers will be reported.'

   if (report_smooth) &
   write (unit = LOG_UNIT, &
      fmt = '(a,a)') &
        NFE_INFO, 'smoothed centers will be reported.'

   if (report_reparam) &
   write (unit = LOG_UNIT, &
      fmt = '(a,a)') &
        NFE_INFO, 'reparametrized centers will be reported.'

   write (unit = LOG_UNIT, fmt = '(a)') NFE_INFO
   do n = 1, ncolvars
      write (unit = LOG_UNIT, &
         fmt = '(a,a,'//pfmt(n)//',a,'//pfmt(a_position(n), 6)//',a,'//pfmt(a_strength(n), 6)//',a)') &
         NFE_INFO, 'CV #', n, ' << anchor : position = ', a_position(n), &
         ', strength = ', a_strength(n), ' >>'
      call colvar_print(cv(n), LOG_UNIT)
      write (unit = LOG_UNIT, fmt = '(a)') NFE_INFO
      if (colvar_is_quaternion(cv(n))) then
          write (unit = OUT_UNIT, fmt = '(a,a,I3)') NFE_INFO, &
          ' <> q_index = ', cv(n)%q_index
      end if
      if (colvar_has_axis(cv(n))) then
           write (unit = OUT_UNIT, fmt = '(a,a,f8.4,a,f8.4,a,f8.4,a)') NFE_INFO, &
           ' <> axis = [',cv(n)%axis(1),', ', cv(n)%axis(2),', ', cv(n)%axis(3),']'
      end if

   end do

   write (unit = LOG_UNIT, fmt = '(a,a/)') NFE_INFO, &
      '~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~'

   NFE_MASTER_ONLY_END

   return
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &stsm namelist!'
    call terminate()

end subroutine on_sander_init

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_sander_exit()

   use nfe_utils, only : close_UNIT

   implicit none

#  include "nfe-mpi.h"

   NFE_MASTER_ONLY_BEGIN
   call close_UNIT(stsm_UNIT)
   NFE_MASTER_ONLY_END

end subroutine on_sander_exit

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_force(x, f, pot)

   NFE_USE_AFAILED

   use nfe_utils
   use nfe_colvar
   use nfe_constants
   use nfe_sander_proxy
   use nfe_colvar_type

   implicit none

   NFE_REAL, intent(in) :: x(*)

   NFE_REAL, intent(inout) :: f(*)
   NFE_REAL, intent(inout) :: pot

   integer :: n, q
   integer, DIMENSION(4) :: cv_q = (/COLVAR_QUATERNION0, COLVAR_QUATERNION1, &
                                    COLVAR_QUATERNION2, COLVAR_QUATERNION3/)
   NFE_REAL :: norm4(100), cv_N(4, 100)

#  ifdef MPI
   integer :: m
   integer :: mm
   integer :: copy_index

   NFE_REAL :: curlen
   NFE_REAL :: ratio
   NFE_REAL, allocatable :: totlen(:)
   integer, allocatable :: totcopies(:)
#  endif /* MPI */

#  ifdef MPI
#     include "nfe-mpi.h"
   integer :: error
#  endif /* MPI */

   if (ncolvars.eq.0) &
      return

   do n = 1, ncolvars
      cv_inst(n) = colvar_value(cv(n), x)
   end do

   NFE_MASTER_ONLY_BEGIN
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n))) then
        do q = 1, 4
         if (cv(n)%type == cv_q(q)) then
            cv_N(q, cv(n)%q_index) = cv_inst(n)
         end if
        end do
      end if
   end do
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n))) then
       norm4(cv(n)%q_index) = sqrt(cv_N(1,cv(n)%q_index)**2 &
                                 + cv_N(2,cv(n)%q_index)**2 &
                                 + cv_N(3,cv(n)%q_index)**2 &
                                 + cv_N(4,cv(n)%q_index)**2)
      end if
   end do
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n))) then
         cv_inst(n) = cv_inst(n) / norm4(cv(n)%q_index)
      else
         cv_inst(n) = cv_inst(n)
      end if
   end do

   if (mod(mdstep,run_freq).eq.0) then
      if (mod(mdstep/run_freq,num_repeats).eq.0) then
         write (unit = OUT_UNIT, fmt = '(/a,a)') &
            NFE_INFO,'#   new restraint:'
      else 
         write (unit = OUT_UNIT, fmt = '(/a,a)') &
            NFE_INFO,'#   restoring restraint:'
      end if
      do n = 1, ncolvars
         write (unit = OUT_UNIT, fmt = &
            '(a,a,'//pfmt(n)//',a,'//pfmt(a_position(n), 6)//',a)') &
            NFE_INFO,'#   << colvar(', n,') = ',a_position(n),' >>'
      end do
      write (unit = OUT_UNIT, fmt = '(a,a/)') &
         NFE_INFO,'#   equilibration begins...'
   end if
   
   if (mod(mdstep,run_freq).lt.eq_freq) then
      do n = 1, ncolvars
         f_cv(n) = &
            - a_strength(n)*colvar_difference(cv(n), cv_inst(n), a_position(n))
         pot = pot + a_strength(n)*colvar_difference(cv(n), cv_inst(n), a_position(n))**2/2
      end do
   else
      f_cv(1:ncolvars) = 0.0
   end if
   NFE_MASTER_ONLY_END

#  ifdef MPI
   call mpi_bcast(f_cv, ncolvars, &
      MPI_DOUBLE_PRECISION, 0, commsander, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(pot, 1, &
      MPI_DOUBLE_PRECISION, 0, commsander, error)
   nfe_assert(error.eq.0)
#  endif /* MPI */

   ! FIXME: virial
   do n = 1, ncolvars
      call colvar_force(cv(n), x, f_cv(n), f)
   end do

   NFE_MASTER_ONLY_BEGIN
   if (nfe_real_mdstep) then

      if (mod(mdstep, output_freq).eq.0) then
         write (unit = stsm_UNIT, fmt = output_fmt) &
            sander_mdtime(), cv_inst(1:ncolvars)
         call flush_UNIT(stsm_UNIT)
      end if

#ifdef MPI
      NFE_MASTER_ONLY_BEGIN
      if (mod(mdstep,run_freq).eq.run_freq-1) then
         copy_index = (mdstep+1)/run_freq
         copy_index = mod(copy_index-1,num_repeats)+1
         if (copy_index.gt.1) then
            do m = 1, ncolvars
               cv_copy(m) = colvar_interpolate( cv(m)     &
                  , ONE/copy_index, cv_inst(m)            &
                  , ONE-(ONE/copy_index), cv_copy(m) )
            end do
         else
            allocate(cv_copy(ncolvars))
            cv_copy(:) = cv_inst(:)
         end if
         do m = 1, ncolvars
           write (unit = OUT_UNIT, fmt = &
             '(a,a,'//pfmt(n)//',a,'//pfmt(cv_copy(m), 6)//',a,'&
              //pfmt(cv_copy(m), 6)//',a)') &
             NFE_INFO,'#   << colvar(', m,') = ', &
             cv_copy(m),' ',cv_inst(m),' >>'
         end do
      end if
      if (mod((mdstep+1),update_freq).eq.0) then
         if (multisander_numgroup().gt.1) then
            nfe_assert(commmaster.ne.mpi_comm_null)
            ! get all cv_inst
            nfe_assert(allocated(all_cv))
            call mpi_allgather(cv_copy, ncolvars, &
               MPI_DOUBLE_PRECISION, all_cv, ncolvars, &
               MPI_DOUBLE_PRECISION, commmaster,error)
            nfe_assert(error.eq.0)
            call mpi_barrier(commmaster,error)
!	    <<recalculate the centers>>
            if (masterrank.eq.0) then
!	       Drift
               all_anchor(1:num_images*ncolvars) = 0.0
               allocate(totcopies(num_images))
               totcopies(1:num_images) = 0
               do n = 0, multisander_numgroup()-1
                  totcopies(all_image(n+1)+1) = &
                     totcopies(all_image(n+1)+1) + 1
                  do m = 1, ncolvars
                     all_anchor(all_image(n+1)*ncolvars+m) &
                     = colvar_interpolate(cv(m)          &
                     ,(totcopies(all_image(n+1)+1)-1)*ONE/ &
                       totcopies(all_image(n+1)+1)         &
                     ,all_anchor(all_image(n+1)*ncolvars+m)&
                     ,ONE/totcopies(all_image(n+1)+1)      &
                       , all_cv(n*ncolvars+m))
                  end do
               end do
               deallocate(totcopies)
               if (report_drift) then
                  do n = 0, num_images-1
                     write (unit = OUT_UNIT, fmt = &
                       '(a,a,'//pfmt(n+1)//',a)',advance="no") &
                       NFE_INFO,'#  drifted center of image ',n+1,' : '
                     write (unit = OUT_UNIT, fmt = output_fmt_centers) &
                       (mdstep+1)/run_freq, all_anchor(n*ncolvars+1:(n+1)*ncolvars)
                     call flush_UNIT(OUT_UNIT)
                  end do
               end if
!              Smoothing
               all_cv(1:ncolvars)=all_anchor(1:ncolvars)
               if (num_images.gt.2) then
                  do n = 1, num_images-2
                     do m = 1, ncolvars
                        all_cv(n*ncolvars+m) &
                        = colvar_interpolate(cv(m) &
                        , 1.0-smooth &
                        , all_anchor(n*ncolvars+m) &
                        , 0.5 * smooth &
                        , all_anchor((n-1)*ncolvars+m) &
                          +all_anchor((n+1)*ncolvars+m))
                     end do
                  end do
               end if
               all_cv((num_images-1)*ncolvars:num_images*ncolvars) &
               =all_anchor((num_images-1)*ncolvars:num_images*ncolvars)
               if (report_smooth) then
                  do n = 0, num_images-1
                     write (unit = OUT_UNIT, fmt = &
                       '(a,a,'//pfmt(n+1)//',a)',advance="no") &
                       NFE_INFO,'#  smoothed center of image ',n+1,' : '
                     write (unit = OUT_UNIT, fmt = output_fmt_centers) &
                       (mdstep+1)/run_freq, all_cv(n*ncolvars+1:(n+1)*ncolvars)
                     call flush_UNIT(OUT_UNIT)
                  end do
               end if
!	       Reparametrizing
               if (num_images.gt.2) then
                  allocate(totlen(num_images))
                  totlen(1) = 0.0
                  do n = 1, num_images-1
                     totlen(n+1) = 0.0
                     do m = 1, ncolvars
                        totlen(n+1) = totlen(n+1) &
                           + colvar_difference(cv(m) &
                           , all_cv(n*ncolvars+m) &
                           , all_cv((n-1)*ncolvars+m))**2
                     end do
                     totlen(n+1) = totlen(n) + sqrt(totlen(n+1))
                  end do
                  n = 0
                  m = 0
                  do while ( n < num_images-2 )
                     n = n + 1
                     curlen = n*totlen(num_images)/(num_images-1)
                     do while (m<num_images-1.and.curlen>totlen(m+1))
                        m = m + 1
                     end do
                     ratio = (curlen-totlen(m)) &
                            / (totlen(m+1)-totlen(m))
                     do mm = 1, ncolvars
                        all_anchor(n*ncolvars+mm) &
                        = colvar_interpolate(cv(mm) &
                          ,ONE-ratio,all_cv((m-1)*ncolvars+mm) &
                          ,ratio,all_cv(m*ncolvars+mm))
                     end do
                  end do
                  deallocate(totlen)
               end if
               if (report_reparam) then
                  do n = 0, num_images-1
                     write (unit = OUT_UNIT, fmt = &
                       '(a,a,'//pfmt(n+1)//',a)',advance="no") &
                       NFE_INFO,'#  reparametrized center of image ',n+1,' : '
                     write (unit = OUT_UNIT, fmt = output_fmt_centers) &
                       (mdstep+1)/run_freq, all_anchor(n*ncolvars+1:(n+1)*ncolvars)
                     call flush_UNIT(OUT_UNIT)
                  end do
               end if
            end if
!	    <<broadcast the new centers>>
            nfe_assert(allocated(all_anchor))
            call mpi_bcast(all_anchor, &
               ncolvars*num_images, &
               MPI_DOUBLE_PRECISION, 0, commmaster,error)
            nfe_assert(error.eq.0)
            a_position(1:ncolvars) = &
              all_anchor(all_image(masterrank+1)*ncolvars+1&
              : (all_image(masterrank+1)+1)*ncolvars)
         else
            a_position(:) = cv_copy(:)
         end if
         nfe_assert(allocated(cv_copy))
         deallocate(cv_copy)
      end if
      NFE_MASTER_ONLY_END
#endif /* MPI */

      mdstep = mdstep + 1

   end if ! nfe_real_mdstep
  
   NFE_MASTER_ONLY_END

end subroutine on_force

end module nfe_stsm_hooks
