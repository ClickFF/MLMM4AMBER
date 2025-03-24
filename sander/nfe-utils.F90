#include "nfe-utils.h"
#include "nfe-config.h"

module nfe_utils

implicit none

private

public :: fatal
public :: out_of_memory

public :: close_UNIT

#ifdef MPI
public :: cpus_enter
public :: cpus_leave
#endif /* MPI */

#ifndef NFE_DISABLE_ASSERT
public :: afailed
#endif /* NFE_DISABLE_ASSERT */

interface swap
   module procedure swap_i, swap_r4, swap_r8, swap_cSL
end interface swap

private :: swap_i, swap_r4, swap_r8, swap_cSL

public :: swap

! pfmt - print format is an interface to functions taking numbers
! and emitting character strings that are format edit descriptors,
! eg, i5 or f7.3.

interface pfmt
   module procedure i0_format_1, f0_format_1
end interface

private :: i0_format_1, f0_format_1, n_digits_i, n_digits_r

public :: pfmt ! this is needed because "I0" and "F0.X" editing is not portable

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

subroutine fatal(message)

   use nfe_constants, only : ERR_UNIT
   use nfe_sander_proxy, only : terminate, is_master

   implicit none

   character(len = *), intent (in) :: message

   if (is_master()) &
      write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, message

   call terminate()

end subroutine fatal

!-----------------------------------------------------------------------------

subroutine out_of_memory(filename, lineno)

   use nfe_constants, only : ERR_UNIT
   use nfe_sander_proxy, only : terminate

   implicit none

   character(len = *), intent(in) :: filename
   integer,            intent(in) :: lineno

   write (unit = ERR_UNIT, fmt = '(a,a,a,a,'//pfmt(lineno)//')') &
      NFE_ERROR, 'memory allocation failed at ', filename, ':', lineno

   call terminate()

end subroutine out_of_memory

!-----------------------------------------------------------------------------

subroutine close_UNIT(u)

   implicit none

   integer, intent(in) :: u

   logical :: o

   inquire (unit = u, opened = o)
   if (o) close (unit = u)

end subroutine close_UNIT

!-----------------------------------------------------------------------------

#ifndef NFE_DISABLE_ASSERT
subroutine afailed(filename, lineno)

   use nfe_constants, only : ERR_UNIT
   use nfe_sander_proxy, only : terminate, flush_UNIT

   implicit none

   character(len = *), intent(in) :: filename
   integer,            intent(in) :: lineno

   write(unit = ERR_UNIT, fmt = '(/a,a,a,'//pfmt(lineno)//',a/)') &
         NFE_ERROR, filename, ':', lineno, ': nfe_assert() failed'
   call flush_UNIT(ERR_UNIT)
   call terminate()

end subroutine afailed
#endif /* NFE_DISABLE_ASSERT */

!-----------------------------------------------------------------------------

#ifdef MPI

!
! stolen from mpb-1.4.2
!

subroutine cpus_enter(comm, tag)

   implicit none

   integer, intent(in) :: comm, tag

#  include "nfe-mpi.h"

   integer :: commrank, commsize
   integer :: recv_tag, recv_status(MPI_STATUS_SIZE), error

   call mpi_comm_rank(comm, commrank, error)
   nfe_assert(error.eq.0)

   call mpi_comm_size(comm, commsize, error)
   nfe_assert(error.eq.0)

   nfe_assert(commrank.ge.0)
   nfe_assert(commsize.gt.0)
   nfe_assert(commrank.lt.commsize)

   if (commrank.gt.0) then
      recv_tag = tag - 1
      call mpi_recv(recv_tag, 1, MPI_INTEGER, &
         commrank - 1, tag, comm, recv_status, error)
      nfe_assert(error.eq.0)
      nfe_assert(recv_tag.eq.tag)
   end if

end subroutine cpus_enter

!-----------------------------------------------------------------------------

subroutine cpus_leave(comm, tag)

   implicit none

   integer, intent(in) :: comm, tag

#include "nfe-mpi.h"

   integer :: error, commsize, commrank

   call mpi_comm_rank(comm, commrank, error)
   nfe_assert(error.eq.0)

   call mpi_comm_size(comm, commsize, error)
   nfe_assert(error.eq.0)

   nfe_assert(commrank.ge.0)
   nfe_assert(commsize.gt.0)
   nfe_assert(commrank.lt.commsize)

   if (commrank.ne.(commsize - 1)) then
      call mpi_send(tag, 1, MPI_INTEGER, &
         commrank + 1, tag, comm, error)
      nfe_assert(error.eq.0)
   end if

end subroutine cpus_leave
#endif /* MPI */

!-----------------------------------------------------------------------------

subroutine swap_i(a, b)

   implicit none

   integer, intent(inout) :: a, b

   integer :: tmp

   tmp = a
   a = b
   b = tmp

end subroutine swap_i

!-----------------------------------------------------------------------------

subroutine swap_r4(a, b)

   implicit none

   real(4), intent(inout) :: a, b

   real(4) :: tmp

   tmp = a
   a = b
   b = tmp

end subroutine swap_r4

!-----------------------------------------------------------------------------

subroutine swap_r8(a, b)

   implicit none

   real(8), intent(inout) :: a, b

   real(8) :: tmp

   tmp = a
   a = b
   b = tmp

end subroutine swap_r8

!-----------------------------------------------------------------------------

subroutine swap_cSL(a, b)

   use nfe_constants, only : STRING_LENGTH

   implicit none

   character(STRING_LENGTH), intent(inout) :: a, b

   character(STRING_LENGTH) :: tmp

   tmp = a
   a = b
   b = tmp

end subroutine swap_cSL

!-----------------------------------------------------------------------------
! emit an integer edit descriptor of width the number of digits in x, eg, i4

function i0_format_1(x) result(f)

   implicit none

   character(8) :: f
   integer, intent(in) :: x

   integer :: n

   n = n_digits_i(x)

   if (n.le.9) then
      write (unit = f, fmt = '(a,i1)') 'i', n
   else
      write (unit = f, fmt = '(a,i2)') 'i', n
   end if

end function i0_format_1

!-----------------------------------------------------------------------------
! emit a fixed-real edit descriptor of width the number of digits in x
! and with y digits after the decimal point, eg, f8.4

function f0_format_1(x, y) result(f)

   implicit none

   character(8) :: f
   real(8), intent(in) :: x
   integer, intent(in) :: y

   integer :: n

   n = n_digits_r(x)

   if ((n + y + 1).le.9) then
      write (unit = f, fmt = '(a,i1,a,i1)') 'f', (n + y + 1), '.', y
   else
      write (unit = f, fmt = '(a,i2,a,i1)') 'f', (n + y + 1), '.', y
   end if

end function f0_format_1

!-----------------------------------------------------------------------------

pure function n_digits_i(i) result(n)

   implicit none

   integer :: n

   integer, intent(in)  :: i

   n = n_digits_r(dble(i))

end function n_digits_i

!-----------------------------------------------------------------------------

pure function n_digits_r(r) result(n)

   implicit none

   integer :: n

   real(8), intent(in)  :: r

   if (r.gt.1.0D0) then
      n = 1 + int(floor(log10(r)))
   else if (r.ge.0.0D0) then
      n = 1
   else if (r.ge.-1.0D0) then
      n = 2
   else
      n = 2 + int(floor(log10(-r)))
   end if

end function n_digits_r

!-----------------------------------------------------------------------------

end module nfe_utils
