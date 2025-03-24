!F90 ISO_C_BINDING wrapper for socket communication.

!Copyright (C) 2013, Michele Ceriotti

!Permission is hereby granted, free of charge, to any person obtaining
!a copy of this software and associated documentation files (the
!"Software"), to deal in the Software without restriction, including
!without limitation the rights to use, copy, modify, merge, publish,
!distribute, sublicense, and/or sell copies of the Software, and to
!permit persons to whom the Software is furnished to do so, subject to
!the following conditions:

!The above copyright notice and this permission notice shall be included
!in all copies or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!Contains both the functions that transmit data to the socket and read the data
!back out again once finished, and the function which opens the socket initially.

!Functions:
!   open_socket: Opens a socket with the required host server, socket type and
!      port number.
!   write_buffer: Writes a string to the socket.
!   read_buffer: Reads data from the socket.

module F90SOCKETS
  use ISO_C_BINDING

  implicit none

  integer,parameter  :: AF_INET = 0, AF_UNIX = 1
  integer,parameter  :: AI_PASSIVE = 1
  integer,parameter  :: SOCK_STREAM = 1
  integer,parameter  :: UNIX_PATH_MAX = 108

  type, bind(C) :: sockaddr_un
     integer(KIND=C_SHORT) sun_family
     character(LEN=1, KIND=C_CHAR) sun_path(UNIX_PATH_MAX)
  end type sockaddr_un

  type, bind(C) :: addrinfo
     integer(KIND=C_INT) :: ai_flags, ai_family, ai_socktype,&
          & ai_protocol
     integer(KIND=C_SIZE_T) :: ai_addrlen
     type(C_PTR) :: ai_addr, ai_canonname, ai_next
  end type addrinfo

  interface
     function getaddrinfo(node, service, hints, res) bind(C)
       use ISO_C_BINDING
       type(C_PTR), intent(IN), value :: node, service, hints
       type(C_PTR) :: RES
       integer(KIND=C_INT) :: getaddrinfo
     end function getaddrinfo
  end interface

  interface
     subroutine freeaddrinfo( res) bind(C)
       use ISO_C_BINDING
       type(C_PTR), intent(IN), value :: RES
     end subroutine freeaddrinfo
  end interface

  interface
     function socket_make(af, type, protocol) bind(C, name="socket")
       use ISO_C_BINDING
       integer(KIND=C_INT), intent(IN), value   :: af, type,&
            & protocol
       integer(KIND=C_INT)  :: socket_make
     end function socket_make
  end interface

  interface
     function socket_connect(sockfd, addr, addrlen) bind(C, name&
          &="connect")
       use ISO_C_BINDING
       integer(KIND=C_INT), value  :: sockfd
       type(C_PTR), intent(IN), value :: addr
       integer(KIND=C_SIZE_T), intent(IN), value  :: addrlen
       integer(KIND=C_INT)  :: socket_connect
     end function socket_connect
  end interface

  interface
     function socket_write(sid, data, dcount) bind(C, name="write")
       use ISO_C_BINDING
       integer(KIND=C_INT), intent(IN), value     :: sid
       type(C_PTR), intent(IN), value     :: data
       integer(KIND=C_SIZE_T), intent(IN), value  :: dcount
       integer(KIND=C_SIZE_T)     :: socket_write
     end function socket_write
  end interface

  interface 
     function socket_read(sid, data, dcount) bind(C, name="read")
       use ISO_C_BINDING
       integer(KIND=C_INT), intent(IN), value     :: sid
       type(C_PTR), intent(IN), value     :: data
       integer(KIND=C_SIZE_T), intent(IN), value  :: dcount
       integer(KIND=C_SIZE_T)     :: socket_read
     end function socket_read
  end interface

  interface
     subroutine memcpy(dout, din, dcount) bind(C, name="memcpy")
       use ISO_C_BINDING
       type(C_PTR), intent(IN), value     :: din, dout
       integer(KIND=C_SIZE_T), intent(IN), value  :: dcount
     end subroutine memcpy
  end interface

  interface
     type(C_PTR) function memset(s, c, n) bind(C, name="memset")
       use ISO_C_BINDING
       type(C_PTR), intent(IN), value     :: s
       integer(KIND=C_INT), intent(IN), value :: c
       integer(KIND=C_SIZE_T), intent(IN), value  :: n
     end function memset
  end interface


  ! writebuffer interfaces
  interface writebuffer
     module procedure writebuffer_s, writebuffer_d, writebuffer_dv,&
          & writebuffer_i

  end interface writebuffer

  ! readbuffer interfaces
  interface readbuffer
     module procedure readbuffer_s, readbuffer_dv, readbuffer_d,&
          & readbuffer_i

  end interface readbuffer
contains

  subroutine fstr2cstr(fstr, cstr, plen)
    implicit none
    character(LEN=*), intent(IN) :: fstr
    character(LEN=1,KIND=C_CHAR), intent(OUT) :: cstr(:)
    integer, intent(IN), optional :: plen

    integer i,n
    if (present(plen)) then
       n = plen
       do i=1,n
          cstr(i) = fstr(i:i)
       enddo
    else
       n = LEN_trim(fstr)
       do i=1,n
          cstr(i) = fstr(i:i)
       enddo
       cstr(n+1) = C_NULL_CHAR
    end if
  end subroutine fstr2cstr

  subroutine open_socket(psockfd, inet, port, host)      
    implicit none
    integer, intent(IN) :: inet, port
    integer, intent(OUT) :: psockfd
    character(LEN=1024), intent(IN) :: host

    integer :: ai_err
    character(LEN=256) :: service
    character(LEN=1,KIND=C_CHAR), target :: cservice(256)
    character(LEN=1,KIND=C_CHAR), target :: chost(1024)

    type(addrinfo), target :: hints
    type(addrinfo), target :: res
    type(sockaddr_un), target :: addrun
    type(C_PTR) :: ptr



    call fstr2cstr(host, chost)
    if (INET>0) then
       ! creates an internet socket

       ! fetches information on the host

       ptr = memset(c_loc(hints), 0, SIZEOF(hints))
       hints%ai_socktype = SOCK_STREAM
       hints%ai_family = AF_INET
       hints%ai_flags = AI_PASSIVE

       write(service,'(I10)') port
       service=adjustl(trim(service))
       call fstr2cstr(service,cservice)

       ai_err = getaddrinfo(c_loc(chost(1)), c_loc(cservice(1)),&
            & c_loc(hints), ptr)
       if (ai_err < 0) then
          write(6,*) "Error fetching host data. Wrong host name?"
          stop " ENDED "
       endif

       call memcpy(c_loc(res), ptr, sizeof(res))
       !WRITE(6,*) "pointer", res
       psockfd = socket_make(res%ai_family, res%ai_socktype, res&
            &%ai_protocol)
       if (psockfd < 0)  then 
          write(6,*) "Error opening socket"
          stop " ENDED "
       endif

       ai_err = socket_connect(psockfd, res%ai_addr, res%ai_addrlen)
       if (ai_err < 0) then
          write(6,*) "Error opening INET socket: wrong port or&
               & server unreachable"
          stop " ENDED "
       endif

       call freeaddrinfo(ptr)
    else
       ! creates an unix socket
       ptr = memset(c_loc(addrun), 0, SIZEOF(addrun))

       addrun%sun_family = AF_UNIX
       call fstr2cstr("/tmp/ipi_"//host, addrun%sun_path) 

       psockfd = socket_make(AF_UNIX, SOCK_STREAM, 0)

       ai_err = socket_connect(psockfd, c_loc(addrun),&
            & sizeof(addrun))
       if (ai_err < 0) then
          write(6,*) "Could not open UNIX socket. Non-existing&
               & path?"
          stop " ENDED "
       endif
    end if
  end subroutine open_socket


  ! Set of wrappers to socket_write. Write data to a socket.
  ! Args:
  !  psockfd: The id of the socket that will be written to.
  !   data: The data to be written to the socket (different kinds)
  !   plen: The length of the data (not for single element writes)

  subroutine writebuffer_s(psockfd, fstring, plen)
    implicit none
    integer, intent(IN) :: psockfd
    integer, intent(IN) ::  plen
    character(LEN=*), target, intent(IN)  :: fstring
    integer(8) :: nwrite, nlen
    character(LEN=1,KIND=C_CHAR), target :: cstring(plen)
    integer i,n

    n = plen
    do i = 1,n
       cstring(i) = fstring(i:i)
    enddo
    nlen = plen
    nwrite = socket_write(psockfd, c_loc(cstring(1)), nlen)
    if (nwrite/=nlen) then
       write(6,*) "Error in writing to socket buffer"
       stop " ENDED "
    endif
  end subroutine writebuffer_s

  subroutine writebuffer_d(psockfd, fdata)
    implicit none
    integer, intent(IN) :: psockfd
    double precision, target, intent(IN)  :: fdata
    integer(8) :: nwrite, nlen

    nlen = 8
    nwrite = socket_write(psockfd, c_loc(fdata), nlen)
    if (nwrite/=nlen) then
       write(6,*) "Error in writing to socket buffer"
       stop " ENDED "
    endif
  end subroutine writebuffer_d

  subroutine writebuffer_i(psockfd, fdata)
    implicit none
    integer, intent(IN) :: psockfd
    integer, target, intent(IN)  :: fdata
    integer(8) :: nwrite, nlen

    nlen = 4
    nwrite = socket_write(psockfd, c_loc(fdata), nlen)
    if (nwrite/=nlen) then
       write(6,*) "Error in writing to socket buffer"
       stop " ENDED "
    endif
  end subroutine writebuffer_i

  subroutine writebuffer_dv(psockfd, fdata, plen)
    implicit none
    integer, intent(IN) :: psockfd, plen
    double precision, target, intent(IN)  :: fdata(plen)
    integer(8) :: nwrite, nlen

    nlen = 8*plen
    nwrite = socket_write(psockfd, c_loc(fdata(1)), nlen)
    if (nwrite/=nlen) then
       write(6,*) "Error in writing to socket buffer"
       stop " ENDED "
    endif
  end subroutine writebuffer_dv


  ! Set of wrappers to socket_read. Read data from a socket.
  ! Args:
  !  psockfd: The id of the socket we will read from.
  !   data: The data to be read from the socket (different kinds)
  !   plen: The length of the data (not for single element reads)
  ! NB: we always need to read to a c_str, since socket_read can
  ! return
  ! before the read has completed. Then, the c_str data must be
  ! copied to
  ! the output variable

  subroutine readbuffer_cstr(psockfd, cstring, plen)
    implicit none
    integer, intent(IN) :: psockfd
    integer, intent(IN) ::  plen
    character(LEN=1,KIND=C_CHAR), intent(OUT), target ::&
         & cstring(plen)
    integer(8) :: nread, nlen, n

    nlen = plen

    nread = socket_read(psockfd, c_loc(cstring(1)), nlen)
    n = nread
    do while(nread>0 .and. n<nlen)
       nread = socket_read(psockfd, c_loc(cstring(n+1)), nlen-n)
       n = n + nread
    enddo

    if (n<nlen) then
       write(6,*) "Error in reading from socket"
       stop " ENDED "
    endif
  end subroutine readbuffer_cstr

  subroutine readbuffer_s(psockfd, fstring, plen)
    implicit none
    integer, intent(IN) :: psockfd
    integer, intent(IN) ::  plen
    character(LEN=*), target, intent(OUT)  :: fstring
    integer(8) :: n, i
    character(LEN=1,KIND=C_CHAR), target :: cstring(plen)

    call readbuffer_cstr(psockfd, cstring, plen)

    n = plen
    do i = 1,n
       fstring(i:i) = cstring(i)
    enddo
  end subroutine readbuffer_s

  subroutine readbuffer_dv(psockfd, fdata, plen)
    implicit none
    integer, intent(IN) :: psockfd, plen
    double precision, target, intent(OUT)  :: fdata(plen)
    integer(8) :: n
    character(LEN=1,KIND=C_CHAR), target :: cstring(plen*8)

    call readbuffer_cstr(psockfd, cstring, plen*8)
    n = plen*8
    call memcpy(c_loc(fdata(1)), c_loc(cstring(1)), n)
  end subroutine readbuffer_dv

  subroutine readbuffer_d(psockfd, fdata)
    implicit none
    integer, intent(IN) :: psockfd
    double precision, target, intent(OUT)  :: fdata
    integer(8) :: n
    character(LEN=1,KIND=C_CHAR), target :: cstring(8)

    call readbuffer_cstr(psockfd, cstring, 8)
    n = 8
    call memcpy(c_loc(fdata), c_loc(cstring(1)), n)
  end subroutine readbuffer_d

  subroutine readbuffer_i(psockfd, fdata)
    implicit none
    integer, intent(IN) :: psockfd
    integer, target, intent(OUT)  :: fdata
    integer(8) :: n
    character(LEN=1,KIND=C_CHAR), target :: cstring(4)

    call readbuffer_cstr(psockfd, cstring, 4)
    n = 4
    call memcpy(c_loc(fdata), c_loc(cstring(1)), n)
  end subroutine readbuffer_i
end module F90SOCKETS
