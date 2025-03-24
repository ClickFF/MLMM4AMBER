#include "../include/dprec.fh"


module ipimod

  use state
  
  implicit none

  private

  public :: init_ipi
  public :: loop_ipi

  character(len=1024),public,save :: ipi_hostname = "127.0.0.1"
  integer,public,save :: ipi_port = 31415
  integer,private,save :: socket
  logical,private,save :: isinit = .false.
  logical,private,save :: hasdata = .false.
  integer,private,save :: verbose = 0
  integer,private,save :: nstep = 0
  _REAL_,private,save :: oldbox(3)
  _REAL_,pointer,save :: msgbuffer(:) => NULL()
  type(state_rec),save :: myener   ! energy values per time step

  _REAL_,private,parameter :: AMBERELE = 18.2223d0
  _REAL_,private,parameter :: BOHRS_TO_A = 0.529177249d0
  _REAL_,private,parameter :: AU_PER_AMBER_CHARGE = 1.d0 / amberele
  _REAL_,private,parameter :: AU_PER_AMBER_ANGSTROM = 1.d0 / bohrs_to_a
  _REAL_,private,parameter :: AU_PER_AMBER_KCAL_PER_MOL = AU_PER_AMBER_CHARGE &
       & * AU_PER_AMBER_CHARGE / AU_PER_AMBER_ANGSTROM
  _REAL_,private,parameter :: AU_PER_AMBER_FORCE = AU_PER_AMBER_KCAL_PER_MOL &
       & / AU_PER_AMBER_ANGSTROM

  _REAL_,private,parameter :: AU_PER_IPI_ANGSTROM = 1.8897261d0
  
contains

  
  subroutine init_ipi()
    use F90SOCKETS
    implicit none
    integer :: inet
    inet=1
    socket=0
    call open_socket(socket, inet, ipi_port, ipi_hostname)
    isinit = .true.
    hasdata=.false.
  end subroutine init_ipi




  
  subroutine loop_ipi(ipairs)
    use F90SOCKETS
    use constants, only : DEG_TO_RAD, one, half
    use nblist, only : alpha, beta, gamma, ucell, recip
    use memory_module, only : x,ix,ih
    implicit none


    
#  include "../include/memory.h"
#  include "box.h"
#  include "ew_pme_recip.h"
#  include "ew_frc.h"
#  include "extra.h"
#  include "extra_pts.h"
#  include "../include/md.h"
#  include "nmr.h"
#ifdef MPI
#  include "mpif.h"
#  include "parallel.h"
#endif
    
    integer,intent(inout) :: ipairs(*)
    integer, parameter :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
    character(len=12) :: header
    double precision :: mtxbuf(9)
    integer :: newnfft1, newnfft2, newnfft3
    integer :: nat,i,k,ierr
    logical :: qsetup,do_list_update
    integer :: cbuf,rid
    character(len=4096) :: initbuffer
    double precision :: onefac(3)

    onefac = 0.d0

#ifdef MPI
    notdone = 1
#endif

    qsetup=.true.
    do_list_update=.true.
       
#ifdef MPI
    if ( master ) then
#endif

       call readbuffer(socket, header, MSGLEN)
       if (verbose > 0) WRITE(*,*) " Message from server: ", trim(header)
       
       if (trim(header) == "STATUS") THEN
          
          
          ! The wrapper is inquiring on what we are doing
          if (isinit) THEN
             call writebuffer(socket,"NEEDINIT    ",MSGLEN)  ! Signals that we need initialization data
             if (verbose > 1) WRITE(*,*) "    !write!=> ", "NEEDINIT    "
          else if (hasdata) THEN
             call writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces
             if (verbose > 1) WRITE(*,*) "    !write!=> ", "HAVEDATA    "
          else
             call writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
             if (verbose > 1) WRITE(*,*) "    !write!=> ", "READY       "
          endif

       
       else if (trim(header) == "INIT") THEN     ! The driver is kindly providing a string for initialization
          myener       = null_state_rec
          
          oldbox(1:3) = (/ box(1), box(2), box(3) /)
          nstep = 0
          call readbuffer(socket, rid)
          if (verbose > 0) WRITE(*,*) "    !read!=> RID: ", rid
          call readbuffer(socket, cbuf)
          if (verbose > 0) WRITE(*,*) "    !read!=> init_length: ", cbuf
          call readbuffer(socket, initbuffer, cbuf)
          if (verbose > 0) WRITE(*,*) "    !read!=> init_string: ", cbuf
          if (verbose > 10) WRITE(*,*) " Initializing system from wrapper, using ", initbuffer(1:cbuf)
          isinit=.false. ! We actually do nothing with this string, thanks anyway. Could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver

          
       else if (trim(header) == "POSDATA") THEN  ! The driver is sending the positions of the atoms. Here is where we do the calculation!



          nstep = nstep + 1

          
       
          ! Parses the flow of data from the socket
          call readbuffer(socket, mtxbuf, 9)  ! Cell matrix
          if (verbose > 1) write(*,*) "    !read!=> cell: ", mtxbuf
          !cell_h = reshape(mtxbuf, (/3,3/))
          
          mtxbuf = mtxbuf / AU_PER_IPI_ANGSTROM
          oldbox(1) = box(1)
          oldbox(2) = box(2)
          oldbox(3) = box(3)
          box(1) = mtxbuf(1)
          box(2) = mtxbuf(5) / sin(DEG_TO_RAD*gamma)
          !box(3) = sqrt(mtxbuf(9)*mtxbuf(9)+mtxbuf(7)*mtxbuf(7)+mtxbuf(8)*mtxbuf(8))
          box(3) = sqrt(mtxbuf(9)*mtxbuf(9)+mtxbuf(3)*mtxbuf(3)+mtxbuf(6)*mtxbuf(6))
          
          
          call readbuffer(socket, mtxbuf, 9)  ! Inverse of the cell matrix (so we don't have to invert it every time here)
          if (verbose > 1) write(*,*) "    !read!=> cell-1: ", mtxbuf
          !cell_ih = reshape(mtxbuf, (/3,3/))

          !write(6,'(3f12.7)')oldbox
          !write(6,'(3f12.7)')box
          

          ! 3- If the box size has changed, some parameters have to be recalc.
          if (box(1)/=oldbox(1).or.box(2)/=oldbox(2).or.box(3)/=oldbox(3)) then
             ! 3.a - Update PME cell parameters
             call fill_ucell(box(1),box(2),box(3),alpha,beta,gamma,.true.)
             
             ! 3c- Recompute grid sizes
             ! If the grid sizes change the non-bonded energy will probably be
             ! significantly different than it was in the simulation.
             call compute_nfft((box(1) + one)*half   ,newnfft1)
             newnfft1=newnfft1*2
             call compute_nfft(box(2),newnfft2)
             call compute_nfft(box(3),newnfft3)
             if ( (newnfft1/=nfft1) .or. (newnfft2/=nfft2) .or. (newnfft3/=nfft3) ) then
                write(6,'(a)') 'I-PI WARNING: nfft size would change with this box!'
                write (6,'(a,3i6)') "I-PI: New NFFTs: ",newnfft1,newnfft2,newnfft3
             endif
          endif ! box sizes have changed
          

       
          ! The wrapper uses atomic units for everything, and row major storage.
          ! At this stage one should take care that everything is converted in the
          ! units and storage mode used in the driver.
          !cell_h = transpose(cell_h)
          !cell_ih = transpose(cell_ih)
          ! We assume an upper triangular cell-vector matrix
          !volume = cell_h(1,1)*cell_h(2,2)*cell_h(3,3)
          
          call readbuffer(socket, nat)       ! The number of atoms in the cell
          if (verbose > 1) write(*,*) "    !read!=> cbuf: ", nat
          !if (nat < 0) THEN  ! Assumes that the number of atoms does not change throughout a simulation, so only does this once
          !   nat = cbuf
          !   if (verbose > 0) write(*,*) " Allocating buffer and data arrays, with ", nat, " atoms"
          !   allocate(msgbuffer(3*nat))
          !   allocate(atoms(nat,3), datoms(nat,3))
          !   allocate(forces(nat,3))
          !   allocate(friction(3*nat,3*nat))
          !   atoms = 0.0d0
          !   datoms = 0.0d0
          !   forces = 0.0d0
          !   friction = 0.0d0
          !   msgbuffer = 0.0d0
          !end if
          if ( nat /= natom ) then
             write(6,'(A,I7,A,I7,A)')"I-PI sending ",nat, &
                  & " atomic crds, but amber expected ",natom," atoms"
          end if
          if ( associated(msgbuffer) ) then
             if ( size(msgbuffer) /= nat*3 ) then
                deallocate(msgbuffer)
                allocate( msgbuffer(nat*3) )
             end if
          else
             allocate( msgbuffer(nat*3) )
          end if
          msgbuffer=0.d0
          call readbuffer(socket, msgbuffer, nat*3)
          do i = 1, nat
             do k=1,3
                x(lcrd+3*(i-1)+k-1) = msgbuffer(3*(i-1)+k) / AU_PER_IPI_ANGSTROM
             end do
             
          enddo

          !do i=4,6
          !   write(6,'(i7,3f12.7)')i,x(lcrd+3*(i-1)),x(lcrd+3*(i-1)+1),x(lcrd+3*(i-1)+2)
          !end do

          !call minrit(1,nat,1,x(lcrd))
          
          call repair_molwrap(nspm,ix(i70),ucell,recip,natom,x(lcrd))
          !call minrit(2,nat,1,x(lcrd))

          
          call wrap_molecules(nspm,ix(i70),x(lcrd))
          !call minrit(3,nat,1,x(lcrd))

          !stop
          
          !do i = 4,6
          !   write(6,'(i7,3f12.7)')i,x(lcrd+3*(i-1)),x(lcrd+3*(i-1)+1),x(lcrd+3*(i-1)+2)
          !end do
          
          if( numextra > 0 ) call local_to_global(x(lcrd),x,ix)

#ifdef MPI
          notdone = 1
          call mpi_bcast(notdone,1,mpi_integer,0,commsander,ierr)
          call mpi_bcast(x(lcrd), 3*natom, mpi_double_precision, 0, commsander, ierr)

#endif

          
          myener = null_state_rec

          call force(x,ix,ih,ipairs,x(lcrd),x(lforce),myener,myener%vir, &
               x(l96),x(l97),x(l98),x(l99), qsetup, &
               do_list_update,nstep)


          ! write(6,'(a,f18.4)')"tot  ",myener%pot%tot
          ! write(6,'(a,f18.4)')"vdw  ",myener%pot%vdw
          ! write(6,'(a,f18.4)')"elec ",myener%pot%elec
          ! write(6,'(a,f18.4)')"gb   ",myener%pot%gb
          ! write(6,'(a,f18.4)')"bond ",myener%pot%bond
          ! write(6,'(a,f18.4)')"angl ",myener%pot%angle
          ! write(6,'(a,f18.4)')"dihe ",myener%pot%dihedral
          ! write(6,'(a,f18.4)')"vdw14",myener%pot%vdw_14
          ! write(6,'(a,f18.4)')"ele14",myener%pot%elec_14
          ! write(6,'(a,f18.4)')"con  ",myener%pot%constraint
          ! write(6,'(a,f18.4)')"pol  ",myener%pot%polar
          ! write(6,'(a,f18.4)')"hbond",myener%pot%hbond
          ! write(6,'(a,f18.4)')"surf ",myener%pot%surf
          ! write(6,'(a,f18.4)')"scf  ",myener%pot%scf
          ! write(6,'(a,f18.4)')"disp ",myener%pot%disp

          call prntmd(nstep,0.d0,myener,onefac,0,.false.)
          
#ifdef MPI
          !cbuf = 0
          !call fdist(x(lforce), x(lfrctmp), myener%pot, myener%vir, cbuf)
          call xdist(x(lforce), x(lfrctmp), natom)
#endif
          hasdata = .true.
          
       
       else if (trim(header) == "GETFORCE") THEN  ! The driver calculation is finished, it's time to send the results back to the wrapper

          ! Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
          do i = 1, nat
             do k=1,3
                msgbuffer(3*(i-1)+k) = x(lforce+3*(i-1)+k-1) * AU_PER_AMBER_FORCE
                !msgbuffer(3*(i-1)+1:3*i) = frc(:,i) * AU_PER_AMBER_FORCE
             end do
          end do
          
          mtxbuf = reshape( atvir, (/ 9 /) ) * AU_PER_AMBER_FORCE * AU_PER_AMBER_ANGSTROM !* (0.1810565d0)
          !mtxbuf = reshape( molvir, (/ 9 /) ) * AU_PER_AMBER_FORCE * AU_PER_AMBER_ANGSTROM

          
          !do i = 1, 3
          !   write(6,'(a,i7,3f16.4)')"force ",i,msgbuffer(3*(i-1)+1:3*(i-1)+3)/AU_PER_AMBER_FORCE
          !end do

          !write(6,*)
          ! do i=1,natom
          !    write(6,'(i6,3f12.7)')i,x(lcrd-1+(i-1)*3+1:lcrd-1+(i-1)*3+3)
          !    write(6,'(i6,3f12.7)')i,msgbuffer(3*(i-1)+1:3*(i-1)+3)/AU_PER_AMBER_FORCE
          ! end do
          
          
          
          CALL writebuffer(socket,"FORCEREADY  ",MSGLEN)
          IF (verbose > 1) WRITE(*,*) "    !write!=> ", "FORCEREADY  "
          CALL writebuffer(socket,myener%pot%tot * AU_PER_AMBER_KCAL_PER_MOL)  ! Writing the potential
          !IF (verbose > 1) WRITE(*,*) "    !write!=> pot: ", myener%pot%tot
          CALL writebuffer(socket,natom)  ! Writing the number of atoms
          !IF (verbose > 1) WRITE(*,*) "    !write!=> nat:", nat
          CALL writebuffer(socket,msgbuffer,3*nat) ! Writing the forces
          !IF (verbose > 1) WRITE(*,*) "    !write!=> forces:", msgbuffer
          CALL writebuffer(socket,mtxbuf,9)  ! Writing the virial tensor, NOT divided by the volume
          !IF (verbose > -1) WRITE(*,*) "    !write!=> strss: ", mtxbuf
          cbuf=1 ! Size of the "extras" string
          CALL writebuffer(socket,cbuf)
          CALL writebuffer(socket,' ',1)

          
          hasdata = .false.

       else
          
          write(6,*)"Unexpected header ",trim(header)

#ifdef MPI
          notdone=0

          myener = null_state_rec

          call mpi_bcast(notdone,1,mpi_integer,0,commsander,ierr)

          call mexit(6,0)
#else
          
          stop "ENDED"
#endif
          
       end if

#ifdef MPI
    else

       nstep = nstep + 1


       qsetup = .true.
       do_list_update = .true.

       call mpi_bcast(notdone,1,mpi_integer,0,commsander,ierr)
       if ( notdone == 0 ) then
          call mexit(6,0)
       end if
       
       ! oldbox(1) = box(1)
       ! oldbox(2) = box(2)
       ! oldbox(3) = box(3)
       ! call mpi_bcast(box(1),3,mpi_double_precision,0,commsander,ierr)
       ! if (box(1)/=oldbox(1).or.box(2)/=oldbox(2).or.box(3)/=oldbox(3)) then
       !    call fill_ucell(box(1),box(2),box(3),alpha,beta,gamma,.true.)
       ! end if
       
       call mpi_bcast(x(lcrd), 3*natom, mpi_double_precision, 0, commsander, ierr)
       
       
       do i=1,3*natom
          x(lvel-1+i)   = 0.d0
          x(lvel2-1+i)  = 0.d0
          x(lforce-1+i) = 0.d0
       end do

       myener = null_state_rec
          
       call force(x,ix,ih,ipairs,x(lcrd),x(lforce),myener,myener%vir, &
            & x(l96),x(l97),x(l98),x(l99), qsetup, &
            & do_list_update,nstep)

       !write(6,*)"notdone",sanderrank,notdone
       
       if ( notdone /= 1 ) then
          call mexit(6,0)
       end if
       
       !cbuf = 0
       !call fdist(x(lforce), x(lfrctmp), myener%pot, myener%vir, cbuf)
       call xdist(x(lforce), x(lfrctmp), natom)
       
       ! write(6,*)
       ! do i=1,natom
       !    write(6,'(i6,3f12.7)')i,x(lcrd-1+(i-1)*3+1:lcrd-1+(i-1)*3+3)
       !    write(6,'(i6,3f12.7)')i,x(lforce-1+(i-1)*3+1:lforce-1+(i-1)*3+3)
       ! end do
       
       
    end if
#endif
          


  end subroutine loop_ipi






  subroutine repair_molwrap(nmol,natpermol,ucell,recip,nat,crds)
    implicit none

    integer,intent(in)   :: nmol
    integer,intent(in)   :: natpermol(nmol)
    _REAL_,intent(in)    :: ucell(3,3), recip(3,3)
    integer,intent(in)   :: nat
    _REAL_,intent(inout) :: crds(3,nat)

    integer :: imol,iat,jat
    integer,allocatable :: firstatom(:)
    _REAL_ :: r2,r2min,dcrd(3),fcrd(3),dmin(3)

    allocate(firstatom(nmol+1))
    firstatom=1
    do imol=1,nmol
       firstatom(imol+1) = firstatom(imol) + natpermol(imol)
       !write(6,*)firstatom(imol),firstatom(imol+1),firstatom(imol)+1,firstatom(imol+1)-1
    end do
    
  ! double fcrd[3] = { 0,0,0 };
  ! dgemv_( "T", &m, &m, &alpha, recip.data(), &m, crd, &inc, &zero, fcrd, &inc );
  ! fcrd[0] -= anint( fcrd[0] );
  ! fcrd[1] -= anint( fcrd[1] );
  ! fcrd[2] -= anint( fcrd[2] );
  ! dgemv_( "N", &m, &m, &alpha, ucell.data(), &m, fcrd, &inc, &zero, crd, &inc );

    do imol=1,nmol
       !write(6,*)"imol",imol,
       do iat=firstatom(imol)+1,firstatom(imol+1)-1
          r2min = 1.d+30
          do jat=firstatom(imol),iat-1
             dcrd = crds(1:3,iat)-crds(1:3,jat)
             fcrd = matmul( dcrd, recip )
             fcrd = fcrd - anint(fcrd)
             dcrd = matmul( ucell, fcrd )
             r2 = dcrd(1)*dcrd(1) + dcrd(2)*dcrd(2) + dcrd(3)*dcrd(3)
             if ( r2 < r2min ) then
                r2min = r2
                dmin = dcrd + crds(1:3,jat)
             end if
             !write(6,'(3i6,4F12.4)')imol,iat,jat,sqrt(r2),dmin
          end do
          crds(1:3,iat) = dmin
       end do
    end do
       
  end subroutine repair_molwrap








  
end module ipimod
