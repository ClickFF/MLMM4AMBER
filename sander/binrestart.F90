! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

!! Module for generating binary restart-type output in NetCDF format
!! Developed by Dan Roe
!! 2010-01-10

module binrestart
   private

   integer, save :: coordVID, velocityVID, cellAngleVID, cellLengthVID
   integer, save :: timeVID, remd_values_var_id
   integer, save :: remd_indices_var_id, remd_groups_var_id
   integer, save :: remd_dimension_var_id, remd_types_var_id
#ifdef MPI
   integer, save :: repidx_var_id, crdidx_var_id
#endif

   public write_nc_restart, &
          read_nc_restart_box, &
          read_nc_restart, &
          read_nc_restart_extents, &
          read_nc_remd_dimension, &
          read_nc_remd_types
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write Netcdf restart file.
!-------------------------------------------------------------------
!     --- WRITE_NC_RESTART ---
!-------------------------------------------------------------------
!     Write Netcdf Restart file with given filename and title. 
!     owrite indicates overwrite status (N is no overwrite), natom 
!     is the # atoms, ntb>0 indicates presence of box coords, first
!     indicates the file should be created and set-up, Coords and 
!     Velo are the coordinates and velocities, temp0 is the current
!     temperature and Time is the current simulation time. If hasV
!     is false, no velocities will be written (for e.g. during min)
subroutine write_nc_restart(filename,title,owrite,natom,ntb,first,Coords,Velo,&
                            Time,hasVin&
#ifdef MPI
                            , temp0, rem, remd_dimension, remd_types, group_num &
                            , replica_indexes, stagid, remd_repidx, remd_crdidx &
                            , solvph, solve &
#endif
                           )
   use AmberNetcdf_mod
#ifdef BINTRAJ
   use netcdf
   use nblist, only: a,b,c,alpha,beta,gamma
#endif
#ifdef MPI
   use sgld, only: trxsgld
#endif

   implicit none
   ! Input variables
   character(len=*), intent(in)      :: filename
   character(len=*), intent(in)      :: title
   character, intent(in)             :: owrite
   integer, intent(in)               :: natom, ntb
   logical, intent(in)               :: first
   _REAL_, dimension(*), intent(in)  :: Coords, Velo
   _REAL_, intent(in)                :: Time
   logical, intent(in)               :: hasVin
#  ifdef MPI
   _REAL_, intent(in)                :: temp0, solvph, solve
   integer, intent(in)               :: rem, remd_dimension
   integer, intent(in), dimension(:) :: remd_types, group_num, replica_indexes
   integer, intent(in)               :: stagid, remd_repidx, remd_crdidx 
#  endif
#ifdef BINTRAJ
   ! Local vars
   integer :: ncid, natom3
   logical :: has1DRemdValues = .false.
   integer :: frcVID ! dummy variable
#  ifdef MPI
   double precision, dimension(remd_dimension) :: remd_values
   integer :: i
   has1DRemdValues = (rem.gt.0)
#  endif
   if (first) then
      ! If first call, create the file and set up all dimensions and vars
      ! owrite status code: 'N', 'O', 'R', 'U' = new, old, replace, unknown
      ! sander flag -O='R', -A='U', default='N'
      if (NC_create(filename, owrite, .true., natom, .true.,&
                    hasVin, (ntb.gt.0), has1DRemdValues, .true., .false., &
                    title, ncid, timeVID, coordVID, velocityVID, frcVID, &
                    cellLengthVID, cellAngleVID, remd_values_var_id)) call mexit(6,1)
#     ifdef MPI
      ! REMD indices
      if (NC_defineRemdIndices(ncid, remd_dimension, remd_indices_var_id,&
                               repidx_var_id, crdidx_var_id, &
                               remd_types, .true., (rem.ne.0), (rem.eq.-1), &
                               remd_groupsVID=remd_groups_var_id,&
                               remd_valuesVID=remd_values_var_id,&
                               remd_typesVID=remd_types_var_id )) &
         call mexit(6,1)
#     endif
   else
      ! If not the first call, just reopen the existing file
      if (NC_openWrite(filename, ncid)) then
        write (6,'(a)') 'write_nc_restart(): Could not open restart'
        call mexit(6,1)
      endif
   endif

   natom3 = natom * 3
   ! Write time
   call checkNCerror(nf90_put_var(ncid, timeVID, Time), 'write time')
   ! Write coords
   call checkNCerror(nf90_put_var(ncid,coordVID, Coords(1:natom3), &
                  start = (/ 1, 1 /), count = (/ 3, natom /) ), 'write atom coords')
   ! Write velocities TODO: Should this check velocityVID instead?
   if (hasVin) then
      call checkNCerror(nf90_put_var(ncid,velocityVID, Velo(1:natom3), &
                     start = (/ 1, 1 /), count = (/ 3, natom /) ), 'write velocities')
   endif
   ! Write box information
   if (ntb > 0) then
      call checkNCerror(nf90_put_var(ncid,cellLengthVID, &
              (/ a, b, c /), start = (/ 1 /), count = (/ 3 /) ), 'write cell lengths')
      call checkNCerror(nf90_put_var(ncid,cellAngleVID, &
              (/ alpha,beta,gamma /), start = (/ 1 /), count = (/ 3 /) ), &
           'write cell angles')
   endif
#  ifdef MPI
  ! Write replica temperature, solvent pH, redox potential, indices
   if (rem.ne.0) then
      ! Write overall coordinate and replica index
      call checkNCerror(nf90_put_var(ncid, repidx_var_id, remd_repidx), &
                        'write overall replica index')
      call checkNCerror(nf90_put_var(ncid, crdidx_var_id, remd_crdidx), &
                        'write overall coordinate index')
      call checkNCerror(nf90_put_var(ncid, remd_indices_var_id, &
                        replica_indexes(:), &
                        start = (/ 1 /), count = (/ remd_dimension /)), &
                        'write replica index for each dimension')
      ! multi-D remd: Store indices of this replica in each dimension
      if (rem .eq. -1) then
         call checkNCerror(nf90_put_var(ncid, remd_groups_var_id, group_num(:), &
                           start = (/ 1 /), count = (/ remd_dimension /)), &
                           'write replica group for each dimension')
         ! Preparing remd_values vector
         do i = 1, remd_dimension
           if (remd_types(i) == 1) then
             remd_values(i) = temp0
           else if (remd_types(i) == 3) then
             remd_values(i) = replica_indexes(i)
           else if (remd_types(i) == 4) then
             remd_values(i) = solvph
           else if (remd_types(i) == 5) then
             remd_values(i) = solve
           end if
         end do
         call checkNCerror(nf90_put_var(ncid, remd_values_var_id, remd_values(:), &
                           start = (/ 1 /), count = (/ remd_dimension /)), &
                           'write replica values to be restarted for each dimension')
      endif
      if (has1DRemdValues) then
        if (trxsgld) then
            call checkNCerror(nf90_put_var(ncid, remd_values_var_id, REAL(stagid)), &
                              'write SGLD replica index')
        else if (rem.eq.1) then
            call checkNCerror(nf90_put_var(ncid, remd_values_var_id, temp0), 'write temp0')
        else if (rem.eq.4) then
            call checkNCerror(nf90_put_var(ncid, remd_values_var_id, solvph), 'write solvph')
        else if (rem.eq.5) then
            call checkNCerror(nf90_put_var(ncid, remd_values_var_id, solve), 'write solve')
        endif
      endif
   endif
#  endif
   ! Close restart file       
   call NC_close(ncid)
#else
   write(0,*) 'sander has been compiled without netcdf support'
   call mexit(6,1)
#endif
end subroutine write_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read box information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART_BOX ---
!-------------------------------------------------------------------
!     Read box information from the Netcdf Restart file with 
!     specified filename.
!     The box read is called from load_ewald_info() in ew_setup.f
!     and is separate from the coord/velocity read since the box
!     information is needed to set up certain ewald parameters.
subroutine read_nc_restart_box(filename,a,b,c,alpha,beta,gamma)
   use AmberNetcdf_mod 
   implicit none

   character(len=*), intent(in) :: filename
   _REAL_, intent(out) :: a,b,c,alpha,beta,gamma
   ! local
   integer ncid
#ifdef BINTRAJ
   if (NC_openRead(filename, ncid)) call mexit(6,1)
   if (NC_readRestartBox(ncid,a,b,c,alpha,beta,gamma)) call mexit(6,1)
   call NC_close(ncid)
   write(6,'(a)') '| NetCDF restart box info found'
#else
   write(0,*) 'sander has been compiled without netcdf support'
   call mexit(6,1)
#endif
end subroutine read_nc_restart_box

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coord/velocity information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART ---
!-------------------------------------------------------------------
!     Read coordinates and velocities from the Netcdf Restart file 
!     with specified filename. This is called from getcor.f. Title
!     will be read in and set. 
!     ntx specifies whether coords and velo or just coords will be 
!     read. ntx=1 means read coords only, and ntx=5 means read 
!     coords and velocities.
!     parmatoms is the expected number of atoms in the restart. If -1, we assume
!        it is not set (if the prmtop hasn't been read yet, such as for the API).
!        In this case, parmatoms is set to the number of atoms in the inpcrd
!        file before returning.
!     Coords and Velo are the coordinates and 
!     velocities, remd_values contains temperature, pH and/or redox potential
!     (if present) and Time is the time.
!     NOTE: Box info is not read here; it is obtained using 
!     read_nc_restart_box in load_ewald_info.
subroutine read_nc_restart(filename,title,ntx,parmatoms,Coords,Velo,remd_values,&
                           remd_values_dim,Time)
   use AmberNetcdf_mod
#ifdef BINTRAJ
   use netcdf
#endif
   use constants, only : NO_INPUT_VALUE_FLOAT
   implicit none

   character(len=*), intent(in)      :: filename
   character(len=80), intent(out)    :: title
   integer, intent(in)               :: ntx, remd_values_dim
   integer, intent(in out)           :: parmatoms
   _REAL_, dimension(*), intent(out) :: Coords, Velo
   _REAL_, intent(out)               :: Time
   _REAL_, dimension(:), intent(out) :: remd_values
#ifdef BINTRAJ
   integer :: ncid, ncatom, ncatom3

   ! ---=== Open file
   if (NC_openRead(filename, ncid)) call mexit(6,1)
   ! Setup restart: title, coordVID, velocityVID, time, remd_values_var_id
   if (NC_setupRestart(ncid, title, ncatom, coordVID, velocityVID, &
                       remd_values_var_id, Time)) call mexit(6,1)
   ! Check that number of atoms matches
   if (ncatom /= parmatoms .and. parmatoms /= -1) then
      write(6,'(2x,a)') "FATAL: NATOM mismatch in restart and topology files."
      call mexit(6,1)
   else if (parmatoms == -1) then
      parmatoms = ncatom
   endif
   ncatom3 = ncatom * 3
   ! ---=== Get Coords
   if (NC_error(nf90_get_var(ncid, coordVID, Coords(1:ncatom3), &
                             start = (/ 1, 1 /), count = (/ 3, ncatom /)),&
                'reading restart coordinates')) call mexit(6,1)
   ! ---=== Get velocities
   !        ntx=1 No Velocity Read 
   !        ntx=5 Read Velocity
   if (ntx == 5) then
      if (velocityVID .eq. -1) then
        write(6,'(2x,a)') "FATAL: ntx=5 specified but no velocities in INPCRD"
        call mexit(6,1)
      endif
      if (NC_error(nf90_get_var(ncid, velocityVID, Velo(1:ncatom3), &
                                start = (/ 1, 1 /), count = (/ 3, ncatom /)),&
                   'reading restart velocities')) call mexit(6,1)
   endif
   ! ---=== Replica remd_values
   if (remd_values_var_id .ne. -1) then
      if (NC_error(nf90_get_var(ncid, remd_values_var_id, remd_values, &
                                start = (/ 1 /), count = (/ remd_values_dim /)),&
                   "read_nc_restart(): Getting restart remd_values")) call mexit(6,1)
   else
     remd_values(:)=NO_INPUT_VALUE_FLOAT
   endif 

   ! NOTE: TO BE ADDED
   !labelDID;
   !int cell_spatialDID, cell_angularDID;
   !int spatialVID, cell_spatialVID, cell_angularVID;
  
   ! ---=== Close file
   call NC_close(ncid)
#else
   write(0,*) 'sander has been compiled without netcdf support'
   call mexit(6,1)
#endif
end subroutine read_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read remd_dimension information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_REMD_DIMENSION ---
!-------------------------------------------------------------------
!     Read remd_dimension from the Netcdf Restart file with specified 
!     filename. Title will be read in and set.
subroutine read_nc_remd_dimension(filename,title,remd_dimension)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only: NC_openRead, NC_setupRemdDimension, NC_error, NC_close
#endif
  implicit none
  ! Formal Arguments
  character(len=*), intent(in)                :: filename
  character(len=80), intent(out)              :: title
  integer, intent(out) :: remd_dimension
#ifdef BINTRAJ
  ! Local variables
  integer :: ncid

  ! ---=== Open file
  if (NC_openRead(filename, ncid)) then
    write(6,'(a)') "read_nc_remd_dimension(): Could not open restart file."
    call mexit(6,1)
  endif
  ! Get remd_dimension
  remd_dimension_var_id = NC_setupRemdDimension(ncid,remd_dimension)
  if (remd_dimension_var_id.eq.-1) remd_dimension = -1
  ! ---=== Close file
  call NC_close(ncid)
#else
  write(0,*) 'sander has been compiled without netcdf support'
  call mexit(6,1)
#endif
end subroutine read_nc_remd_dimension

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read remd_types information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_REMD_TYPES ---
!-------------------------------------------------------------------
!     Read remd_types from the Netcdf Restart file with specified 
!     filename. Title will be read in and set.
subroutine read_nc_remd_types(filename,title,remd_types,remd_dimension)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only: NC_openRead, NC_setupRemdTypes, NC_error, NC_close
#endif
  implicit none
  ! Formal Arguments
  character(len=*), intent(in)                    :: filename
  character(len=80), intent(out)                  :: title
  integer, intent(in)                             :: remd_dimension
  integer, dimension(remd_dimension), intent(out) :: remd_types
#ifdef BINTRAJ
  ! Local variables
  integer :: ncid

  ! ---=== Open file
  if (NC_openRead(filename, ncid)) then
    write(6,'(a)') "read_nc_remd_types(): Could not open restart file."
    call mexit(6,1)
  endif
  ! Setup restart: remd_types_var_id
  remd_types_var_id = NC_setupRemdTypes(ncid)
  if (remd_types_var_id.eq.-1) then
    remd_types(1) = -1
  else
    ! Get remd_types
    if (NC_error(nf90_get_var(ncid, remd_types_var_id, remd_types, &
                              start=(/ 1 /), count=(/ remd_dimension /)),&
                 "read_nc_remd_types(): Getting remd_types")) call mexit(6,1)
  end if
  ! ---=== Close file
  call NC_close(ncid)  
#else
  write(0,*) 'sander has been compiled without netcdf support'
  call mexit(6,1)
#endif
end subroutine read_nc_remd_types

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read extents information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART_EXTENTS ---
!-------------------------------------------------------------------
!     Read extents from the Netcdf Restart file 
!     with specified filename.
subroutine read_nc_restart_extents(filename, extents)
   use AmberNetcdf_mod
#ifdef BINTRAJ
   use netcdf
#endif
   use constants, only : NO_INPUT_VALUE_FLOAT
   implicit none

   character(len=*), intent(in)        :: filename
   _REAL_, dimension(3,2), intent(out) :: extents
#ifdef BINTRAJ
   _REAL_, dimension(:), allocatable :: Coords
   _REAL_  :: Time
   integer :: ncid, ncatom, ncatom3, ierror, j, i
   character(len=256) :: title

   title = ''
   ! ---=== Open file
   if (NC_openRead(filename, ncid)) call mexit(6,1)
   ! Setup restart: title, coordVID, velocityVID, time, remd_values_var_id
   if (NC_setupRestart(ncid, title, ncatom, coordVID, velocityVID, &
                       remd_values_var_id, Time)) call mexit(6,1)
   ! Check that number of atoms matches
   ncatom3 = ncatom * 3
   ! Allocate the coordinate array
   allocate(Coords(ncatom3), stat=ierror)
   REQUIRE( ierror == 0 )
   ! ---=== Get Coords
   if (NC_error(nf90_get_var(ncid, coordVID, Coords(1:ncatom3), &
                             start = (/ 1, 1 /), count = (/ 3, ncatom /)),&
                'reading restart coordinates')) call mexit(6,1)
   ! ---=== Close file
   call NC_close(ncid)
   ! Now go through and find the max/min of x, y, and z
   extents(1, 1) = Coords(1) ! Min X
   extents(2, 1) = Coords(2) ! Min Y
   extents(3, 1) = Coords(3) ! Min Z
   extents(1, 2) = Coords(1) ! Max X
   extents(2, 2) = Coords(2) ! Max Y
   extents(3, 2) = Coords(3) ! Max Z
   j = 4
   do i = 2, ncatom
      if ( Coords(j) < extents(1, 1) ) extents(1, 1) = Coords(j)
      if ( Coords(j) > extents(1, 2) ) extents(1, 2) = Coords(j)
      j = j + 1
      if ( Coords(j) < extents(2, 1) ) extents(2, 1) = Coords(j)
      if ( Coords(j) > extents(2, 2) ) extents(2, 2) = Coords(j)
      j = j + 1
      if ( Coords(j) < extents(3, 1) ) extents(3, 1) = Coords(j)
      if ( Coords(j) > extents(3, 2) ) extents(3, 2) = Coords(j)
      j = j + 1
   end do
   ! Deallocate our work array
   deallocate(Coords, stat=ierror)
   REQUIRE( ierror == 0 )
#else
   write(0,*) 'sander has been compiled without netcdf support'
   call mexit(6,1)
#endif
end subroutine read_nc_restart_extents

end module binrestart
