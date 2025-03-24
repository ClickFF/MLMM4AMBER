! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

!**************************************************************************************************
!
! Module: constante
!
! Description: This module handles all of the constant Redox potential capabilities in Amber.
!              It was adapted from the constantph module initially written on Sander by John Mongan.
!
!              Adapted here by Vinicius Wilian D. Cruzeiro
!
!**************************************************************************************************
module constante

use random

! Private by default
private

#  include "dyne.h"
integer, save :: STATEINF_FLD_C, TITR_RES_C, TITR_STATES_C, ATOM_CHRG_C

! Public variable declaration
_REAL_, public, allocatable :: chrgdat_e(:)
_REAL_, public  :: ce_intdiel
integer, public :: ce_igb
integer, public :: cefirst_sol
integer, public :: relaxations_e
integer, public :: tot_relax_e
integer, public :: mccycles_e

! Private variable declaration
type(const_e_info), allocatable :: stateinf(:)
type(rand_gen_state), save :: ce_gen

character(len=40), allocatable :: resname(:)

_REAL_, allocatable  :: statene(:)
_REAL_, allocatable  :: eo_corr(:)
integer, allocatable :: resstate(:)
integer, allocatable :: iselres(:)
integer :: iselstat(0:2)
integer, allocatable :: eleccnt(:)
integer :: trescnt
logical :: ce_success
logical :: first

integer, allocatable :: mobile_atoms(:)
_REAL_, allocatable  :: frac_pop(:,:)

!  Variable descriptions
!
!  ce_gen        : random number generator for constant Redox potential MC
!  stateinf       : information about titratable residues. See dyne.h
!  chrgdat_e        : partial atomic charges for every state of every titr. res.
!  statene        : relative energy of every state for each titr. res.
!  eo_corr       : contains relative standard redox potential values
!  resstate       : list of current oxidation states for each titr. res.
!  iselres        : array of selected, titrating residues
!  iselstat       : array of selected, proposed, states for selected titr. res.
!  eleccnt        : number of electrons for each state of each titr. res.
!  trescnt        : number of titrating residues
!  cefirst_sol   : atom number of the first solvent residue
!  ce_igb        : igb model to use for MC evals. in explicit solvent Redox potential MD
!  ce_success    : did the electron exchange attempt succeed or not
!  ce_intdiel    : Internal dielectric to use for MC moves in explicit solvent
!  first          : is this the first pass through? (for CEOUT printing)
!  relaxations_e    : how many times we have run water relaxations (for timing)
!  tot_relax_e      : how many relaxations we've done in total (for timing)
!  mccycles_e       : The number of MC steps that will be attempted at once in
!                   explicit solvent CEMD
!  frac_pop       : Track the populations of each state for each residue when
!                   statistics for a cycle of MC attempts in explicit solvent

! Make necessary subroutines public
public cnsteinit, cnstebeginstep, cnsteendstep,  &
       cnstewriterestart, cnstewrite, cnste_explicitmd, cnsteread, &
       total_reduction, cnste_finalize, &
#ifdef MPI
       cnste_bcast, &
#endif
       cnste_zero, cnstereadlimits

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Allocate constant E variables
subroutine cnsteallocate()
   implicit none

   integer :: alloc_failed

   allocate(chrgdat_e(0:ATOM_CHRG_C-1),    &
            stateinf(0:TITR_RES_C-1),    &
            resname(0:TITR_RES_C),       &
            statene(0:TITR_STATES_C-1),  &
            eo_corr(0:TITR_STATES_C-1),  &
            resstate(0:TITR_RES_C-1),    &
            iselres(0:TITR_RES_C),       &
            eleccnt(0:TITR_STATES_C-1),  &
            stat=alloc_failed )

   if (alloc_failed .ne. 0) then
     write(6,'(a)') 'Error in constant E allocation!'
     call mexit(6, 1)
   end if

end subroutine cnsteallocate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Zero all of the constant Redox potential arrays
subroutine cnste_zero()
   use constants, only : ZERO, ONE
   implicit none

   ! Allocate constant pH variables
   call cnsteallocate()

   ! I zero out all of the arrays at first because when the cerestrt is written,
   ! it simply writes the whole namelist, but unused values would not otherwise
   ! be initialized.

   ! NULL_CE_INFO is a zero-filled const_e_info type. See dyne.h

   stateinf(0:TITR_RES_C-1)                  = NULL_CE_INFO
   chrgdat_e(0:ATOM_CHRG_C-1)                  = 0
   resstate(0:TITR_RES_C-1)                  = 0
   iselres(0:2)                              = 0
   iselstat(0:2)                             = 0
   statene(0:TITR_STATES_C-1)                = ZERO
   eo_corr(0:TITR_STATES_C-1)               = ZERO
   eleccnt(0:TITR_STATES_C-1)                = 0
   resname(0:TITR_RES_C)                     = ' '
   relaxations_e                               = 0
   tot_relax_e                                 = 0
   ce_igb                                   = 0
   ce_intdiel                               = ONE
   cefirst_sol                              = 0

end subroutine cnste_zero

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read the constant E limits for allocating variables from the CEIN
!+ file, if this information is available
subroutine cnstereadlimits()
   use file_io_dat, only : CNSTE_UNIT, cein
   implicit none

   integer :: ntres, ntstates, natchrg, maxh, ifind

   namelist /cnstphe_limits/ ntres, ntstates, natchrg, maxh

   ! Set default values
   ntres = 50
   ntstates = 200
   natchrg = 1000

   ! Open the cein file and read the cnstphe_limits namelist
   call amopen(CNSTE_UNIT, cein, 'O', 'F', 'R')
   call nmlsrc('cnstphe_limits', CNSTE_UNIT, ifind)
   if (ifind .ne. 0) then        ! Namelist found. Read it:
     read (CNSTE_UNIT, nml=cnstphe_limits)
     close (CNSTE_UNIT)
   end if

   ! Setting limits
   STATEINF_FLD_C = 5
   TITR_RES_C = ntres
   TITR_STATES_C = ntstates
   ATOM_CHRG_C = natchrg

end subroutine cnstereadlimits

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read the CEIN file and initialize arrays
subroutine cnsteread(charge)
   use constants, only : AMBER_ELECTROSTATIC, ZERO, ONE
   use file_io_dat, only : CNSTE_UNIT, cein, CE_DUMP_UNIT, ce_dump, owrite, &
                           write_dump_e, target_e
   use memory_module, only : natom
   implicit none
#  include "extra.h"
#  include "../include/md.h"

   ! Variable Descriptions:
   !
   ! Passed:
   !  charge         : partial charge array
   !
   ! Internal:
   !  itres          : residue loop pointer
   !  res            : residue pointer
   !  iatom          : atom pointer
   !  icumstat       : sum of states
   !  icumchrg       : sum of charge
   !
   ! Common memory:
   !  gbgamma        : (md.h) OBC GB parameter
   !  gbbeta         : (md.h) OBC GB parameter
   !  gbalpha        : (md.h) OBC GB parameter

   _REAL_, intent(inout) :: charge(*)
   _REAL_ :: chrgdat(0:ATOM_CHRG_C-1)

   integer  :: itres, res, iatom, icumstat, icumchrg, ier, i

   ! This subroutine reads in and initializes the constant Redox potential data. It
   ! multiplies in the AMBER_ELECTROSTATIC factor and adjusts the initial charge
   ! array for the initial oxidation states.

   namelist /cnste/ stateinf, resstate, eleccnt, chrgdat, statene, eo_corr, &
                     trescnt, resname, cefirst_sol, ce_igb, ce_intdiel

   icumstat = -1
   icumchrg = 0
   res = 0

   write(6,'(a,a)') '|reading charge increments from file: ', cein
   call amopen(CNSTE_UNIT, cein, 'O', 'F', 'R')
   read(CNSTE_UNIT, nml=cnste)
   do iatom = 0, ATOM_CHRG_C-1
      chrgdat(iatom) = chrgdat(iatom) * AMBER_ELECTROSTATIC
   end do
   close(CNSTE_UNIT)

   ! Set initial charges
   do itres = 0, trescnt - 1
      do iatom = 0, stateinf(itres)%num_atoms - 1
         charge(iatom + stateinf(itres)%first_atom) &
             = chrgdat(stateinf(itres)%first_charge + iatom + resstate(itres)* &
               stateinf(itres)%num_atoms)
      end do
   end do

   chrgdat_e = chrgdat

   ! Set proper GB stuff
   if ( ce_igb .eq. 2 ) then
      gbgamma = 2.90912499999d0
      gbbeta  = ZERO
      gbalpha = 0.8d0
   else if ( ce_igb .eq. 5 ) then
      gbgamma = 4.85d0
      gbbeta  = 0.8d0
      gbalpha = ONE
   end if

   ! Overflow checking
   if (trescnt > TITR_RES_C) then
      write(6,*) 'Too many titrating residues. Alter ntres in the &
                 &cnstphe_limits namelist'
      call mexit(6,1)
   end if
   do itres = 0, trescnt - 1
      if (stateinf(itres)%first_state > icumstat) then
         icumstat = stateinf(itres)%first_state
         res = itres
      end if
   end do

   icumstat = stateinf(res)%first_state + stateinf(res)%num_states
   icumchrg = stateinf(res)%first_charge + stateinf(res)%num_atoms * &
              stateinf(res)%num_states

   if (icumstat > TITR_STATES_C) then
      write(6,*) 'Too many titrating states. Alter ntstates in the &
                 &cnstphe_limits namelist'
      call mexit(6,1)
   end if
   if (icumchrg > ATOM_CHRG_C) then
      write(6,*) 'Too much charge data. Alter natchrg in the &
                 &cnstphe_limits namelist'
      call mexit(6,1)
   end if

   if (icnste > 1) then
      write(6,'(a,i4,a,i4,a)') '| Attempting ', trescnt, &
               ' MC oxidation changes every ', ntcnste, ' steps.'
      allocate(mobile_atoms(natom), stat=ier)
      REQUIRE (ier == 0)
      mobile_atoms(1:cefirst_sol-1) = 0
      mobile_atoms(cefirst_sol:natom) = 1
      ! Only dump to separate file if we have more than 1 MC cycle
      if (mccycles_e == 1 .or. .not. master) write_dump_e = .false.
      if (write_dump_e) then
         call setup_pop_arry ! only master needs this
         call amopen(CE_DUMP_UNIT, ce_dump, owrite, 'F', 'W')
         write(CE_DUMP_UNIT, '("Residue Number",6x)', advance='NO')
         do i = 0, max_states() - 1
            write(CE_DUMP_UNIT, '("State",x,i2,2x)', advance='NO') i
         end do
         write(CE_DUMP_UNIT, '()')
         write(CE_DUMP_UNIT, '()')
      end if
   else
      ! Never write all data in implicit solvent (all data is already written to
      ! ceout file)
      write_dump_e = .false.
   end if

   target_e = solve

   ! Set intdiel here. It's not used anywhere except GB, so it doesn't matter if
   ! we set it here

   intdiel = ce_intdiel

end subroutine cnsteread

#ifdef MPI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Broadcasts all of the constant Redox potential information to each thread
subroutine cnste_bcast(ierr)
   use memory_module, only : natom
   implicit none
#include "parallel.h"
#include "../include/md.h"
#include "extra.h"
 include 'mpif.h'
   integer, intent(out) :: ierr

   ! Slave nodes have not been allocated
   call mpi_bcast(STATEINF_FLD_C, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(TITR_RES_C, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(TITR_STATES_C, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ATOM_CHRG_C, 1, MPI_INTEGER, 0, commsander, ierr)
   if (.not. master) then
     call cnsteallocate()
   end if

   ! Broadcast everything
   call mpi_bcast(stateinf, TITR_RES_C*STATEINF_FLD_C, MPI_INTEGER, 0, &
                  commsander, ierr)
   call mpi_bcast(trescnt, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(chrgdat_e, ATOM_CHRG_C, MPI_DOUBLE_PRECISION, 0, &
                  commsander, ierr)
   call mpi_bcast(resstate, TITR_RES_C, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(statene, TITR_STATES_C, MPI_DOUBLE_PRECISION, 0, &
                  commsander, ierr)
   call mpi_bcast(eo_corr, TITR_STATES_C, MPI_DOUBLE_PRECISION, 0, &
                  commsander, ierr)
   call mpi_bcast(eleccnt, TITR_STATES_C, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(resname, 40*(TITR_RES_C+1), MPI_CHARACTER, 0, &
                  commsander, ierr)
   call mpi_bcast(cefirst_sol, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ce_igb, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(mccycles_e, 1, MPI_INTEGER, 0, commsander, ierr)
   ! Now we have to allocate and broadcast our mobile atoms array
   if (icnste > 1) then
      if (.not. allocated(mobile_atoms)) &
         allocate(mobile_atoms(natom), stat=ierr)
      REQUIRE (ierr == 0)
      call mpi_bcast(mobile_atoms, natom, MPI_INTEGER, 0, commsander, ierr)
   end if

end subroutine cnste_bcast
#endif /* MPI */

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Set up an array to track populations of ea state of ea residue for mc cycles
subroutine setup_pop_arry

   implicit none

   integer :: ierror
   integer :: maxstates

   ! Make this subroutine re-usable---deallocate the array if it's already
   ! allocated
   if (allocated(frac_pop)) deallocate(frac_pop)

   maxstates = max_states()

   ! Index from zero, like the rest of our arrays here
   allocate(frac_pop(0:trescnt-1, 0:maxstates-1), stat=ierror)

   REQUIRE( ierror == 0 )

   return

end subroutine setup_pop_arry

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find the most states in any titratable residue
integer function max_states()

   implicit none

   integer i

   max_states = stateinf(0)%num_states

   do i = 1, trescnt - 1
      if (stateinf(i)%num_states > max_states) &
         max_states = stateinf(i)%num_states
   end do

   return

end function max_states

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine cnsteinit(x, ig)
   implicit none

#ifdef MPI
   include 'mpif.h'
#  include "parallel.h"
#endif

   ! Variable Descriptions:
   !
   ! Passed:
   !  x              : global position array
   !  ig             : random seed -- used for initializing ce_gen
   !
   ! Local:
   !  res, state     : descriptive pointers
   !  atom, i        : (descriptive) pointers
   !  used_atoms     : memory for which titr. electrons have been found
   !  is_used        : whether the current electron has been found yet or not


   _REAL_, intent(in)  :: x(*)
   integer, intent(in) :: ig

   integer :: seed
#ifdef MPI
   integer :: ierr
#endif

   ! Every thread in parallel MUST have an identical random number stream.
   ! However, setting ig=-1 in the input file forces desynchronization of the
   ! random numbers. For simplicity's sake, we'll just broadcast the master's
   ! "ig" value (but we don't want to clobber it, so copy it here)

   seed = ig
#ifdef MPI
   call mpi_bcast(seed, 1, mpi_integer, 0, commsander, ierr)
#endif
   call amrset_gen(ce_gen, seed)

   ! this is our first pass through
   first = .true.

end subroutine cnsteinit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Select random new state for random titrating residue, set charges
subroutine cnstebeginstep(dcharge)
   implicit none

   ! Variable descriptions:
   !
   ! Passed:
   !  dcharge        : charge array for the proposed state
   !
   ! Internal:
   !  randval        : random value from random number generator

   _REAL_, intent(inout) :: dcharge(1:*)

   _REAL_   :: randval

   ! This subroutine will randomly select what to titrate. There is a 25% chance
   ! that we will attempt a multi-site move among neighboring residues. This is
   ! done so we can at times observe a electron transfer when that would otherwise
   ! be a very very rare event. The code only currently supports a 2-site move,
   ! but iselres/iselstat and associated loops have been generalized to an
   ! arbitrary number of sites.
   ! First determine if we will do a multisite move. Then randomly select a res
   ! to titrate. If we tried a multisite move, see if our selected site has any
   ! neighbors, then select from among those to titrate also. Update the dcharge
   ! array for the proposed states. If we have no neighbors, don't do multisite.

   iselres(0) = 1

   ! select a random residue
   call amrand_gen(ce_gen, randval)
   iselres(1) = int((randval * 0.9999999d0)*trescnt)

   ! select a different, random state
   call amrand_gen(ce_gen, randval)
   iselstat(1) = int((randval * 0.9999999d0) * &
                     (stateinf(iselres(1))%num_states-1))
   if (iselstat(1) .ge. resstate(iselres(1))) &
      iselstat(1) = iselstat(1)+1
   call cnsteupdatechrg(dcharge, iselres(1), iselstat(1))

end subroutine cnstebeginstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Same as cnstebeginstep, but for explicit solvent, so we can do multi-MC
subroutine excnstebeginstep(dcharge, done_res, num_done)
   implicit none

   ! Variable descriptions:
   !
   ! Passed:
   !  dcharge        : charge array for the proposed state
   !  done_res       : array for which residues have been done
   !  num_done       : how many residues have been finished
   !
   ! Internal:
   !  i              : loop counters
   !  randval        : random value from random number generator
   !  holder_res     : holder variable so we insert selected residue in order

   _REAL_, intent(inout)  :: dcharge(1:*)
   integer, intent(inout) :: done_res(1:TITR_RES_C)
   integer, intent(in)    :: num_done

   integer  :: i, holder_res(1:TITR_RES_C)
   _REAL_   :: randval

   ! This subroutine will randomly select what to titrate, as long as it has.
   ! not yet been titrated. That is the only difference from cnstebeginstep

   iselres(0) = 1

   ! select a random residue
   call amrand_gen(ce_gen, randval)
   iselres(1) = int((randval * 0.9999999d0) * (trescnt-num_done))

   if (num_done .eq. 0) then
      done_res(1) = iselres(1)
   else
      ! now move the selected residue off of those that have been selected already
      ! THIS ARRAY SHOULD ALREADY BE SORTED (see below)
      do i = 1, num_done
         if (iselres(1) .ge. done_res(i)) iselres(1) = iselres(1) + 1
      end do

      ! If the residue we chose was the larger than our last number, add it to
      ! the end
      if (iselres(1) .gt. done_res(num_done)) &
         done_res(num_done + 1) = iselres(1)
   end if ! (num_done .eq. 0)

   ! now add this to the list of residues that we've titrated already, but
   ! keep this list in order!
   do i = 1, num_done
      if (iselres(1) .lt. done_res(i)) then
         holder_res(i:num_done) = done_res(i:num_done)
         done_res(i) = iselres(1)
         done_res(i+1:num_done+1) = holder_res(i:num_done)
         exit
      end if
   end do

   ! select a different, random state
   call amrand_gen(ce_gen, randval)
   iselstat(1) = int((randval * 0.9999999d0) * &
                     (stateinf(iselres(1))%num_states-1))
   if (iselstat(1) .ge. resstate(iselres(1))) &
      iselstat(1) = iselstat(1)+1
   call cnsteupdatechrg(dcharge, iselres(1), iselstat(1))

end subroutine excnstebeginstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Select random new state for random titrating residue, set charges
subroutine cnsteendstep(dcharge, charge, dvdl, temp, solve)
   use file_io_dat, only : CE_DUMP_UNIT
   use constants, only : KB, FARADAY
   implicit none
#include "extra.h"

   ! Variable descriptions
   !
   ! Passed:
   !  dcharge        : charge array for the proposed state
   !  charge         : charge array for the current state
   !  dvdl           : change in energy between the two states
   !  temp           : system temperature
   !  solve         : solvent Redox potential
   !
   ! Internal:
   !  randval        : random value from random number generator
   !  deltae         : Non-electrostatic factors adjusting dvdl
   !  i              : loop counter
   !  statebase      : the first state of the selected res. in state arrays

   _REAL_, intent(in)    :: temp, solve
   _REAL_, intent(inout) :: dcharge(1:*), charge(1:*)
   _REAL_, intent(inout) :: dvdl

   _REAL_   :: randval, deltae
   integer  :: i, statebase
   integer  :: orig_states(2) ! 2 is max # of residues changed in multi-site MC

   ! This subroutine adjusts the delta E passed from electrostatic calculations
   ! for non-EEL factors, and uses the final result to evaluate a MC transition.
   ! If successful, it will call cnsteupdatechrg charge to modify the charge
   ! array to be the same as the dcharge array and update the resstate to be
   ! the newly updated residue state. If it fails, it will reset the
   ! dcharge array to be the same as the charge array. Then it will modify the
   ! ce_success variable based on whether or not we succeeded.

   ! Back up original states
   orig_states(1) = iselstat(1)
   if (iselres(0) == 2) orig_states(2) = iselstat(2)

   deltae = 0.d0
   do i = 1, iselres(0)
      statebase = stateinf(iselres(i))%first_state
      !  deltae = E(proposed state) - E(current state)
      deltae = deltae + statene(iselstat(i)+statebase) - &
                        statene(resstate(iselres(i))+statebase)
      ! Adjust for the Eo term (Eoref * FARADAY)
      deltae = deltae + (eo_corr(iselstat(i)+statebase) - &
                         eo_corr(resstate(iselres(i))+statebase)) * &
                         FARADAY
      !  correct for Redox potential (delta electrons * E * FARADAY)
      deltae = deltae - (eleccnt(iselstat(i)+statebase) - &
                         eleccnt(resstate(iselres(i))+statebase)) * solve * &
                         FARADAY
   end do
   call amrand_gen(ce_gen, randval)

   dvdl = dvdl - deltae

   if ((dvdl .lt. 0.d0) .or. (randval .le. exp(-dvdl/(KB * temp)))) then
      ce_success = .true.
      do i = 1, iselres(0)
         call cnsteupdatechrg(charge, iselres(i), iselstat(i))
         resstate(iselres(i)) = iselstat(i)
      end do
   else
      do i = 1, iselres(0)
         call cnsteupdatechrg(dcharge, iselres(i), resstate(iselres(i)))
      end do
   end if

   return

end subroutine cnsteendstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Update charges to reflect new state
subroutine cnsteupdatechrg(charge, res, inewstate)
   implicit none

   ! Variables
   !
   ! Passed:
   !  charge         : charge array to be updated
   !  res            : the titr. res. to change the charges of
   !  inewstate      : the charge state to change to
   !
   ! Internal:
   !  i              : loop counter

   _REAL_, intent(inout) :: charge(*)
   integer, intent(in)   :: res, inewstate

   integer :: i

   do i = 0, stateinf(res)%num_atoms - 1
      charge(i+stateinf(res)%first_atom) &
            = chrgdat_e(stateinf(res)%first_charge &
              + inewstate * stateinf(res)%num_atoms + i)
   end do

end subroutine cnsteupdatechrg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write the cerestrt file
subroutine cnstewriterestart(nstep,inchrgdat_e)
   use constants, only : AMBER_ELECTROSTATIC
   use file_io_dat, only : CNSTE_UNIT, cerestrt, owrite, CEOUT_UNIT, ntwr, MAX_FN_LEN
   implicit none

   ! Variable descriptions
   !
   ! Passed:
   !  nstep          : current step number
   !  inchrgdat_e      : chrgdat_e array - sent as a copy to multiply EEL factor
   !
   ! Internal:
   !  chrgdat_e        : chgrdat array adjusted by AMBER_ELECTROSTATIC
   !  iatom          : atom counter
   !  stat           : file status
   !  first          : is this the first pass through?

   integer, intent(in)  :: nstep
   _REAL_, intent(in)   :: inchrgdat_e(0:ATOM_CHRG_C-1)

   _REAL_            :: chrgdat(0:ATOM_CHRG_C-1)
   integer           :: iatom, istart, iend
   character(12)     :: num
   character(len=7)  :: stat
   logical, save     :: first_cerestrt = .true.
   character(256), save :: cerestrt2

   integer :: ntres, ntstates, natchrg, ifind

   namelist /cnstphe_limits/ ntres, ntstates, natchrg

   namelist /cnste/ stateinf, resstate, eleccnt, chrgdat, statene, eo_corr, &
                     trescnt, resname, cefirst_sol, ce_igb, ce_intdiel

   ntres = TITR_RES_C
   ntstates = TITR_STATES_C
   natchrg = ATOM_CHRG_C

   do iatom = 0, ATOM_CHRG_C - 1
      chrgdat(iatom) = inchrgdat_e(iatom) / AMBER_ELECTROSTATIC
   end do

   if (first_cerestrt) then
      if (owrite == 'N') then
         stat = 'NEW'
      else if (owrite == 'O') then
         stat = 'OLD'
      else if (owrite == 'R') then
         stat = 'REPLACE'
      else if (owrite == 'U') then
         stat = 'UNKNOWN'
      end if
      open(unit=CNSTE_UNIT, file=cerestrt, status=stat, form='FORMATTED', &
           delim='APOSTROPHE')
      first_cerestrt = .false.
   else
      open(unit=CNSTE_UNIT, file=cerestrt, status='OLD', form='FORMATTED', &
           delim='APOSTROPHE')
   end if

   write(CNSTE_UNIT, nml=cnstphe_limits)
   write(CNSTE_UNIT, nml=cnste)
   close(CNSTE_UNIT)

   ! flush all ceout data
   call amflsh(CEOUT_UNIT)

   ! Consider whether to save 2ndary restrt:

   if (ntwr .ge. 0) return

   do iend = 1, MAX_FN_LEN
      if (cerestrt(iend:iend) .le. ' ') exit
   end do

   iend = iend - 1

   write(num, '(i12)') nstep

   do istart = 1, 12
      if (num(istart:istart) .ne. ' ') exit
   end do

   write(cerestrt2, '(a,a,a)') cerestrt(1:iend), '_', num(istart:12)

   if (owrite == 'N') then
      stat = 'NEW'
   else if (owrite == 'O') then
      stat = 'OLD'
   else if (owrite == 'R') then
      stat = 'REPLACE'
   else if (owrite == 'U') then
      stat = 'UNKNOWN'
   end if
   open(unit=CNSTE_UNIT, file=cerestrt2, status=stat, form='FORMATTED', &
        delim='APOSTROPHE')

   write(CNSTE_UNIT, nml=cnste)
   close(CNSTE_UNIT)

   return

end subroutine cnstewriterestart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write data to the ceout file
subroutine cnstewrite(rem,remd_types,replica_indexes)
   use file_io_dat, only : CEOUT_UNIT, ntwx
   implicit none
#  include "../include/md.h"

   ! Variables
   !
   ! Passed:
   !  rem            : remd method used
   !
   ! Internal:
   !  full           : do we write a full record or not
   !  i              : loop counter
   !
   ! Common memory:
   !  irespa         : (md.h) step counter for updating slow-varying terms
   !  nstlim         : (md.h) limit of how many steps we run
   !  solve         : (md.h) solvent Redox potential
   !  ntcnste       : (md.h) how often MC oxidation jumps are attempted
   !  t              : (md.h) time

   integer, intent(in) :: rem
   integer, dimension(:), intent(in) :: remd_types, replica_indexes

   logical  :: full
   integer  :: i
   integer  :: nstlim_
   character(len=300) :: txt

   nstlim_ = nstlim
   if (numexchg > 0) &
      nstlim_ = nstlim * numexchg

   if (ntwx .gt. 0) then
      full = (first .or. irespa .eq. nstlim_ .or. mod(irespa,ntwx) .eq. 0)
   else
      full = (first .or. irespa .eq. nstlim_)
   end if

   write(txt, '()')
   if (rem>0) then
      if (rem == 1) then
         write(txt, '(a,f7.2,a)') ' T: ', temp0, ' K'
      else if (rem == 3) then
         write(txt, '(a,i3)') ' H: ', replica_indexes(1)
      else if (rem == 4) then
         write(txt, '(a,f7.3)') ' pH: ', solvph
      else if (rem == 5) then
         write(txt, '(a,f12.7,a)') ' E: ', solve, ' V'
      end if
   else if (rem == -1) then
      do i = 1, size(remd_types)
         if (remd_types(i) == 1) then
            write(txt, '(a,a,f7.2,a)') trim(txt), ' T: ', temp0, ' K'
         else if (remd_types(i) == 3) then
            write(txt, '(a,a,i3)') trim(txt), ' H: ', replica_indexes(i)
         else if (remd_types(i) == 4) then
            write(txt, '(a,a,f7.3)') trim(txt), ' pH: ', solvph
         else if (remd_types(i) == 5) then
            write(txt, '(a,a,f12.7,a)') trim(txt), ' E: ', solve, ' V'
         end if
      end do
   end if
   if (full) then
      write(CEOUT_UNIT, '(a,f12.7,a,f7.2,a)') 'Redox potential: ', solve, ' V Temperature: ', temp0, ' K'
      write(CEOUT_UNIT, '(a,i8)')   'Monte Carlo step size: ', ntcnste
      write(CEOUT_UNIT, '(a,i8)')   'Time step: ', irespa
      write(CEOUT_UNIT, '(a,f14.3)') 'Time: ', t
      do i = 0, trescnt - 1
         write(CEOUT_UNIT, '(a,i4,a,i2,a)') &
            'Residue ', i, ' State: ', resstate(i), trim(txt)
      end do
   else
      do i = 1, iselres(0)
         write(CEOUT_UNIT, '(a,i4,a,i2,a)') 'Residue ', iselres(i), &
            ' State: ', resstate(iselres(i)), trim(txt)
      end do
   end if
   write(CEOUT_UNIT, '()')
   first = .false.

end subroutine cnstewrite

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calculates the total number of "active" electrons
subroutine total_reduction(tot_prot)

   implicit none

! Passed Variables
   integer, intent(out) :: tot_prot

! Local Variables
   integer              :: i ! counter

   tot_prot = 0

   do i = 0, trescnt - 1
      tot_prot = tot_prot + eleccnt(stateinf(i)%first_state + resstate(i))
   end do

end subroutine total_reduction

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Sets up explicit constant Redox potential transition and water relaxation
subroutine cnste_explicitmd(xx,ix,ih,ipairs,x,winv,amass,f,v,vold, &
                             xr,xc,conp,skip,nsp,tma,erstop,qsetup, &
                             do_list_update,rem,remd_types,replica_indexes)
   use amd_mod, only : iamd
   use constantph, only : cph_success
   use file_io_dat, only : CE_DUMP_UNIT, on_cpstep, write_dump_e
   use memory_module, only : natom, l15, l190, l96, l97, l98, l99, ibellygp
   use nbips, only: ips
   use state

   implicit none

#include "box.h"
#include "extra.h"
#include "../include/md.h"
#include "nmr.h"

!  Passed Variables
!
!  XX    : Global real array
!  X     : Global coordinate array
!  WINV  : Inverted masses
!  AMASS : Atomic masses
!  F     : Force array
!  V     : Velocity array
!  VOLD  : Old velocity array (for velocity verlet integrator)
!  XR    : Coords relative to COM of molecule
!  XC    : Restraint reference coords
!  CONP  : Bond parameter for SHAKE
!  SKIP  : Skip array for shake with logical value for each atom
!  NSP   : Submolecule index array (?)
!  TMA   : Submolecular weight array (?)
!  ERSTOP: Should we stop in error (?)
!  QSETUP: Not quite sure...
!  DO_LIST_SETUP : logical that indicates if we should rebuild pairlist (?)

!  Local Variables
!
!  holder_cut     : temporarily store the &cntrl CUT value while we do GB
!  holder_ntb     : temporarily store the &cntrl NTB value while we do GB
!  holder_natom   : temporarily store the true number of atoms
!  holder_ips     : temporarily store the ips value
!  holder_iamd    : temporarily store the iamd value
!  last_solute    : last atom of the solute (first atom of solvent - 1)
!  ce_ener       : record for storing energies
!  iselres        : constant Redox potential residue selection (change prot state)
!  iselstat       : constant Redox potential state selection for chosen residues
!  vtemp          : temporary array holder for the velocities
!  voldtemp       : temporary array holder for the old velocities
!  natom3         : natom * 3
!  selres_holder  : residues that have been chosen to titrate

!  Other used variables from common memory
!
!  igb         : (md.h) -- GB solvent model to use
!  ntb         : (box.h) -- PBC (1 -> V, 2 -> P, 0 -> none)
!  ntp         : (md.h) -- constant pressure control
!  cut         : (box.h) -- value of cutoff
!  master      : (extra.h) -- True if this is PID 0, False otherwise
!  ntt         : (md.h) -- thermostat
!
!  NOTE: All pointer variables into IX and XX are defined in memory.h and
!        described in locmem.f, where they're used during the allocation
!        of the global integer, real, and hollerith arrays.

   ! Passed arguments
   _REAL_  :: xx(*), x(*), winv(*), amass(*), f(*), v(*), vold(*), &
              xr(*), xc(*), conp(*), tma(*)
   logical :: skip(*), erstop, qsetup, do_list_update
   integer :: ix(*), ipairs(*), nsp(*), rem
   integer, dimension(:), intent(in) :: remd_types, replica_indexes
   character (len=4) :: ih(*)

   ! local variable declarations
   integer          :: holder_ntb, holder_natom, j,  &
                       i, k, natom3, selres_holder(1:TITR_RES_C),&
                       holder_ntp, holder_ips, holder_iamd
   type (state_rec) :: ce_ener
   _REAL_, dimension(natom*3) :: vtemp, voldtemp
   _REAL_           :: holder_cut

   ! This subroutine is the main driver for running constant Redox potential MD in explicit
   ! solvent. The first thing it does is call cnstebeginstep to set up the MC.
   ! Then it removes the periodic boundary conditions (PBC), selects a GB model,
   ! then calls force with nstep = 0 which will force oncpstep to be .true., so
   ! we'll get dvdl back. This is then passed to cnsteendstep to evaluate the
   ! transition. If it fails, then we restore the PBC, turn off the GB model,
   ! restore the CUToff, and return to the calling routine. If it's successful,
   ! we still return the above variables, but then we also turn on belly and fix
   ! the protein while we call runmd to relax the solvent for ntrelaxe steps.
   ! Some important things we do are:
   !  o  turn off T-coupling, since ntt=3 doesn't work with belly
   !  o  turn off trajectory/energy printing
   !  o  zero-out the v and vold arrays for the solute
   !  o  turn on belly and assign belly arrays

   ! Local variable assignments
   natom3          = natom * 3
   holder_cut      = cut
   holder_natom    = natom
   holder_ntb      = ntb
   holder_ntp      = ntp
   holder_ips      = ips

   cut   = 9801.0d0    ! 99 * 99
   natom = cefirst_sol - 1
   ntb   = 0
   igb   = ce_igb

   ce_success = .false.

   if (write_dump_e) &
      frac_pop(:,:) = 0

   do j = 1, mccycles_e
      ! Zero-out the holder for the selected residues
      selres_holder(1:trescnt) = 0

      ! Initiate trial moves for oxidation states
      do i = 1, trescnt
         call excnstebeginstep(xx(l190), selres_holder, i-1)

         ! call force to get energies
         call force(xx,ix,ih,ipairs,x,f,ce_ener,ce_ener%vir, &
                 xx(l96),xx(l97),xx(l98),xx(l99),qsetup, &
                 do_list_update,0)

         ! end constant Redox potential stuff
         call cnsteendstep(xx(l190), xx(l15), ce_ener%pot%dvdl, temp0, solve)

         ! We do _not_ want this force call to count toward our nmropt counter,
         ! so we decrement it here if nmropt is on
         if (nmropt /= 0) call nmrdcp()
      end do ! i = 1, trescnt

      if (write_dump_e) then
         do i = 0, trescnt - 1
            frac_pop(i, resstate(i)) = frac_pop(i, resstate(i)) + 1
         end do
      end if

   end do ! j = 1, mccycles_e

   if (write_dump_e) then
      write(CE_DUMP_UNIT, '(29("="),x,"E",f8.3,x,29("="))') solve
      do i = 0, trescnt - 1

         write(CE_DUMP_UNIT, '(("Residue "),i3,9x)', advance='NO') i

         do k = 0, stateinf(i)%num_states - 1
            write(CE_DUMP_UNIT, '(1f8.6,2x)', advance='NO') &
               dble(frac_pop(i, k)) / dble(mccycles_e)
         end do

         ! now add the newline
         write(CE_DUMP_UNIT, '()')

      end do
   end if

   ! Restore PBC variables
   cut   = holder_cut
   natom = holder_natom
   ntb   = holder_ntb
   igb   = 0
   ntp   = holder_ntp
   ips   = holder_ips

   ! set iselres so that every residue is printed out
   iselres(0) = trescnt
   do i = 1, trescnt
      iselres(i) = i - 1
   end do

   ! write out the results to the ceout file
   if (master) call cnstewrite(rem,remd_types,replica_indexes)

   ! if no exchange attempts succeeded, just return and continue dynamics
   if ( .not. ((on_cpstep .and. cph_success) .or. ce_success) ) return

   ! increase the relaxation count
   relaxations_e = relaxations_e + 1
   tot_relax_e = tot_relax_e + 1

   ! Disable AMD for relaxation dynamics
   holder_iamd = iamd
   iamd = 0

   do i = 1, cefirst_sol * 3 - 3
      vtemp(i) = v(i)
      voldtemp(i) = vold(i)
      v(i)     = 0
   end do

   ! call runmd to relax the waters
   if (on_cpstep .and. cph_success) then
     i = MAX(ntrelax,ntrelaxe)
   else
     i = ntrelaxe
   end if
   call relaxmd(xx,ix,ih,ipairs,x,winv,amass,f, &
      v,vold,xr,xc,conp,skip,nsp,tma,erstop, qsetup, &
      i, mobile_atoms, .false.)

   ! restore the original values
   do i = 1, cefirst_sol * 3 - 3
      v(i) = vtemp(i)
      vold(i) = voldtemp(i)
   end do

   ! restore iamd
   iamd = holder_iamd

end subroutine cnste_explicitmd

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Closes all relevant file units and does necessary deallocations
subroutine cnste_finalize

   use file_io_dat, only : CE_DUMP_UNIT, write_dump_e
   implicit none
#  include "extra.h"

   if (allocated(chrgdat_e)) deallocate(chrgdat_e)
   if (allocated(stateinf)) deallocate(stateinf)
   if (allocated(resname)) deallocate(resname)
   if (allocated(statene)) deallocate(statene)
   if (allocated(eo_corr)) deallocate(eo_corr)
   if (allocated(resstate)) deallocate(resstate)
   if (allocated(iselres)) deallocate(iselres)
   if (allocated(eleccnt)) deallocate(eleccnt)

   if (allocated(mobile_atoms)) deallocate(mobile_atoms)
   if (allocated(frac_pop)) deallocate(frac_pop)

   if (write_dump_e) close(CE_DUMP_UNIT)

end subroutine cnste_finalize

end module constante
