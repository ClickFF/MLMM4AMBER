!<compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
!------------------------------------------------------------------------------
! relaxmd: A stripped-down runmd routine for running relaxation dynamics on a
!          given mask.
!
! Units:
!   Runmd operates in kcal/mol units for energies, amu for masses, and
!   Angstroms for distances.  To convert the input time parameters from
!   picoseconds to internal units, multiply by 10.0*sqrt(4.184) = 20.455.
!
! Arguments:
!   xx:        global real array. See locmem.F90 for structure/pointers.
!   ix:        global integer array. See locmem.F90 for structure/pointers.
!   ih:        global Hollerith array holding atom names, residues names, atom
!              types, and other information (see the indexing fully described
!              in locmem.F90)
!   ipairs:    pairlist of nonbonded interactions
!   x:         global position array *
!   winv:      array with inverse masses *
!   amass:     mass array *
!   f:         force array, used to hold old coordinates temporarily, too
!   v:         velocity array
!   vold:      old velocity array, from the previous step
!   xr:        coordinates with respect to COM of molecule
!   xc:        array of reals, matching the size of x itself, used for scratch
!              space in various subroutine calls
!   conp:      bond parameters for SHAKE
!   skip:      logical skip array for SHAKE (and QM/MM too, I think)
!   nsp:       submolecule index array (?)
!   tma:       submolecular weight array (?)
!   erstop:    should we stop in error (?)
!   qsetup:    Flag to activate setup of multiple components, .false. on
!              first call
!
! Additional arguments in relaxmd not common to the main runmd subroutine:
!   mobile_atoms:     bellymask-style array with 1s for moving atoms and 0s for
!                     frozen atoms
!   relax_nstlim:     number of relaxation dynamics steps to run
!   increment_nmropt: flag to signal whether the nmropt counter will increment
!------------------------------------------------------------------------------
subroutine relaxmd(xx, ix, ih, ipairs, x, winv, amass, f, v, vold, xr, xc, &
                   conp, skip, nsp, tma, erstop, qsetup, relax_nstlim, &
                   mobile_atoms, increment_nmropt)

  use bintraj, only: end_binary_frame
  use barostats, only: mcbar_trial
  use constants, only: third, ten_to_minus3
  use fastwt
  use file_io_dat
  use nblist,only: fill_tranvec,volume,oldrecip,ucell
  use qmmm_module, only: qmmm_struct
#ifdef RELAXATION_TRAJ
  use qmmm_module, only: qmmm_nml
#endif /* RELAXATION_TRAJ */
  use random
  use stack
  use state
  use trace

  ! Local variables
  !  factt       : degree-of-freedom correction factor for temperature scaling
  !  nr          : local copy of nrp, number of atoms
  !  nr3         : 3 * nr, used for runtime efficiency
  !
  ! Common memory variables
  !  nrp         : number of atoms, adjusted for LES copies

  implicit none
  integer   ipairs(*), ix(*), relax_nstlim
  integer, intent(in) :: mobile_atoms(*)
  logical, intent(in) :: increment_nmropt
  _REAL_ xx(*)
  character(len=4) ih(*)

#ifdef MPI
#  include "parallel.h"
  include 'mpif.h'
  _REAL_ mpitmp(8) !Use for temporary packing of mpi messages.
  integer ierr
#else
  ! mdloop and REM is always 0 in serial
  integer, parameter :: mdloop = 0, rem = 0
#endif

#include "../include/md.h"
#include "box.h"
#include "nmr.h"
#include "../include/memory.h"
#include "extra.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "ew_mpole.h"
#include "def_time.h"
#include "extra_pts.h"

  _REAL_ sysx, sysy, sysz, sysrange(3,2)

  integer m
   
  logical qspatial

  logical resetvelo
  _REAL_ etot_save,ekpbs
   
  logical do_list_update
  logical skip(*), lout, loutfm, erstop, vlim, onstep
  _REAL_ x(*), winv(*), amass(*), f(*), v(*), vold(*), &
         xr(*), xc(*), conp(*)
  type(state_rec) :: ener
  _REAL_ rmu(3), fac(3), onefac(3)
  _REAL_ tma(*)

  _REAL_ fln,scaltp
  _REAL_ vmax
  _REAL_ aamass, rterm, ekmh, ekph, wfac, rsd
  _REAL_ fit, fiti, fit2
  logical is_langevin  ! Is this a Langevin dynamics simulation
  _REAL_ gammai, c_implic, c_explic, c_ave, sdfac, ekins0
  _REAL_ dtx, dtxinv, dt5, factt, ekin0, ekinp0, dtcp, dttp
  _REAL_ rndf, rndfs, rndfp, boltz2, pconv

  ! Variables and parameters for constant surface tension:
  ! ten_conv converts dyne/cm to bar angstroms
  _REAL_, parameter :: ten_conv = 100.0d0
  _REAL_ :: pres0x
  _REAL_ :: pres0y
  _REAL_ :: pres0z
  _REAL_ :: gamma_ten_int
  _REAL_ :: press_tan_ave

  integer nsp(*)
  integer idumar(4)
#ifdef RELAXATION_TRAJ
  ! DEBUG -- RELAXATION_TRAJ will print the relaxation dynamics
  ! to a trajectory file
  character(kind=1,len=7) :: routine="relaxmd"
  integer l_temp
#endif
  integer i, j, im, i3, nitp, nits
  integer nstep, nrep, nrek, iend, istart3, iend3
  integer nrx, nr, nr3, ntcmt, izero, istart
  logical qsetup
   
  integer nvalid, nvalidi
  _REAL_ eke
  _REAL_ small
  data small/1.0d-7/
  _REAL_, parameter :: pressure_constant = 6.85695d+4

  call trace_enter( 'relaxmd' )
  vlim = vlimit > small
  ntcmt = 0
  izero = 0
  lout = .true.
  loutfm = ioutfm <= 0
  nr = natom
  nr3 = 3*nr
  ekmh = 0.d0
  onstep = .true.
  etot_save = 0.d0
  pres0x = 0.d0
  pres0y = 0.d0
  pres0z = 0.d0
  gamma_ten_int = 0.d0
  dtcp = 0.d0
  dttp = 0.d0

  do_list_update=.false.
#ifdef MPI
  istart = iparpt(mytaskid) + 1
  iend = iparpt(mytaskid+1)
#else
  istart = 1
  iend = nr
#endif
  istart3 = 3*istart -2
  iend3 = 3*iend

  ! If ntwprt is not 0, print only the atoms up to this value
  nrx = nr3
  if (ntwprt > 0) nrx = ntwprt*3
   
  ! Determine system degrees of freedom (for T scaling, reporting).
  ! Call DEGCNT to get the actual number of degrees of freedom for the
  ! solute and solvent. This call returns the correct numbers for belly
  ! simulations and simulations with separate solute/solvent scaling -- dap
  ! "IDUMAR" is dummy array. Used since this routine was also used w/ GIBBS.
  call degcnt(1, nr, mobile_atoms, nsolut, nbonh, nbona, 0, ix(iibh), &
              ix(ijbh), ix(iiba), ix(ijba), idumar, idumar, ntc, idumar, &
              0, 0, 0, idumar, rndfp, rndfs)
   
  ! RNDFP = # degrees of freedom for solute
  ! RNDFS = # degrees of freedom for solvent
  ! RNDF = total number of degrees of freedom.

  ! qtw - substract the number of overlapping noshake QM atoms in noshakemask
  rndfp = rndfp - qmmm_struct%noshake_overlap

  ! Modify RNDFP to reflect NDFMIN (set in mdread) and num_noshake
  rndfp = rndfp - ndfmin + num_noshake
  rndf = rndfp + rndfs

  ! Correct the degree of freedom count for extra points; this
  ! ends to work of computing degrees of freedom.
  call fix_degree_count(rndf)
   
  ! Begin unit conversion.  pconv eventually becomes a factor to convert
  ! pressure in kcal/mole-A^3 to bar.
  boltz2 = 8.31441d-3 * 0.5d0
  pconv = 1.6604345d+04
  boltz2 = boltz2/4.184d0
  dtx = dt*20.455d+00
  dtxinv = 1.0d0 / dtx
  dt5 = dtx * 0.5d0
  pconv = pconv*4.184d0
   
  ! The values in the array fac() are degrees of freedom * kboltz / 2.
  ! Multiply by T to get expected kinetic energy.
  ! fac(1) applies to the whole system.
  fac(1) = boltz2*rndf
  fac(2) = boltz2*rndfp

  if (rndfp < 0.1d0) then
    fac(2) = 1.d-6
  end if

  fac(3) = boltz2*rndfs
  if (rndfs < 0.1d0) then
    fac(3) = 1.d-6
  end if
  onefac(1) = 1.0d0/fac(1)
  onefac(2) = 1.0d0/fac(2)
  onefac(3) = 1.0d0/fac(3)
  factt = rndf / (rndf + ndfmin)
   
  ! These are "desired" kinetic energies based on the
  ! number of degrees freedom and target temperature.
  ! They will be used for calculating the velocity scaling factor.
  ekinp0 = fac(2)*temp0
  ekins0 = fac(3)*temp0
  ekin0  = fac(1)*temp0

  ! Langevin dynamics setup:
  is_langevin = (gamma_ln > 0.0d0)
  gammai = gamma_ln / 20.455d0
  c_implic = 1.d0 / (1.d0 + gammai*dt5)
  c_explic = 1.d0 - gammai*dt5
  c_ave = 1.d0 + gammai*dt5
  sdfac = sqrt(4.d0 * gammai * boltz2 * temp0 / dtx)
  if (is_langevin .and. ifbox == 0) then
    call get_position(nr, x, sysx, sysy, sysz, sysrange, 0)
  end if
  if (ntt == 1) then
    dttp = dt / tautp
  end if
  if (ntp > 0) then
    dtcp = comp * 1.0d-06 * dt / taup
  end if

  ! Constant surface tension setup:
  if (csurften > 0) then

    ! Set pres0 in direction of surface tension.
    ! The reference pressure is held constant in on direction dependent
    ! on what the surface tension direction is set to.
    if (csurften .eq. 1) then

      ! pres0 in the x direction
      pres0x = pres0
    else if (csurften .eq. 2) then

      ! pres0 in the y direction
      pres0y = pres0

    else

      ! pres0 in the z direction
      pres0z = pres0
    end if

    ! Multiply surface tension by the number of interfaces
    gamma_ten_int = dble(ninterface) * gamma_ten
  end if

  nrek = 4
  nrep = 15
  
  nvalid = 0
  nvalidi = 0
  nstep = 0
  fit = 0.d0
  fiti = 0.d0
  fit2 = 0.d0
  ener = null_state_rec
  ener%kin%pres_scale_solt = 1.d0
  ener%kin%pres_scale_solv = 1.d0
  ener%box(1:3) = box(1:3)
  ener%cmt(1:4) = 0.d0
  nitp = 0
  nits = 0
  ekmh = 0.0d0

  i3 = 0
  do j = 1, nrp
    aamass = amass(j)
    do m = 1, 3
      i3 = i3 + 1
      rterm = v(i3) * v(i3) * aamass
      ekmh = ekmh + rterm
    end do
  end do

  do im = 1, iscale
    ekmh = ekmh + scalm * v(nr3+im) * v(nr3+im)
  end do
  ekmh = ekmh * 0.5d0
  do i = 1, nr3+iscale
    vold(i) = v(i)
  end do
   
  ! The main loop for performing dynamics steps: at this point, the
  ! coordinates are a half-step "ahead" of the velocities; the variable
  ! EKMH holds the kinetic energy at these "-1/2" velocities, which are
  ! stored in the array vold.
  260 continue

  ! Step 1a: do some setup for pressure calculations
  if (ntp > 0) then
    ener%cmt(1:3) = 0.d0
    xr(1:nr3) = x(1:nr3)
      
    ! Calculate, for each molecule, the center of mass, kinetic energy
    ! of the center of mass, and the coordinates of the molecule
    ! relative to that center of mass.
    call timer_start(TIME_EKCMR)
    call ekcmr(nspm, nsp, tma, ener%cmt, xr, v, amass, istart, iend)
#ifdef MPI
    call trace_mpi('mpi_allreduce', 3, 'MPI_DOUBLE_PRECISION', mpi_sum)
#  ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE, ener%cmt, 3, MPI_DOUBLE_PRECISION, &
                       mpi_sum, commsander, ierr)
#  else
    call mpi_allreduce(ener%cmt, mpitmp, 3, MPI_DOUBLE_PRECISION, mpi_sum, &
                       commsander, ierr)
    ener%cmt(1:3) = mpitmp(1:3)
#  endif
#endif
    call timer_stop(TIME_EKCMR)
  end if

  ! If we're using the MC barostat, go ahead and do the trial move now
  if (ntp > 0 .and. barostat == 2 .and. mod(nstep+1, mcbarint) == 0) then
    call mcbar_trial(xx, ix, ih, ipairs, x, xc, f, ener%vir, xx(l96), &
                     xx(l97), xx(l98), xx(l99), qsetup, do_list_update, &
                     nstep, nsp, amass)
  end if

  ! Step 1b: Get the forces for the current coordinates
  iprint = 0
  if (nstep == 0 .or. nstep+1 == relax_nstlim) then
    iprint = 1
  end if
  call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
             xx(l98), xx(l99), qsetup, do_list_update, nstep)

  ! If we don't want to increment the NMROPT counter, decrement it here.
  if (.not. increment_nmropt) then
    call nmrdcp
  end if

  ! Reset quantities depending on TEMP0 and TAUTP (which may have been
  ! changed by MODWT during FORCE call).
  ekinp0 = fac(2)*temp0
  ekins0 = fac(3)*temp0
  ekin0 = fac(1)*temp0
  if (ntt == 1) then
    dttp = dt / tautp
  end if
  if (ntp > 0) then
    ener%volume = volume
    ener%density = tmass / (0.602204d0*volume)
    ener%cmt(4) = 0.d0
    ener%vir(4) = 0.d0
    ener%pres(4) = 0.d0
    do m = 1, 3
      ener%cmt(m)  = ener%cmt(m) * 0.5d0
      ener%cmt(4)  = ener%cmt(4) + ener%cmt(m)
      ener%vir(4)  = ener%vir(4) + ener%vir(m)
      ener%pres(m) = (pconv+pconv) * (ener%cmt(m) - ener%vir(m)) / volume
      ener%pres(4) = ener%pres(4) + ener%pres(m)
    end do
    ener%pres(4) = ener%pres(4)/3.d0

    ! Constant surface tension output
    if (csurften > 0) then
      if (csurften == 1) then

        ! Surface tension in the x direction
        ener%surface_ten = box(1) * &
          (ener%pres(1) - 0.5d0 * (ener%pres(2) + ener%pres(3))) / &
          (ninterface * ten_conv)
      else if (csurften .eq. 2) then

        ! Surface tension in the y direction
        ener%surface_ten = box(2) * &
          (ener%pres(2) - 0.5d0 * (ener%pres(1) + ener%pres(3))) / &
          (ninterface * ten_conv)
      else
 
        ! Surface tension in the z direction
        ener%surface_ten = box(3) * &
          (ener%pres(3) - 0.5d0 * (ener%pres(1) + ener%pres(2))) / &
          (ninterface * ten_conv)
      end if
    end if
  end if

  ! Step 1c: do randomization of velocities, if needed.
  ! Assign new random velocities every Vrand steps for Andersen thermostating
  resetvelo = .false.
  if (vrand .ne. 0 .and. ntt == 2) then
    if (mod((nstep+1), vrand) == 0) then
      resetvelo = .true.
    end if
  end if
  if (resetvelo) then
    ! DAN ROE: Why are only the masters doing this?  Even if the velocities 
    ! are broadcast to the child processes, the wont the different # of random
    ! calls put the randomg num generators out of sync, or do we not care?
    if (master) then
      call setvel(nr, v, winv, temp0*factt, iscale, scalm)
    end if

#ifdef MPI
    call trace_mpi('mpi_bcast', 3*natom, 'MPI_DOUBLE_PRECISION', 0)
    call mpi_bcast(v, 3*natom, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
#endif /* MPI */
    ! At this point in the code, the velocities lag the positions by half
    ! a timestep.  If we intend for the velocities to be drawn from a 
    ! Maxwell distribution at the timepoint where the positions and 
    ! velocities are synchronized, we have to correct these newly redrawn
    ! velocities by backing them up half a step using the current force.
    ! Note that this fix only works for Newtonian dynamics.
    if (gammai == 0.d0) then
      i3 = 3*(istart - 1)
      do j = istart, iend
        wfac = winv(j) * dt5
        v(i3+1) = v(i3+1) - f(i3+1)*wfac
        v(i3+2) = v(i3+2) - f(i3+2)*wfac
        v(i3+3) = v(i3+3) - f(i3+3)*wfac
        i3 = i3+3
      end do
    end if
  end if
  ! End contingency for resetting velocities based on Andersen thermocoupling

  

  ! Step 2: Do the velocity update.
  ! Step 2a: apply quenched MD if needed.
  ! This is useful in Nudged Elastic Band simulations.
  call timer_start(TIME_VERLET)
  if (vv == 1) then
    call quench(f, v)
  end if
  if (gammai == 0.d0) then

    ! Newtonian dynamics
    i3 = 3*(istart - 1)
    do j = istart, iend
      wfac = winv(j) * dtx
      v(i3+1) = v(i3+1) + f(i3+1)*wfac
      v(i3+2) = v(i3+2) + f(i3+2)*wfac
      v(i3+3) = v(i3+3) + f(i3+3)*wfac
      i3 = i3+3
    end do
  else

    ! gamma_ln .ne. 0, which implies Langevin dynamics (ntt == 3)
    ! (see mdread.f).  The simple model for Langevin dynamics, is
    ! taken from Loncharich, Brooks and Pastor, Biopolymers
    ! 32:523-535 (1992), Eq. 11. (Note that the first term on the
    ! rhs of Eq. 11b should not be there.).  Update the Langevin
    ! parameters, since temp0 might have changed:
    sdfac = sqrt(4.d0 * gammai * boltz2 * temp0 / dtx)
    i3 = 3*(istart - 1)
    if (no_ntt3_sync == 1) then

      ! We don't worry about synchronizing the random number stream
      ! across processors.
      do j = istart, iend
        wfac = winv(j) * dtx
        aamass = amass(j)         
        rsd = sdfac*sqrt(aamass)
        call gauss( 0.d0, rsd, fln )
        v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
        call gauss( 0.d0, rsd, fln )
        v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
        call gauss( 0.d0, rsd, fln )
        v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
        i3 = i3+3
      end do
    else
      do j = 1, nr
        if (j < istart .or. j > iend) then

          ! In order to generate the same sequence of pseudorandom
          ! numbers that you would using a single processor you have
          ! to go through the atoms in order.  The unused results
          ! are thrown away
          call gauss(0.d0, 1.d0, fln)
          call gauss(0.d0, 1.d0, fln)
          call gauss(0.d0, 1.d0, fln)
          cycle
        end if
        wfac = winv(j) * dtx
        aamass = amass(j)         
        rsd = sdfac*sqrt(aamass)
        call gauss(0.d0, rsd, fln)
        v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
        call gauss(0.d0, rsd, fln)
        v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
        call gauss(0.d0, rsd, fln)
        v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
        i3 = i3+3
      end do
    end if
    ! End branch for Langevin dynamics sychronization (based on the
    ! global integer no_ntt3_sync coming in from ../include/md.h)

  end if
  ! End the branch based on whether Langevin dynamics is in effect

  if (vlim) then
    vmax = 0.0d0
    do i=istart3,iend3
      vmax = max(vmax, abs(v(i)))
      v(i) = sign(min(abs(v(i)), vlimit), v(i))
    end do

    ! Only violations on the master node are actually reported
    ! to avoid both MPI communication and non-master writes.
    if (vmax > vlimit) then
      if (master) then
        write(6, '(a,i6,a,f10.4)') 'vlimit exceeded for step ', nstep, &
              '; vmax = ',vmax
      end if
    end if
  end if
   
  do im = 1, iscale
    v(nr3+im) = (v(nr3+im) + f(nr3+im)*dtx/scalm)
  end do

  ! Step 3: update the positions, putting the
  ! "old" positions into the array f
  i = istart - 1
  do i3 = istart3, iend3, 3
    f(i3) = x(i3)
    f(i3+1) = x(i3+1)
    f(i3+2) = x(i3+2)
    if (mobile_atoms(i) == 1) then
      x(i3) = x(i3) + v(i3)*dtx
      x(i3+1) = x(i3+1) + v(i3+1)*dtx
      x(i3+2) = x(i3+2) + v(i3+2)*dtx
    end if
    i = i + 1
  end do
  do i = 1,iscale
    f(nr3+i) = x(nr3+i)
    x(nr3+i) = x(nr3+i) + v(nr3+i)*dtx
  end do
  call timer_stop(TIME_VERLET)
  ! End of this segment of Velocity Verlet work

  if (ntc .ne. 1) then

    ! Step 4a: if shake is being used, update the new positions to fix
    !          the bond lengths.
    call timer_start(TIME_SHAKE)
    qspatial = .false.
    call shake(nrp, nbonh, nbona, 0, ix(iibh), ix(ijbh), ix(ibellygp), &
               winv, conp, skip, f, x, nitp, .false., ix(iifstwt), &
               ix(noshake), qspatial)
    call quick3(f, x, ix(iifstwr), natom, nres, ix(i02))
    if (nitp == 0) then
      erstop = .true.
      goto 480
    end if

    ! Step 4b: fix the velocities, calculate KE, and re-estimate
    !          the velocities from differences in positions.
    v(istart3:iend3) = (x(istart3:iend3)-f(istart3:iend3)) * dtxinv
    call timer_stop(TIME_SHAKE)
  end if
  call timer_start(TIME_VERLET)
  if (ntt == 1 .or. onstep) then

    ! Step 4c: get the KE, either for averaging or for Berendsen:
    eke = 0.d0
    ekph = 0.d0
    ekpbs = 0.d0
    if (gammai == 0.0d0) then
      i3 = 3*(istart - 1)
      do j = istart, iend
        aamass = amass(j)
        do m = 1, 3
          i3 = i3 + 1
          eke = eke + aamass*0.25d0*(v(i3) + vold(i3))**2

          ! Try pseudo KE from Eq. 4.7b of Pastor, Brooks & Szabo,
          ! Mol. Phys. 65, 1409-1419 (1988):
          ekpbs = ekpbs + aamass*v(i3)*vold(i3)
          ekph = ekph + aamass*v(i3)**2
        end do
      end do
    else
      i3 = 3*(istart - 1)
      do j = istart, iend
        aamass = amass(j)
        do m = 1, 3
          i3 = i3 + 1
          eke = eke + aamass*0.25d0*c_ave*(v(i3) + vold(i3))**2
        end do
      end do
    end if
    ! End the branch based on whether Langevin dynamics is in effect

#ifdef MPI
    ! Sum up the partial kinetic energies:
    if (numtasks > 1) then
      call trace_mpi('mpi_allreduce', 1, 'MPI_DOUBLE_PRECISION', mpi_sum)
      mpitmp(1) = eke
      mpitmp(2) = ekph
      mpitmp(3) = ekpbs
#  ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE,mpitmp, 3, MPI_DOUBLE_PRECISION, &
                       mpi_sum, commsander, ierr)
    eke = mpitmp(1)
    ekph = mpitmp(2)
    ekpbs = mpitmp(3)
#  else
      call mpi_allreduce(mpitmp,mpitmp(4), 3, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commsander, ierr)
      eke = mpitmp(4)
      ekph = mpitmp(5)
      ekpbs = mpitmp(6)
#  endif /* USE_MPI_IN_PLACE */
    end if
#endif /* MPI */
      
    ! All processors handle the "extra" variables:
    do im = 1, iscale
      eke = eke + scalm*0.25d0*(v(nr3+im) + vold(nr3+im))**2
      ekpbs = ekpbs + scalm*v(nr3+im)*vold(nr3+im)
      ekph = ekph + scalm*v(nr3+im)**2
    end do
    eke = eke * 0.5d0
    ekph = ekph * 0.5d0
    ekpbs = ekpbs * 0.5d0
    if (ntt == 1) then
         
      ! The following is from T.E. Cheatham, III and B.R. Brooks,
      ! Theor. Chem. Acc. 99:279, 1998.
      scaltp = sqrt(1.d0 + 2.d0*dttp*(ekin0 - eke) / (ekmh + ekph))
      do j = istart, iend
        i3 = (j - 1)*3 + 1
        v(i3) = v(i3) * scaltp
        v(i3+1) = v(i3+1) * scaltp
        v(i3+2) = v(i3+2) * scaltp
      end do
      do im = 1, iscale
        v(nr3+im) = v(nr3+im) * scaltp
      end do
    end if
    ! End contingency for Berendsen thermocoupling
  end if
  ! End of step 4c: a contingency for kinetic energy computation when
  ! we are either on a reportable step or doing Berendsen thermocoupling

  ! Step 5: several tasks related to dumping of trajectory information
  !
  ! Determine if trajectory, velocity, or restart writing is
  ! imminent, or if the center of mass motion will be removed.
  ! These requires distribution of velocities or dipoles to
  ! all processes by the subroutine xdist in parallel runs.
  !
  ! Modified so that when running REMD, writing can occur less often
  ! than exchanges (e.g. ntwx > nstlim).  Two new variables, total_nstep
  ! and total_nstlim were added.  For non-REMD runs,
  ! total_nstep = nstep + 1 and total_nstlim = nstlim as before.
  !
  ! For REMD runs, total_nstep = (mdloop-1)*nstlim + nstep + 1, where
  ! mdloop is the current exchange - this is the current
  ! replica exchange MD step.  total_nstlim = numexchg*nstlim, which
  ! is the maximum number of REMD steps.
#ifdef MPI
  ! Distribute the coordinates, dipoles, and velocities as necessary
  call timer_barrier(commsander)
  call timer_stop_start(TIME_VERLET, TIME_DISTCRD)
  if (numtasks > 1) then
    call xdist(x, xx(lfrctmp), natom)
  end if
  call timer_stop(TIME_DISTCRD)
#endif  /* MPI */

  ! Fix lone pair positions
  if (numextra > 0) then
    call local_to_global(x, xx, ix)
  end if

#ifdef MPI
  ! Additional Velocity Verlet work
  call timer_start(TIME_VERLET)
#endif  /* MPI */
   
  ! Step 6: zero COM velocity if requested.  This is used for preventing
  ! the "block of ice flying thru space" phenomenon (caused by Berendsen
  ! thermocoupling, or by coarse Ewald approximations).  It also
  ! prevents accumulations of rotational momentum in vacuum simulations.
  vold(istart3:iend3) = v(istart3:iend3)
  do im = 1, iscale
    vold(nr3+im) = v(nr3+im)
  end do

  ! Step 7: scale coordinates if constant pressure run:
  if (ntp == 1) then
    rmu(1) = (1.d0 - dtcp*(pres0 - ener%pres(4)))**third
    rmu(2) = rmu(1)
    rmu(3) = rmu(1)
  else if (ntp == 2) then
    if (csurften > 0) then

      ! Constant surface tension adjusts the tangential pressures
      ! See Zhang, Feller, Brooks, Pastor. J. Chem. Phys. 1995
      if (csurften == 1) then

        ! For surface tension in the x direction
        pres0y = pres0x - gamma_ten_int * ten_conv / box(1)
        pres0z = pres0y
      else if (csurften == 2) then

        ! For surface tension in the y direction
        pres0x = pres0y - gamma_ten_int * ten_conv / box(2)
        pres0z = pres0x
      else

        ! For surface tension in the z !direction
        pres0x = pres0z - gamma_ten_int * ten_conv / box(3)
        pres0y = pres0x
      end if
      rmu(1) = (1.d0 - dtcp * (pres0x - ener%pres(1)))**third
      rmu(2) = (1.d0 - dtcp * (pres0y - ener%pres(2)))**third
      rmu(3) = (1.d0 - dtcp * (pres0z - ener%pres(3)))**third
    else
      rmu(1) = (1.d0-dtcp*(pres0-ener%pres(1)))**third
      rmu(2) = (1.d0-dtcp*(pres0-ener%pres(2)))**third
      rmu(3) = (1.d0-dtcp*(pres0-ener%pres(3)))**third
    end if
  else

    ! This means ntp = 3, semi-isotropic pressure coupling.
    ! Currently this only works with csurften > 0, constant
    ! surface tension.  Semi-isotropic pressure scaling in
    ! any direction with no constant surface tension has yet
    ! to be implemented.
    if (csurften > 0) then
      if (csurften == 1) then

        ! For surface tension in the x direction
        pres0y = pres0x - gamma_ten_int * ten_conv / box(1)
        pres0z = pres0y
        press_tan_ave = (ener%pres(2) + ener%pres(3))/2
        rmu(1) = (1.d0 - dtcp * (pres0x - ener%pres(1)))**third
        rmu(2) = (1.d0 - dtcp * (pres0y - press_tan_ave))**third
        rmu(3) = (1.d0 - dtcp * (pres0z - press_tan_ave))**third
      else if (csurften == 2) then

        ! For surface tension in the y direction
        pres0x = pres0y - gamma_ten_int * ten_conv / box(2)
        pres0z = pres0x
        press_tan_ave = (ener%pres(1) + ener%pres(3))/2
        rmu(1) = (1.d0 - dtcp * (pres0x - press_tan_ave))**third
        rmu(2) = (1.d0 - dtcp * (pres0y - ener%pres(2)))**third
        rmu(3) = (1.d0 - dtcp * (pres0z - press_tan_ave))**third
      else

        ! For surface tension in the z direction
        pres0x = pres0z - gamma_ten_int * ten_conv / box(3)
        pres0y = pres0x
        press_tan_ave = (ener%pres(1) + ener%pres(2))/2
        rmu(1) = (1.d0 - dtcp * (pres0x - press_tan_ave))**third
        rmu(2) = (1.d0 - dtcp * (pres0y - press_tan_ave))**third
        rmu(3) = (1.d0 - dtcp * (pres0z - ener%pres(3)))**third
      end if
    end if
    ! End contingency for semi-isotropic or anisotropic pressure scaling
  end if
  if (ntp > 0) then
    box(1:3) = box(1:3)*rmu(1:3)
    ener%box(1:3) = box(1:3)
         
    ! WARNING!!   This is not correct for non-orthogonal boxes if
    ! NTP > 1 (i.e. non-isotropic scaling).  Currently, general cell
    ! updates which allow cell angles to change are not implemented.
    ! The viral tensor computed for ewald is the general Nose Klein,
    ! however the cell response needs a more general treatment.
    call redo_ucell(rmu)

    ! Keep tranvec up to date, rather than recomputing each MD step.
    ! Tranvec is dependent on only ucell
    call fill_tranvec()
    call ew_pscale(natom,x,amass,nspm,nsp,npscal)
    if (ntr > 0 .and. nrc > 0) then
      call ew_pscale(natom, xc, amass, nspm, nsp, npscal)
    end if
  endif

  ! Pastor, Brooks, Szabo conserved quantity
  ! for harmonic oscillator: Eq. 4.7b of Mol.
  ! Phys. 65:1409-1419, 1988
  ener%kin%solv = ekpbs + ener%pot%tot  
  ener%kin%solt = eke
  ener%kin%tot  = ener%kin%solt
  if (ntt == 1 .and. onstep) then
    ekmh = max(ekph, fac(1)*10.d0)
  end if

  ! If velocities were reset, the KE is not accurate; fudge it
  ! here to keep the same total energy as on the previous step.
  ! Note that this only affects printout and averages for Etot
  ! and KE -- it has no effect on the trajectory, or on any averages
  ! of potential energy terms.  Total energy is the sum of KE and PE.
  if (resetvelo) then
    ener%kin%tot = etot_save - ener%pot%tot
  end if
  ener%tot  = ener%kin%tot + ener%pot%tot
  etot_save = ener%tot
 
  ! Step 8:  update the step counter and the integration time:
  nstep = nstep + 1

  ! Full energies are only calculated every nrespa steps.
  ! nvalid is the number of steps where all energies are calculated.
  ntnb = 0
  if (mod(nstep,nsnb) == 0) then
    ntnb = 1
  end if

#ifdef RELAXATION_TRAJ
  ! DEBUG code -- this will print out every frame of the
  ! relaxation dynamics to the trajectory if uncommented
  if (master) then
   
    ! Coordinate archiving will always occur in this case
    if (iwrap == 0) then
      call corpac(x, 1, nrx, MDCRD_UNIT, loutfm)
      if (ntb > 0) then
        call corpac(box, 1, 3, MDCRD_UNIT, loutfm)
      end if
    else if (iwrap == 1) then
      call get_stack(l_temp, nr3, routine)
      if (.not. rstack_ok) then
        deallocate(r_stack)
        allocate(r_stack(1:lastrst), stat=alloc_ier)
        call reassign_rstack(routine)
      endif
      REQUIRE(rstack_ok)
      do m = 1, nr3
        r_stack(l_temp+m-1) = x(m)
      end do
      call wrap_molecules(nspm, nsp, r_stack(l_temp))
      if (ifbox == 2) then
        call wrap_to(nspm, nsp, r_stack(l_temp), box)
      end if
      call corpac(r_stack(l_temp), 1, nrx, MDCRD_UNIT, loutfm)
      call corpac(box, 1, 3, MDCRD_UNIT, loutfm)
      call free_stack(l_temp, routine)
    end if
    ! End branch based on molecule wrapping

    ! If using variable QM solvent, try to write a new pdb file
    ! with the QM coordinates for this step. This is done here
    ! to keep the PDB file in sync with the mdcrd file, which
    ! makes it easier to check later.
    if (qmmm_nml%vsolv > 0 .and. qmmm_nml%verbosity == 0) then
      call qm_print_coords(nstep,.false.)
    end if
  end if
  ! End output writing operations on the master process

  if (ioutfm > 0) then
    if (.true.) call end_binary_frame(MDCRD_UNIT)
  end if
#endif /* RELAXATION_TRAJ */

  ! Major cycle back to new step unless we have reached our limit:
  call trace_integer( 'end of step', nstep )
  call trace_output_mpi_tally( )
  call timer_stop(TIME_VERLET)
  if (nstep < relax_nstlim) then
    goto 260
  end if
  480 continue

end subroutine relaxmd

