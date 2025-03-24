!<compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
#include "nfe-config.h"

!  In vim, :set foldmethod=marker, then zo,zc, to open/close single folds
!          zO,ZC to open close folds recursively
!          zM closes all folds; zR opens all folds

!------------------------------------------------------------------------------
! runmd: main driver routine for molecular dynamics.  This is over 4000 lines
!        long and incorporates virtually every method that sander offers.
!        Prior to calling this routine, the sander subroutine itself will take
!        user input, read system specifications from coordinates and topology
!        files, and lay out a series of arrays to hold the information that
!        runmd then operates with.
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
!------------------------------------------------------------------------------
subroutine runmd(xx, ix, ih, ipairs, x, winv, amass, f, v, vold, xr, xc, &
                 conp, skip, nsp, tma, erstop, qsetup)

!------------------------------------------------------------------------------
! modules used:  {{{

  use state

#if !defined(DISABLE_NFE) && defined(NFE_ENABLE_BBMD)
  use nfe_sander_hooks, only : nfe_on_mdstep => on_mdstep
#endif
  use nfe_sander_proxy, only : infe, nfe_real_mdstep, nfe_prt

  use commandline_module, only: cpein_specified
  use molecule, only: n_iwrap_mask_atoms, iwrap_mask_atoms
  use cmd_vars, only: activate, file_pos_cmd, file_vel_cmd, nstep_cmd, &
                      t_cmd, eq_cmd, restart_cmd, etot_cmd, eke_cmd, temp_cmd
  use pimd_vars, only: ipimd, natomCL, bnd_vir, Eimp_virial, equal_part, &
                       Epot_deriv, tau_vol, Epot_spring, NMPIMD, CMD, &
                       cartpos, cartvel
  use neb_vars, only: ineb, neb_nbead
  use lscivr_vars, only: ilscivr, ndof_lsc, natom_lsc, mass_lsc, v2_lsc, &
                         ilsc, x_lsc, f_lsc, dx_lsc
  use md_scheme, only: ischeme, thermostat_step
  use nose_hoover_module, only: thermo_lnv, x_lnv, x_lnv_old, v_lnv, &
                                f_lnv_p, f_lnv_v, c2_lnv, mass_lnv, &
                                Thermostat_init
#ifdef RISMSANDER
  use sander_rism_interface, only: rismprm, RISM_NONE, RISM_FULL, &
                                   RISM_INTERP, rism_calc_type, &
                                   rism_solvdist_thermo_calc, mylcm
#endif /* RISMSANDER */

  use full_pimd_vars, only: totener,totenert,totenert2,mybeadid,xall
  use qmmm_module, only: qmmm_nml,qmmm_struct, qmmm_mpi, qm2_struct
  use ani_module, only: ani_struct,ani_nml
  use call_ani, only: ani_get_mass_center
#ifdef MPI
  use qmmm_module, only: qmmm_vsolv
#endif /* MPI */

  use file_io_dat
  use constants, only: third, ten_to_minus3
  use trace
  use stack

#ifdef MPI
  use decomp, only: decpr, jgroup, indx, irespw, nrs, collect_dec, &
                    checkdec, printdec
#else
  use decomp, only: decpr, jgroup, indx, irespw, checkdec, printdec
#endif /* MPI */

  use fastwt
  use bintraj, only: end_binary_frame
  use nblist,only: fill_tranvec,volume,oldrecip,ucell,olducell,bcshift

  use nose_hoover_module, only: Thermostat_switch, Thermostat_integrate_1, &
                                Thermostat_integrate_2, &
                                Thermostat_hamiltonian, &
                                Adaptive_Thermostat_integrate, &
                                Adaptive_Thermostat_hamiltonian, &
                                file_nhc, nchain, thermo, nthermo, Econserved
  use sgld, only: isgld, sgenergy, sgldw, sgmdw
  use resamplekin_mod, only: resamplekin

#ifdef LES
  ! Self-Guided molecular/Langevin Dynamics (SGLD)
  use les_data, only: cnum, temp0les, tempsules, sdfacles, scaltles, ekeles, &
                      ekinles0, ekmhles, ekphles, rndfles
#else
  use sgld, only: isgsta,isgend,sg_fix_degree_count
  use pimd_vars, only: nbead, itimass
#endif

#ifdef MPI
  use evb_parm,  only: evb_dyn, nbias
  use evb_data,  only: evb_frc, evb_vel0, evb_bias, evb_nrg, evb_nrg_ave, &
                       evb_nrg_rms, evb_nrg_tmp, evb_nrg_old, evb_nrg_tmp2, &
                       evb_nrg_old2
  use remd, only: rem, mdloop, remd_ekmh, repnum, stagid, my_remd_data, &
                  hybrid_remd_ene, next_rem_method, remd_types, replica_indexes, &
                  my_pressure, my_volume
#  ifdef LES
  use evb_pimd,  only: evb_pimd_dealloc
  use miller,    only: i_qi
  use pimd_vars, only : real_mass
#  endif /* LES */
  use softcore, only: ifsc, sc_dvdl, sc_tot_dvdl, sc_tot_dvdl_partner, &
                      sc_dvdl_ee, sc_tot_dvdl_ee, sc_tot_dvdl_partner_ee, &
                      extra_atoms, mix_temp_scaling, sc_pscale, &
                      adj_dvdl_stat, sc_mix_velocities, sc_nomix_frc, &
                      sc_sync_x, sc_print_energies, calc_softcore_ekin, &
                      sc_belly_pscale, bellyf_softcore, &
                      sc_ener, sc_ener_ave, sc_ener_rms, sc_lngdyn, &
                      sc_ener_tmp, sc_ener_tmp2, sc_ener_old, sc_ener_old2, &
                      sc_mix_position, sc_print_dvdl_values, &
                      sc_degrees_o_freedom, dynlmb, sc_change_clambda, &
                      ti_ene_cnt, sc_compare
  use mbar, only : ifmbar, bar_intervall, calc_mbar_energies, &
                   bar_collect_cont, do_mbar
#endif /* MPI */

  use amoeba_mdin, only: iamoeba
  use amoeba_runmd, only: AM_RUNMD_scale_cell
  use constantph, only: cnstphinit, cnstphwrite, cnstphupdatepairs, &
                        cnstphbeginstep, cnstphendstep, chrgdat, &
                        cnstph_explicitmd, cnstphwriterestart, cphfirst_sol
  use constante, only: cnsteinit, cnstewrite, &
                       cnstebeginstep, cnsteendstep, chrgdat_e, &
                       cnste_explicitmd, cnstewriterestart, cefirst_sol
  use emap, only:temap,emap_move
  use barostats, only : mcbar_trial, mcbar_summary

#ifdef EMIL
  use emil_mod, only : emil_do_calc, emil_init, emil_step
#endif /* EMIL */

  use memory_module, only: mass
  use random

#ifdef MPI
  use sgld, only : trxsgld
#endif

  ! Andreas Goetz's adaptive QM/MM
  use qmmm_adaptive_module, only: adaptive_qmmm
  use crg_reloc, only: ifcr, crprintcharges, cr_print_charge
  use abfqmmm_module, only: abfqmmm_param, abfqmmm_combine_forces

  ! Accelerated Molecular Dynamics (aMD)
  use amd_mod

  ! scaledMD
  use scaledMD_mod

  ! SEBOMD: Semi-Empirical Born-Oppenheimer MD
  use sebomd_module, only: sebomd_obj, sebomd_gradient_write, &
                           sebomd_hessian_compute

  ! SINR (Stochastic Isokinetic Nose-Hoover RESPA integrator)
  use sinr_t

  ! }}}

  ! Local variables  {{{
  !  factt       : degree-of-freedom correction factor for temperature scaling
  !  nr          : local copy of nrp, number of atoms
  !  nr3         : 3 * nr, used for runtime efficiency
  !
  ! Common memory variables
  !  nrp         : number of atoms, adjusted for LES copies

!------------------------------------------------------------------------------
! local variables:

  implicit none
  integer ani_i
  double precision ani_x, ani_y, ani_z
  character(kind=1,len=5) :: routine="runmd"
  integer   ipairs(*), ix(*)
  _REAL_ xx(*)
  character(len=4) ih(*)
#ifdef MPI
#  include "parallel.h"
  include 'mpif.h'
  _REAL_ mpitmp(8) !Use for temporary packing of mpi messages.
  integer ist(MPI_STATUS_SIZE), partner
#else
  ! mdloop and REM is always 0 in serial
  integer, parameter :: mdloop = 0, rem = 0, remd_types(1) = 0, replica_indexes(1) = 0
#endif

  ! The following variables are needed since nstep and nstlim behave
  ! differently in a REMD run.  In certain places where output is required,
  ! total_nstep and total_nstlim take the place of nstep and nstlim. This
  ! allows replica runs to output less often than every exchange.  They are
  ! the absolute step # of the REMD or MD simulation.
  integer total_nstep, total_nstlim

  ! New variables for the generalized canonical-isokinetic algorithm:...
  _REAL_ s1(natom), s2(natom), s3(natom), glcl(natom), etlc(natom)
  _REAL_ sinsh2, sinsh4, sinsh8, tt, erlxt, erlxt2, erlixt2, tktk
  _REAL_ st, st1, aaa, bbb, esi, esim, erlst2, hkin, hhin, tkkk
  _REAL_ boltz, etlci, davalev, clfs, tkik, rndbim(3)
#ifdef MPI
  _REAL_ sdavalev(numtasks), sclfs, davalevs
  _REAL_ s1a(natom), s2a(natom), s3a(natom)
  _REAL_ s1s(natom), s2s(natom), s3s(natom)
  _REAL_ sv(3*natom), svs(3*natom)
#endif
  integer iii, kija, iseed(4)

  ! And related distribution functions:

  integer iconfs, iconfd, lep, kves, jl
#ifdef MPI
  integer inum, stxtr
#endif
  _REAL_ txtr, dargvs, distr, distra, distrmaxwel, pif1
  parameter(lep=250)
  parameter(pif1=3.141592653589793d0)
  _REAL_ argvel(-lep:lep), distrs(-lep:lep), distrh1(lep), distrh2(-lep:lep), &
         distrh3(-lep:lep)
#ifdef MPI
  _REAL_ sdistrs(-lep:lep), sdistrh1(lep), sdistrh2(-lep:lep), &
         sdistrh3(-lep:lep)
#endif
  LOGICAL :: exstf1,exstf2

  ! Stochastic Isokinetic Nose-Hoover RESPA integrator (SINR)
  type(sinr) :: sinrdata ! Variables for SINR
  ! Used to save averaged velocity information for on-step output
  _REAL_ onstepvel(3*natom)
  ! Used to save coordinates at time=t for on-step output
  _REAL_ onstepcrd(3*natom)

#include "../include/md.h"
#include "box.h"
#include "nmr.h"
#include "tgtmd.h"
#include "multitmd.h"
#include "../include/memory.h"
#include "extra.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "ew_mpole.h"
#include "def_time.h"
#include "extra_pts.h"
#include "../pbsa/pb_md.h"

#ifdef LES
  ! additional variables for PIMD and LES
  _REAL_  :: xcmd(3*natomCL), vcmd(3*natomCL), vscalt
  integer :: ncmd
#else
  _REAL_ sgsta_rndfp, sgend_rndfp, ignore_solvent
#endif
  ! for const press PIMD
  _REAL_ atomvir
  _REAL_ sysx, sysy, sysz, sysrange(3,2)
  logical mv_flag

  integer iatom, ierr
  ! Workaround for gfortran optimization bug that retains support for
  ! gfortran 4.2 and lower
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
  integer, volatile :: m
#else
  integer :: m
#endif

  _REAL_  Ekin2_tot,tmp
  integer :: idim
  _REAL_ :: E_nhc, exp1, exp2

  logical ivscm
  logical qspatial
  logical resetvelo
  _REAL_ etot_save,ekpbs

  logical do_list_update
  logical skip(*), belly, lout, loutfm, erstop, vlim, onstep
  ! Fortran does not guarantee short circuit logical expressions:
  logical is_remainder0
  _REAL_ x(*), winv(*), amass(*), f(*), v(*), vold(*), xr(*), xc(*), conp(*)
  type(state_rec) :: ener   ! energy values per time step
  type(state_rec) :: enert  ! energy values tallied over the time steps
  type(state_rec) :: enert2 ! energy values squared tallied over the time steps
  type(state_rec) :: enert_old, enert2_old
  type(state_rec) :: enert_tmp, enert2_tmp
  type(state_rec) :: edvdl
  type(state_rec) :: edvdl_r

#ifdef MPI
  type(state_rec) :: ecopy
  _REAL_ :: clfac, tmpvir(3,3)
#endif
  _REAL_ rmu(3), fac(3), onefac(3), etot_start
  _REAL_ tma(*)
  _REAL_ tspan, atempdrop, fln, scaltp
  _REAL_ vel, vel2, vcmx, vcmy, vcmz, vmax
  _REAL_ winf, aamass, rterm, ekmh, ekph, wfac, rsd
  _REAL_ fit, fiti, fit2

  ! Variables to control a Langevin dynamics simulation
  logical is_langevin
  _REAL_ gammai, c_implic, c_explic, c_ave, sdfac, ekins0
  _REAL_ dtx, dtxinv, dt5, factt, ekin0, ekinp0, dtcp, dttp
  _REAL_ rndf, rndfs, rndfp, boltz2, pconv, tempsu
  _REAL_ xcm(3), acm(3), ocm(3), vcm(3), ekcm, ekrot
  _REAL_ emtmd

  ! Variables and parameters for constant surface tension:
  ! ten_conv converts dyne/cm to bar angstroms
  _REAL_, parameter :: ten_conv = 100.0d0
  _REAL_  :: pres0x
  _REAL_  :: pres0y
  _REAL_  :: pres0z
  _REAL_  :: gamma_ten_int
  _REAL_  :: press_tan_ave

  integer nsp(*)
  integer idumar(4)
  integer l_temp
  integer i, j, im, i3, nitp, nits, iskip_start, iskip_end
  integer nstep, nrep, nrek, iend, istart3, iend3
  integer nrx, nr, nr3, ntcmt, izero, istart
  logical ixdump, ivdump, itdump, ifdump
  logical qsetup
  _REAL_, allocatable, dimension(:) :: f_or
#ifdef RISMSANDER
  logical irismdump
#  ifdef RISM_DEBUG
  _REAL_ r(3),cm(3),angvel(3),erot,moi,proj(3),rxv(3)
#  endif
#endif

  integer nvalid, nvalidi
  _REAL_ eke

  _REAL_, allocatable, dimension(:) :: frcti

  _REAL_ small
  data small/1.0d-7/

  !--- VARIABLES FOR DIPOLE PRINTING ---
  integer prndipngrp
  integer prndipfind
  character(len=4) prndiptest

  _REAL_,parameter :: pressure_constant = 6.85695d+4
  ! variables used in constant pressure PIMD
  _REAL_ :: Nkt,centvir,pressure, aa, arg2, poly, e2, e4, e6, e8
  _REAL_ :: box_center(3)

  ! variables used in middle scheme
  _REAL_ :: xold(3*natom)       ! for constrained MD
  _REAL_ :: x_centroid(3,natom) ! for PIMD
  ! Ek(t) is written in mden or mdout, as the default scheme does
  ! Ek(t+dt/2) is written in file_ek, as defined below
  integer,parameter :: fh_ek = 1024
  character(len=*),parameter :: file_ek = "eke.dat"
  _REAL_ :: ekhf, ekhf2

#ifdef MPI
  ! variable used in CMD
  _REAL_ :: tmp_eke_cmd ! Use for temporary packing of mpi messages.

  ! for adaptive qm/mm runs
  _REAL_ :: adqmmm_first_energy, etotcorr, tadc
  _REAL_ :: corrected_energy
  integer :: nstepadc
  logical :: flag_first_energy = .true.
#endif

  _REAL_ :: kinetic_E_save(2)
  integer :: aqmmm_flag

  ! PLUMED related variables
  _REAL_ :: plumed_box(3,3), plumed_virial(3,3), plumed_kbt
  integer :: plumed_version, plumed_stopflag
  _REAL_ :: plumed_energyUnits, plumed_timeUnits, plumed_lengthUnits
  _REAL_ :: plumed_chargeUnits
  ! }}}


  double precision      :: target_ekin               !< Target EKIN, used in Bussi's thermostat
#ifdef LES
  double precision      :: target_ekin_les           !< Target EKIN, used in Bussi's thermostat
#endif
  integer               :: target_ekin_update_nstep  !< Target EKIN update rate, used in Bussi's thermostat

  logical               :: update_bussi_target_kin_energy_on_current_step
  logical               :: update_kin_energy_on_current_step
!------------------------------------------------------------------------------
!  execution/initialization begins here:

  call trace_enter( 'runmd' )

  ! Initialize some variables {{{
  ixdump = .false.
  ivdump = .false.
  itdump = .false.
  ifdump = .false.
  onstepvel(:) = 0.d0
  onstepcrd(:) = 0.d0
#ifdef MPI
  if (master) then
    ! In Replica Exchange Molecular Dynamics (REMD), runmd will be called many
    ! times, so we dont want to open files every time.  For normal md, mdloop
    ! will just be 0.
    if (mdloop .eq. 0) then
      call amopen(7, mdinfo, 'U', 'F', facc)
    end if
  end if

  adqmmm_first_energy = 0.d0
#else
  call amopen(7, mdinfo, 'U', 'F', 'W')
#endif

  vlim = vlimit > small
  ntcmt = 0
  izero = 0
  belly = (ibelly > 0)
  lout = .true.
  loutfm = (ioutfm <= 0)
  nr = nrp
  nr3 = 3*nr
  ekmh = 0.d0

  aqmmm_flag = 0
  pressure = 0.d0
  etot_save = 0.d0
  E_nhc = 0.d0
  tktk = 0.d0
  sinsh8 = 0.d0
  sinsh4 = 0.d0
  sinsh2 = 0.d0
  erlxt2 = 0.d0
  erlxt = 0.d0
  erlixt2 = 0.d0
  pres0x = 0.d0
  pres0y = 0.d0
  pres0z = 0.d0
  gamma_ten_int = 0.d0
  dtcp = 0.d0
  dttp = 0.d0
  ekph = 0.d0
  ekpbs = 0.d0
  eke = 0.d0

#ifdef LES
  aa = 0.d0
  poly = 0.d0
  ekmhles = 0.d0
#endif
  do_list_update = .false.
#ifdef MPI
  if (mpi_orig) then
    istart = 1
    iend = natom
  else
    istart = iparpt(mytaskid) + 1
    iend = iparpt(mytaskid+1)
  end if
#else
  istart = 1
  iend = nr
#endif
  istart3 = 3*istart -2
  iend3 = 3*iend

#ifdef MPI
  if (icfe /= 0) then
    allocate(frcti(nr3 + 3*extra_atoms), stat=ierr)
    REQUIRE(ierr == 0)
  end if
#endif

  ! If ntwprt.NE.0, only print the atoms up to this value
  nrx = nr3
  if (ntwprt > 0) nrx = ntwprt*3
  if (.not. allocated(f_or)) allocate(f_or(nr3))
  if (abfqmmm_param%abfqmmm == 1) then
#ifdef MPI
    call xdist(v, xx(lfrctmp), natom)
#endif
    if (abfqmmm_param%system == 1) then
      if (abfqmmm_param%qmstep == 1) then
        abfqmmm_param%v(1:nr3+iscale) = v(1:nr3+iscale)
      end if
      v(1:nr3+iscale) = 0.d0
      t = t+dt
      if (abfqmmm_param%maxqmstep == 0) then
        t = 0
      end if
    else
      v(1:nr3+iscale) = abfqmmm_param%v(1:nr3+iscale)
    endif
  endif

  ! Cleanup the velocity if belly run
  if (belly) then
     call bellyf(nr,ix(ibellygp),v)
#ifdef MPI /* SOFT CORE */
     if ( ifsc /= 0 ) then
        call bellyf_softcore(v,vold,istart,iend)
     end if
#endif
  end if


!------------------------------------------------------------------------------
  ! Determine system degrees of freedom (for T scaling, reporting)
#   include "degcnt.inc"
!   }}}

!------------------------------------------------------------------------------
!    Pressure/temp units  {{{
  ! Begin unit conversion.  pconv eventually becomes a factor to convert
  ! pressure in kcal/mole-A^3 to bar.
  boltz2 = 8.31441d-3 * 0.5d0
  pconv = 1.6604345d+04
  boltz2 = boltz2/4.184d0
  dtx = dt*20.455d+00
  dtxinv = 1.0d0 / dtx
  dt5 = dtx * 0.5d0
  pconv = pconv*4.184d0

  ! fac() are #deg freedom * kboltz / 2.  Multiply by T to get the expected
  ! kinetic energy.  fac(1) is for the entire system.
  fac(1) = boltz2*rndf
  fac(2) = boltz2*rndfp
  if (rndfp < 0.1d0) fac(2) = 1.d-6

#ifdef LES
  ! Replaced solvent variables with LES ones
  ! since separate solvent coupling no longer used
  ! ASSUME SAME COUPLING CONSTANT FOR BOTH BATHS, just different target T

  ! will also have to accumulate LES and non-LES kinetic energies separately
  if (temp0les < 0.d0) then
    fac(3) = boltz2*rndfs
    if (rndfs < 0.1d0) fac(3) = 1.d-6
  else
    fac(3) = boltz2*rndfles
    if (rndfles < 0.1d0) fac(3) = 1.d-6
  end if
#else
  fac(3) = boltz2*rndfs
  if (rndfs < 0.1d0) fac(3) = 1.d-6
#endif
  if (ipimd == CMD) then
    if (eq_cmd) then
      fac(1) = boltz2 * dble(3*natomCL)
    else
      fac(1) = boltz2 * dble(3*(natomCL-1))
    endif
  endif
  onefac(1) = 1.0d0 / fac(1)
  onefac(2) = 1.0d0 / fac(2)
  onefac(3) = 1.0d0 / fac(3)
  factt = rndf/(rndf+ndfmin)

  ! These are "desired" kinetic energies based on the number of
  ! degrees of freedom and target temperature.  They will be used
  ! for calculating the velocity scaling factor
  ekinp0 = fac(2)*temp0
#ifdef LES

  ! Modified for LES temperature
  ekins0=0.d0
  ekinles0=0.d0
  if (temp0les < 0.d0) then
    ekins0 = fac(3) * temp0
    ekin0  = fac(1) * temp0
    if (master) &
      write (6,*) "Single temperature bath for LES and non-LES"
  else
    ekinles0 = fac(3)*temp0les
    ekin0  = ekinp0 + ekinles0
    if (master) then
      write (6,*) "LES particles coupled to separate bath"
      write (6,'(a,f8.2)')"    LES target temperature:    ",temp0les
      write (6,'(a,f8.2)')"    LES target kinetic energy: ",ekinles0
      write (6,'(a,f8.2)')"non-LES target temperature:    ",temp0
      write (6,'(a,f8.2)')"non-LES target kinetic energy: ",ekinp0
    end if
  end if
  target_ekin_les = ekinles0
#else
  ekins0 = fac(3)*temp0
  ekin0  = fac(1)*temp0
#endif
  target_ekin = ekin0
  ! }}}

  if (ischeme > 0 .and. rem == 0 .and. master) then
     open(fh_ek, file=file_ek)
     write(fh_ek, *) '#     NSTEP   EKIN(kcal/mol)    TEMP(K)'
  endif

!------------------------------------------------------------------------------
!    Nose-Hoover initialization {{{
#ifdef LES
  if (abfqmmm_param%abfqmmm /= 1) then
    if (ntt == 4) call nose_hoover_init_LES(amass,v,f)
  else
    if (ntt == 4 .and. abfqmmm_param%qmstep == 1 .and. &
        abfqmmm_param%system == 1) &
      call nose_hoover_init_LES(amass,abfqmmm_param%v,abfqmmm_param%f)
  end if
#else
  if (abfqmmm_param%abfqmmm .ne. 1) then
    if (ntt >= 4 .and. ntt <= 8) call nose_hoover_init(amass, v, f)
  else
    if (ntt >= 4 .and. ntt <= 8 .and. abfqmmm_param%qmstep == 1 .and. &
        abfqmmm_param%system == 1) &
      call nose_hoover_init(amass, abfqmmm_param%v, abfqmmm_param%f)
  endif
#endif /* LES */
!  }}}

!------------------------------------------------------------------------------
  !    Langevin dynamics setup  {{{
  is_langevin = (gamma_ln > 0.0d0)
  gammai = gamma_ln / 20.455d0
  c_implic = 1.d0 / (1.d0 + gammai*dt5)
  c_explic = 1.d0 - gammai*dt5
  c_ave = 1.d0 + gammai*dt5
  sdfac = sqrt(4.d0 * gammai * boltz2 * temp0 / dtx)
#ifdef LES
  if (temp0les < 0.d0) then
    sdfacles = sqrt( 4.d0*gammai*boltz2*temp0/dtx )
  else
    sdfacles = sqrt( 4.d0*gammai*boltz2*temp0les/dtx )
  endif
#endif /* LES */
  if (ischeme == 0 .and. is_langevin .and. ifbox==0) then
    call get_position(nr, x, sysx, sysy, sysz, sysrange, 0)
#ifdef MPI
    ! Soft core position mixing
    if (ifsc == 1) call sc_mix_position(sysx, sysy, sysz, clambda)
#endif /* MPI */
  end if
  ! }}}
  ! middle scheme when nscm is enabled
  if (nscm > 0 .and. ischeme > 0) then
     if (ifbox == 0) then
        if (ipimd == 0) then
           call get_position(nr, x, sysx, sysy, sysz, sysrange, 0)
#ifdef MPI
    ! Soft core position mixing
           if (ifsc == 1) call sc_mix_position(sysx, sysy, sysz, clambda)
#endif /* MPI */
        else if (ipimd > 0) then
#ifdef MPI
#ifndef LES
           REQUIRE( allocated(xall) )
           if (sanderrank.eq.0) then
              call MPI_Allgather(x   ,3*natom,MPI_DOUBLE_PRECISION,&
                                 xall,3*natom,MPI_DOUBLE_PRECISION,&
                                 commmaster,ierr) 
           endif
           
           call MPI_Bcast(xall,3*natom*nbead,MPI_DOUBLE_PRECISION,&
              0, commsander,ierr)
           x_centroid(1:3,1:natom)=0.0
           do i = 1, nbead
              do j = istart, iend
                 x_centroid(1:3,j)=x_centroid(1:3,j)+xall(1:3,j,i)
              end do
           end do
           x_centroid(1:3,1:natom) = x_centroid(1:3,1:natom)/dble(nbead)
           ! position of mass center
           sysx = 0.d0
           sysy = 0.d0
           sysz = 0.d0
           tmp = 0.d0
           do i = 1, natom
             aamass = amass(i)
             sysx = sysx + aamass*x_centroid(1,i)
             sysy = sysy + aamass*x_centroid(2,i)
             sysz = sysz + aamass*x_centroid(3,i)
           end do
           sysx = sysx * tmassinv * dble(nbead)
           sysy = sysy * tmassinv * dble(nbead)
           sysz = sysz * tmassinv * dble(nbead)
#endif
#endif 
        end if ! pimd or not
     end if ! ifbox==0: non-periodic
  end if ! middle scheme when nscm is enabled

!------------------------------------------------------------------------------
  ! Constant pH and constant Redox potential setup: {{{
  if ((icnstph /= 0 .or. (icnste /= 0 .and. cpein_specified)) .and. mdloop .eq. 0) call cnstphinit(x, ig)
  if (icnste /= 0 .and. .not. cpein_specified .and. mdloop .eq. 0) call cnsteinit(x, ig)
  if (ntt == 1 .or. ntt == 11) dttp = dt/tautp
  if (ntp > 0) dtcp = comp * 1.0d-06 * dt / taup
  ! }}}

!------------------------------------------------------------------------------
  !    Constant surface tension setup: {{{
  if (csurften > 0) then

    ! Set pres0 in direction of surface tension.
    ! The reference pressure is held constant in on direction dependent
    ! on what the surface tension direction is set to.
    if (csurften .eq. 1) then           ! pres0 in the x direction
      pres0x = pres0
    else if (csurften .eq. 2) then      ! pres0 in the y direction
      pres0y = pres0
    else                                ! pres0 in the z direction
      pres0z = pres0
    end if

    ! Multiply surface tension by the number of interfaces
    gamma_ten_int = dble(ninterface) * gamma_ten
  else if (csurften < 0) then

    ! Set pres0 in out-of-plane direction
    if (csurften .eq. 1) then           ! pres0 in the x direction
      pres0x = pres0
    else if (csurften .eq. 2) then      ! pres0 in the y direction
      pres0y = pres0
    else                                ! pres0 in the z direction
      pres0z = pres0
    end if

  end if
  ! }}}

  !   General initialization:  {{{
  nrek = 4
  nrep = 15

  nvalid = 0
  nvalidi = 0
  nstep = 0
  total_nstep = 0
#ifdef MPI
  ! For REMD, total_nstep is the number of steps * the number of exchanges
  ! we've already attempted
  if (rem .ne. 0) total_nstep = (mdloop - 1) * nstlim
#endif /* MPI */
  fit = 0.d0
  fiti = 0.d0
  fit2 = 0.d0

!------------------------------------------------------------------------------
  ! Zero all elements of these sequence types
  ener       = null_state_rec
  enert      = null_state_rec
  enert2     = null_state_rec
  enert_old  = null_state_rec
  enert2_old = null_state_rec
  edvdl      = null_state_rec
  edvdl_r    = null_state_rec

!------------------------------------------------------------------------------
  ! For Path Integral Molecular Dynamics (PIMD), Normal Mode PIMD, Centroid
  ! Centroid Molecular Dynamics (CMD), or Ring Polymer Molecular Dynamics
  ! (RPMD):
  totenert   = null_state_rec
  totenert2  = null_state_rec

  ener%kin%pres_scale_solt = 1.d0
  ener%kin%pres_scale_solv = 1.d0
  ener%box(1:3) = box(1:3)
  ener%cmt(1:4) = 0.d0
  nitp = 0
  nits = 0
  if (ischeme > 0) then
     ekhf = 0.0d0
     ekhf2 = 0.0d0
  endif
  ! }}}

!------------------------------------------------------------------------------
  ! PLUMED initialization.  PLUMED is an open-source plugin that
  ! confers the functionality of a number of enhanced sampling methods.
  if (plumed == 1) then
#   include "Plumed_init.inc"
  endif

!------------------------------------------------------------------------------
  ! Transform cartesian positions into normal mode positions for
  ! (Normal Mode) Path Integral Molecular Dynamics (PIMD)
  if (ipimd == NMPIMD .or. ipimd == CMD) call trans_pos_cart_to_nmode(x)

!------------------------------------------------------------------------------
  ! Make a first dynamics step. {{{
  ! init = 3: general startup if not continuing a previous run
  if (init == 3 .or. nstlim == 0 .or. &
      (abfqmmm_param%abfqmmm == 1 .and. abfqmmm_param%system == 1)) then

    ! Constant-pressure MD, when we are not also dealing with AMOEBA or PIMD:
    ! Calculate the center of mass for each molecule, kinetic energy of the
    ! molecule's center of mass, and coordinates of each molecule relative to
    ! its center of mass in the simulation.
    if (ntp > 0 .and. iamoeba == 0 .and. ipimd == 0) then
      xr(1:nr3) = x(1:nr3)
      call ekcmr(nspm, nsp, tma, ener%cmt, xr, v, amass, 1, nr)
    end if

    ! Calculate the force.  Set irespa to get full
    ! energies calculated on step "0":
    irespa = 0
    iprint = 1

!------------------------------------------------------------------------------
    ! for path-integral MD:

    if (ipimd == NMPIMD .or. ipimd == CMD) then
#if defined(MPI) || defined(LES)
      call trans_pos_nmode_to_cart(x, cartpos)
#else
      call trans_pos_nmode_to_cart(cartpos)
#endif /* MPI or LES */
      call force(xx, ix, ih, ipairs, cartpos, f, ener, ener%vir, xx(l96), &
                 xx(l97), xx(l98), xx(l99), qsetup, do_list_update,nstep)
#if defined(MPI) && defined(LES)
      if (ievb == 1 .and. i_qi > 0) then
        call evb_umb(f, cartpos, real_mass, natom, istart3, iend3)
        if (i_qi == 2) then
          call qi_corrf_les(cartpos, real_mass)
        end if
        evb_nrg(1) = evb_frc%evb_nrg
        evb_nrg(2) = evb_vel0%evb_nrg
        if (nbias > 0) evb_nrg(3) = sum(evb_bias%nrg_bias(:))
      end if
#endif /* MPI and LES */
      call trans_frc_cart_to_nmode(f)
      i3 = 3*(istart-1)
#if defined(MPI) && defined(LES)
      if (ievb /= 0 .and. i_qi == 0) then
        call evb_umb(f, x, real_mass, natom, istart3, iend3)
        evb_nrg(1) = evb_frc%evb_nrg
        evb_nrg(2) = evb_vel0%evb_nrg
        if (nbias > 0) evb_nrg(3) = sum(evb_bias%nrg_bias(:))
      end if
#endif /* MPI and LES */

!------------------------------------------------------------------------------
    ! for LSC-IVR:

    else if (ilscivr == 1) then

      ! Prepare the Hessian Matrix of the potential for the Linearized
      ! Semi-Classical Initial Value Representation (LSC-IVR).  At this
      ! point, x is the position of a bead at equilibrium.  Initialize
      ! the LSC-IVR variables.
      natom_lsc = natom
      ndof_lsc = natom * 3
      call lsc_init
      do ilsc = 1, natom_lsc
        mass_lsc(3*ilsc-2) = amass(ilsc)
        mass_lsc(3*ilsc-1) = amass(ilsc)
        mass_lsc(3*ilsc  ) = amass(ilsc)
      end do
      v2_lsc = 0.0d0
      do ilsc = 1, ndof_lsc

        ! ith vector of the Hessian matrix
        x_lsc = 0.0d0
        x_lsc(1:ndof_lsc) = x(1:ndof_lsc)
        x_lsc(ilsc) = x(ilsc) + dx_lsc
        call force(xx, ix, ih, ipairs, x_lsc, f_lsc, ener, ener%vir, xx(l96), &
                   xx(l97), xx(l98), xx(l99), qsetup, do_list_update, nstep)
#ifdef MPI
        call xdist( f_lsc, xx(lfrctmp), natom )
#endif /* MPI */
        v2_lsc(1:ndof_lsc,ilsc) = f_lsc(1:ndof_lsc)
      enddo

      call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
                 xx(l98), xx(l99), qsetup, do_list_update,nstep)
#ifdef MPI
      call xdist(f, xx(lfrctmp), natom)
#endif /* MPI */

      ! Second derivative of the potential
      do ilsc = 1, ndof_lsc
        v2_lsc(1:ndof_lsc,ilsc) = &
          (f(1:ndof_lsc) - v2_lsc(1:ndof_lsc,ilsc))/dx_lsc
      end do

      ! Get the initial position of the momentum
      call lsc_xp(x, v)

    else

!------------------------------------------------------------------------------
      ! Thermodynamic Integration (TI) decomposition
      if (idecomp > 0 .and. ntpr > 0) then
        decpr = .false.
        if (mod(nstep+1, ntpr) == 0) decpr = .true.
      end if

!------------------------------------------------------------------------------
      call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
                 xx(l98), xx(l99), qsetup, do_list_update,nstep)
#ifdef MPI
      if (ievb /= 0) then
#  ifdef LES
        call evb_umb_primitive(f, x, real_mass, natom, istart, iend)
#  else
        call evb_umb_primitive(f, x, amass, natom, istart, iend)
#  endif /* LES */
        evb_nrg(1) = evb_frc%evb_nrg
        evb_nrg(2) = evb_vel0%evb_nrg
        if (nbias > 0) evb_nrg(3) = sum(evb_bias%nrg_bias(:))
      endif
#endif /* MPI */
    endif
    ! End branches into Normal Mode Path Integral MD, Centroid MD, Linearized
    ! Semi-Classical Initial Value Representation

!------------------------------------------------------------------------------
    ! Semi-Empirical Born-Oppenheimer Molecular Dynamics:
    if (sebomd_obj%do_sebomd) then

      ! Computes the Hessian matrix if necessary
      if (sebomd_obj%ntwh /= 0) then

        ! Don't output atomic charges and bond orders
        sebomd_obj%iflagch_old = sebomd_obj%iflagch
        sebomd_obj%iflagch = 0
        sebomd_obj%iflagbo_old = sebomd_obj%iflagbo
        sebomd_obj%iflagbo = 0
        call sebomd_gradient_write(f, 3*natom)
        call sebomd_hessian_compute(xx, ix, ih, ipairs, x, ener, &
                                    qsetup, do_list_update, nstep)
        sebomd_obj%iflagch = sebomd_obj%iflagch_old
        sebomd_obj%iflagbo = sebomd_obj%iflagbo_old
      endif
    endif
    ! End SEBOMD

!------------------------------------------------------------------------------
    ! Write information concerning constant PH molecular dynamics
    if ((icnstph /= 0 .or. (icnste /= 0 .and. cpein_specified)) .and. master .and. &
      ((rem /= 0 .and. mdloop > 0) .or. rem == 0)) then
      call cnstphwrite(rem,remd_types,replica_indexes)
    end if

!------------------------------------------------------------------------------
    ! Write information concerning constant Redox potential molecular dynamics
    if (icnste /= 0 .and. .not. cpein_specified .and. master .and. &
      ((rem /= 0 .and. mdloop > 0) .or. rem == 0)) then
      call cnstewrite(rem,remd_types,replica_indexes)
    end if

!------------------------------------------------------------------------------
    ! Needed for Adaptive Biasing-Force mixing QM/MM
    f_or(1:nr3) = f(1:nr3)
#ifdef MPI
    call xdist(f_or, xx(lfrctmp), natom)
#endif /* MPI */
    if (abfqmmm_param%abfqmmm == 1) then
      if (abfqmmm_param%system == 1) then
        abfqmmm_param%f1(1:nr3) = f_or(1:nr3)
      end if
      if (abfqmmm_param%system == 2) then
        abfqmmm_param%f2(1:nr3) = f_or(1:nr3)
        call abfqmmm_combine_forces()
        f_or(1:nr3) = abfqmmm_param%f(1:nr3)
        f(1:nr3) = abfqmmm_param%f(1:nr3)
      end if
    end if

    ! This force call does not count as a "step". CALL NMRDCP to decrement
    ! local NMR step counter and MTMDUNSTEP to decrease the local MTMD step
    ! counter
    call nmrdcp
    call mtmdunstep

    ! PLUMED force is added in this routine.
    plumed_stopflag=0
    if (plumed == 1) then
#     include "Plumed_force.inc"
    end if

#ifdef MPI /* SOFT CORE */
!------------------------------------------------------------------------------
    ! If softcore potentials are used, collect their dvdl contributions:
    if (ifsc /= 0) then
      call mpi_reduce(sc_dvdl, sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0, commsander, ierr)

      ! Zero dV/dLambda for the next stetp
      sc_dvdl=0.0d0
      call mpi_reduce(sc_dvdl_ee, sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0, commsander, ierr)

      ! Zero for the next step
      sc_dvdl_ee=0.0d0
      call mpi_reduce(sc_ener, sc_ener_tmp, ti_ene_cnt, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0, commsander, ierr)
      sc_ener(1:ti_ene_cnt) = sc_ener_tmp(1:ti_ene_cnt)
    end if
    if (ifsc == 2) then

      ! If this is a perturb to nothing run, scale forces and calculate dvdl
      call sc_nomix_frc(f, nr3, ener)
      if (numtasks > 1) then
        call mpi_bcast(f, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(ener, state_rec_len, MPI_DOUBLE_PRECISION, 0, &
                       commsander, ierr)
      end if
    end if

    if (icfe /= 0) then
!------------------------------------------------------------------------------
      ! Free energies using thermodynamic integration (icfe /= 0)
      if (master) then

        ! First, partner threads exchange forces and energies
        partner = ieor(masterrank, 1)
        call mpi_sendrecv(f, nr3, MPI_DOUBLE_PRECISION, partner, 5, &
                          frcti, nr3+3*extra_atoms, MPI_DOUBLE_PRECISION, &
                          partner, 5, commmaster, ist, ierr )
        call mpi_sendrecv(ener, state_rec_len, MPI_DOUBLE_PRECISION, partner, &
                          5, ecopy, state_rec_len, MPI_DOUBLE_PRECISION, &
                          partner, 5, commmaster, ist, ierr)

        ! Exchange sc-dvdl contributions between masters
        call mpi_sendrecv(sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                          sc_tot_dvdl_partner, 1, MPI_DOUBLE_PRECISION, &
                          partner, 5, commmaster, ist, ierr)
        call mpi_sendrecv(sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, partner, &
                          5, sc_tot_dvdl_partner_ee, 1, MPI_DOUBLE_PRECISION, &
                          partner, 5, commmaster, ist, ierr )
        if (masterrank == 0) then
          call mix_frcti(frcti, ecopy, f, ener, nr3, clambda, klambda)
        else
          call mix_frcti(f, ener, frcti, ecopy, nr3, clambda, klambda)
        end if
      end if

      if (numtasks > 1) then
        call mpi_bcast(f, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(ener, state_rec_len, MPI_DOUBLE_PRECISION, 0, &
                       commsander, ierr)
      end if
    end if
#endif /* MPI SOFT CORE */
    irespa = 1

!------------------------------------------------------------------------------
    ! Reset quantities depending on temp0 and tautp (which may have been
    ! changed by modwt during the call to force()).  Recalculate target
    ! kinetic energies.
    ekinp0 = fac(2) * temp0

#ifdef LES

    ! Modified for LES temperature, not solvent
    ekins0 = 0.d0
    ekinles0 = 0.d0
    if (temp0les < 0.d0) then
      ekins0 = fac(3) * temp0
      ekin0 = fac(1) * temp0
    else
      ekinles0 = fac(3) * temp0les
      ekin0 = ekinp0 + ekinles0
    end if
    
    target_ekin_les = ekinles0
#else
    ekins0 = fac(3) * temp0
    ekin0 = fac(1) * temp0
#endif
    target_ekin = ekin0

    if (ntt == 1 .or. ntt == 11) dttp = dt / tautp
    if (ntp > 0) then
      ener%volume = volume
      ener%density = tmass / (0.602204d0*volume)
      if (iamoeba == 0) then
        ener%cmt(4) = 0.d0
        ener%vir(4) = 0.d0
        ener%pres(4) = 0.d0
        do m = 1,3
          ener%cmt(m)  = ener%cmt(m) * 0.5d0
          ener%cmt(4)  = ener%cmt(4) + ener%cmt(m)
          ener%vir(4)  = ener%vir(4) + ener%vir(m)
          ener%pres(m) = (pconv+pconv) * (ener%cmt(m)-ener%vir(m)) / volume
          ener%pres(4) = ener%pres(4) + ener%pres(m)
        end do
        ener%pres(4) = ener%pres(4) / 3.d0
      end if
    end if
    ntnb = 0
    i3 = 0
    tempsu = 0.0d0

#ifdef LES
    ! Added LES tempsu (actual LES sum of m*v**2 )
    tempsules = 0.0d0
#endif
    eke_cmd = 0.d0
    do j = 1,nrp
      winf = winv(j) * dt5
      aamass = amass(j)
      do m = 1,3
        i3 = i3+1
        rterm = v(i3)*v(i3) * aamass
#ifdef LES
        if (temp0les < 0.d0) then
          tempsu = tempsu + rterm
          if (ipimd == CMD .and. (cnum(j) == 0 .or. cnum(j) == 1)) then
            eke_cmd = eke_cmd + aamass*v(i3)*v(i3)
          endif
        else
          if (cnum(j) == 0) then
            tempsu = tempsu + rterm
          else
            tempsules = tempsules + rterm
          end if
        end if
#else
        if (ipimd == CMD .and. mybeadid == 1) &
          eke_cmd = eke_cmd + aamass*v(i3)*v(i3)
        tempsu = tempsu + rterm
#endif
        if (ipimd .ne. NMPIMD .and. ipimd .ne. CMD) &
          v(i3) = v(i3) - f(i3) * winf
        if (vlim) v(i3) = sign(min(abs(v(i3)), vlimit), v(i3))
      end do
    end do

#ifdef MPI /* SOFT CORE */
    if (ifsc /= 0) then
      call calc_softcore_ekin(amass, v, v, istart, iend)
      sc_ener(13) = sc_ener(6) + sc_ener(12)
    end if
#endif

    do im = 1, iscale
      tempsu = tempsu + scalm * v(nr3+im)*v(nr3+im)
      v(nr3+im) = v(nr3+im) - f(nr3+im) * dt5 / scalm
    end do
    ener%kin%solt = tempsu * 0.5d0

    ! middle scheme for constrain MD
    if (ischeme == 1 .and. ntc /= 1) then
       qspatial = .false.
       ! RATTLE-V, correct velocities
       call rattlev(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
       winv,conp,skip,x,v,nitp,belly,ix(iifstwt),ix(noshake), qspatial)
 
       ! use SETTLE to deal with water model
       call quick3v(x, v, ix(iifstwr), natom, nres, ix(i02))
    end if

#ifdef LES
    ! Added for LES temperature using old solvent variable for ener(4)
    if (temp0les < 0.d0) then
      ener%kin%solv = 0.d0
      ener%kin%tot = ener%kin%solt
      ! For CMD: "virial" estimate of KE
      if (ipimd > 0) then
        ener%kin%solv = equal_part + Epot_deriv
        ener%tot = ener%kin%solv + ener%pot%tot
      else
        ener%tot = ener%kin%tot + ener%pot%tot
      endif
      if (ipimd == CMD) then
        ener%kin%tot  = eke_cmd*0.5d0
        ener%kin%solv = ener%kin%tot
      endif
    else
      ener%kin%solv = tempsules * 0.5d0
      ener%kin%tot = ener%kin%solt + ener%kin%solv
    end if
#else
    ! For better output for parallel PIMD/NMPIM/CMD/RPMD
    if (ipimd > 0) then
      ener%tot = 0.d0
      ener%kin%tot = 0.d0
      ener%kin%solt = 0.d0
      ener%kin%solv = 0.d0
      ener%volume = 0.d0
    endif
    ener%kin%tot = ener%kin%solt
    ener%tot = ener%kin%tot + ener%pot%tot
#endif /* LES */

    if (ntt == 1) then
#ifdef LES
      if (temp0les >= 0.d0) then
        ekmh = max(ener%kin%solt, fac(2)*10.d0)
        ekmhles = max(ener%kin%solv, fac(3)*10.d0)
      else
        ekmh = max(ener%kin%solt,fac(1)*10.d0)
      end if
#else
      ekmh = max(ener%kin%solt,fac(1)*10.d0)
#endif /* LES */
    end if
  end if
  ! This ends a HUGE conditional branch in which init == 3, general startup
  ! when not continuing a previous dynamics run.  }}}

!------------------------------------------------------------------------------
  ! Continue init=3 or start init=4: {{{
  ! What follows applies to the case of init = 4: continuation of a previous
  ! trajectory.  This will also done for init = 3 (that code is simply run as
  ! startup and effectively performs a dynamics step to prime the pump).
  !
  ! Note: if the last printed energy from the previous trajectory was at time
  ! "t", then the restrt file has velocities at time t + 0.5dt and
  ! coordinates at time t + dt.
  ekmh = 0.0d0
#ifdef LES
  ekmhles = 0.0d0
#endif /* LES */

  i3 = 0
  do j = 1, nrp
    aamass = amass(j)
    do m = 1, 3
      i3 = i3+1
      rterm = v(i3) * v(i3) * aamass
#  ifdef LES
      ! use copy number, not solute/solvent
      if (temp0les < 0.d0) then
        ! 1 bath
        ekmh = ekmh + rterm
      else
        if (cnum(j) == 0) then
          ekmh = ekmh + rterm
        else
          ekmhles = ekmhles + rterm
        end if
      end if
#  else
      ekmh = ekmh + rterm
#  endif /* LES */
    end do
  end do

#ifdef MPI /* SOFT CORE */
  if (ifsc /= 0) then
    call calc_softcore_ekin(amass, v, v, istart, iend)
    sc_ener(13) = sc_ener(6) + sc_ener(12)
  end if
#endif /* MPI */

   do im = 1, iscale
      ekmh = ekmh + scalm*v(nr3+im)*v(nr3+im)
   end do
   ekmh = ekmh * 0.5d0
#ifdef LES
   ekmhles = ekmhles * 0.5d0
#endif /* LES */

  vold(1:nr3+iscale) = v(1:nr3+iscale)

#ifdef EMIL
  ! Setup the emil calculation if required.  EMIL is a
  ! sort of thermodynamic integration tool.
  if (emil_do_calc .gt. 0) &
    call emil_init(natom, 1.0/(temp0 * 2 * boltz2 ), &
                   mass, xx(lcrd), f, v, ener%box)
#endif /* EMIL */

  ! Adjust the step count if Adaptive Buffered Force QM/MM is in effect
  if (abfqmmm_param%abfqmmm == 1) then
    nstep=abfqmmm_param%qmstep
    if (abfqmmm_param%maxqmstep == 0) nstep = 0
  end if
  ! }}}

!------------------------------------------------------------------------------
  ! If init .ne. 4, or only one step, or ABF QM/MM {{{
  if (init .ne. 4 .or. nstlim == 0 .or. &
      (abfqmmm_param%abfqmmm == 1 .and. abfqmmm_param%system == 1)) then

    ! Print the initial energies and temperatures
#ifdef RISMSANDER

    if (rismprm%rism == 1 .and. rismprm%write_thermo==1 .and. &
        nstep <= 0 .and. facc .ne. 'A') then
      if (rism_calc_type(0) == RISM_FULL) then
        if (nstlim == 0) then
          call rism_solvdist_thermo_calc(.true., 0)
        else
          call rism_solvdist_thermo_calc(.false., 0)
        end if
      end if
    end if
#endif /* RISMSANDER */
    is_remainder0 = .false.
    if ( ntpr > 0 ) is_remainder0 = mod(abfqmmm_param%qmstep,ntpr) == 0
    if ( (nstep <= 0 .and. master .and. facc .ne. 'A') .or. &
        (master .and. abfqmmm_param%abfqmmm == 1 .and. is_remainder0) ) then
      if (isgld > 0) call sgenergy(ener)
      rewind(7)

#ifdef LES
      if (.not. ipimd > 0) ener%tot = ener%kin%tot+ener%pot%tot
#endif /* LES */

      if (abfqmmm_param%abfqmmm /= 1 .or. &
          abfqmmm_param%system == 1 .or. nstep == 0) then
        call prntmd(nstep, t, ener, onefac, 7, .false.)
      end if

#ifdef MPI /* SOFT CORE */
      if (ifsc .ne. 0) call sc_print_energies(6, sc_ener)
      if (ifsc .ne. 0) call sc_print_energies(7, sc_ener)
#endif

      if (ifcr > 0 .and. crprintcharges > 0) &
        call cr_print_charge(xx(l15), nstep)

      ! Begin dipole printing code  {{{
      call nmlsrc('dipoles', 5, prndipfind)
      if (prndipfind /= 0) then
        write(6,*) '------------------------------- DIPOLE &
                    &INFO ----------------------------------'
        write(6,9018) nstep, t
        read (5,'(a)') prndiptest
        call rgroup(natom, natc, nres, prndipngrp, ix(i02), ih(m02), &
                    ih(m04), ih(m06), ih(m08), ix(icnstrgp), jgroup, indx, &
                    irespw, npdec, xx(l60), 0, 0, 0, idecomp, 5, .false.)
        rewind(5)
        if (prndipngrp > 0) then
          call printdip(prndipngrp, ix(icnstrgp), xx(lcrd), &
                        xx(l15), xx(linddip), xx(Lmass), natom)
        end if
        write(6,*) '----------------------------- END DIPOLE &
                    &INFO --------------------------------'
      end if
      !--- END DIPOLE PRINTING CODE --- }}}

      if (nmropt > 0) call nmrptx(6)
      if (infe == 1) call nfe_prt(6)
      call amflsh(7)
    end if
    if (abfqmmm_param%abfqmmm == 1 .and. abfqmmm_param%system == 1) then
      deallocate(f_or, stat=ierr)
      REQUIRE(ierr == 0)
      return
    end if
    if (nstlim == 0) then
      if (abfqmmm_param%abfqmmm == 1) then
        v(1:nr3) = abfqmmm_param%v(1:nr3)
      end if
#ifdef MPI
      call xdist(x, xx(lfrctmp), natom)
      call xdist(v, xx(lfrctmp), natom)
#endif
      if (master) then
        if ( nstlim > 0 ) then
           call mdwrit(nstep, nr, ntxo, ntb, x, v, t, temp0,solvph,solve)
        end if
        if (ntwx>0) call corpac(x, 1, nrx, MDCRD_UNIT, loutfm)
        if (ntwv>0) call corpac(v, 1, nrx, MDVEL_UNIT, loutfm)
        if (ntwf>0) call corpac(f_or, 1, nrx, MDFRC_UNIT, loutfm)
        if (ntwe>0) call mdeng(15, nstep, t, ener, onefac, ntp, csurften)
      end if
      return
    end if
    init = 4
  end if
  ! End of contingencies primarily related to init not equal to 4 }}}

  if (ntp > 0 .and. ipimd > 0) then
    REQUIRE(ipimd == NMPIMD)
#ifdef LES
    call part_setup_cnst_press_pimd(Nkt, tau_vol)
#else
    call full_setup_cnst_press_pimd(Nkt, tau_vol)
#endif
    e2 = 1.0 / (2.0*3.0)
    e4 = e2 / (4.0*5.0)
    e6 = e4 / (6.0*7.0)
    e8 = e6 / (8.0*9.0)
    x_lnv = log(box(1)*box(2)*box(3)) / 3
  end if

!------------------------------------------------------------------------------
  ! For Centroid MD {{{
  if (ipimd == CMD) then
    if (.not. eq_cmd) then

      ! De-activate thermostat for path-centroid.
#ifdef LES
      do iatom = 1, natom
        do idim = 1, 3
          if (cnum(iatom) == 0 .or. cnum(iatom) == 1) then
            activate = .false.
          else
            activate = .true.
          end if
          call Thermostat_switch(thermo(idim,iatom), activate)
        enddo
      enddo
      if (.not. restart_cmd) then

        ! Scale path-centroid velocity and set total momentum equal to zero.
        call part_scale_vel_centroid(v, amass, istart, iend)
        nstep_cmd = 0
        t_cmd = 0.d0
      else
        t_cmd = t
        nstep_cmd = int( t / dt )
      end if
#else
      if (mybeadid == 1) then
        activate = .false.
      else
        activate = .true.
      end if
      do iatom = 1, natom
        do idim  = 1, 3
          call Thermostat_switch(thermo(idim, iatom), activate)
        enddo
      enddo
      if (.not. restart_cmd) then

        ! Scale path-centroid velocity and set total momentum equal to zero.
        call full_scale_vel_centroid(v,amass,istart,iend)
        nstep_cmd = 0
        t_cmd = 0.d0
      else
        nstep_cmd = nstep
        t_cmd = t
      end if
#endif /* LES */
    else
      nstep_cmd = nstep
      t_cmd = t
    end if
  end if
  ! End of contingency for Centroid MD with constant pressure conditions }}}

#ifdef MPI
  ! If this is a replica run and we are on exchange > 1, restore the
  ! old ekmh value since it was reset after we left runmd last time.
  ! DAN ROE: Only for ntt==1??
  if (rem /= 0 .and. mdloop >= 1) ekmh = remd_ekmh
#endif

!------------------------------------------------------------------------------
  if(ntt.eq.9) then
  !   initialize the OIN integrator: {{{
  !  Initialization of variables for the canonical-isokinetic algorithm:

      boltz=2.d0*boltz2
      erlxt=1.d0/gammai
      erlxt2=erlxt**2
      erlixt2=1.d0/erlxt2
      sinsh2=dtx/2.d0/nkija
      sinsh4=sinsh2/2.d0
      sinsh8=sinsh2/4.d0
      tktk=temp0*boltz
      davalev=0.d0
      clfs=0.d0

!     Read or generate the atomic velocities and chain-thermostat variables:

      if(irest.eq.1) then
         INQUIRE(FILE='vfreez.rst',EXIST=exstf1)
         INQUIRE(FILE='tfreez.rst',EXIST=exstf2)
      end if

!     Reading the atomic velocities and chain-thermostat variables:

#ifdef MPI
      if(master.and.irest.eq.1.and.exstf1.and.exstf2) then
         write(6,576) 'Using the existing files vfreez.rst and tfreez.rst', &
                'to read initial velocities and thermostat variables'
  576       format(/,a,/,a,/)

         open(578,file='vfreez.rst')
         iii = 0
         do j=1,natom
            read(578,*) svs(iii+1),svs(iii+2),svs(iii+3)
            iii = iii+3
         end do
         close(578)

         open(579,file='tfreez.rst')
         do j=1,natom
            read(579,*) s1s(j),s2s(j),s3s(j)
         end do
         close(579)

      end if

      if(irest.eq.1.and.exstf1.and.exstf2) then

         v(1:3*natom)=0.d0

         call mpi_allreduce(svs,v,3*natom,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)

         s1=0.d0
         s2=0.d0
         s3=0.d0

         call mpi_allreduce(s1s,s1,natom,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
         call mpi_allreduce(s2s,s2,natom,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
         call mpi_allreduce(s3s,s3,natom,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)

      end if

#else

      if(irest.eq.1.and.exstf1.and.exstf2) then
         write(6,576) 'Using the existing files vfreez.rst and tfreez.rst', &
                'to read initial velocities and thermostat variables'
  576       format(/,a,/,a,/)

         open(578,file='vfreez.rst')
         iii = 0
         do j=1,natom
            read(578,*) v(iii+1),v(iii+2),v(iii+3)
            iii = iii+3
         end do
         close(578)

         open(579,file='tfreez.rst')
         do j=1,natom
            read(579,*) s1(j),s2(j),s3(j)
         end do
         close(579)

      end if

#endif

      if(irest.eq.0.or..not.exstf1.or..not.exstf2) then

         if(master) write(6,577) 'Using random generation of initial &
                                 &velocities and thermostat variables'
  577       format(/,a,/)

!        At the very beginning of the simulations the velocities of atoms
!        are generated randomly to satisfy the kinetic energy constraint

         iseed(1)=0
         iseed(2)=1
         iseed(3)=2
         iseed(4)=3

         iii=0

         do j=1,natom

!           Random generation within the allowed interval

            call dlarnv(2,iseed,3,rndbim)

            v(iii+1)=rndbim(1)*dsqrt(tktk/amass(j))
            v(iii+2)=rndbim(2)*dsqrt(tktk/amass(j))
            v(iii+3)=rndbim(3)*dsqrt(tktk/amass(j))

!           Scaling to exactly satisfy the required initial temperature:

            if(tempi.ge.0.d0) then

               tkik=(amass(j)*(v(iii+1)**2+v(iii+2)**2+v(iii+3)**2))/3.d0/boltz*4.d0/3.d0

               v(iii+1)=v(iii+1)*dsqrt(tempi/tkik)
               v(iii+2)=v(iii+2)*dsqrt(tempi/tkik)
               v(iii+3)=v(iii+3)*dsqrt(tempi/tkik)

               tkik=(amass(j)*(v(iii+1)**2+v(iii+2)**2+v(iii+3)**2))/3.d0/boltz*4.d0/3.d0

            end if

!           Satisfying the kinetic energy constraint:

            hhin=amass(j)*(v(iii+1)**2+v(iii+2)**2+v(iii+3)**2)
            s1(j)=dsqrt(4.d0*(1.d0-hhin/(3.d0*tktk)))/erlxt

!           while the chain-thermostat variables are put to zero:

            s2(j)=0.d0
            s3(j)=0.d0
            iii=iii+3
         end do
      end if

!     Initialization for local conservation laws:

      do j=1,natom
         glcl(j)=0.d0
         etlc(j)=(s2(j)**2+s3(j)**2)*erlxt2/2.d0+glcl(j)
      end do
      iconfs=0

!        Initialization for distribution of dynamical variables:

      if(idistr.ne.0) then

         iconfd=0
         dargvs=5.d0/lep
         do jl=1,lep
            argvel(jl)=0.5d0*dargvs+(jl-1)*dargvs
            argvel(-jl)=-argvel(jl)
            distrs(jl)=0.d0
            distrs(-jl)=0.d0
            distrh1(jl)=0.d0
            distrh2(jl)=0.d0
            distrh2(-jl)=0.d0
            distrh3(jl)=0.d0
            distrh3(-jl)=0.d0
         end do

      end if

   ! }}}
  end if
  if (ntt == 10) then
     !   initialize the SINR integrator: {{{

#ifdef MPI
    call sinr_init(natom, nkija, dtx, boltz2, temp0, gammai, sinrtau, &
                   sinrdata, commsander)
#else
    call sinr_init(natom, nkija, dtx, boltz2, temp0, gammai, sinrtau, &
                   sinrdata)
#endif
    if (irest == 1) then
      call sinr_read_vels(v, sinrtau, sinrdata)
    else
      call init_sinr_vels(v, amass, sinrdata)
    endif

#ifdef MPI
    if (numtasks > 1) then
      call sinr_mpi_init(v,numtasks,sinrdata)
    end if
#endif
    call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, &
               xx(l96), xx(l97), xx(l98), xx(l99), qsetup, &
               do_list_update,nstep)
    call iLndt(v, amass, istart, iend, sinrdata)
    call iLvdt(v, f, amass, istart, iend, sinrdata)
    call iLudt(x, v, istart, iend, sinrdata)
    ! }}}
  endif

!------------------------------------------------------------------------------
  ! The main loop for performing dynamics steps: at this point, the
  ! coordinates are a half-step "ahead" of the velocities; the variable
  ! EKMH holds the kinetic energy at these "-1/2" velocities, which are
  ! stored in the array vold.
  260 continue
  onstep = mod(irespa,nrespa) == 0

  ! Deciding if in the current step we are gonna perform a constant pH titration
  ! attempt and/or a constant Redox potential titration attempt
  if (icnstph .gt. 0 .or. (icnste .gt. 0 .and. cpein_specified)) then
    on_cpstep = mod(irespa + nstlim*mdloop, ntcnstph) == 0
  end if
  if (icnste .gt. 0 .and. .not. cpein_specified) then
    on_cestep = mod(irespa + nstlim*mdloop, ntcnste) == 0
  end if

  ! Constant pH setup  {{{
  if ((icnstph /= 0 .or. (icnste /= 0 .and. cpein_specified)) .and. &
      ((rem /= 0 .and. mdloop > 0) .or. rem == 0)) then
    if (ntnb == 1) then ! rebuild pairlist
      call cnstphupdatepairs(x)
    end if
    if (on_cpstep) then
      if (icnstph .eq. 1 .or. (icnste .eq. 1 .and. cpein_specified)) then
        call cnstphbeginstep(xx(l190))
      else
        call cnstph_explicitmd(xx, ix, ih, ipairs, x, winv, amass, f, v, &
                               vold, xr, xc, conp, skip, nsp, tma, erstop, &
                               qsetup, do_list_update,rem,remd_types,replica_indexes)
      end if
    end if
  end if
  ! }}}

  ! If constant pH and constant Redox potential titration attempts are performed in the
  ! same step for implicit solvent, then we need to compute gb_pot_ene here and finish
  ! constant pH
  if (on_cpstep .and. icnstph == 1 .and. &
      on_cestep .and. icnste == 1 .and. .not. cpein_specified) then
    call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
               xx(l98), xx(l99), qsetup, do_list_update,nstep)

    call cnstphendstep(xx(l190), xx(l15), ener%pot%dvdl, temp0, solvph, solve)
    if (master) call cnstphwrite(rem,remd_types,replica_indexes)
  end if

  ! Constant Redox potential setup {{{
  if (icnste /= 0 .and. .not. cpein_specified .and. &
      ((rem /= 0 .and. mdloop > 0) .or. rem == 0)) then
    if (on_cestep) then
      if (icnste .eq. 1) then
        call cnstebeginstep(xx(l190))
      else
        call cnste_explicitmd(xx, ix, ih, ipairs, x, winv, amass, f, v, &
                              vold, xr, xc, conp, skip, nsp, tma, erstop, &
                              qsetup, do_list_update,rem,remd_types,replica_indexes)
      end if
    end if
  end if
  ! }}}

  ! EVB reactive flux: {{{
  ! driver for coordinating backward and
  ! forward propagation as well as for enforcing stopping criteria
#if defined(MPI)
   if (ievb .ne. 0 .and. trim(adjustl(evb_dyn)) == "react_flux") then
     REQUIRE( ipimd.eq.0 .or. ipimd.eq.NMPIMD )
     call react_flux(x, v, f, winv, tempi*factt, dt5, dtx, nr, nstep, nstlim)
   endif
#endif
   ! }}}

!------------------------------------------------------------------------------
  ! Step 1a: do some setup for pressure calculations: {{{
  if (ntp > 0 .and. iamoeba == 0 .and. ipimd==0) then
    ener%cmt(1:3) = 0.d0
    xr(1:nr3) = x(1:nr3)

    ! Calculate, for each molecule, the center of mass, kinetic energy
    ! of the center of mass, and the coordinates of the molecule
    ! relative to that center of mass.
    call timer_start(TIME_EKCMR)
    call ekcmr(nspm, nsp, tma, ener%cmt, xr, v, amass, istart, iend)
#ifdef MPI
    call trace_mpi('mpi_allreduce', 3, 'MPI_DOUBLE_PRECISION', mpi_sum)
# ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE, ener%cmt, 3, MPI_DOUBLE_PRECISION, &
                       mpi_sum, commsander, ierr)
# else
    call mpi_allreduce(ener%cmt, mpitmp, 3, MPI_DOUBLE_PRECISION, mpi_sum, &
                       commsander, ierr)
    ener%cmt(1:3) = mpitmp(1:3)
# endif
#endif
    call timer_stop(TIME_EKCMR)
  end if

  ! If we're using the MC barostat, go ahead and do the trial move now
  if (ntp > 0 .and. barostat == 2 .and. mod(total_nstep+1, mcbarint) == 0) &
    call mcbar_trial(xx, ix, ih, ipairs, x, xc, f, ener%vir, xx(l96), &
                     xx(l97), xx(l98), xx(l99), qsetup, do_list_update, &
                     nstep, nsp, amass)
  ! }}}

!------------------------------------------------------------------------------
  ! Step 1b: prepare to get the forces for the system's current coordinates {{{
  npbstep = nstep
  iprint = 0
  if (nstep == 0 .or. nstep+1 == nstlim) iprint = 1
  if (rem .eq. 0 .or. mdloop .gt. 0) nfe_real_mdstep = .True.
  if (sebomd_obj%do_sebomd) then

    ! Write down atomic charges and density matrix if needed
    sebomd_obj%iflagch = 0
    if (sebomd_obj%ntwc .ne. 0) then
      if (mod(nstep+1,sebomd_obj%ntwc) == 0) sebomd_obj%iflagch = 1
    end if
    sebomd_obj%iflagbo = 0
    if (sebomd_obj%ntwb .ne. 0) then
      if (mod(nstep+1, sebomd_obj%ntwb) == 0) sebomd_obj%iflagbo = 1
    end if
  end if

#ifdef MPI
  ! Set do_mbar for the force contributions
  if (ifmbar /= 0) then
    do_mbar = .false.
    if (mod(nstep+1, bar_intervall) == 0) do_mbar = .true.
  end if
#endif
  ! }}}

!------------------------------------------------------------------------------
  ! Step 1b': actually get the forces for the system's current coordinates {{{

  if (ipimd == NMPIMD .or. ipimd == CMD) then
  ! Contingency for Normal Mode Path Integral MD or Centroid MD {{{
    call trans_pos_nmode_to_cart(x, cartpos)
    call force(xx, ix, ih, ipairs, cartpos, f, ener, ener%vir, xx(l96), &
               xx(l97), xx(l98), xx(l99), qsetup, do_list_update,nstep)
#if defined(MPI) && defined(LES)
    if (ievb == 1 .and. i_qi > 0) then
      call evb_umb(f, cartpos, real_mass, natom, istart3, iend3)
      if (i_qi == 2) then
        call qi_corrf_les(cartpos, real_mass)
      end if
      evb_nrg(1) = evb_frc%evb_nrg
      evb_nrg(2) = evb_vel0%evb_nrg
      if (nbias > 0) evb_nrg(3) = sum(evb_bias%nrg_bias(:))
    end if
#endif
    call trans_frc_cart_to_nmode(f)
#if defined(MPI) && defined(LES)
    if (ievb .ne. 0 .and. i_qi == 0) then
      call evb_umb (f, x, real_mass, natom, istart3, iend3)
      evb_nrg(1) = evb_frc%evb_nrg
      evb_nrg(2) = evb_vel0%evb_nrg
      if (nbias > 0) evb_nrg(3) = sum(evb_bias%nrg_bias(:))
    end if
#endif /* MPI and LES */
  ! }}}
  else
!------------------------------------------------------------------------------
    ! Thermodynamic Integration decomposition {{{
    if (idecomp > 0 .and. ntpr > 0) then
      decpr = .false.
      if (mod(nstep+1, ntpr) == 0) decpr = .true.
    end if
    ! }}}

!------------------------------------------------------------------------------
    ! This(!) is where the force() routine mainly gets called:

    ! determine mass center
#ifdef MPI
    if (sanderrank.eq.0) then
#endif
    if (nstep == 0) then
      write(6,*) 'Determining mass center'
      call ani_get_mass_center(ani_struct%ani_mass_center, ani_struct%ani_coordinates, &
                               ani_struct%ani_mass, ani_nml%ani_natom)
      write(6,*) 'Mass center: ', ani_struct%ani_mass_center
    end if
#ifdef MPI

    end if
#endif

    call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
               xx(l98), xx(l99), qsetup, do_list_update,nstep)

    !   EVB umbrellas:  {{{
#if defined(MPI)
    if (ievb .ne. 0) then
#  ifdef LES
      call evb_umb_primitive(f, x, real_mass, natom, istart, iend)
#  else
      call evb_umb_primitive(f, x, amass, natom, istart, iend)
#  endif /* LES */
      evb_nrg(1) = evb_frc%evb_nrg
      evb_nrg(2) = evb_vel0%evb_nrg
      if (nbias > 0) evb_nrg(3) = sum(evb_bias%nrg_bias(:))
    endif
#endif
    ! }}}
  endif

!------------------------------------------------------------------------------
  ! Semi-Empirical Born-Oppenheimer MD {{{
  if (sebomd_obj%do_sebomd) then

    ! Compute the hessian matrix if necessary
    is_remainder0 = .false.
    if (sebomd_obj%ntwh /= 0) is_remainder0 = mod(nstep+1, sebomd_obj%ntwh) == 0
    if ( is_remainder0 ) then
      ! don't output atomic charges and bond orders
      sebomd_obj%iflagch_old = sebomd_obj%iflagch
      sebomd_obj%iflagch = 0
      sebomd_obj%iflagbo_old = sebomd_obj%iflagbo
      sebomd_obj%iflagbo = 0
      call sebomd_gradient_write(f, 3*natom)
      call sebomd_hessian_compute(xx, ix, ih, ipairs, x, ener, qsetup, &
                                  do_list_update, nstep)
      sebomd_obj%iflagch = sebomd_obj%iflagch_old
      sebomd_obj%iflagbo = sebomd_obj%iflagbo_old
    endif
  endif
  ! }}}

!------------------------------------------------------------------------------
  ! distribute the forces:  (dac: what does f_or stand for?)
  f_or(1:nr3) = f(1:nr3)

#ifdef MPI
  call xdist(f_or, xx(lfrctmp), natom)
#endif

  if (abfqmmm_param%abfqmmm == 1) then
    abfqmmm_param%f2(1:nr3) = f_or(1:nr3)
    call abfqmmm_combine_forces()
#ifdef MPI
    call mpi_bcast(abfqmmm_param%f, 3*natom, mpi_double_precision, 0, &
                   commsander, ierr)
#endif
    f_or(1:nr3) = abfqmmm_param%f(1:nr3)
    f(1:nr3) = abfqmmm_param%f(1:nr3)
  end if

!------------------------------------------------------------------------------
  ! Constant pH transition evaluation for GB CpHMD (not explicit CpHMD)
  if ((icnstph == 1 .or. (icnste == 1 .and. cpein_specified)) .and. on_cpstep .and. &
      .not. on_cestep) then
    call cnstphendstep(xx(l190), xx(l15), ener%pot%dvdl, temp0, solvph, solve)
    if (master) call cnstphwrite(rem,remd_types,replica_indexes)
  end if

  ! Constant Redox potential transition evaluation for GB CEMD (not explicit CEMD)
  if (icnste == 1 .and. .not. cpein_specified .and. on_cestep) then
    call cnsteendstep(xx(l190), xx(l15), ener%pot%dvdl, temp0, solve)
    if (master) call cnstewrite(rem,remd_types,replica_indexes)
  end if

!------------------------------------------------------------------------------
  ! PLUMED force added
  if (plumed == 1) then
#     include "Plumed_force.inc"
  end if

#ifdef MPI
!------------------------------------------------------------------------------
  ! If softcore potentials are used, collect their dvdl contributions: {{{
  if (ifsc .ne. 0) then
    call mpi_reduce(sc_dvdl, sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, commsander, ierr)
    sc_dvdl=0.0d0 ! zero for next step
    call mpi_reduce(sc_dvdl_ee, sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, commsander, ierr)
    sc_dvdl_ee = 0.0d0 ! zero for next step
    call mpi_reduce(sc_ener, sc_ener_tmp, ti_ene_cnt, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, commsander, ierr)
    sc_ener(1:ti_ene_cnt) = sc_ener_tmp(1:ti_ene_cnt)
  end if
  if (ifsc == 2) then

    ! If this is a perturb to nothing run, scale forces and calculate dvdl
    call sc_nomix_frc(f,nr3,ener)
    if (numtasks > 1) then
      call mpi_bcast(f, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call mpi_bcast(ener, state_rec_len, MPI_DOUBLE_PRECISION, 0, &
                     commsander, ierr)
    end if
  end if
  ! }}}

!------------------------------------------------------------------------------
  ! Multi-state Bennet Acceptance Ratio upkeep
  if (ifmbar .ne. 0 .and. do_mbar) call bar_collect_cont()

!------------------------------------------------------------------------------
  ! Free energies using thermodynamic integration: {{{
  if (icfe .ne. 0) then

    ! First, partners exchange forces, energies, and the virial:
    if (master) then
      partner = ieor(masterrank, 1)
      call mpi_sendrecv(f, nr3, MPI_DOUBLE_PRECISION, partner, 5, frcti, &
                        nr3 + 3*extra_atoms, MPI_DOUBLE_PRECISION, partner, &
                        5, commmaster, ist, ierr )
      call mpi_sendrecv(ener, state_rec_len, MPI_DOUBLE_PRECISION, partner, &
                        5, ecopy, state_rec_len, MPI_DOUBLE_PRECISION, &
                        partner, 5, commmaster, ist, ierr)

      ! Exchange sc-dvdl contributions between masters:
      call mpi_sendrecv(sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                        sc_tot_dvdl_partner, 1, MPI_DOUBLE_PRECISION, &
                        partner, 5, commmaster, ist, ierr)
      call mpi_sendrecv(sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                        sc_tot_dvdl_partner_ee, 1, MPI_DOUBLE_PRECISION, &
                        partner, 5, commmaster, ist, ierr )

      ! Collect statistics for free energy calculations
      if (onstep) then
        if (masterrank == 0) then
          if (klambda == 1) then
            edvdl = edvdl - ener + ecopy
            edvdl_r = edvdl_r - ener + ecopy
          else
            clfac = klambda*(1.d0 - clambda)**(klambda-1)
            edvdl = edvdl - (ener - ecopy)*clfac
            edvdl_r = edvdl_r - (ener - ecopy)*clfac
          end if
        else
          if (klambda == 1) then
            edvdl = edvdl + ener - ecopy
            edvdl_r = edvdl_r + ener - ecopy
          else
            clfac = klambda*(1.d0 - clambda)**(klambda-1)
            edvdl = edvdl + (ener - ecopy)*clfac
            edvdl_r = edvdl_r + (ener - ecopy)*clfac
          end if
        end if

        ! This includes the sc-dvdl contribution into the vdw-part
        ! and potential energy parts of the dvdl-statistics
        if (ifsc == 1) call adj_dvdl_stat(edvdl, edvdl_r)
      end if

      ! Do energy collection for MBAR FEP runs
      if (ifmbar .ne. 0 .and. do_mbar) &
        call calc_mbar_energies(ener%pot%tot, ecopy%pot%tot)
      if (masterrank == 0) then
        call mix_frcti(frcti, ecopy, f, ener, nr3, clambda, klambda)
      else
        call mix_frcti(f, ener, frcti, ecopy, nr3, clambda, klambda)
      endif
    endif

    if (numtasks > 1) then
      call mpi_bcast(f, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call mpi_bcast(ener, state_rec_len, MPI_DOUBLE_PRECISION, 0, &
                     commsander, ierr)
    end if
  end if
  ! End contingency for free energies by Thermodynamic Integration }}}
#endif /* MPI */

#ifdef EMIL
  ! Call the EMIL absolute free energy calculation.
  if (emil_do_calc .gt. 0) &
    call emil_step(natom, nstep, 1.0 / (temp0*2*boltz2), &
                   xx(lcrd), f, v, ener%pot, ener%pot, ener%box)
#endif

!------------------------------------------------------------------------------
  ! Reset quantities depending on TEMP0 and TAUTP  {{{
  ! (which may have been changed by MODWT during FORCE call).
  ekinp0 = fac(2)*temp0

  target_ekin_update_nstep = tautp / dt
  update_bussi_target_kin_energy_on_current_step = (ntt==11 .and. target_ekin_update_nstep /=0 .and. &
      mod(total_nstep, target_ekin_update_nstep)==0)
  update_kin_energy_on_current_step = ntt==1 .or. onstep .or. update_bussi_target_kin_energy_on_current_step


#ifdef LES
  ! TEMP0LES may have changed too
  ekinles0=0.d0
  ekins0=0.d0
  if (temp0les >= 0.d0) then
    ekinles0 = fac(3)*temp0les
    ekin0 = ekinp0 + ekinles0
  else
    ekins0 = fac(3)*temp0
    ekin0 = fac(1)*temp0
  end if
  target_ekin_les = ekinles0
#else
  ekins0 = fac(3)*temp0
  ekin0 = fac(1)*temp0
#endif /* LES */
  
  target_ekin = ekin0
  
  if (ntt == 1 .or. ntt == 11) dttp = dt/tautp
  ! }}}

!------------------------------------------------------------------------------
  ! Pressure coupling for PIMD:  {{{
  if (ntp > 0 .and. ipimd > 0) then
    REQUIRE(ipimd == NMPIMD)
    centvir = 0.0
#ifdef LES
    do iatom=istart,iend
      if (cnum(iatom) == 0 .or. cnum(iatom) == 1) then
        centvir = centvir - x(3*iatom-2)*f(3*iatom-2)
        centvir = centvir - x(3*iatom-1)*f(3*iatom-1)
        centvir = centvir - x(3*iatom)*f(3*iatom)
      end if
    end do
#else
    if (mybeadid == 1) then
      do iatom = istart, iend
        centvir = centvir-x(3*iatom-2)*f(3*iatom-2)
        centvir = centvir-x(3*iatom-1)*f(3*iatom-1)
        centvir = centvir-x(3*iatom  )*f(3*iatom)
      end do
    end if
#endif /* LES */

    if (iamoeba == 1) then
      atomvir = sum(ener%vir(1:3))
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
      call mpi_allreduce(MPI_IN_PLACE, centvir, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      call mpi_allreduce(MPI_IN_PLACE, atomvir, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
#  else
      call mpi_allreduce(centvir, mpitmp, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      centvir = mpitmp(1)
      tmp = 0.0
      call mpi_allreduce(atomvir, tmp, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      atomvir=tmp
#  endif /* USE_MPI_IN_PLACE */
#endif /* MPI */
    else
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
      call mpi_allreduce(MPI_IN_PLACE, centvir, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      call mpi_allreduce(MPI_IN_PLACE, bnd_vir, 9, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      call mpi_allreduce(MPI_IN_PLACE, e14vir, 9, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
#    ifndef LES
      if (master) &
        call mpi_allreduce(MPI_IN_PLACE, atvir, 9, MPI_DOUBLE_PRECISION, &
                           mpi_sum, commmaster, ierr)
#    endif /* LES */
#  else
      call mpi_allreduce(centvir, tmp, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      centvir = tmp
      tmpvir = 0.0
      call mpi_allreduce(bnd_vir, tmpvir, 9, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      bnd_vir = tmpvir
#    ifndef LES
      if (master) then
        tmpvir = 0.0
        call mpi_allreduce(e14vir, tmpvir, 9, MPI_DOUBLE_PRECISION, &
                           mpi_sum, commmaster, ierr)
        e14vir = tmpvir
        tmpvir = 0.0
        call mpi_allreduce(atvir, tmpvir, 9, MPI_DOUBLE_PRECISION, &
                           mpi_sum, commmaster, ierr)
        atvir=tmpvir
      endif
#    else
      tmpvir = 0.0
      call mpi_allreduce(e14vir, tmpvir, 9, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      e14vir = tmpvir
#    endif /* LES */
#  endif /* USE_MPI_IN_PLACE */
      call mpi_bcast(atvir, 9, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call mpi_bcast(e14vir, 9, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
#endif /* MPI */
      atomvir = 0.0
      atomvir = atomvir + atvir(1,1) + bnd_vir(1,1) + e14vir(1,1)
      atomvir = atomvir + atvir(2,2) + bnd_vir(2,2) + e14vir(2,2)
      atomvir = atomvir + atvir(3,3) + bnd_vir(3,3) + e14vir(3,3)
    end if
    pressure = (Nkt*3.0 - centvir - (atomvir-Eimp_virial)) / (3.0*volume)
    f_lnv_p = (pressure - pres0/pconv) * volume * 3.0
  end if
  ! }}}

  ! Constant pressure conditions: {{{
  if (ntp > 0) then
    ener%volume = volume
    ener%density = tmass / (0.602204d0*volume)
    if (iamoeba == 0 .and. ipimd == 0) then
      ener%cmt(4) = 0.d0
      ener%vir(4) = 0.d0
      ener%pres(4) = 0.d0
      do m = 1,3
        ener%cmt(m)  = ener%cmt(m)*0.5d0
        ener%cmt(4)  = ener%cmt(4) + ener%cmt(m)
        ener%vir(4)  = ener%vir(4) + ener%vir(m)
        ener%pres(m) = (pconv + pconv) * (ener%cmt(m) - ener%vir(m)) / volume
        ener%pres(4) = ener%pres(4) + ener%pres(m)
      end do
      ener%pres(4) = ener%pres(4) / 3.d0

      ! Constant surface tension output:
      if (csurften > 0) then
        if (csurften == 1) then

          ! Surface tension in the x direction
          ener%surface_ten = &
            box(1) * (ener%pres(1) - 0.5d0 * &
                      (ener%pres(2) + ener%pres(3))) / (ninterface * ten_conv)
        else if (csurften .eq. 2) then

          ! Surface tension in the y direction
          ener%surface_ten = &
            box(2) * (ener%pres(2) - 0.5d0 * &
                      (ener%pres(1) + ener%pres(3))) / (ninterface * ten_conv)
        else

          ! Surface tension in the z direction
          ener%surface_ten = &
            box(3) * (ener%pres(3) - 0.5d0 * &
                      (ener%pres(1) + ener%pres(2))) / (ninterface * ten_conv)
        end if
      end if
    end if
  end if
  ! End contingency for constant pressure conditions }}}

!------------------------------------------------------------------------------
#ifdef MPI
  ! Replica Exchange Molecular Dynamics: if rem /= 0 and mdloop == 0, this is
  ! the first sander call and we don't want to actually do any MD or change the
  ! initial coordinates.  Exit here since we only wanted to get the potential
  ! energy for the first subrem exchange probability calc.
  if (rem /= 0 .and. mdloop == 0) then
#  ifdef VERBOSE_REMD
    if (master) &
      write (6,'(a,i3)') '| REMD: Exiting runmd after getting initial energies &
                          &for replica', repnum
#  endif /* VERBOSE_REMD */

    ! This diverts all the way to the end of the runmd loop
    goto 480
  endif
  ! End contingency for the first sander call in REMD
  ! (rem is not 0, and mdloop == 0)

  ! REB Do adaptive QMMM
  if ( qmmm_nml%vsolv > 1 ) then
    ! Mix forces for adaptive QM/MM and calculate adaptive energy if requested.
    ! Note: nstep is zero during first call; this is the energy/force
    ! calculation with the starting geometry / velocities.
    call adaptive_qmmm(nstep, natom, x, f, ener%pot%tot, ntpr, ntwx, xx, ix, &
                       ih, ipairs, qsetup, do_list_update, corrected_energy, &
                       aqmmm_flag)
  endif
#endif /* MPI */
  ! }}}

!------------------------------------------------------------------------------
  ! Step 1c: do Andersen randomization of velocities: {{{
  !          (ntt == 2).  Assign new random velocities every Vrand steps.
  resetvelo = .false.
  if (vrand .ne. 0 .and. ntt == 2) then
    resetvelo = mod((nstep+1), vrand) == 0
  end if
  if (resetvelo) then

    ! DAN ROE: Why are only the masters doing this?  Even if the velocities
    ! are broadcast to the child processes, the wont the different # of random
    ! calls put the randomg num generators out of sync, or do we not care?
    if (master) then
      write (6,'(a,i8)') 'Setting new random velocities at step ', nstep + 1
      call setvel(nr, v, winv, temp0*factt, iscale, scalm)

#ifdef MPI /* SOFT CORE */
      ! Make sure all common atoms have the same v (that of V0) in TI runs:
      if (icfe .ne. 0 .and. ifsc .ne. 0) then
        call sc_sync_x(v, nr3)
      end if
#endif /* MPI */

#ifdef LES
      ! newvel call is fixed for the dual target temperatures
      if (temp0les >= 0.d0 .and. temp0 .ne. temp0les) then
        vscalt = sqrt(temp0les/temp0)
        do j = 1, natom
          if (cnum(j) > 0) then
            i3 = 3*(j-1)
            v(i3+1) = v(i3+1) * vscalt
            v(i3+2) = v(i3+2) * vscalt
            v(i3+3) = v(i3+3) * vscalt
          end if
        end do
      end if
#endif /* LES */
      if (ibelly > 0) then
         call bellyf(nr, ix(ibellygp), v)
#ifdef MPI /* SOFT CORE */
         if ( ifsc /= 0 ) then
            call bellyf_softcore(v,vold,istart,iend)
         end if
#endif
      end if
    end if
#ifdef MPI
    call trace_mpi('mpi_bcast', 3*natom, 'MPI_DOUBLE_PRECISION', 0)
    call mpi_bcast(v, 3*natom, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
#endif /* MPI */

    ! At this point in the code, the velocities lag the positions by
    ! half a timestep.  If we intend for the velocities to be drawn
    ! from a Maxwell distribution at the timepoint where the positions
    ! and velocities are synchronized, we must correct the new
    ! velocities by backing them up half a step using the current forces.
    ! Note that this fix only works for Newtonian dynamics.
    if (gammai == 0.d0 .and. (ipimd .ne. NMPIMD .or. ipimd .ne. CMD)) then
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
  ! End contingency for velocity reset in Andersen temperature coupling }}}

!------------------------------------------------------------------------------
  ! Step 2: Do the velocity update: {{{
  ! Step 2a: apply quenched MD if needed.  This is useful in NEB>0
  call timer_start(TIME_VERLET)
  if (vv == 1) call quench(f, v)

  ! Car-Parrinello on dipoles: note that the (small?) kinetic energy
  ! of the dipoles is included in the epol energy.
  if (induced > 0 .and. indmeth == 3) call cp_dips(natom,xx(lpol),xx,dt)

!------------------------------------------------------------------------------
  if (ischeme == 1) then
  ! leap-frog middle scheme {{{
  ! the 1st step for updating p
     i3 = 3*(istart - 1)
     do j = istart, iend
       wfac = winv(j) * dtx
       v(i3+1) = v(i3+1) + f(i3+1)*wfac
       v(i3+2) = v(i3+2) + f(i3+2)*wfac
       v(i3+3) = v(i3+3) + f(i3+3)*wfac
       i3 = i3+3
     end do

     ! for constrained MD
     if (ntc /= 1) then
       ! @TODO what if tsgld?
       qspatial = .false.
       ! RATTLE-V, correct velocities
       call rattlev(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
       winv,conp,skip,x,v,nitp,belly,ix(iifstwt),ix(noshake), qspatial)

       ! use SETTLE to deal with water model
       call quick3v(x, v, ix(iifstwr), natom, nres, ix(i02))
     end if
     !}}}
!------------------------------------------------------------------------------
  else if (ntt == 9) then
  ! The generalized canonical-isokinetic algorithm {{{

    ! Propagation of the physical and virtual velocities as well as chain
    ! variables related to the generalized canonical-isokinetic ensemble
    ! using a specified set of atoms for each processor
    iii = 3*(istart-1)
    do j = istart, iend

      ! The contribution of thermostat forces (first half step):
      do kija = 1, nkija
        glcl(j) = glcl(j) + (s2(j) + s3(j) - (s1(j)**2)*s2(j)*erlxt2)*sinsh4
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        if (s3(j) == 0.d0) then
          s2(j) = s2(j) + (s1(j)**2 - erlixt2)*sinsh4
        else
          s2(j) = s2(j) * dexp(-s3(j)*sinsh4) + &
                  (s1(j)**2 - erlixt2)*(1.d0 - dexp(-s3(j)*sinsh4))/s3(j)
        end if
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        hkin = amass(j)*(v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2) + &
               3.d0 * tktk * erlxt2 * (s1(j)**2) / 4.d0
        erlst2 = erlxt2 * 3.d0 * tktk / (hkin * 2.d0)
        esi = dexp(-s2(j) * sinsh2)
        esim = dsqrt((s1(j)**2) * erlst2 * (esi**2 - 1.d0)/2.d0 + 1.d0)
        v(iii+1) = v(iii+1) / esim
        v(iii+2) = v(iii+2) / esim
        v(iii+3) = v(iii+3) / esim
        s1(j) = s1(j) * esi / esim
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        if (s3(j) == 0.d0) then
          s2(j) = s2(j) + (s1(j)**2 - erlixt2)*sinsh4
        else
          s2(j) = s2(j) * dexp(-s3(j)*sinsh4) + &
                  (s1(j)**2 - erlixt2)*(1.d0 - dexp(-s3(j)*sinsh4))/s3(j)
        end if
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        glcl(j) = glcl(j) + (s2(j) + s3(j) - (s1(j)**2)*s2(j)*erlxt2)*sinsh4
      end do

      ! The contribution of physical forces (full step)
      hkin = amass(j) * (v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2) + &
             3.d0 * tktk * erlxt2 * (s1(j)**2) / 4.d0
      aaa = f(iii+1)*v(iii+1) + f(iii+2)*v(iii+2) + f(iii+3)*v(iii+3)
      bbb = (f(iii+1)**2 + f(iii+2)**2 + f(iii+3)**2)/amass(j)
      aaa = aaa / hkin
      bbb = bbb / hkin
      tt = dtx * dsqrt(bbb)
      st = aaa * (dcosh(tt) - 1.d0)/bbb + dsinh(tt)/dsqrt(bbb)
      st1 = dcosh(tt) + aaa*dsinh(tt)/dsqrt(bbb)
      s1(j) = s1(j) / st1
      v(iii+1) = (v(iii+1) + f(iii+1)/amass(j)*st) / st1
      v(iii+2) = (v(iii+2) + f(iii+2)/amass(j)*st) / st1
      v(iii+3) = (v(iii+3) + f(iii+3)/amass(j)*st) / st1

      ! The contribution of thermostat forces (second half step)
      do kija = 1, nkija
        glcl(j) = glcl(j) + (s2(j) + s3(j) - (s1(j)**2)*s2(j)*erlxt2)*sinsh4
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        if (s3(j) == 0.d0) then
          s2(j) = s2(j) + (s1(j)**2 - erlixt2)*sinsh4
        else
          s2(j) = s2(j) * dexp(-s3(j) * sinsh4) + &
                  (s1(j)**2 - erlixt2)*(1.d0 - dexp(-s3(j)*sinsh4))/s3(j)
        end if
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        hkin = amass(j) * (v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2) + &
               3.d0 * tktk * erlxt2 * (s1(j)**2) / 4.d0
        erlst2 = erlxt2 * 3.d0 * tktk / (hkin * 2.d0)
        esi = dexp(-s2(j) * sinsh2)
        esim = dsqrt((s1(j)**2) * erlst2 * (esi**2 - 1.d0)/2.d0 + 1.d0)
        v(iii+1) = v(iii+1) / esim
        v(iii+2) = v(iii+2) / esim
        v(iii+3) = v(iii+3) / esim
        s1(j) = s1(j) * esi / esim
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        if (s3(j) == 0.d0) then
          s2(j) = s2(j) + (s1(j)**2 - erlixt2)*sinsh4
        else
          s2(j) = s2(j) * dexp(-s3(j) * sinsh4) + &
                  (s1(j)**2 - erlixt2)*(1.d0 - dexp(-s3(j)*sinsh4))/s3(j)
        end if
        s3(j) = s3(j) + (s2(j)**2 - erlixt2)*sinsh8
        glcl(j) = glcl(j) + (s2(j) + s3(j) - (s1(j)**2)*s2(j)*erlxt2)*sinsh4
      end do
      iii = iii + 3
    end do

    ! Calculation of some quantities using parallel processing
    iii = 3*(istart - 1)
    do j = istart, iend

      ! Measuring the maximal kinetic constraints deviations
      tkkk = (amass(j) * (v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2) + &
              3.d0 * tktk * erlxt2 * (s1(j)**2) / 4.d0) / (3.d0 * boltz)
      if (dabs(tkkk - temp0) > davalev) then
        davalev = dabs(tkkk - temp0)
      end if

      ! Measuring the conservation laws fluctuations:
      etlci = (s2(j)**2 + s3(j)**2)*erlxt2/2.d0 + glcl(j)
      clfs = clfs + (etlci - etlc(j))**2
      iii = iii+3
    end do

#ifdef MPI
    ! Extracting the quantities from each processor
    ! and collecting them in the root
    call MPI_Gather(davalev, 1, MPI_DOUBLE_PRECISION, sdavalev, 1, &
                    MPI_DOUBLE_PRECISION, 0, commsander, ierr);

    ! Summing up over all processors and sending the result to the master root
    call mpi_reduce(clfs, sclfs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                    commsander, ierr)
    if (master) then
      davalevs = 0.d0
      do inum = 1, numtasks
        if (sdavalev(inum) > davalevs) then
          davalevs = sdavalev(inum)
        end if
      end do
    end if
#endif /* MPI */
    iconfs = iconfs + 1

    ! Writing the results only once using the master processor:
#ifdef MPI
    if (master) then
#endif /* MPI */
      if ((nstep+1) == ntpr*((nstep+1)/ntpr)) then
#ifdef MPI
        write(6,256) '| Max deviations in the kinetic constraints:', &
          davalevs / temp0 * 100.d0, ' %'
        write(6,266) '| Mean fluctuations in the conservation laws:', &
          dsqrt(sclfs / (iconfs*natom))
#else
        write(6,256) '| Max deviations in the kinetic constraints:', &
          davalev / temp0 * 100.d0, ' %'
        write(6,266) '| Mean fluctuations in the conservation laws:', &
          dsqrt(clfs / (iconfs*natom))
#endif /* MPI */
  256   format(a44,e12.3,a2)
  266   format(a45,e12.3)
      end if
#ifdef MPI
    end if
#endif /* MPI */

    ! Calculation of velocity and thermostat distribution functions
    if (idistr .ne. 0) then
      if ((nstep+1) == idistr * ((nstep+1)/idistr)) then
        txtr = 0.d0
        iii = 3*(istart - 1)
        do j = istart, iend
          txtr = txtr + amass(j)*(v(iii+1)**2 + v(iii+2)**2 + v(iii+3)**2)
          iii = iii+3
        end do
#ifdef MPI

        ! Summing up the partial terms over all processors and
        ! sending the result to all
        call mpi_allreduce(txtr, stxtr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           commsander, ierr)
        txtr = stxtr / (natom * 3.d0)
#else
        txtr = txtr / (natom * 3.d0)
#endif
        iii = 3*(istart - 1)
        do j = istart, iend

          ! Distribution of atomic velocities:
          kves = idint(dabs(dsqrt(amass(j))*v(iii+1))/(dsqrt(txtr)*dargvs)) + 1
          if (v(iii+1) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrs(kves) = distrs(kves) + 1.d0
          end if
          kves = idint(dabs(dsqrt(amass(j))*v(iii+2))/(dsqrt(txtr)*dargvs)) + 1
          if (v(iii+2) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrs(kves) = distrs(kves) + 1.d0
          end if
          kves = idint(dabs(dsqrt(amass(j))*v(iii+3))/(dsqrt(txtr)*dargvs)) + 1
          if (v(iii+3) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrs(kves) = distrs(kves) + 1.d0
          end if

          ! And their thermostat counterparts
          kves = idint(s1(j) * erlxt / dargvs) + 1
          if (kves <= lep) then
            distrh1(kves) = distrh1(kves) + 1.d0
          end if

          ! Include chain variables:
          kves = idint(dabs(s2(j)*erlxt)/dargvs) + 1
          if (s2(j) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrh2(kves) = distrh2(kves) + 1.d0
          end if
          kves = idint(dabs(s3(j)*erlxt)/dargvs) + 1
          if (s3(j) < 0.d0) then
            kves = -kves
          end if
          if (iabs(kves) <= lep) then
            distrh3(kves) = distrh3(kves) + 1.d0
          end if
          iii = iii+3
        end do
        iconfd = iconfd + 1
      end if
    end if

#ifdef MPI
    ! Summing up over all processors and sending the results to the master root
    call mpi_reduce(distrs, sdistrs, 2*lep+1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, commsander, ierr)
    call mpi_reduce(distrh1, sdistrh1, lep, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, commsander, ierr)
    call mpi_reduce(distrh2, sdistrh2, 2*lep+1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, commsander, ierr)
    call mpi_reduce(distrh3, sdistrh3, 2*lep+1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, commsander, ierr)

    ! Writing the results only once using the master processor
    if (master) then
#endif /* MPI */

    ! Writing the velocity distribution functions:
    if (iconfd .ne. 0 .and. idistr .ne. 0) then
      if ((nstep+1) == ntpr*((nstep+1)/ntpr)) then
        open(26, file='velocds.dat')
        do jl = -lep, lep
          if (jl .ne. 0) then
#ifdef MPI
            distr = sdistrs(jl) / (iconfd * natom * dargvs * 3.d0)
#else
            distr = distrs(jl) / (iconfd * natom * dargvs * 3.d0)
#endif /* MPI */
            distrmaxwel = ((1.d0 / (2.d0*pif1))**(1.d0/2.d0)) * &
                           dexp(-0.5d0 * argvel(jl)**2)
            write(26, 604) argvel(jl), distr, distrmaxwel
  604       format(1x, f10.5, 2(2x, f10.5))
          end if
        end do
        close(26)
        ! Writing the accompanied distribution function:
        open(27, file='velocds1.dat')
        do jl = 1, lep
#ifdef MPI
          distr = sdistrh1(jl) / (iconfd * natom * dargvs)
#else
          distr = distrh1(jl) / (iconfd * natom * dargvs)
#endif /* MPI */
          write(27, 605) argvel(jl), distr
  605     format(1x, f10.5, 2x, f10.5)
        end do
        close(27)
        ! Writing the chain distributions:
        open(27, file='velocd23.dat')
        do jl = -lep, lep
          if (jl .ne. 0) then
#ifdef MPI
            distr = sdistrh2(jl) / (iconfd * natom * dargvs)
            distra = sdistrh3(jl) / (iconfd * natom * dargvs)
#else
            distr = distrh2(jl) / (iconfd * natom * dargvs)
            distra = distrh3(jl) / (iconfd * natom * dargvs)
#endif /* MPI */
            write(27, 604) argvel(jl), distr, distra
          end if
        end do
        close(27)
      end if
    end if
#ifdef MPI
    end if
    ! End of contingency for writing out velocity distributions,
    ! which occurs only on the master process.  The presence of
    ! this contingency keeps us at two levels of indentation.
#endif /* MPI */

    ! Writing the atomic velocities and chain-thermostat variables
    ! every ntwr steps
    if ((nstep+1) == ntwr*((nstep+1)/ntwr) .or. nstep+1 == nstlim) then
#ifdef MPI
      sv = 0.d0
      s1a = 0.d0
      s2a = 0.d0
      s3a = 0.d0
      iii = 3*(istart-1)
      do j = istart, iend
        sv(iii+1) = v(iii+1)
        sv(iii+2) = v(iii+2)
        sv(iii+3) = v(iii+3)
        iii = iii + 3
        s1a(j) = s1(j)
        s2a(j) = s2(j)
        s3a(j) = s3(j)
      end do
      svs = 0.d0
      call mpi_reduce(sv, svs, 3*natom, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                      commsander, ierr)
      s1s = 0.d0
      s2s = 0.d0
      s3s = 0.d0
      call mpi_reduce(s1a, s1s, natom, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                      commsander, ierr)
      call mpi_reduce(s2a, s2s, natom, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                      commsander, ierr)
      call mpi_reduce(s3a, s3s, natom, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                      commsander, ierr)
      if (master) then
        open(578, file='vfreez.rst')
        iii = 0
        do j = 1, natom
          write(578, *) svs(iii+1), svs(iii+2), svs(iii+3)
          iii = iii+3
        end do
        close(578)
        open(579, file='tfreez.rst')
        do j = 1, natom
          write(579, *) s1s(j), s2s(j), s3s(j)
        end do
        close(579)
      end if
#else
      open(578, file='vfreez.rst')
      iii = 0
      do j = 1, natom
        write(578, *) v(iii+1), v(iii+2), v(iii+3)
        iii = iii+3
      end do
      close(578)
      open(579, file='tfreez.rst')
      do j = 1, natom
        write(579,*) s1(j),s2(j),s3(j)
      end do
      close(579)
#endif /* MPI */
    end if
    ! }}}

!------------------------------------------------------------------------------
  else if (ntt == 10) then
  ! Stochastic Isokinetic Nose-Hoover RESPA (SINR) integrator) {{{
    call iLvdt(v, f, amass, istart, iend, sinrdata)
    call iLndt(v, amass, istart, iend, sinrdata)
    call iLndt(v, amass, istart, iend, sinrdata)
    call iLvdt(v, f, amass, istart, iend, sinrdata)
    call iLudt(x, v, istart, iend, sinrdata)
    if (mod(nstep+1, ntpr*nrespa) == 0) then
      call sinr_temp(v, amass, istart3, iend3, sinrdata)
    end if
    if (mod(nstep+1, ntwr) == 0 .or. nstep+1 == nstlim) then
      call sinr_write_vels(v, nstep, istart3, iend3, sinrtau, sinrdata)
    end if
    ! }}}

  else if (ntt == 4) then
  ! Nose'-Hoover thermostat (1st step). {{{
    Ekin2_tot = 0.d0
    i3 = 3*(istart - 1)
    do j = istart, iend
      wfac = dtx / amass(j)
      do idim = 1, 3
#ifdef LES
        if (ntp > 0 .and. ipimd == NMPIMD .and. &
            (cnum(j) .eq. 0 .or. cnum(j) .eq. 1)) then
#else
        if (ntp > 0 .and. ipimd == NMPIMD .and. mybeadid == 1) then
#endif /* LES */
          exp1 = exp(-dt5*thermo(idim,j)%v(1) - dt5*v_lnv*c2_lnv)
          Ekin2_tot = Ekin2_tot + amass(j) * v(i3+idim) * v(i3+idim)
        else
          exp1 = exp(-dt5 * thermo(idim,j)%v(1))
        end if
        exp2 = exp1 * exp1
        vold(i3+idim) = v(i3+idim)
        v(i3+idim) = v(i3+idim) * exp2 + f(i3+idim) * wfac * exp1
      end do
      i3 = i3 + 3
    end do
    if (ntp > 0 .and. ipimd == NMPIMD) then
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
      call mpi_allreduce(MPI_IN_PLACE, Ekin2_tot, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
#  else
      call mpi_allreduce(Ekin2_tot, mpitmp, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      Ekin2_tot = mpitmp(1)
#  endif /* USE_MPI_IN_PLACE */
#endif /* MPI */
      f_lnv_v = Ekin2_tot*(c2_lnv - 1)
      tmp = exp(-dt5 * thermo_lnv%v(1))
      v_lnv = tmp*(tmp*v_lnv + dtx*(f_lnv_v + f_lnv_p)/mass_lnv)
    end if
    call Thermostat_integrate_1(nchain, thermo, nthermo, dtx, ntp)
    ! }}}

  else if (ntt > 4 .and. ntt <= 8) then
    !   Other thermostats?: {{{
    Ekin2_tot = 0.d0
    i3 = 3*(istart - 1)
    do j = istart, iend
      wfac = dtx / amass(j)
      do idim = 1, 3
#ifdef LES
        if (ntp > 0 .and. ipimd == NMPIMD .and. &
            (cnum(j) == 0 .or. cnum(j) == 1)) then
#else
        if (ntp > 0 .and. ipimd == NMPIMD .and. mybeadid == 1) then
#endif
          Ekin2_tot = Ekin2_tot + amass(j)*v(i3+idim)*v(i3+idim)
          exp1 = exp(-dt5 * v_lnv * c2_lnv)
        else
          exp1 = 1.d0
        end if
        exp2 = exp1 * exp1
        vold(i3+idim) = v(i3+idim)
        v(i3+idim) = v(i3+idim) * exp2
        f(i3+idim) = f(i3+idim) * exp1
      end do
      i3 = i3 + 3
    end do
    if (ntp > 0 .and. ipimd == NMPIMD) then
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
      call mpi_allreduce(MPI_IN_PLACE, Ekin2_tot, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
#  else
      call mpi_allreduce(Ekin2_tot, mpitmp, 1, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commworld, ierr)
      Ekin2_tot = mpitmp(1)
#  endif /* USE_MPI_IN_PLACE */
#endif /* MPI */
      f_lnv_v = Ekin2_tot*(c2_lnv - 1)
    end if
    if (abfqmmm_param%abfqmmm == 1) then
#ifdef MPI
      call xdist(v, xx(lfrctmp), natom)
      call xdist(f, xx(lfrctmp), natom)
#endif /* MPI */
      abfqmmm_param%v(1:nr3+iscale) = v(1:nr3+iscale)
      abfqmmm_param%f(1:nr3+iscale) = f(1:nr3+iscale)
    end if
    call Adaptive_Thermostat_integrate(nchain, thermo, nthermo, dtx, ntp, 1)
    if (abfqmmm_param%abfqmmm == 1) then
      v(1:nr3+iscale) = abfqmmm_param%v(1:nr3+iscale)
#ifdef MPI
      call xdist(v, xx(lfrctmp), natom)
#endif
      abfqmmm_param%v(1:nr3+iscale) = v(1:nr3+iscale)
    end if
    ! }}}

!------------------------------------------------------------------------------
  else if (gammai == 0.d0) then
    ! Newtonian dynamics: possible with guiding force: {{{
    if (isgld > 0) &
      call sgmdw(natom, istart, iend,  dtx, ener, amass, winv, x, f, v)
    i3 = 3*(istart - 1)
    do j = istart, iend
      wfac = winv(j) * dtx
      v(i3+1) = v(i3+1) + f(i3+1)*wfac
      v(i3+2) = v(i3+2) + f(i3+2)*wfac
      v(i3+3) = v(i3+3) + f(i3+3)*wfac
      i3 = i3 + 3
    end do
    ! }}}

!------------------------------------------------------------------------------
  else if (isgld > 0) then
    call sgldw(natom, istart, iend, dtx, temp0, ener, amass, winv, &
               x, f, v)

!------------------------------------------------------------------------------
  else
    !   Langevin dynamics: {{{
    ! gamma_ln .ne. 0, which also implies ntt=3 (see mdread.f)
    ! Simple model for Langevin dynamics, basically taken from
    ! Loncharich, Brooks and Pastor, Biopolymers 32:523-535 (1992),
    ! Eq. 11. (Note that the first term on the rhs of Eq. 11b
    ! should not be there.) Update Langevin parameters, since temp0
    ! might have changed:
    sdfac = sqrt(4.d0 * gammai * boltz2 * temp0 / dtx)
#ifdef LES
    sdfacles = sqrt(4.d0 * gammai * boltz2 * temp0les / dtx)
#endif /* LES */

#ifdef MPI /* SOFT CORE */
    if (ifsc == 1) then
      call sc_lngdyn(winv, amass, v, f, sdfac, c_explic, c_implic, &
                     istart, iend, nr, dtx)
      if (ibelly > 0) then
         call bellyf_softcore(v,vold,istart,iend)
         call bellyf(nr,ix(ibellygp),v)
      end if
    else
#endif
    if (no_ntt3_sync == 1) then

      ! We don't worry about synchronizing the random number stream
      ! across processors.
      iskip_start = 0
      iskip_end = 0
    else

      ! In order to generate the same sequence of pseudorandom numbers that
      ! you would using a single processor you have to go through the atoms
      ! in order: skip those that have are being used on other processors
      iskip_start = 3*(istart-1)
      iskip_end = 3*(nr-iend)
#ifndef LES
      ! Always sync random number stream for PIMD
      ! (AWG: not sure if this is required)
      if (ipimd > 0) then
        iskip_start = iskip_start + 3*nr*(mybeadid-1)
        iskip_end = iskip_end + 3*nr*(nbead-mybeadid)
      end if
#endif
    endif
      do j = 1, iskip_start

        ! Skip some random numbers
        call gauss(0.d0, 1.d0, fln)
      end do

      ! Do Langevin step
      i3 = 3*(istart - 1)
      do j = istart, iend
        wfac = winv(j) * dtx
        aamass = amass(j)
#ifdef LES
        if (temp0les >= 0 .and. temp0 .ne. temp0les .and. cnum(j) .ne. 0) then
          rsd = sdfacles * sqrt(aamass)
        else
          rsd = sdfac * sqrt(aamass)
        endif
#else
        rsd = sdfac*sqrt(aamass)
#endif /* LES */
        call gauss(0.d0, rsd, fln)
        v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
        call gauss(0.d0, rsd, fln)
        v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
        call gauss(0.d0, rsd, fln)
        v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
        i3 = i3 + 3
      end do
      if (ibelly > 0) then
         call bellyf(nr,ix(ibellygp),v)
#ifdef MPI /* SOFT CORE */
         if ( ifsc /= 0 ) then
            call bellyf_softcore(v,vold,istart,iend)
         end if
#endif
      end if
      do j = 1, iskip_end

        ! Skip some random numbers
        call gauss(0.d0, 1.d0, fln)
      end do
#ifdef MPI /* SOFT CORE */
    end if ! for (ifsc==1) call sc_lngdyn
#endif
  end if  ! ( gammai == 0.d0 )
  !  }}}
  ! End case switch for various thermostats

!------------------------------------------------------------------------------
  ! Update EMAP rigid domains
  if (temap) call emap_move()

  ! Consider vlimit {{{
  if (vlim .and. ipimd == 0) then
    vmax = 0.0d0
    do i = istart3, iend3
      vmax = max(vmax, abs(v(i)))
      v(i) = sign(min(abs(v(i)), vlimit), v(i))
    end do

    ! Only violations on the master node are actually reported
    ! to avoid both MPI communication and non-master writes.
    if (vmax > vlimit) then
      if (master) &
        write(6,'(a,i6,a,f10.4)') 'vlimit exceeded for step ', nstep, &
              '; vmax = ', vmax
    end if
  end if
  ! }}}

  !  Simple Newtonian dynamics on the "extra" variables  (why?)
  do im = 1, iscale
    v(nr3+im) = (v(nr3+im) + f(nr3+im)*dtx/scalm)
  end do

!------------------------------------------------------------------------------
  ! Dump forces: {{{
  ! We do the force dump here if requested, since the "old"
  ! positions are about to be dumped into the force array...
  if (master) then

    ! Test if forces are to be written this step
    ifdump = .false.
    if (ntwf > 0) ifdump = (mod(total_nstep+1,ntwf) == 0)
    if (ntwf == -1 .and. ntwx /= 0) ifdump = mod(total_nstep+1,ntwx) == 0
    ! ABF QM/MM has a combined coordinate and force file:
    ! never write forces by themselves for this type of dynamics
    if (abfqmmm_param%abfqmmm == 1) ifdump = .false.
#ifdef MPI
    ! For adaptive QM/MM, only the master does a dump.
    if (qmmm_nml%vsolv > 1) then
      if (nodeid .ne. 0) ifdump = .false.
    end if
    if (ifdump) call xdist(f, xx(lfrctmp), natom)
#endif
    ! Force archive:
    if (ifdump) then
#ifdef MPI
      ! Write out current replica#, exchange#, step#, and mytargettemp
      ! If mdloop==0 this is a normal md run (since REMD never calls corpac
      !  when mdloop==0) and we don't want the REMD header.
      if (mdloop > 0 .and. loutfm) then
        if (trxsgld) then
          write (MDFRC_UNIT,'(a,4(1x,i8))') "RXSGLD ", repnum, mdloop, &
                total_nstep+1, stagid
        else
          write (MDFRC_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, &
                mdloop, total_nstep+1, my_remd_data%mytargettemp
        end if
      end if
#endif
      ! ipimd forces will probably not be right if some type of
      ! transformation is necessary. This is from the vel dump code -- keep
      ! it here as a holder in case somebody wants to fix it.
      call corpac(f,1,nrx,MDFRC_UNIT,loutfm)
    end if
  else

    ! Slave processes need to participate in force distribution,
    ! but do not print forces themselves.
    ifdump = .false.
    if (ntwf > 0) ifdump = (mod(total_nstep+1,ntwf) == 0)
    if (ntwf == -1 .and. ntwx /= 0) ifdump = mod(total_nstep+1,ntwx) == 0
    if (abfqmmm_param%abfqmmm == 1) ifdump = .false.
#ifdef MPI
    if (ifdump) call xdist(f, xx(lfrctmp), natom)
#endif /* MPI */
  end if
  ! End branch for master and slave processes in force printing.
  ! The "old" positions can now be dumped into the force array. }}}
  ! }}}
  ! Determine if we need to save the on-step coordinates
  if (isynctraj.ne.0) then
    ! Decision to write trajectory coords
    if (ntwx > 0 .and. (mod(total_nstep+1, ntwx) == 0)) then !onstepcrd(1:3*natom) = x(1:3*natom)
      do i3 = istart3, iend3
        onstepcrd(i3) = x(i3)
      end do
    end if
  endif

!------------------------------------------------------------------------------
  ! Step 3: update the positions, putting the "old" positions into F: {{{
  if (ischeme == 1) then
    ! middle scheme {{{
    ! the step for updating x-T-x
    do i3 = istart3, iend3
       f(i3) = x(i3)
       x(i3) = x(i3) + v(i3)*dt5
    end do
    ! for random number, controlled by ig, see ntt part
    if (no_ntt3_sync == 1) then
      iskip_start = 0
      iskip_end = 0
    else
      iskip_start = 3*(istart-1)
      iskip_end = 3*(nr-iend)
#ifndef LES
      ! Always sync random number stream for PIMD
      ! (AWG: not sure if this is required)
      if (ipimd > 0) then
        iskip_start = iskip_start + 3*nr*(mybeadid-1)
        iskip_end = iskip_end + 3*nr*(nbead-mybeadid)
      end if
#endif
    endif
    call thermostat_step(v, winv, dtx, istart, iend, iskip_start,iskip_end)
    do i3 = istart3, iend3
       x(i3) = x(i3) + v(i3)*dt5
    end do
    !}}}
  else
    ! default AMBER scheme {{{
#ifdef LES
  if (ntp > 0 .and. ipimd == NMPIMD) then
    aa = exp(dt5 * v_lnv)
    arg2 = v_lnv * dt5 * v_lnv * dt5
    poly = 1.0d0 + arg2*(e2 + arg2*(e4 + arg2*(e6 + arg2*e8)))
  endif
  i3 = 3*(istart-1)
  do j = istart, iend
    if (ntp > 0 .and. ipimd == NMPIMD .and. &
        (cnum(j) == 0 .or. cnum(j) == 1)) then
      do idim = 1, 3
        f(i3+idim) = x(i3+idim)
        x(i3+idim) = aa*(x(i3+idim)*aa + v(i3+idim)*poly*dtx)
      enddo
    else
      do idim = 1, 3
        f(i3+idim) = x(i3+idim)
        x(i3+idim) = x(i3+idim)+v(i3+idim)*dtx
      enddo
    endif
    i3 = i3 + 3
  enddo
#else
  if (ntp > 0 .and. ipimd == NMPIMD .and. mybeadid == 1) then
    aa = exp(dt5 * v_lnv)
    arg2 = v_lnv * dt5 * v_lnv * dt5
    poly = 1.0d0 + arg2*(e2 + arg2*(e4 + arg2*(e6 + arg2*e8)))
    do i3 = istart3, iend3
      f(i3) = x(i3)
      x(i3) = aa*(x(i3)*aa + v(i3)*poly*dtx)
    end do
  else if (ntt == 10) then
    do i3 = istart3, iend3
      f(i3) = x(i3)
    end do
  else
    do i3 = istart3, iend3
      f(i3) = x(i3)
      x(i3) = x(i3) + v(i3)*dtx
    end do
  end if
#endif /* LES */
   !}}}
  end if
  !}}}
  ! End case branching for updating positions

!------------------------------------------------------------------------------
  !Nose'-Hoover thermostat (2nd step).{{{
  if (ntt == 4) then
    call Thermostat_integrate_2(nchain, thermo, nthermo, dtx, ntp)
    E_nhc = Thermostat_hamiltonian(nchain, thermo, nthermo)
  else if (ntt >= 4 .and. ntt <= 8) then
    if (abfqmmm_param%abfqmmm == 1) then
#ifdef MPI
      call xdist(v, xx(lfrctmp), natom)
#endif
      abfqmmm_param%v(1:nr3+iscale)=v(1:nr3+iscale)
    end if
    call Adaptive_Thermostat_integrate(nchain,thermo,nthermo,dtx,ntp,2)
    if (abfqmmm_param%abfqmmm == 1) then
      v(1:nr3+iscale)=abfqmmm_param%v(1:nr3+iscale)
#ifdef MPI
      call xdist(v, xx(lfrctmp), natom)
#endif
      abfqmmm_param%v(1:nr3+iscale)=v(1:nr3+iscale)
    end if
    E_nhc = Adaptive_Thermostat_hamiltonian(nchain,thermo,nthermo)
  end if

  ! position update for the "extra" variables"
  do i = 1,iscale
    f(nr3+i) = x(nr3+i)
    x(nr3+i) = x(nr3+i) + v(nr3+i)*dtx
  end do

  call timer_stop(TIME_VERLET)
  ! }}}

!------------------------------------------------------------------------------
  if (ntc .ne. 1) then
    ! Step 4a: if shake is being used, update the positions {{{
    call timer_start(TIME_SHAKE)
    if (ischeme == 1) xold(istart3:iend3) = x(istart3:iend3)
   qspatial = .false.
    call shake(nrp, nbonh, nbona, 0, ix(iibh), ix(ijbh), ix(ibellygp), &
               winv, conp, skip, f, x, nitp, belly, ix(iifstwt), &
               ix(noshake), qspatial)
    call quick3(f, x, ix(iifstwr), natom, nres, ix(i02))
    if (nitp == 0) then
      erstop = .true.
      goto 480
    end if


    ! Need to synchronize coordinates for linearly scaled atoms after shake
#ifdef MPI
    if (icfe .ne. 0) then
      call timer_barrier( commsander )
      call timer_stop_start(TIME_SHAKE,TIME_DISTCRD)
      if (.not. mpi_orig .and. numtasks > 1) &
        call xdist(x, xx(lfrctmp), natom)

      ! In dual-topology this is done within softcore.f
      if (ifsc .ne. 1) then
        if (master) call mpi_bcast(x, nr3, MPI_DOUBLE_PRECISION, &
                         0, commmaster, ierr)
      else
        if (master) call sc_sync_x(x, nr3)
      end if
      if (numtasks > 1) &
        call mpi_bcast(x, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call timer_stop_start(TIME_DISTCRD, TIME_SHAKE)
    end if
#endif  /* MPI */

!------------------------------------------------------------------------------
    ! Step 4b: Now fix the velocities and calculate KE.
    ! Re-estimate the velocities from differences in positions.
    if (.not. (ipimd == NMPIMD .and. ipimd == CMD .and. mybeadid .ne. 1)) then
      if (ischeme == 0) then
        v(istart3:iend3) = (x(istart3:iend3) - f(istart3:iend3))*dtxinv
      else
        v(istart3:iend3) = v(istart3:iend3) &
          + (x(istart3:iend3) - xold(istart3:iend3))*dtxinv
      end if
    end if
    call timer_stop(TIME_SHAKE)
    ! }}}
  end if
  call timer_start(TIME_VERLET)

  if (ischeme == 1) then
    if (ntc /= 1) then
      qspatial = .false.
      ! RATTLE-V, correct velocities
      call rattlev(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
      winv,conp,skip,x,v,nitp,belly,ix(iifstwt),ix(noshake), qspatial)

      ! use SETTLE to deal with water model
      call quick3v(x, v, ix(iifstwr), natom, nres, ix(i02))
    end if
  end if

  ! NEB: remove velocities but ONLY for the end beads so V doesn't
  ! accumulate if there are high forces
  if (ineb > 0 .and. (mybeadid == 1 .or. mybeadid == neb_nbead)) then
    x(1:3*natom) = f(1:3*natom)
    v(1:3*natom) = 0.d0
  end if

!------------------------------------------------------------------------------
  ! NOTE: The decision to write velocities is made here to decide if we need
  !       to calculate the on-step velocities. The decision to write coordinates
  !       is made as well to determine if this is a combined coords/velocities
  !       file.
  ! Decision to write trajectory coords
  itdump = .false.
  if (ntwx > 0) itdump = (mod(total_nstep+1, ntwx) == 0)
  ! Decision to write velocities
  ivdump = .false.
  if (ntwv > 0) ivdump = (mod(total_nstep+1,ntwv) == 0)
  ! Combined coordinate and velocity file writing
  if (ntwv == -1 .and. itdump) ivdump = .true.
  ! Determine whether we want to save the on-step velocities.
  if ((isynctraj.ne.0).and.ivdump) then
    i3 = (3*(istart-1)) + 1
    do j = istart, iend
      onstepvel(i3  ) = (v(i3  ) + vold(i3  )) / 2.d0
      onstepvel(i3+1) = (v(i3+1) + vold(i3+1)) / 2.d0
      onstepvel(i3+2) = (v(i3+2) + vold(i3+2)) / 2.d0
      i3 = i3 + 3
    enddo
  endif

  if (update_kin_energy_on_current_step) then
    ! Step 4c: get the KE, either for averaging or for Berendsen:  {{{
    eke = 0.d0
    ekph = 0.d0
    ekpbs = 0.d0
#ifdef LES
    ekeles = 0.d0
    ekphles = 0.d0
#endif
    eke_cmd = 0.d0

    if (ischeme == 1) then
    ! LF-Middle: use velocity v(t) for KE calculation {{{
      i3 = 3*(istart-1)
      do j = istart, iend
        aamass = amass(j)
        do m = 1, 3
          i3 = i3 + 1
#ifdef LES
          if (temp0les < 0.d0) then
            eke = eke + aamass*0.25d0*(v(i3) + vold(i3))**2
            ekph = ekph + aamass*v(i3)**2
          else
            if (cnum(j) == 0) then
              eke = eke + aamass*0.25d0*(v(i3) + vold(i3))**2
              ekph = ekph + aamass*v(i3)**2
            else
              ekeles = ekeles + aamass*0.25d0*(v(i3) + vold(i3))**2
              ekphles = ekphles + aamass*v(i3)**2
            end if
          end if
#else
          eke = eke + aamass*0.25d0*(v(i3) + vold(i3))**2
          ! Try pseudo KE from Eq. 4.7b of Pastor, Brooks & Szabo,
          ! Mol. Phys. 65, 1409-1419 (1988):
          ekpbs = ekpbs + aamass*v(i3)*vold(i3)
          ekph = ekph + aamass*v(i3)**2
#endif
        end do
      end do

      !}}}
    else if (gammai == 0.0d0) then
      ! No Langevin forces:  {{{

      i3 = 3*(istart - 1)
      do j = istart, iend
        aamass = amass(j)
        do m = 1, 3
          i3 = i3 + 1
#ifdef LES
          if (temp0les < 0.d0) then
            eke = eke + aamass*0.25d0*(v(i3) + vold(i3))**2
            ekph = ekph + aamass*v(i3)**2
            if (ipimd == CMD .and. (cnum(j) == 0 .or. cnum(j) == 1)) then
              eke_cmd = eke_cmd + aamass*0.25d0*(v(i3)+vold(i3))**2
            endif
          else
            if (cnum(j) == 0) then
              eke = eke + aamass*0.25d0*(v(i3) + vold(i3))**2
              ekph = ekph + aamass*v(i3)**2
            else
              ekeles = ekeles + aamass*0.25d0*(v(i3) + vold(i3))**2
              ekphles = ekphles + aamass*v(i3)**2
            end if
          end if
#else
          eke = eke + aamass*0.25d0*(v(i3) + vold(i3))**2
          if (mybeadid == 1) &
            eke_cmd = eke_cmd + aamass*0.25d0*(v(i3) + vold(i3))**2

          ! Try pseudo KE from Eq. 4.7b of Pastor, Brooks & Szabo,
          ! Mol. Phys. 65, 1409-1419 (1988):
          ekpbs = ekpbs + aamass*v(i3)*vold(i3)
          ekph = ekph + aamass*v(i3)**2
#endif
        end do
      end do
      ! }}}
    else
      ! Langevin integrator KE estimation: {{{
      i3 = 3*(istart - 1)
      do j = istart, iend
        aamass = amass(j)
        do m = 1, 3
          i3 = i3 + 1
#ifdef LES
          if (temp0les < 0.d0) then
            eke = eke + aamass*0.25d0*c_ave*(v(i3) + vold(i3))**2
          else
            if (cnum(j) == 0) then
              eke = eke + aamass*0.25d0*c_ave*(v(i3) + vold(i3))**2
            else
              ekeles = ekeles + aamass*0.25d0*c_ave*(v(i3) + vold(i3))**2
            end if
          end if
#else
          if (ntt == 9) then
            eke = eke + (aamass*v(i3)**2)*4.d0/3.d0
          else
            eke = eke + aamass*0.25d0*c_ave*(v(i3) + vold(i3))**2
          end if
#endif
        end do
      end do
    end if
    ! End branch based on gammai }}}

    ! Sum up the partial kinetic energies: {{{
#ifdef MPI
    if (ipimd == CMD) then
      call mpi_reduce(eke_cmd, tmp_eke_cmd, 1, MPI_DOUBLE_PRECISION, &
                      mpi_sum, 0, commsander, ierr)
      eke_cmd = tmp_eke_cmd
    endif
#  ifdef LES
    if (.not. mpi_orig .and. numtasks > 1) then
      if (temp0les < 0) then
        mpitmp(1) = eke
        mpitmp(2) = ekph
#    ifdef USE_MPI_IN_PLACE
        call mpi_allreduce(MPI_IN_PLACE, mpitmp, 2, MPI_DOUBLE_PRECISION, &
                           mpi_sum, commsander, ierr)
        eke = mpitmp(1)
        ekph = mpitmp(2)
#    else
        call mpi_allreduce(mpitmp, mpitmp(3), 2, MPI_DOUBLE_PRECISION, &
                           mpi_sum, commsander, ierr)
        eke = mpitmp(3)
        ekph = mpitmp(4)
#    endif
      else
        mpitmp(1) = eke
        mpitmp(2) = ekph
        mpitmp(3) = ekeles
        mpitmp(4) = ekphles
#    ifdef USE_MPI_IN_PLACE
        call mpi_allreduce(MPI_IN_PLACE, mpitmp, 4, MPI_DOUBLE_PRECISION, &
                           mpi_sum, commsander, ierr)
        eke = mpitmp(1)
        ekph = mpitmp(2)
        ekeles = mpitmp(3)
        ekphles = mpitmp(4)
#    else
        call mpi_allreduce(mpitmp, mpitmp(5), 4, MPI_DOUBLE_PRECISION, &
                           mpi_sum, commsander, ierr)
        eke = mpitmp(5)
        ekph = mpitmp(6)
        ekeles = mpitmp(7)
        ekphles = mpitmp(8)
#    endif
      endif
    end if
#  else
    if (.not. mpi_orig .and. numtasks > 1) then
      call trace_mpi('mpi_allreduce', 1, 'MPI_DOUBLE_PRECISION', mpi_sum)
      mpitmp(1) = eke
      mpitmp(2) = ekph
      mpitmp(3) = ekpbs
#    ifdef USE_MPI_IN_PLACE
      call mpi_allreduce(MPI_IN_PLACE, mpitmp, 3, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commsander, ierr)
      eke = mpitmp(1)
      ekph = mpitmp(2)
      ekpbs = mpitmp(3)
#    else
      call mpi_allreduce(mpitmp, mpitmp(4), 3, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commsander, ierr)
      eke = mpitmp(4)
      ekph = mpitmp(5)
      ekpbs = mpitmp(6)
#    endif
    end if
#  endif

    ! Calculate Ekin of the softcore part of the system
    if (ifsc .ne. 0) then
      call calc_softcore_ekin(amass, v, vold, istart, iend)
      sc_ener(13) = sc_ener(6) + sc_ener(12)
    end if
#endif

    ! All processors handle the "extra" variables:
    do im = 1, iscale
      eke = eke + scalm*0.25d0*(v(nr3+im) + vold(nr3+im))**2
      ekpbs = ekpbs + scalm*v(nr3+im)*vold(nr3+im)
      ekph = ekph + scalm*v(nr3+im)**2
    end do
    eke = eke * 0.5d0
    ekph = ekph * 0.5d0
    ekpbs = ekpbs * 0.5d0
#ifdef LES
    ekeles = ekeles * 0.5d0
    ekphles = ekphles * 0.5d0
#endif /* LES */
    if (ntt == 1 .or. update_bussi_target_kin_energy_on_current_step) then
      if (ntt == 1) then
#ifdef LES
        if (temp0les < 0.d0) then
          scaltp = sqrt(1.d0 + 2.d0*dttp*(ekin0-eke)/(ekmh+ekph))
        else
          scaltp = sqrt(1.d0 + 2.d0*dttp*(ekinp0-eke)/(ekmh+ekph))
          scaltles = sqrt(1.d0 + 2.d0*dttp*(ekinles0-ekeles)/(ekmhles+ekphles))
        end if
#else
        ! Following is from T.E. Cheatham, III and B.R. Brooks,
        ! Theor. Chem. Acc. 99:279, 1998.
        scaltp = sqrt(1.d0 + 2.d0*dttp*(ekin0-eke)/(ekmh+ekph))
#endif /* LES */
      ! }}}
      else if (update_bussi_target_kin_energy_on_current_step) then
#ifdef LES
        if (temp0les < 0.d0) then
          target_ekin  = resamplekin(eke, ekin0, int(rndf, 4), dttp)
          scaltp = sqrt(target_ekin/eke)
        else
          target_ekin  = resamplekin(eke, ekinp0, int(rndf, 4), dttp)
          scaltp = sqrt(target_ekin/eke)
          
          target_ekin_les  = resamplekin(ekeles, ekinles0, int(rndf, 4), dttp)
          scaltles = sqrt(target_ekin/eke)
        end if
#else
        target_ekin  = resamplekin(eke, ekin0, int(rndf, 4), dttp)
        scaltp = sqrt(target_ekin/eke)
#endif
      endif
#ifdef MPI /* SOFT CORE */
      if (icfe .ne. 0) then
        if (ifsc == 1) then
          if (master) then

            ! Linearly combine the scaling factors from both processes
            ! the combined factor is broadcast to all nodes
            ! the subroutine also correctly scales the softcore atom v's
            call mix_temp_scaling(scaltp, v)
          end if
          call mpi_bcast(scaltp, 1, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        end if
      end if
#endif /* MPI */
      do j = istart, iend
        i3 = (j-1)*3 + 1
#ifdef LES
        if (temp0les > 0.d0 .and. cnum(j) /= 0 ) then
          v(i3) = v(i3)*scaltles
          v(i3+1) = v(i3+1)*scaltles
          v(i3+2) = v(i3+2)*scaltles
        else
          v(i3) = v(i3)*scaltp
          v(i3+1) = v(i3+1)*scaltp
          v(i3+2) = v(i3+2)*scaltp
        end if
#else
        v(i3) = v(i3)*scaltp
        v(i3+1) = v(i3+1)*scaltp
        v(i3+2) = v(i3+2)*scaltp
#endif
      end do
      do im=1,iscale
        v(nr3+im) = v(nr3+im)*scaltp
      end do
    end if
    ! End contingency for Berendsen thermostat (ntt == 1)
  end if
  ! End contingency for Berendsen thermostat or onstep; end of step 4c
  !    (also: end of big section devoted to estimating the kinetic energy) }}}

!------------------------------------------------------------------------------
  ! Step 5: several tasks related to dumping of trajectory information  {{{
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
  total_nstep = nstep + 1
  total_nstlim = nstlim
  if (abfqmmm_param%abfqmmm == 1) then
    total_nstep = abfqmmm_param%qmstep
  end if

#ifdef MPI
  if (rem .ne. 0) then
    total_nstep = (mdloop - 1) * nstlim + nstep + 1
    total_nstlim = nstlim * numexchg
  endif
#endif

  ! NOTE: The decisions to write coords (itdump)/velocities (ivdump) have
  !       been moved to right before the KE calculation in order to have
  !       the option to write on-step velocities (i.e. determine if we need
  !       to calculate the on-step velocities).
  !       Similarly, the decision to write forces (ifdump) has been moved to
  !       right before step 3 (positions update) because the old positions
  !       get dumped into the force array at that point. However, we retain
  !       the additional logic to write forces here as well for QM/MM forces.

  ! Decision to write forces
  ifdump = .false.
  if (ntwf > 0) ifdump = (mod(total_nstep,ntwf) == 0)

  ! Decision to write a restart file, or the final restart file
  ixdump = .false.
  if ( ntwr /= 0) ixdump = mod(total_nstep, ntwr ) == 0
  if (total_nstep >= total_nstlim) ixdump = .true.
  if ( nstlim == 0 ) ixdump = .false.

  ! Decision to remove velocity of the system center of mass
  ivscm  = .false.
  if (nscm > 0) ivscm = mod(total_nstep,nscm) == 0

#ifdef MPI
  ! Adaptive QM/MM via multisander: all groups have identical
  ! coords and velocities only master of first group needs to dump results
  ! We have to leave the dump values for all threads in the group, though
  ! since for dumping the coords, these are broadcast within the group
  ! (see call to xdist() below)
  if (qmmm_nml%vsolv > 1) then
    if (nodeid .ne. 0) then
      ixdump = .false.
      itdump = .false.
      ivdump = .false.
    end if
  end if
#endif

#ifdef RISMSANDER
  ! Write RISM files this step?
  irismdump = .false.
  if (rismprm%rism == 1) then
    if (rismprm%ntwrism > 0) then
      irismdump = (mod(nstep+1,rismprm%ntwrism) == 0)
      if (nstep + 1 >= nstlim) irismdump = .true.
    end if
  end if
#endif

#ifdef MPI
!------------------------------------------------------------------------------
  ! Distribute the coordinates, dipoles, and velocities as necessary
  call timer_barrier(commsander)
  call timer_stop_start(TIME_VERLET,TIME_DISTCRD)
  if (.not. mpi_orig .and. numtasks > 1) then
    !write(6,'(a)') 'DEBUG: Calling XDIST for coordinates.'
    call xdist(x, xx(lfrctmp), natom)
    if (isynctraj.ne.0.and.itdump) call xdist(onstepcrd, xx(lfrctmp), natom) 
  endif

  ! DAC/knut change: force the coordinates to be the same on both masters.
  ! For certain compilers, addition may not be strictly commutative, so
  ! the forces on group 0 may be different by roundoff from the forces on
  ! group 1.  This can lead to divergent trajectories.  The interval at
  ! which they are resynchronized is hard-wired here to 20, which seems to
  ! work fine in our tests.
  !
  ! jwk change: coordinates are synchronized when shake is enabled above
  if (icfe .ne. 0 .and. mod(nstep+1, 20) == 0 .and. ntc == 1) then

    ! In dual-topology this is done within softcore.f
    if (ifsc .ne. 1) then
      if (master) &
        call mpi_bcast(x, nr3, MPI_DOUBLE_PRECISION, 0, commmaster, ierr)
    else
      if (master) then

        ! First, check if coordinates have desynced,
        ! then do the same for velocities
        call sc_compare(x, nr3, 'CRD')
        if (numtasks == 1) call sc_compare(v,nr3,'VEL')

        ! Resync coordinates and velocities
        call sc_sync_x(x,nr3)
      end if
    end if
    if (numtasks > 1) &
      call mpi_bcast(x, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
  end if
  call timer_stop(TIME_DISTCRD)
#endif  /* MPI */

!------------------------------------------------------------------------------
  ! Fix lone pair positions
  if (numextra > 0) call local_to_global(x, xx, ix)

#ifdef MPI
  if (.not. mpi_orig .and. numtasks > 1) then
    call timer_start(TIME_DISTCRD)

!------------------------------------------------------------------------------
    ! Here we provide every processor a full copy of the velocities
    ! for removal of center of mass motion, or for archiving.
    ! (Note: this is actually over-kill: for example, only the master
    ! node really needs the velocities for archiving.  But the extra
    ! overhead of doing it this way is probably small in most cases.)
    if (ivdump .or. ivscm .or. ixdump) then
      call xdist(v, xx(lfrctmp), natom)
      call xdist(onstepvel, xx(lfrctmp), natom)
    endif
    if (ixdump .and. (induced > 0 .and. indmeth == 3)) then
      call xdist(xx(ldipvel), xx(lfrctmp), natom)
      call xdist(xx(linddip), xx(lfrctmp), natom)
    end if
    call timer_stop(TIME_DISTCRD)
  end if
  call timer_start(TIME_VERLET)
  ! This is the end of major broadcasting operations for parallel runs
#endif  /* MPI */
  ! }}}

!------------------------------------------------------------------------------
  ! Step 6: zero COM velocity if requested.  {{{
  ! This is used for preventing
  ! the "block of ice flying thru space" phenomenon (caused by Berendsen
  ! thermocoupling, or by coarse Ewald approximations).  It also
  ! prevents accumulations of rotational momentum in vacuum simulations.
  if (ivscm) then
    if ( nsnb /= 0 ) then
      if (mod(nstep,nsnb) == 0) ntnb = 1
    end if
    if (ifbox == 0) then
      if ((ischeme==0.and.is_langevin) .or. (ischeme>0.and.ipimd==0)) then

        ! Get current center of the system
        call get_position(nr, x, vcmx, vcmy, vcmz, sysrange, 0)
#ifdef MPI /* SOFT CORE */
        if (ifsc == 1) call sc_mix_position(vcmx, vcmy, vcmz, clambda)
#endif /* MPI */
        ! Center the system to the original center
        call re_position(nr, ntr, x, xc, vcmx, vcmy, vcmz, sysx, sysy, &
                         sysz, sysrange, mv_flag, 0)

      elseif (ischeme>0 .and. ipimd>0) then
#ifdef MPI
#ifndef LES
         if (sanderrank.eq.0) then
            call MPI_Allgather(x   ,3*natom,MPI_DOUBLE_PRECISION,&
                               xall,3*natom,MPI_DOUBLE_PRECISION,&
                               commmaster,ierr) 
         endif
         
         call MPI_Bcast(xall,3*natom*nbead,MPI_DOUBLE_PRECISION,&
            0, commsander,ierr)
         x_centroid(1:3,1:natom)=0.0
         do i = 1, nbead
            do j = istart, iend
               x_centroid(1:3,j)=x_centroid(1:3,j)+xall(1:3,j,i)
            end do
         end do
         x_centroid(1:3,1:natom)=x_centroid(1:3,1:natom)/dble(nbead)
         ! position of mass center
         vcmx = 0.d0
         vcmy = 0.d0
         vcmz = 0.d0
         do i = 1, natom
           aamass = amass(i)
           vcmx = vcmx + aamass*x_centroid(1,i)
           vcmy = vcmy + aamass*x_centroid(2,i)
           vcmz = vcmz + aamass*x_centroid(3,i)
         end do
         vcmx = vcmx * tmassinv * dble(nbead) 
         vcmy = vcmy * tmassinv * dble(nbead)
         vcmz = vcmz * tmassinv * dble(nbead)
         call re_position(nr, ntr, x, xc, vcmx, vcmy, vcmz, sysx, sysy, &
                          sysz, sysrange, mv_flag, 0)
#endif
#endif
      else
        ! Non-periodic simulation: remove both translation and rotation.
        ! Back the coords up 1/2 step, so that they correspond to the
        ! velocities; temporarily store in the F() array:
        f(1:nr3) = x(1:nr3) - v(1:nr3)*dt5

        ! Now compute the com motion, remove it, and recompute (just
        ! to check that it is really gone.....)
        call cenmas(nr, f, v, amass, ekcm, xcm, vcm, acm, ekrot, ocm, 4)
        call stopcm(nr, f, v, xcm, vcm, ocm, .true.)
        call cenmas(nr, f, v, amass, ekcm, xcm, vcm, acm, ekrot, ocm, 4)
      end if
    else

      ! This must be checked, as the branch is based on ifbox
      if (ischeme==0 .and. (.not. is_langevin)) then

        ! Periodic simulation: just remove the translational velocity:
        vcmx = 0.d0
        vcmy = 0.d0
        vcmz = 0.d0
        j = 1
        do i = 1, 3*natom,3
          aamass = amass(j)
          vcmx = vcmx + aamass*v(i)
          vcmy = vcmy + aamass*v(i+1)
          vcmz = vcmz + aamass*v(i+2)
          j = j + 1
        end do
        vcmx = vcmx * tmassinv
        vcmy = vcmy * tmassinv
        vcmz = vcmz * tmassinv
        vel2 = vcmx*vcmx + vcmy*vcmy + vcmz*vcmz

        ! The array onefac is the inverse of the array fac
        atempdrop = 0.5d0 * tmass * vel2 * onefac(1)
        vel = sqrt(vel2)
        if (master) &
          write (6, '(a,f15.6,f9.2,a)') 'check COM velocity, temp: ', vel, &
                atempdrop, '(Removed)'
        do i = 1, 3*natom, 3
          v(i) = v(i) - vcmx
          v(i+1) = v(i+1) - vcmy
          v(i+2) = v(i+2) - vcmz
        end do

#ifdef MPI /* SOFT CORE */
        if (icfe == 1) then
          if (ifsc == 1) then
            if (master) call sc_mix_velocities(v, nr3)
            call mpi_bcast(v, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
          end if
        end if
#endif /* MPI */
      end if
      ! End of contingency for thermostats other than Langevin
    end if
    ! End of contingency for isolated, non-periodic systems (ifbox == 0)
  end if
  ! End of contingency for zeroing system center of mass velocity  }}}


!------------------------------------------------------------------------------
  !  Also zero out the non-moving velocities if a belly is active:
  if (belly) then
     call bellyf(nr, ix(ibellygp), v)
#ifdef MPI /* SOFT CORE */
     if ( ifsc /= 0 ) then
        call bellyf_softcore(v,vold,istart,iend)
     end if
#endif
  end if

!------------------------------------------------------------------------------
  ! Put current velocities into vold

  ! ANI P validation
    !ANI P validation
  !call cenmas(nr, f, v, amass, ekcm, xcm, vcm, acm, ekrot, ocm, 4)
  !ani_struct%ani_mom_p(2) = 0.0d0
  !ani_struct%ani_mom_p(3) = 0.0d0
  !do i = 1,natom
  !  aamass = amass(i)
  !  ani_struct%ani_mom_p(1) = ani_struct%ani_mom_p(1) + aamass*v((i-1)*3+1)
  !  ani_struct%ani_mom_p(2) = ani_struct%ani_mom_p(2) + aamass*v((i-1)*3+2)
  !  ani_struct%ani_mom_p(3) = ani_struct%ani_mom_p(3) + aamass*v((i-1)*3+3)
  !end do

  vold(istart3:iend3) = v(istart3:iend3)
  do im = 1, iscale
    vold(nr3+im) = v(nr3+im)
  end do


  
  !write(6,*) 'ANI P validation: ', ani_struct%ani_mom_p(1)

!------------------------------------------------------------------------------
  ! Step 7: scale coordinates if NPT with Berendsen barostat:
  !  First, for ipimd: {{{
  if (ntp > 0 .and. ipimd > 0 .and. barostat == 1) then
    x_lnv_old = x_lnv
    x_lnv = x_lnv_old + v_lnv * dtx
    rmu(1:3) = exp(x_lnv - x_lnv_old)
    box(1:3) = box(1:3) * rmu(1:3)
    volume = box(1) * box(2) * box(3)
    ener%box(1:3) = box(1:3)
    ! Only for NMPIMD in sander.LES
    ! (in sander.MPI volume, pressure and density printed in pimdout)
#ifdef LES
    ener%volume = volume
#else
    ener%volume = 0.
    totener%volume = volume
#endif
    call redo_ucell(rmu)
    call fill_tranvec()
    call ew_pscale(natom,x,amass,nspm,nsp,2)
  end if
  ! }}}
  if (iamoeba == 0 .and. barostat == 1) then
    !  next, for non-amoeba: {{{
    if (ntp == 1) then

      ! Isotropic pressure coupling
      rmu(1) = (1.d0-dtcp*(pres0 - ener%pres(4)))**third
      rmu(2) = rmu(1)
      rmu(3) = rmu(1)
    else if (ntp == 2) then

      ! Anisotropic pressure scaling
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

          ! For surface tension in the z direction
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
      ! End branch for anisotropic pressure scaling

    else

      ! This means ntp = 3, semi-isotropic pressure coupling.
      ! For csurften > 0, constant surface tension.
      ! For csurften < 0, semi-isotropic pressure scaling in
      !      any direction with no constant surface tension.
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

      else if (csurften < 0) then
        ! Semiisotropic pressure scaling in any direction with no constant 
        ! surface tension here
          if (csurften == -1) then        ! For pressure control in x
            rmu(1) = (1.d0 - dtcp * (pres0x - ener%pres(1)))**third
            rmu(2) = 1.0d0
            rmu(3) = 1.0d0

          else if (csurften == -2) then   ! For pressure control in y
            rmu(1) = 1.0d0
            rmu(2) = (1.d0 - dtcp * (pres0y - ener%pres(2)))**third
            rmu(3) = 1.0d0

          !else if (csurften == 3) then   ! For pressure control in z
          else
            rmu(1) = 1.0d0
            rmu(2) = 1.0d0
            rmu(3) = (1.d0 - dtcp * (pres0z - ener%pres(3)))**third
          end if

      end if
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

#ifdef MPI /* SOFT CORE */
      ! If softcore potentials and the dual topology approach are used
      ! C.O.M. scaling has to be changed to account for different masses
      ! of the same molecule in V0 and V1. This is quite inefficient and is
      ! therefore done in a separate routine in softcore.f
      ! only both masters actually do the computation for ifsc==1
      ! the scaled coordinates are then broadcast to the nodes
      if (icfe .ne. 0 .and. ifsc == 1) then
        if (master) then
           if (ibelly > 0) then 
              call sc_belly_pscale(x,amass,nspm,natom,nsp,oldrecip,ucell,ix(ibellygp),bcshift)
           else
              call sc_pscale(x,amass,nspm,nsp,oldrecip,ucell)
           endif
        end if
        if (ibelly > 0) then 
        call mpi_bcast(bcshift,3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
        end if
        call mpi_bcast(x,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
      else
#endif /* MPI */
      if (ibelly > 0) then 
        call belly_pscale(natom,x,amass,nspm,nsp,ix(ibellygp))
      else
      call ew_pscale(natom,x,amass,nspm,nsp,npscal)
      endif
#ifdef MPI /* SOFT CORE */
      end if
#endif /* MPI */
      if (ntr > 0 .and. nrc > 0) &
        call ew_pscale(natom, xc, amass, nspm, nsp, npscal)
    end if
    if (ipimd == NMPIMD .and. ntp > 0) then
      ener%cmt(4) = 0.d0
      ener%vir(4) = 0.d0
      ener%pres(4) = pressure*pconv
    end if
    ! }}}
  else if (barostat == 1) then
    ! finally, for amoeba: {{{
    if (ntp > 0) then
      if (ipimd == 0) then     ! for classical AMOEBA
        ener%cmt(4) = eke      !  for printing in prntmd()
        ener%vir(4) = ener%vir(1) + ener%vir(2) + ener%vir(3)
        ener%pres(4) = (pressure_constant/volume) * &
                       (2.d0*eke - ener%vir(4)) / 3.d0
      elseif (ipimd == NMPIMD) then     ! for NMPIMD AMOEBA
        ener%cmt(4)  = 0.d0
        ener%vir(4)  = 0.d0
        ener%pres(4) = pressure*pconv
      endif
      call AM_RUNMD_scale_cell(natom, ener%pres(4), dt, pres0, taup, x)
      call fill_tranvec()
    end if
    ! }}}
  end if

#ifdef LES
  ener%kin%solt = eke
  ener%kin%solv = ekeles
  ener%kin%tot  = ener%kin%solt + ener%kin%solv
  if (ntt == 1 .and. onstep) then
    if (temp0les < 0) then
      ekmh = max(ekph, fac(1)*10.d0)
    else
      ekmh = max(ekph,fac(2)*10.d0)
      ekmhles = max(ekphles,fac(3)*10.d0)
    endif
  end if

  if (ipimd > 0) then
    ener%kin%solv = equal_part + Epot_deriv  ! "virial" estimate of KE
    ener%tot = ener%kin%solv + ener%pot%tot
  endif
#else
  if (ipimd > 0) then

    ! Use a "virial" estimator for the KE, rather
    ! than one derived from the  bead velocities
    totener%kin%solv = equal_part + Epot_deriv
  else

    ! Pastor, Brooks, Szabo conserved quantity
    ! for harmonic oscillator: Eq. 4.7b of Mol.
    ! Phys. 65:1409-1419, 1988
    ener%kin%solv = ekpbs + ener%pot%tot
  endif
  ener%kin%solt = eke
  ener%kin%tot  = ener%kin%solt
  if (ntt == 1 .and. onstep) ekmh = max(ekph,fac(1)*10.d0)
#endif /* LES */

  ! If velocities were reset, the KE is not accurate; fudge it
  ! here to keep the same total energy as on the previous step.
  ! Note that this only affects printout and averages for Etot
  ! and KE -- it has no effect on the trajectory, or on any
  ! averages of potential energy terms.
  if (resetvelo) ener%kin%tot = etot_save - ener%pot%tot

!------------------------------------------------------------------------------
  ! Total energy is sum of KE + PE:
  if (ipimd > 0) then
    totener%tot = totener%kin%solv + totener%pot%tot
    etot_save   = totener%kin%tot  + totener%pot%tot
    if (ipimd == CMD) then
      etot_cmd = eke_cmd*0.5 + ener%pot%tot
      totener%tot = etot_cmd
      ener%tot = etot_cmd
      ener%kin%tot  = eke_cmd*0.5
      ener%kin%solv = ener%kin%tot
    end if
  else
    ener%tot = ener%kin%tot + ener%pot%tot
    etot_save = ener%tot
  end if

!------------------------------------------------------------------------------
  ! Step 8: update the step counter and the integration time: {{{
  if (abfqmmm_param%abfqmmm .ne. 1) then
    nstep = nstep + 1
    t = t + dt
  end if

  ! For Centroid MD
  if (ipimd == CMD) then
    nstep_cmd = nstep_cmd + 1
    t_cmd = t_cmd + dt
  end if

  ! Full energies are only calculated every nrespa steps.
  ! nvalid is the number of steps where all energies are calculated.
  if (onstep .or. aqmmm_flag > 0) then
    nvalid = nvalid + 1

    ! Update all elements of these sequence types
    enert  = enert + ener
    enert2 = enert2 + (ener*ener)
    if (ischeme > 0) then
       ekhf = ekhf + ekph
       ekhf2 = ekhf2 + ekph*ekph
    endif
#ifdef MPI
    if (ievb .ne. 0) then
      evb_nrg_ave(:) = evb_nrg_ave(:) + evb_nrg(:)
      evb_nrg_rms(:) = evb_nrg_rms(:) + evb_nrg(:)**2
    endif
    if (ifsc .ne. 0) then
      sc_ener_ave(1:ti_ene_cnt) = sc_ener_ave(1:ti_ene_cnt) + &
                                  sc_ener(1:ti_ene_cnt)
      sc_ener_rms(1:ti_ene_cnt) = sc_ener_rms(1:ti_ene_cnt) + &
                                  sc_ener(1:ti_ene_cnt)**2
    end if
#endif /* MPI */
    if (nvalid == 1) etot_start = ener%tot

#ifndef LES
#  ifdef MPI
    if (master .and. (ipimd > 0 .or. ineb > 0)) &
      call mpi_reduce(ener%kin%tot, totener%kin%tot, 1, &
                      MPI_DOUBLE_PRECISION, mpi_sum, 0, commmaster, ierr)
#  endif /* MPI */

!------------------------------------------------------------------------------
    ! Passing of dvdl=dV/dl for Thermodynamic Integration w.r.t. mass
    ! Note that ener(39) (in runmd and mix_frcti) =
    !       = ener(17) = ene(21) (in force). All denote dvdl.
    ! Note, ener() is now historical, MJW Feb 2010
    if (ipimd>0 .and. itimass>0) totener%pot%dvdl = ener%pot%dvdl

    if (ipimd == NMPIMD .and. ntp > 0) then
      totener%pres(4) = pressure * pconv
      totener%density = tmass / (0.602204d0*volume)
    endif
    if (ipimd == CMD) then
      totener%kin%tot = eke_cmd * 0.5d0
      totener%kin%solv = totener%kin%tot
      totener%tot = totener%kin%tot + totener%pot%tot
    endif
    totenert  = totenert + totener
    totenert2 = totenert2 + (totener*totener)
#endif /* LES is not defined */
    kinetic_E_save(2) = kinetic_E_save(1)
    kinetic_E_save(1) = ener%kin%tot
  end if
  ! End contingency to calculate energies when on a reportable step

  ! Added for rbornstat
  ! !FIX: TL - do we need to put in rismnrespa here?
  if (mod(irespa,nrespai) == 0 .or. irespa < 2) nvalidi = nvalidi + 1
  ntnb = 0
  if ( nsnb /= 0 ) then
    if (mod(nstep,nsnb) == 0) ntnb = 1
  end if

  ! Since nstep has been incremented, total_nstep is now equal to
  ! (mdloop-1)*nstlim+nstep for REMD and nstep for MD.
  lout = .false.
  if (ntpr > 0) lout = (mod(total_nstep,ntpr) == 0 .and. onstep)
  irespa = irespa + 1

  ! Reset flags related to Poisson-Boltzmann electrostatics
#ifdef MPI
  if (mytaskid == 0) then
#endif /* MPI */
  if (igb == 10 .or. ipb .ne. 0) then
    if (mod(nstep,npbgrid) == 0 .and. nstep .ne. nstlim) pbgrid = .true.
    pbprint = nstep == nstlim
    if ( ntpr > 0 ) pbprint = pbprint .or. mod(nstep,ntpr) == 0
    if (mod(nstep,nsnbr) == 0 .and. nstep /= nstlim) ntnbr = 1
    if (mod(nstep,nsnba) == 0 .and. nstep /= nstlim) ntnba = 1
  end if
#ifdef MPI
  end if
#endif /* MPI */
  ! }}}

!------------------------------------------------------------------------------
  ! Step 9: output from this step if required: {{{
  !    RISM dumping: {{{
#ifdef RISMSANDER
  ! Some 3D-RISM files require all processes to participate in output
  ! due to the distributed memory.  RISM archive:
  if (rismprm%rism == 1) then
    ! Combined thermodynamics and distribution output.
    ! Execute if we need to do either.
    if (irismdump .or. (rism_calc_type(nstep) == RISM_FULL .and. &
        rismprm%write_thermo == 1 .and. lout)) &
      call rism_solvdist_thermo_calc(irismdump, nstep)
  end if
#endif
   ! }}}
  !    coordinate dumps: {{{
  if (ixdump) then
    if (ipimd .eq. NMPIMD .or. ipimd == CMD) then
#if defined(LES) || defined(MPI)
      call trans_pos_nmode_to_cart(x, cartpos)
      call trans_vel_nmode_to_cart(v, cartvel)
#else
      call trans_pos_nmode_to_cart(cartpos)
      call trans_vel_nmode_to_cart(cartvel)
#endif /* LES or MPI */
    endif
  endif
  ! }}}
  !    some non-standard dumps: {{{
  if (itdump) then
    if (ipimd == NMPIMD .or. ipimd == CMD) then
#if defined(LES) || defined(MPI)
      call trans_pos_nmode_to_cart(x, cartpos)
#else
      call trans_pos_nmode_to_cart(cartpos)
#endif /* LES or MPI */
    endif
    ! Accelerated MD: Flush amdlog file
    if (iamd > 0) then
#ifdef MPI
      if (worldrank == 0) then
#endif /* MPI */
      call write_amd_weights(ntwx,total_nstep)
#ifdef MPI
      end if
#endif /* MPI */
    end if

    ! ScaledMD: Flush scaledMDlog file
    if (scaledMD > 0) then
#ifdef MPI
      if (worldrank == 0) then
#endif /* MPI */
      call write_scaledMD_log(ntwx,total_nstep)
#ifdef MPI
      end if
#endif /* MPI */
    end if
  end if
  ! }}}
  !    velocity dumps: {{{
  if (ivdump) then
    if (ipimd == NMPIMD .or. ipimd == CMD) then
#if defined(LES) || defined(MPI)
      call trans_vel_nmode_to_cart(v, cartvel)
#else
      call trans_vel_nmode_to_cart(cartvel)
#endif /* LES or MPI */
    endif
  endif
  ! }}}

  ! Begin writing output on the master process.
  if (master) then

    !    Restart file writing {{{
    if (ixdump) then

      ! NOTE - This assumes that if numextra > 0, then velocities are
      !        found in the array v...
      if (numextra > 0) call zero_extra_pnts_vec(v, ix)

!------------------------------------------------------------------------------
      if (iwrap == 0) then
        nr = nrp
#ifdef LES
        if (ipimd == NMPIMD .or. ipimd == CMD) then
          call mdwrit(nstep, nr, ntxo, ntb, cartpos, cartvel, t, temp0,solvph,solve)
        else
          call mdwrit(nstep,nr,ntxo,ntb,x,v,t,temp0les,solvph,solve)
        endif
#else
        if (ipimd == NMPIMD .or. ipimd == CMD) then
          call mdwrit(nstep, nr, ntxo, ntb, cartpos, cartvel, t, temp0,solvph,solve)
        else
          call mdwrit(nstep, nr, ntxo, ntb, x, v, t, temp0,solvph,solve)
        endif
#endif /* LES */
!------------------------------------------------------------------------------
      else if (iwrap == 1) then

        ! Use a temporary array to hold coordinates so that the master's
        ! values are always identical to those on all other nodes:
        call get_stack(l_temp, nr3, routine)
        if (.not. rstack_ok) then
          deallocate(r_stack)
          allocate(r_stack(1:lastrst), stat=alloc_ier)
          call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)
        if (ipimd == NMPIMD .or. ipimd == CMD) then
          do iatom = 1, natom
            do m = 1, 3
              r_stack(l_temp+3*(iatom-1)+m-1) = cartpos(m,iatom)
            end do
          end do
        else
          do m = 1, nr3
            r_stack(l_temp+m-1) = x(m)
          end do
        end if
        call wrap_molecules(nspm, nsp, r_stack(l_temp))
        if (ifbox == 2) then
          call wrap_to(nspm, nsp, r_stack(l_temp), box)
        end if
        nr = nrp
#ifdef LES
        call mdwrit(nstep, nr, ntxo, ntb, r_stack(l_temp), v, t, temp0les,solvph,solve)
#else
        call mdwrit(nstep, nr, ntxo, ntb, r_stack(l_temp), v, t, temp0,solvph,solve)
#endif
        call free_stack(l_temp, routine)
!------------------------------------------------------------------------------
      else if (iwrap == 2) then

        ! Wrapping around a pre-determined mask: center the system on
        ! the mask center of mass first, then wrap it normally as it
        ! happens on the iwrap = 1 case.
        call get_stack(l_temp, nr3, routine)
        if (.not. rstack_ok) then
          deallocate(r_stack)
          allocate(r_stack(1:lastrst), stat=alloc_ier)
          call reassign_rstack(routine)
        endif
        REQUIRE(rstack_ok)
        if (ipimd == NMPIMD .or. ipimd == CMD) then
          do iatom = 1, natom
            do m = 1, 3
              r_stack(l_temp+3*(iatom-1)+m-1) = cartpos(m,iatom)
            end do
          end do
        else
          do m = 1, nr3
            r_stack(l_temp+m-1) = x(m)
          end do
        end if
        nr = nrp

        ! Now, wrap the coordinates around the iwrap_mask:
        call iwrap2(n_iwrap_mask_atoms, iwrap_mask_atoms, r_stack(l_temp), &
                    box_center)
#ifdef LES
        call mdwrit(nstep, nr, ntxo, ntb, r_stack(l_temp), v, t, temp0les,solvph,solve)
#else
        call mdwrit(nstep, nr, ntxo, ntb, r_stack(l_temp), v, t, temp0,solvph,solve)
#endif /* LES */
        call free_stack(l_temp, routine)
      end if
      ! End branches for molecule wrapping in periodic boundary conditions

      if (igb == 0 .and. ipb == 0 .and. induced > 0 .and. indmeth == 3) then
        call wrt_dips(xx(linddip), xx(ldipvel), nr, t, title)
      end if
      if ((icnstph .ne. 0 .or. (icnste .ne. 0 .and. cpein_specified)) .and. &
          ((rem .ne. 0 .and. mdloop > 0) .or. rem == 0)) then
        call cnstphwriterestart(total_nstep,chrgdat)
      end if
      if (icnste .ne. 0 .and. .not. cpein_specified .and. &
          ((rem .ne. 0 .and. mdloop > 0) .or. rem == 0)) then
        call cnstewriterestart(total_nstep,chrgdat_e)
      end if
    end if
    ! End decision process for restart file writing (ixdump flag) }}}

    !    Coordinate archive: {{{
    !   for formatted writes and replica exchange, write out a header line.
    if (itdump) then
#ifdef MPI
      ! Write out current replica#, exchange#, step#, and mytargettemp
      ! If mdloop==0 this is a normal md run (since REMD never calls
      ! corpac when mdloop==0) and we don't want the REMD header.
      ! total_nstep is set in step 5.
      if (mdloop > 0 .and. loutfm) then
        if (trxsgld) then
          write (MDCRD_UNIT,'(a,4(1x,i8))') "RXSGLD ", repnum, mdloop, &
                total_nstep, stagid
        else
          write (MDCRD_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, mdloop, &
                total_nstep, my_remd_data%mytargettemp
        end if
      end if
#endif /* MPI */
      if (iwrap == 0) then
        if (ipimd == NMPIMD .or. ipimd == CMD) then
          call corpac(cartpos, 1, nrx, MDCRD_UNIT, loutfm)
        else if (isynctraj.ne.0) then
          call corpac(onstepcrd, 1, nrx, MDCRD_UNIT, loutfm)
        else
          call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
        endif
        if (ntb > 0) call corpac(box, 1, 3, MDCRD_UNIT, loutfm)
      else if (iwrap == 1) then
        call get_stack(l_temp, nr3, routine)
        if (.not. rstack_ok) then
          deallocate(r_stack)
          allocate(r_stack(1:lastrst), stat=alloc_ier)
          call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)
        if (ipimd == NMPIMD .or. ipimd == CMD) then
          do iatom = 1, natom
            do m = 1, 3
              r_stack(l_temp+3*(iatom-1)+m-1) = cartpos(m,iatom)
            end do
          end do
        else if (isynctraj.ne.0) then
          do m = 1, nr3
            r_stack(l_temp+m-1) = onstepcrd(m)
          end do
        else
          do m = 1, nr3
            r_stack(l_temp+m-1) = x(m)
          end do
        endif
        call wrap_molecules(nspm, nsp, r_stack(l_temp))
        if (ifbox == 2) call wrap_to(nspm, nsp, r_stack(l_temp), box)
        call corpac(r_stack(l_temp), 1, nrx, MDCRD_UNIT, loutfm)
        call corpac(box, 1, 3, MDCRD_UNIT, loutfm)
        call free_stack(l_temp, routine)

      else if (iwrap == 2) then

        ! Wrapping around a pre-determined mask: center the system on
        ! the mask center of mass first, then wrap it normally as it
        ! happens on the iwrap = 1 case.
        call get_stack(l_temp, nr3, routine)
        if (.not. rstack_ok) then
          deallocate(r_stack)
          allocate(r_stack(1:lastrst), stat=alloc_ier)
          call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)
        if (ipimd == NMPIMD .or. ipimd == CMD) then
          do iatom = 1, natom
            do m = 1, 3
              r_stack(l_temp+3*(iatom-1)+m-1) = cartpos(m,iatom)
            end do
          end do
        else if (isynctraj.ne.0) then
          do m = 1, nr3
            r_stack(l_temp+m-1) = onstepcrd(m)
          end do
        else
          do m = 1, nr3
            r_stack(l_temp+m-1) = x(m)
          end do
        endif
        call iwrap2(n_iwrap_mask_atoms, iwrap_mask_atoms, r_stack(l_temp), &
                    box_center)
        call corpac(r_stack(l_temp), 1, nrx, MDCRD_UNIT, loutfm)
        call corpac(box, 1, 3, MDCRD_UNIT, loutfm)
        call free_stack(l_temp, routine)
      end if
      ! End branch for wrapping molecules (iwrap)

      ! If using variable QM solvent, try to write a new pdb file
      ! with the QM coordinates for this step. This is done here
      ! to keep the PDB file in sync with the mdcrd file, which
      ! makes it easier to check later.
      if (qmmm_nml%vsolv > 0 .and. qmmm_nml%verbosity == 0) &
        call qm_print_coords(nstep,.false.)
    end if
    ! End contingency for trajectory writing on the master process (itdump) }}}

    !    Velocity archive: {{{
    if (ivdump) then

      ! NOTE - This assumes that if numextra > 0, then velocities are
      !        found in the array v...
      if (numextra > 0) call zero_extra_pnts_vec(v, ix)
#ifdef MPI
      ! Write out current replica#, exchange#, step#, and mytargettemp
      ! If mdloop==0 this is a normal md run (since REMD never calls corpac
      ! when mdloop==0) and we don't want the REMD header.
      if (mdloop > 0 .and. loutfm) then
        if (trxsgld) then
          write (MDVEL_UNIT,'(a,4(1x,i8))') "RXSGLD ", repnum, mdloop, &
                total_nstep, stagid
        else
          write (MDVEL_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, mdloop, &
                total_nstep, my_remd_data%mytargettemp
        end if
      end if
#endif /* MPI */

      if (ipimd == NMPIMD .or. ipimd == CMD) then
        call corpac(cartvel, 1, nrx, MDVEL_UNIT, loutfm)
      else if (isynctraj.ne.0) then
        call corpac(onstepvel, 1, nrx, MDVEL_UNIT, loutfm)
      else
        call corpac(v, 1, nrx, MDVEL_UNIT, loutfm)
      endif
    end if
    ! }}}

!------------------------------------------------------------------------------
    ! Force archive lam81
    if (ifdump .and. (abfqmmm_param%abfqmmm == 1)) &
      call corpac(f_or,1,nrx,MDFRC_UNIT,loutfm)

    ! Energy archive: (total_nstep set in Step 5.)
    if ( ntwe > 0 ) then
      if (mod(total_nstep,ntwe) == 0 .and. onstep) &
        call mdeng(15,nstep,t,ener,onefac,ntp,csurften)
    end if

    if (ioutfm > 0) then
      if (itdump) call end_binary_frame(MDCRD_UNIT)
      if (ivdump .and. ntwv > 0) call end_binary_frame(MDVEL_UNIT)
      if (ifdump .and. ntwf > 0) call end_binary_frame(MDFRC_UNIT)
    end if

#ifdef MPI
    if (ievb .ne. 0) call out_evb(nstep)
#endif /* MPI */

!------------------------------------------------------------------------------
    !    General printed output:  {{{
    if (lout) then
      if (facc .ne. 'A') rewind(7)

      ! Conserved quantity for Nose'-Hoover based thermostats.
      if (ipimd  == 0 .and. ntt > 4 .and. ntt <= 8) then
        Econserved = ener%kin%tot + ener%pot%tot + E_nhc
        if (ntp > 0) Econserved = Econserved + pres0 / pconv * volume
#ifdef MPI
        if (worldrank.eq.0) then
#endif /* MPI */
        write(file_nhc,'(I10,F14.4)') nstep, Econserved
#ifdef MPI
        end if
#endif /* MPI */
      endif
#ifdef LES
      if (ipimd > 0 .and. ntt == 4) then
        Econserved = ener%kin%tot + ener%pot%tot + E_nhc
        Econserved = Econserved   + Epot_spring
        if (ntp > 0) Econserved = Econserved + pres0 / pconv * volume
        write(file_nhc, '(I10,F14.4)') nstep, Econserved
      endif
      if (ipimd == CMD) then
        ener%kin%tot  = eke_cmd*0.5d0
        ener%kin%solv = ener%kin%tot
        ener%tot   = ener%kin%tot + ener%pot%tot
      end if
#else
      if (ipimd > 0) then
        ener%tot = 0.d0
        ener%kin%tot = 0.d0

        ! Conserved quantity for Nose'-Hoover thermostat.
        if (ntt == 4) then
          Econserved = totener%kin%tot + totener%pot%tot + E_nhc
          Econserved = Econserved + Epot_spring
          if (ntp > 0) Econserved = Econserved + volume*(pres0/pconv)
#  ifdef MPI
          if (worldrank == 0) then
#  endif /* MPI */
          write(file_nhc,'(I10,F14.4)') nstep, Econserved
#  ifdef MPI
          end if
#  endif /* MPI */
        endif
#  ifdef MPI
        if (worldrank == 0) then
#  endif /* MPI */
        call pimd_report(nstep, t, pimd_unit, totener, onefac)
#  ifdef MPI
        end if
#  endif /* MPI */
      end if
#endif /* LES */
      call prntmd(total_nstep, t, ener, onefac, 7, .false.)
      if (ischeme > 0 .and. rem == 0 .and. master) then
        write(fh_ek, '(I10,2(1x,F14.4))') nstep, ekph, ekph*onefac(1)
      endif

#ifdef MPI
      ! AWG FIXME - this should be in a subroutine
      ! Print corrected energy for adaptive qm/mm runs.
      ! Note: nstep has already been increased here
      !       (it was not increased when adaptive_qmmm() was called above)
      if (qmmm_nml%vsolv > 1) then
        if (masterrank == 0) then
          if (aqmmm_flag > 0 .and. nstep > aqmmm_flag) then
            etotcorr = corrected_energy + kinetic_E_save(aqmmm_flag)
            nstepadc = nstep - aqmmm_flag + 1
            tadc = t - dt * (dble(aqmmm_flag - 1))
            write(6, '(a)') ' Adaptive QM/MM energies:'
            write(6, '(x,a,i5,x,a,f11.4,x,2(a,f15.4,x))') 'adQMMM STEP=', &
                  nstepadc, 'TIME(PS)=', tadc, 'ETC=', etotcorr, &
                  'EPC=', corrected_energy

            ! print total energy for adaptive qm/mm into a separate file
            ! when qmmm_vsolv%verbosity > 0
            ! set reference energy to zero only for energy dumping purposes
            if (flag_first_energy) then
              flag_first_energy = .false.
              adqmmm_first_energy = etotcorr
              etotcorr = 0.0d0
            else
              etotcorr = etotcorr - adqmmm_first_energy
            end if
            if (qmmm_vsolv%verbosity > 0) then
              open(80,file='adqmmm_tot_energy.dat',position='append')
              write(80,'(i9,5x,f11.4,5x,f15.4)') nstepadc, tadc, etotcorr
              close(80)
            end if
          end if
        end if
      end if

      ! ANI energy
      if (ani_nml%ifmlp) then
        if (masterrank == 0) then
          write(6, '(a)') ' ANI Properties:'
          write(6, '(2(f16.6))') ani_struct%ani_energy
        end if
      end if
#endif

#ifdef MPI /* SOFT CORE */
      if (ifsc .ne. 0) call sc_print_energies(6, sc_ener)
      if (ifsc .ne. 0) call sc_print_energies(7, sc_ener)
#endif
      if (ifcr > 0 .and. crprintcharges > 0) &
        call cr_print_charge(xx(l15), total_nstep)

!------------------------------------------------------------------------------
      ! Output for Centroid MD
#ifdef LES
      if (ipimd == CMD) then
        ncmd = 0
        do iatom = 1, natom
          if (cnum(iatom) == 0 .or. cnum(iatom) == 1) then
            xcmd(ncmd+1) = x(3*iatom-2)
            xcmd(ncmd+2) = x(3*iatom-1)
            xcmd(ncmd+3) = x(3*iatom)
            vcmd(ncmd+1) = v(3*iatom-2)
            vcmd(ncmd+2) = v(3*iatom-1)
            vcmd(ncmd+3) = v(3*iatom)
            ncmd = ncmd+3
          endif
        enddo
        write(file_pos_cmd, '(10f8.3)') xcmd(1:ncmd)
        write(file_vel_cmd, '(10f8.3)') vcmd(1:ncmd)
        write(file_pos_cmd, '(10f8.3)') box(1:3)
        eke_cmd = eke_cmd * 0.5d0
        etot_cmd = eke_cmd + ener%pot%tot
        if (eq_cmd) then
          temp_cmd = eke_cmd / (boltz2 * dble(3*natomCL))
        else
          temp_cmd = eke_cmd / (boltz2 * dble(3*(natomCL-1)))
        endif
      endif
#else
      if (ipimd == CMD .and. mybeadid == 1) then
        write(file_pos_cmd, '(10f8.3)') x(1:3*natom)
        write(file_vel_cmd, '(10f8.3)') v(1:3*natom)
        write(file_pos_cmd, '(10f8.3)') box(1:3)
        eke_cmd = eke_cmd * 0.5d0
        etot_cmd = eke_cmd + totener%pot%tot
        if (eq_cmd) then
          temp_cmd = eke_cmd / (boltz2*dble(3*natom))
        else
          temp_cmd = eke_cmd / (boltz2*dble(3*(natom-1)))
        endif
      end if
#endif /* LES */

!------------------------------------------------------------------------------
      ! Print QMMM Muliken Charges if needed
      if (qmmm_nml%ifqnt) then
        if (qmmm_nml%printcharges .and. qmmm_mpi%commqmmm_master) then
          call qm2_print_charges(nstep, qmmm_nml%dftb_chg, &
                                 qmmm_struct%nquant_nlink, &
                                 qm2_struct%scf_mchg, &
                                 qmmm_struct%iqm_atomic_numbers)
        end if
      end if
      if (qmmm_nml%printdipole .ne. 0) &
        call qmmm_dipole(x, xx(Lmass), ix(i02), ih(m02), nres)

!------------------------------------------------------------------------------
      ! Begin dipole printing.  Also output dipole information if
      ! the dipoles namelist has been specified and corresponding
      ! groups defined.  Check input unit 5 for namelist &dipoles.
      ! We expect to find &dipoles followed by a group specification
      ! of the dipoles to print.
      call nmlsrc('dipoles', 5, prndipfind)
      if (prndipfind .ne. 0) then
        ! We calculate the dipoles
        write(6,*) '------------------------------- DIPOLE INFO -&
                    &---------------------------------'
        write(6,9018) nstep, t
        9018 format(/1x, 'NSTEP =', i7, 1x, 'TIME(PS) =', f10.3)

        ! Get the groups for the dipoles - Ideally we only really want
        ! to call this the once but for the time being I will call it
        ! every time
        read (5, '(a)') prndiptest
        call rgroup(natom, natc, nres, prndipngrp, ix(i02), ih(m02), &
                    ih(m04), ih(m06), ih(m08), ix(icnstrgp), jgroup, indx, &
                    irespw, npdec, xx(l60), 0, 0, 0, idecomp, 5, .false.)

        ! Need to rewind input file after rgroup
        ! so it is available when we next loop through
        rewind(5)
        if (prndipngrp > 0) then
          ! prndipngrp - holds number of groups specified + 1
          ! ix(icnstrgp) - holds map of group membership for each atom
          ! x(lcrd) - X,Y,Z coords of atoms - (3,*)
          ! x(l15) - Partial Charges
          ! x(linddip) - induced dipoles X,Y,Z for each atom (3,*)
          ! x(Lmass) - Mass of each atom
          call printdip(prndipngrp, ix(icnstrgp), xx(lcrd), &
                        xx(l15), xx(linddip), xx(Lmass), natom)
        end if
        write(6,*) '----------------------------- END DIPOLE INFO -&
                    &-------------------------------'
      end if
      ! End of dipole printing

      if (nmropt > 0) call nmrptx(6)
      if (infe > 0) call nfe_prt(6)
      if (itgtmd == 2) then
        emtmd = 0.0d0
        call mtmdcall(emtmd, xx(lmtmd01), ix(imtmd02), x, f, ih(m04), &
                      ih(m02), ix(i02), ih(m06), xx(lmass), natom, &
                      nres, 'PRNT')
      end if
      call amflsh(7)
    end if
    ! end of giant "if (lout)" contingency related to data output }}}

!------------------------------------------------------------------------------
    !    Output running averages: {{{
    ! total_nstep = Total nstep REMD/MD, set in step 5
    if (ntave > 0) then
      if (mod(total_nstep,ntave) == 0 .and. onstep) then
        write(6, 542)
#ifdef RISMSANDER
        if (rismprm%rism == 1) then
          tspan = ntave / mylcm(nrespa, rismprm%rismnrespa)
        else
          tspan = ntave / nrespa
        end if
#else
        tspan = ntave / nrespa
#endif /* RISMSANDER */

        ! Update all elements of these sequence types
        enert_tmp  = enert - enert_old
        enert2_tmp = enert2 - enert2_old
        enert_old  = enert
        enert2_old = enert2
        enert_tmp  = enert_tmp/tspan
        enert2_tmp = enert2_tmp/tspan - enert_tmp*enert_tmp
        call zero_neg_values_state(enert2_tmp)
        enert2_tmp = sqrt(enert2_tmp)
#ifdef MPI
        if (ievb /= 0) then
          evb_nrg_tmp (:) = evb_nrg_ave(:) - evb_nrg_old (:)
          evb_nrg_tmp2(:) = evb_nrg_rms(:) - evb_nrg_old2(:)
          evb_nrg_old (:) = evb_nrg_ave(:)
          evb_nrg_old2(:) = evb_nrg_rms(:)
          evb_nrg_tmp (:) = evb_nrg_tmp (:) / tspan
          evb_nrg_tmp2(:) = evb_nrg_tmp2(:) / tspan - evb_nrg_tmp(:)**2
          evb_nrg_tmp2(:) = max(evb_nrg_tmp2(:), 0.0d0)
          evb_nrg_tmp2(:) = sqrt(evb_nrg_tmp2(:))
        endif
        if (ifsc .ne. 0) then
          do m = 1,ti_ene_cnt
            sc_ener_tmp(m) = sc_ener_ave(m) - sc_ener_old(m)
            sc_ener_tmp2(m) = sc_ener_rms(m) - sc_ener_old2(m)
            sc_ener_old(m) = sc_ener_ave(m)
            sc_ener_old2(m) = sc_ener_rms(m)
            sc_ener_tmp(m) = sc_ener_tmp(m) / tspan
            sc_ener_tmp2(m) = sc_ener_tmp2(m)/tspan - sc_ener_tmp(m)**2
            if (sc_ener_tmp2(m) < 0.0d0) then
              sc_ener_tmp2(m) = 0.0d0
            end if
            sc_ener_tmp2(m) = sqrt(sc_ener_tmp2(m))
          end do
        end if
        if (ievb .ne. 0) evb_frc%evb_ave = .true.
#endif
#ifdef RISMSANDER
        if (rismprm%rism == 1) then
          write(6, 540) ntave / mylcm(nrespa, rismprm%rismnrespa)
        else
          write(6, 540) ntave/nrespa
        end if
#else
        write(6, 540) ntave/nrespa
#endif /* RISMSANDER */
        call prntmd(total_nstep, t, enert_tmp, onefac, 0, .false.)
#ifdef MPI
        if (ifsc .ne. 0) call sc_print_energies(6, sc_ener_tmp)
        if (ievb .ne. 0) evb_frc%evb_rms = .true.
#endif /* MPI */
        write(6, 550)
        call prntmd(total_nstep, t, enert2_tmp, onefac, 0, .true.)
#ifdef MPI /* SOFT CORE */
        if (ifsc .ne. 0) call sc_print_energies(6, sc_ener_tmp2)
#endif /* MPI */
        if (icfe > 0) then
#ifdef RISMSANDER
          if (rismprm%rism == 1) then
            write (6, 541) ntave / mylcm(nrespa, rismprm%rismnrespa)
          else
            write (6, 541) ntave/nrespa
          end if
#else
          write(6,541) ntave/nrespa
#endif /* RISMSANDER */
          edvdl_r = edvdl_r/tspan
          edvdl_r%pot%dvdl = enert_tmp%pot%dvdl  ! fix for DV/DL output
          edvdl_r%virvsene = 0.d0 ! virvsene should not but included here
          call prntmd(total_nstep, t, edvdl_r, onefac, 0, .false.)
          edvdl_r = null_state_rec
        end if
        write(6,542)
      end if
    end if
    ! End contingency to output running averages (ntave)      }}}
  end if
  ! End output work on the master process

  ! Soft core output:  {{{
#ifdef MPI /* SOFT CORE */
  if (ntave > 0 .and. icfe > 0 .and. dynlmb > 0) then
    if (mod(nstep,ntave) == 0 .and. onstep) then

!------------------------------------------------------------------------------
      ! For runs with dynamically changing lambda, raise lambda here
      ! and flush all buffers for the next averages
      clambda = clambda + dynlmb
      call sc_change_clambda(clambda)
      if (master) then
        sc_ener(1:ti_ene_cnt) = 0.0d0
        sc_ener_ave(1:ti_ene_cnt) = 0.0d0
        sc_ener_rms(1:ti_ene_cnt) = 0.0d0
        sc_ener_old(1:ti_ene_cnt) = 0.0d0
        sc_ener_old2(1:ti_ene_cnt) = 0.0d0
        enert = null_state_rec
        enert2 = null_state_rec
        enert_old = null_state_rec
        enert2_old = null_state_rec
        write (6, *)
        write (6, '(a,f12.4,a,f12.4)') &
              'Dynamically changing lambda: Increased clambda by ', &
              dynlmb, ' to ', clambda
        write (6,*)
      end if
    end if
  end if
#endif
  ! }}}

  ! }}}
!------------------------------------------------------------------------------
  ! Major cycle back to new step unless we have reached our limit:
  call trace_integer( 'end of step', nstep )
  call trace_output_mpi_tally( )
  call timer_stop(TIME_VERLET)
#if !defined(DISABLE_NFE) && defined(NFE_ENABLE_BBMD)
  if (infe == 1) then
    call xdist(x, xx(lfrctmp), natom)
    call xdist(v, xx(lfrctmp), natom)
    call nfe_on_mdstep(ener%pot%tot, x, v, ekmh)
    vold(istart3:iend3) = v(istart3:iend3)
    do im = 1, iscale
       vold(nr3+im) = v(nr3+im)
    end do
  end if
#endif /* DISABLE_NFE is NOT defined, but NFE_ENABLE_BBMD is */

  !  RISM debug: {{{
#if defined(RISMSANDER) && defined(RISM_DEBUG)
  if (rismprm%rism == 1) then
    angvel=0
    do m = 1, natom
      r = x((m-1)*3+1:(m-1)*3+3)-cm
      call cross(r,v((m-1)*3+1:(m-1)*3+3),rxv)
      angvel = angvel + rxv/sum(r**2)
    end do
    moi = 0
    erot = 0
    do m = 1, natom
      r = x((m-1)*3+1:(m-1)*3+3)-cm
      call cross(r, v((m-1)*3+1:(m-1)*3+3), rxv)
      proj = angvel * (sum(r*angvel) / sum(angvel**2))
      moi=moi+amass(m)*sum((r-proj)**2)
      erot = erot + .5*amass(m)*sum((r-proj)**2)*sum((rxv/sum(r**2))**2)
    end do
  end if
#endif /*RISMSANDER && RISM_DEBUG*/
  ! }}}

  if (abfqmmm_param%abfqmmm == 1) then
#ifdef MPI
    call xdist(v, xx(lfrctmp), natom)
#endif
    abfqmmm_param%v(1:nr3+iscale) = v(1:nr3+iscale)
    deallocate(f_or, stat=ierr)
    return
  end if

  if (plumed .ne. 0 .and. plumed_stopflag .ne. 0) goto 480

  ! This is where we actually cycle back to the next MD step
  if (nstep < nstlim) goto 260

  480 continue

#ifdef MPI
!------------------------------------------------------------------------------
  ! Replica Exchange MD post-dynamics work {{{
  if (next_rem_method == 1) then
    remd_ekmh = ekmh

    ! Hybrid REMD
    if (numwatkeep >= 0) then
      ! This is a hybrid REMD run. Get energy of
      ! stripped system for next exchange.
      call hybrid_remd_ene(xx, ix, ih, ipairs, qsetup, numwatkeep, hybridgb, &
                           igb, ntr, nspm, t, temp0, ntb, cut, ener, &
                           do_list_update, nstep, onefac)
    else

      ! The positions are currently one step ahead of the energy ener%pot%tot,
      ! since force was called prior to the position propagation. Thus, call
      ! force one more time to update ener%pot%tot to reflect the current
      ! coordinates.
      call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
                 xx(l98), xx(l99), qsetup, do_list_update,nstep)
    end if
    ! End branch for additional energy computations in Hybrid REMD

    ! Set myeptot, mytemp, and mytargettemp
    my_remd_data%mytemp = ener%kin%tot * onefac(1)
    my_remd_data%myeptot = ener%pot%tot
    my_pressure = pres0
    my_volume = ener%volume

    my_remd_data%mytargettemp = temp0
#  ifdef VERBOSE_REMD
    if (master) write(6, '(a,f15.4,2(a,f6.2))') "| REMD: myEptot= ", &
           my_remd_data%myeptot, " myTargetTemp= ", &
           my_remd_data%mytargettemp, " mytemp= ", my_remd_data%mytemp
#  endif /* VERBOSE_REMD */
#  ifdef LES
   else if(next_rem_method == 2 ) then
    my_remd_data%mytemp       = ener%kin%solv * onefac(3)
    my_remd_data%myeptot      = ener%eptot
    my_remd_data%mytargettemp = temp0les
    my_pressure = pres0
    my_volume = ener%volume
#  endif /* LES */
  else if (next_rem_method == 3) then
    remd_ekmh = ekmh
    if (mdloop > 0) then
      my_remd_data%mytemp = ener%kin%tot * onefac(1)
    end if
    my_remd_data%mytargettemp = temp0

    ! This call to force will bring all energies up-to-date
    call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
               xx(l98), xx(l99), qsetup, do_list_update,nstep)
    my_remd_data%myeptot = ener%pot%tot
    my_pressure = pres0
    my_volume = ener%volume
    ! Call nmrdcp to decrement the NMR counter, since this should not count as
    ! a real step (JMS 2/12). This is OK, since the counter got incremented at
    ! the very end of nmrcal, so we haven't already printed an unwanted value.
    if (nmropt /= 0) then
      call nmrdcp
    end if

    ! Call xdist such that master has all the velocities
    call xdist(v, xx(lfrctmp), natom)
  else if (next_rem_method == 4 .or. next_rem_method == 5) then
    remd_ekmh = ekmh
  endif
#endif /* MPI */
  ! End Replica Exchange MD post-dynamics work }}}

  ! Stochastic Isokinetic Nose-Hoover RESPA (SINR) integrator clean-up
  if (ntt==10) call sinr_cleanup(sinrdata)

  ! Print averages {{{
#ifdef MPI
  ! Thermodynamic Integration decomposition
  if (icfe .ne. 0 .and. idecomp .ne. 0) then
    if (idecomp == 1 .or. idecomp == 2) call collect_dec(nrs)
  end if

  ! Turn off avg. for REMD. and explicit solvent CpHMD, since it's not
  ! accumulated correctly in that case for each compiler
  if (master .and. rem == 0) then
#else
  if (master) then
#endif /*MPI*/
    tspan = nvalid
    if (nvalid > 0) then

      ! Update all elements of these sequence types
      enert  = enert/tspan
      enert2 = enert2/tspan - enert*enert
      call zero_neg_values_state(enert2)
      enert2 = sqrt(enert2)
      edvdl = edvdl/tspan

      ! For PIMD/NMPIMD/CMD/RPMD averages
      if (ipimd>0) then
        totenert = totenert/tspan
        totenert2 = totenert2/tspan - (totenert*totenert)
        call zero_neg_values_state(totenert2)
        totenert2 =  sqrt(totenert2)
      end if

      if (ischeme > 0 .and. rem == 0 .and. master) then
        ekhf = ekhf / tspan
        ekhf2 = ekhf2/tspan - ekhf*ekhf
        ekhf2 = sqrt(ekhf2)
        write(fh_ek, '(a,i12,a)') ' # average of EKIN and TEMP over ', nvalid, &
                                  '  steps'
        write(fh_ek, '(2(1x,F14.4))') ekhf, ekhf*onefac(1)
        write(fh_ek, '(a,i12,a)') ' # rms of EKIN and TEMP over ', nvalid, &
                                  '  steps'
        write(fh_ek, '(2(1x,F14.4))') ekhf2, ekhf2*onefac(1)
      endif
#ifdef MPI
      if (ievb .ne. 0) then
        evb_nrg_ave(:) = evb_nrg_ave(:) / tspan
        evb_nrg_rms(:) = evb_nrg_rms(:) / tspan - evb_nrg_ave(:)**2
        evb_nrg_rms(:) = max( evb_nrg_rms(:), 0.0d0 )
        evb_nrg_rms(:) = sqrt( evb_nrg_rms(:) )
      endif
      if (ifsc .ne. 0) then
        do m = 1, ti_ene_cnt
          sc_ener_ave(m) = sc_ener_ave(m)/tspan
          sc_ener_rms(m) = sc_ener_rms(m)/tspan - sc_ener_ave(m)**2
          if (sc_ener_rms(m) < 0.0d0) then
            sc_ener_rms(m) = 0.0d0
          end if
          sc_ener_rms(m) = sqrt(sc_ener_rms(m))
        end do
      end if
      if (ievb .ne. 0) evb_frc%evb_ave = .true.
#endif
      write(6, 540) nvalid
      call prntmd(total_nstep, t, enert, onefac, 0, .false.)
#ifdef MPI /* SOFT CORE */
      if (ifsc .ne. 0) call sc_print_energies(6, sc_ener_ave)
      if (ievb .ne. 0) evb_frc%evb_rms = .true.
      if (ipimd > 0 .and. worldrank == 0) then
        write(pimd_unit, 540) nvalid
        call pimd_report(nstep, t, pimd_unit, totenert, onefac)
        write(pimd_unit, 550)
        call pimd_report(nstep, t, pimd_unit, totenert2, onefac)
      endif
#endif
      if (nmropt > 0) call nmrptx(6)
      write(6, 550)
      call prntmd(total_nstep, t, enert2, onefac, 0, .true.)
#ifdef MPI
      if (ifsc .ne. 0) call sc_print_energies(6, sc_ener_rms)
      if (ifsc .ne. 0) call sc_print_dvdl_values()
      if (icfe > 0) then
        write(6, 541) nvalid
        edvdl%pot%dvdl = enert%pot%dvdl  ! fix for DV/DL output
        edvdl%virvsene = 0.d0 ! virvsene should not but included here
        call prntmd(total_nstep, t, edvdl, onefac, 0, .false.)

        ! Thermodynamic Integration decomposition
        if (worldrank == 0 .and. idecomp .ne. 0) then
          call checkdec(idecomp)
          if (idecomp == 1 .or. idecomp == 2) then
            call printdec(ix)
          end if
        end if
      end if
#endif /* MPI */
      if (nmropt >= 1) then
        write(6, 500)
        if (iredir(7) .ne. 0) then
          call pcshift(-1, x, f)
        end if
        call ndvptx(x, f, ih(m04), ih(m02), ix(i02), nres, xx(l95), &
                    natom, xx(lwinv), xx(lnmr01), ix(inmr02), 6)
      end if

!------------------------------------------------------------------------------
      ! Print Born radii statistics
      if (rbornstat == 1 .and. (igb .ne. 0 .or. ipb .ne. 0)) then

        ! Born radii stats collected every nrespai step not nrespa step
        tspan = nvalidi
        write(6, 580) nstep
        write(6, 590)
        do m = 1, natom
          xx(l188-1+m) = xx(l188-1+m)/tspan
          xx(l189-1+m) = xx(l189-1+m)/tspan - xx(l188-1+m)*xx(l188-1+m)
          xx(l189-1+m) = sqrt(xx(l189-1+m))
          write(6, 600) m, xx(l186-1+m), xx(l187-1+m), xx(l188-1+m), &
                        xx(l189-1+m)
        end do
      end if

      enert%kin%tot = enert%kin%tot*onefac(1)
      enert2%kin%tot = enert2%kin%tot*onefac(1)
      enert%kin%solt = enert%kin%solt*onefac(2)
      enert2%kin%solt = enert2%kin%solt*onefac(2)
      enert%kin%solv = enert%kin%solv*onefac(3)
      enert2%kin%solv = enert2%kin%solv*onefac(3)
      temp = enert%kin%tot
    end if
    ! End contingency for nvalid > 0, signifying that
    ! all energies must be calculated
    if (ischeme > 0 .and. rem == 0 .and. master) close(fh_ek)
    if (ntp > 0 .and. barostat == 2) call mcbar_summary
  end if
  ! End of contingency for work on the master process;
  ! this started about 120 lines ago and the way it started
  ! depends on whether MPI is part of the compilation. }}}

#ifdef MPI
  if (ievb .ne. 0) then
    call evb_dealloc
#  if defined(LES)
    if (master) call evb_pimd_dealloc
#  endif /* LES */
  endif
#endif /* MPI */
  if (icfe .ne. 0) then
    deallocate( frcti, stat = ierr )
    REQUIRE( ierr == 0 )
  end if
  if (plumed .ne. 0) call plumed_f_gfinalize()

  500 format(/,' NMR restraints on final step:'/)
  540 format(/5x,' A V E R A G E S   O V E R ',i7,' S T E P S',/)
  541 format(/5x,' DV/DL, AVERAGES OVER ',i7,' STEPS',/)
  542 format('|',79('='))
  550 format(/5x,' R M S  F L U C T U A T I O N S',/)
  580 format('STATISTICS OF EFFECTIVE BORN RADII OVER ',i7,' STEPS')
  590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
  600 format(i4,2x,4f12.4)
  call trace_exit( 'runmd' )
  return
end subroutine runmd
