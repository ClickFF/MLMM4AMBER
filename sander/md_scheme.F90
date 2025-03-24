! include for _REAL_ definition, compatible with runmd.F90
#include "../include/dprec.fh"

! module for "middle" scheme for MD / PIMD
!
! References:
! 1. Zhang#, Z.; Liu#, X.; Chen, Z.; Zheng, H.; Yan, K.; Liu*, J.
!   "A unified thermostat scheme for efficient configurational sampling
!    for classical/quantum canonical ensembles via molecular dynamics",
!    The Journal of Chemical Physics 147, 034109 (2017) 
! 2. Zhang, Z.; Yan, K.; Liu, X.; Liu*, J.
!   "A leap-frog algorithm-based efficient unified thermostat scheme 
!    for molecular dynamics",
!    Chinese Science Bulletin 63(33), 3467-3483 (2018)
!
module md_scheme
  implicit none
  ! the scheme for integration algorithm
  integer :: ischeme     
  integer,parameter :: &
    ENUM_NO_SCHEME       = 0, & ! default AMBER
    ENUM_LFMIDDLE_SCHEME = 1    ! leap-frog middle scheme

  ! thermostat options
  integer :: ithermostat ! thermostat method
  integer,parameter :: &
    ENUM_NO_THERM   = 0, & ! no thermostatting
    ENUM_LGV_THERM  = 1, & ! Langevin thermostat
    ENUM_ADS_THERM  = 2    ! Andersen thermostat

  ! thermostat parameter, unit is 1/ps
  !   = gamma (friction coefficient) for Langevin dynamics
  !   = nu    (collision frequency)  for Andersen thermostat
  _REAL_ :: therm_par    

#include "../include/md.h"
#ifdef MPI
#  include "parallel.h"
   include"mpif.h"
#endif

contains
  !-- integration step for the thermostat
  ! @param (in AMBER internal unit)
  !   [inout] v, velocity array (for each DOF)
  !   [in] winv, inverse of mass array (for each atom)
  !   [in] dtx, time step, see runmd.F90
  !   [in] istart, start index of atom, see runmd.F90
  !   [in] iend, end index of atom, see runmd.F90
  ! @external
  !   amrand: random number in uniform dist., in ../lib/random.F90
  !   gauss:  random number in gaussian dist., in ../lib/random.F90
  subroutine thermostat_step(v, winv, dtx, istart, iend, iskip_start, iskip_end)
    !use constants, only: KB
    use random
#ifdef LES
    use les_data, only: cnum, temp0les
#endif    
    implicit none

    character(kind=1,len=*),parameter :: routine="thermostat_step"
    _REAL_ :: v(3*nrp), winv(nrp), dtx
    integer :: istart,iend, iskip_start,iskip_end

    integer :: i, j, k, idof
    ! rand          : random number
    ! dt            : time interval for this step, in ps
    ! lgv_c1, lgv_c2: parameters used in Langevin dynamics
    ! ads_prob      : collision probability for Andersen thermostat
    ! rtKT          : square root of kB*Temperature
    ! stdvel        : standard deviation of velocities within Maxwell distribution
    _REAL_ :: rand, dt, lgv_c1, lgv_c2, ads_prob, rtkT, stdvel

#ifdef LES
    _REAL_ :: vscalt
#endif    
#ifdef MPI
    integer :: ierr
#endif 

    double precision :: boltz

    boltz = 8.31441d0 / 4.184d3

    dt = dtx / 20.455d0
    select case (ithermostat)
    case (ENUM_LGV_THERM) ! Langevin thermostat
      ! therm_par in ps^-1, dt in ps
      lgv_c1 = exp(-therm_par*dt) 
      lgv_c2 = sqrt(1.0d0 - lgv_c1*lgv_c1)
      rtkT = sqrt(boltz*temp0)

      do i = 1, iskip_start
        call gauss(0.0d0, 1.0d0, rand)
      end do
      do i = istart, iend
        ! the unit of stdvel is Angstrom per 1/20.455 ps, same as vel
        ! 1 Angstrom per 1/20.455 ps = 1 sqrt(kcal/mol / amu)
        stdvel = rtkT * sqrt(winv(i))
#ifdef LES
        if (temp0les >= 0 .and. cnum(i) > 0) then
          stdvel = stdvel * sqrt(temp0les/temp0)
        end if
#endif
        do j = 1, 3
          call gauss(0.0d0, stdvel, rand) 
          idof = 3*(i-1) + j
          v(idof) = lgv_c1 * v(idof) + lgv_c2 * rand
        end do
      end do
      do i = 1, iskip_end
        call gauss(0.0d0, 1.0d0, rand)
      end do
    case (ENUM_ADS_THERM) ! Andersen thermostat
      ads_prob = 1.0d0 - exp(-therm_par*dt)
      rtkT = sqrt(boltz*temp0)
      call amrand(rand)
#ifdef MPI      
      ! one rand to rule them all
      call MPI_Bcast(rand, 1, MPI_DOUBLE_PRECISION, 0, commsander, ierr) 
#endif      
      if (rand < ads_prob) then
        do i = 1, iskip_start
          call gauss(0.0d0, 1.0d0, rand)
        end do

        do i = istart, iend
          stdvel = rtkT * sqrt(winv(i))
          do j = 1, 3
            idof = 3*(i-1) + j
            call gauss(0.0d0, stdvel, rand) 
            v(idof) = rand
          end do
        end do

        do i = 1, iskip_end
          call gauss(0.0d0, 1.0d0, rand)
        end do
#ifdef LES
        if (temp0les >= 0.0 .and. temp0 .ne. temp0les) then
          vscalt = sqrt(temp0les/temp0)
          do i = istart, iend
            if (cnum(i) > 0) then
              do j = 1,3
                idof = 3*(i-1)+j
                v(idof) = v(idof) * vscalt
              end do
            end if
          end do
        end if 
#endif      
      end if
    case default
      write(*,*) 'Error in '//routine//': unknown ithermostat'
      stop
    end select
  end subroutine

end module
