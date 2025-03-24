#include "../include/dprec.fh"
#include "../include/assert.fh"

module pol_gauss_mdin
  implicit none
  private
  integer,save :: ipgm=0,do_pol_gauss_nonbond=1
  _REAL_,save :: dipole_scf_tol = 0.01d0
  integer,save :: use_average = 0
  integer,save :: dipole_scf_init = 3
  integer,save :: dipole_scf_init_order = 3
  integer,save :: dipole_scf_init_step = 2
  integer,save :: scf_solv_opt = 3
  integer,save :: scf_sor_niter = 100
  _REAL_,save :: scf_sor_coefficient = 0.65d0
  integer,save :: scf_cg_niter = 50
  integer,save :: scf_local_niter = 3
  _REAL_,save :: scf_local_cut=4.0d0
  _REAL_,save :: ee_dsum_cut=9.0d0
  _REAL_,save :: ee_damped_cut=4.5d0
  integer, save :: pol_gauss_verbose = 0
  _REAL_,save :: erfc_tol = 1.0e-5
  _REAL_,save :: ee_gauss_cut = 3.12342d0 ! corresponding to max error of 10^-5 in erfc

  public POL_GAUSS_read_mdin,ipgm, &
         do_pol_gauss_nonbond,pol_gauss_verbose, &
         dipole_scf_tol,use_average,dipole_scf_init,dipole_scf_init_order,dipole_scf_init_step,&
         scf_solv_opt,scf_cg_niter,scf_local_niter, &
         scf_sor_coefficient,scf_sor_niter, &
         ee_dsum_cut,ee_damped_cut,scf_local_cut
  contains
!-------------------------------------------------------------------------------
subroutine POL_GAUSS_read_mdin(nf)
  use file_io_dat
  implicit none

  integer,intent(in) :: nf

  integer :: do_recip=1,do_direct=1,do_self=1,do_vdw=1,do_induced=1

  namelist/pol_gauss/do_pol_gauss_nonbond, &
                     do_recip,do_direct,do_self, &
                     do_vdw,do_induced,pol_gauss_verbose,use_average, &
                     dipole_scf_tol,dipole_scf_init,dipole_scf_init_order,dipole_scf_init_step, &
                     scf_solv_opt,scf_sor_coefficient,scf_sor_niter, &
                     scf_cg_niter,scf_local_niter,scf_local_cut, &
                     ee_dsum_cut,ee_damped_cut, &
                     erfc_tol,ee_gauss_cut

  read(nf,nml=pol_gauss)
  call pGM_NONBOND_set_user_bit(do_recip,do_direct,do_self,do_vdw,do_induced)
end subroutine POL_GAUSS_read_mdin
!-------------------------------------------------------------------------------
end module pol_gauss_mdin
!-------------------------------------------------------------------------------
