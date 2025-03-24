#include "../include/dprec.fh"
#include "../include/assert.fh"

!-------------------------------------------------------------------------------
module pol_gauss_interface
  implicit none
  private

  public POL_GAUSS_readparm,POL_GAUSS_deallocate,pGM_NonBond_eval

contains
!-------------------------------------------------------------------------------
subroutine POL_GAUSS_readparm(nf,numatoms)
   integer,intent(in) :: nf,numatoms

   call POL_GAUSS_check_parm_legal(nf)
   call pGM_NONBOND_readparm(nf,numatoms)
end subroutine POL_GAUSS_readparm
!-------------------------------------------------------------------------------
subroutine pGM_NONBOND_readparm(nf,natom)
  use pol_gauss_multipoles, only : pGM_MPOLE_readparm
  use pol_gauss_induced, only : pGM_INDUCED_readparm
  use pol_gauss_recip, only : pGM_RECIP_allocate

  integer,intent(in) :: nf,natom

  integer :: mpole_valid,vdw_valid,polar_valid,adjust_valid

  mpole_valid = pGM_MPOLE_readparm(nf,natom)
  polar_valid = pGM_INDUCED_readparm(nf,natom)
  ! alocate to initialize recip
  call pGM_RECIP_allocate(natom)
end subroutine pGM_NONBOND_readparm
!-------------------------------------------------------------------------------
subroutine POL_GAUSS_check_parm_legal(nf)
  use pol_gauss_mdin, only : ipgm
  integer,intent(in) :: nf

  integer :: iok,ionerr,indicator
  character(len=80) :: fmt,fmtin,dtype

  fmtin = '(I5)'
  dtype = 'POL_GAUSS_FORCEFIELD'
  ionerr = 1 ! not fatal if missing
  call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
  if ( iok == 0 )then !this data type found in prmtop
    read(nf,fmt)indicator
    if ( ipgm /= 1 )then
       write(6,*)'POL_GAUSS style prmtop, but pgm NOT set to one'
       call mexit(6,1)
    endif
    return ! pGM prmtop, pGM set to 1
  else
    if ( ipgm == 1 )then
       write(6,*)'NOT an POL_GAUSS style prmtop, but pgm IS set to one'
       call mexit(6,1)
    endif
    return ! NON pGM prmtop, pGM NOT set to 1
  endif
end subroutine POL_GAUSS_check_parm_legal
!-------------------------------------------------------------------------------
subroutine POL_GAUSS_deallocate()

  call pGM_NONBOND_deallocate()
end subroutine POL_GAUSS_deallocate
!-------------------------------------------------------------------------------
subroutine pGM_NONBOND_deallocate()
  use pol_gauss_multipoles,only : pGM_MPOLE_deallocate

  call pGM_MPOLE_deallocate()
end subroutine pGM_NONBOND_deallocate
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_eval(numatoms,crd,frc,sander_vir,x,nummols,molsiz,amass,ipairs, &
                       ntypes,iac,ico,cn1,cn2,& ! R. Luo: vdw
                       evdw,eelt,epolar,evdw_14,eelt_14, &
                       diprms,dipiter)
  use pol_gauss_mdin, only : do_pol_gauss_nonbond,pol_gauss_verbose
  use pol_gauss_multipoles,only : pGM_MPOLE_local_to_global, &
                               global_multipole, coulomb_const_kcal_per_mole
  use pol_gauss_induced,only : pGM_INDUCED_eval
  use nblist, only: recip,adjust_imagcrds, map_coords

  integer,intent(in) :: numatoms
  _REAL_ :: crd(3,*), amass(*)
  _REAL_,intent(inout) :: frc(3,*),sander_vir(4)
  _REAL_,intent(in) :: x(*)
  integer nummols,molsiz(*)
  integer, intent(in) :: ipairs(*)
  integer,intent(in) :: ntypes,iac(*),ico(*) ! R. Luo: vdw
  _REAL_,intent(in) :: cn1(*),cn2(*) ! R. Luo: vdw
  _REAL_,intent(out) :: evdw,eelt,epolar,evdw_14,eelt_14,dipiter,diprms

# include "box.h"

  _REAL_ :: pgm_virial(3,3), vdw_virial(3,3)

  evdw = 0.d0
  eelt = 0.d0
  epolar = 0.d0
  evdw_14 = 0.d0
  eelt_14 = 0.d0

  pgm_virial = 0.0
  vdw_virial = 0.0

  if ( do_pol_gauss_nonbond == 1 )then

    call zero_array(global_multipole,10*numatoms)

    ! update the imaged crds
    call map_coords(crd,numatoms,recip)
    call adjust_imagcrds(crd,numatoms)

    ! R. Luo: molecular moments can be read in the global frame for development
    ! the conversion to local is done here after crd is obtained
    !call pGM_MPOLE_global_to_local(crd)
    call pGM_MPOLE_local_to_global(crd)
    call pGM_INDUCED_eval(numatoms,crd,x,ipairs,diprms,dipiter)
    call pGM_NonBond_ene_frc(numatoms,crd,x,nummols,molsiz,amass,ipairs,ntypes,iac,ico,cn1,cn2,& ! R. Luo: vdw
                     eelt,epolar,evdw,evdw_14,frc,pgm_virial,vdw_virial)

    if ( ntb > 0 ) then
      if ( pol_gauss_verbose == 2 )then
        write(6,'(a,3(1x,g16.8))') &
                  ' pGM nonbond vir = ',pgm_virial(1,1),pgm_virial(1,2),pgm_virial(1,3)
        write(6,'(a,3(1x,g16.8))') &
                  ' pGM nonbond vir = ',pgm_virial(2,1),pgm_virial(2,2),pgm_virial(2,3)
        write(6,'(a,3(1x,g16.8))') &
                  ' pGM nonbond vir = ',pgm_virial(3,1),pgm_virial(3,2),pgm_virial(3,3)
      endif
    endif

  endif

  if ( ntb > 0 ) then
    ! the factor 0.5 is due to the fact that sander uses molecular kinetic energy as the
    ! kinetic part of virial, which is actually half of the real virial
    pgm_virial = 0.5 * coulomb_const_kcal_per_mole * pgm_virial
    vdw_virial = 0.5 * vdw_virial

    sander_vir(1) = sander_vir(1) + pgm_virial(1,1) + vdw_virial(1,1)
    sander_vir(2) = sander_vir(2) + pgm_virial(2,2) + vdw_virial(2,2)
    sander_vir(3) = sander_vir(3) + pgm_virial(3,3) + vdw_virial(3,3)
    sander_vir(4) = sander_vir(4) + pgm_virial(1,1) + pgm_virial(2,2) + pgm_virial(3,3) +&
                                    vdw_virial(1,1) + vdw_virial(2,2) + vdw_virial(3,3)
  end if

end subroutine pGM_NonBond_eval
!-------------------------------------------------------------------------------
end module pol_gauss_interface
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_perm_fields(numatoms,is_polarizable,crd,x,ipairs, &
                                   dip_field,cart_dipole_field)
  use pol_gauss_direct, only : pGM_DIRECT_perm_field
  use pol_gauss_recip, only : pGM_RECIP_perm_field
  use pol_gauss_self, only : pGM_SELF_perm_field
  integer,intent(in) :: numatoms
  logical,intent(in) :: is_polarizable(*)
  _REAL_,intent(in) :: crd(3,*)
  _REAL_,intent(in) :: x(*)
  integer, intent(in) :: ipairs(*)
  _REAL_,intent(out) :: dip_field(3,*),cart_dipole_field(3,*)

# include "ew_erfc_spline.h"
#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
   logical master
#endif

  call zero_array(dip_field,3*numatoms)
  call zero_array(cart_dipole_field,3*numatoms)

#ifdef MPI
  call mpi_comm_size(recip_comm,numtasks,ierr)
  call mpi_comm_rank(recip_comm,mytaskid,ierr)
  master = mytaskid.eq.0
#endif

  ! All computed fields are actually grad phi
  if (eedmeth.ne.4) call pGM_RECIP_perm_field(numatoms,crd,cart_dipole_field,x)

#ifdef MPI
  call mpi_comm_rank(commsander,mytaskid,ierr)
  call mpi_comm_size(commsander,numtasks,ierr)
  master = mytaskid.eq.0
#endif

  call pGM_DIRECT_perm_field(ipairs,x,cart_dipole_field)
  if (eedmeth.ne.4) call pGM_SELF_perm_field(numatoms,dip_field)

  ! R. Luo: I think this can be removed by using dip_field or
  ! card_dipole_field only
  dip_field(1:3,1:numatoms) = dip_field(1:3,1:numatoms) + cart_dipole_field(1:3,1:numatoms)

end subroutine pGM_NonBond_perm_fields
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_dip_dip_fields_short(numatoms,x,ind_dip,dip_field)
  use pol_gauss_direct, only : pGM_DIRECT_dip_dip_field_short

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: x(*)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(out) :: dip_field(3,*)

#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
  logical master
#endif

  call zero_array(dip_field,3*numatoms)

#ifdef MPI
  call mpi_comm_size(recip_comm,numtasks,ierr)
  call mpi_comm_rank(recip_comm,mytaskid,ierr)
  master = mytaskid.eq.0
#endif

#ifdef MPI
  call mpi_comm_rank(commsander,mytaskid,ierr)
  call mpi_comm_size(commsander,numtasks,ierr)
  master = mytaskid.eq.0
#endif
 
  ! All computed fields are actually grad phi
  call pGM_DIRECT_dip_dip_field_short(ind_dip,dip_field)

end subroutine pGM_NonBond_dip_dip_fields_short
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_dip_dip_fields(numatoms,x, &
                       ind_dip,dip_field)
  use pol_gauss_recip, only : pGM_RECIP_dipole_field
  use pol_gauss_direct, only : pGM_DIRECT_dip_dip_field
  use pol_gauss_self, only : pGM_SELF_dipole_field

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: x(*)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(out) :: dip_field(3,*)

# include "ew_erfc_spline.h"
#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
  logical master
#endif

  call zero_array(dip_field,3*numatoms)

#ifdef MPI
  call mpi_comm_size(recip_comm,numtasks,ierr)
  call mpi_comm_rank(recip_comm,mytaskid,ierr)
  master = mytaskid.eq.0
#endif

  ! All computed fields are actually grad phi
  if (eedmeth.ne.4) call pGM_RECIP_dipole_field(numatoms,x,ind_dip,dip_field)

#ifdef MPI
  call mpi_comm_rank(commsander,mytaskid,ierr)
  call mpi_comm_size(commsander,numtasks,ierr)
  master = mytaskid.eq.0
#endif
 
  call pGM_DIRECT_dip_dip_field(ind_dip,dip_field)
  if (eedmeth.ne.4) call pGM_SELF_dipole_field(numatoms,ind_dip,dip_field)
end subroutine pGM_NonBond_dip_dip_fields
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_ene_frc(numatoms,crd,x,nummols,molsiz,amass,ipairs, &
                       ntypes,iac,ico,cn1,cn2, & ! R. Luo: vdw
                       ene_perm,ene_ind,ene_vdw, &
                       ene_vdw_14,frc,pgm_virial,vdw_virial)
  use pol_gauss_recip, only : pGM_RECIP_ene_frc
  use pol_gauss_direct, only : pGM_DIRECT_ene_frc
  use pol_gauss_self, only : pGM_SELF_ene_frc
  use pol_gauss_induced, only : ind_dip
  use pol_gauss_mdin, only : pol_gauss_verbose
  use pol_gauss_multipoles, only : coulomb_const_kcal_per_mole

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: crd(3,*),x(*),amass(*)
  integer nummols,molsiz(*)
  integer,intent(in) :: ipairs(*)
  integer,intent(in) :: ntypes,iac(*),ico(*) ! R. Luo: vdw
  _REAL_,intent(in) :: cn1(*),cn2(*) ! R. Luo: vdw
  _REAL_,intent(out) :: ene_perm,ene_ind,ene_vdw,ene_vdw_14
  _REAL_,intent(inout) :: frc(3,numatoms),pgm_virial(3,3),vdw_virial(3,3)

# include "box.h"
# include "ew_erfc_spline.h"
#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
  logical master
#endif

  _REAL_ :: e_rec_perm,e_rec_ind,e_dir_perm,e_dir_ind,e_adj_perm, &
            e_adj_ind,e_self_perm,e_self_ind,e_dir_vdw,e_adj_vdw,e_rec_vdw
  _REAL_, allocatable :: pGM_force(:,:)
  _REAL_, allocatable :: phi(:,:), phi_rec(:,:), phi_dir(:,:), phi_self(:,:)
  _REAL_, allocatable :: displacement(:,:)
  _REAL_ :: mass, massmol, xmol, ymol, zmol
  integer :: imol, num, isave, iatom, i

  ! electric potentials and derivatives
  allocate(pGM_force(3,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(phi(10,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(phi_rec(10,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(phi_dir(10,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(phi_self(10,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(displacement(3,numatoms), stat=ier)
  REQUIRE(ier==0)

  e_dir_perm = 0.0d0
  e_adj_perm = 0.0d0
  e_rec_perm = 0.0d0
  e_dir_ind = 0.0d0
  e_adj_ind = 0.0d0
  e_rec_ind = 0.0d0
  e_dir_vdw = 0.0d0
  e_adj_vdw = 0.0d0
  e_rec_vdw = 0.0d0
  e_self_perm = 0.0d0
  e_self_ind = 0.0d0

  phi = 0.0d0
  phi_rec = 0.0d0
  phi_dir = 0.0d0
  phi_self = 0.0d0

#ifdef MPI
  call mpi_comm_size(recip_comm,numtasks,ierr)
  call mpi_comm_rank(recip_comm,mytaskid,ierr)
  master = mytaskid.eq.0
#endif

  if (eedmeth .ne. 4) call pGM_RECIP_ene_frc(numatoms,crd,x,ind_dip, &
                 e_rec_perm,e_rec_ind,frc,pgm_virial,phi_rec)

#ifdef MPI
  call mpi_comm_rank(commsander,mytaskid,ierr)
  call mpi_comm_size(commsander,numtasks,ierr)
  master = mytaskid.eq.0
#endif

  i = 0
  do imol = 1,nummols
    massmol = 0.d0
    xmol = 0.d0
    ymol = 0.d0
    zmol = 0.d0
    num = molsiz(imol)
    isave = i
     
    ! get c.o.m. of molecule - TODO: use precalced masses
     
    do iatom = 1,num
      i = i + 1
      mass = amass(i)
      massmol = massmol + mass
      xmol = xmol + mass*crd(1,i)
      ymol = ymol + mass*crd(2,i)
      zmol = zmol + mass*crd(3,i)
    end do
    xmol = xmol / massmol
    ymol = ymol / massmol
    zmol = zmol / massmol

    i = isave
    do iatom = 1,num
      i = i + 1
      displacement(1,i) = crd(1,i) - xmol
      displacement(2,i) = crd(2,i) - ymol
      displacement(3,i) = crd(3,i) - zmol
    end do
  end do

  call pGM_DIRECT_ene_frc(ipairs,ntypes,iac,ico,cn1,cn2,crd,x,ind_dip, &
                  e_dir_perm,e_dir_ind,e_dir_vdw,frc,pgm_virial,vdw_virial,phi_dir,numatoms,displacement)

  if (eedmeth .ne. 4) call pGM_SELF_ene_frc(numatoms,ind_dip,e_self_perm,e_self_ind,phi_self)

  phi = phi_rec + phi_dir + phi_self

  pGM_force = 0d0

  call pGM_energy_force_virial(numatoms,displacement,crd,nummols,molsiz,amass,&
                               phi,ind_dip,pGM_force,ene_ind,phi_rec,phi_dir,phi_self,pgm_virial)

  frc = frc + coulomb_const_kcal_per_mole * pGM_force
  ene_ind = coulomb_const_kcal_per_mole * ene_ind

  if ( pol_gauss_verbose == 2 )then
    write(6,'(a,/,4(1x,f14.4))') &
            ' e_rec_perm,e_dir_perm,e_adj_perm,e_self_perm = ', &
              e_rec_perm,e_dir_perm,e_adj_perm,e_self_perm
    write(6,'(a,/,4(1x,f14.4))') &
            ' e_rec_ind,e_dir_ind,e_adj_ind,e_self_ind = ', &
              e_rec_ind,e_dir_ind,e_adj_ind,e_self_ind
    write(6,'(a,/,3(1x,f14.4))') &
            ' e_dir_vdw,e_adj_vdw,e_rec_vdw = ', &
             e_dir_vdw,e_adj_vdw,e_rec_vdw
  endif

  deallocate(pGM_force, stat=ier)
  REQUIRE(ier==0)
  deallocate(phi, stat=ier)
  REQUIRE(ier==0)
  deallocate(phi_rec, stat=ier)
  REQUIRE(ier==0)
  deallocate(phi_dir, stat=ier)
  REQUIRE(ier==0)
  deallocate(phi_self, stat=ier)
  REQUIRE(ier==0)
  deallocate(displacement, stat=ier)
  REQUIRE(ier==0)

  ene_vdw = e_dir_vdw

end subroutine pGM_NonBond_ene_frc
!-------------------------------------------------------------------------------
subroutine pGM_energy_force_virial(numatoms,displacement,crd,nummols,molsiz,amass,&
                                   phi,induced_dipole,pGM_force,pGM_energy,phi_rec,phi_dir,phi_self,pgm_virial)
  use pol_gauss_multipoles, only : global_multipole,covalent_ptr,covalent_atm,covalent_dip
  use nblist, only: ucell,recip
  implicit none

# include "box.h"
# include "pol_gauss_mpole_index.h"

  integer :: numatoms
  _REAL_ :: crd(3,numatoms),amass(*),displacement(3,numatoms)
  integer nummols,molsiz(*)
  _REAL_ :: phi(10,numatoms),induced_dipole(3,numatoms), pGM_force(3,numatoms)
  _REAL_ :: phi_rec(10,numatoms), phi_dir(10,numatoms), phi_self(10,numatoms), pgm_virial(3,3)
  _REAL_ :: pGM_energy

  integer :: i, j, k, atom_first, atom_last, covalent_atom
  _REAL_ :: distance, deltr_dot_field, deltr_dot_field2!, deltr_dot_dipole

  _REAL_ :: virial_force(3,numatoms)

  ! Haixin: convert grad of phi to electric field and its derivatives
  phi(2:10,:) = -phi(2:10,:)
  phi_rec(2:10,:) = -phi_rec(2:10,:)
  phi_dir(2:10,:) = -phi_dir(2:10,:)
  phi_self(2:10,:) = -phi_self(2:10,:)

  virial_force = 0.0

  do i = 1, numatoms
    ! energy can be expressed as 1/2*q*phi and -1/2*mu*E
    pGM_energy = pGM_energy + 0.5d0*global_multipole(1,i)*phi(1,i)&
                            - 0.5d0*global_multipole(2,i)*phi(2,i)&
                            - 0.5d0*global_multipole(3,i)*phi(3,i)&
                            - 0.5d0*global_multipole(4,i)*phi(4,i)!&

    ! forces due to covalent monopoles
    pGM_force(1:3,i) = pGM_force(1:3,i) + global_multipole(1,i)*phi(2:4,i)

    ! forces due to total dipoles (covalent dipoles + induced dipoles)
    pGM_force(1,i) = pGM_force(1,i) + (global_multipole(2,i)+induced_dipole(1,i))*phi(Ind_200,i)&
                                    + (global_multipole(3,i)+induced_dipole(2,i))*phi(Ind_110,i)&
                                    + (global_multipole(4,i)+induced_dipole(3,i))*phi(Ind_101,i)
    pGM_force(2,i) = pGM_force(2,i) + (global_multipole(2,i)+induced_dipole(1,i))*phi(Ind_110,i)&
                                    + (global_multipole(3,i)+induced_dipole(2,i))*phi(Ind_020,i)&
                                    + (global_multipole(4,i)+induced_dipole(3,i))*phi(Ind_011,i)
    pGM_force(3,i) = pGM_force(3,i) + (global_multipole(2,i)+induced_dipole(1,i))*phi(Ind_101,i)&
                                    + (global_multipole(3,i)+induced_dipole(2,i))*phi(Ind_011,i)&
                                    + (global_multipole(4,i)+induced_dipole(3,i))*phi(Ind_002,i)

    ! forces due to changes in covalent dipoles 
    atom_last = covalent_ptr(i)
    atom_first = covalent_ptr(i-1) + 1
    do j = atom_first, atom_last
      covalent_atom = covalent_atm(j)

      distance = dot_product((crd(1:3,covalent_atom)-crd(1:3,i)),(crd(1:3,covalent_atom)-crd(1:3,i)))
      distance = sqrt(distance)
      if ( distance .gt. 3.0d0 ) then
        write(6,*) 'pGM Fatal Error: covalent bond too long!'
        call mexit(6,1)
      endif
      deltr_dot_field = dot_product((crd(1:3,covalent_atom) - crd(1:3,i)),phi(2:4,i))
      do k = covalent_ptr(covalent_atom-1)+1, covalent_ptr(covalent_atom)
        if (covalent_atm(k) == i) exit
      end do
      deltr_dot_field2 = dot_product((crd(1:3,i) - crd(1:3,covalent_atom)),phi(2:4,covalent_atom))

      pGM_force(1,i) = pGM_force(1,i) - &
                       covalent_dip(j)*(&
                                                                phi(2,i)/distance - &
                       (crd(1,covalent_atom) - crd(1,i))*deltr_dot_field/distance**3&
                       ) + &
                       covalent_dip(k)*(&
                                                     phi(2,covalent_atom)/distance - &
                       (crd(1,i) - crd(1,covalent_atom))*deltr_dot_field2/distance**3&
                       )
      pGM_force(2,i) = pGM_force(2,i) - &
                       covalent_dip(j)*(&
                                                                phi(3,i)/distance - &
                       (crd(2,covalent_atom) - crd(2,i))*deltr_dot_field/distance**3&
                       ) + &
                       covalent_dip(k)*(&
                                                     phi(3,covalent_atom)/distance - &
                       (crd(2,i) - crd(2,covalent_atom))*deltr_dot_field2/distance**3&
                       )
      pGM_force(3,i) = pGM_force(3,i) - &
                       covalent_dip(j)*(&
                                                                phi(4,i)/distance - &
                       (crd(3,covalent_atom) - crd(3,i))*deltr_dot_field/distance**3&
                       ) + &
                       covalent_dip(k)*(&
                                                     phi(4,covalent_atom)/distance - &
                       (crd(3,i) - crd(3,covalent_atom))*deltr_dot_field2/distance**3&
                       )

    end do
  end do

  if ( ntb > 0 ) then
!  phi = phi_dir + phi_self ! now the phi array is contaminated, this is only for virial calculation
!
!  do i = 1, numatoms
!
!    atom_last = covalent_ptr(i)
!    atom_first = covalent_ptr(i-1) + 1
!    do j = atom_first, atom_last
!      covalent_atom = covalent_atm(j)
!      distance = dot_product((crd(1:3,covalent_atom)-crd(1:3,i)),(crd(1:3,covalent_atom)-crd(1:3,i)))
!      distance = sqrt(distance)
!      if ( distance .gt. 3.0d0 ) then
!        write(6,*) 'pGM Fatal Error: covalent bond too long!'
!        call mexit(6,1)
!      endif
!      deltr_dot_field = dot_product((crd(1:3,covalent_atom) - crd(1:3,i)),phi(2:4,i))
!      do k = covalent_ptr(covalent_atom-1)+1, covalent_ptr(covalent_atom)
!        if (covalent_atm(k) == i) exit
!      end do
!      deltr_dot_field2 = dot_product((crd(1:3,i) - crd(1:3,covalent_atom)),phi(2:4,covalent_atom))
!
!      virial_force(1,i) = virial_force(1,i) - &
!                       covalent_dip(j)*(&
!                                                                phi(2,i)/distance - &
!                       (crd(1,covalent_atom) - crd(1,i))*deltr_dot_field/distance**3&
!                       ) + &
!                       covalent_dip(k)*(&
!                                                     phi(2,covalent_atom)/distance - &
!                       (crd(1,i) - crd(1,covalent_atom))*deltr_dot_field2/distance**3&
!                       )
!      virial_force(2,i) = virial_force(2,i) - &
!                       covalent_dip(j)*(&
!                                                                phi(3,i)/distance - &
!                       (crd(2,covalent_atom) - crd(2,i))*deltr_dot_field/distance**3&
!                       ) + &
!                       covalent_dip(k)*(&
!                                                     phi(3,covalent_atom)/distance - &
!                       (crd(2,i) - crd(2,covalent_atom))*deltr_dot_field2/distance**3&
!                       )
!      virial_force(3,i) = virial_force(3,i) - &
!                       covalent_dip(j)*(&
!                                                                phi(4,i)/distance - &
!                       (crd(3,covalent_atom) - crd(3,i))*deltr_dot_field/distance**3&
!                       ) + &
!                       covalent_dip(k)*(&
!                                                     phi(4,covalent_atom)/distance - &
!                       (crd(3,i) - crd(3,covalent_atom))*deltr_dot_field2/distance**3&
!                       )
!
!    end do
!  end do
!
!  do i = 1, numatoms
!
!    pgm_virial(1,1) = pgm_virial(1,1) - &
!                      (virial_force(1,i)*crd(1,i) - phi_rec(2,i)*induced_dipole(1,i))
!    pgm_virial(1,2) = pgm_virial(1,2) - &
!                      (virial_force(1,i)*crd(2,i) - phi_rec(2,i)*induced_dipole(2,i))
!    pgm_virial(1,3) = pgm_virial(1,3) - &
!                      (virial_force(1,i)*crd(3,i) - phi_rec(2,i)*induced_dipole(3,i))
!    pgm_virial(2,1) = pgm_virial(2,1) - &
!                      (virial_force(2,i)*crd(1,i) - phi_rec(3,i)*induced_dipole(1,i))
!    pgm_virial(2,2) = pgm_virial(2,2) - &
!                      (virial_force(2,i)*crd(2,i) - phi_rec(3,i)*induced_dipole(2,i))
!    pgm_virial(2,3) = pgm_virial(2,3) - &
!                      (virial_force(2,i)*crd(3,i) - phi_rec(3,i)*induced_dipole(3,i))
!    pgm_virial(3,1) = pgm_virial(3,1) - &
!                      (virial_force(3,i)*crd(1,i) - phi_rec(4,i)*induced_dipole(1,i))
!    pgm_virial(3,2) = pgm_virial(3,2) - &
!                      (virial_force(3,i)*crd(2,i) - phi_rec(4,i)*induced_dipole(2,i))
!    pgm_virial(3,3) = pgm_virial(3,3) - &
!                      (virial_force(3,i)*crd(3,i) - phi_rec(4,i)*induced_dipole(3,i))
!    ! this is a different way of collecting part of direct and self related
!    ! virial, different from below which is one obvious way of implimentation of
!    ! the formula in the paper. If interested, it can be a good exercise to
!    ! prove this two ways of collecting are equal.
!
!    atom_last = covalent_ptr(i)
!    atom_first = covalent_ptr(i-1) + 1
!    do j = atom_first, atom_last
!      covalent_atom = covalent_atm(j)
!      distance = dot_product((crd(1:3,covalent_atom)-crd(1:3,i)),(crd(1:3,covalent_atom)-crd(1:3,i)))
!      distance = sqrt(distance)
!      if ( distance .gt. 3.0d0 ) then
!        write(6,*) 'pGM Fatal Error: covalent bond too long!'
!        call mexit(6,1)
!      endif
!      deltr_dot_field = dot_product((crd(1:3,covalent_atom) - crd(1:3,i)),phi_rec(2:4,i))
!
!      pgm_virial(1,1) = pgm_virial(1,1) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(1,covalent_atom)-crd(1,i))*(crd(1,covalent_atom)-crd(1,i))
!      pgm_virial(1,2) = pgm_virial(1,2) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(1,covalent_atom)-crd(1,i))*(crd(2,covalent_atom)-crd(2,i))
!      pgm_virial(1,3) = pgm_virial(1,3) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(1,covalent_atom)-crd(1,i))*(crd(3,covalent_atom)-crd(3,i))
!      pgm_virial(2,1) = pgm_virial(2,1) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(2,covalent_atom)-crd(2,i))*(crd(1,covalent_atom)-crd(1,i))
!      pgm_virial(2,2) = pgm_virial(2,2) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(2,covalent_atom)-crd(2,i))*(crd(2,covalent_atom)-crd(2,i))
!      pgm_virial(2,3) = pgm_virial(2,3) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(2,covalent_atom)-crd(2,i))*(crd(3,covalent_atom)-crd(3,i))
!      pgm_virial(3,1) = pgm_virial(3,1) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(3,covalent_atom)-crd(3,i))*(crd(1,covalent_atom)-crd(1,i))
!      pgm_virial(3,2) = pgm_virial(3,2) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(3,covalent_atom)-crd(3,i))*(crd(2,covalent_atom)-crd(2,i))
!      pgm_virial(3,3) = pgm_virial(3,3) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(3,covalent_atom)-crd(3,i))*(crd(3,covalent_atom)-crd(3,i))
!    end do
!
!  end do

    do i = 1,numatoms

       pgm_virial(1,1) = pgm_virial(1,1) - &
                         (-phi_rec(2,i)*&
                          (induced_dipole(1,i)+global_multipole(2,i)+global_multipole(1,i)*displacement(1,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_200,i)*displacement(1,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_110,i)*displacement(1,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_101,i)*displacement(1,i)&
                         )

       pgm_virial(1,2) = pgm_virial(1,2) - &
                         (-phi_rec(2,i)*&
                          (induced_dipole(2,i)+global_multipole(3,i)+global_multipole(1,i)*displacement(2,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_200,i)*displacement(2,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_110,i)*displacement(2,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_101,i)*displacement(2,i)&
                         )

       pgm_virial(1,3) = pgm_virial(1,3) - &
                         (-phi_rec(2,i)*&
                          (induced_dipole(3,i)+global_multipole(4,i)+global_multipole(1,i)*displacement(3,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_200,i)*displacement(3,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_110,i)*displacement(3,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_101,i)*displacement(3,i)&
                         )

       pgm_virial(2,1) = pgm_virial(2,1) - &
                         (-phi_rec(3,i)*&
                          (induced_dipole(1,i)+global_multipole(2,i)+global_multipole(1,i)*displacement(1,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_110,i)*displacement(1,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_020,i)*displacement(1,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_011,i)*displacement(1,i)&
                         )

       pgm_virial(2,2) = pgm_virial(2,2) - &
                         (-phi_rec(3,i)*&
                          (induced_dipole(2,i)+global_multipole(3,i)+global_multipole(1,i)*displacement(2,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_110,i)*displacement(2,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_020,i)*displacement(2,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_011,i)*displacement(2,i)&
                         )

       pgm_virial(2,3) = pgm_virial(2,3) - &
                         (-phi_rec(3,i)*&
                          (induced_dipole(3,i)+global_multipole(4,i)+global_multipole(1,i)*displacement(3,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_110,i)*displacement(3,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_020,i)*displacement(3,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_011,i)*displacement(3,i)&
                         )

       pgm_virial(3,1) = pgm_virial(3,1) - &
                         (-phi_rec(4,i)*&
                          (induced_dipole(1,i)+global_multipole(2,i)+global_multipole(1,i)*displacement(1,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_101,i)*displacement(1,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_011,i)*displacement(1,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_002,i)*displacement(1,i)&
                         )

       pgm_virial(3,2) = pgm_virial(3,2) - &
                         (-phi_rec(4,i)*&
                          (induced_dipole(2,i)+global_multipole(3,i)+global_multipole(1,i)*displacement(2,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_101,i)*displacement(2,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_011,i)*displacement(2,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_002,i)*displacement(2,i)&
                         )

       pgm_virial(3,3) = pgm_virial(3,3) - &
                         (-phi_rec(4,i)*&
                          (induced_dipole(3,i)+global_multipole(4,i)+global_multipole(1,i)*displacement(3,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_101,i)*displacement(3,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_011,i)*displacement(3,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_002,i)*displacement(3,i)&
                         )

    end do
  end if

end subroutine pGM_energy_force_virial
!-------------------------------------------------------------
subroutine pGM_NONBOND_set_user_bit(do_recip,do_direct,do_self, &
                       do_vdw,do_induce)
  use pol_gauss_recip, only : pGM_RECIP_set_user_bit
  use pol_gauss_direct, only : pGM_DIRECT_set_user_bit
  use pol_gauss_self, only : pGM_SELF_set_user_bit
  use pol_gauss_induced, only : pGM_INDUCED_set_user_bit
  implicit none

  integer, intent(in) :: do_recip,do_direct,do_self, &
                         do_vdw,do_induce

  call pGM_RECIP_set_user_bit(do_recip)
  call pGM_DIRECT_set_user_bit(do_direct)
  call pGM_SELF_set_user_bit(do_self)
  call pGM_INDUCED_set_user_bit(do_induce)
end subroutine pGM_NONBOND_set_user_bit
!-------------------------------------------------------------------------------
subroutine POL_GAUSS_get_numlist(header,nf,num_list)
  implicit none

  character(len=*), intent(in) :: header
  integer, intent(in) :: nf
  integer, intent(out) :: num_list

  integer :: iok,ionerr
  character(len=80) :: fmt
  character(len=80) :: fmtin,dtype

  fmtin = '(10I8)'
  dtype = header//'NUM_LIST'
  ionerr = 1 ! not fatal if missing
  call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
  if ( iok == 0 )then !this data type found in prmtop
    read(nf,fmt)num_list
  else !either old style prmtop or data not found
    num_list = 0 ! upon return this will invalidate valence_term
  endif
end subroutine POL_GAUSS_get_numlist
!------------------------------------------------------------------------
subroutine POL_GAUSS_read_list_data(header,nf,dim1,num_list,list)
  implicit none

  character(len=*), intent(in) :: header
  integer, intent(in) :: nf,dim1,num_list
  integer, intent(out) :: list(dim1,num_list)

  integer :: iok,ionerr,j,k
  character(len=80) :: fmt
  character(len=80) :: fmtin,dtype

  ionerr = 0 !fatal if missing
  fmtin = '(10I8)'
  dtype = header//'LIST'
  call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
  read(nf,fmt)((list(j,k),j=1,dim1),k=1,num_list)
end subroutine POL_GAUSS_read_list_data
!----------------------------------------------------------
subroutine POL_GAUSS_read_real_list_data(header,nf,dim1,num_list,list)
  implicit none

  character(len=*), intent(in) :: header
  integer, intent(in) :: nf,dim1,num_list
  _REAL_, intent(out) :: list(dim1,num_list)

  integer :: iok,ionerr,j,k
  character(len=80) :: fmt
  character(len=80) :: fmtin,dtype

  ionerr = 0 !fatal if missing
  fmtin = '(5E16.8)'
  dtype = header//'LIST'
  call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
  read(nf,fmt)((list(j,k),j=1,dim1),k=1,num_list)
end subroutine POL_GAUSS_read_real_list_data
!--------------------------------------------------------------------
subroutine POL_GAUSS_read_real_scalar(flag,nf,scalar_value)
  implicit none

  character(len=*), intent(in) :: flag
  integer, intent(in) :: nf
  _REAL_, intent(out) :: scalar_value

  integer :: iok,ionerr
  character(len=80) :: fmt
  character(len=80) :: fmtin

  ionerr = 0 !fatal if missing
  fmtin = '(E16.8)'
  call nxtsec(nf,  6,  ionerr,fmtin,  flag,  fmt,  iok)
  read(nf,fmt)scalar_value
end subroutine POL_GAUSS_read_real_scalar
!--------------------------------------------------------------------
subroutine POL_GAUSS_get_startlist_endlist(num_list,startlist,endlist,siztask)
  implicit none

  integer,intent(in) :: num_list
  integer,intent(out) :: startlist,endlist,siztask

  integer piece
  integer numtasks,mytaskid

  numtasks=1  
  mytaskid=0
  if (numtasks > 1) then
    piece = num_list/numtasks
    startlist = mytaskid*piece+1
    endlist   = mytaskid*piece+piece
    if (mytaskid == (numtasks-1)) endlist = num_list
    siztask = endlist-startlist+1
  else
    startlist=1
    endlist = num_list
    siztask = endlist-startlist+1
  end if
end subroutine POL_GAUSS_get_startlist_endlist
!-------------------------------------------------------------
subroutine pGM_array_add(a,b,num)
  _REAL_,intent(inout) :: a(*)
  _REAL_,intent(in) :: b(*)
  integer,intent(in) :: num

  integer :: n
  do n = 1,num
    a(n) = a(n) + b(n)
  enddo
end subroutine pGM_array_add
!-------------------------------------------------------------
subroutine  pGM_dump_dipoles(dip,nsites,nf)
  implicit none

  _REAL_ dip(3,*)
  integer nsites,nf
  integer j,n

  write(nf,*) 'num of dipoles', nsites ! R. Luo
  do n = 1,nsites
    write(nf,'(3f20.12)')(dip(j,n),j=1,3)
  enddo
  call amflsh(nf)
end subroutine  pGM_dump_dipoles
