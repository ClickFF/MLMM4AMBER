#include "../include/dprec.fh"
#include "../include/assert.fh"

!-----------------------------------------------------------
module pol_gauss_direct
  implicit none
  private

# include "pol_gauss_mpole_index.h" 

  integer, save :: do_flag
  integer,save ::  num_pairs_in_ee_cut = 0,size_dipole_dipole_list = -1
  integer,save ::  num_pairs_in_ee_cut_short = 0,size_dipole_dipole_list_short = -1
  integer,save :: num_tensor, num_tensor_short
  _REAL_,save,allocatable :: dipole_dipole_tensor(:,:)
  _REAL_,save,allocatable :: dipole_dipole_tensor_short(:,:)
  integer,save,allocatable :: dipole_dipole_list(:,:)
  integer,save,allocatable :: dipole_dipole_list_short(:,:)
  _REAL_,parameter :: safety = 1.25d0, checklist = 0.9d0

  public  pGM_DIRECT_perm_field,pGM_DIRECT_dip_dip_field,pGM_DIRECT_ene_frc, &
          pGM_DIRECT_set_user_bit,pGM_DIRECT_dip_dip_field_short
# ifdef MPI
  public  pGM_DIRECT_bcast
# endif

 contains

#ifdef MPI
!-----------------------------------------------------------
subroutine pGM_DIRECT_bcast
  implicit none
  integer ierr

  include 'mpif.h'
# include "parallel.h"
  call mpi_bcast(do_flag,1,MPI_INTEGER,0,commsander,ierr)
end subroutine pGM_DIRECT_bcast
#endif
!-----------------------------------------------------------
subroutine pGM_DIRECT_set_user_bit(do_this)
  integer,intent(in) :: do_this
#include "do_flag.h"

  ! set the valid bit---this part always since no parmread needed
  do_flag = ibset(do_flag,VALID_BIT)

  if ( do_this == 1 )then ! do in all cases
    do_flag = ibset(do_flag,USER_INDUCE_BIT)
    do_flag = ibset(do_flag,USER_POSTINDUCE_BIT)
  elseif ( do_this == 2 )then ! do the induction, not the post-induction
    do_flag = ibset(do_flag,USER_INDUCE_BIT)
    do_flag = ibclr(do_flag,USER_POSTINDUCE_BIT)
  elseif ( do_this == 3 )then ! do the post-induction, not the induction
    do_flag = ibclr(do_flag,USER_INDUCE_BIT)
    do_flag = ibset(do_flag,USER_POSTINDUCE_BIT)
  elseif ( do_this == 0 )then 
    do_flag = ibclr(do_flag,USER_INDUCE_BIT)
    do_flag = ibclr(do_flag,USER_POSTINDUCE_BIT)
  else
    write(6,*)'pGM_DIRECT_set_user_bit: bad value of user do_this'
    call mexit(6,1)
  endif
end subroutine pGM_DIRECT_set_user_bit
!-----------------------------------------------------------
subroutine pGM_DIRECT_perm_field(ipairs,x,cart_dipole_field)
  use nblist, only: imagcrds,bckptr,nlogrid,nhigrid,numvdw,numhbnd,numexc,&
                    myindexlo,myindexhi,numimg,listtot_kp
  use pol_gauss_multipoles, only : global_multipole
  use pol_gauss_induced, only : is_polarizable
  use pol_gauss_mdin, only : ee_dsum_cut,ee_damped_cut,scf_local_cut

  integer,intent(in) :: ipairs(*)
  _REAL_,intent(in) :: x(*)
  _REAL_,intent(inout) :: cart_dipole_field(3,*)

# include "extra.h"
# include "def_time.h"
# include "ew_erfc_spline.h"
# include "ew_pme_recip.h"

  integer :: i,k,numpack,index,ncell_lo,ncell_hi,ntot,ier
  _REAL_ xk,yk,zk
!YDEBUG
  integer, save :: size_dipole_dipole_list_kp = 0
!ENDYDEBUG

  call timer_start(TIME_SHORT_ENE)

  ! check if dipole_dipole_list big enough
!YDEBUG
  if (eedmeth .eq. 4) then

  if (size_dipole_dipole_list_kp .ne. listtot_kp) then
    size_dipole_dipole_list = listtot_kp
    size_dipole_dipole_list_kp = listtot_kp
    size_dipole_dipole_list_short = listtot_kp
    if ( allocated(dipole_dipole_tensor) ) deallocate(dipole_dipole_tensor)
    if ( allocated(dipole_dipole_list) ) deallocate(dipole_dipole_list)
    if ( allocated(dipole_dipole_tensor_short) ) deallocate(dipole_dipole_tensor_short)
    if ( allocated(dipole_dipole_list_short) ) deallocate(dipole_dipole_list_short)
    if (master) then
      write(6,'(a,i10,i10)') &
            ' num_pairs_in_ee_cut,size_dipole_dipole_list = ', &
            num_pairs_in_ee_cut,size_dipole_dipole_list
      write(6,'(a,i10,i10)') &
            ' num_pairs_in_ee_cut_short,size_dipole_dipole_list_short = ', &
            num_pairs_in_ee_cut_short,size_dipole_dipole_list_short
    end if
    allocate(dipole_dipole_tensor(6,size_dipole_dipole_list), &
             dipole_dipole_list(2,size_dipole_dipole_list), &
             dipole_dipole_tensor_short(6,size_dipole_dipole_list_short), &
             dipole_dipole_list_short(2,size_dipole_dipole_list_short), &
             stat=ier)
    REQUIRE( ier==0 )
  end if

  else
!ENDYDEBUG

  if ( num_pairs_in_ee_cut > checklist*size_dipole_dipole_list )then
    if ( allocated(dipole_dipole_tensor) ) deallocate(dipole_dipole_tensor)
    if ( allocated(dipole_dipole_list) ) deallocate(dipole_dipole_list)
    call pGM_DIRECT_count_num_ee_pairs(ipairs,ee_dsum_cut,num_pairs_in_ee_cut)
    size_dipole_dipole_list = int(safety * num_pairs_in_ee_cut)
    if (master) then
      write(6,'(a,i10,i10)') &
              ' num_pairs_in_ee_cut,size_dipole_dipole_list = ', &
              num_pairs_in_ee_cut,size_dipole_dipole_list
    end if
    allocate(dipole_dipole_tensor(6,size_dipole_dipole_list),stat=ier)
    REQUIRE( ier==0 )
    allocate(dipole_dipole_list(2,size_dipole_dipole_list),stat=ier)
    REQUIRE( ier==0 )
  end if

  if ( num_pairs_in_ee_cut_short > checklist*size_dipole_dipole_list_short )then
    if ( allocated(dipole_dipole_tensor_short) ) deallocate(dipole_dipole_tensor_short)
    if ( allocated(dipole_dipole_list_short) ) deallocate(dipole_dipole_list_short)
    call pGM_DIRECT_count_num_ee_pairs(ipairs,scf_local_cut,num_pairs_in_ee_cut_short)
    size_dipole_dipole_list_short = int(safety * num_pairs_in_ee_cut_short)
    if (master) then
      write(6,'(a,i10,i10)') &
              ' num_pairs_in_ee_cut_short,size_dipole_dipole_list_short = ', &
              num_pairs_in_ee_cut_short,size_dipole_dipole_list_short
    end if
    allocate(dipole_dipole_tensor_short(6,size_dipole_dipole_list_short),stat=ier)
    REQUIRE( ier==0 )
    allocate(dipole_dipole_list_short(2,size_dipole_dipole_list_short),stat=ier)
    REQUIRE( ier==0 )
    dipole_dipole_tensor_short = 0d0
  end if

!YDEBUG
  end if
!ENDYDEBUG

  numpack = 1
  num_tensor = 0
  num_tensor_short = 0
  do index = myindexlo,myindexhi
    if ( numimg(index) > 0 )then
      ncell_lo = nlogrid(index)
      ncell_hi = nhigrid(index)
      do k = ncell_lo,ncell_hi
        i = bckptr(k)
        xk = imagcrds(1,k)
        yk = imagcrds(2,k)
        zk = imagcrds(3,k)
        ntot = numvdw(i) + numhbnd(i) + numexc(i)
        if ( ntot > 0 )then
          call pGM_DIRECT_perm_field_i(i,ipairs(numpack),ntot, &
                   xk,yk,zk,ew_coeff,eedtbdns,x(leed_cub),x(leed_lin), &
                   ee_type,eedmeth,dxdr,ee_dsum_cut, &
                   ee_damped_cut,scf_local_cut, &
                   is_polarizable, &
                   dipole_dipole_tensor,dipole_dipole_list, &
                   size_dipole_dipole_list,num_tensor, &
                   dipole_dipole_tensor_short,dipole_dipole_list_short, &
                   size_dipole_dipole_list_short,num_tensor_short, &
                   global_multipole,cart_dipole_field)
          numpack = numpack + ntot
        end if  ! ( ntot > 0 )
      end do  !  k = ncell_lo,ncell_hi
    end if  ! ( numimg(k) > 0 )
  end do  !  index = myindexlo,myindexhi
  call timer_stop(TIME_SHORT_ENE)
  return
end subroutine pGM_DIRECT_perm_field
!-------------------------------------------------------
subroutine pGM_DIRECT_ene_frc(ipairs,ntypes,iac,ico,cn1,cn2,crd,x,ind_dip, &
               ene_perm,ene_ind,ene_vdw,frc,virial,vdw_virial,phi,numatoms,displacement)

  use nblist, only: imagcrds,bckptr,nlogrid,nhigrid,numvdw,numhbnd,numexc, &
                    myindexlo,myindexhi,numimg
  use pol_gauss_multipoles, only : global_multipole
  use pol_gauss_induced, only : is_polarizable
  use pol_gauss_mdin, only : ee_dsum_cut,ee_damped_cut

  integer,intent(in) :: ipairs(*)
  integer,intent(in) :: ntypes,iac(*),ico(*)
  _REAL_,intent(in) :: cn1(*),cn2(*)
  _REAL_,intent(in) :: crd(3,*),x(*)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(inout) :: ene_perm,ene_ind,ene_vdw,frc(3,*),virial(3,3),vdw_virial(3,3)
  integer :: numatoms
  _REAL_ :: phi(10,numatoms), displacement(3,numatoms)

# include "parallel.h"
# include "def_time.h"
# include "ew_erfc_spline.h"
# include "ew_pme_recip.h"

  integer :: i,k,numpack,index,ncell_lo,ncell_hi,ntot
  _REAL_ xk,yk,zk

  numpack = 1
  ene_perm = 0.d0
  ene_ind = 0.d0
  ene_vdw = 0.d0
  call timer_start(TIME_SHORT_ENE)
  do index = myindexlo,myindexhi
    if ( numimg(index) > 0 )then
      ncell_lo = nlogrid(index)
      ncell_hi = nhigrid(index)
      do k = ncell_lo,ncell_hi
        i = bckptr(k)
        xk = imagcrds(1,k)
        yk = imagcrds(2,k)
        zk = imagcrds(3,k)
        ntot = numvdw(i) + numhbnd(i) + numexc(i)
        if ( ntot > 0 )then
          call pGM_DIRECT_ene_force_i(i,ipairs(numpack),ntot,numvdw(i),xk,yk,zk, &
                   ew_coeff,eedtbdns, &
                   x(leed_cub),x(leed_lin), &
                   ee_type,eedmeth,dxdr,ee_dsum_cut, &
                   ee_damped_cut, &
                   is_polarizable, &
                   ind_dip,global_multipole,&
                   ntypes,iac,ico,cn1,cn2,&
                   ene_perm,ene_ind,ene_vdw,frc,virial,vdw_virial,phi,numatoms,displacement)
          numpack = numpack + ntot
        end if  ! ( ntot > 0 )
      end do  !  k = ncell_lo,ncell_hi
    end if  ! ( numimg(k) > 0 )
  end do  !  index = myindexlo,myindexhi
  call timer_stop(TIME_SHORT_ENE)
  return
end subroutine pGM_DIRECT_ene_frc
!-------------------------------------------------------
subroutine pGM_DIRECT_count_num_ee_pairs(ipairs,ee_dsum_cut,num_pairs_in_ee_cut)
  use nblist, only: imagcrds,bckptr,nlogrid,nhigrid,numvdw,numhbnd,numexc, &
                    myindexlo,myindexhi,numimg
  integer, intent(in) :: ipairs(*)
  _REAL_,intent(in) :: ee_dsum_cut
  integer,intent(out) :: num_pairs_in_ee_cut

  integer :: i,k,numpack,index,ncell_lo,ncell_hi,ntot
  _REAL_ xk,yk,zk

  numpack = 1
  num_pairs_in_ee_cut = 0
  do index = myindexlo,myindexhi
    if ( numimg(index) > 0 )then
      ncell_lo = nlogrid(index)
      ncell_hi = nhigrid(index)
      do k = ncell_lo,ncell_hi
        i = bckptr(k)
        xk = imagcrds(1,k)
        yk = imagcrds(2,k)
        zk = imagcrds(3,k)
        ntot = numvdw(i) + numhbnd(i) + numexc(i)
        if ( ntot > 0 )then
          call pGM_DIRECT_increment_ee_pairs( &
                   ipairs(numpack),ntot, &
                   xk,yk,zk,ee_dsum_cut,  &
                   num_pairs_in_ee_cut)
          numpack = numpack + ntot
        end if  ! ( ntot > 0 )
      end do  !  k = ncell_lo,ncell_hi
    end if  ! ( numimg(k) > 0 )
  end do  !  index = myindexlo,myindexhi
end subroutine pGM_DIRECT_count_num_ee_pairs
!-------------------------------------------------------
subroutine pGM_DIRECT_increment_ee_pairs(ipairs,numtot,xk,yk,zk,ee_dsum_cut,  &
               num_pairs_in_ee_cut)
  use nblist, only: imagcrds,tranvec

  integer,intent(in) :: ipairs(*),numtot
  _REAL_,intent(in) :: xk,yk,zk,ee_dsum_cut
  integer,intent(inout) :: num_pairs_in_ee_cut

  _REAL_ :: ee_dsum_cut2
  _REAL_ :: xktran(3,18)
  integer :: mask27
  integer :: m,np,itran
  _REAL_ :: delx,dely,delz,delr2

  mask27 = 2**27 - 1
  ee_dsum_cut2 = ee_dsum_cut*ee_dsum_cut
  do m=1,18
    xktran(1,m) = tranvec(1,m) - xk
    xktran(2,m) = tranvec(2,m) - yk
    xktran(3,m) = tranvec(3,m) - zk
  end do
  do m = 1,numtot
    np=ipairs(m)
    itran=ishft(np,-27)
    np = iand(np,mask27)
    delx = imagcrds(1,np) + xktran(1,itran)
    dely = imagcrds(2,np) + xktran(2,itran)
    delz = imagcrds(3,np) + xktran(3,itran)
    delr2 = delx*delx + dely*dely+delz*delz
    if ( delr2 < ee_dsum_cut2 )then
      num_pairs_in_ee_cut = num_pairs_in_ee_cut + 1
    end if
  end do
end subroutine pGM_DIRECT_increment_ee_pairs
!-------------------------------------------------------
subroutine pGM_DIRECT_perm_field_i(i,ipairs,numtot,xk,yk,zk, &
               ewaldcof,eedtbdns,eed_cub,eed_lin, &
               ee_type,eedmeth,beta,ee_dsum_cut, &
               ee_damped_cut,scf_local_cut, &
               is_polarizable, &
               dipole_dipole_tensor,dipole_dipole_list, &
               size_dipole_dipole_list, num_tensor, &
               dipole_dipole_tensor_short,dipole_dipole_list_short, &
               size_dipole_dipole_list_short,num_tensor_short, &
               global_multipole,gradphi)
  use nblist, only: bckptr,imagcrds,tranvec
  use constants, only : zero,third,half,one,two,three,five
  use pol_gauss_multipoles, only : r_gauss ! R. Luo

  integer,intent(in) :: i,ipairs(*),numtot
  _REAL_,intent(in) :: xk,yk,zk,ewaldcof,eedtbdns,eed_cub(4,*),eed_lin(2,*)
  integer,intent(in) :: ee_type,eedmeth
  _REAL_,intent(in) :: beta,ee_dsum_cut,ee_damped_cut,scf_local_cut
  logical,intent(in) :: is_polarizable(*)
  _REAL_,intent(out) :: dipole_dipole_tensor(6,*)
  integer,intent(out) :: dipole_dipole_list(2,*)
  integer,intent(in) :: size_dipole_dipole_list
  integer,intent(inout) :: num_tensor
  _REAL_,intent(out) :: dipole_dipole_tensor_short(6,*)
  integer,intent(out) :: dipole_dipole_list_short(2,*)
  integer,intent(in) :: size_dipole_dipole_list_short
  integer,intent(inout) :: num_tensor_short
  _REAL_,intent(in) :: global_multipole(10,*)
  _REAL_,intent(inout) :: gradphi(3,*)

  _REAL_ :: ee_dsum_cut2,ee_damped_cut2,scf_local_cut2
  _REAL_ :: xktran(3,18)
  integer :: mask27
  integer :: m,n,np,itran,j,ind
  _REAL_ :: delx,dely,delz,delr2,delr,delr2inv,x,dx
  _REAL_ :: switch,d_switch_dx
  _REAL_ :: B(0:3),fac,fact,del,gphi_i(3),gphi_j(3)
  _REAL_ :: Rn(1),Rn_1(4),Rn_2(10),Rn_3(20)

  ! R. Luo: gaussian variables
  _REAL_ :: betaij, facij
  _REAL_ :: switchij,d_switchij_dx,y
  _REAL_ :: Bij(0:3)

# include "do_flag.h"

  if ( iand(do_flag,PROCEED_INDUCE) /= PROCEED_INDUCE ) return

  mask27 = 2**27 - 1
  del = one / eedtbdns
  ee_dsum_cut2 = ee_dsum_cut*ee_dsum_cut
  ee_damped_cut2 = ee_damped_cut*ee_damped_cut
  scf_local_cut2 = scf_local_cut*scf_local_cut
 
  do m = 1, 18
    xktran(1,m) = tranvec(1,m) - xk
    xktran(2,m) = tranvec(2,m) - yk
    xktran(3,m) = tranvec(3,m) - zk
  end do
  do m = 1, numtot
    np = ipairs(m)
    itran = ishft(np,-27)
    np = iand(np,mask27)
    j = bckptr(np)
    delx = imagcrds(1,np) + xktran(1,itran)
    dely = imagcrds(2,np) + xktran(2,itran)
    delz = imagcrds(3,np) + xktran(3,itran)
    delr2 = delx*delx + dely*dely+delz*delz

    if ( (delr2 < ee_dsum_cut2) ) then
      !---------------------------------------------------------
      ! McMurchie-Davidson recursion for interaction tensor:
      ! interaction at point charge level given by complementary Boys
      ! B(0) = int_0^1 exp(-pr^2^t^2)dt
      ! complementary Boys is BC(0) = 1/r - B(0)
      ! (d/dr) B(0) = (-2p)*r int_0^1 t^2*exp(-pr^2t^2)dt = (-2p)*r B(1)
      ! and so if R(0,0,0,n) = (-2p)^n B(n) then we have
      ! (1/r)*(d/dr) R(0,0,0,n) = R(0,0,0,n+1)
      ! Now let R(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
      ! Then e.g. R(t+1,u,v,n) = t*R(t-1,u,v,n+1) + x*R(t,u,v,n)
      ! proof:  
      ! R(t+1,u,v,n) = (d/dx)^t+1 (d/dy)^u (d/dz)^v R(0,0,0,n)
      !              = (d/dx)^t (d/dy)^u (d/dz)^v (d/dx) R(0,0,0,n)
      !              = (d/dx)^t (d/dy)^u (d/dz)^v x*(1/r)(d/dr)R(0,0,0,n)
      !              = (d/dx)^t (d/dy)^u (d/dz)^v [x*R(0,0,0,n+1)]
      !              = (d/dx)^t [x*R(0,u,v,n+1)]
      !              = t*R(t-1,u,v,n+1) + x*R(t,u,v,n+1) (Leibniz)
      ! similar recursions hold for R(t,u+1,v,n),R(t,u,v+1,n)
      ! R(t,u+1,v,n) = u*R(t,u-1,v,n+1) + y*R(t,u,v,n+1)
      ! R(t,u,v+1,n) = v*R(t,u,v-1,n+1) + z*R(t,u,v,n+1)
      ! below array is packed---hence use of Ind_tuv
      ! Rn(Ind_tuv) denotes R(t,u,v,n)
      ! Due to its form--we recur downwards in n
      !---------------------------------------------------------
      ! top n is 3 for dipole fields
      ! get boys and R(0,0,0,n), n=0,1,2,3
      if ( eedmeth /= 1 .and. eedmeth /=4 ) then
        write(6,'(/a)') ' pGM_DIRECT_perm_field_i: Only eedmeth=1/4 are supported'
        call mexit(6,1)
      end if

      delr = sqrt(delr2)
      delr2inv = one/delr2

      ! for free-space coulombic calculation
      if ( eedmeth == 4 ) then
        ! r_gauss and r_ewald are saved as 2*radi^2 for speed
        betaij  = one/sqrt(r_gauss(i) + r_gauss(j))
        y = betaij*delr

        !if ( x >= ee_gauss_cut ) then ! too far away for gaussian
        !   switch = one
        !   d_switch_dx = zero

        !   B(0) = switch*delr*delr2inv ! 1/r
        !   fact = d_switch_dx*beta ! zero
        !   B(1) = (B(0) - fact)*delr2inv ! 1/r^3
        !   fact = fac*fact ! zero
        !   B(2) = (three*B(1) - fact)*delr2inv ! 3/r^5
        !end if

        !if ( x < ee_gauss_cut ) then ! gaussian function coulomb potential
        facij = two*betaij*betaij
        ind = int(eedtbdns*y) + 1
        dx = y - (dble(ind)-one)*del
        switchij = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
                   dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*third)*half)
        d_switchij_dx = eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                   dx*eed_cub(4,ind)*half)

        !call get_ee_func(x,switchij,d_switchij_dx,1)

        Bij(0) = (one - switchij)*delr*delr2inv ! erf/r
        fact = -d_switchij_dx*betaij
        Bij(1) = (Bij(0) - fact)*delr2inv
        fact = facij*fact
        Bij(2) = (three*Bij(1) - fact)*delr2inv
        !end if

        B(0:2) = Bij(0:2)
      end if ! eedmeth == 4

      ! for PME calculation, only do the cubic spline option
      if ( eedmeth == 1 ) then
        ! compute erfc=one-erf function and its derivative for PME gaussians.
        x = beta*delr
        fac = two*ewaldcof*ewaldcof
        ind = int(eedtbdns*x) + 1
        dx = x - (dble(ind)-one)*del
        switch = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
                 dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*third)*half)
        d_switch_dx = eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                 dx*eed_cub(4,ind)*half)

        !call get_ee_func(x,switch,d_switch_dx,1)

        ! compute erfc=one-erf function and its derivative for pGM gaussians.
        betaij  = one/sqrt(r_gauss(i) + r_gauss(j))
        y = betaij*delr
        !if ( y < ee_gauss_cut ) then 
        facij = two*betaij*betaij
        ind = int(eedtbdns*y) + 1
        dx = y - (dble(ind)-one)*del
        switchij = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
                   dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*third)*half)
        d_switchij_dx = eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                   dx*eed_cub(4,ind)*half)

        !call get_ee_func(y,switchij,d_switchij_dx,1)

        !end if

        ! compute Boys series
        ! TD: Got the idea for B_l from Walter Smith's CCP5 article 1982
        ! Ewald for point multipoles
        ! B_l satisfies grad_i B_l(|r_j - r_i|) = (r_j - r_i)B_{l+1}(|r_j-r_i|)
        ! grad_j B_l(|r_j - r_i|) = -grad_i B_l(|r_j - r_i|)
        !if ( y >= ee_gauss_cut) ) then ! too far away for gaussian
        !  B(0) = switch*delr*delr2inv
        !  fact = d_switch_dx*beta
        !  B(1) = (B(0) - fact)*delr2inv
        !  fact = fac*fact
        !  B(2) = (three*B(1) - fact)*delr2inv
        !end if
        !if ( y < ee_gauss_cut ) then ! gaussian function coulomb potential
        B(0) = (one - switch)*delr*delr2inv
        fact = -d_switch_dx*beta
        B(1) = (B(0) - fact)*delr2inv
        fact = fac*fact
        B(2) = (three*B(1) - fact)*delr2inv

        Bij(0) = (one - switchij)*delr*delr2inv
        fact = -d_switchij_dx*betaij
        Bij(1) = (Bij(0) - fact)*delr2inv
        fact = facij*fact
        Bij(2) = (three*Bij(1) - fact)*delr2inv

        ! erf_pgm/r - erf_pme/r
        B(0) = Bij(0) - B(0)
        B(1) = Bij(1) - B(1)
        B(2) = Bij(2) - B(2)
        !end if
      end if ! eedmeth == 1

      ! negate the odd order boys factors
      B(1) = -B(1)

      ! multipdelz*Rn_1(Ind_000)ole interaction tensor elements
      Rn_1(Ind_000) = B(2)
      Rn_2(Ind_000) = B(1)
      Rn_2(Ind_100) = delx*Rn_1(Ind_000)
      Rn_2(Ind_010) = dely*Rn_1(Ind_000)
      Rn_2(Ind_001) = delz*Rn_1(Ind_000)
      Rn_3(Ind_100) = delx*Rn_2(Ind_000)
      Rn_3(Ind_010) = dely*Rn_2(Ind_000)
      Rn_3(Ind_001) = delz*Rn_2(Ind_000)
      Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
      Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
      Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
      Rn_3(Ind_110) = delx*Rn_2(Ind_010)
      Rn_3(Ind_101) = delx*Rn_2(Ind_001)
      Rn_3(Ind_011) = dely*Rn_2(Ind_001)

      ! store the dipole-dipole component
      num_tensor = num_tensor+1
      if ( num_tensor > size_dipole_dipole_list )then
        write(6,'(/a)')' pGM_DIRECT_perm_field_i: Too many dipole_dipole pairs'
        call mexit(6,1)
      end if
      dipole_dipole_list(1,num_tensor) = i
      dipole_dipole_list(2,num_tensor) = j
      dipole_dipole_tensor(1,num_tensor) = Rn_3(Ind_200)
      dipole_dipole_tensor(2,num_tensor) = Rn_3(Ind_110)
      dipole_dipole_tensor(3,num_tensor) = Rn_3(Ind_101)
      dipole_dipole_tensor(4,num_tensor) = Rn_3(Ind_020)
      dipole_dipole_tensor(5,num_tensor) = Rn_3(Ind_011)
      dipole_dipole_tensor(6,num_tensor) = Rn_3(Ind_002)

      ! store the short dipole-dipole component
      if ( delr2 <= scf_local_cut2 ) then
      num_tensor_short = num_tensor_short+1
      if ( num_tensor_short > size_dipole_dipole_list_short )then
        write(6,'(/a)')' pGM_DIRECT_perm_field_i: Too many short dipole_dipole pairs'
        call mexit(6,1)
      end if
      dipole_dipole_list_short(1,num_tensor_short) = i
      dipole_dipole_list_short(2,num_tensor_short) = j
      dipole_dipole_tensor_short(1,num_tensor_short) = -Bij(1) + delx*delx*Bij(2) != Rn_3(Ind_200)
      dipole_dipole_tensor_short(2,num_tensor_short) = delx*dely*Bij(2)           != Rn_3(Ind_110)
      dipole_dipole_tensor_short(3,num_tensor_short) = delx*delz*Bij(2)           != Rn_3(Ind_101)
      dipole_dipole_tensor_short(4,num_tensor_short) = -Bij(1) + dely*dely*Bij(2) != Rn_3(Ind_020)
      dipole_dipole_tensor_short(5,num_tensor_short) = dely*delz*Bij(2)           != Rn_3(Ind_011)
      dipole_dipole_tensor_short(6,num_tensor_short) = -Bij(1) + delz*delz*Bij(2) != Rn_3(Ind_002)
      end if

      ! next compute the grad of phi, not field, so pay attention to sign
      gphi_i(1) = Rn_3(Ind_100)*global_multipole(Ind_000,j) + &
                  Rn_3(Ind_200)*global_multipole(Ind_100,j) + &
                  Rn_3(Ind_110)*global_multipole(Ind_010,j) + &
                  Rn_3(Ind_101)*global_multipole(Ind_001,j)
      gphi_i(2) = Rn_3(Ind_010)*global_multipole(Ind_000,j) + &
                  Rn_3(Ind_110)*global_multipole(Ind_100,j) + &
                  Rn_3(Ind_020)*global_multipole(Ind_010,j) + &
                  Rn_3(Ind_011)*global_multipole(Ind_001,j)
      gphi_i(3) = Rn_3(Ind_001)*global_multipole(Ind_000,j) + &
                  Rn_3(Ind_101)*global_multipole(Ind_100,j) + &
                  Rn_3(Ind_011)*global_multipole(Ind_010,j) + &
                  Rn_3(Ind_002)*global_multipole(Ind_001,j)
      gradphi(1,i) = gradphi(1,i) - gphi_i(1) ! R. Luo:
      gradphi(2,i) = gradphi(2,i) - gphi_i(2) ! We are returning to the convention
      gradphi(3,i) = gradphi(3,i) - gphi_i(3) ! gradphi is grad of phi, not field

      gphi_j(1) = -Rn_3(Ind_100)*global_multipole(Ind_000,i) + &
                   Rn_3(Ind_200)*global_multipole(Ind_100,i) + &
                   Rn_3(Ind_110)*global_multipole(Ind_010,i) + &
                   Rn_3(Ind_101)*global_multipole(Ind_001,i)
      gphi_j(2) = -Rn_3(Ind_010)*global_multipole(Ind_000,i) + &
                   Rn_3(Ind_110)*global_multipole(Ind_100,i) + &
                   Rn_3(Ind_020)*global_multipole(Ind_010,i) + &
                   Rn_3(Ind_011)*global_multipole(Ind_001,i)
      gphi_j(3) = -Rn_3(Ind_001)*global_multipole(Ind_000,i) + &
                   Rn_3(Ind_101)*global_multipole(Ind_100,i) + &
                   Rn_3(Ind_011)*global_multipole(Ind_010,i) + &
                   Rn_3(Ind_002)*global_multipole(Ind_001,i)
      gradphi(1,j) = gradphi(1,j) - gphi_j(1) ! R.Luo
      gradphi(2,j) = gradphi(2,j) - gphi_j(2) ! we are returning to the convention
      gradphi(3,j) = gradphi(3,j) - gphi_j(3) ! gradphi is grad of phi, not field
    end if !( delr2 < ee_dsum_cut2 )then
  end do !m = 1,numlist
end subroutine pGM_DIRECT_perm_field_i
!-------------------------------------------------------
subroutine pGM_DIRECT_dip_dip_field(ind_dip,dip_field)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(inout) :: dip_field(3,*)
# include "do_flag.h"

  if ( iand(do_flag,PROCEED_INDUCE) /= PROCEED_INDUCE )return

  call timer_start(TIME_SHORT_ENE)
  call pGM_DIRECT_calc_dipdip_field(num_tensor, &
           dipole_dipole_list,dipole_dipole_tensor, &
           ind_dip,dip_field)
  call timer_stop(TIME_SHORT_ENE)
end subroutine pGM_DIRECT_dip_dip_field
!------------------------------------------------------- 
subroutine pGM_DIRECT_dip_dip_field_short(ind_dip,dip_field)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(inout) :: dip_field(3,*)
# include "do_flag.h"

  integer n, i, j

  if ( iand(do_flag,PROCEED_INDUCE) /= PROCEED_INDUCE )return

  call timer_start(TIME_SHORT_ENE)
  call pGM_DIRECT_calc_dipdip_field(num_tensor_short, &
           dipole_dipole_list_short,dipole_dipole_tensor_short, &
           ind_dip,dip_field)
  call timer_stop(TIME_SHORT_ENE)
end subroutine pGM_DIRECT_dip_dip_field_short
!-------------------------------------------------------
subroutine pGM_DIRECT_calc_dipdip_field(num_tensor, &
               dipole_dipole_list,dipole_dipole_tensor, &
               ind_dip,dip_field)
  integer,intent(in) :: num_tensor,dipole_dipole_list(2,*)
  _REAL_,intent(in) :: dipole_dipole_tensor(6,*)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(inout) :: dip_field(3,*)

  integer :: i,j,n

  ! again we are computing grad of phi, not field
  ! minus signs due to deriv wrt position of i
  do n = 1, num_tensor
    i = dipole_dipole_list(1,n)
    j = dipole_dipole_list(2,n)
    dip_field(1,i) = dip_field(1,i) -  &
         ( dipole_dipole_tensor(1,n)*ind_dip(1,j) + &
           dipole_dipole_tensor(2,n)*ind_dip(2,j) + &
           dipole_dipole_tensor(3,n)*ind_dip(3,j) )
    dip_field(2,i) = dip_field(2,i) -  &
         ( dipole_dipole_tensor(2,n)*ind_dip(1,j) + &
           dipole_dipole_tensor(4,n)*ind_dip(2,j) + &
           dipole_dipole_tensor(5,n)*ind_dip(3,j) )
    dip_field(3,i) = dip_field(3,i) -  &
         ( dipole_dipole_tensor(3,n)*ind_dip(1,j) + &
           dipole_dipole_tensor(5,n)*ind_dip(2,j) + &
           dipole_dipole_tensor(6,n)*ind_dip(3,j) )
    dip_field(1,j) = dip_field(1,j) -  &
         ( dipole_dipole_tensor(1,n)*ind_dip(1,i) + &
           dipole_dipole_tensor(2,n)*ind_dip(2,i) + &
           dipole_dipole_tensor(3,n)*ind_dip(3,i) )
    dip_field(2,j) = dip_field(2,j) -  &
         ( dipole_dipole_tensor(2,n)*ind_dip(1,i) + &
           dipole_dipole_tensor(4,n)*ind_dip(2,i) + &
           dipole_dipole_tensor(5,n)*ind_dip(3,i) )
    dip_field(3,j) = dip_field(3,j) -  &
         ( dipole_dipole_tensor(3,n)*ind_dip(1,i) + &
           dipole_dipole_tensor(5,n)*ind_dip(2,i) + &
           dipole_dipole_tensor(6,n)*ind_dip(3,i) )
   end do
end subroutine pGM_DIRECT_calc_dipdip_field
!-------------------------------------------------------
subroutine pGM_DIRECT_ene_force_i(i,ipairs,numtot,numvdw,xk,yk,zk, & 
               ewaldcof,eedtbdns,eed_cub,eed_lin, &
               ee_type,eedmeth,dxdr,ee_dsum_cut, &
               ee_damped_cut, &
               is_polarizable, &
               ind_dip,global_multipole,&
               ntypes,iac,ico,cn1,cn2,&
               ene_perm,ene_ind,ene_vdw,frc,virial,vdw_virial,phi,numatoms,displacement) 
  use nblist, only: bckptr,imagcrds,tranvec,cutoffnb
  use pol_gauss_multipoles, only : coulomb_const_kcal_per_mole,r_gauss
  use constants, only : zero,third,half,one,two,three,five,seven,nine

  integer,intent(in) :: i,ipairs(*),numtot,numvdw
  _REAL_,intent(in) :: xk,yk,zk,ewaldcof,eedtbdns,eed_cub(4,*),eed_lin(2,*)
  integer,intent(in) :: ee_type,eedmeth
  _REAL_,intent(in) :: dxdr,ee_dsum_cut,ee_damped_cut
  logical,intent(in) :: is_polarizable(*)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(in) :: global_multipole(10,*)
  integer,intent(in) :: ntypes,iac(*),ico(*) 
  _REAL_,intent(in) :: cn1(*),cn2(*) 
  _REAL_,intent(inout) :: ene_perm,ene_ind,ene_vdw,frc(3,*),virial(3,3),vdw_virial(3,3)
  integer,intent(in) :: numatoms
  _REAL_,intent(inout) :: phi(10,numatoms),displacement(3,numatoms)

# include "box.h"
# include "do_flag.h"

  ! local variables
  _REAL_ :: cutoffnb2
  _REAL_ :: ee_dsum_cut2,ee_damped_cut2
  _REAL_ :: xktran(3,18)
  integer :: mask27
  integer :: m,n,np,itran,j,ind,jj
  _REAL_ :: delx,dely,delz,delr2,delr,delr2inv,x,dx,switch,d_switch_dx
  integer :: ic,iaci 
  _REAL_ :: delr6inv,delr12inv,f6,f12,df
  _REAL_ :: dfx, dfy, dfz
  _REAL_ :: B(0:5),BD(5),fac,fact,del,gmj(10),gmi(10),tmi(10),tmj(10),phi_i(10),phi_j(10)
  _REAL_ :: e_pp,e_ind,g_pp(3),g_ind(3),  &
            i_di(3),i_pi(3),i_dj(3),i_pj(3),i_mi(3),i_mj(3)
  _REAL_ :: Rn(1),Rn_1(4),Rn_2(10),Rn_3(20),Rn_4(35),Rn_5(56)
  _REAL_,parameter :: const1 = 0.6d0,const2 = 18.d0/35.d0,const3 = 9.d0/35.d0

  ! R. Luo: gaussian variables
  _REAL_ :: beta, betaij, facij
  _REAL_ :: switchij,d_switchij_dx,y
  _REAL_ :: Bij(0:3)

  if ( iand(do_flag,PROCEED_POSTINDUCE) /= PROCEED_POSTINDUCE )return

  beta = dxdr

  mask27 = 2**27 - 1
  fac = two*ewaldcof*ewaldcof
  del = one / eedtbdns
  ee_dsum_cut2 = ee_dsum_cut*ee_dsum_cut
  cutoffnb2 = cutoffnb*cutoffnb  
  ee_damped_cut2 = ee_damped_cut*ee_damped_cut
   
  do m=1,18
    xktran(1,m) = tranvec(1,m) - xk
    xktran(2,m) = tranvec(2,m) - yk
    xktran(3,m) = tranvec(3,m) - zk
  end do
  do m = 1,numtot
    np=ipairs(m)
    itran=ishft(np,-27)
    np = iand(np,mask27)
    j = bckptr(np)
    delx = imagcrds(1,np) + xktran(1,itran)
    dely = imagcrds(2,np) + xktran(2,itran)
    delz = imagcrds(3,np) + xktran(3,itran)
    delr2 = delx*delx + dely*dely+delz*delz

    if ( delr2 < cutoffnb2 ) then
      delr = sqrt(delr2)
      delr2inv = one/delr2

      ! do vdw first
      if ( m <= numvdw ) then
        iaci = ntypes*(iac(i)-1)
        ic = ico(iaci+iac(j))
        delr6inv = delr2inv*delr2inv*delr2inv
        delr12inv = delr6inv*delr6inv
        f6 = cn2(ic)*delr6inv
        f12 = cn1(ic)*delr12inv
        df = (12.d0*f12 - 6.d0*f6)*delr2inv
        ene_vdw = ene_vdw + f12 - f6
        dfx = df*delx
        dfy = df*dely
        dfz = df*delz
        frc(1,i) = frc(1,i) - dfx
        frc(2,i) = frc(2,i) - dfy
        frc(3,i) = frc(3,i) - dfz
        frc(1,j) = frc(1,j) + dfx
        frc(2,j) = frc(2,j) + dfy
        frc(3,j) = frc(3,j) + dfz

        if ( ntb > 0 ) then
        vdw_virial(1,1) = vdw_virial(1,1) - dfx*(delx-(displacement(1,j)-displacement(1,i)))
        vdw_virial(1,2) = vdw_virial(1,2) - dfx*(dely-(displacement(2,j)-displacement(2,i)))
        vdw_virial(1,3) = vdw_virial(1,3) - dfx*(delz-(displacement(3,j)-displacement(3,i)))
        vdw_virial(2,1) = vdw_virial(2,1) - dfy*(delx-(displacement(1,j)-displacement(1,i)))
        vdw_virial(2,2) = vdw_virial(2,2) - dfy*(dely-(displacement(2,j)-displacement(2,i)))
        vdw_virial(2,3) = vdw_virial(2,3) - dfy*(delz-(displacement(3,j)-displacement(3,i)))
        vdw_virial(3,1) = vdw_virial(3,1) - dfz*(delx-(displacement(1,j)-displacement(1,i)))
        vdw_virial(3,2) = vdw_virial(3,2) - dfz*(dely-(displacement(2,j)-displacement(2,i)))
        vdw_virial(3,3) = vdw_virial(3,3) - dfz*(delz-(displacement(3,j)-displacement(3,i)))

        !vdw_virial(1,1) = vdw_virial(1,1) - dfx*delx
        !vdw_virial(1,2) = vdw_virial(1,2) - dfx*dely
        !vdw_virial(1,3) = vdw_virial(1,3) - dfx*delz
        !vdw_virial(2,1) = vdw_virial(2,1) - dfy*delx
        !vdw_virial(2,2) = vdw_virial(2,2) - dfy*dely
        !vdw_virial(2,3) = vdw_virial(2,3) - dfy*delz
        !vdw_virial(3,1) = vdw_virial(3,1) - dfz*delx
        !vdw_virial(3,2) = vdw_virial(3,2) - dfz*dely
        !vdw_virial(3,3) = vdw_virial(3,3) - dfz*delz
        end if ! ( ntb > 0) then

      end if ! ( m <= numvdw ) then
    end if ! ( delr2 < cutoffnb2 ) then

    ! do eel, only grad of phi is computed here.
    if ( delr2 < ee_dsum_cut2 ) then
      !---------------------------------------------------------
      ! McMurchie-Davidson recursion for interaction tensor:
      ! interaction at point charge level given by complementary Boys
      ! B(0) = int_0^1 exp(-pr^2^t^2)dt
      ! complementary Boys is BC(0) = 1/r - B(0)
      ! (d/dr) B(0) = (-2p)*r int_0^1 t^2*exp(-pr^2t^2)dt = (-2p)*r B(1)
      ! and so if R(0,0,0,n) = (-2p)^n B(n) then we have
      ! (1/r)*(d/dr) R(0,0,0,n) = R(0,0,0,n+1)
      ! Now let R(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
      ! Then e.g. R(t+1,u,v,n) = t*R(t-1,u,v,n+1) + x*R(t,u,v,n)
      ! proof:  
      ! R(t+1,u,v,n) = (d/dx)^t+1 (d/dy)^u (d/dz)^v R(0,0,0,n)
      !              = (d/dx)^t (d/dy)^u (d/dz)^v (d/dx) R(0,0,0,n)
      !              = (d/dx)^t (d/dy)^u (d/dz)^v x*(1/r)(d/dr)R(0,0,0,n)
      !              = (d/dx)^t (d/dy)^u (d/dz)^v [x*R(0,0,0,n+1)]
      !              = (d/dx)^t [x*R(0,u,v,n+1)]
      !              = t*R(t-1,u,v,n+1) + x*R(t,u,v,n+1) (Leibniz)
      ! similar recursions hold for R(t,u+1,v,n),R(t,u,v+1,n)
      ! R(t,u+1,v,n) = u*R(t,u-1,v,n+1) + y*R(t,u,v,n+1)
      ! R(t,u,v+1,n) = v*R(t,u,v-1,n+1) + z*R(t,u,v,n+1)
      ! below array is packed---hence use of Ind_tuv
      ! Rn(Ind_tuv) denotes R(t,u,v,n)
      ! Due to its form--we recur downwards in n
      !---------------------------------------------------------
      ! top n is 5 for energy and forces
      ! get boys and R(0,0,0,n), n=0,...,5
      if ( eedmeth /= 1 .and. eedmeth /=4 ) then
        write(6,'(a)') 'This version only supports eedmeth=1 and 4'
        call mexit(6,1)
      end if

      ! for free-space coulombic calculation
      if ( eedmeth == 4 ) then
        ! r_gauss and r_ewald are saved as 2*radi^2 for speed
        betaij  = one/sqrt(r_gauss(i) + r_gauss(j))
        x = betaij*delr
        !if ( x < ee_gauss_cut ) then ! gaussian function coulomb potential
        facij = two*betaij*betaij
        ind = int(eedtbdns*x) + 1
        dx = x - (dble(ind)-one)*del
        switchij = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
               dx*(eed_cub(3,ind)+dx* eed_cub(4,ind)*third)*half)
        d_switchij_dx = eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                     dx*eed_cub(4,ind)*half)

        !call get_ee_func(x,switchij,d_switchij_dx,1)

        B(0) = (one - switchij)*delr*delr2inv ! erf/r
        fact = -d_switchij_dx*betaij
        B(1) = (B(0) - fact)*delr2inv
        fact = facij*fact
        B(2) = (three*B(1) - fact)*delr2inv
        fact = facij*fact
        B(3) = (five*B(2) - fact)*delr2inv
        !end if
      end if ! eedmeth == 4

      ! for PME calculation, only do the cubic spline option
      if ( eedmeth == 1 ) then
        ! compute erfc=one-erf function and its derivative for PME gaussians.
        x = beta*delr
        fac = two*ewaldcof*ewaldcof
        ind = int(eedtbdns*x) + 1
        dx = x - (dble(ind)-one)*del
        switch = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
             dx*(eed_cub(3,ind)+dx* eed_cub(4,ind)*third)*half)
        d_switch_dx = eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                   dx*eed_cub(4,ind)*half)

        !call get_ee_func(x,switch,d_switch_dx,1)

        ! compute erfc=one-erf function and its derivative for pGM gaussians.
        betaij  = one/sqrt(r_gauss(i) + r_gauss(j))
        y = betaij*delr
        !if ( y < ee_gauss_cut ) then 
        facij = two*betaij*betaij
        ind = int(eedtbdns*y) + 1
        dx = y - (dble(ind)-one)*del
        switchij = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
               dx*(eed_cub(3,ind)+dx* eed_cub(4,ind)*third)*half)
        d_switchij_dx = eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                dx*eed_cub(4,ind)*half)

        !call get_ee_func(y,switchij,d_switchij_dx,1)

        !end if

        ! compute Boys series
        ! TD: Got the idea for B_l from Walter Smith's CCP5 article 1982
        ! Ewald for point multipoles
        ! B_l satisfies grad_i B_l(|r_j - r_i|) = (r_j - r_i)B_{l+1}(|r_j-r_i|)
        ! grad_j B_l(|r_j - r_i|) = -grad_i B_l(|r_j - r_i|)
        ! Ruxi: Revised to compute Boys from switch (=erfc())
        !if ( y < ee_gauss_cut ) then ! gaussian function coulomb potential
        B(0) = (one - switch)*delr*delr2inv
        fact = -d_switch_dx*beta
        B(1) = (B(0) - fact)*delr2inv
        fact = fac*fact
        B(2) = (three*B(1) - fact)*delr2inv
        fact = fac*fact
        B(3) = (five*B(2) - fact)*delr2inv

        Bij(0) = (one - switchij)*delr*delr2inv
        fact = -d_switchij_dx*betaij
        Bij(1) = (Bij(0) - fact)*delr2inv
        fact = facij*fact
        Bij(2) = (three*Bij(1) - fact)*delr2inv
        fact = facij*fact
        Bij(3) = (five*Bij(2) - fact)*delr2inv

        ! erf_pgm/r - erf_pme/r
        B(0) = Bij(0) - B(0)
        B(1) = Bij(1) - B(1)
        B(2) = Bij(2) - B(2)
        B(3) = Bij(3) - B(3)
        !end if
      end if ! eedmeth == 1

      ! negate the odd order boys factors
      B(1) = -B(1)
      B(3) = -B(3)

      n = 5
      Rn_2(Ind_000) = B(n-2)
      Rn_3(Ind_000) = B(n-3)
      Rn_3(Ind_100) = delx*Rn_2(Ind_000)
      Rn_3(Ind_010) = dely*Rn_2(Ind_000)
      Rn_3(Ind_001) = delz*Rn_2(Ind_000)
      Rn_4(Ind_000) = B(n-4)
      Rn_4(Ind_100) = delx*Rn_3(Ind_000)
      Rn_4(Ind_010) = dely*Rn_3(Ind_000)
      Rn_4(Ind_001) = delz*Rn_3(Ind_000)
      Rn_4(Ind_200) = Rn_3(Ind_000) + delx*Rn_3(Ind_100)
      Rn_4(Ind_020) = Rn_3(Ind_000) + dely*Rn_3(Ind_010)
      Rn_4(Ind_002) = Rn_3(Ind_000) + delz*Rn_3(Ind_001)
      Rn_5(Ind_000) = B(n-5)
      Rn_5(Ind_100) = delx*Rn_4(Ind_000)
      Rn_5(Ind_010) = dely*Rn_4(Ind_000)
      Rn_5(Ind_001) = delz*Rn_4(Ind_000)
      Rn_5(Ind_200) = Rn_4(Ind_000) + delx*Rn_4(Ind_100)
      Rn_5(Ind_020) = Rn_4(Ind_000) + dely*Rn_4(Ind_010)
      Rn_5(Ind_002) = Rn_4(Ind_000) + delz*Rn_4(Ind_001)
      Rn_5(Ind_110) = delx*Rn_4(Ind_010)
      Rn_5(Ind_101) = delx*Rn_4(Ind_001)
      Rn_5(Ind_011) = dely*Rn_4(Ind_001)
      Rn_5(Ind_300) = 2.d0*Rn_4(Ind_100) + delx*Rn_4(Ind_200)
      Rn_5(Ind_030) = 2.d0*Rn_4(Ind_010) + dely*Rn_4(Ind_020)
      Rn_5(Ind_003) = 2.d0*Rn_4(Ind_001) + delz*Rn_4(Ind_002)
      Rn_5(Ind_210) = dely*Rn_4(Ind_200)
      Rn_5(Ind_201) = delz*Rn_4(Ind_200)
      Rn_5(Ind_120) = delx*Rn_4(Ind_020)
      Rn_5(Ind_021) = delz*Rn_4(Ind_020)
      Rn_5(Ind_102) = delx*Rn_4(Ind_002)
      Rn_5(Ind_012) = dely*Rn_4(Ind_002)
      Rn_5(Ind_111) = delx*dely*delz*Rn_2(Ind_000)

      do jj = 1,4!10
        gmi(jj) = global_multipole(jj,i) ! local perm moments: charges/dipoles
        gmj(jj) = global_multipole(jj,j) ! local perm moments: charges/dipoles
      enddo
      do jj = 1,3
        g_ind(jj) = 0.d0 ! induced force components
        i_di(jj) = ind_dip(jj,i) ! local induced dipoles
        i_dj(jj) = ind_dip(jj,j) ! local induced dipoles
        i_mi(jj) = ind_dip(jj,i) + gmi(jj+1) ! local total dipoles
        i_mj(jj) = ind_dip(jj,j) + gmj(jj+1) ! local total dipoles
      enddo

      ! electrostatic potential at i due to permanent + induced mpoles at j
      ! and derivatives of that with respect to r_i
      tmj(1) = gmj(1)
      tmj(2:4) = i_mj(1:3)

      phi_i(Ind_000)=Rn_5(Ind_000)*tmj(Ind_000)+Rn_5(Ind_100)*tmj(Ind_100)+ &
                     Rn_5(Ind_010)*tmj(Ind_010)+Rn_5(Ind_001)*tmj(Ind_001)
      phi_i(Ind_100)=-(Rn_5(Ind_100)*tmj(Ind_000)+Rn_5(Ind_200)*tmj(Ind_100)+ &
                       Rn_5(Ind_110)*tmj(Ind_010)+Rn_5(Ind_101)*tmj(Ind_001))
      phi_i(Ind_010)=-(Rn_5(Ind_010)*tmj(Ind_000)+Rn_5(Ind_110)*tmj(Ind_100)+ &
                       Rn_5(Ind_020)*tmj(Ind_010)+Rn_5(Ind_011)*tmj(Ind_001))
      phi_i(Ind_001)=-(Rn_5(Ind_001)*tmj(Ind_000)+Rn_5(Ind_101)*tmj(Ind_100)+ &
                       Rn_5(Ind_011)*tmj(Ind_010)+Rn_5(Ind_002)*tmj(Ind_001))
      phi_i(Ind_200)=Rn_5(Ind_200)*tmj(Ind_000)+Rn_5(Ind_300)*tmj(Ind_100)+ &
                     Rn_5(Ind_210)*tmj(Ind_010)+Rn_5(Ind_201)*tmj(Ind_001)
      phi_i(Ind_020)=Rn_5(Ind_020)*tmj(Ind_000)+Rn_5(Ind_120)*tmj(Ind_100)+ &
                     Rn_5(Ind_030)*tmj(Ind_010)+Rn_5(Ind_021)*tmj(Ind_001)
      phi_i(Ind_002)=Rn_5(Ind_002)*tmj(Ind_000)+Rn_5(Ind_102)*tmj(Ind_100)+ &
                     Rn_5(Ind_012)*tmj(Ind_010)+Rn_5(Ind_003)*tmj(Ind_001)
      phi_i(Ind_110)=Rn_5(Ind_110)*tmj(Ind_000)+Rn_5(Ind_210)*tmj(Ind_100)+ &
                     Rn_5(Ind_120)*tmj(Ind_010)+Rn_5(Ind_111)*tmj(Ind_001)
      phi_i(Ind_101)=Rn_5(Ind_101)*tmj(Ind_000)+Rn_5(Ind_201)*tmj(Ind_100)+ &
                     Rn_5(Ind_111)*tmj(Ind_010)+Rn_5(Ind_102)*tmj(Ind_001)
      phi_i(Ind_011)=Rn_5(Ind_011)*tmj(Ind_000)+Rn_5(Ind_111)*tmj(Ind_100)+ &
                     Rn_5(Ind_021)*tmj(Ind_010)+Rn_5(Ind_012)*tmj(Ind_001)

      phi(Ind_000,i)=phi(Ind_000,i) + phi_i(Ind_000)
      phi(Ind_100,i)=phi(Ind_100,i) + phi_i(Ind_100)
      phi(Ind_010,i)=phi(Ind_010,i) + phi_i(Ind_010)
      phi(Ind_001,i)=phi(Ind_001,i) + phi_i(Ind_001)
      phi(Ind_200,i)=phi(Ind_200,i) + phi_i(Ind_200)
      phi(Ind_020,i)=phi(Ind_020,i) + phi_i(Ind_020)
      phi(Ind_002,i)=phi(Ind_002,i) + phi_i(Ind_002)
      phi(Ind_110,i)=phi(Ind_110,i) + phi_i(Ind_110)
      phi(Ind_101,i)=phi(Ind_101,i) + phi_i(Ind_101)
      phi(Ind_011,i)=phi(Ind_011,i) + phi_i(Ind_011)

      tmi(1) = gmi(1)
      tmi(2:4) = i_mi(1:3)
      Rn_5(2:4) = -Rn_5(2:4) ! now these should be derivatives against i
      Rn_5(11:20) = -Rn_5(11:20) ! now these should be derivatives against i

      phi_j(Ind_000)=Rn_5(Ind_000)*tmi(Ind_000)+Rn_5(Ind_100)*tmi(Ind_100)+ &
                     Rn_5(Ind_010)*tmi(Ind_010)+Rn_5(Ind_001)*tmi(Ind_001)!+ &
      phi_j(Ind_100)=-(Rn_5(Ind_100)*tmi(Ind_000)+Rn_5(Ind_200)*tmi(Ind_100)+ &
                       Rn_5(Ind_110)*tmi(Ind_010)+Rn_5(Ind_101)*tmi(Ind_001))!+ &
      phi_j(Ind_010)=-(Rn_5(Ind_010)*tmi(Ind_000)+Rn_5(Ind_110)*tmi(Ind_100)+ &
                       Rn_5(Ind_020)*tmi(Ind_010)+Rn_5(Ind_011)*tmi(Ind_001))!+ &
      phi_j(Ind_001)=-(Rn_5(Ind_001)*tmi(Ind_000)+Rn_5(Ind_101)*tmi(Ind_100)+ &
                       Rn_5(Ind_011)*tmi(Ind_010)+Rn_5(Ind_002)*tmi(Ind_001))!+ &
      phi_j(Ind_200)=Rn_5(Ind_200)*tmi(Ind_000)+Rn_5(Ind_300)*tmi(Ind_100)+ &
                     Rn_5(Ind_210)*tmi(Ind_010)+Rn_5(Ind_201)*tmi(Ind_001)!+ &
      phi_j(Ind_020)=Rn_5(Ind_020)*tmi(Ind_000)+Rn_5(Ind_120)*tmi(Ind_100)+ &
                     Rn_5(Ind_030)*tmi(Ind_010)+Rn_5(Ind_021)*tmi(Ind_001)!+ &
      phi_j(Ind_002)=Rn_5(Ind_002)*tmi(Ind_000)+Rn_5(Ind_102)*tmi(Ind_100)+ &
                     Rn_5(Ind_012)*tmi(Ind_010)+Rn_5(Ind_003)*tmi(Ind_001)!+ &
      phi_j(Ind_110)=Rn_5(Ind_110)*tmi(Ind_000)+Rn_5(Ind_210)*tmi(Ind_100)+ &
                     Rn_5(Ind_120)*tmi(Ind_010)+Rn_5(Ind_111)*tmi(Ind_001)!+ &
      phi_j(Ind_101)=Rn_5(Ind_101)*tmi(Ind_000)+Rn_5(Ind_201)*tmi(Ind_100)+ &
                     Rn_5(Ind_111)*tmi(Ind_010)+Rn_5(Ind_102)*tmi(Ind_001)!+ &
      phi_j(Ind_011)=Rn_5(Ind_011)*tmi(Ind_000)+Rn_5(Ind_111)*tmi(Ind_100)+ &
                     Rn_5(Ind_021)*tmi(Ind_010)+Rn_5(Ind_012)*tmi(Ind_001)!+ &

      phi(Ind_000,j)=phi(Ind_000,j) + phi_j(Ind_000)
      phi(Ind_100,j)=phi(Ind_100,j) + phi_j(Ind_100)
      phi(Ind_010,j)=phi(Ind_010,j) + phi_j(Ind_010)
      phi(Ind_001,j)=phi(Ind_001,j) + phi_j(Ind_001)
      phi(Ind_200,j)=phi(Ind_200,j) + phi_j(Ind_200)
      phi(Ind_020,j)=phi(Ind_020,j) + phi_j(Ind_020)
      phi(Ind_002,j)=phi(Ind_002,j) + phi_j(Ind_002)
      phi(Ind_110,j)=phi(Ind_110,j) + phi_j(Ind_110)
      phi(Ind_101,j)=phi(Ind_101,j) + phi_j(Ind_101)
      phi(Ind_011,j)=phi(Ind_011,j) + phi_j(Ind_011)

      if ( ntb > 0 ) then
      virial(1,1) = virial(1,1) - &
                    half*( gmi(1)*(-phi_i(Ind_100))*(-(delx-(displacement(1,j)-displacement(1,i)))) + &
                          (i_mi(1)*(-phi_i(Ind_200))+i_mi(2)*(-phi_i(Ind_110))+i_mi(3)*(-phi_i(Ind_101)))&
                         *(-(delx-(displacement(1,j)-displacement(1,i)))) ) - &
                    half*( gmj(1)*(-phi_j(Ind_100))*(delx-(displacement(1,j)-displacement(1,i))) + &
                          (i_mj(1)*(-phi_j(Ind_200))+i_mj(2)*(-phi_j(Ind_110))+i_mj(3)*(-phi_j(Ind_101)))&
                         *(delx-(displacement(1,j)-displacement(1,i))) )
    
      virial(1,2) = virial(1,2) - &
                    half*( gmi(1)*(-phi_i(Ind_100))*(-(dely-(displacement(2,j)-displacement(2,i)))) + &
                          (i_mi(1)*(-phi_i(Ind_200))+i_mi(2)*(-phi_i(Ind_110))+i_mi(3)*(-phi_i(Ind_101)))&
                         *(-(dely-(displacement(2,j)-displacement(2,i)))) ) - &
                    half*( gmj(1)*(-phi_j(Ind_100))*(dely-(displacement(2,j)-displacement(2,i))) + &
                          (i_mj(1)*(-phi_j(Ind_200))+i_mj(2)*(-phi_j(Ind_110))+i_mj(3)*(-phi_j(Ind_101)))&
                         *(dely-(displacement(2,j)-displacement(2,i))) )
    
      virial(1,3) = virial(1,3) - &
                    half*( gmi(1)*(-phi_i(Ind_100))*(-(delz-(displacement(3,j)-displacement(3,i)))) + &
                          (i_mi(1)*(-phi_i(Ind_200))+i_mi(2)*(-phi_i(Ind_110))+i_mi(3)*(-phi_i(Ind_101)))&
                         *(-(delz-(displacement(3,j)-displacement(3,i)))) ) - &
                    half*( gmj(1)*(-phi_j(Ind_100))*(delz-(displacement(3,j)-displacement(3,i))) + &
                          (i_mj(1)*(-phi_j(Ind_200))+i_mj(2)*(-phi_j(Ind_110))+i_mj(3)*(-phi_j(Ind_101)))&
                         *(delz-(displacement(3,j)-displacement(3,i))) )

      virial(2,1) = virial(2,1) - &
                    half*( gmi(1)*(-phi_i(Ind_010))*(-(delx-(displacement(1,j)-displacement(1,i)))) + &
                          (i_mi(1)*(-phi_i(Ind_110))+i_mi(2)*(-phi_i(Ind_020))+i_mi(3)*(-phi_i(Ind_011)))&
                         *(-(delx-(displacement(1,j)-displacement(1,i)))) ) - &
                    half*( gmj(1)*(-phi_j(Ind_010))*(delx-(displacement(1,j)-displacement(1,i))) + &
                          (i_mj(1)*(-phi_j(Ind_110))+i_mj(2)*(-phi_j(Ind_020))+i_mj(3)*(-phi_j(Ind_011)))&
                         *(delx-(displacement(1,j)-displacement(1,i))) )

      virial(2,2) = virial(2,2) - &
                    half*( gmi(1)*(-phi_i(Ind_010))*(-(dely-(displacement(2,j)-displacement(2,i)))) + &
                          (i_mi(1)*(-phi_i(Ind_110))+i_mi(2)*(-phi_i(Ind_020))+i_mi(3)*(-phi_i(Ind_011)))&
                         *(-(dely-(displacement(2,j)-displacement(2,i)))) ) - &
                    half*( gmj(1)*(-phi_j(Ind_010))*(dely-(displacement(2,j)-displacement(2,i))) + &
                          (i_mj(1)*(-phi_j(Ind_110))+i_mj(2)*(-phi_j(Ind_020))+i_mj(3)*(-phi_j(Ind_011)))&
                         *(dely-(displacement(2,j)-displacement(2,i))) )

      virial(2,3) = virial(2,3) - &
                    half*( gmi(1)*(-phi_i(Ind_010))*(-(delz-(displacement(3,j)-displacement(3,i)))) + &
                          (i_mi(1)*(-phi_i(Ind_110))+i_mi(2)*(-phi_i(Ind_020))+i_mi(3)*(-phi_i(Ind_011)))&
                         *(-(delz-(displacement(3,j)-displacement(3,i)))) ) - &
                    half*( gmj(1)*(-phi_j(Ind_010))*(delz-(displacement(3,j)-displacement(3,i))) + &
                          (i_mj(1)*(-phi_j(Ind_110))+i_mj(2)*(-phi_j(Ind_020))+i_mj(3)*(-phi_j(Ind_011)))&
                         *(delz-(displacement(3,j)-displacement(3,i))) )

      virial(3,1) = virial(3,1) - &
                    half*( gmi(1)*(-phi_i(Ind_001))*(-(delx-(displacement(1,j)-displacement(1,i)))) + &
                          (i_mi(1)*(-phi_i(Ind_101))+i_mi(2)*(-phi_i(Ind_011))+i_mi(3)*(-phi_i(Ind_002)))&
                         *(-(delx-(displacement(1,j)-displacement(1,i)))) ) - &
                    half*( gmj(1)*(-phi_j(Ind_001))*(delx-(displacement(1,j)-displacement(1,i))) + &
                          (i_mj(1)*(-phi_j(Ind_101))+i_mj(2)*(-phi_j(Ind_011))+i_mj(3)*(-phi_j(Ind_002)))&
                         *(delx-(displacement(1,j)-displacement(1,i))) )

      virial(3,2) = virial(3,2) - &
                    half*( gmi(1)*(-phi_i(Ind_001))*(-(dely-(displacement(2,j)-displacement(2,i)))) + &
                          (i_mi(1)*(-phi_i(Ind_101))+i_mi(2)*(-phi_i(Ind_011))+i_mi(3)*(-phi_i(Ind_002)))&
                         *(-(dely-(displacement(2,j)-displacement(2,i)))) ) - &
                    half*( gmj(1)*(-phi_j(Ind_001))*(dely-(displacement(2,j)-displacement(2,i))) + &
                          (i_mj(1)*(-phi_j(Ind_101))+i_mj(2)*(-phi_j(Ind_011))+i_mj(3)*(-phi_j(Ind_002)))&
                         *(dely-(displacement(2,j)-displacement(2,i))) )

      virial(3,3) = virial(3,3) - &
                    half*( gmi(1)*(-phi_i(Ind_001))*(-(delz-(displacement(3,j)-displacement(3,i)))) + &
                          (i_mi(1)*(-phi_i(Ind_101))+i_mi(2)*(-phi_i(Ind_011))+i_mi(3)*(-phi_i(Ind_002)))&
                         *(-(delz-(displacement(3,j)-displacement(3,i)))) ) - &
                    half*( gmj(1)*(-phi_j(Ind_001))*(delz-(displacement(3,j)-displacement(3,i))) + &
                          (i_mj(1)*(-phi_j(Ind_101))+i_mj(2)*(-phi_j(Ind_011))+i_mj(3)*(-phi_j(Ind_002)))&
                         *(delz-(displacement(3,j)-displacement(3,i))) )

!      virial(1,1) = virial(1,1) - &
!                    half*( gmi(1)*(-phi_i(Ind_100))*(-delx) + &
!                          (i_mi(1)*(-phi_i(Ind_200))+i_mi(2)*(-phi_i(Ind_110))+i_mi(3)*(-phi_i(Ind_101)))&
!                         *(-delx) ) - &
!                    half*( gmj(1)*(-phi_j(Ind_100))*delx + &
!                          (i_mj(1)*(-phi_j(Ind_200))+i_mj(2)*(-phi_j(Ind_110))+i_mj(3)*(-phi_j(Ind_101)))&
!                         *delx )
!    
!      virial(1,2) = virial(1,2) - &
!                    half*( gmi(1)*(-phi_i(Ind_100))*(-dely) + &
!                          (i_mi(1)*(-phi_i(Ind_200))+i_mi(2)*(-phi_i(Ind_110))+i_mi(3)*(-phi_i(Ind_101)))&
!                         *(-dely) ) - &
!                    half*( gmj(1)*(-phi_j(Ind_100))*dely + &
!                          (i_mj(1)*(-phi_j(Ind_200))+i_mj(2)*(-phi_j(Ind_110))+i_mj(3)*(-phi_j(Ind_101)))&
!                         *dely )
!    
!      virial(1,3) = virial(1,3) - &
!                    half*( gmi(1)*(-phi_i(Ind_100))*(-delz) + &
!                          (i_mi(1)*(-phi_i(Ind_200))+i_mi(2)*(-phi_i(Ind_110))+i_mi(3)*(-phi_i(Ind_101)))&
!                         *(-delz) ) - &
!                    half*( gmj(1)*(-phi_j(Ind_100))*delz + &
!                          (i_mj(1)*(-phi_j(Ind_200))+i_mj(2)*(-phi_j(Ind_110))+i_mj(3)*(-phi_j(Ind_101)))&
!                         *delz )
!
!      virial(2,1) = virial(2,1) - &
!                    half*( gmi(1)*(-phi_i(Ind_010))*(-delx) + &
!                          (i_mi(1)*(-phi_i(Ind_110))+i_mi(2)*(-phi_i(Ind_020))+i_mi(3)*(-phi_i(Ind_011)))&
!                         *(-delx) ) - &
!                    half*( gmj(1)*(-phi_j(Ind_010))*delx + &
!                          (i_mj(1)*(-phi_j(Ind_110))+i_mj(2)*(-phi_j(Ind_020))+i_mj(3)*(-phi_j(Ind_011)))&
!                         *delx )
!
!      virial(2,2) = virial(2,2) - &
!                    half*( gmi(1)*(-phi_i(Ind_010))*(-dely) + &
!                          (i_mi(1)*(-phi_i(Ind_110))+i_mi(2)*(-phi_i(Ind_020))+i_mi(3)*(-phi_i(Ind_011)))&
!                         *(-dely) ) - &
!                    half*( gmj(1)*(-phi_j(Ind_010))*dely + &
!                          (i_mj(1)*(-phi_j(Ind_110))+i_mj(2)*(-phi_j(Ind_020))+i_mj(3)*(-phi_j(Ind_011)))&
!                         *dely )
!
!      virial(2,3) = virial(2,3) - &
!                    half*( gmi(1)*(-phi_i(Ind_010))*(-delz) + &
!                          (i_mi(1)*(-phi_i(Ind_110))+i_mi(2)*(-phi_i(Ind_020))+i_mi(3)*(-phi_i(Ind_011)))&
!                         *(-delz) ) - &
!                    half*( gmj(1)*(-phi_j(Ind_010))*delz + &
!                          (i_mj(1)*(-phi_j(Ind_110))+i_mj(2)*(-phi_j(Ind_020))+i_mj(3)*(-phi_j(Ind_011)))&
!                         *delz )
!
!      virial(3,1) = virial(3,1) - &
!                    half*( gmi(1)*(-phi_i(Ind_001))*(-delx) + &
!                          (i_mi(1)*(-phi_i(Ind_101))+i_mi(2)*(-phi_i(Ind_011))+i_mi(3)*(-phi_i(Ind_002)))&
!                         *(-delx) ) - &
!                    half*( gmj(1)*(-phi_j(Ind_001))*delx + &
!                          (i_mj(1)*(-phi_j(Ind_101))+i_mj(2)*(-phi_j(Ind_011))+i_mj(3)*(-phi_j(Ind_002)))&
!                         *delx )
!
!      virial(3,2) = virial(3,2) - &
!                    half*( gmi(1)*(-phi_i(Ind_001))*(-dely) + &
!                          (i_mi(1)*(-phi_i(Ind_101))+i_mi(2)*(-phi_i(Ind_011))+i_mi(3)*(-phi_i(Ind_002)))&
!                         *(-dely) ) - &
!                    half*( gmj(1)*(-phi_j(Ind_001))*dely + &
!                          (i_mj(1)*(-phi_j(Ind_101))+i_mj(2)*(-phi_j(Ind_011))+i_mj(3)*(-phi_j(Ind_002)))&
!                         *dely )
!
!      virial(3,3) = virial(3,3) - &
!                    half*( gmi(1)*(-phi_i(Ind_001))*(-delz) + &
!                          (i_mi(1)*(-phi_i(Ind_101))+i_mi(2)*(-phi_i(Ind_011))+i_mi(3)*(-phi_i(Ind_002)))&
!                         *(-delz) ) - &
!                    half*( gmj(1)*(-phi_j(Ind_001))*delz + &
!                          (i_mj(1)*(-phi_j(Ind_101))+i_mj(2)*(-phi_j(Ind_011))+i_mj(3)*(-phi_j(Ind_002)))&
!                         *delz )
      end if ! ( ntb > 0 ) then

    end if ! ( delr2 < ee_dsum_cut2 )then
  end do ! m = 1,numlist

end subroutine pGM_DIRECT_ene_force_i
!-------------------------------------------------------
end module pol_gauss_direct
