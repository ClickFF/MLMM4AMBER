#include "../include/dprec.fh"
#include "../include/assert.fh"

!------------------------------------------------------------
module pol_gauss_self
  implicit none
  private

# include "pol_gauss_mpole_index.h"

  integer,save :: do_flag

  public pGM_SELF_perm_field,pGM_SELF_dipole_field,pGM_SELF_ene_frc, &
         pGM_SELF_set_user_bit
#ifdef MPI
  public pGM_SELF_bcast
#endif

  contains
!-------------------------------------------------------
#ifdef MPI
subroutine pGM_SELF_bcast
  implicit none
  integer ierr

  include 'mpif.h'
# include "parallel.h"
  call mpi_bcast(do_flag,1,MPI_INTEGER,0,commsander,ierr)
end subroutine pGM_SELF_bcast
#endif
!-------------------------------------------------------
subroutine pGM_SELF_set_user_bit(do_this)
  integer,intent(in) :: do_this
# include "do_flag.h"

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
    write(6,*)'pGM_SELF_set_user_bit: bad value of user do_this'
    call mexit(6,1)
  endif
end subroutine pGM_SELF_set_user_bit
!-------------------------------------------------------
subroutine pGM_SELF_perm_field(numatoms,direct_field)
  use pol_gauss_multipoles, only : global_multipole
  use constants, only : three,four,pi

  integer,intent(in) :: numatoms
  _REAL_,intent(inout) :: direct_field(3,*)

# include "ew_pme_recip.h"
# include "do_flag.h"

  _REAL_ :: factor
  integer n

  if ( iand(do_flag,PROCEED_INDUCE) /= PROCEED_INDUCE )return

  factor = four * ew_coeff**3 / (three*sqrt(pi))
  do n = 1, numatoms
    direct_field(1,n) = direct_field(1,n)-factor*global_multipole(Ind_100,n)
    direct_field(2,n) = direct_field(2,n)-factor*global_multipole(Ind_010,n)
    direct_field(3,n) = direct_field(3,n)-factor*global_multipole(Ind_001,n)
  end do
end subroutine pGM_SELF_perm_field
!-------------------------------------------------------
subroutine pGM_SELF_dipole_field(numatoms,ind_dip,dip_field)
  use constants, only : three,four,pi

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(inout) :: dip_field(3,*)

# include "do_flag.h"
# include "ew_pme_recip.h"

  _REAL_ :: factor
  integer n

  if ( iand(do_flag,PROCEED_INDUCE) /= PROCEED_INDUCE ) return

  factor = four * ew_coeff**3 / (three*sqrt(pi))
  do n = 1, numatoms
    dip_field(1,n) = dip_field(1,n) - factor*ind_dip(1,n)
    dip_field(2,n) = dip_field(2,n) - factor*ind_dip(2,n)
    dip_field(3,n) = dip_field(3,n) - factor*ind_dip(3,n)
  end do
end subroutine pGM_SELF_dipole_field
!-------------------------------------------------------
subroutine pGM_SELF_ene_frc(numatoms,ind_dip,ene_perm,ene_ind,phi)
  use pol_gauss_multipoles, only : global_multipole, &
                                coulomb_const_kcal_per_mole
  use constants, only : zero,half,two,pi

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(inout) :: ene_perm,ene_ind
  _REAL_,intent(inout) :: phi(10,numatoms)

  _REAL_ :: delx,dely,delz,B(0:4),fac,fact,gmi(10),gphi(10),i_di(3), &
            i_mi(3),e_pp,e_ind
  _REAL_ :: Rn(1),Rn_1(4),Rn_2(10),Rn_3(20),Rn_4(35)
  integer i,j,n

# include "do_flag.h"
# include "ew_pme_recip.h"

  if ( iand(do_flag,PROCEED_POSTINDUCE) /= PROCEED_POSTINDUCE )return

  fact = two*ew_coeff / sqrt(PI)
  fac = -two*ew_coeff*ew_coeff
  do j = 0,4
    B(j) = fact/(two*j+1)
    fact = fac*fact
  end do
  delx = zero
  dely = zero
  delz = zero

  n = 4
  Rn(Ind_000) = B(n)
  Rn_1(Ind_000) = B(n-1)
  Rn_1(Ind_100) = delx*Rn(Ind_000)
  Rn_1(Ind_010) = dely*Rn(Ind_000)
  Rn_1(Ind_001) = delz*Rn(Ind_000)
  Rn_2(Ind_000) = B(n-2)
  Rn_2(Ind_100) = delx*Rn_1(Ind_000)
  Rn_2(Ind_010) = dely*Rn_1(Ind_000)
  Rn_2(Ind_001) = delz*Rn_1(Ind_000)
  Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
  Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
  Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
  Rn_2(Ind_110) = delx*Rn_1(Ind_010)
  Rn_2(Ind_101) = delx*Rn_1(Ind_001)
  Rn_2(Ind_011) = dely*Rn_1(Ind_001)
  Rn_3(Ind_000) = B(n-3) 
  Rn_3(Ind_100) = delx*Rn_2(Ind_000)
  Rn_3(Ind_010) = dely*Rn_2(Ind_000)
  Rn_3(Ind_001) = delz*Rn_2(Ind_000)
  Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
  Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
  Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
  Rn_3(Ind_110) = delx*Rn_2(Ind_010)
  Rn_3(Ind_101) = delx*Rn_2(Ind_001)
  Rn_3(Ind_011) = dely*Rn_2(Ind_001)
  Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
  Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
  Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
  Rn_3(Ind_210) = dely*Rn_2(Ind_200)
  Rn_3(Ind_201) = delz*Rn_2(Ind_200)
  Rn_3(Ind_120) = delx*Rn_2(Ind_020)
  Rn_3(Ind_021) = delz*Rn_2(Ind_020)
  Rn_3(Ind_102) = delx*Rn_2(Ind_002)
  Rn_3(Ind_012) = dely*Rn_2(Ind_002)
  Rn_3(Ind_111) = delx*Rn_2(Ind_011)
  Rn_4(Ind_000) = B(n-4) 
  Rn_4(Ind_100) = delx*Rn_3(Ind_000)
  Rn_4(Ind_010) = dely*Rn_3(Ind_000)
  Rn_4(Ind_001) = delz*Rn_3(Ind_000)
  Rn_4(Ind_200) = Rn_3(Ind_000) + delx*Rn_3(Ind_100)
  Rn_4(Ind_020) = Rn_3(Ind_000) + dely*Rn_3(Ind_010)
  Rn_4(Ind_002) = Rn_3(Ind_000) + delz*Rn_3(Ind_001)
  Rn_4(Ind_110) = delx*Rn_3(Ind_010)
  Rn_4(Ind_101) = delx*Rn_3(Ind_001)
  Rn_4(Ind_011) = dely*Rn_3(Ind_001)
  Rn_4(Ind_300) = 2.d0*Rn_3(Ind_100) + delx*Rn_3(Ind_200)
  Rn_4(Ind_030) = 2.d0*Rn_3(Ind_010) + dely*Rn_3(Ind_020)
  Rn_4(Ind_003) = 2.d0*Rn_3(Ind_001) + delz*Rn_3(Ind_002)
  Rn_4(Ind_210) = dely*Rn_3(Ind_200)
  Rn_4(Ind_201) = delz*Rn_3(Ind_200)
  Rn_4(Ind_120) = delx*Rn_3(Ind_020)
  Rn_4(Ind_021) = delz*Rn_3(Ind_020)
  Rn_4(Ind_102) = delx*Rn_3(Ind_002)
  Rn_4(Ind_012) = dely*Rn_3(Ind_002)
  Rn_4(Ind_111) = delx*Rn_3(Ind_011)
  Rn_4(Ind_400) = 3.d0*Rn_3(Ind_200) + delx*Rn_3(Ind_300)
  Rn_4(Ind_040) = 3.d0*Rn_3(Ind_020) + dely*Rn_3(Ind_030)
  Rn_4(Ind_004) = 3.d0*Rn_3(Ind_002) + delz*Rn_3(Ind_003)
  Rn_4(Ind_310) = dely*Rn_3(Ind_300)
  Rn_4(Ind_301) = delz*Rn_3(Ind_300)
  Rn_4(Ind_130) = delx*Rn_3(Ind_030)
  Rn_4(Ind_031) = delz*Rn_3(Ind_030)
  Rn_4(Ind_103) = delx*Rn_3(Ind_003)
  Rn_4(Ind_013) = dely*Rn_3(Ind_003)
  Rn_4(Ind_220) = Rn_3(Ind_020) + delx*Rn_3(Ind_120)
  Rn_4(Ind_202) = Rn_3(Ind_002) + delx*Rn_3(Ind_102)
  Rn_4(Ind_022) = Rn_3(Ind_002) + dely*Rn_3(Ind_012)
  Rn_4(Ind_211) = dely*Rn_3(Ind_201)
  Rn_4(Ind_121) = delx*Rn_3(Ind_021)
  Rn_4(Ind_112) = delx*Rn_3(Ind_012)
  do i = 1,numatoms
    do j = 1,10
      gmi(j) = global_multipole(j,i)
    end do
    do j = 1,3
      i_di(j) = ind_dip(j,i)
    end do
    gmi(2:4) = gmi(2:4) + i_di(1:3)

    ! self-field due to permanent mpoles at i and derivs wrt r_j-r_i (at 0)
    gphi(Ind_000)=  Rn_4(Ind_000)*gmi(Ind_000)+Rn_4(Ind_100)*gmi(Ind_100)+ &
                    Rn_4(Ind_010)*gmi(Ind_010)+Rn_4(Ind_001)*gmi(Ind_001)+ &
                    Rn_4(Ind_200)*gmi(Ind_200)+Rn_4(Ind_020)*gmi(Ind_020)+ &
                    Rn_4(Ind_002)*gmi(Ind_002)+Rn_4(Ind_110)*gmi(Ind_110)+ &
                    Rn_4(Ind_101)*gmi(Ind_101)+Rn_4(Ind_011)*gmi(Ind_011)
    gphi(Ind_100)=-(Rn_4(Ind_100)*gmi(Ind_000)+Rn_4(Ind_200)*gmi(Ind_100)+ &
                    Rn_4(Ind_110)*gmi(Ind_010)+Rn_4(Ind_101)*gmi(Ind_001)+ &
                    Rn_4(Ind_300)*gmi(Ind_200)+Rn_4(Ind_120)*gmi(Ind_020)+ &
                    Rn_4(Ind_102)*gmi(Ind_002)+Rn_4(Ind_210)*gmi(Ind_110)+ &
                    Rn_4(Ind_201)*gmi(Ind_101)+Rn_4(Ind_111)*gmi(Ind_011))
    gphi(Ind_010)=-(Rn_4(Ind_010)*gmi(Ind_000)+Rn_4(Ind_110)*gmi(Ind_100)+ &
                    Rn_4(Ind_020)*gmi(Ind_010)+Rn_4(Ind_011)*gmi(Ind_001)+ &
                    Rn_4(Ind_210)*gmi(Ind_200)+Rn_4(Ind_030)*gmi(Ind_020)+ &
                    Rn_4(Ind_012)*gmi(Ind_002)+Rn_4(Ind_120)*gmi(Ind_110)+ &
                    Rn_4(Ind_111)*gmi(Ind_101)+Rn_4(Ind_021)*gmi(Ind_011))
    gphi(Ind_001)=-(Rn_4(Ind_001)*gmi(Ind_000)+Rn_4(Ind_101)*gmi(Ind_100)+ &
                    Rn_4(Ind_011)*gmi(Ind_010)+Rn_4(Ind_002)*gmi(Ind_001)+ &
                    Rn_4(Ind_201)*gmi(Ind_200)+Rn_4(Ind_021)*gmi(Ind_020)+ &
                    Rn_4(Ind_003)*gmi(Ind_002)+Rn_4(Ind_111)*gmi(Ind_110)+ &
                    Rn_4(Ind_102)*gmi(Ind_101)+Rn_4(Ind_012)*gmi(Ind_011))
    gphi(Ind_200)=  Rn_4(Ind_200)*gmi(Ind_000)+Rn_4(Ind_300)*gmi(Ind_100)+ &
                    Rn_4(Ind_210)*gmi(Ind_010)+Rn_4(Ind_201)*gmi(Ind_001)+ &
                    Rn_4(Ind_400)*gmi(Ind_200)+Rn_4(Ind_220)*gmi(Ind_020)+ &
                    Rn_4(Ind_202)*gmi(Ind_002)+Rn_4(Ind_310)*gmi(Ind_110)+ &
                    Rn_4(Ind_301)*gmi(Ind_101)+Rn_4(Ind_211)*gmi(Ind_011)
    gphi(Ind_020)=  Rn_4(Ind_020)*gmi(Ind_000)+Rn_4(Ind_120)*gmi(Ind_100)+ &
                    Rn_4(Ind_030)*gmi(Ind_010)+Rn_4(Ind_021)*gmi(Ind_001)+ &
                    Rn_4(Ind_220)*gmi(Ind_200)+Rn_4(Ind_040)*gmi(Ind_020)+ &
                    Rn_4(Ind_022)*gmi(Ind_002)+Rn_4(Ind_130)*gmi(Ind_110)+ &
                    Rn_4(Ind_121)*gmi(Ind_101)+Rn_4(Ind_031)*gmi(Ind_011)
    gphi(Ind_002)=  Rn_4(Ind_002)*gmi(Ind_000)+Rn_4(Ind_102)*gmi(Ind_100)+ &
                    Rn_4(Ind_012)*gmi(Ind_010)+Rn_4(Ind_003)*gmi(Ind_001)+ &
                    Rn_4(Ind_202)*gmi(Ind_200)+Rn_4(Ind_022)*gmi(Ind_020)+ &
                    Rn_4(Ind_004)*gmi(Ind_002)+Rn_4(Ind_112)*gmi(Ind_110)+ &
                    Rn_4(Ind_103)*gmi(Ind_101)+Rn_4(Ind_013)*gmi(Ind_011)
    gphi(Ind_110)=  Rn_4(Ind_110)*gmi(Ind_000)+Rn_4(Ind_210)*gmi(Ind_100)+ &
                    Rn_4(Ind_120)*gmi(Ind_010)+Rn_4(Ind_111)*gmi(Ind_001)+ &
                    Rn_4(Ind_310)*gmi(Ind_200)+Rn_4(Ind_130)*gmi(Ind_020)+ &
                    Rn_4(Ind_112)*gmi(Ind_002)+Rn_4(Ind_220)*gmi(Ind_110)+ &
                    Rn_4(Ind_211)*gmi(Ind_101)+Rn_4(Ind_121)*gmi(Ind_011)
    gphi(Ind_101)=  Rn_4(Ind_101)*gmi(Ind_000)+Rn_4(Ind_201)*gmi(Ind_100)+ &
                    Rn_4(Ind_111)*gmi(Ind_010)+Rn_4(Ind_102)*gmi(Ind_001)+ &
                    Rn_4(Ind_301)*gmi(Ind_200)+Rn_4(Ind_121)*gmi(Ind_020)+ &
                    Rn_4(Ind_103)*gmi(Ind_002)+Rn_4(Ind_211)*gmi(Ind_110)+ &
                    Rn_4(Ind_202)*gmi(Ind_101)+Rn_4(Ind_112)*gmi(Ind_011)
    gphi(Ind_011)=  Rn_4(Ind_011)*gmi(Ind_000)+Rn_4(Ind_111)*gmi(Ind_100)+ &
                    Rn_4(Ind_021)*gmi(Ind_010)+Rn_4(Ind_012)*gmi(Ind_001)+ &
                    Rn_4(Ind_211)*gmi(Ind_200)+Rn_4(Ind_031)*gmi(Ind_020)+ &
                    Rn_4(Ind_013)*gmi(Ind_002)+Rn_4(Ind_121)*gmi(Ind_110)+ &
                    Rn_4(Ind_112)*gmi(Ind_101)+Rn_4(Ind_022)*gmi(Ind_011)
    phi(1:10,i) = phi(1:10,i) - gphi(1:10)
  end do
end subroutine pGM_SELF_ene_frc
!-------------------------------------------------------
end module pol_gauss_self
