! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

module pol_gauss

implicit none

! for Damped Gaussian Model, use UNIT 19 to read in exponents from parm
! and 20 to debug
integer     POLG_UNIT
parameter ( POLG_UNIT = 19 )
integer     POLG_DEBUG_UNIT
parameter ( POLG_DEBUG_UNIT = 20 )

!File number to read in md input
integer :: MDIN_UNIT
   parameter ( MDIN_UNIT = 5 )

!MD options
integer :: ipolg
   _REAL_  :: chg_chg_14, chg_pol_14
   _REAL_  :: pol_pol_12, pol_pol_13, pol_pol_14 
   _REAL_  :: erfc_tol

!Debug Variable
   logical :: POLG_DEBUG
   parameter ( POLG_DEBUG = .TRUE. )


!Pt Charge and Gaussian Polarizability Exponent
   _REAL_, dimension(:), allocatable  :: Q_tot, b_polar

!Cuttoff radii for erfc(x) r = 1/b_exp*erfc_inverse(x)
   _REAL_, dimension(:), allocatable  :: r_polar
   
!Interaction Lists
   integer :: n_list_12, n_list_13, n_list_14
   integer, dimension(:, :), allocatable  :: list_12, list_13, list_14

!Constant
   _REAL_  :: PI, SQRT_PI, EEDT, EEDT_INV, third, half, erfc_inv_tol

   PUBLIC :: ipolg, chg_chg_14, chg_pol_14, &
             pol_pol_12, pol_pol_13, pol_pol_14, &
             erfc_tol 

   PRIVATE ::POLG_DEBUG,b_polar, &
n_list_12, n_list_13, n_list_14, &
list_12, list_13, list_14, &
PI, SQRT_PI, EEDT, EEDT_INV, third, half, &
            erfc_inv_tol

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!I/O and INITIALIZATION FUNCTIONS!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_pol_gauss()
      ! read in md options 

         use file_io_dat, only : mdin_pol_gauss
         implicit none
!#        include "files.h"
         
         namelist /pol_gauss/ &
chg_chg_14, chg_pol_14, &
            pol_pol_12, pol_pol_13, pol_pol_14, &
            erfc_tol 
        
	!set up default values for scales
         chg_chg_14 = 0.8d0
         chg_pol_14 = 0.8d0

         pol_pol_12 = 1.0d0
         pol_pol_13 = 1.0d0
         pol_pol_14 = 1.0d0

         erfc_tol   = 1.0E-6
         if ( mdin_pol_gauss ) read(MDIN_UNIT,nml=pol_gauss)

      end subroutine read_pol_gauss

      subroutine setup_pol_gauss()
      ! read in exponent, Z effective nucleus, 1-2, 1-3, 1-4 interaction lists
         use file_io_dat, only : parm
         implicit none
!#  include "files.h"
#  include "ew_erfc_spline.h"
         integer :: i
         integer :: iostatv, alloc_status
         character(1000) :: line, ref_line
         character(1000) :: parm_filename

         !local variables from parm file
         integer :: natom_loc
         integer :: i_temp(1:10)
         _REAL_ :: erfc_x

         !write to debugging file
         if (POLG_DEBUG) &
            call amopen(POLG_DEBUG_UNIT,'Debug_Gauss.txt','R','F','W')
  
 
         !open up parm top file
         parm_filename = parm
         call amopen(POLG_UNIT,parm,'O','F','R')
         

      !read in local variables from parm file
         ref_line(1:14) = "%FLAG POINTERS"
         call Go_To_Line (POLG_UNIT, ref_line, 14, parm_filename)
         read(POLG_UNIT, '(A1000)', iostat=iostatv) line
         read(POLG_UNIT, '(10I8)', iostat=iostatv) i_temp(1:10)
         natom_loc = i_temp(1)
         read(POLG_UNIT, '(10I8)', iostat=iostatv) i_temp(1:10)
         read(POLG_UNIT, '(10I8)', iostat=iostatv) i_temp(1:10)
         read(POLG_UNIT, '(10I8)', iostat=iostatv) i_temp(1:4)
         n_list_12 = i_temp(2)
         n_list_13 = i_temp(3)
         n_list_14 = i_temp(4)

         if (POLG_DEBUG) then
            write(POLG_DEBUG_UNIT,'(1(A20, I4))') 'natom_loc = ', natom_loc
            write(POLG_DEBUG_UNIT,'(1(A20, I4))') 'n_list_12 = ', n_list_12
            write(POLG_DEBUG_UNIT,'(1(A20, I4))') 'n_list_13 = ', n_list_13
            write(POLG_DEBUG_UNIT,'(1(A20, I4))') 'n_list_14 = ', n_list_14
            write(POLG_DEBUG_UNIT, '(///)')
         end if

      !read in Gaussian Polarizability Exponents
          if (ipolg .eq. 1) then
            allocate(b_polar(1:natom_loc), stat=alloc_status)
            if (alloc_status /= 0) then
               print *, 'Error in pol_gauss.f allocating b_polar'
               STOP
            end if
            ref_line(1:29) = "%FLAG POLARIZABILITY EXPONENT"
            call Go_To_Line (POLG_UNIT, ref_line, 29, parm_filename)
            read(POLG_UNIT, '(A1000)', iostat=iostatv) line
            call read_real_data (POLG_UNIT, 5, natom_loc, b_polar)
            if (POLG_DEBUG) then
               write(POLG_DEBUG_UNIT,*) 'Finding b_polar'
               call write_real_data (POLG_DEBUG_UNIT, natom_loc, b_polar)
               write(POLG_DEBUG_UNIT, '(///)')
            end if
         end if

      !read in total charge, Q_tot
         allocate(Q_tot(1:natom_loc), stat=alloc_status)
         if (alloc_status /= 0) then
            print *, 'Error in pol_gauss.f allocating Q_tot'
            STOP
         end if
         ref_line(1:12) = "%FLAG CHARGE"
         call Go_To_Line (POLG_UNIT, ref_line, 12, parm_filename)
         read(POLG_UNIT, '(A1000)', iostat=iostatv) line
         call read_real_data (POLG_UNIT, 5, natom_loc, Q_tot)
         if (POLG_DEBUG) then
            write(POLG_DEBUG_UNIT,*) 'Finding Q_tot'
            call write_real_data (POLG_DEBUG_UNIT, natom_loc, Q_tot)
            write(POLG_DEBUG_UNIT, '(///)')
         end if


      !read in interaction lists
         allocate(list_12(1:2, 1:n_list_12), stat=alloc_status)
         if (alloc_status /= 0) then
            print *, 'Error in pol_gauss.f allocating list_12'
            STOP
         end if
         allocate(list_13(1:2, 1:n_list_13), stat=alloc_status)
         if (alloc_status /= 0) then
            print *, 'Error in pol_gauss.f allocating list_13'
            STOP
         end if
         allocate(list_14(1:2, 1:n_list_14), stat=alloc_status)
         if (alloc_status /= 0) then
            print *, 'Error in pol_gauss.f allocating list_14'
            STOP
         end if

      !1-2 list
         ref_line(1:14) = "%FLAG 1-2 LIST"
         call Go_To_Line (POLG_UNIT, ref_line, 14, parm_filename)
         read(POLG_UNIT, '(A1000)', iostat=iostatv) line
         call read_interact_list (POLG_UNIT, 10, n_list_12, list_12)
         if (POLG_DEBUG) then
            write(POLG_DEBUG_UNIT,*) 'Finding list_12'
            call write_interact_list (POLG_DEBUG_UNIT, n_list_12, list_12)
            write(POLG_DEBUG_UNIT, '(///)')
         end if
      !1-3 list
         ref_line(1:14) = "%FLAG 1-3 LIST"
         call Go_To_Line (POLG_UNIT, ref_line, 14, parm_filename)
         read(POLG_UNIT, '(A1000)', iostat=iostatv) line
         call read_interact_list (POLG_UNIT, 10, n_list_13, list_13)
         if (POLG_DEBUG) then
            write(POLG_DEBUG_UNIT,*) 'Finding list_13'
            call write_interact_list (POLG_DEBUG_UNIT, n_list_13, list_13)
            write(POLG_DEBUG_UNIT, '(///)')
         end if
      !1-4 list
         ref_line(1:14) = "%FLAG 1-4 LIST"
         call Go_To_Line (POLG_UNIT, ref_line, 14, parm_filename)
         read(POLG_UNIT, '(A1000)', iostat=iostatv) line
         call read_interact_list (POLG_UNIT, 10, n_list_14, list_14)
         if (POLG_DEBUG) then
            write(POLG_DEBUG_UNIT,*) 'Finding list_14'
            call write_interact_list (POLG_DEBUG_UNIT, n_list_14, list_14)
            write(POLG_DEBUG_UNIT, '(///)')
         end if

      !evaluate constants
         PI             = 3.1415926535897932384626433832795d0
         SQRT_PI        = sqrt(PI)
         EEDT           = eedtbdns
         EEDT_INV       = 1.0d0/eedtbdns
         third          = 1.0d0/3.0d0
         half           = 1.0d0/2.0d0
         erfc_inv_tol   = erfc_inv(erfc_tol)

         if (POLG_DEBUG) then
            write(POLG_DEBUG_UNIT,*) 'PI           = ', PI
            write(POLG_DEBUG_UNIT,*) 'SQRT_PI      = ', SQRT_PI
            write(POLG_DEBUG_UNIT,*) 'EEDT         = ', EEDT
            write(POLG_DEBUG_UNIT,*) 'EEDT_INV     = ', EEDT_INV
            write(POLG_DEBUG_UNIT,*) 'third        = ', third
            write(POLG_DEBUG_UNIT,*) 'half         = ', half
            write(POLG_DEBUG_UNIT,*) 'erfc_inv_tol = ', erfc_inv_tol
            call erfcfun(erfc_inv_tol, erfc_x)
            write(POLG_DEBUG_UNIT,*) 'tol          = ', erfc_x
         end if


      !cutoff radii for erfc(x)
      !r = 1/b_exp*erfc_inverse(x)
         if (ipolg .eq. 1) then
            allocate(r_polar(1:natom_loc), stat=alloc_status)
            if (alloc_status /= 0) then
               print *, 'Error in pol_gauss.f allocating r_polar'
               STOP
            end if

            do i=1, natom_loc
               r_polar(i) = 1.0d0/b_polar(i)*erfc_inv_tol
            end do
            if (POLG_DEBUG) then
               write(POLG_DEBUG_UNIT,*) 'r_polar'
               call write_real_data (POLG_DEBUG_UNIT, natom_loc, r_polar)
               write(POLG_DEBUG_UNIT, '(///)')
            end if
         end if 
      end subroutine setup_pol_gauss



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  1-2, 1-3, 1-4 INTRAMOLECULAR        !
!  ELECTROSTATIC +    GAUSSIAN         !
!  POLARIZATION CONTRIBUTION           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine Get_Intra_Gauss_Dip (crd,frc,vir,eelg, epolg,dipole,field)
     !charge and nuclei contribution to 1-2, 1-3, 1-4 interactions
         _REAL_ crd(3,*),frc(3,*),vir(3,3)
         _REAL_ dipole(3,*),field(3,*)
         _REAL_, intent(INOUT) :: eelg, epolg

         !1-2 Interactions
         call Get_1_X_Inter_Dip (0.0d0, 0.0d0, pol_pol_12, &
            n_list_12, list_12, crd, frc, vir, eelg, epolg, dipole, field)

         !1-3 Interactions
         call Get_1_X_Inter_Dip (0.0d0, 0.0d0, pol_pol_13, &
            n_list_13, list_13, crd, frc, vir, eelg, epolg, dipole, field)

         !1-4 Interactions
         call Get_1_X_Inter_Dip (chg_chg_14, chg_pol_14, pol_pol_14, &
            n_list_14, list_14, crd, frc, vir, eelg, epolg, dipole, field)

      end subroutine Get_Intra_Gauss_Dip
      
      subroutine Get_1_X_Inter_Dip (chg_chg_sc, chg_pol_sc, pol_pol_sc, &
            n_inter, inter_list, crd, frc, vir, eelg, epolg, dipole, field)
         
         _REAL_ crd(3,*),frc(3,*), dipole(3,*),field(3,*),vir(3,3)
         _REAL_, intent(INOUT)   :: eelg, epolg
         _REAL_, intent(IN)      :: chg_chg_sc, chg_pol_sc, pol_pol_sc
         integer, dimension(:, :), intent(IN) :: inter_list
         integer, intent(IN)     :: n_inter
         
         integer :: i, j, k, l
         _REAL_  :: delx, dely, delz, delr2, delr, dfx,dfy,dfz
         _REAL_  :: qq, c0, c1, c2, c3, b_exp
         _REAL_  :: dotir, dotjr, dotij, term, termi, termj, termr
         _REAL_  :: f_loc_r, f_loc_i, f_loc_j, efield_i_r, efield_j_r
         _REAL_  :: efield_i_j, efield_j_i

         do k = 1,3
            do l = 1,3
               vir(k,l) = 0.d0
            end do
         end do

         !df*delx = x component of frc on j
         do k=1, n_inter
            i = inter_list(1,k)
            j = inter_list(2,k)

            delx  = crd(1, j) - crd(1, i)
            dely  = crd(2, j) - crd(2, i)
            delz  = crd(3, j) - crd(3, i)
            delr2 = delx*delx + dely*dely + delz*delz
            delr  = sqrt(delr2)
            dotir = dipole(1,i)*delx + dipole(2,i)*dely + dipole(3,i)*delz
            dotjr = dipole(1,j)*delx + dipole(2,j)*dely + dipole(3,j)*delz
            dotij = dipole(1,i)*dipole(1,j) + dipole(2,i)*dipole(2,j) &
                  + dipole(3,i)*dipole(3,j)

         !force on particle j = f_loc_r*r + f_loc_i*dipole(i) + f_loc_j*dipole(j)
            f_loc_r = 0.0d0  
            f_loc_i = 0.0d0   
            f_loc_j = 0.0d0            
            efield_i_r = 0.0d0      
            efield_j_r = 0.0d0

      ! permanent charge - permanent charge
            if (chg_chg_sc > 0.0d0) then
            !Q_tot(j) - Q_tot(i)
               c0      = 1.0d0/delr
               c1      = 1.0d0/(delr2*delr)
               qq      = Q_tot(j)*Q_tot(i)
               eelg    = eelg    + qq*c0*chg_chg_sc
               f_loc_r = f_loc_r + qq*c1*chg_chg_sc
            end if

      ! permanent charge - induced interactions
            if (chg_pol_sc > 0.0d0) then
            !Point Charges and Gaussian Dipoles

            !Q_tot(j) - dip(i)
               b_exp = b_polar(i)
               call Calc_c0_c1_c2 (delr, b_exp, c0, c1, c2)
               term  = Q_tot(j)*dotir
               termi = -Q_tot(j)*c1
               termr = term*c2
               epolg      = epolg + term*c1*chg_pol_sc
               f_loc_r    = f_loc_r + termr*chg_pol_sc
               f_loc_i    = f_loc_i + termi*chg_pol_sc
               efield_i_r = efield_i_r + termi*chg_pol_sc
   
            !dip(j) - Q_tot(i)
               b_exp = b_polar(j)
               call Calc_c0_c1_c2 (delr, b_exp, c0, c1, c2)
               term  = -Q_tot(i)*dotjr
               termj = Q_tot(i)*c1
               termr = term*c2
               epolg      = epolg + term*c1*chg_pol_sc
               f_loc_r    = f_loc_r + termr*chg_pol_sc
               f_loc_j    = f_loc_j + termj*chg_pol_sc
               efield_j_r = efield_j_r + termj*chg_pol_sc
            end if

      ! induced interactions - induced interactions
         ! Gaussian correction only possible
            b_exp = b_polar(j)*b_polar(i)/ &
                     SQRT(b_polar(j)*b_polar(j) + b_polar(i)*b_polar(i) )
            call Calc_c0_c1_c2_c3 (delr, b_exp, c0, c1, c2, c3)

            term  = -dotir*dotjr
            termi = dotjr*c2
            termj = dotir*c2
            termr = dotij*c2 + term*c3

            epolg      = epolg + (dotij*c1 + term*c2)*pol_pol_sc
            f_loc_r    = f_loc_r + (dotij*c2 + term*c3)*pol_pol_sc
            f_loc_i    = f_loc_i + termi*pol_pol_sc
            f_loc_j    = f_loc_j + termj*pol_pol_sc
            efield_i_r = efield_i_r + termi*pol_pol_sc
            efield_j_r = efield_j_r + termj*pol_pol_sc
            efield_i_j = -c1*pol_pol_sc
            efield_j_i = -c1*pol_pol_sc

         !collect frc terms on j
            dfx       = f_loc_r*delx + f_loc_i*dipole(1,i) + f_loc_j*dipole(1,j)
            dfy       = f_loc_r*dely + f_loc_i*dipole(2,i) + f_loc_j*dipole(2,j)
            dfz       = f_loc_r*delz + f_loc_i*dipole(3,i) + f_loc_j*dipole(3,j)
            frc(1, j) = frc(1, j) + dfx
            frc(2, j) = frc(2, j) + dfy
            frc(3, j) = frc(3, j) + dfz
            frc(1, i) = frc(1, i) - dfx
            frc(2, i) = frc(2, i) - dfy
            frc(3, i) = frc(3, i) - dfz

         !collect field terms j
            field(1,j) = field(1,j) + efield_j_r*delx + efield_j_i*dipole(1,i)
            field(2,j) = field(2,j) + efield_j_r*dely + efield_j_i*dipole(2,i)
            field(3,j) = field(3,j) + efield_j_r*delz + efield_j_i*dipole(3,i)
     
         !collect field terms i
            field(1,i) = field(1,i) + efield_i_r*delx + efield_i_j*dipole(1,j)
            field(2,i) = field(2,i) + efield_i_r*dely + efield_i_j*dipole(2,j)
            field(3,i) = field(3,i) + efield_i_r*delz + efield_i_j*dipole(3,j)

            vir(1,1) = vir(1,1) - dfx*delx
            vir(1,2) = vir(1,2) - dfx*dely
            vir(1,3) = vir(1,3) - dfx*delz
            vir(2,1) = vir(2,1) - dfy*delx
            vir(2,2) = vir(2,2) - dfy*dely
            vir(2,3) = vir(2,3) - dfy*delz
            vir(3,1) = vir(3,1) - dfz*delx
            vir(3,2) = vir(3,2) - dfz*dely
            vir(3,3) = vir(3,3) - dfz*delz
          
         end do
      end subroutine Get_1_X_Inter_Dip



      subroutine Calc_c0_c1 (delr, b_exp, c0, c1)
      !c0 = erf(b_exp*r)/r and cn+1 = -1/r dcn/dr
      !(similiar to bn except erf(x) instead of erfc(x) )
      !cn = 1/r2{ (2n-1)*cn-1 - 
      !     (2*b_exp*b_exp)^n/b_exp/sqrt(Pi)*exp(-b_exp*b_exp*r*r) }
         implicit none         

         _REAL_, intent(IN)  :: delr, b_exp
         _REAL_, intent(OUT) :: c0, c1
         _REAL_ :: x, switch, d_switch_dx, fact, delr2inv
         _REAL_ :: erf_x, erfc_x

         x        = delr*b_exp
         delr2inv = 1.0d0/(delr*delr)

         !later change to spline fit for switch and d_switch_dx
         call erfcfun(x, erfc_x)
         erf_x       = 1.0d0 - erfc_x 
         switch      = erf_x
         d_switch_dx = 2.0d0/SQRT_PI*exp(-x*x)

         c0   = switch*delr*delr2inv
         fact = d_switch_dx*b_exp
         c1   = (c0 - fact)*delr2inv
   

      end subroutine Calc_c0_c1

      subroutine Calc_c0_c1_c2 (delr, b_exp, c0, c1, c2)
      !c0 = erf(b_exp*r)/r and cn+1 = -1/r dcn/dr
      !(similiar to bn except erf(x) instead of erfc(x) )
      !cn = 1/r2{ (2n-1)*cn-1 - 
      !     (2*b_exp*b_exp)^n/b_exp/sqrt(Pi)*exp(-b_exp*b_exp*r*r) }
         implicit none         

         _REAL_, intent(IN)  :: delr, b_exp
         _REAL_, intent(OUT) :: c0, c1, c2
         _REAL_ :: x, switch, d_switch_dx, fac, fact, delr2inv
         _REAL_ :: erf_x, erfc_x

         x        = delr*b_exp
         delr2inv = 1.0d0/(delr*delr)
         fac      = 2.0d0*b_exp*b_exp

         !later change to spline fit for switch and d_switch_dx
         call erfcfun(x, erfc_x)
         erf_x       = 1.0d0 - erfc_x
         switch      = erf_x
         d_switch_dx = 2.0d0/SQRT_PI*exp(-x*x)

         c0   = switch*delr*delr2inv
         fact = d_switch_dx*b_exp
         c1   = (c0 - fact)*delr2inv
         fact = fac*fact
         c2   = (3.0d0*c1 - fact)*delr2inv

      end subroutine Calc_c0_c1_c2

      subroutine Calc_c0_c1_c2_c3 (delr, b_exp, c0, c1, c2, c3)
      !c0 = erf(b_exp*r)/r and cn+1 = -1/r dcn/dr
      !(similiar to bn except erf(x) instead of erfc(x) )
      !cn = 1/r2{ (2n-1)*cn-1 - 
      !     (2*b_exp*b_exp)^n/b_exp/sqrt(Pi)*exp(-b_exp*b_exp*r*r) }
         implicit none         

         _REAL_, intent(IN)  :: delr, b_exp
         _REAL_, intent(OUT) :: c0, c1, c2, c3
         _REAL_ :: x, switch, d_switch_dx, fac, fact, delr2inv
         _REAL_ :: erf_x, erfc_x

         x        = delr*b_exp
         delr2inv = 1.0d0/(delr*delr)
         fac      = 2.0d0*b_exp*b_exp

         !later change to spline fit for switch and d_switch_dx
         call erfcfun(x, erfc_x)
         erf_x       = 1.0d0 - erfc_x
         switch      = erf_x
         d_switch_dx = 2.0d0/SQRT_PI*exp(-x*x)

         c0   = switch*delr*delr2inv
         fact = d_switch_dx*b_exp
         c1   = (c0 - fact)*delr2inv
         fact = fac*fact
         c2   = (3.0d0*c1 - fact)*delr2inv
         fact = fac*fact
         c3   = (5.0d0*c2 - fact)*delr2inv

      end subroutine Calc_c0_c1_c2_c3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!GAUSSIAN ELECTROSTATIC CORRECTION!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!the Gaussian correction to point charges -chgi*chgj*B0
!where B0 = erfc(Bij*Rij)/Rij where Bij^2 is the combined 
!Gaussian product exponent


      subroutine Calc_Dip_Gauss_Correc (i, j, delx, dely, delz, delr2inv, &
               delr, delr2, eelt, epol, dfx, dfy, dfz, dotjr, dotir, dotij, &
               dipole_jx, dipole_jy, dipole_jz, &
               dipole_ix, dipole_iy, dipole_iz, &
               dphii_dx_cor, dphii_dy_cor, dphii_dz_cor, &
               dphij_dx_cor, dphij_dy_cor, dphij_dz_cor, eed_cub)


         !del = r(j) - r(i) and df is force on j
         implicit none

         integer, intent(IN)     :: i, j
         _REAL_,  intent(IN)     :: delx, dely, delz, delr2inv, delr, delr2
         _REAL_,  intent(IN)     :: dotjr, dotir, dotij
         _REAL_,  intent(IN)     :: dipole_jx, dipole_jy, dipole_jz
         _REAL_,  intent(IN)     :: dipole_ix, dipole_iy, dipole_iz
         _REAL_,  intent(IN)     :: eed_cub(4,*)
         _REAL_,  intent(INOUT)  :: eelt, epol, dfx, dfy, dfz
         _REAL_,  intent(INOUT)  :: dphii_dx_cor, dphii_dy_cor, dphii_dz_cor
         _REAL_,  intent(INOUT)  :: dphij_dx_cor, dphij_dy_cor, dphij_dz_cor

         _REAL_ :: b_exp, b0, b1, b2, b3, term, termi, termj, termr


         dphii_dx_cor = 0.0d0
         dphii_dy_cor = 0.0d0
         dphii_dz_cor = 0.0d0
         dphij_dx_cor = 0.0d0
         dphij_dy_cor = 0.0d0
         dphij_dz_cor = 0.0d0

      !Pt Charge + Gaussian Dipoles
         !Q_tot(j) - dip(i)
         if ( r_polar(i) > delr) then
            b_exp = b_polar(i)
            call Calc_b0_b1_b2 (delr, b_exp, b0, b1, b2)
            term  = -Q_tot(j)*dotir
            termi = Q_tot(j)*b1
            termr = term*b2
            epol  = epol + term*b1
            dfx   = dfx + termi*dipole_ix + termr*delx
            dfy   = dfy + termi*dipole_iy + termr*dely
            dfz   = dfz + termi*dipole_iz + termr*delz

            !contribution to the field
            dphii_dx_cor = dphii_dx_cor - termi*delx
            dphii_dy_cor = dphii_dy_cor - termi*dely
            dphii_dz_cor = dphii_dz_cor - termi*delz
         end if

         !dip(j) - Q_tot(i)
         if ( r_polar(j) > delr) then   
            b_exp = b_polar(j)
            call Calc_b0_b1_b2 (delr, b_exp, b0, b1, b2)
            term  = Q_tot(i)*dotjr
            termj = -Q_tot(i)*b1
            termr = term*b2
            epol  = epol + term*b1
            dfx   = dfx + termj*dipole_jx + termr*delx
            dfy   = dfy + termj*dipole_jy + termr*dely
            dfz   = dfz + termj*dipole_jz + termr*delz
   
            !contribution to the field
            dphij_dx_cor = dphij_dx_cor - termj*delx
            dphij_dy_cor = dphij_dy_cor - termj*dely
            dphij_dz_cor = dphij_dz_cor - termj*delz
         end if         

         !dip-dip
         if ( r_polar(j) + r_polar(i) > delr) then   
            b_exp = b_polar(j)*b_polar(i)/ &
            SQRT(b_polar(j)*b_polar(j) + b_polar(i)*b_polar(i) )
            call Calc_b0_b1_b2_b3 (delr, b_exp, b0, b1, b2, b3)
            epol  = epol - dotij*b1 + dotjr*dotir*b2
            termi = -dotjr*b2
            termj = -dotir*b2
            termr = -dotij*b2 + dotjr*dotir*b3
            dfx   = dfx + termi*dipole_ix + termj*dipole_jx + termr*delx
            dfy   = dfy + termi*dipole_iy + termj*dipole_jy + termr*dely
            dfz   = dfz + termi*dipole_iz + termj*dipole_jz + termr*delz
   
            !contribution to the field
            termr  = dotjr*b2
            dphii_dx_cor = dphii_dx_cor - b1*dipole_jx + termr*delx
            dphii_dy_cor = dphii_dy_cor - b1*dipole_jy + termr*dely
            dphii_dz_cor = dphii_dz_cor - b1*dipole_jz + termr*delz

            !contribution to the field
            termr  = dotir*b2
            dphij_dx_cor = dphij_dx_cor - b1*dipole_ix + termr*delx
            dphij_dy_cor = dphij_dy_cor - b1*dipole_iy + termr*dely
            dphij_dz_cor = dphij_dz_cor - b1*dipole_iz + termr*delz
         end if   
      end subroutine Calc_Dip_Gauss_Correc




      subroutine Calc_b0_b1 (delr, b_exp, b0, b1, eed_cub)
      !b0 = erfc(b_exp*r)/r and bn+1 = -1/r dbn/dr
      !(same definition of b0 in direct sum)
      !bn = 1/r2{ (2n-1)*bn-1 + 
      !     (2*b_exp*b_exp)^n/b_exp/sqrt(Pi)*exp(-b_exp*b_exp*r*r) }
      !using same code as in direct sum
         implicit none         

         _REAL_, intent(IN)  :: delr, b_exp
         _REAL_, intent(IN)  :: eed_cub(4,*)
         _REAL_, intent(OUT) :: b0, b1
         _REAL_ :: x, switch, d_switch_dx, fact, delr2inv
         _REAL_ :: erfc_x

         x        = delr*b_exp
         delr2inv = 1.0d0/(delr*delr)

         !later change to spline fit for switch and d_switch_dx
         call erfcfun(x, erfc_x)
         switch      = erfc_x
         d_switch_dx = -2.0d0/SQRT_PI*exp(-x*x)

         b0   = switch*delr*delr2inv
         fact = d_switch_dx*b_exp
         b1   = (b0 - fact)*delr2inv
   

      end subroutine Calc_b0_b1

      subroutine Calc_b0_b1_b2 (delr, b_exp, b0, b1, b2)
      !b0 = erfc(b_exp*r)/r and bn+1 = -1/r dbn/dr
      !(same definition of b0 in direct sum)
      !bn = 1/r2{ (2n-1)*bn-1 + 
      !     (2*b_exp*b_exp)^n/b_exp/sqrt(Pi)*exp(-b_exp*b_exp*r*r) }
      !using same code as in direct sum
         implicit none         

         _REAL_, intent(IN)  :: delr, b_exp
         _REAL_, intent(OUT) :: b0, b1, b2
         _REAL_ :: x, switch, d_switch_dx, fact, fac, delr2inv
         _REAL_ :: erfc_x

         x        = delr*b_exp
         delr2inv = 1.0d0/(delr*delr)
         fac      = 2.0d0*b_exp*b_exp

         !later change to spline fit for switch and d_switch_dx
         call erfcfun(x, erfc_x)
         switch      = erfc_x
         d_switch_dx = -2.0d0/SQRT_PI*exp(-x*x)

         b0   = switch*delr*delr2inv
         fact = d_switch_dx*b_exp
         b1   = (b0 - fact)*delr2inv
         fact = fac*fact
         b2   = (3.0d0*b1 - fact)*delr2inv
   

      end subroutine Calc_b0_b1_b2

      subroutine Calc_b0_b1_b2_b3 (delr, b_exp, b0, b1, b2, b3)
      !b0 = erfc(b_exp*r)/r and bn+1 = -1/r dbn/dr
      !(same definition of b0 in direct sum)
      !bn = 1/r2{ (2n-1)*bn-1 + 
      !     (2*b_exp*b_exp)^n/b_exp/sqrt(Pi)*exp(-b_exp*b_exp*r*r) }
      !using same code as in direct sum
         implicit none         

         _REAL_, intent(IN)  :: delr, b_exp
         _REAL_, intent(OUT) :: b0, b1, b2, b3
         _REAL_ :: x, switch, d_switch_dx, fact, fac, delr2inv
         _REAL_ :: erfc_x

         x        = delr*b_exp
         delr2inv = 1.0d0/(delr*delr)
         fac      = 2.0d0*b_exp*b_exp

         !later change to spline fit for switch and d_switch_dx
         call erfcfun(x, erfc_x)
         switch      = erfc_x
         d_switch_dx = -2.0d0/SQRT_PI*exp(-x*x)

         b0   = switch*delr*delr2inv
         fact = d_switch_dx*b_exp
         b1   = (b0 - fact)*delr2inv
         fact = fac*fact
         b2   = (3.0d0*b1 - fact)*delr2inv
         fact = fac*fact
         b3   = (5.0d0*b2 - fact)*delr2inv
      end subroutine Calc_b0_b1_b2_b3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!MISC. I/O FUNCTIONS!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine Go_To_Line (unit_num, ref_line, n_char, filename)
         implicit none         

         integer, intent(IN) :: unit_num
         integer, intent(IN) :: n_char
         character(1000), intent(IN) :: ref_line
         character(1000), intent(IN) :: filename

         integer :: iostatv
         character(1000) :: line
         rewind unit_num

         do 
            read(unit_num, '(A1000)', iostat=iostatv) line
            if (line(1:n_char) .eq. ref_line(1:n_char) ) exit
            if (iostatv < 0) then
               print *, 'Missing: ', ref_line(1:24), ' in ', filename(1:30)
               STOP
            end if
         end do
      end subroutine Go_To_Line

      subroutine read_real_data (unit_num, n_per_line, n_tot, r_data)
         implicit none
         integer, intent(IN) :: unit_num, n_tot, n_per_line
         _REAL_, dimension(:), intent(INOUT) :: r_data

         integer :: i, r, a

         r = mod (n_tot, n_per_line)
         a = (n_tot-r)/n_per_line

         do i=1,a
            read(unit_num, *) r_data(n_per_line*(i-1)+1 : n_per_line*i)
         end do
         if (r > 0) then
            read(unit_num, *) r_data(n_per_line*a+1 : n_tot)
         end if   
      end subroutine read_real_data

      subroutine read_integer_data (unit_num, n_per_line, n_tot, i_data)
         implicit none
         integer, intent(IN) :: unit_num, n_tot, n_per_line
         integer, dimension(:), intent(INOUT) :: i_data

         integer :: i, r, a

         r = mod (n_tot, n_per_line)
         a = (n_tot-r)/n_per_line

         do i=1,a
            read(unit_num, *) i_data(n_per_line*(i-1)+1 : n_per_line*i)
         end do
         if (r > 0) then
            read(unit_num, *) i_data(n_per_line*a+1 : n_tot)
         end if   
      end subroutine read_integer_data

      subroutine read_interact_list (unit_num, n_per_line, n_tot, i_data)
         implicit none
         
         integer, intent(IN) :: unit_num, n_tot, n_per_line
         integer, dimension(:, :), intent(INOUT) :: i_data

         integer :: i, alloc_status
         integer, dimension(:), allocatable :: i_temp

         !read in 1D array from parm top file
         allocate(i_temp(1:2*n_tot), stat=alloc_status)
         if (alloc_status /= 0) then
            print *, 'Error in read_interact_list in pol_gauss.f'
            print *, 'allocating i_temp'
            STOP
         end if
         call read_integer_data (unit_num, n_per_line, 2*n_tot, i_temp)

         !convert to 2D array interaction list
         do i=1,n_tot
            i_data(1,i) = abs(i_temp(2*(i-1)+1))/3+1
            i_data(2,i) = abs(i_temp(2*(i-1)+2))/3+1
         end do
         deallocate(i_temp)
         
      end subroutine read_interact_list

      subroutine write_real_data (unit_num, n_tot, r_data)
         implicit none
         integer, intent(IN) :: unit_num, n_tot
         _REAL_, dimension(:), intent(INOUT) :: r_data

         integer :: i, r, a, n_per_line

         n_per_line = 5
         r = mod (n_tot, n_per_line)
         a = (n_tot-r)/n_per_line

         do i=1,a
            write(unit_num, '(5F10.4)') r_data(n_per_line*(i-1)+1 : n_per_line*i)
         end do
         if (r > 0) then
            do i = n_per_line*a+1, n_tot
               write(unit_num, '(F10.4)', advance='NO') r_data(i)
            end do
            write(unit_num, '(//)')
         end if   
      end subroutine write_real_data

      subroutine write_interact_list (unit_num, n_tot, i_data)
         implicit none
         integer, intent(IN) :: unit_num, n_tot
         integer, dimension(:, :), intent(INOUT) :: i_data
         integer :: i

         do i=1,n_tot
            write(unit_num, '(2I8,4X)', advance='NO') &
               i_data(1,i), i_data(2,i)
            if (mod (i, 3) .eq. 0) write(unit_num, *)
         end do
         write(unit_num, '(//)')
      end subroutine write_interact_list


      function erfc_inv(x)
      !crude asymptotic approximation to the inverse of erfc
      !based on erfc(x) ~ exp(-x^2)/sqrt(PI)/x 
         implicit none
         _REAL_, INTENT(IN) :: x
         _REAL_   :: erfc_inv, x0
   
         x0       = SQRT( -log(x) - 0.5d0*log(PI) )
         erfc_inv = SQRT( -log(x) - 0.5d0*log(PI) - log(x0))
      end function erfc_inv

end module pol_gauss

