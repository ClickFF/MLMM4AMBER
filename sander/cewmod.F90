#ifdef CEW
#include "../include/dprec.fh"
module cewmod

  use, intrinsic :: iso_c_binding
  implicit none

  public :: cew_prescf_wrapper
  public :: cew_postscf_wrapper
  public :: cew_free
  public :: cew_call_qmmm
  

  !
  ! Amber-consistent conversion factors
  ! Parts of the Ewald are performed by sander, cew, and quick;
  ! they need to use consistent conversions in order for real
  ! and reciprocal space interactions to properly cancel.
  ! The conversion factors are based on the conversions
  ! used within sander
  !
  _REAL_,parameter,private :: AMBERELE = 18.2223d0
  _REAL_,parameter,private :: BOHRS_TO_A = 0.529177249d0
  _REAL_,parameter,public  :: AU_PER_AMBER_CHARGE = 1.d0 / amberele
  _REAL_,parameter,public  :: AU_PER_AMBER_ANGSTROM = 1.d0 / bohrs_to_a
  _REAL_,parameter,public  :: AU_PER_AMBER_KCAL_PER_MOL = &
       & AU_PER_AMBER_CHARGE *  AU_PER_AMBER_CHARGE / AU_PER_AMBER_ANGSTROM
  _REAL_,parameter,public  :: AU_PER_AMBER_FORCE = &
       & AU_PER_AMBER_KCAL_PER_MOL / AU_PER_AMBER_ANGSTROM


  _REAL_,pointer,public :: cew_mmcharges(:) => NULL()
  _REAL_,pointer,public :: cew_mmcrd(:,:) => NULL()
  _REAL_,pointer,public :: cew_mmgrd(:,:) => NULL()
  _REAL_,pointer,public :: cew_qmcharges(:) => NULL()
  _REAL_,pointer,public :: cew_qmcrd(:,:) => NULL()
  _REAL_,pointer,public :: cew_qmgrd(:,:) => NULL()
  integer,public,save :: cew_nmm = 0

  logical,public,save :: use_cew = .false.
  

  !
  ! fortran interfaces to the c-functions exposed by the cew library 
  !
  
  interface
     
     subroutine cew_init( &
#ifdef MPI
          & commgroup, & 
#endif
          & p_natom, &
          & p_nquant, &
          & p_nquant_nlink, &
          & p_iqmatoms, &
          & p_charges, &
          & p_qm_charges, &
          & p_nucq, &
          & p_nfft1, &
          & p_nfft2, &
          & p_nfft3, &
          & p_order, &
          & p_ew_coeff, &
          & p_cut &
          & ) bind(c,name="cew_init_")
       use, intrinsic :: iso_c_binding
       implicit none
       
#ifdef MPI
       integer(kind=c_int),intent(in) :: commgroup 
#endif
       integer(kind=c_int),intent(in) :: p_natom
       integer(kind=c_int),intent(in) :: p_nquant
       integer(kind=c_int),intent(in) :: p_nquant_nlink
       integer(kind=c_int),intent(in) :: p_iqmatoms
       real(kind=c_double),intent(in) :: p_charges
       real(kind=c_double),intent(in) :: p_qm_charges
       real(kind=c_double),intent(in) :: p_nucq
       integer(kind=c_int),intent(in) :: p_nfft1
       integer(kind=c_int),intent(in) :: p_nfft2
       integer(kind=c_int),intent(in) :: p_nfft3
       integer(kind=c_int),intent(in) :: p_order
       real(kind=c_double),intent(in) :: p_ew_coeff
       real(kind=c_double),intent(in) :: p_cut

     end subroutine cew_init


 
     subroutine cew_prescf(crds,ucell,qmc,nmm,mmq,mmc) &
          & bind(c,name="cew_prescf_")
       use, intrinsic :: iso_c_binding
       implicit none

       real(kind=c_double),intent(in) :: crds
       real(kind=c_double),intent(in) :: ucell
       real(kind=c_double),intent(out) :: qmc
       integer(kind=c_int),intent(out) :: nmm
       real(kind=c_double),intent(out) :: mmq
       real(kind=c_double),intent(out) :: mmc

     end subroutine cew_prescf
     

     subroutine cew_postscf(qmg,mmg,e,frc) &
          & bind(c,name="cew_postscf_")
       use, intrinsic :: iso_c_binding
       implicit none

       real(kind=c_double),intent(in) :: qmg
       real(kind=c_double),intent(in) :: mmg
       real(kind=c_double),intent(inout) :: e
       real(kind=c_double),intent(inout) :: frc

     end subroutine cew_postscf
     
  end interface

  
  private


  
contains


  
  subroutine cew_initialize()
    
    use qmmm_module, only : qmmm_nml, qmmm_struct, qmmm_mpi
    use memory_module, only : x 
    implicit none

#include "extra.h"
    
#include "../include/memory.h"
    ! memory.h provides
    ! integer :: natom, l15

#include "ew_pme_recip.h"
    ! ew_pme_recip.h provides
    ! integer :: order,nfft1,nfft2,nfft3
    ! _REAL_ :: ew_coeff

    
    
    integer, pointer :: iqmatoms(:) => NULL()
    _REAL_, pointer :: nucq(:) => NULL()

    integer :: i,j
    logical :: found
    _REAL_ :: delta,beta,cutoff
    _REAL_ :: netcharge,netlinkcharge,netqmcharge
    _REAL_ :: netmmcharge
    
    if ( associated( cew_mmcharges ) ) then
       deallocate( cew_mmcharges )
    end if
    allocate( cew_mmcharges( natom ) )
    cew_mmcharges = 0.d0

    
    if ( associated( cew_mmcrd ) ) then
       deallocate( cew_mmcrd )
    end if
    allocate( cew_mmcrd( 3, natom ) )
    cew_mmcrd = 0.d0


    if ( associated( cew_mmgrd ) ) then
       deallocate( cew_mmgrd )
    end if
    allocate( cew_mmgrd( 3, natom ) )
    cew_mmgrd = 0.d0

    
    cew_nmm = 0

    
    if ( associated( cew_qmcharges ) ) then
       deallocate( cew_qmcharges )
    end if
    allocate( cew_qmcharges( qmmm_struct%nquant_nlink ) )
    cew_qmcharges = 0.d0

    
    if ( associated( cew_qmcrd ) ) then
       deallocate( cew_qmcrd )
    end if
    allocate( cew_qmcrd( 3, qmmm_struct%nquant_nlink ) )
    cew_qmcrd = 0.d0

    
    if ( associated( cew_qmgrd ) ) then
       deallocate( cew_qmgrd )
    end if
    allocate( cew_qmgrd( 3, qmmm_struct%nquant_nlink ) )
    cew_qmgrd = 0.d0

    
    allocate( iqmatoms( qmmm_struct%nquant_nlink ) )
    iqmatoms = 0

    
    allocate( nucq( qmmm_struct%nquant_nlink ) )
    nucq = 0.d0


    do i=1,qmmm_struct%nquant_nlink
       nucq(i) = qmmm_struct%iqm_atomic_numbers(i)
    end do
    
    do i=1,qmmm_struct%nquant
       iqmatoms(i) = qmmm_struct%iqmatoms(i)-1
    end do
    do i=1,qmmm_struct%nlink
       iqmatoms(i + qmmm_struct%nquant) = qmmm_struct%link_pairs(1,i)-1
    end do

    netqmcharge = 0.d0
    do i=1,qmmm_struct%nquant
       netqmcharge = netqmcharge + qmmm_struct%qm_resp_charges(i) &
            & * AU_PER_AMBER_CHARGE
       cew_qmcharges(i) = qmmm_struct%qm_resp_charges(i) &
            & * AU_PER_AMBER_CHARGE
    end do
    
    netlinkcharge = 0.d0
    do i=1,qmmm_struct%nlink
        netlinkcharge = netlinkcharge + &
             & qmmm_struct%mm_link_pair_resp_charges(i) * AU_PER_AMBER_CHARGE
        cew_qmcharges(qmmm_struct%nquant + i) = 0.d0
    end do

    netcharge = 0.d0
    do i=1,natom
       netcharge = netcharge + qmmm_struct%scaled_mm_charges(i)
    end do
    netmmcharge = netcharge - netlinkcharge - netqmcharge
    
#ifdef MPI
    if ( master ) then
#endif
       write(6,'(2A)')" cew: Checking partial charges.", &
            & " Modifications will be made if necessary."
       write(6,'(A)')" cew: On input, the following charge",&
            & " information was found."
       write(6,'(A,F10.5)')" cew: Sum of all partial charges:       ",&
            & netcharge
       write(6,'(A,F10.5)')" cew: Sum of qmmask partial charges:    ",&
            & netqmcharge
       write(6,'(A,F10.5)')" cew: Sum of link atom partial charges: ",&
            & netlinkcharge
       write(6,'(A,F10.5)')" cew: Sum of MM partial charges:        ", &
            & netcharge-netqmcharge-netlinkcharge
       write(6,'(A,I4)')" cew: Expected sum of qmmask charges:   ", &
            & qmmm_nml%qmcharge
       write(6,'(A,F10.5)')" cew: Expected sum of link atom charges:",0.d0


    
       if ( abs(qmmm_nml%qmcharge - netqmcharge) > 0.5d0 ) then
          write(6,'(A,I3,A,F9.5)')" cew: qmcharge = ",qmmm_nml%qmcharge, &
               & " but MM charge of QM region is ",netqmcharge
          write(6,'(2A)')" cew: Either there is an error in qmmask or", &
               & " you should modify the", &
               & " charges in the parmeter file to make them agree"
          call mexit(6,1)
       end if
#ifdef MPI
    end if
#endif
    
    delta = ( qmmm_nml%qmcharge - netqmcharge ) / qmmm_struct%nquant


     
#ifdef MPI
    if ( master ) then
#endif
       write(6,'(A,F9.5,A)')" cew: Each QM atom charge shifted by ",delta,&
            & " so the underlying MM charges match the qmcharge" 
#ifdef MPI
    end if
#endif

    
    do i=1,qmmm_struct%nquant
       cew_qmcharges(i) = cew_qmcharges(i) + delta
    end do

    do i=1,qmmm_struct%nquant_nlink
       qmmm_struct%scaled_mm_charges( iqmatoms(i)+1 ) = cew_qmcharges(i)
    end do


    
    delta = (netlinkcharge+netqmcharge-qmmm_nml%qmcharge) &
         & / (natom - qmmm_struct%nquant_nlink)

!     if ( qmmm_struct%nlink > 0 ) then
       
     
! #ifdef MPI
!        if ( master ) then
! #endif
!           write(6,'(A)')" cew: Setting link atom charges to zero"
!           write(6,'(A,I3,A,F9.5)')" cew: The sum of the ",qmmm_struct%nlink, &
!                & " link atom charges was ",netlinkcharge
!           write(6,'(A,I7,A,F9.5,A)')"The ",natom-qmmm_struct%nquant_nlink, &
!                & " MM charges (non-link, non-QM) will be shifted by ", &
!                & delta, " due to the zeroing of the link atom charges"
! #ifdef MPI
!        end if
! #endif

!     end if

    ! netmmcharge = 0.d0
    ! do i=1,natom
    !    found=.false.
    !    do j=1,qmmm_struct%nquant_nlink
    !       if ( i == iqmatoms(j)+1 ) then
    !          found = .true.
    !       end if
    !    end do
    !    if ( .not. found ) then
    !       !qmmm_struct%scaled_mm_charges(i) = qmmm_struct%scaled_mm_charges(i) &
    !       !     & +  delta
    !       netmmcharge = netmmcharge + qmmm_struct%scaled_mm_charges(i)
    !    end if
    ! end do


    
    do i=1,natom
       x(l15-1+i) = qmmm_struct%scaled_mm_charges(i) / AU_PER_AMBER_CHARGE
    end do

    netcharge = 0.d0
    do i=1,natom
       netcharge = netcharge + qmmm_struct%scaled_mm_charges(i)
    end do
    
    netqmcharge = 0.d0
    do i=1,qmmm_struct%nquant
       netqmcharge = netqmcharge + cew_qmcharges(i)
    end do
    
    netlinkcharge = 0.d0
    do i=1,qmmm_struct%nlink
       netlinkcharge = netlinkcharge + cew_qmcharges(qmmm_struct%nquant+i)
    end do

    
    !write(6,'(a,4f20.10)')"netcharge",netcharge,netqmcharge,netlinkcharge,&
    !     & netmmcharge
    

#ifdef MPI
    if ( master ) then
#endif
       
       write(6,'(A)')" cew: Charge information after modifications."
       write(6,'(A,F10.5)')" cew: Sum of all partial charges:       ",&
            & netcharge
       write(6,'(A,F10.5)')" cew: Sum of qmmask partial charges:    ",&
            & netqmcharge
       write(6,'(A,F10.5)')" cew: Sum of link atom partial charges: ",&
            & netlinkcharge
       write(6,'(A,F10.5)')" cew: Sum of MM partial charges:        ", &
            & netcharge-netqmcharge-netlinkcharge

#ifdef MPI
    end if
#endif


    beta = ew_coeff / AU_PER_AMBER_ANGSTROM
    cutoff = qmmm_nml%qmcut * AU_PER_AMBER_ANGSTROM
    
    call cew_init( &
#ifdef MPI
          & qmmm_mpi%commqmmm, & 
#endif
          & natom, &
          & qmmm_struct%nquant, &
          & qmmm_struct%nquant_nlink, &
          & iqmatoms(1), &
          & qmmm_struct%scaled_mm_charges(1), &
          & cew_qmcharges(1), &
          & nucq(1), &
          & nfft1, &
          & nfft2, &
          & nfft3, &
          & order, &
          & beta, &
          & cutoff )


    
    
    do i=1,qmmm_struct%nquant_nlink
       j = iqmatoms(i)+1
       x(l15-1+j) = 0.d0
    end do

    deallocate( iqmatoms )
    deallocate( nucq )
    
  end subroutine cew_initialize


  
  subroutine cew_free
    implicit none
    if ( associated( cew_mmcharges ) ) then
       deallocate( cew_mmcharges )
    end if
    if ( associated( cew_qmcharges ) ) then
       deallocate( cew_qmcharges )
    end if
    if ( associated( cew_mmcrd ) ) then
       deallocate( cew_mmcrd )
    end if
    if ( associated( cew_mmgrd ) ) then
       deallocate( cew_mmgrd )
    end if
    if ( associated( cew_qmcrd ) ) then
       deallocate( cew_qmcrd )
    end if
    if ( associated( cew_qmgrd ) ) then
       deallocate( cew_qmgrd )
    end if
    cew_nmm = 0
  end subroutine cew_free

  
  subroutine cew_prescf_wrapper()
    use memory_module, only : x
    implicit none
    
#include "../include/memory.h"
    ! integer :: natom, l15

    call cew_prescf_worker( natom, x(lcrd) )

  end subroutine cew_prescf_wrapper



  subroutine cew_postscf_wrapper( escf )
    use memory_module, only : x
    implicit none
    
#include "../include/memory.h"
    ! integer :: natom, l15
    _REAL_,intent(inout) :: escf

    call cew_postscf_worker( natom, x(lcrd), x(lforce), escf )

  end subroutine cew_postscf_wrapper



  

  subroutine cew_prescf_worker( nat, crd )

    use constants, only : zero
    use qmmm_module, only : qmmm_struct, qm2_struct, qmmm_nml, qmewald, &
         & qmmm_scratch, qmmm_mpi
    use nblist, only: ucell
    use memory_module, only : x

    use quick_cew_module, only : new_quick_cew
    
    implicit none
    
#include "../include/memory.h"
    ! memory.h provides
    ! integer :: natom, l15

#include "ew_pme_recip.h"
    ! ew_pme_recip.h provides
    ! integer :: order,nfft1,nfft2,nfft3
    ! _REAL_ :: ew_coeff

#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#include "extra.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
   integer ierr
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

   integer, intent(in) :: nat
   _REAL_, intent(inout)  :: crd(3,nat)
   
   integer :: i,qm_temp_count
   _REAL_ :: beta
!   _REAL_,pointer,save :: savecrd(:,:) => NULL()
!   integer :: k
   

   if ( qmmm_struct%qm_mm_first_call ) then

!      allocate(savecrd(3,nat))
!      savecrd = crd
      
      qmmm_struct%qm_mm_first_call = .false.
      
      call cew_initialize()
      
      if ( associated( qmmm_scratch%qm_real_scratch ) ) then
         deallocate( qmmm_scratch%qm_real_scratch )
      end if
      allocate( qmmm_scratch%qm_real_scratch(3*natom) )
      qmmm_scratch%qm_real_scratch = zero


      if ( qmmm_nml%qmtheory%ISQUICK ) then
         beta = ew_coeff / AU_PER_AMBER_ANGSTROM
         call new_quick_cew(beta,size(cew_qmcharges),cew_qmcharges)
      end if
      
      !do qm_temp_count = 1, qmmm_struct%nquant_nlink
      !   charges(qmmm_struct%iqmatoms(qm_temp_count)) = zero
      !end do

   end if

!   do i=1,nat
!      do k=1,3
!         if ( abs(crd(k,i)-savecrd(k,i)) > 5.d-7 ) then
!            write(6,'(a,i6,i4,3f15.7)')"DIFF",i,k, &
!                 & crd(k,i),savecrd(k,i),crd(k,i)-savecrd(k,i)
!         end if
!      end do
!   end do

   
   !Replace the MM atom's coordinates 
   call adj_mm_link_pair_crd(crd) 
   
   qmmm_scratch%qm_real_scratch = zero
   
   cew_mmcharges = zero
   cew_mmcrd = zero
   cew_qmcrd = zero
   cew_mmgrd = zero
   cew_qmgrd = zero
   
   call cew_prescf(crd(1,1),ucell(1,1), &
        & cew_qmcrd(1,1), &
        & cew_nmm,cew_mmcharges(1),cew_mmcrd(1,1))
   
   ! now call qm/mm energy and force
   
   
 end subroutine cew_prescf_worker


  
 subroutine cew_postscf_worker( nat, crd, frc, escf )

    use constants, only : zero
    use qmmm_module, only : qmmm_struct, qm2_struct, qmmm_nml, qmewald, &
         & qmmm_scratch, qmmm_mpi
    use nblist, only: ucell
    use memory_module, only : x 

    implicit none
    
#include "../include/memory.h"
    ! memory.h provides
    ! integer :: natom, l15

#include "ew_pme_recip.h"
    ! ew_pme_recip.h provides
    ! integer :: order,nfft1,nfft2,nfft3
    ! _REAL_ :: ew_coeff

#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#include "extra.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
   integer ierr
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif


    integer, intent(in) :: nat
    _REAL_, intent(inout) :: crd(3,nat)
    _REAL_, intent(inout) :: frc(3,nat)
    _REAL_, intent(inout) :: escf
    
    integer :: qm_temp_count,i,mm_no,qm_no,lnk_no

    _REAL_ :: forcemod(3)



   qmmm_scratch%qm_real_scratch = zero

   call cew_postscf( cew_qmgrd(1,1), cew_mmgrd(1,1), escf, &
        & qmmm_scratch%qm_real_scratch(1) )

   

   !do i=1,16
   !   write(6,'(a,i3,3f9.4,3f9.4)')"prelink ",i,crd(1:3,i),qmmm_scratch%qm_real_scratch(3*i-2:3*i)
   !end do
   
   
   forcemod = 0.d0
   
   !Restore the MM atom's coordinates 
   call rst_mm_link_pair_crd(crd)

   do i=1,qmmm_struct%nlink
      mm_no = qmmm_struct%link_pairs(1,i)  !location of atom in x array
      lnk_no = qmmm_struct%link_pairs(2,i) !Nquant number of QM atom bound to link atom
      qm_no = qmmm_struct%iqmatoms(lnk_no)
      !Note this routine uses the flink in the form -flink.
      call distribute_lnk_f(forcemod, &
           qmmm_scratch%qm_real_scratch(3*mm_no-2:3*mm_no), &
           crd(1,mm_no), crd(1,qm_no), qmmm_nml%lnk_dis)
      
      frc(1:3,mm_no) = frc(1:3,mm_no) + forcemod(1:3)
      frc(1:3,qm_no) = frc(1:3,qm_no) &
           + qmmm_scratch%qm_real_scratch(3*mm_no-2:3*mm_no) &
           - forcemod(1:3)
      !Zero out the link atom force in the scratch array.
      qmmm_scratch%qm_real_scratch(3*mm_no-2) = zero
      qmmm_scratch%qm_real_scratch(3*mm_no-1) = zero
      qmmm_scratch%qm_real_scratch(3*mm_no) = zero
   end do
   
   do i = 1, nat
      frc(1,i)=frc(1,i)+qmmm_scratch%qm_real_scratch(3*i-2)
      frc(2,i)=frc(2,i)+qmmm_scratch%qm_real_scratch(3*i-1)
      frc(3,i)=frc(3,i)+qmmm_scratch%qm_real_scratch(3*i)
   end do


   !do i=1,16
   !   write(6,'(a,i3,3f9.4,3f9.4)')"postlink",i,crd(1:3,i),frc(1:3,i)
   !end do
   
   
 end subroutine cew_postscf_worker

   
  subroutine cew_call_qmmm( escf )
    
    use qmmm_module, only : qmmm_nml, qmmm_struct
    use quick_module, only: get_quick_qmmm_forces
    
    implicit none

    _REAL_,intent(out) :: escf

    _REAL_,pointer :: mmdata(:,:) => NULL()
    integer :: i

    escf = 0.d0
    
    if (qmmm_nml%qmtheory%ISQUICK) then

       allocate( mmdata(4,cew_nmm) )
       mmdata = 0.d0
       do i=1,cew_nmm
          mmdata(1:3,i) = cew_mmcrd(1:3,i)
          mmdata(4,i) = cew_mmcharges(i)
       end do
       
       call get_quick_qmmm_forces(qmmm_struct%nquant_nlink, &
                                  cew_qmcrd, &
                                  qmmm_struct%iqm_atomic_numbers, &
                                  cew_nmm, mmdata, escf, &
                                  cew_qmgrd, cew_mmgrd, use_cew)

       deallocate( mmdata )
       
    else
       
       write(6,*)"ERROR in cewmod.F90 (cew_call_qmmm) : unimplemented QM theory"
       call mexit(6,1)
       
    end if
    
  end subroutine cew_call_qmmm

    
end module cewmod
#endif