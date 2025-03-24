! <compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

module nfe_colvar_type

implicit none

private

integer, public, parameter :: COLVAR_ANGLE            = 1
integer, public, parameter :: COLVAR_TORSION          = 2
integer, public, parameter :: COLVAR_DISTANCE         = 3
integer, public, parameter :: COLVAR_MULTI_RMSD       = 4
integer, public, parameter :: COLVAR_R_OF_GYRATION    = 5
integer, public, parameter :: COLVAR_HANDEDNESS       = 6
integer, public, parameter :: COLVAR_N_OF_BONDS       = 7
integer, public, parameter :: COLVAR_N_OF_STRUCTURES  = 8
integer, public, parameter :: COLVAR_LCOD             = 9
integer, public, parameter :: COLVAR_COS_OF_DIHEDRAL  = 10
integer, public, parameter :: COLVAR_COM_ANGLE        = 11
integer, public, parameter :: COLVAR_COM_TORSION      = 12
integer, public, parameter :: COLVAR_COM_DISTANCE     = 13
integer, public, parameter :: COLVAR_PCA              = 14
integer, public, parameter :: COLVAR_SIN_OF_DIHEDRAL  = 15
integer, public, parameter :: COLVAR_PAIR_DIHEDRAL    = 16
integer, public, parameter :: COLVAR_PATTERN_DIHEDRAL = 17
integer, public, parameter :: COLVAR_DF_COM_DISTANCE  = 18
integer, public, parameter :: COLVAR_ORIENTATION_ANGLE= 19
integer, public, parameter :: COLVAR_ORIENTATION_PROJ = 20
integer, public, parameter :: COLVAR_SPINANGLE        = 21
integer, public, parameter :: COLVAR_TILT             = 22
integer, public, parameter :: COLVAR_QUATERNION0      = 23
integer, public, parameter :: COLVAR_QUATERNION1      = 24
integer, public, parameter :: COLVAR_QUATERNION2      = 25
integer, public, parameter :: COLVAR_QUATERNION3      = 26



type, public :: colvar_t

   integer :: type = -1

   integer,   pointer :: i(:) => null()
   NFE_REAL, pointer :: r(:) => null()
   
   integer :: tag ! (see nfe-cv-priv.*)

   ! avgcrd : average crd of the trajectory 
   ! r      : reference crd 
   ! evec   : eigenvector from PCA 
   NFE_REAL, pointer  :: avgcrd(:) => null() 
   NFE_REAL, pointer  :: evec(:) => null()
   NFE_REAL, pointer  :: axis(:) => null()
   integer,  pointer  :: q_index => null()
   ! state(:) stores the sate of reference part of ref.crd
   
   integer,  pointer :: state_ref(:) => null()
   integer,  pointer :: state_pca(:) => null()
   integer,  pointer :: ipca_to_i(:) => null()
 
end type colvar_t

double precision, public                 :: cv_min,cv_max,resolution
character(len = 256), public             :: cv_type
integer,dimension(20000),public          :: cv_i 
double precision,dimension(20000),public :: cv_r
integer, public                          :: cv_ni, cv_nr, refcrd_len

double precision,dimension(20000),public :: path, harm
integer, public                          :: npath, nharm, q_index
character(len = 256), public             :: path_mode, harm_mode, refcrd_file


double precision,dimension(4), public    :: anchor_position 
double precision,dimension(2), public    :: anchor_strength
double precision,dimension(3), public    :: axis


public  :: colvar
namelist / colvar /      cv_min, cv_max, resolution, cv_type, &
                         cv_i, cv_r, cv_ni, cv_nr, &
                         path, harm, npath, nharm, path_mode, harm_mode, &
                         anchor_position, anchor_strength, axis, q_index, &
                         refcrd_file

end module nfe_colvar_type
