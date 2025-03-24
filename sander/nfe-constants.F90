#include "nfe-config.h"

module nfe_constants

implicit none

private

! Modified by M Moradi
#ifdef NFE_REAL_IS_DOUBLE
NFE_REAL, public, parameter :: ZERO  = 0.d0
NFE_REAL, public, parameter :: ONE   = 1.d0
NFE_REAL, public, parameter :: TWO   = 2.d0
NFE_REAL, public, parameter :: THREE = 3.d0
NFE_REAL, public, parameter :: FOUR  = 4.d0
NFE_REAL, public, parameter :: kB = 1.9872041d-3  !Boltzmann's constant in (kcal/mol)/K
NFE_REAL, public, parameter :: PI    = 4 * atan (1.0_8)
#else
NFE_REAL, public, parameter :: ZERO  = 0.0
NFE_REAL, public, parameter :: ONE   = 1.0
NFE_REAL, public, parameter :: TWO   = 2.0
NFE_REAL, public, parameter :: THREE = 3.0
NFE_REAL, public, parameter :: FOUR  = 4.0
NFE_REAL, public, parameter :: kB = 1.9872041d-3  !Boltzmann's constant in (kcal/mol)/K
#endif /* NFE_REAL_IS_DOUBLE */
! Moradi end

integer, public, parameter :: STRING_LENGTH = 256

integer, public, parameter :: ERR_UNIT = SANDER_STDERR_UNIT
integer, public, parameter :: OUT_UNIT = SANDER_STDOUT_UNIT

integer, public, parameter :: ABMD_MONITOR_UNIT = SANDER_LAST_UNIT + 1
integer, public, parameter :: SMD_OUTPUT_UNIT = SANDER_LAST_UNIT + 2
integer, public, parameter :: PMD_OUTPUT_UNIT = SANDER_LAST_UNIT + 3
integer, public, parameter :: STSM_OUTPUT_UNIT = SANDER_LAST_UNIT + 4

#ifdef MPI
integer, public, parameter :: REM_MDIN_UNIT = SANDER_LAST_UNIT + 5
integer, public, parameter :: PMD_REMLOG_UNIT = SANDER_LAST_UNIT + 6
integer, public, parameter :: ABMD_REMLOG_UNIT = SANDER_LAST_UNIT + 7
integer, public, parameter :: BBMD_MONITOR_UNIT = SANDER_LAST_UNIT + 8
integer, public, parameter :: BBMD_LOG_UNIT = SANDER_LAST_UNIT + 9
integer, public, parameter :: STSM_REMLOG_UNIT = SANDER_LAST_UNIT + 10
#endif /* MPI */

integer, public, parameter :: EVEC_UNIT1 = SANDER_LAST_UNIT + 11
integer, public, parameter :: CRD_UNIT1  = SANDER_LAST_UNIT + 12
integer, public, parameter :: REF_UNIT1 = SANDER_LAST_UNIT + 13
integer, public, parameter :: IDX_UNIT1 = SANDER_LAST_UNIT + 14
integer, public, parameter :: SMD_CV_UNIT = SANDER_LAST_UNIT + 15
integer, public, parameter :: ABMD_CV_UNIT = SANDER_LAST_UNIT + 16
integer, public, parameter :: PMD_CV_UNIT = SANDER_LAST_UNIT + 17
integer, public, parameter :: BBMD_CV_UNIT = SANDER_LAST_UNIT + 18
integer, public, parameter :: STSM_CV_UNIT = SANDER_LAST_UNIT + 19

end module nfe_constants
