#ifndef NFE_CONFIG_H
#define NFE_CONFIG_H

#if defined(LES) || defined(QMMM)
#  define DISABLE_NFE yes
#endif

#if defined(MPI) && defined(BINTRAJ) && !defined(DISABLE_NFE)
#  define NFE_ENABLE_BBMD sure,whynot
#endif

#ifndef _REAL_
#  include "../include/dprec.fh"
#  define NFE_REAL _REAL_
#  ifdef DPREC
#    define NFE_REAL_IS_DOUBLE indeed
#  endif
#endif /* _REAL_ */

#ifdef NFE_REAL_IS_DOUBLE
#  define NFE_TO_REAL(x) dble(x)
#else
#  define NFE_TO_REAL(x) real(x)
#endif /* NFE_REAL_IS_DOUBLE */

#define SANDER_STDERR_UNIT 6
#define SANDER_STDOUT_UNIT 6

!   /* EVB uses 75th (see files.h) */
#define SANDER_LAST_UNIT 77

#ifndef BINTRAJ
#  define NFE_NO_NETCDF entirely
#endif /* BINTRAJ */

#endif /* NFE_CONFIG_H */
