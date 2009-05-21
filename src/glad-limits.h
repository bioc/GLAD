#ifdef IS_MAC_OS
#include <limits.h>
#else
#include <values.h>
#endif

#ifndef MAXDOUBLE
#include <float.h>
#define MAXDOUBLE DBL_MAX
#endif


#ifndef MAXINT
#define MAXINT INT_MAX
#endif
