#if !defined WIN64
#include <limits.h>
#ifndef IS_MAC_OS
#include <values.h>
#endif
#endif

#ifndef MAXDOUBLE
#include <float.h>
#define MAXDOUBLE DBL_MAX
#endif


#ifndef MAXINT
#define MAXINT INT_MAX
#endif
