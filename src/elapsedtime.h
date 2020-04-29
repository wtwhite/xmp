#ifdef _WIN32
#include "elapsedtime_win.h"
#else	// not _WIN32
#ifdef HAS_GETTIMEOFDAY
#include "elapsedtime_posix.h"
#else	// not HAS_GETTIMEOFDAY
#error Could not find an appropriate elapsedtime_XXX.h file for your system.
#endif	// not HAS_GETTIMEOFDAY
#endif	// not _WIN32
