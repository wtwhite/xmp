//#ifdef _WIN32
#ifdef _WIN32
#include "timer_win.h"
#else	// not _WIN32
#ifdef HAS_POSIX_TIMER
#include "timer_posix.h"
#else	// not HAS_POSIX_TIMER
//#error Could not find an appropriate timer_XXX.h file for your system.
#include "timer_dummy.h"
#endif	// not HAS_POSIX_TIMER
#endif	// not _WIN32
