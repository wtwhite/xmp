#include <windows.h>
#include <mmsystem.h>		// Needed for timeGetTime()
//HACK: WTJW 11/6/2004: To get this damn thing to compile on the free BCC compiler, I need to
// work around Windows' settings of these macros...
//#ifndef __BORLANDC__
#if defined(__BORLANDC__) || defined(_MSC_VER)
#pragma message ("This file has an awful hack built into it so that it will compile with BCC and MS VC++.  Fix it!")
#undef LITTLEENDIAN
#undef BIGENDIAN
#define LITTLEENDIAN 1
#endif
#include "common.h"
#include "elapsedtime_win.h"

// In a sense, this doesn't really belong here since it deals with determining differences in
// time rather than setting up timers that will cause callbacks to run in the future, but it's
// easier to lump it in here than to create a new set of source and header files.
// The return value in this case will be the number of milliseconds since Windows started, although
// the caller should not depend on that, instead only depending on the result of calling TimeDifference()
// on 2 TimePoints.  Note this will overflow after around 49.7 days, but see TimeDifference().
TimePoint GetTimeNow(void) {
	return timeGetTime();
}

// Returns the difference between two time points in seconds.
// This will overflow if the times differ by more than 49.7 days or so, though it's not important
// whether the time counter wraps in between (twos-complement subtraction does the right thing).
double TimeDifference(TimePoint start, TimePoint end) {
	return (double) (end - start) / 1000;		// Since Windows timepoints are in ms
}
