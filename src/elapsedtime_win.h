#ifndef __ELAPSEDTIME_WIN_H
#define __ELAPSEDTIME_WIN_H
// It's pretty horrible that including just windef.h generates a bunch of spurious errors, forcing
// us to include the enormous windows.h just to get the definition for DWORD!  We can't even #define
// WIN32_LEAN_AND_MEAN, since the source file including us might not want that!
//#include <windef.h>			// SHOULD BE all that's needed for DWORD
#include <windows.h>			// Needed for DWORD; YUCK

typedef DWORD TimePoint;

TimePoint GetTimeNow(void);
double TimeDifference(TimePoint start, TimePoint end);
#endif	// #include guard
