#include <signal.h>
#include "elapsedtime_posix.h"

TimePoint GetTimeNow(void) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv;
}

double TimeDifference(TimePoint start, TimePoint end) {
	//HACK: We rely on suseconds_t being a signed integral type, so that if start.tv_usec > end.tv_usec,
	// we will get a negative double result, rather than a large positive double result.
	// Unfortunately I see no easy way of verifying this at compile time.
	return difftime(end.tv_sec, start.tv_sec) + (double) (end.tv_usec - start.tv_usec) / 1000000;
}
