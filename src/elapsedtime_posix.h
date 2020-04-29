#ifndef __ELAPSEDTIME_POSIX_H
#define __ELAPSEDTIME_POSIX_H
#include <sys/time.h>
#include <time.h>

typedef struct timeval TimePoint;

TimePoint GetTimeNow(void);
double TimeDifference(TimePoint start, TimePoint end);
#endif	// #include guard
