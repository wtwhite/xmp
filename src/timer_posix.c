#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include "timer_posix.h"

struct TimerData {
	struct sigaction oldSigAct;
	unsigned *flagPtr;
	unsigned secs;
	unsigned repeating;
	unsigned inUse;
};

// There can be only one
static struct TimerData gTimer;

static void AlarmHandler(int sig);

void InitTimer(void) {
	gTimer.inUse = 0;
}

// It is safe to call this more than once.
void FinaliseTimer(void) {
	StopTimer(0);		// 0 is a dummy value
}

// It is safe to call this more than once for any given timer.
void StopTimer(TimerHandle th) {
	if (gTimer.inUse) {
		alarm(0);		// Cancels any active alarm.
		
		if (sigaction(SIGALRM, &gTimer.oldSigAct, NULL)) {
			fprintf(stderr, "sigaction() failed, aborting.\n");
			exit(1);
		}
		
		gTimer.inUse = 0;
	}
}

TimerHandle StartTimer(unsigned *flagPtr, unsigned secs, unsigned repeating) {
	struct sigaction newSigAct;
	
	if (gTimer.inUse) {
		fprintf(stderr, "Attempt to allocate more than one timer, aborting.\n");
		exit(1);
	}
	
	gTimer.inUse = 1;
	gTimer.flagPtr = flagPtr;
	gTimer.secs = secs;
	gTimer.repeating = repeating;
	newSigAct.sa_handler = AlarmHandler;
//	newSigAct.sa_mask = (sigset_t) 0;
	sigemptyset(&newSigAct.sa_mask);
	newSigAct.sa_flags = 0;
	
	if (sigaction(SIGALRM, &newSigAct, &gTimer.oldSigAct)) {
		fprintf(stderr, "sigaction() failed, aborting.\n");
		exit(1);
	}
	
	alarm(secs);
	
	return 0;		// Dummy value -- there's only one timer
}

void AlarmHandler(int sig) {
	*(gTimer.flagPtr) = 1;
	
	if (gTimer.repeating) {
		alarm(gTimer.secs);
	}
}
