#include <stdlib.h>
#include <windows.h>
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
#include "timer_win.h"

//HACK: use a fixed-size global array of timers for easy management.
#define NUM_TIMERS 10

struct TimerData {
	HANDLE threadHandle;
	DWORD threadId;
	HANDLE threadTermEventHandle;
	unsigned *flagPtr;
	unsigned secs;
	unsigned repeating;
	unsigned inUse;
};

static struct TimerData gTimers[NUM_TIMERS];

static DWORD WINAPI TimerThreadFunc(LPVOID param);

void InitTimer(void) {
	unsigned i;
	
	for (i = 0; i < NUM_TIMERS; ++i) {
		gTimers[i].inUse = 0;
	}
	
	// Ensure that all threads are stopped before the program is exited.
	atexit(FinaliseTimer);
}

// This can safely be called more than once, so it is safe to pass it to atexit().
void __cdecl FinaliseTimer(void) {
	unsigned i;
	
	// Stop any timers that are still going.
	for (i = 0; i < NUM_TIMERS; ++i) {
		StopTimer(gTimers + i);		// Safe to call more than once
	}
}

TimerHandle StartTimer(unsigned *flagPtr, unsigned secs, unsigned repeating) {
	unsigned i;
	
	// Find an unused timer slot
	for (i = 0; i < NUM_TIMERS; ++i) {
		if (!gTimers[i].inUse) break;
	}
	
	if (i == NUM_TIMERS) {
		fprintf(stderr, "Attempt to allocate more than NUM_TIMERS=%u timers, aborting.\n", NUM_TIMERS);
		exit(1);
	}
	
	gTimers[i].flagPtr = flagPtr;
	gTimers[i].secs = secs;
	gTimers[i].repeating = repeating;
	gTimers[i].inUse = 1;
	gTimers[i].threadTermEventHandle = CreateEvent(
		NULL,	// pointer to security attributes  
		FALSE,	// flag for manual-reset event 
		FALSE,	// flag for initial state 
		NULL 	// pointer to event-object name  
		);
	gTimers[i].threadHandle = CreateThread(
		NULL,	// pointer to thread security attributes  
		1024,	// 1Kb stack size should be sufficient (but it will grow if necessary so don't worry)
		TimerThreadFunc,	// pointer to thread function 
		gTimers + i,	// argument for new thread 
		0,	// creation flags 
		&(gTimers[i].threadId) 	// pointer to returned thread identifier 
		);
	
	return gTimers + i;
}

// You MUST call StopTimer() on a particular timer once for every time you
// call StartTimer() to avoid running out of timers.  This function can safely
// be called more than once for any given timer.
void StopTimer(TimerHandle th) {
	if (th->inUse) {
		// Ask the timer thread to stop.  Could use TerminateThread() but this is cleaner.
		SetEvent(th->threadTermEventHandle);
		
		// Wait for the timer thread to stop.
		WaitForSingleObject(th->threadHandle, INFINITE);
		
		// Clean up the TimerData object
		CloseHandle(th->threadTermEventHandle);
		CloseHandle(th->threadHandle);
		th->inUse = 0;
	}
}

DWORD WINAPI TimerThreadFunc(LPVOID param) {
	struct TimerData *timer = (struct TimerData *) param;
	DWORD waitHandle;
	
	do {
		waitHandle = WaitForSingleObject(timer->threadTermEventHandle, timer->secs * 1000);
		
		if (waitHandle == WAIT_OBJECT_0) {
			// StopTimer() called: thread should now terminate.
			break;
		}
		
		// This is what we're here for, folks
		*(timer->flagPtr) = 1;
	} while (timer->repeating);
	
	return 0;
}
