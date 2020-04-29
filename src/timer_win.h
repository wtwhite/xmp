#ifndef __TIMER_WIN_H
#define __TIMER_WIN_H
struct TimerData;			// Partial definition so that we can talk about pointer types
typedef struct TimerData *TimerHandle;

void InitTimer(void);
void __cdecl FinaliseTimer(void);
TimerHandle StartTimer(unsigned *flagPtr, unsigned secs, unsigned repeating);
void StopTimer(TimerHandle th);
#endif	// #include guard
