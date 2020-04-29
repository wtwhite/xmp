// There is only one timer, so the TimerHandle value is unimportant and can be of any convenient type.
typedef unsigned TimerHandle;

void InitTimer(void);
void FinaliseTimer(void);
TimerHandle StartTimer(unsigned *flagPtr, unsigned secs, unsigned repeating);
void StopTimer(TimerHandle th);
