#ifndef __MPIWORKER_H
#define __MPIWORKER_H
#include "common.h"
void SendNewBound(struct TreeData *td);
int WorkerMain(struct TreeData *td, struct tree *root, int idx, MPI_Status *mpiStatus, int waitingForWork);
void Worker(struct tree *root, struct TreeData *td);
#endif	// #include guard
