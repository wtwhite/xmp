#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "common.h"
#include <dbgprint.h>
#include "mpiworker.h"
#include "mpiboss.h"
#include "yanbader.h"
#include "measure.h"
#include "maxmin.h"
#ifdef TARGETMULTI
#include "mpi.h"
#endif	// TARGETMULTI

#ifdef TARGETMULTI
void NewBoundReceived(struct TreeData *td) {
	DBGPRINT3("W%04d: Received bound update of %u.\n", td->rank, *td->receiveBuf);
	if (*td->receiveBuf < td->bound) {
		DBGPRINT4("W%04d: Newly received bound of %u improves on our previous bound of %u!\n", td->rank, *td->receiveBuf, td->bound);
		td->bound = *td->receiveBuf;

		DeleteOldTrees(td);
		td->pendingBoundBroadcast = 0;	// The received UB is better than our own, even after improvements since we last broadcast it
	}
}

// This function assumes that it is safe to call MPI_CHOSEN_Isend(), i.e. there is not already a MSG_NEWBOUND being
// sent by us that has not completed yet.
void SendNewBound(struct TreeData *td) {
	static unsigned ubToSend;		// Must be static!
	assert(td->pendingBoundBroadcast);
	DBGPRINT3("W%04d: Sending new bound of %u.\n", td->rank, td->bound);
	assert(td->mpiRequests[WCR_SND_NEWBOUND] == MPI_REQUEST_NULL);
	// Length-1 message indicates NEWBOUND.
	// It's safer to copy td->bound into a static variable and use that as the buffer --
	// otherwise it's (remotely) possible that the send will occur just as the bound is being updated.
	ubToSend = td->bound;
	MPI_CHOSEN_Isend(&ubToSend, 1, MPI_UNSIGNED, 0, MSG_WORKERREQUEST, MPI_COMM_WORLD, td->mpiRequests + WCR_SND_NEWBOUND);
	td->pendingBoundBroadcast = 0;
}

// This function assumes that it is safe to call MPI_CHOSEN_Isend(), i.e. there is not already a MSG_NEWBOUND being
// sent by us that has not completed yet.
void AskForNewJob(struct TreeData *td) {
	static unsigned dummy = -1;		// Must be static!
	assert(td->pendingAskForJob);
	DBGPRINT2("W%04d: Asking for new job.\n", td->rank);
	// Length-0 message indicates ASKFORJOB.
	MPI_CHOSEN_Isend(&dummy, 0, MPI_UNSIGNED, 0, MSG_WORKERREQUEST, MPI_COMM_WORLD, td->mpiRequests + WCR_SND_ASKFORJOB);
	td->pendingAskForJob = 0;
}

// This is called in one of the following situations:
// (1) We have run out of work.
// (2) A request for work sent by us has completed.
// (3) A new bound we sent has completed.
//
// In all cases the action is now actually the same: initiate the next thing
// on our lists.  This must always be sending another UB update if that is pending; it's
// important that we only ask for a new job once sends for all UBs have been initiated.
// This is needed to guarantee that we never initiate a UB send when we already have an
// outstanding request for work -- since if the boss runs out of work, it will then respond
// to the work request with a "show's over" message and never receive the UB, leading to a hang.
// And because the only way for us to finish is by sending a request for a new job and having
// it turned down, this means that our final UB is always sent to the boss before we finish.
//
// All this requires thinking about things more like a finite state machine.  We can only initiate requests for
// work once all UB sends have been initiated.  (Since the non-overtaking property of MPI guarantees that these
// requests will be received in order of *initiation*.  And even though the boss's MPI_Waitsome() might complete both
// receives at the same, UBs are always processed first there, so that's OK.)
//
// MPI will reset the correct MPI_Request to MPI_REQUEST_NULL automatically.
// We can only have one outstanding bound update at a time, so it's possible that our
// UB has improved between the last improvement for which a MSG_NEWBOUND was sent by us
// and receipt of this event.
static void SendNewBoundOrAskForWork(struct TreeData *td) {
	if (td->pendingBoundBroadcast) {
		if (td->mpiRequests[WCR_SND_NEWBOUND] == MPI_REQUEST_NULL) {
			// Initiate send of new UB
			SendNewBound(td);
		}
	} else {
		// Not waiting on any UB-related stuff, but we might need to ask for work.
		if (td->mpiRequests[WCR_SND_ASKFORJOB] == MPI_REQUEST_NULL && td->pendingAskForJob) {
			// It's safe to make the request now.
			AskForNewJob(td);
		}
	}
}

// Returns 0 to indicate that the boss has told us that there is no more work.  Otherwise returns 1.
int WorkerMain(struct TreeData *td, struct tree *root, int idx, MPI_Status *mpiStatus, int waitingForWork) {
	int count;
	int i, j;
	static unsigned noWork;		// This dummy value is not actually used, but setting it to all-bits-on may simplify debugging
	static unsigned dummy;		// Used as a buffer for 0-length messages

	DBGPRINT4("W%04d: Received instructions.  idx=%d.  waitingForWork=%d\n", td->rank, idx, waitingForWork);

	// Because we now just use a single MPI_ANY_TAG request for all receives, we need to check the
	// tag instead of the index.
	if (idx == WCR_RCV_ALL && mpiStatus->MPI_TAG == MSG_NEWJOB) {
		assert(waitingForWork);		// It should be impossible to receive this message during BranchAndBound().

		MPI_Get_count(mpiStatus, MPI_UNSIGNED, &count);
		if (count == 1) {
			// We're finished.
			DBGPRINT2("W%04d: Boss says we're finished.\n", td->rank);
			return 0;
		} else {
			// More work to do!
			DBGPRINT3("W%04d: Received new job of size %d.\n", td->rank, count / 2);

			// Now that we no longer receive directly into td->workStack, need to copy the data there.
			memcpy(td->workStack + 3, (unsigned (*)[2]) td->receiveBuf, count * sizeof (unsigned));

			// Repost the receive.
			//MPI_Irecv(jobBuf, td->numTaxa * 2, MPI_UNSIGNED, 0, MSG_NEWJOB, MPI_COMM_WORLD, td->mpiRequests + WCR_RCV_NEWJOB);
			//MPI_Irecv(td->workStack + 3, td->numTaxa * 2, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, td->mpiRequests + WCR_RCV_ALL);
			MPI_Irecv(td->receiveBuf, td->numTaxa * 2, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, td->mpiRequests + WCR_RCV_ALL);

			for (j = 0; j < count / 2; ++j) {
				DBGPRINT5("W%04d: New job: Level %d: %d..%d\n", td->rank, j, td->workStack[j + 3][0], td->workStack[j + 3][1]);
			}

			//TODO: Actually call BranchAndBound() here.  Need to somehow pass in all MPI_Request state.
			// Skip over the 1st two elements of workStack -- i.e. an empty job means "start at the root".
			td->jobTreeSize = 2 + count / 2;
			//memcpy(td->workStack + 3, jobBuf, count * sizeof (unsigned));
			BranchAndBound(root, 3, td);

			// Record the fact that we need to ask for more work.
			td->pendingAskForJob = 1;

			// Use identical logic to handle all three events: (1) running out of work, (2) send of new UB completes, (3) send of request for new work completes.
			SendNewBoundOrAskForWork(td);
		}
	} else if (idx == WCR_RCV_ALL && mpiStatus->MPI_TAG == MSG_NEWBOUND) {
		NewBoundReceived(td);

		// Repost the receive.
		MPI_Irecv(td->receiveBuf, td->numTaxa * 2, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, td->mpiRequests + WCR_RCV_ALL);
	} else if (idx == WCR_RCV_ALL && mpiStatus->MPI_TAG == MSG_STEALJOB) {
		// Find the first available post slot and use that.
		for (i = 0; td->mpiRequests[i + MAX_WCR] != MPI_REQUEST_NULL && i < td->nWorkers; ++i) {
			// Do nothing
		}

		// We should have found a slot, because it should not be possible to receive more than td->nWorkers
		// steal requests before the boss waits for a response from us.
		assert(i < td->nWorkers);
		assert(!td->stealBuf[i]);

		// Need to keep track of the number of outstanding steal responses so that we don't bite off more than we can chew.
		++td->nStealsFromUs;
		assert(td->nStealsFromUs <= td->nWorkers);

		if (waitingForWork) {
			DBGPRINT2("W%04d: Received steal request while waiting for work; will deny it.\n", td->rank);
			MPI_CHOSEN_Isend(&noWork, 1, MPI_UNSIGNED, 0, MSG_NEWJOB, MPI_COMM_WORLD, td->mpiRequests + i + MAX_WCR);
		} else {
			DBGPRINT2("W%04d: Received steal request while working; will satisfy it.\n", td->rank);
			td->stealBuf[i] = (unsigned (*)[2]) malloc((td->numTaxa + 1) * 2 * sizeof (unsigned));

			// Find a job ;)
			for (j = 3; td->workStack[j][0] + 1 == td->workStack[j][1]; ++j) {
				td->stealBuf[i][j - 3][0] = td->workStack[j][0];
				td->stealBuf[i][j - 3][1] = td->workStack[j][1];
				DBGPRINT5("W%04d: Job to be stolen: Level %d: %d..%d\n", td->rank, j - 3, td->stealBuf[i][j - 3][0], td->stealBuf[i][j - 3][1]);
			}

			td->stealBuf[i][j - 3][0] = td->workStack[j][1] - 1;
			td->stealBuf[i][j - 3][1] = td->workStack[j][1];
			DBGPRINT5("W%04d: Job to be stolen: Level %d: %d..%d\n", td->rank, j - 3, td->stealBuf[i][j - 3][0], td->stealBuf[i][j - 3][1]);

			DBGPRINT4("W%04d: Found a job of size %d to send back.  Will use slot %d to hold it.\n", td->rank, j - 2, i);
			--td->workStack[j][1];
			MPI_CHOSEN_Isend(td->stealBuf[i], (j - 2) * 2, MPI_UNSIGNED, 0, MSG_NEWJOB, MPI_COMM_WORLD, td->mpiRequests + i + MAX_WCR);
		}

		// Repost the receive.
		MPI_Irecv(td->receiveBuf, td->numTaxa * 2, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, td->mpiRequests + WCR_RCV_ALL);
	} else if (idx == WCR_SND_ASKFORJOB || idx == WCR_SND_NEWBOUND) {
		// See the comments at the top of WorkerTransition().
		SendNewBoundOrAskForWork(td);
	} else if (idx >= MAX_WCR) {
		// A steal response we sent has completed.  Don't need to do anything --
		// MPI will reset the correct MPI_Request to MPI_REQUEST_NULL automatically.
		assert(td->mpiRequests[idx] == MPI_REQUEST_NULL);
		if (td->stealBuf[idx - MAX_WCR]) {
			// This resulted from a successful steal from us: deallocate the buffer.
			// (A failed steal request does not allocate any buffer.)
			free(td->stealBuf[idx - MAX_WCR]);
			td->stealBuf[idx - MAX_WCR] = NULL;
		}
		--td->nStealsFromUs;
	} else {
		DBGPRINT5("W%04d: Unexpected event! MPI_ERROR=%d, MPI_SOURCE=%d, MPI_TAG=%d\n", td->rank,
			mpiStatus->MPI_ERROR,
			mpiStatus->MPI_SOURCE,
			mpiStatus->MPI_TAG
		);
		assert(!"Unexpected message tag received!");
	}

	return 1;
}

//HACK: Because we need to wait on many different requests simultaneously, we need to keep them all in a single
// array.  Currently we have specific requests positioned at specific indices, which is very fragile.
void Worker(struct tree *root, struct TreeData *td) {
	unsigned dummy;		// Used as a buffer for 0-length messages
	int idx;
	MPI_Status mpiStatus;
	int i;
	MPI_Status *mpiStatuses;		// Needed for MPI_Waitall() at the end.
	char *treeBuf;
	int nBytesRead;
	unsigned nTreesRemaining;			// Only used if td->maxTrees != 0
	unsigned nTreesToSend;
	unsigned nTreesSent = 0;			// Only used if td->maxTrees != 0

	DBGPRINT3("W%04d: I'm a worker, with rank %d.\n", td->rank, td->rank);

	td->nMpiRequests = MAX_WCR + td->nWorkers;
	td->mpiRequests = (MPI_Request *) malloc(td->nMpiRequests * sizeof (MPI_Request));
	td->nStealsFromUs = 0;

	// This is only necessary for assert() checks.
	for (i = 0; i < td->nMpiRequests; ++i) {
		td->mpiRequests[i] = MPI_REQUEST_NULL;
	}

	td->receiveBuf = (unsigned *) malloc(td->numTaxa * 2 * sizeof (unsigned));
	td->stealBuf = (unsigned (**)[2]) malloc(td->nWorkers * 2 * sizeof (unsigned (*)[2]));
	memset(td->stealBuf, 0, td->nWorkers * 2 * sizeof (unsigned (*)[2]));

	// Set up outstanding requests.
	// This wastes some space on request types that we don't wait for, but makes the code easier to read and write.
	// Instead of 3 separate receives, try posting a single one with MPI_ANY_TAG -- this should guarantee that messages sent
	// arrive in order, which is necessary so we know whether a given MSG_NEWBOUND was sent before or after the final length-1
	// MSG_NEWJOB.  Note that we still don't care about the actual bound sent with MSG_NEWBOUND -- we use a 0-length message here for that.
	MPI_Irecv(td->receiveBuf, td->numTaxa * 2, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, td->mpiRequests + WCR_RCV_ALL);

	// Ask for work.
	//MPI_CHOSEN_Isend(&dummy, 0, MPI_UNSIGNED, 0, MSG_ASKFORJOB, MPI_COMM_WORLD, td->mpiRequests + WCR_SND_ASKFORJOB);
	td->pendingAskForJob = 1;
	AskForNewJob(td);

	// Initialise slot for completion of sending a bound update.  There can be at most 1 outstanding
	// bound update at a time.
	td->mpiRequests[WCR_SND_NEWBOUND] = MPI_REQUEST_NULL;

	// Initialise slots for handling completion of send requests (e.g. responses to steal attempts)
	for (i = 0; i < td->nWorkers; ++i) {
		td->mpiRequests[i + MAX_WCR] = MPI_REQUEST_NULL;			// Will be used for sending back stolen jobs later.
	}

	while (1) {
		//HACK: We could make the following slightly more efficient by remembering what the highest-numbered in-use
		// request is, and only waiting on that many requests.
		DBGPRINT2("W%04d: Waiting for instructions.\n", td->rank);
		// There is a tricky race condition whereby if everyone has stolen from us, if we're not careful we can
		// accidentally start processing another steal request from the next round of steals before waiting on
		// our sends of stolen jobs to complete, meaning there won't be any slots free!  If we have any outstanding
		// sends of stolen jobs, we can wait on them without deadlocking, so let's do that.
		if (td->nStealsFromUs < td->nWorkers) {
			MPI_Waitany(td->nMpiRequests, td->mpiRequests, &idx, &mpiStatus);
		} else {
			// Hold your horses!  We need to wait for a send to complete before doing anything else.
			MPI_Waitany(td->nMpiRequests - MAX_WCR, td->mpiRequests + MAX_WCR, &idx, &mpiStatus);
			idx += MAX_WCR;
		}

		if (!WorkerMain(td, root, idx, &mpiStatus, 1)) {
			break;
		}
	}

	// Although we could omit the special-case waiting for the final bound and steal requests since
	// they will be waited on by the final MPI_Waitall() that waits for outstanding sends to complete,
	// it may simplify debugging to wait on them separately, and the performance penalty is tiny.
	DBGPRINT2("W%04d: Waiting for final bound update.\n", td->rank);
	MPI_Irecv(td->receiveBuf, td->numTaxa * 2, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, td->mpiRequests + WCR_RCV_ALL);
	MPI_Wait(td->mpiRequests + WCR_RCV_ALL, &mpiStatus);
	DBGPRINT3("W%04d: mpiStatus.tag=%d\n", td->rank, mpiStatus.MPI_TAG);		/*DEBUG*/
	assert(mpiStatus.MPI_TAG == MSG_NEWBOUND);
	NewBoundReceived(td);

	DBGPRINT2("W%04d: Waiting for final (dummy) steal request.\n", td->rank);
	MPI_Irecv(td->receiveBuf, td->numTaxa * 2, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, td->mpiRequests + WCR_RCV_ALL);
	MPI_Wait(td->mpiRequests + WCR_RCV_ALL, &mpiStatus);

	DBGPRINT2("W%04d: Waiting for any outstanding sends to complete.\n", td->rank);
	mpiStatuses = (MPI_Status *) malloc(td->nMpiRequests * sizeof (MPI_Status));
	MPI_Waitall(td->nMpiRequests, td->mpiRequests, mpiStatuses);
	free(mpiStatuses);

	// Iff a maximum number of trees was requested, only send back as many as will be needed.  This entails
	// communicating with prior and subsequent workers.
	if (td->maxTrees) {
		// There is a limit, so we need to communicate with neighbouring workers.
		// Remember td->numTrees counts even trees we didn't store, uMin(td->numTrees, td->maxTrees) is the number of trees we actually stored.
		if (td->rank == 1) {
			nTreesRemaining = td->maxTrees;
		} else {
			MPI_Recv(&nTreesRemaining, 1, MPI_UNSIGNED, td->rank - 1, MSG_NTREESLEFT, MPI_COMM_WORLD, &mpiStatus);
		}

		nTreesToSend = uMin(td->numTrees, nTreesRemaining);		// No need to include td->maxTrees here, since nTreesRemaining is always <= it.
		DBGPRINT6("W%04d: maxTrees=%u; nTreesRemaining=%u; actual number of trees we have stored = %u; nTreesToSend=%u.\n", td->rank, td->maxTrees, nTreesRemaining, uMin(td->numTrees, td->maxTrees), nTreesToSend);

		// Tell the next guy right away
		if (td->rank < td->nWorkers) {
			nTreesRemaining -= nTreesToSend;
			DBGPRINT3("W%04d: nTreesRemaining for next guy = %u.\n", td->rank, nTreesRemaining);
			MPI_Send(&nTreesRemaining, 1, MPI_UNSIGNED, td->rank + 1, MSG_NTREESLEFT, MPI_COMM_WORLD);
		}
	} else {
		nTreesToSend = td->numTrees;
	}

	DBGPRINT3("W%04d: Sending back %u trees.\n", td->rank, nTreesToSend);
	if (td->tempTreeFile) {
		// Send our trees in chunks of TREEBUFSIZE.
		treeBuf = malloc(TREEBUFSIZE);
		rewind(td->tempTreeFile);
		do {
			if (td->maxTrees && nTreesSent == nTreesToSend) {
				// This case occurs in 2 ways:
				// 1. nTreesToSend == 0, but we have some trees stored here.
				// 2. The last byte of the previous block contained the last byte of the final tree to be sent.
				// In either case, we need a zero-length send to signal the end of the data.
				nBytesRead = 0;
			} else {
				nBytesRead = fread(treeBuf, 1, TREEBUFSIZE, td->tempTreeFile);
			}

			if (td->maxTrees) {
				// We need to scan through this buffer to make sure we don't send more trees than we should.
				// I.e. all control of the number of trees sent is done at our end, not the boss's end -- this way
				// we never transfer any unnecessary tree data.
				for (i = 0; i < nBytesRead; ++i) {
					if (treeBuf[i] == '\n') {
						++nTreesSent;
						if (nTreesSent == nTreesToSend) {
							nBytesRead = i + 1;				// "Truncate" the data to be sent
							break;
						}
					}
				}
			}

			DBGPRINT3("W%04d: Sending back %d bytes worth of trees to the boss.\n", td->rank, nBytesRead);
			MPI_Send(treeBuf, nBytesRead, MPI_CHAR, 0, MSG_TREES, MPI_COMM_WORLD);
		} while (nBytesRead == TREEBUFSIZE);
		free(treeBuf);
	} else {
		// No trees to send.
		DBGPRINT2("W%04d: No temporary tree file, so no trees to send back the boss.\n", td->rank);
		MPI_Send(&dummy, 0, MPI_CHAR, 0, MSG_TREES, MPI_COMM_WORLD);
	}

	// We no longer need to (nor should we) send a dummy UB update, now that all requests from workers use MSG_WORKERREQUEST.

	DBGPRINT2("W%04d: Done!\n", td->rank);

	free(td->mpiRequests);
	free(td->stealBuf);
	free(td->receiveBuf);
	td->mpiRequests = NULL;

#ifdef MEASURE
	DumpMeasurements(td);			// This involves waiting our turn.
#endif	// MEASURE
}
#endif	// TARGETMULTI
