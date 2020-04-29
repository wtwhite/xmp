#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#include "common.h"
#include <dbgprint.h>
#include "mpiboss.h"
#include "measure.h"
#ifdef TARGETMULTI
#include "mpi.h"
#endif	// TARGETMULTI

#ifdef TARGETMULTI
enum worker_state {
	HASWORK,				// To the best of our knowledge, this worker has work
	NOWORKTOGIVE,			// This worker responded to a steal request in the negative, but hasn't asked for work yet
	OUTOFWORK,				// This worker has requested more work
	MAXWORKERSTATE
};

//#define MAX_MPI_REQUESTS_AT_ONCE 50

//HACK: Doesn't belong here...
int RandomIntLessThan(int n) {
	return rand() % n;			//TODO: This gives poor randomness.  Replace with something better.
}

// All MPI_Issend()s can probably be changed the following to MPI_Isend() to use buffering if it's available.
void Boss(struct TreeData *td, FILE *outFile) {
	int stateCount[MAXWORKERSTATE];		// The number of workers in each state
	enum worker_state *state = (enum worker_state *) malloc(td->nWorkers * sizeof (enum worker_state));
	MPI_Status *mpiStatuses = (MPI_Status *) malloc(td->nWorkers * sizeof (MPI_Status));
	MPI_Request *workerListenRequests = (MPI_Request *) malloc(td->nWorkers * sizeof (MPI_Request));
	unsigned *workerMsgBufs = (unsigned *) malloc(td->nWorkers * sizeof (unsigned));
	int nRequests;						// Number of workers currently asking or telling us something
	int *indices = (int *) malloc(td->nWorkers * sizeof (int));
	int i, j, k;
	int foundBetterUB;
	MPI_Request *boundBroadcastRequests = (MPI_Request *) malloc(td->nWorkers * sizeof (MPI_Request));
	MPI_Status *boundBroadcastStatuses = (MPI_Status *) malloc(td->nWorkers * sizeof (MPI_Status));
	MPI_Request *stealRequestsAndResponses = (MPI_Request *) malloc(td->nWorkers * 2 * sizeof (MPI_Request));
	MPI_Status stealRequestOrResponseStatus;
	int nStealRequests;
	int victim;
	unsigned (*jobBuf)[2] = (unsigned (*)[2]) malloc(td->nWorkers * (td->numTaxa + 1) * 2 * sizeof (unsigned));
	unsigned *thieves = (unsigned *) malloc(td->nWorkers * sizeof (unsigned));
	MPI_Request *forwardRequests = (MPI_Request *) malloc(td->nWorkers * sizeof (MPI_Request));
	MPI_Status *forwardStatuses = (MPI_Status *) malloc(td->nWorkers * sizeof (MPI_Status));
	int nForwardRequests;
	int nFinished = 0;
	unsigned char *finished = (unsigned char *) malloc(td->nWorkers);
	int count;
	int nJobRequests;

	td->branchAndBoundStartTime = GetTimeNow();

	DBGPRINT1("B: I'm the boss.\n");

	// Initialisation
	stateCount[HASWORK] = 0;
	stateCount[NOWORKTOGIVE] = 0;
	stateCount[OUTOFWORK] = td->nWorkers;
	for (i = 0; i < td->nWorkers; ++i) {
		state[i] = OUTOFWORK;
		finished[i] = 0;
	}

	DBGPRINT1("B: Waiting for the 1st request...\n");
	MPI_Recv(jobBuf, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, MSG_WORKERREQUEST, MPI_COMM_WORLD, mpiStatuses);		// Use jobBuf arbitrarily
	DBGPRINT2("B: Received 1st request from W%d.\n", mpiStatuses[0].MPI_SOURCE);
	MPI_Get_count(mpiStatuses, MPI_UNSIGNED, &count);
	assert(count == 0);		// The first request must be for work (obviously it can't be a new UB...)

	// Give the entire problem to the 1st worker that asks for it!
	// The entire problem looks like 0-length list of edge ranges.
	MPI_Send(jobBuf, 0, MPI_UNSIGNED, mpiStatuses[0].MPI_SOURCE, MSG_NEWJOB, MPI_COMM_WORLD);		// Buffer is unimportant
	state[mpiStatuses[0].MPI_SOURCE - 1] = HASWORK;
	--stateCount[OUTOFWORK];
	++stateCount[HASWORK];

	// Wait for requests from all workers.  Each request has at most 1 unsigned int of data.
	for (i = 0; i < td->nWorkers; ++i) {
		// Tricky: Need to wait separately on each type of message from workers, since if we use MPI_ANY_TAG
		// we may inadvertently grab a MSG_NEWJOB after asking the worker for a job later on.
		// Another way to do this would be to have a single MSG_WORKERREQUEST tag, and supply the type of request
		// (either a new UB or a request for more work) as data.  That would save on MPI resources.  But having
		// separate comm. objects for each tag type is more flexible so let's do that for now ;)
		// Waiting on separate objects does mean that at completion time, each worker needs to send in a
		// dummy message of each type.
		// We now just use a single tag for receiving both new UBs and requests for work from workers.  As
		// explained at the top, this is needed to guarantee that these messages arrive in order, which in turn is
		// needed to avoid deadlock.
		MPI_Irecv(workerMsgBufs + i, 1, MPI_UNSIGNED, i + 1, MSG_WORKERREQUEST, MPI_COMM_WORLD, workerListenRequests + i);
	}

	// Now wait for further work requests to come in.
	//while (stateCount[HASWORK]) {
	while (nFinished < td->nWorkers) {
		// Use MPI_Waitsome() for fairness.
		//int MPI_Waitsome(int incount, MPI_Request *array_of_requests, int *outcount, int *array_of_indices, MPI_Status *array_of_statuses)
		MPI_Waitsome(td->nWorkers, workerListenRequests, &nRequests, indices, mpiStatuses);
		DBGPRINT2("B: %d requests have arrived.\n", nRequests);
		INCMEASURE(bossMainLoopCount);

		// Gather up any new upper bounds, and mark all work-requesting processes as OUTOFWORK
		foundBetterUB = 0;
		nJobRequests = 0;
		for (i = 0; i < nRequests; ++i) {
			// Length == 0 => ASKFORJOB; length == 1 => NEWBOUND.
			MPI_Get_count(mpiStatuses + i, MPI_UNSIGNED, &count);
			DBGPRINT6("B: Request %d is from W%d and has tag %d.  Count = %d.  Index = %d.\n", i, mpiStatuses[i].MPI_SOURCE, mpiStatuses[i].MPI_TAG, count, indices[i]);
			assert(mpiStatuses[i].MPI_TAG == MSG_WORKERREQUEST);		// Seems a bit pointless now...

			// Actually, we don't want to do this here -- we want to wait and see whether this is an
			// unsuccessful steal request.  We only post a new pending receive if it is not.
			//// First, regardless of what happens: post a new pending receive to replace the one just used up.
			//MPI_Irecv(workerMsgBufs + indices[i], 1, MPI_UNSIGNED, indices[i] + 1, mpiStatuses[i].MPI_TAG, MPI_COMM_WORLD, workerListenRequests + i);

			// Is this a UB update?
			if (count == 1) {
				INCMEASURE(newBoundReportedCount);
				if (workerMsgBufs[indices[i]] < td->bound) {
					DBGPRINT3("B: W%d reported a better bound of %u!\n", mpiStatuses[i].MPI_SOURCE, workerMsgBufs[indices[i]]);
					td->bound = workerMsgBufs[indices[i]];
					foundBetterUB = 1;
				} else {
					DBGPRINT3("B: W%d reported a bound of %u, but we're as good or better already.\n", mpiStatuses[i].MPI_SOURCE, workerMsgBufs[indices[i]]);
				}
			}

			// Is this a request for more work?  If so, we need to record it now in a first pass, before making
			// steal attempts, to avoid trying to steal from a worker that is requesting work.
			if (count == 0) {
				DBGPRINT2("B: W%d asked for more work.\n", mpiStatuses[i].MPI_SOURCE);
				// The following assert() is incorrect -- every worker except the one who grabs the very first
				// job initially has status OUTOFWORK.
				//assert(state[indices[i]] != OUTOFWORK);

				--stateCount[state[indices[i]]];
				state[indices[i]] = OUTOFWORK;
				++stateCount[state[indices[i]]];
				++nJobRequests;
				INCMEASURE(jobRequestCount);
			}
		}

		DBGPRINT4("B: After 1st pass, worker states are: %d in HASWORK, %d in NOWORKTOGIVE, %d in OUTOFWORK.\n",
			stateCount[HASWORK],
			stateCount[NOWORKTOGIVE],
			stateCount[OUTOFWORK]
		);

		// If a better UB was found, tell everyone who hasn't finished yet.
		//HACK: Could avoid telling workers who we know have at least as good a bound, but who cares...
		if (foundBetterUB) {
			DBGPRINT2("B: A better bound of %u was reported this round.\n", td->bound);
			INCMEASURE(newBoundBroadcastCount);
			j = 0;
			for (i = 0; i < td->nWorkers; ++i) {
				if (!finished[i]) {
					MPI_CHOSEN_Isend(&td->bound, 1, MPI_UNSIGNED, i + 1, MSG_NEWBOUND, MPI_COMM_WORLD, boundBroadcastRequests + j++);
				}
			}

			assert(j == td->nWorkers - nFinished);		// This caught a bug in the Promela model!

			// The following MPI_Waitall() implies that unfinished workers must be able to receive a new bound at any
			// time -- otherwise we'll deadlock.
			DBGPRINT2("B: Waiting for new bound to be sent to all %d unfinished workers.\n", td->nWorkers - nFinished);
			MPI_Waitall(td->nWorkers - nFinished, boundBroadcastRequests, boundBroadcastStatuses);
		}

		// Handle all other requests, which should be for more work.
		// Because we already updated the states of work-requesting workers in a previous pass,
		// we will never make a steal request to a worker which is itself making a steal request this round.
		// It also simplifies thinking about the logic to know that the set of all thieves is disjoint from
		// the set of all victims.
		nStealRequests = 0;
		for (i = 0; i < nRequests; ++i) {
			MPI_Get_count(mpiStatuses + i, MPI_UNSIGNED, &count);
			if (count == 0) {
				assert(state[indices[i]] == OUTOFWORK);			// First pass establishes this.

				if (!stateCount[HASWORK]) {
					// We're done!
					break;
				}

				// Pick a random worker and try to steal a job from them.
				// Although we only pick from among workers with status HASWORK, it could well be that
				// the worker chosen was just about to run out of work and now has, so the steal could fail.
				j = RandomIntLessThan(stateCount[HASWORK]) + 1;
				for (victim = 0; j; ++victim) {
					if (state[victim] == HASWORK) {
						--j;
					}
				}
				--victim;

				// Workers that we're attempting to steal from may try to send ASKFORJOB requests or NEWBOUNDs
				// of their own, but we're only interested in the response to the STEALJOB attempt.  Because
				// MPI doesn't allow receiving of partial wildcards, that means we need to ask for only MSG_NEWJOB
				// responses, and encode "can't satisfy the steal" in this somehow, rather than use a different
				// message type to indicate this.
				// Tricky: jobBuf elements have type unsigned [2], so pointer arithmetic implicitly multiplies by 2.
				// We use stealRequestsAndResponses to hold statuses for both the steal request sent to the victim and
				// the response.  By always putting the latter first, we can quickly check whether a completed
				// communication refers to a send or receive with a simple index comparison, and in the important case
				// of a receive, thieves[index] is the relevant thief.
				DBGPRINT3("B: Attempting to steal from W%d to give to W%d.\n", victim + 1, indices[i] + 1);
				MPI_CHOSEN_Isend(jobBuf, 0, MPI_UNSIGNED, victim + 1, MSG_STEALJOB, MPI_COMM_WORLD, stealRequestsAndResponses + nJobRequests + nStealRequests);		// The buffer isn't used
				MPI_Irecv(jobBuf + indices[i] * td->numTaxa, td->numTaxa * 2, MPI_UNSIGNED, victim + 1, MSG_NEWJOB, MPI_COMM_WORLD, stealRequestsAndResponses + nStealRequests);
				thieves[nStealRequests] = indices[i];
				++nStealRequests;
				INCMEASURE(stealCount);
			}
		}

		DBGPRINT2("B: Initial number of steal requests: %d\n", nStealRequests);

		// Repeatedly wait for the outcome of steal requests, and make new ones for any that failed.
		while (nStealRequests) {
			// At this point, there should be nStealRequests outgoing MSG_STEALJOBs pending, and the same number of
			// MSG_NEWJOB receives pending.

			// See whether our steal requests worked.

			// First pass: record all failed steals.  This avoids wastefully requesting multiple further steals
			// from workers that have already run out of work.
			// We also perform forwarding here for successful steals.  This is important as the successful thieves go back
			// to having HASWORK state, meaning they themselves can be stolen from in the next round of re-attempts.  If
			// we didn't do this, workers may starve.
			// This loop waits for both the sends and receives, hence the "* 2".
			nForwardRequests = 0;
			for (i = 0; i < nStealRequests * 2; ++i) {
				MPI_Waitany(nJobRequests * 2, stealRequestsAndResponses, &k, &stealRequestOrResponseStatus);

				if (k < nJobRequests) {		// Otherwise, a steal request has completed: big deal.
					// A steal response has been received.  In this case, thieves[k] will be the relevant thief.
					MPI_Get_count(&stealRequestOrResponseStatus, MPI_UNSIGNED, &count);
					if (count == 1) {
						// Steal request can't be satisfied.  Try another one.
						DBGPRINT3("B: Steal from W%d to give to W%d FAILED.  Will re-attempt if some workers have work left.\n", stealRequestOrResponseStatus.MPI_SOURCE, thieves[k] + 1);
						INCMEASURE(stealFailCount);

						//HACK: Some code near-duplication follows.
						// Although we only asked workers with status HASWORK, it might be that > 1 steal request went to
						// the same worker, and one of these requests could not be satisfied.  So it's possible that
						// the victim's state is already NOWORKTOGIVE.
						assert(state[stealRequestOrResponseStatus.MPI_SOURCE - 1] == HASWORK || state[stealRequestOrResponseStatus.MPI_SOURCE - 1] == NOWORKTOGIVE);
						--stateCount[state[stealRequestOrResponseStatus.MPI_SOURCE - 1]];
						state[stealRequestOrResponseStatus.MPI_SOURCE - 1] = NOWORKTOGIVE;
						++stateCount[state[stealRequestOrResponseStatus.MPI_SOURCE - 1]];
					} else {
						// The steal request can be satisfied.  Forward the job to the thief.
						DBGPRINT4("B: Steal from W%d to give to W%d SUCCEEDED (job tree size = %d).  Forwarding.\n", stealRequestOrResponseStatus.MPI_SOURCE, thieves[k] + 1, count / 2);
						MPI_CHOSEN_Isend(jobBuf + thieves[k] * td->numTaxa, count, MPI_UNSIGNED, thieves[k] + 1, MSG_NEWJOB, MPI_COMM_WORLD, forwardRequests + nForwardRequests++);

						--stateCount[state[thieves[k]]];
						state[thieves[k]] = HASWORK;
						++stateCount[state[thieves[k]]];
					}
				}
			}

			// Wait for all stolen jobs to be forwarded to successful thieves.  This must be done before allowing
			// attempts to steal from these successful ex-thieves, since otherwise such a steal request may fail
			// due to the forwarding operation not having completed yet.  The victim would then be marked as NOWORKTOGIVE,
			// meaning that the job it had just stolen would effectively "disappear off the radar".  This would not
			// make the algorithm incorrect, but risks overloading some workers while starving others.

			// The following MPI_Waitall() implies that all workers who asked for more work must be able to receive
			// those jobs without blocking.
			DBGPRINT2("B: Waiting for %d forwarded jobs to be sent to workers.\n", nForwardRequests);
			MPI_Waitall(nForwardRequests, forwardRequests, forwardStatuses);

			if (!stateCount[HASWORK]) {
				DBGPRINT1("B: No more workers with work to give!\n");
				break;
			}

			// Second pass: attempt re-stealing of all failed steals.  Again, the previous pass
			// ensures that we have that thieves are disjoint from victims.  We also know that there is at least
			// 1 worker in state HASWORK who it is worth attempting to steal from, since we exit the loop above
			// if there are none left, and state counts are not changed by the loop below.
			// Note that all steals could fail even though computation is not complete -- it may be that some workers
			// with state HASWORK were not chosen as victims in the previous round of steal attempts.
			k = 0;					// k is the location to write the next re-steal request
			for (i = 0; i < nStealRequests; ++i) {
				// "Can't satisfy the steal" is encoded as a length-1 response.  Normally, all jobs have even length
				// as they are an array of pairs.
				if (state[thieves[i]] == OUTOFWORK) {
					// Steal request couldn't be satisfied.  Try another one.

					// Pick a random worker and try to steal a job from them.
					// Although we only pick from among workers with status HASWORK, it could well be that
					// the worker chosen was just about to run out of work and now has, so the steal could fail.
					j = RandomIntLessThan(stateCount[HASWORK]) + 1;
					for (victim = 0; j; ++victim) {
						if (state[victim] == HASWORK) {
							--j;
						}
					}
					--victim;

					// Tricky: in order to be able to reuse the job buffers for forwarding the jobs to the workers
					// if and when they are successfully stolen, we need to always use the same job buffer for each
					// thief.  That means using the thieves[i]-th buffer, not the ith -- otherwise later iters of
					// this loop can overwrite received jobs that are waiting to be forwarded.
					DBGPRINT3("B: Re-attempting to steal from W%d to give to W%d.\n", victim + 1, thieves[i] + 1);
					MPI_CHOSEN_Isend(jobBuf, 0, MPI_UNSIGNED, victim + 1, MSG_STEALJOB, MPI_COMM_WORLD, stealRequestsAndResponses + nJobRequests + k);		// The buffer isn't used
					MPI_Irecv(jobBuf + thieves[i] * td->numTaxa, td->numTaxa * 2, MPI_UNSIGNED, victim + 1, MSG_NEWJOB, MPI_COMM_WORLD, stealRequestsAndResponses + k);
					thieves[k] = thieves[i];
					++k;
					INCMEASURE(stealCount);
				}
			}

			nStealRequests = k;
		}

		// For each thief, the possible state transitions during the above "while (nStealRequests) { ... }" loop are:
		//
		// OUTOFWORK (only possible if all steal attempts failed and there are no more HASWORK workers to steal from;
		//    in this case, the thief will still be waiting for a response to its ASKFORJOB request)
		// OUTOFWORK->HASWORK (a steal attempt succeeded; job forwarding has been completed; possibly other thieves later
		//    stole successfully from us)
		// OUTOFWORK->HASWORK->NOWORKTOGIVE (a steal attempt succeeded; job forwarding has been completed;
		//    >= 1 other thieves later stole unsuccessfully from us)
		//
		// In all cases, stateCount[HASWORK] == 0 iff some thief could not be satisfied.  This is the termination condition.
		// An equivalent termination condition is if any of the workers whose requests were processed in this round
		// now has state OUTOFWORK.  (Other workers may have OUTOFWORK status -- e.g. near the start when not
		// all initial requests have been processed.)
		//TODO: find a way to assert() this.  The assert() below is wrong because it doesn't account for other workers
		// who didn't have requests processed during this round.
//		assert((stateCount[HASWORK] == 0) == (stateCount[OUTOFWORK] > 0)) -------  THIS IS WRONG

		DBGPRINT4("B: After 2nd pass, worker states are: %d in HASWORK, %d in NOWORKTOGIVE, %d in OUTOFWORK.\n",
			stateCount[HASWORK],
			stateCount[NOWORKTOGIVE],
			stateCount[OUTOFWORK]
		);

		// 3rd and final pass, telling unsuccessful thieves that it's all over.
		// We also reinstate receives for any other requests that we just "used up".
		nForwardRequests = 0;	// Re-use to count the number of "finished" messages sent
		for (i = 0; i < nRequests; ++i) {
			// Determine which worker we are talking about.  Remember we have two
			// requests pending for each.
			j = indices[i] + 1;		// Don't need "%" now that we only have NWORKERS requests.  TODO: Get rid of this completely.
			assert(j == mpiStatuses[i].MPI_SOURCE);		// Added from the Promela model
			MPI_Get_count(mpiStatuses + i, MPI_UNSIGNED, &count);
			DBGPRINT7("B: 3rd pass: looking at request %d, index %d, from W%d (in state %d), of tag %d, length %d.\n", i, indices[i], j, state[j - 1], mpiStatuses[i].MPI_TAG, count);
			if (count == 0 && state[j - 1] == OUTOFWORK) {
				// We were unable to perform a steal for this thief.  Tell them we're finished.
				// This will be the last message sent to this worker during this loop.  This is an
				// important property to guarantee, since otherwise the worker doesn't know whether to
				// wait for more NEWBOUND or STEALJOB requests.
				assert(stateCount[HASWORK] == 0);
				DBGPRINT2("B: Unable to perform a steal for W%d, and there are no more potential victims.  Telling thief we're finished.\n", j); 
				MPI_CHOSEN_Isend(jobBuf + j - 1, 1, MPI_UNSIGNED, mpiStatuses[i].MPI_SOURCE, MSG_NEWJOB, MPI_COMM_WORLD, forwardRequests + nForwardRequests++);
				finished[mpiStatuses[i].MPI_SOURCE - 1] = 1;
			} else {
				// Reinstate this pending receive.
				MPI_Irecv(workerMsgBufs + indices[i], 1, MPI_UNSIGNED, j, MSG_WORKERREQUEST, MPI_COMM_WORLD, workerListenRequests + indices[i]);
			}
		}

		MPI_Waitall(nForwardRequests, forwardRequests, forwardStatuses);
		nFinished += nForwardRequests;
	}

	DBGPRINT1("B: Main loop has finished.\n");

	// All workers have been sent a "finished" message.  They are still all waiting on NEWBOUND and STEALJOB
	// receives, however.
	
	// Send the final bound to each worker (this may well have changed since the time that worker was sent "finished").
	// Also send a dummy STEALJOB request.  This is just to consume each worker's pending receive on this message type, meaning
	// we can avoid having to call MPI_Cancel().
	for (i = 0; i < td->nWorkers; ++i) {
		// I think we could just use MPI_Isend() here.
		//HACK: We should put both requests in a single array so we can use a single MPI_Waitall() call below.
		MPI_CHOSEN_Isend(&td->bound, 1, MPI_UNSIGNED, i + 1, MSG_NEWBOUND, MPI_COMM_WORLD, boundBroadcastRequests + i);
		MPI_CHOSEN_Isend(jobBuf, 0, MPI_UNSIGNED, i + 1, MSG_STEALJOB, MPI_COMM_WORLD, stealRequestsAndResponses + i);		// Buffer is unimportant
	}

	DBGPRINT2("B: Waiting for all final UB updates to be sent.  Final UB = %d.\n", td->bound);
	MPI_Waitall(td->nWorkers, boundBroadcastRequests, boundBroadcastStatuses);

	DBGPRINT1("B: Waiting for all dummy steal requests to be sent.\n");
	for (i = 0; i < td->nWorkers; ++i) {
		MPI_Waitany(td->nWorkers, stealRequestsAndResponses, &k, &stealRequestOrResponseStatus);		// We ignore k
	}

	// We want to exclude the time it takes to send the actual tree data via the network.
	td->outputResultsStartTime = GetTimeNow();

	OutputResults(td, outFile);

	// Now that we have combined the ASKFORJOB and NEWBOUND receives into a single receive, it's impossible to receive
	// a new UB after a worker has been sent a "show's over" message, meaning it no longer makes sense to wait for receipt
	// of a "dummy" UB update.

	DBGPRINT1("B: Done!\n");
	td->endTime = GetTimeNow();

#ifdef MEASURE
	DumpMeasurements(td);		// This involves waiting our turn.  We'll write our stuff last in fact.
#endif	// MEASURE
	ReportTiming(td);

	free(finished);
	free(state);
	free(mpiStatuses);
	free(workerListenRequests);
	free(workerMsgBufs);
	free(indices);
	free(boundBroadcastRequests);
	free(boundBroadcastStatuses);
	free(stealRequestsAndResponses);
	free(jobBuf);
	free(thieves);
	free(forwardRequests);
	free(forwardStatuses);
}
#endif	 // TARGETMULTI
