#include "common.h"
#include <dbgprint.h>

int fitchCount = 0;
int fitch2SeqsCount = 0;
int newNodeCount = 0;
int insertNodeCount = 0;
int branchAndBoundCount = 0;
int branchAndBoundScanEdgeCount = 0;
int branchAndBoundScanEdgeCountByDepth[MAXTAXA] = { 0 };
int minBoundSavingsCount = 0;
int colBoundSavingsCount = 0;
int colBoundSavingsCountByDepth[MAXTAXA] = { 0 };	//HACK
int treeList_AddTreeCount = 0;
int treeList_DestroyCount = 0;
int saveTreeListCount = 0;
int copyTreeCount = 0;
int destroyTreeCount = 0;
int zeroAdditionalScoreCount = 0;
int bestEqualCounts[MAXTAXA] = { 0 };
int fitchBasesCount = 0;
int scoreFitchCount = 0;
int postorderFitchCallsFitchBasesCount = 0;
int preorderFitchCallsFitchBasesCount = 0;
int updateTreePointerCopiesCount = 0;
#ifdef DETECTREDUNDANTFITCH
int redundantFitchCountByDepth[MAXTAXA + 5] = { 0 };		//HACK
#endif	// DETECTREDUNDANTFITCH
int incFitchRootCount = 0;
int fullTreesConsideredCount = 0;
#ifdef TARGETMULTI
int bossMainLoopCount = 0;
int jobRequestCount = 0;
int stealCount = 0;
int stealFailCount = 0;
int newBoundReportedCount = 0;
int newBoundBroadcastCount = 0;
#endif	// TARGETMULTI

// This seems to be obsolete -- these variables are not updated anywhere.
#ifdef FITCHEXTMEASURE
int unnecessaryFitchCount = 0;
int unnecessaryFitchCountHalf1 = 0;
int unnecessaryFitchCountHalf2 = 0;
#endif	// FITCHEXTMEASURE

void DumpMeasurements(struct TreeData *td) {
	unsigned i;
	FILE *measureFile = stderr;		//HACK: should really be in struct TreeData, but I can't be bothered...
	
#ifdef TARGETMULTI
	MPI_Status mpiStatus;

	// Write measurement output in the order W1, W2, ..., B.  (Put the boss at the end because that's most
	// likely to contain the most useful info, making it easy to find with tail.)
	if (!td->rank) {
		// Kick things off.  Empirically, flushing doesn't guarantee that this output won't be overlapped
		// with other output, unfortunately.
		fflush(measureFile);
		MPI_Ssend(NULL, 0, MPI_CHAR, 1, 0, MPI_COMM_WORLD);		// Don't care about the buffer or tag
	}

	// Wait for the previous guy to finish writing his output
	i = td->rank ? td->rank - 1 : td->nWorkers;		// The boss waits for the last worker
	MPI_Recv(NULL, 0, MPI_CHAR, i, 0, MPI_COMM_WORLD, &mpiStatus);		// Don't care about the buffer or tag

	if (td->rank) {
		fprintf(measureFile, "W%04d measurements:\n", td->rank);
	} else {
		fprintf(measureFile, "B measurements:\n");
	}
#endif	// TARGETMULTI

	fprintf(measureFile, "minBoundSavingsCount =                 %10d\n", minBoundSavingsCount);
	fprintf(measureFile, "BranchAndBound() executed              %10d times\n", branchAndBoundCount);
	fprintf(measureFile, "BAndBScanEdge() executed               %10d times\n", branchAndBoundScanEdgeCount);
	
	for (i = 0; i < td->numTaxa; ++i) {
		fprintf(measureFile, "bAndBScanEdgeCountByDepth[%2d] =        %10d times\n", i, branchAndBoundScanEdgeCountByDepth[i]);
	}
	fprintf(measureFile, "zeroAdditionalScoreCount =             %10d times\n", zeroAdditionalScoreCount);
	
	fprintf(measureFile, "InsertNode() executed                  %10d times\n", insertNodeCount);
	fprintf(measureFile, "Fitch() executed                       %10d times\n", fitchCount);
	fprintf(measureFile, "Fitch2Seqs() executed                  %10d times\n", fitch2SeqsCount);
	fprintf(measureFile, "treeList_AddTree executed              %10d times\n", treeList_AddTreeCount);
	fprintf(measureFile, "treeList_Destroy executed              %10d times\n", treeList_DestroyCount);
	fprintf(measureFile, "saveTreeList executed                  %10d times\n", saveTreeListCount);
	fprintf(measureFile, "copyTree executed                      %10d times\n", copyTreeCount);
	fprintf(measureFile, "destroyTree executed                   %10d times\n", destroyTreeCount);
#ifdef FITCHEXTMEASURE
	fprintf(measureFile, "Fitch() executed unnecessarily         %10d times\n", unnecessaryFitchCount);
	fprintf(measureFile, "Fitch() executed unnecessarily (0-0.5) %10d times\n", unnecessaryFitchCountHalf1);
	fprintf(measureFile, "Fitch() executed unnecessarily (0.5-1) %10d times\n", unnecessaryFitchCountHalf2);
#endif	// FITCHEXTMEASURE
	fprintf(measureFile, "FitchBases() executed                  %10d times\n", fitchBasesCount);
	fprintf(measureFile, "ScoreFitch() executed                  %10d times\n", scoreFitchCount);
	fprintf(measureFile, "FitchBases() called by Postorder...()  %10d times\n", postorderFitchCallsFitchBasesCount);
	fprintf(measureFile, "FitchBases() called by Preorder...()   %10d times\n", preorderFitchCallsFitchBasesCount);
	fprintf(measureFile, "Fast updates in UpdateTree()           %10d times\n", updateTreePointerCopiesCount);
	fprintf(measureFile, "Number of trees in search space:       %10u * number of subtrees for each tree on %u taxa\n", td->proportionComplete, COUNT_UNIT);
	fprintf(measureFile, "Number of trees actually examined:     %10u * number of subtrees for each tree on %u taxa\n", td->treesActuallyExamined, COUNT_UNIT);
	fprintf(measureFile, "Number of full trees considered:       %10u\n", fullTreesConsideredCount);
#ifdef TARGETMULTI
	if (!td->rank) {
		fprintf(measureFile, "Boss main loop iterations:             %10u\n", bossMainLoopCount);
		fprintf(measureFile, "Jobs requested:                        %10u\n", jobRequestCount);
		fprintf(measureFile, "Steal requests made:                   %10u\n", stealCount);
		fprintf(measureFile, "Steal requests failed:                 %10u\n", stealFailCount);
		fprintf(measureFile, "Better LBs reported by workers:        %10u\n", newBoundReportedCount);
		fprintf(measureFile, "Actually-better LBs broadcast back:    %10u\n", newBoundBroadcastCount);
	} else {
		// Tell the next guy to write his output
		fflush(measureFile);
		i = td->rank == td->nWorkers ? 0 : td->rank + 1;		// The last worker tells the boss
		MPI_Ssend(NULL, 0, MPI_CHAR, i, 0, MPI_COMM_WORLD);		// Don't care about the buffer or tag
	}
#endif	// TARGETMULTI
}

void ShowDepthsForPeriod(struct TreeData *td) {
	unsigned i;
	unsigned long total = 0;
	unsigned highestBestEqual;

	for (i = 0; i < td->numTaxa; ++i) {
		total += branchAndBoundScanEdgeCountByDepth[i];
	}
	
	for (i = 0; i < td->numTaxa; ++i) {
		DBGPRINT4("bAndBScanEdgeCountByDepth[%2d] =        %10d times (%.2f%%)\n", i, branchAndBoundScanEdgeCountByDepth[i], ((double) branchAndBoundScanEdgeCountByDepth[i]) / total * 100);
		branchAndBoundScanEdgeCountByDepth[i] = 0;		//HACK
	}
	
#ifdef DETECTREDUNDANTFITCH
	for (i = 0; i < td->numTaxa; ++i) {
		DBGPRINT3("redundantFitchCountByDepth[%2d] =       %10d times\n", i, redundantFitchCountByDepth[i]);
		redundantFitchCountByDepth[i] = 0;				//HACK
	}
#endif	// DETECTREDUNDANTFITCH
	
	DBGPRINT2("zeroAdditionalScoreCount =             %10d times\n", zeroAdditionalScoreCount);
	zeroAdditionalScoreCount = 0;
	
//	for (i = 0; i < td->numTaxa; ++i) {
//		DBGPRINT3("bestEqualCounts[%2d] =                  %10d times\n", i, bestEqualCounts[i]);
//		branchAndBoundScanEdgeCountByDepth[i] = 0;		//HACK
//	}
	for (i = td->numTaxa - 1; i > 0; --i) {
		if (bestEqualCounts[i]) {
			highestBestEqual = i;
			break;
		}
	}
	
	DBGPRINT3("highest bestEqualCount = %d (happened %d times)\n", highestBestEqual, bestEqualCounts[highestBestEqual]);
	memset(bestEqualCounts, 0, sizeof bestEqualCounts);
}
