#ifndef __MEASURE_H
#define	__MEASURE_H
#include "common.h"

// INCMEASURE(x) becomes a no-op if MEASURE is not defined.
#ifdef MEASURE
#define INCMEASURE(x) ++(x)
#else
#define INCMEASURE(x)
#endif

#ifdef MEASURE
extern int fitchCount;
extern int fitch2SeqsCount;
extern int newNodeCount;
extern int insertNodeCount;
extern int branchAndBoundCount;
extern int branchAndBoundScanEdgeCount;
extern int branchAndBoundScanEdgeCountByDepth[MAXTAXA];
extern int minBoundSavingsCount;
extern int colBoundSavingsCount;
extern int colBoundSavingsCountByDepth[MAXTAXA];	//HACK
extern int treeList_AddTreeCount;
extern int treeList_DestroyCount;
extern int saveTreeListCount;
extern int copyTreeCount;
extern int destroyTreeCount;
extern int zeroAdditionalScoreCount;
extern int bestEqualCounts[MAXTAXA];
extern int fitchBasesCount;
extern int scoreFitchCount;
extern int postorderFitchCallsFitchBasesCount;
extern int preorderFitchCallsFitchBasesCount;
extern int updateTreePointerCopiesCount;
#ifdef DETECTREDUNDANTFITCH
extern int redundantFitchCountByDepth[MAXTAXA + 5];		//HACK
#endif	// DETECTREDUNDANTFITCH
extern int incFitchRootCount;
extern int fullTreesConsideredCount;
#ifdef TARGETMULTI
extern int bossMainLoopCount;
extern int jobRequestCount;
extern int stealCount;
extern int stealFailCount;
extern int newBoundReportedCount;
extern int newBoundBroadcastCount;
#endif	// TARGETMULTI

// This seems to be obsolete -- these variables are not updated anywhere.
#ifdef FITCHEXTMEASURE
extern int unnecessaryFitchCount;
extern int unnecessaryFitchCountHalf1;
extern int unnecessaryFitchCountHalf2;
#endif	// FITCHEXTMEASURE

void DumpMeasurements(struct TreeData *td);
void ShowDepthsForPeriod(struct TreeData *td);
#endif	// MEASURE
#endif	// #include guard
