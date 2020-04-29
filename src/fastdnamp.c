#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
//#include <mem.h>
#include <string.h>
#include <assert.h>
#include "common.h"
#include "seq.h"
#include "yanbader.h"
#include "measure.h"
#include "timer.h"
#include "elapsedtime.h"
#include <dbgprint.h>		// My DBGPRINTn() macros
#ifdef TARGETMULTI
#include "mpi.h"
#include "mpiboss.h"
#include "mpiworker.h"
#endif	// TARGETMULTI

int main(int argc, char **argv) {
	struct tree *root;
	char *inFileName, *outFileName;
	FILE *inFile, *outFile = stdout;
	struct CharData *cd = &gCD;
	struct TreeData *td = &gTD;
#ifdef TARGETMULTI
	int rank, nProcs;
#endif	// TARGETMULTI

	td->startTime = GetTimeNow();

#ifdef TARGETMULTI
	// MPI initialisation
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
#endif	// TARGETMULTI

	//dbgfile = stdout;		// Just use stdout for output for SP version
	dbgfile = stderr;		// stderr is easier for debugging MP version -- assertion failures etc. will be in sequence
	InitTimer();
	
#ifdef DEBUG
#ifdef TARGETMULTI
	// Only the master outputs settings.
	if (!rank) {
#endif	// TARGETMULTI
	DebugPrintSwitches();
#ifdef TARGETMULTI
	}
#endif	// TARGETMULTI
#endif	// DEBUG
	
	SetupDefaults(td);
	
	if (!ProcessCmdLineArgs(argc, argv, td, &inFileName, &inFile, &outFileName, &outFile)) {
		return 0;
	}
	
#ifdef TARGETMULTI
	// Only the master outputs settings.
	if (!rank) {
#endif	// TARGETMULTI
	ReportSettings(td, inFileName, inFile, outFileName, outFile);
#ifdef TARGETMULTI
	}
#endif	// TARGETMULTI
	
	// Read a PHYLIP-format sequence alignment into a character matrix
	ReadPhylipAlignment(inFile, cd);
#ifdef TARGETMULTI
	// Only the master outputs messages.
	if (!rank) {
#endif	// TARGETMULTI
	if (!(td->options & OPT_QUIET)) {
		fprintf(stderr, "Nucleotide data for %u taxa x %u sites read in.\n", cd->numTaxa, cd->seqLen);
	}
#ifdef TARGETMULTI
	}
#endif	// TARGETMULTI
	
	if (inFile != stdin) {
		fclose(inFile);
	}
	
	// Set up the data structure that contains all information relevant to
	// building a tree on this character matrix
	InitTreeData(cd, td);
	
	if (!(td->options & OPT_ONLYCOMPUTELOWERBOUND)) {
		root = ConstructBaseTree(td, 1, 2, 3);
		
#ifdef DEBUG
		branchAndBoundMinNumTaxa = branchAndBoundMaxNumTaxa = 3;
#endif	// DEBUG
		
		// Perform search
#ifdef TARGETMULTI
		// Only the master outputs messages.
		if (!rank) {
#endif	// TARGETMULTI
		if (!(td->options & OPT_QUIET)) {
			fprintf(stderr, "Beginning branch and bound search...\n");
		}
#ifdef TARGETMULTI
		}
#endif	// TARGETMULTI

#ifdef TARGETMULTI
		td->rank = rank;
		td->nWorkers = nProcs - 1;
		if (!rank) {
			Boss(td, outFile);
		} else {
			Worker(root, td);
		}
#else	// not TARGETMULTI
		if (td->updateIntervalSecs) {
			td->timerHandle = StartTimer(&td->updateProgressNow, td->updateIntervalSecs, 1);
			DBGPRINT2("Just started timer for the first time.  td->timerHandle=%u.\n", td->timerHandle);
		}
		
		td->branchAndBoundStartTime = GetTimeNow();
		
		BranchAndBound(root, 3, td);
		
		EMMS;	// Needed for MMXASMFITCH
		
		td->outputResultsStartTime = GetTimeNow();
		
		// Clean up
		if (td->updateIntervalSecs) {
			StopTimer(td->timerHandle);
			DBGPRINT2("Just stopped timer for the last time.  td->timerHandle=%u.\n", td->timerHandle);
		}

		OutputResults(td, outFile);
		td->endTime = GetTimeNow();

#ifdef MEASURE
		DumpMeasurements(td);
#endif	// MEASURE
		ReportTiming(td);
#endif	// not TARGETMULTI
	}

	FinaliseTimer();
	DestroyCharData(cd);
	DestroyTreeData(td);
	if (dbgfile != stdout && dbgfile != stderr) {
		fclose(dbgfile);
	}
	
#ifdef TARGETMULTI
	MPI_Finalize();
#endif	// TARGETMULTI

	return 0;
}
