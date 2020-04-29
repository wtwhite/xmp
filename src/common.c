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
#include <procargs.h>		// Part of my own library
#include <dbgprint.h>		// My DBGPRINTn() macros
#include "bound.h"
#include "ordertaxa.h"
#include "maketempfile.h"
#include "measure.h"
#include "alignedalloc.h"
#include "elapsedtime.h"
#include "greedytree.h"
#include "tbr.h"

//HACK: Linux has a random() function which takes no args, unlike Windows' version which takes one.  This isn't the right place to do this but hey...
// It turns out that MS VC++ doesn't have this function either!  Bad Borland.
//#ifdef __linux__
#ifndef __BORLANDC__
#define random(x) (rand() % (x))
#endif	// not __BORLANDC__

void FillTaxonSequence(struct TreeData *td);
//unsigned PadToBlockSize(unsigned char *fromBuf, unsigned *fromWeights, unsigned **pToBuf, unsigned **pToWeights, unsigned seqLen, unsigned numTaxa, unsigned effNumTaxa);
unsigned PadToBlockSize(unsigned char *fromBuf, unsigned *fromWeights, unsigned char **pToBuf, unsigned **pToWeights, unsigned seqLen, unsigned numTaxa, unsigned effNumTaxa);

//HACK: limits.h is supposed to #define PATH_MAX, but it doesn't on Windows and according to the Linux man pages,
// systems that do define it may give it an unreasonably large value or -1 to indicate "no limit".  So instead
// of dynamically resizing until things fit, just pick this random size as something that should be large enough
// without breaking the bank.
#define MY_PATH_MAX 4096

// Global Variables
struct CharData gCD;
struct TreeData gTD;
FILE *dbgfile;
unsigned _ColumnCompareWidth;			// Used by several qsort() callbacks in this and other modules.

#ifdef DEBUG
unsigned branchAndBoundMinNumTaxa = 0;
unsigned branchAndBoundMaxNumTaxa = 0;
#endif	// DEBUG

void SetupDefaults(struct TreeData *td) {
	td->options =
		OPT_NEXUSFMT |			// By default, look for all most parsimonious trees and output in NEXUS format
		OPT_DISPLAYNEWBOUNDS |
		OPT_DISPLAYNEWTREES |
		OPT_USEDISTBASESBOUNDREST |
		OPT_IMPROVEINITUPPERBOUND |
		OPT_TBR;
	td->bound = INT_MAX;
	td->maxTrees = INT_MAX;
	td->startCol = 0;
	td->endCol = INT_MAX;
	td->restBoundInputFname = NULL;
	td->partBoundStrategy = PBS_CLIPPEDANDADDITIVE;
	td->partBoundMaxItersNoImprovement = 2;
	td->updateIntervalSecs = 2;
	td->rngSeed = 0;			// 0 means "Use a 'random' seed"
	td->rand1 = 0;				// See ProcessCmdLineArgs() for an explanation.
	td->rand2 = 0;
	td->rand3 = 0;
}

int ProcessCmdLineArgs(int argc, char **argv, struct TreeData *td, char **pInFileName, FILE **pInFile, char **pOutFileName, FILE **pOutFile) {
	int i;
	size_t j;
	unsigned do_generateSitePairCostsFile = 0;			//HACK: This is just used for regenerating the sitepaircosts.c file.
	
	if (argc < 2) {
		PrintUsage();
		return 0;
	}
	
	*pInFile = NULL;
	*pOutFile = stdout;
	
	// Read command-line arguments
	for (i = 1; i < argc; ++i) {
		if (argv[i][0] == '-') {
			// Process switches
			switch (argv[i][1]) {
			case 'b':		// -b <initial bound>
				td->bound = (unsigned) atoi(argv[++i]);
				if (td->bound == 0) td->bound = INT_MAX;	// Interpret zero as "no bound"
				break;
			
			case 'o':		// -o <output file>
				if ((*pOutFile = fopen(argv[++i], "w")) == NULL) {
					fprintf(stderr, "File '%s' could not be opened for output, aborting.\n", argv[i]);
					exit(1);
				} else {
					*pOutFileName = argv[i];
				}
				break;
			
			case 'm':		// -m <no of trees to keep>
				td->maxTrees = (unsigned) atoi(argv[++i]);
				if (td->maxTrees == 0) td->maxTrees = INT_MAX;	// Interpret zero as "no limit"
				break;
			
			case 'r':		// Output trees in raw (Newick) format, rather than NEXUS format
				td->options &= ~OPT_NEXUSFMT;
				break;
			
			case 'd':		// Discard sites containing gaps or ambiguous characters
				td->options |= OPT_DISCARDMISSAMBIG;
				break;
			
			case 'q':		// Quiet mode
				td->options |= OPT_QUIET;
				td->options &= ~(OPT_DISPLAYNEWBOUNDS | OPT_DISPLAYNEWTREES);
				td->updateIntervalSecs = 0;
				break;
			
			case 'B':		//HACK: should really make this an extended option.  Added 12/6/2004.
				// 28/2/2005: Now make this a 2-char switch to choose all bounding options.
				td->options &= ~(
					OPT_USELOADRESTBOUND |
					OPT_USEDISTBASESBOUNDREST |
					OPT_USEINCOMPATBOUNDREST |
					OPT_USEMSTBOUNDREST |
					OPT_USEKICKASSBOUNDREST |
					OPT_USEKA2BOUNDREST |
					OPT_USEGREEDYHITTINGSETBOUND |
					OPT_USEPARTBOUND
				);
				
				for (j = 2; j < strlen(argv[i]); ++j) {
					switch (argv[i][j]) {
						case 'd': td->options |= OPT_USEDISTBASESBOUNDREST; break;
						case 'i': td->options |= OPT_USEINCOMPATBOUNDREST; break;
						case 'm': td->options |= OPT_USEMSTBOUNDREST; break;
						case 'k': td->options |= OPT_USEKICKASSBOUNDREST; break;
						case 'K': td->options |= OPT_USEKA2BOUNDREST; break;
						case 'h': td->options |= OPT_USEGREEDYHITTINGSETBOUND; break;
						case 'p': td->options |= OPT_USEPARTBOUND; break;
						case 'C':
							if (td->options & OPT_USEPARTBOUND) {
								td->partBoundStrategy = PBS_CLIPPED;
							} else {
								fprintf(stderr, "To specify the 'C' option for clipped scores only, you must have previously specified -Bp.\n");
								exit(1);
							}
							break;

						case 'L':
							if (td->options & OPT_USEPARTBOUND) {
								td->partBoundStrategy = PBS_LB;
							} else {
								fprintf(stderr, "To specify the 'L' option for lower bound scores only, you must have previously specified -Bp.\n");
								exit(1);
							}
							break;

						case 'f':
							td->options |= OPT_USELOADRESTBOUND;
							td->restBoundInputFname = argv[i + 1];
							break;
						
						default:
							fprintf(stderr, "Unknown bound type '%c'.\n", argv[i][j]);
							exit(1);
					}
				}
				
				if (td->options & OPT_USELOADRESTBOUND) ++i;		// Advance past filename
				break;
				
			case '-':		// Extended option
				if (ProcessBooleanExtOption(argv[i] + 2, "displaynewbounds", &td->options, OPT_DISPLAYNEWBOUNDS)) break;
				if (ProcessBooleanExtOption(argv[i] + 2, "displaynewtrees", &td->options, OPT_DISPLAYNEWTREES)) break;
				if (ProcessUnsignedExtOption(argv[i] + 2, "updatesecs", &td->updateIntervalSecs)) break;
				if (ProcessUnsignedExtOption(argv[i] + 2, "startcol", &td->startCol)) {
					--td->startCol;
					break;
				}
				if (ProcessUnsignedExtOption(argv[i] + 2, "endcol", &td->endCol)) {
					--td->endCol;
					break;
				}
				if (ProcessBooleanExtOption(argv[i] + 2, "improveinitupperbound", &td->options, OPT_IMPROVEINITUPPERBOUND)) break;
				if (ProcessUnsignedExtOption(argv[i] + 2, "seed", &td->rngSeed)) break;
				if (ProcessBooleanExtOption(argv[i] + 2, "generatesitepaircostsfile", &do_generateSitePairCostsFile, 1)) break;
				if (ProcessBooleanExtOption(argv[i] + 2, "onlycomputelowerbound", &td->options, OPT_ONLYCOMPUTELOWERBOUND)) break;
				if (ProcessBooleanExtOption(argv[i] + 2, "keeptempfiles", &td->options, OPT_KEEPTEMPFILES)) break;
				if (ProcessUnsignedExtOption(argv[i] + 2, "pbmaxiters", &td->partBoundMaxItersNoImprovement)) break;
				if (ProcessBooleanExtOption(argv[i] + 2, "savereordereddataset", &td->options, OPT_SAVEREORDEREDDATASET)) break;
				if (ProcessBooleanExtOption(argv[i] + 2, "tbr", &td->options, OPT_TBR)) break;
				
				// Fall through to here if the argument did not match any extended option
				fprintf(stderr, "Unrecognised option '%s', aborting.\n", argv[i]);
				exit(1);
			
			case '?':
			case 'h':
				PrintUsage();
				return 0;
			
			case 0:					// i.e. just a "-"
				*pInFile = stdin;
				break;
			
			default:
				fprintf(stderr, "Unrecognised option '-%c', aborting.\n", argv[i][1]);
				exit(1);
			}
		} else {
			// Process data file specification
			*pInFileName = argv[i];
			
			if ((*pInFile = fopen(argv[i], "r")) == NULL) {
				fprintf(stderr, "File '%s' could not be opened for input, aborting.\n", argv[i]);
				exit(1);
			}
		}
	}
	
	if (!td->rngSeed) {
		// Choose a "random" random number generator seed ;)
		// Just time(NULL) is unsafe as > 1 run may be begun in < 1s.
		// Also "*" is much safer than "^", since getpid() might increment every 1s or so.
		td->rngSeed = time(NULL) * (getpid() + 1);
	}

	srand(td->rngSeed);
	td->rand1 = rand();		// These apparently pointless values are printed out in ReportSettings() to make it easier to
	td->rand2 = rand();		// catch cases where a run is performed on two different machines/compilers using the same
	td->rand3 = rand();		// seed and we have mistakenly assumed that the results will be the same.

	//HACK: This should not be called more than once, ever.
	if (do_generateSitePairCostsFile) {
		generateSitePairCostsFile();
	}

	if (!*pInFile) {
		fprintf(stderr, "Must specify PHYLIP-format file to open on the command line (or '-' to read data from standard input).\n");
		exit(1);
	}
	
	return 1;
}

void ReportSettings(struct TreeData *td, char *inFileName, FILE *inFile, char *outFileName, FILE *outFile) {
#ifdef TARGETMULTI
	int nProcs;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
#endif	// TARGETMULTI

	if (!(td->options & OPT_QUIET)) {
		fprintf(stderr, "%s\n\n", PROGLONGNAME);

#ifdef TARGETMULTI
		fprintf(stderr, "Parallel (MPI) version running %d processes.\n", nProcs);
#ifdef SYNCHRONOUSSENDS
		fprintf(stderr, "Synchronous-mode MPI sends will be used.\n");
#else	// not SYNCHRONOUSSENDS
		fprintf(stderr, "Standard-mode MPI sends will be used.\n");
#endif	// not SYNCHRONOUSSENDS
#else	// not TARGETMULTI
		fprintf(stderr, "Single-processor version.\n");
#endif	// not TARGETMULTI
		
		if (inFile == stdin) {
			fprintf(stderr, "Reading data from standard input.\n");
		} else {
			fprintf(stderr, "Reading data from file '%s'.\n", inFileName);
		}
		
		if (outFile == stdout) {
			fprintf(stderr, "Writing result trees to standard output.\n");
		} else {
			fprintf(stderr, "Writing result trees to file '%s'.\n", outFileName);
		}
		
		if (td->bound == INT_MAX) {
			fprintf(stderr, "No initial bound specified.\n");
		} else {
			fprintf(stderr, "Initial bound = %u.\n", td->bound);
		}
		
		if (td->maxTrees == INT_MAX) {
			fprintf(stderr, "All optimal trees will be saved.\n");
		} else {
			fprintf(stderr, "Up to %u optimal trees will be saved.\n", td->maxTrees);
		}
		
		if (td->startCol == 0 && td->endCol == INT_MAX) {
			fprintf(stderr, "All sites will be processed.\n");
		} else {
			fprintf(stderr, "Sites %u-%u will be processed.\n", td->startCol + 1, td->endCol + 1);
		}
		
		if (td->options & OPT_DISCARDMISSAMBIG) {
			fprintf(stderr, "Sites containing gaps or ambiguous characters will be discarded.\n");
		} else {
			fprintf(stderr, "Sites containing gaps or ambiguous characters will be retained.\n");
		}
		
		if (td->options & OPT_NEXUSFMT) {
			fprintf(stderr, "Trees will be output in NEXUS format.\n");
		} else {
			fprintf(stderr, "Trees will be output in raw (Newick) format.\n");
		}
		
		if (!(td->options & OPT_DISPLAYNEWBOUNDS)) {
			fprintf(stderr, "New bounds will not be reported as they occur.\n");
		}
		
		if (!(td->options & OPT_DISPLAYNEWTREES)) {
			fprintf(stderr, "New trees will not be reported as they occur.\n");
		}
		
		// WTJW 12/6/2004 (& 25/2/2005)
		if (td->options & OPT_USELOADRESTBOUND) {
			fprintf(stderr, "The restBound[] array will be loaded from the file '%s'.\n", td->restBoundInputFname);
		}
		if (td->options & OPT_USEDISTBASESBOUNDREST) {
			fprintf(stderr, "The distinct-bases-per-site bound will be used.\n");
		}
		if (td->options & OPT_USEINCOMPATBOUNDREST) {
			fprintf(stderr, "The pairs-of-incompatible-sites bound will be used.\n");
		}
		if (td->options & OPT_USEMSTBOUNDREST) {
			fprintf(stderr, "The minimum spanning tree bound will be used.\n");
		}
		if (td->options & OPT_USEKICKASSBOUNDREST) {
			fprintf(stderr, "The KICKASS bound will be used.\n");
		}
		if (td->options & OPT_USEKA2BOUNDREST) {
			fprintf(stderr, "The KA2 bound will be used.\n");
		}
		if (td->options & OPT_USEPARTBOUND) {
			fprintf(stderr, "The PARTBOUND bound will be used.  %s will be maximised.\n",
				td->partBoundStrategy == PBS_CLIPPEDANDADDITIVE ? "Final scores" :
				td->partBoundStrategy == PBS_CLIPPED ? "Clipped scores" :
				"Lower bounds"
			);
			if (td->partBoundMaxItersNoImprovement) {
				fprintf(stderr, "PARTBOUND iterations will stop after %d iterations with no improvement.\n", td->partBoundMaxItersNoImprovement);
			}
		}
		if (td->options & OPT_USEGREEDYHITTINGSETBOUND) {
			fprintf(stderr, "The greedy hitting set bound will be used.\n");
		}

		if (td->options & OPT_ONLYCOMPUTELOWERBOUND) {
			fprintf(stderr, "The program will exit after computing lower bounds.\n");
		}

		if (td->options & OPT_SAVEREORDEREDDATASET) {
			fprintf(stderr, "The dataset produced by reordering taxa and stripping constant and non-PI sites will be saved.\n");
		}

		if (td->options & OPT_KEEPTEMPFILES) {
			fprintf(stderr, "Temporary intermediate files will be kept.\n");
		}

		fprintf(stderr, "Using random number seed %u.  First 3 random numbers: %d, %d, %d.\n", td->rngSeed, td->rand1, td->rand2, td->rand3);
	}
}

void OutputResults(struct TreeData *td, FILE *outFile) {
	int i, j, k;
	char *treeBuf;
	int startOfTree;
	int count;
	int noTreesOutput = 1;
#ifdef TARGETMULTI
	MPI_Status mpiStatus;
#endif	// TARGETMULTI

	if (td->options & OPT_NEXUSFMT) {
		fprintf(outFile,
			"#NEXUS\n"
			"[ Trees produced by " PROGLONGNAME "\n"
			"> Each tree has score %u. ]\n"
			"BEGIN TREES;\n"
			"\tTRANSLATE\n", td->bound);
		
		for (i = 0; i < td->numTaxa; ++i) {
			fprintf(outFile, "\t\t%u\t%s%c\n", i + 1, td->cd->labels[i], i < td->numTaxa - 1 ? ',' : ';');
		}
		
		j = 1;
	}

	treeBuf = malloc(TREEBUFSIZE);
#ifdef TARGETMULTI
	// Retrieve tree data from workers in TREEBUFSIZE chunks.
	// We need to serialise workers here in order to get all chunks from each worker in order, but the performance
	// hit is too small to matter.
	for (i = 0; i < td->nWorkers; ++i) {
		DBGPRINT2("B: Waiting for W%04d to send us its tree data.\n", i + 1);
#else	// not TARGETMULTI
	if (td->tempTreeFile) {
		rewind(td->tempTreeFile);
#endif	// not TARGETMULTI
		startOfTree = 1;
		do {
#ifdef TARGETMULTI
			// Multi-processor: read a block from the ith worker.
			MPI_Recv(treeBuf, TREEBUFSIZE, MPI_CHAR, i + 1, MSG_TREES, MPI_COMM_WORLD, &mpiStatus);
			MPI_Get_count(&mpiStatus, MPI_CHAR, &count);
			DBGPRINT3("B: Received %d bytes of tree data from W%04d.\n", count, i + 1);
#else	// not TARGETMULTI
			// Single-processor: just read a block from our own temporary tree file.
			// We could special-case a rename operation when plain Newick format has been requestsed and tempTreeFile
			// lives on the same device as the final output file, but that's unnecessary optimisation.
			count = fread(treeBuf, 1, TREEBUFSIZE, td->tempTreeFile);
#endif	// not TARGETMULTI

			if (td->options & OPT_NEXUSFMT) {
				// Need to put some identifying stuff at the start of each line.
				for (k = 0; k < count; ++k) {
					if (startOfTree) {
						fprintf(outFile, "\tTREE FASTDNAMP_%u = [&U] ", j++);
						startOfTree = 0;
					}

					fputc(treeBuf[k], outFile);
					if (treeBuf[k] == '\n') {
						startOfTree = 1;
					}
				}
			} else {
				// Can just dump the data verbatim.
				fwrite(treeBuf, 1, count, outFile);
			}

			if (count) {
				noTreesOutput = 0;
			}
		} while (count == TREEBUFSIZE);
	}
	
	if (td->options & OPT_NEXUSFMT) {
		fprintf(outFile, "END;\n");
	}

	free(treeBuf);

	if (!(td->options & OPT_QUIET) && noTreesOutput) {
		fprintf(stderr, "NO BEST TREE FOUND!  (Have you specified an initial bound that is too low?)\n");
	}
}

// Only called by ReportTiming().  The fprintf() field widths have been chosen to agree with the label
// strings passed in from this function.
static void ReportOneTime(TimePoint start, TimePoint end, char *label) {
	char buf[200];			//HACK: I know, possible buffer overflow if I use too-long strings in future...
	double secs;
	
	secs = TimeDifference(start, end);
	if (secs < 60) {
		fprintf(stderr, "%-45s %11.2lfs\n", label, secs);
	} else {
		DurationToString(secs, 1, buf);
		fprintf(stderr, "%-45s %11.2lfs (%s)\n", label, secs, buf);
	}
}

void ReportTiming(struct TreeData *td) {
	if (!(td->options & OPT_QUIET)) {
		ReportOneTime(td->startTime, td->orderTaxaStartTime, "Time to read input:");
		ReportOneTime(td->orderTaxaStartTime, td->upperBoundStartTime, "Time to order taxa:");
		ReportOneTime(td->upperBoundStartTime, td->lowerBoundStartTime, "Time to determine initial upper bound:");
		ReportOneTime(td->lowerBoundStartTime, td->branchAndBoundStartTime, "Time to determine lower bounds:");
		ReportOneTime(td->branchAndBoundStartTime, td->outputResultsStartTime, "Time to perform branch and bound search:");
		ReportOneTime(td->outputResultsStartTime, td->endTime, "Time to write output:");
		ReportOneTime(td->startTime, td->endTime, "Total time taken:");
	}
}

// If precise = 0, the output string will only give a few "sig figs", e.g.  if
// the time is more than one day, no minutes or hours will be displayed.
// buf must be preallocated by the caller and must contain enough room for
// any output string.  80 bytes is definitely an upper bound on this length.
// Returns a pointer to the *end* of the string produced.
char *DurationToString(double secs, int precise, char *buf) {
	unsigned intSecs, days, hours, minutes;
	
	intSecs = (unsigned) secs;
	secs -= intSecs;

	days = intSecs / (24 * 60 * 60);
	intSecs -= days * 24 * 60 * 60;
	hours = intSecs / (60 * 60);
	intSecs -= hours * 60 * 60;
	minutes = intSecs / 60;
	intSecs -= minutes * 60;
	secs += intSecs;			// Now 0.0 <= secs < 60.0

	if (days > 0) {
		buf += sprintf(buf, "%u day%s, ", days, days == 1 ? "" : "s");
	}
	
	if (days > 0 || hours > 0) {
		buf += sprintf(buf, "%u hour%s%s", hours, hours == 1 ? "" : "s", days == 0 ? ", " : ".");
	}
	
	if (precise || days == 0) {
		if (hours > 0 || minutes > 0) {
			buf += sprintf(buf, "%u minute%s%s", minutes, minutes == 1 ? "" : "s", precise || hours == 0 ? ", " : ".");
		}
	}
	
	if (precise || days == 0) {
		if (precise || hours == 0) {
			buf += sprintf(buf, "%.2lf seconds.", secs);		// Yeah, we'll say "1.00 seconds".
		}
	}
	
	return buf;
}

void EstimateTimeRemaining(struct TreeData *td) {
	TimePoint now = GetTimeNow();
	double elapsedSecs = TimeDifference(td->branchAndBoundStartTime, now);
	double remainingSecs;
	char buf[80];
	
	fprintf(stderr, "Estimated time remaining: ");
	
	if (elapsedSecs < 10 || td->proportionComplete == 0) {
		fprintf(stderr, "(too early to tell)\n");
	} else {
		remainingSecs = ((COMPLETE_UNITS - td->proportionComplete) * elapsedSecs) / td->proportionComplete;
		
		DurationToString(remainingSecs, 0, buf);
		fprintf(stderr, "%s\n", buf);
	}
}

struct tree *ConstructBaseTree(struct TreeData *td, int t1, int t2, int t3) {
	struct tree *root;
	
	// Construct a base tree to work with
	root = NewNode(t1, 3, td);
	AddToPendantNode(root, LEFT, NewNode(-1, 3, td));
	AddToPendantNode(root->nodes[LEFT], LEFT, NewNode(t2, 3, td));
	AddToPendantNode(root->nodes[LEFT], RIGHT, NewNode(t3, 3, td));
#ifdef SMARTFITCH
#error "I haven't tested this since making changes on 5/10/2009."
	root->nodes[LEFT]->cumScore = FastFitch2Seqs(root->nodes[LEFT]->nodes[LEFT]->sequence, root->nodes[LEFT]->nodes[RIGHT]->sequence, root->nodes[LEFT]->sequence, td->weights, td->memSeqLen, td->heavySeqLen, td->byteWeights);
	root->cumScore = FastFitch2SeqsCost(root->nodes[LEFT]->sequence, root->sequence, td->weights, td->memSeqLen, td->heavySeqLen, td->byteWeights) + root->nodes[LEFT]->cumScore;
	DBGPRINT2("root->nodes[LEFT]->cumScore = <%d>\n", root->nodes[LEFT]->cumScore);		//DEBUG
	DBGPRINT2("root->cumScore = <%d>\n", root->cumScore);		//DEBUG
#endif	// SMARTFITCH
	{
		// Calculate internal node's UP sequence from its childrens' UP sequences.
		root->nodes[LEFT]->seqsByTreeSize[3].seqs[UP] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
		INCMEASURE(preorderFitchCallsFitchBasesCount);
		CHOSEN_FitchBases(
			root->nodes[LEFT]->nodes[LEFT]->seqsByTreeSize[3].seqs[UP],
			root->nodes[LEFT]->nodes[RIGHT]->seqsByTreeSize[3].seqs[UP],
			root->nodes[LEFT]->seqsByTreeSize[3].seqs[UP],
			td->seqLenBlocks
		);

		// Internal node's DOWN sequence is the same the root's UP sequence: just shallow-copy.
		root->nodes[LEFT]->seqsByTreeSize[3].seqs[DOWN] = root->seqsByTreeSize[3].seqs[UP];

		// Calculate internal node's MIDPOINT sequence.
		root->nodes[LEFT]->seqsByTreeSize[3].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
		INCMEASURE(preorderFitchCallsFitchBasesCount);
		CHOSEN_FitchBases(
			root->nodes[LEFT]->seqsByTreeSize[3].seqs[UP],
			root->nodes[LEFT]->seqsByTreeSize[3].seqs[DOWN],
			root->nodes[LEFT]->seqsByTreeSize[3].seqs[MIDPOINT],
			td->seqLenBlocks
		);

		// Calculate the left leaf's DOWN sequence from its parent's DOWN and its sibling's UP.
		root->nodes[LEFT]->nodes[LEFT]->seqsByTreeSize[3].seqs[DOWN] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
		INCMEASURE(preorderFitchCallsFitchBasesCount);
		CHOSEN_FitchBases(
			root->nodes[LEFT]->seqsByTreeSize[3].seqs[DOWN],
			root->nodes[LEFT]->nodes[RIGHT]->seqsByTreeSize[3].seqs[UP],
			root->nodes[LEFT]->nodes[LEFT]->seqsByTreeSize[3].seqs[DOWN],
			td->seqLenBlocks
		);

		// Calculate the left leaf's MIDPOINT sequence.
		root->nodes[LEFT]->nodes[LEFT]->seqsByTreeSize[3].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
		INCMEASURE(preorderFitchCallsFitchBasesCount);
		CHOSEN_FitchBases(
			root->nodes[LEFT]->nodes[LEFT]->seqsByTreeSize[3].seqs[UP],
			root->nodes[LEFT]->nodes[LEFT]->seqsByTreeSize[3].seqs[DOWN],
			root->nodes[LEFT]->nodes[LEFT]->seqsByTreeSize[3].seqs[MIDPOINT],
			td->seqLenBlocks
		);

		// Calculate the right leaf's DOWN sequence from its parent's DOWN and its sibling's UP.
		root->nodes[LEFT]->nodes[RIGHT]->seqsByTreeSize[3].seqs[DOWN] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
		INCMEASURE(preorderFitchCallsFitchBasesCount);
		CHOSEN_FitchBases(
			root->nodes[LEFT]->seqsByTreeSize[3].seqs[DOWN],
			root->nodes[LEFT]->nodes[LEFT]->seqsByTreeSize[3].seqs[UP],
			root->nodes[LEFT]->nodes[RIGHT]->seqsByTreeSize[3].seqs[DOWN],
			td->seqLenBlocks
		);

		// Calculate the right leaf's MIDPOINT sequence.
		root->nodes[LEFT]->nodes[RIGHT]->seqsByTreeSize[3].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
		INCMEASURE(preorderFitchCallsFitchBasesCount);
		CHOSEN_FitchBases(
			root->nodes[LEFT]->nodes[RIGHT]->seqsByTreeSize[3].seqs[UP],
			root->nodes[LEFT]->nodes[RIGHT]->seqsByTreeSize[3].seqs[DOWN],
			root->nodes[LEFT]->nodes[RIGHT]->seqsByTreeSize[3].seqs[MIDPOINT],
			td->seqLenBlocks
		);

		// Calculate scores.
		td->score = CHOSEN_FitchScore(
			root->nodes[LEFT]->nodes[LEFT]->seqsByTreeSize[3].seqs[UP],
			root->nodes[LEFT]->nodes[RIGHT]->seqsByTreeSize[3].seqs[UP],
			td->weights,
			td->seqLenBlocks,
			INT_MAX
		);
		DBGPRINT2("Score from children of root's left child = <%d>\n", td->score);		//DEBUG
		td->score += CHOSEN_FitchScore(
			root->nodes[LEFT]->seqsByTreeSize[3].seqs[UP],
			root->seqsByTreeSize[3].seqs[UP],
			td->weights,
			td->seqLenBlocks,
			INT_MAX
		);
		DBGPRINT2("Score from entire initial 3-taxon tree = <%d>\n", td->score);		//DEBUG
	}
	
	// Nowadays we also add the score of all sites that are non-PI on the full tree in here.
	// That means that td->score is always actually >= the true score of the tree so far --
	// but that true score can always be recovered by subtracting td->constantWeight.
	td->score += td->constantWeight;
	DBGPRINT5("Score from entire initial 3-taxon tree, containing taxa at positions %d, %d and %d, including weight from all sites that are non-PI on the full tree = <%d>\n", t1, t2, t3, td->score);

	DBGPRINT2("td->seqLenBlocks = <%d>\n", td->seqLenBlocks);		//DEBUG
	DBGPRINT2("td->seqLen = <%d>\n", td->seqLen);		//DEBUG
	return root;
}

void PrintUsage(void) {
	fprintf(stderr,
		"Syntax: " PROGNAME " [<options>] <phylipfile>\n"
		"  -o <fname>   File to write results to, instead of standard output\n"
		"  -b <nn>      Specify an initial upper bound (or 0 for no bound)\n"
		"  -m <nn>      Maximum number of trees to save (or 0 for no limit)\n"
		"  -r           Output trees in raw (Newick) format (default: NEXUS format)\n"
		"  -d           Discard sites containing gaps/ambiguous states\n"
		"                 (default: like GAPMODE=MISSING in PAUP*)\n"
		"  -q           Quiet mode: only the final set of trees is output\n"
		"  -B <fname>   Load the restBound[] array from <fname>\n"
		"  --startcol=<nn>           Starting site number (first site is site 1)\n"
		"  --endcol=<nn>             Ending site number\n"
		"  --displaynewbounds={Y|N}  Control printing of updated bounds while running\n"
		"  --displaynewtrees={Y|N}   Control printing of best trees while running\n"
		"  --updatesecs=<nn>         Update progress every <nn> seconds (use 0 to\n"
		"                              turn off periodic progress updates)\n"
		"  --seed=<nn>               Use <nn> as the random number generation seed --\n"
		"                              0 (the default) means 'pick a random seed'\n"
	);
}

struct tree *NewNode(int label, int nTaxaOnTree, struct TreeData *td) {
	struct tree *t;
	
	INCMEASURE(newNodeCount);
	
	t = malloc(sizeof (struct tree));
	t->nodes[LEFT] = NULL;
	t->nodes[RIGHT] = NULL;
	t->contents = label;
#ifdef SMARTFITCH
	t->cumScore = 0;		// All leaves have cumulative score 0
#endif	// SMARTFITCH
	t->seqsByTreeSize = ArenaAllocate(&td->alignedMem, UROUNDUP(sizeof (struct edge_seqs) * (td->numTaxa + 1), BYTESPERBLOCK));
	if (label == -1) {
		//ALIGNED_MALLOC(t->sequence, td->offsets[td->numTaxa] * 3 * sizeof (unsigned));
		//ALIGNED_MALLOC(t->sequence, td->offsets[td->numTaxa] * BYTESPERBLOCK * 3);
		// We don't actually allocate any sequence buffer space here -- that will be done as necessary
		// by ConstructBaseTree() and UpdateTreeAfterInsertion() (although arguably the latter should
		// actually be merged into this function).
		//t->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK * 3);
		//t->seqsByTreeSize[nTaxaOnTree].seqs[UP] = t->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] + td->seqLenBlocks * BYTESPERBLOCK;
		//t->seqsByTreeSize[nTaxaOnTree].seqs[DOWN] = t->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] + td->seqLenBlocks * BYTESPERBLOCK * 2;
	} else {
		//t->sequence = td->taxonSequence[label - 1];
		//t->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK * 2);
		t->seqsByTreeSize[nTaxaOnTree].seqs[UP] = td->charMat + (label - 1) * td->seqLenBlocks * BYTESPERBLOCK;		// Will never change
		//t->seqsByTreeSize[nTaxaOnTree].seqs[DOWN] = t->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] + td->seqLenBlocks * BYTESPERBLOCK;
	}
#ifdef BACKPOINTERS
	t->nodes[PARENT] = NULL;
#endif	// BACKPOINTERS
	return t;
}

void AddToPendantNode(struct tree *parentNode, int child, struct tree *childNode) {
	if (parentNode->nodes[child] == NULL) {
		parentNode->nodes[child] = childNode;
#ifdef BACKPOINTERS
		childNode->nodes[PARENT] = parentNode;
#endif	// BACKPOINTERS
	}
	else
		DBGPRINT1("Error in AddToPendantNode\n");
	return;
}

// Inserts a taxon node in an edge.  The edge is specified by the parent node
// and a child number.
struct tree * FASTCALL InsertNode(struct tree *parent, int child, int label, int nTaxaOnTree, struct TreeData *td) {
	struct tree *newInternal;
	struct tree *newExternal;
	
	INCMEASURE(insertNodeCount);
	
	//DBGPRINT3("Adding node to parent at address %p (contents = %d).\n", parent, parent->contents);	//WTJW
	newExternal = malloc(sizeof(struct tree));
	newExternal->nodes[LEFT] = NULL;
	newExternal->nodes[RIGHT] = NULL;
	newExternal->contents = label;

	newInternal = malloc(sizeof(struct tree));
	newInternal->nodes[LEFT] = newExternal;		//Arbitrarily always add 
	newInternal->nodes[RIGHT] = parent->nodes[child]; //new node on left
	newInternal->contents = -1;
#ifdef BACKPOINTERS
	newInternal->nodes[PARENT] = parent;			// 7/2/2003: Fixed bug
	newExternal->nodes[PARENT] = newInternal;
	parent->nodes[child]->nodes[PARENT] = newInternal;	// 31/3/2003: Fixed another bug!
#endif	// BACKPOINTERS
	//HACK: not yet deallocating these properly in DestroyTree().  Only internal nodes should be deallocated.
	//ALIGNED_MALLOC(newInternal->sequence, td->offsets[td->numTaxa] * 3 * sizeof (unsigned));
	//ALIGNED_MALLOC(newInternal->sequence, td->offsets[td->numTaxa] * BYTESPERBLOCK * 3);
	// We don't actually allocate any sequence buffer space here -- that will be done as necessary
	// by ConstructBaseTree() and UpdateTreeAfterInsertion() (although arguably the latter should
	// actually be merged into this function).
	newInternal->seqsByTreeSize = ArenaAllocate(&td->alignedMem, UROUNDUP(sizeof (struct edge_seqs) * (td->numTaxa + 1), BYTESPERBLOCK));
	//newInternal->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK * 3);
	//newInternal->seqsByTreeSize[nTaxaOnTree].seqs[UP] = newInternal->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] + td->seqLenBlocks * BYTESPERBLOCK;
	//newInternal->seqsByTreeSize[nTaxaOnTree].seqs[DOWN] = newInternal->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] + td->seqLenBlocks * BYTESPERBLOCK * 2;
	////newExternal->sequence = td->taxonSequence[label - 1];
	newExternal->seqsByTreeSize = ArenaAllocate(&td->alignedMem, UROUNDUP(sizeof (struct edge_seqs) * (td->numTaxa + 1), BYTESPERBLOCK));
	//newExternal->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK * 2);
	// We need UP sequences for both the "old" and "new" trees, since when pushing changes down
	// through the tree, UpdateTree() looks at the sibling's *old* UP sequence.
	newExternal->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP] =
		newExternal->seqsByTreeSize[nTaxaOnTree].seqs[UP] = td->charMat + (label - 1) * td->seqLenBlocks * BYTESPERBLOCK;		// Never changes
	//newExternal->seqsByTreeSize[nTaxaOnTree].seqs[DOWN] = newExternal->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] + td->seqLenBlocks * BYTESPERBLOCK;
#ifdef SMARTFITCH
	newExternal->cumScore = 0;			// All leaves have cumulative score 0
	newInternal->cumScore = 0;			// This should be recomputed ASAP
#endif	// SMARTFITCH
	parent->nodes[child] = newInternal;
#ifdef DYNAMICORDER
	UpdateStates(label - 1, 1, td);
#endif	// DYNAMICORDER
	return newInternal;		//WTJW: for convenience
}

#ifdef DYNAMICORDER
// Updates the treeStateCounts and restStateCounts arrays for computing the
// number-of-distinct-bases-per-site-minus-one bound when DYNAMICORDER is on.
//HACK: for the time being, we completely ignore sites containing ambiguous
// base sets.
void UpdateStates(unsigned taxon, int dir, struct TreeData *td) {
	static unsigned map[] = {
		0,		// 0 = n/a
		0,		// 1 = A
		1,		// 2 = C
		0,		// 3 = n/a
		2,		// 4 = G
		0,		// 5 = n/a
		0,		// 6 = n/a
		0,		// 7 = n/a
		3		// 8 = T
	};
	
	unsigned i, idx;
	
	for (i = 0; i < td->seqLen; ++i) {
		if (!td->safeSite[i]) continue;			//HACK: Ignore sites containing ambiguous base sets
		
		idx = map[getBaseAt(taxon, i, td)];
		treeStateCounts[idx][i] += dir;
		restStateCounts[idx][i] -= dir;
	}
}

//HACK: get a better name
unsigned ComputeDynamicBound(struct TreeData *td) {
	unsigned weight, bound = 0;
	
	for (i = 0; i < td->seqLen; ++i) {
		if (!td->safeSite[i]) continue;			//HACK: Ignore sites containing ambiguous base sets
		
		weight =
			(restStateCounts[0][i] && !treeStateCounts[0][i]) + 
			(restStateCounts[1][i] && !treeStateCounts[1][i]) + 
			(restStateCounts[2][i] && !treeStateCounts[2][i]) + 
			(restStateCounts[3][i] && !treeStateCounts[3][i]);
		
		bound += weight * td->weights[i];
	}
	
	return bound;
}
#endif	// DYNAMICORDER

void FASTCALL DeleteNode(struct tree *parent, int child) {
	struct tree *deadInternal;
	struct tree *deadExternal;
	
	deadInternal = parent->nodes[child];
	deadExternal = deadInternal->nodes[LEFT];	//arbitrarily delete left
												//see InsertNode

	parent->nodes[child] = deadInternal->nodes[RIGHT];
#ifdef BACKPOINTERS
	deadInternal->nodes[RIGHT]->nodes[PARENT] = parent;		// Bugfix 31/3/2003
#endif	// BACKPOINTERS
	// Sequences are allocated using the arena and freed all at once.
	//ALIGNED_FREE(deadInternal->sequence);
	
	free(deadInternal);
	free(deadExternal);
	
	return;
}

#ifdef BACKPOINTERS
int VerifyTree(struct tree *root, struct tree *parent) {
	int ok = 1;
	
	if (root->nodes[PARENT] != parent) {
		return 0;
	}
	
	if (root->nodes[LEFT] != NULL) {
		ok &= VerifyTree(root->nodes[LEFT], root);
	}
	
	if (root->nodes[RIGHT] != NULL) {
		ok &= VerifyTree(root->nodes[RIGHT], root);
	}
	
	return ok;
}
#endif	// BACKPOINTERS

// Can't handle masks that contain zero or more than one base.
char GetBaseFromMask(unsigned char mask) {
	switch (mask) {
	case BM_A: return 'A';
	case BM_C: return 'C';
	case BM_G: return 'G';
	case BM_T: return 'T';		// May have originally been a 'U' but who cares
#ifdef SHOWSHORTBASESETS
	case BM_A | BM_G:               return 'R';
	case BM_C | BM_T:               return 'Y';
	case BM_C | BM_G:               return 'S';
	case BM_A | BM_T:               return 'W';
	case BM_G | BM_T:               return 'K';
	case BM_A | BM_C:               return 'M';
	case BM_C | BM_G | BM_T:        return 'B';
	case BM_A | BM_G | BM_T:        return 'D';
	case BM_A | BM_C | BM_T:        return 'H';
	case BM_A | BM_C | BM_G:        return 'V';
	case BM_A | BM_C | BM_G | BM_T: return 'N';
#endif	// SHOWSHORTBASESETS
	case 0: return '0';		//DEBUG
	default:
		DBGPRINT2("Error! GetBaseFromMask() called with a mask (0x%x) that contains zero or more than one base, aborting.\n", mask);
		exit(1);
		return -1;		// Never runs, exists to satisfy the compiler
	}
}

// Returns a bit mask with a 1-bit in the position corresponding to base.
unsigned char GetMaskFromBase(char base) {
	switch (toupper(base)) {
	case 'A': return BM_A;
	case 'C': return BM_C;
	case 'G': return BM_G;
	case 'U': // Fall through
	case 'T': return BM_T;
	case 'R': return BM_A | BM_G;
	case 'Y': return BM_C | BM_T;
	case 'S': return BM_C | BM_G;
	case 'W': return BM_A | BM_T;
	case 'K': return BM_G | BM_T;
	case 'M': return BM_A | BM_C;
	case 'B': return BM_C | BM_G | BM_T;
	case 'D': return BM_A | BM_G | BM_T;
	case 'H': return BM_A | BM_C | BM_T;
	case 'V': return BM_A | BM_C | BM_G;
	case '-':	// PAUP* treat a gap the same as 'N' when GAPMODE=MISSING, so we will do the same.
	case '?':	// Always treat this as 'N'.
	case 'N': return BM_A | BM_C | BM_G | BM_T;
	default:
		DBGPRINT2("Error! GetMaskFromBase() called with invalid character '%c', aborting.\n", base);
		exit(1);
		return -1;		// Never runs, exists to satisfy the compiler
	}
}

// Returns a code from 0 to 3, inclusive.
int GetCodeFromBase(char base) {
	switch (toupper(base)) {
	case 'A': return 0;
	case 'C': return 1;
	case 'G': return 2;
	case 'U': // Fall through
	case 'T': return 3;
	default:
		DBGPRINT2("Error! GetCodeFromBase() called with invalid character '%c', aborting.\n", base);
		exit(1);
		return -1;		// Never runs, exists to satisfy the compiler
	}
}

void PrintTree(struct tree *root, struct TreeData *td, FILE *f) {
	//DBGPRINT1("Call to PrintTree\n");					//DEBUG

	if ((root->nodes[LEFT] != NULL) && (root->nodes[RIGHT] == NULL)) {	//ie this is the root node
		fprintf(f, "(");
#ifdef DEBUG_SHOWUNRENAMEDTAXONLABELS
		fprintf(f, "%d", root->contents);
#else	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		fprintf(f, "%d", td->taxonMap[root->contents - 1] + 1);
#endif	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		fprintf(f, ",");
		PrintTree(root->nodes[LEFT], td, f);
		fprintf(f, ")");	
	} else if ((root->nodes[LEFT] == NULL) && (root->nodes[RIGHT] != NULL))	{	//ie this is the root node
		fprintf(f, "(");
#ifdef DEBUG_SHOWUNRENAMEDTAXONLABELS
		fprintf(f, "%d", root->contents);
#else	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		fprintf(f, "%d", td->taxonMap[root->contents - 1] + 1);
#endif	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		fprintf(f, ",");
		PrintTree(root->nodes[RIGHT], td, f);
		fprintf(f, ")");	
	} else {	//This is a binary node
		if ((root->nodes[LEFT] == NULL) && (root->nodes[RIGHT] == NULL)) {
#ifdef DEBUG_SHOWUNRENAMEDTAXONLABELS
			fprintf(f, "%d", root->contents);
#else	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
			fprintf(f, "%d", td->taxonMap[root->contents - 1] + 1);
#endif	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		} else {
			fprintf(f, "(");
			if (root->nodes[LEFT] != NULL)
				PrintTree(root->nodes[LEFT], td, f);
			fprintf(f, ",");
			if (root->nodes[RIGHT] != NULL)
				PrintTree(root->nodes[RIGHT], td, f);
			fprintf(f, ")");
		}
	}
	return;	
}

// The caller must provide a buffer that is known to be large enough.
// Returns a pointer to the trailing nul character.
// Converted to use sprintf() instead of (nonstandard?) itoa().
//HACK: maybe should implement PrintTree() in terms of this?
char *TreeToString(struct tree *root, char *buf, struct TreeData *td) {
	//DBGPRINT1("Call to PrintTree\n");					//DEBUG

	if ((root->nodes[LEFT] != NULL) && (root->nodes[RIGHT] == NULL)) {	//ie this is the root node
		*buf++ = '(';
#ifdef DEBUG_SHOWUNRENAMEDTAXONLABELS
		buf += sprintf(buf, "%d", root->contents);
#else	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		buf += sprintf(buf, "%d", td->taxonMap[root->contents - 1] + 1);
#endif	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		*buf++ = ',';
		buf = TreeToString(root->nodes[LEFT], buf, td);
		*buf++ = ')';
	} else if ((root->nodes[LEFT] == NULL) && (root->nodes[RIGHT] != NULL))	{	//ie this is the root node
		*buf++ = '(';
#ifdef DEBUG_SHOWUNRENAMEDTAXONLABELS
		buf += sprintf(buf, "%d", root->contents);
#else	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		buf += sprintf(buf, "%d", td->taxonMap[root->contents - 1] + 1);
#endif	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		*buf++ = ',';
		buf = TreeToString(root->nodes[RIGHT], buf, td);
		*buf++ = ')';
	} else {	//This is a binary node
		if ((root->nodes[LEFT] == NULL) && (root->nodes[RIGHT] == NULL)) {
#ifdef DEBUG_SHOWUNRENAMEDTAXONLABELS
			buf += sprintf(buf, "%d", root->contents);
#else	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
			buf += sprintf(buf, "%d", td->taxonMap[root->contents - 1] + 1);
#endif	// not DEBUG_SHOWUNRENAMEDTAXONLABELS
		} else {
			*buf++ = '(';
			if (root->nodes[LEFT] != NULL)
				buf = TreeToString(root->nodes[LEFT], buf, td);
			*buf++ = ',';
			if (root->nodes[RIGHT] != NULL)
				buf = TreeToString(root->nodes[RIGHT], buf, td);
			*buf++ = ')';
		}
	}
	*buf = 0;		// Don't advance pointer
	return buf;
}

// Takes a non-transposed matrix and produces an un-transposed matrix whose
// sequences are all padded to a multiple of BYTESPERBLOCK bytes.
// The weights are also adjusted.  fromBuf and fromWeights must be different
// from *pToBuf and *pToWeights.  Returns the length of a sequence in
// 32-bit words.
// NO LONGER TRUE: 3 buffers are allocated per sequence, for use with the Yan-Bader algorithm.
// 2/10/2009: This is where the logic that converts "normal" 1-base-per-byte sequences
// into the 2-bases-per-byte "SQUEEZEBASES" format.  The original SqueezeBases() routine is
// actually no longer called, since it is only called from inside OrderColumns() which is
// no longer called from anywhere.
static void BlockifySequence_from_b1_w1(unsigned char *fromBuf, unsigned char *toBuf, unsigned seqLen) {
	unsigned s;
	unsigned char base;
	unsigned newSeqBlocks = GetNumBlocksFromSeqLen(seqLen);		// Probably computed by our caller too, but who cares.
	unsigned newSeqLen = newSeqBlocks * BYTESPERBLOCK * SITESPERBYTE;
	
	// Pad bases
	for (s = 0; s < newSeqLen; ++s) {
		if (s < seqLen) {
			base = fromBuf[s];
		} else {
			base = 0x0F;			// Important to use non-zero base filler for FASTWEIGHT1FITCH to work correctly
		}
		
#if SITESPERBYTE == 1
		toBuf[s] = base;
#elif SITESPERBYTE == 2
		if (s & 1) {
			toBuf[s / 2] |= base << 4;
		} else {
			toBuf[s / 2] = base;
		}
#else	// SITESPERBYTE != 1 && SITESPERBYTE != 2
#error SITESPERBYTE must be 1 or 2.
#endif	// SITESPERBYTE != 1 && SITESPERBYTE != 2
	}
}

// Convert the td->nTaxa-by-td->seqLen-byte rectangle of _b1_w1 taxon sequence data in td->charMat to a
// td->nTaxa by (td->seqLenBlocks * BYTESPERBLOCK)-byte rectangle in the format determined
// by config parameters (e.g. _b2_w4).
// We set td->seqLenBlocks here; it should not be accessed before this function is called.
// We also create the memory arena td->alignedMem.
// NOTE: td->charMat and td->weights, which were previously allocated with malloc(), are free()ed
// and reallocated using the arena.
static void BlockifyTaxa(struct TreeData *td) {
	unsigned newSeqBytes;
	unsigned newSeqLen;
	unsigned char *newMat;
	unsigned *newWeights;
	unsigned i;

	//td->seqLenBlocks = (!td->seqLen - 1) & ((td->seqLen - 1) / (SITESPERBYTE * BYTESPERBLOCK) + 1);		// Safely round up
	td->seqLenBlocks = GetNumBlocksFromSeqLen(td->seqLen);
	newSeqBytes = td->seqLenBlocks * BYTESPERBLOCK;
	newSeqLen = newSeqBytes * SITESPERBYTE;
	
	// Figure out how much memory our arena needs.
	ArenaCreate(&td->alignedMem,
		(2 * td->numTaxa - 3 + 1) *		// Number of edges, +1 for 2-part edge at root
		(
			UROUNDUP(sizeof (struct edge_seqs) * (td->numTaxa + 1), BYTESPERBLOCK) +	// Per-edge array of triples of pointers
			td->seqLenBlocks * BYTESPERBLOCK * 3 * (td->numTaxa + 1)	// History of UP, DOWN, MID sequences for this edge
		)
	);
	
	// Pad sequences
	newMat = ArenaAllocate(&td->alignedMem, td->numTaxa * td->seqLenBlocks * BYTESPERBLOCK);
	for (i = 0; i < td->numTaxa; ++i) {
		BlockifySequence_from_b1_w1(td->charMat + i * td->seqLen, newMat + i * td->seqLenBlocks * BYTESPERBLOCK, td->seqLen);
	}

	free(td->charMat);
	td->charMat = newMat;

	// Pad weights
	// Yes, this memory block is aligned (since newSeqLen is a multiple of BYTESPERBLOCK), which is necessary
	// since e.g. SSE2ASMFITCH will access the weights 128 bits at a time.
	newWeights = (unsigned *) ArenaAllocate(&td->alignedMem, newSeqLen * sizeof (unsigned));
	for (i = 0; i < newSeqLen; ++i) {
		newWeights[i] = i < td->seqLen ? td->weights[i] : 1;
	}

	free(td->weights);
	td->weights = newWeights;
}

// The following fields should be populated before calling InitTreeData():
// options, bound, maxTrees.
void InitTreeData(struct CharData *cd, struct TreeData *td) {
	unsigned char *buf;
	unsigned nConstantSites, nNonPiSites;
	
	if (cd->numTaxa < 4) {
		fprintf(stderr, "You must have four or more taxa to perform this analysis.\n");
		exit(1);
	}
	
	// Initialise
	td->cd = cd;
	td->bestList = NULL;
	td->numTrees = 0;
	td->optimalTreesMissed = 0;
	td->numTaxa = cd->numTaxa;
	td->seqLen  = cd->seqLen;
	td->seqLenBlocks = UINT_MAX;		// Should not be accessed until BlockifyTaxa() has been called
	td->proportionComplete = 0;
	td->treesActuallyExamined = 0;
	td->taxonMap = NULL;
#ifdef TARGETMULTI
	td->rank = -1;			// Indicate it's not initialised
	td->nWorkers = 0;
	td->mpiRequests = NULL;
	td->nMpiRequests = 0;
	td->pendingBoundBroadcast = 0;
	td->receiveBuf = NULL;
	td->stealBuf = NULL;
#endif	// TARGETMULTI
#ifdef LOGTREES
	td->treeLogFile = fopen("treelog.txt", "wb");		//HACK: Binary mode just because it will speed things up hopefully...
	fprintf(td->treeLogFile, "NTaxa\tScore\trestBound\tTotal\tBound\tTree\n");
#endif	// LOGTREES
	td->workStack = (unsigned (*)[2]) malloc((td->numTaxa + 1) * 2 * sizeof (unsigned));
	td->jobTreeSize = 2;		// This is appropriate for the single-CPU setting: no fast-forwarding will occur.
	td->tempTreeFName = malloc(MY_PATH_MAX);		//HACK: Ugh.  See the "explanation" at the top.
	// We *decide* the filename here, but only create it when we first need it -- hopefully this will reduce
	// filesystem contention on large-scale systems using a shared mount.
	GenerateTempFileName(NULL, "fastdnamp", td->tempTreeFName);
	DBGPRINT2("Temporary filename for trees: %s\n", td->tempTreeFName);
	td->tempTreeFile = NULL;
	
	buf = (unsigned char *) malloc(td->numTaxa * td->seqLen);
	td->seqLen = RecodeMatrix_b1t(cd->charMat, buf, td->seqLen, td->numTaxa, td->startCol, td->endCol, &td->weights, td->options & OPT_DISCARDMISSAMBIG);
	DBGPRINT2("After RecodeMatrix_b1t(): td->seqLen = <%d>\n", td->seqLen);		//DEBUG
	
	SortSitesAndWeights_b1t(buf, buf, td->weights, td->weights, td->seqLen, td->numTaxa, ColumnCompareIncludingWeights);		// Necessary before any call to CombineIdenticalSites()
	td->seqLen = CombineIdenticalSites_b1t(buf, td->weights, td->seqLen, td->numTaxa, td->numTaxa);
	DBGPRINT2("After CombineIdenticalSites_b1t(): td->seqLen = <%d>\n", td->seqLen);		//DEBUG
	
	SortSitesAndWeights_b1t(buf, buf, td->weights, td->weights, td->seqLen, td->numTaxa, ColumnCompareWeights);
	
	// Get rid of sites that are parsimoniously uninformative over the entire
	// tree.  Do this before OrderTaxa, so that the taxon order is as good
	// as possible.
	td->seqLen = RemoveUninfSites_b1t(buf, td->weights, td->seqLen, td->numTaxa, td->numTaxa, &nConstantSites, &nNonPiSites, &td->constantWeight);		//HACK: the name "td->constantWeight" is misleading -- this in fact refers to the constant amount of weight added to the tree by the non-PI sites.
	DBGPRINT2("After RemoveUninfSites_b1t(): td->seqLen = <%d>\n", td->seqLen);		//DEBUG
#ifdef TARGETMULTI
	if (!(td->options & OPT_QUIET) && !td->rank) {
#else	// not TARGETMULTI
	if (!(td->options & OPT_QUIET)) {
#endif	// not TARGETMULTI
		fprintf(stderr, "%d constant sites.\n", nConstantSites);
		fprintf(stderr, "%d non-parsimoniously-informative sites, having a total weight of %d.\n", nNonPiSites, td->constantWeight);
		fprintf(stderr, "%d parsimoniously-informative sites.\n", td->seqLen);
	}
	
	td->charMat = (unsigned char *) malloc(td->numTaxa * td->seqLen);
	memcpy(td->charMat, buf, td->numTaxa * td->seqLen);
	Transpose_b1(td->charMat, td->numTaxa, td->seqLen);
	
	td->orderTaxaStartTime = GetTimeNow();
	OrderTaxa(td);
	
#ifdef DEBUG
	WritePhylipAlignment(td->charMat, td->seqLen, td->numTaxa, td->weights, td->taxonMap, stderr, WPA_PRINTWEIGHTS);
#endif	// DEBUG
	if (td->options & OPT_SAVEREORDEREDDATASET) {
		FILE *reorderedDatasetFile = fopen("reordered.phy", "w");
		if (!reorderedDatasetFile) {
			perror("Could not open 'reordered.phy' for output");
			exit(1);
		}

		// Produce an actual legitimate PHYLIP file that other apps can read.
		// Unfortunately this means replicating weighted sites, since PHYLIP format only allows for site weights
		// up to 9 and I suspect that many apps can't handle even that.
		WritePhylipAlignment(td->charMat, td->seqLen, td->numTaxa, td->weights, td->taxonMap, reorderedDatasetFile, 0);
		fclose(reorderedDatasetFile);
	}
	
	td->upperBoundStartTime = GetTimeNow();
	if (td->options & OPT_IMPROVEINITUPPERBOUND) {
		ImproveInitialUpperBound(td);
	}
	
	// Establish the restBound[] array.
	td->lowerBoundStartTime = GetTimeNow();
	ComputeBounds(td);

	//FillTaxonSequence(td);
	BlockifyTaxa(td);

	//TODO: The time spent here is currently allocated to LB computation time, not UB computation time.
	// Should probably move the call to plain old ImproveInitialUpperBound() here and swap the calls to
	// GetTimeNow() around.
	if (td->options & OPT_TBR) {
		ImproveInitialUpperBoundViaTbr(td);
	}

	DBGPRINT2("td->seqLenBlocks = <%d>\n", td->seqLenBlocks);		//DEBUG
}

// Options is a bitmask.  WPA_PRINTWEIGHTS causes the weights to be written
// out underneath, instead of replicating columns as necessary.
// WPA_ALIGN causes character columns to be aligned with their weights.
// To produce actual Phylip-format output, set options to zero.
// Also, weights may be NULL, in which case all sites are considered to have
// weight 1.  Also, taxonMap may be NULL, in which case taxon names begin
// with a U and simply number up from 1.
void WritePhylipAlignment(unsigned char *buf, unsigned seqLen, unsigned numTaxa, unsigned *weights, unsigned *taxonMap, FILE *outFile, unsigned options) {
	unsigned i, j, k;
	unsigned totalSeqLen = seqLen;
	char ch;
	unsigned weight;
	
	if (!(options & WPA_PRINTWEIGHTS)) {
		totalSeqLen = 0;
		for (i = 0; i < seqLen; ++i) {
			totalSeqLen += weights[i];
		}
	}
	
	fprintf(outFile, "%d %d\n", numTaxa, totalSeqLen);
	for (i = 0; i < numTaxa; ++i) {
		if (taxonMap) {
			fprintf(outFile, "T%-9d", taxonMap[i] + 1);
		} else {
			fprintf(outFile, "U%-9d", i + 1);			// Indicate that the taxon number is meaningless
		}
		
		for (j = 0; j < seqLen; ++j) {
			ch = GetBaseFromMask(buf[i * seqLen + j]);
			
			if (options & WPA_PRINTWEIGHTS) {
				if (options & WPA_ALIGN) {
					fprintf(outFile, "  %c ", ch);		// Line up columns with weights
				} else {
					fputc(ch, outFile);
				}
			} else {
				if (weights) {
					weight = weights[j];
				} else {
					weight = 1;					// No weights provided
				}
				
				for (k = 0; k < weight; ++k) {
					fputc(ch, outFile);
				}
			}
		}
		
		fputc('\n', outFile);
	}
	
	if (weights && (options & WPA_PRINTWEIGHTS)) {
		for (i = 0; i < seqLen; ++i) {
			fprintf(outFile, "%3d ", weights[i]);
		}
		
		fputc('\n', outFile);
	}
}

// Recodes a character matrix, transposing it in the process (this is necessary
// for gap-containing columns to be removed efficiently).  fromBuf must not be
// equal to toBuf.  Returns the new number of sites.
// Also allocates a weight vector, which is initially all 1.
// The _b1t suffix is used because the *output* of the function is in this form.
int RecodeMatrix_b1t(unsigned char *fromBuf, unsigned char *toBuf, unsigned seqLen, unsigned numTaxa, unsigned startCol, unsigned endCol, unsigned **pWeights, unsigned discardMissAmbig) {
	// Recode bases
	unsigned newSeqLen = 0;				// This will grow as necessary
	unsigned cList[4];
	int highest;
	int map[4];
	unsigned ambiguous;
	unsigned discardedSites = 0;
	unsigned char intersection;
	unsigned keepSite;
	unsigned i, j;
	unsigned site, taxon;
	char ch;
	unsigned c;
	
	if (endCol > seqLen) {
		endCol = seqLen;
	}
	
	for (site = startCol; site < endCol; ++site) {
		// Recode a column of bases
		// WTJW 20/5/2003: Need to generalise to allow IUPAC codes.  This will
		// be slightly tricky, since e.g. these two transposed columns should map to
		// the same recoded bases:
		// A,[CGT],C
		// A,[CGT],G
		// If the worst comes to the worst, we could just declare that the
		// problem is "too hard" and any columns containing non-singleton base
		// sets will never be combined together into a single column.
		highest = 0;
		map[0] = map[1] = map[2] = map[3] = -1;
		keepSite = 1;
		ambiguous = 0;
		intersection = 0x0F;
		for (taxon = 0; taxon < numTaxa; ++taxon) {
			ch = toupper(fromBuf[taxon * seqLen + site]);
			c = GetMaskFromBase(ch);
			
			// Determine the list of bases at this site (usually 1, but can be up to 4)
			i = 0;
			if (c & BM_A) cList[i++] = 0;
			if (c & BM_C) cList[i++] = 1;
			if (c & BM_G) cList[i++] = 2;
			if (c & BM_T) cList[i++] = 3;
			
			if (i > 1) ambiguous = 1;
			
			// Discard sites containing gaps or ambiguous characters if told to
			if (i > 1 && discardMissAmbig) {
				keepSite = 0;
				++discardedSites;
				break;
			}
			
			// Determine the recoded bases that each base at this site should
			// map to.
			//HACK: this technique will always be correct, but not always optimal,
			// insofar as there may be sites which have the same pattern but
			// which are not condensed into the same column.
			for (j = 0; j < i; ++j) {
				if (map[cList[j]] == -1) {
					map[cList[j]] = highest++;
				}
			}
			
			// Recode this site
			c = 0;
			for (j = 0; j < i; ++j) {
				c |= 1 << map[cList[j]];
			}
			
			//HACK: What I think you need to do to handle all cases is:
			// 1) Determine bases from sets in such a way as to minimise the number of
			//    distinct bases in the column.  This minimises the number-of-
			//    distinct-bases-minus-one bound, which is a lower bound on
			//    the weight added by this column.  So a column produced like
			//    this will provide a minimum score that must be added by
			//    this column.  Call the number of distinct bases produced x.
			// 2) Go back to the start and determine bases from sets in such a
			//    way as to maximi
			
			intersection &= c;
			toBuf[numTaxa * newSeqLen + taxon] = (unsigned char) c;	// Transpose into toBuf
		}
		
		newSeqLen += keepSite;		// If this site should not be retained, it will be overwritten with the next
	}
	
	// Now set up the weights
	*pWeights = (unsigned *) malloc(newSeqLen * sizeof (unsigned));
	for (i = 0; i < newSeqLen; ++i) {
		(*pWeights)[i] = 1;
	}
	
	return newSeqLen;
}

//// Returns a copy of the input matrix which has been lexically sorted by base.
//// fromBuf and toBuf may be the same or different; the same for fromWeights
//// and toWeights.
//// 28/2/2005: Fixed bug: Need to sort weights too, OBVIOUSLY!
//void SortSites_b1t(unsigned char *fromBuf, unsigned char *toBuf, unsigned *fromWeights, unsigned *toWeights, unsigned seqLen, unsigned numTaxa) {
//	unsigned i;
//	unsigned char *augBuf = (unsigned char *) malloc(seqLen * (numTaxa + sizeof (unsigned)));
//	
//	// Build the temporary "augmented" array
//	for (i = 0; i < seqLen; ++i) {
//		*((unsigned *) (augBuf + (numTaxa + sizeof (unsigned)) * i)) = fromWeights[i];
//		memcpy(augBuf + (numTaxa + sizeof (unsigned)) * i + sizeof (unsigned), fromBuf + numTaxa * i, numTaxa);
//	}
//	
//	// Sort
//	_ColumnCompareWidth = numTaxa;					//HACK: using global variable
//	qsort(augBuf, seqLen, numTaxa + sizeof (unsigned), ColumnCompareIncludingWeights);
//	
//	// "Unpack" the augmented array
//	for (i = 0; i < seqLen; ++i) {
//		toWeights[i] = *((unsigned *) (augBuf + (numTaxa + sizeof (unsigned)) * i));
//		memcpy(toBuf + numTaxa * i, augBuf + (numTaxa + sizeof (unsigned)) * i + sizeof (unsigned), numTaxa);
//	}
//}

// Performs in-place compaction of a transposed character matrix.
// Assumes that the transposed matrix is already sorted by row.
// Returns the new number of sites.
int CombineIdenticalSites_b1t(unsigned char *buf, unsigned *weights, unsigned seqLen, unsigned numTaxa, unsigned effNumTaxa) {
	unsigned i;
	unsigned j = 0;
	
	for (i = 1; i < seqLen; ++i) {
		// Is site i equivalent to site j?
		if (!memcmp(buf + numTaxa * j, buf + numTaxa * i, effNumTaxa)) {
			// Sites are equivalent on first effNumTaxa taxa: combine weights
			weights[j] += weights[i];
		} else {
			// Sites are different
			++j;
			if (i != j) {
				memcpy(buf + numTaxa * j, buf + numTaxa * i, numTaxa);		// Copy whole site, even though we only looked at the first effNumTaxa taxa
				weights[j] = weights[i];
			}
		}
	}
	
	return j + 1;
}

// Performs in-place compaction of a transposed character matrix.
// Returns the new number of sites, and optionally records the number
// of constant sites removed and the total weight of non-PI sites removed.
// (pNonPiSites effectively counts sites, while pNonPiWeight counts the total
// weight added by those sites.)
int RemoveUninfSites_b1t(unsigned char *buf, unsigned *weights, unsigned seqLen, unsigned numTaxa, unsigned effNumTaxa, unsigned *pConstantWeight, unsigned *pNonPiSites, unsigned *pNonPiWeight) {
	unsigned i;
	unsigned j = 0;
	unsigned constantWeight = 0;
	unsigned nonPiWeight = 0;
	unsigned nonPiSites = 0;
	unsigned char intersection;
	unsigned keepSite;
	int weight;		// Return value from ParsimoniouslyUninformativeWeight_b1t() may be negative, indicating a PI site
	unsigned t;
	
	for (i = 0; i < seqLen; ++i) {
		keepSite = 1;
		
		// Is site i constant?
		// Iff the same base can be chosen from all base sets, the site is
		// effectively constant.  If the site is in fact a constant site,
		// this will of course still be true.
		intersection = 0x0F;
		for (t = 0; t < effNumTaxa && intersection; ++t) {
			intersection &= buf[numTaxa * i + t];
		}
		
		if (intersection) {
			// The site is constant
			keepSite = 0;
			constantWeight += weights[i];
		} else {
			// The site is not constant, but perhaps it is parsimoniously uninformative?
			weight = ParsimoniouslyUninformativeWeight_b1t(buf + numTaxa * i, effNumTaxa);	// I.e. number-of-distinct-bases-minus-one
			
			if (weight != -1) {
				keepSite = 0;
				nonPiWeight += weight * weights[i];
				nonPiSites += weights[i];
			}
		}
		
		if (keepSite) {
			if (i != j) {
				memcpy(buf + numTaxa * j, buf + numTaxa * i, numTaxa);		// Copy whole site, even though we only looked at the first effNumTaxa taxa
				weights[j] = weights[i];
			}
			++j;
		}
	}
	
	if (pConstantWeight) *pConstantWeight = constantWeight;
	if (pNonPiWeight) *pNonPiWeight = nonPiWeight;
	if (pNonPiSites) *pNonPiSites = nonPiSites;
	return j;
}

// Returns a copy of the input matrix which has been sorted by site weight.
// fromBuf and toBuf may be the same or different; the same for fromWeights
// and toWeights.
void SortSitesAndWeights_b1t(unsigned char *fromBuf, unsigned char *toBuf, unsigned *fromWeights, unsigned *toWeights, unsigned seqLen, unsigned numTaxa, int (*comparator)(const void *, const void *)) {
	unsigned i;
	unsigned char *augBuf = (unsigned char *) malloc(seqLen * (numTaxa + sizeof (unsigned)));
	
	// Build the temporary "augmented" array
	for (i = 0; i < seqLen; ++i) {
		*((unsigned *) (augBuf + (numTaxa + sizeof (unsigned)) * i)) = fromWeights[i];
		memcpy(augBuf + (numTaxa + sizeof (unsigned)) * i + sizeof (unsigned), fromBuf + numTaxa * i, numTaxa);
	}
	
	// Sort by the comparator function
	_ColumnCompareWidth = numTaxa;			//HACK: using global variable
	qsort(augBuf, seqLen, numTaxa + sizeof (unsigned), comparator);
	
	// "Unpack" the augmented array
	for (i = 0; i < seqLen; ++i) {
		toWeights[i] = *((unsigned *) (augBuf + (numTaxa + sizeof (unsigned)) * i));
		memcpy(toBuf + numTaxa * i, augBuf + (numTaxa + sizeof (unsigned)) * i + sizeof (unsigned), numTaxa);
	}
}

struct column {
	unsigned weight;				// Actually, this will hold log2(weight)
	unsigned char *data;
};

// Load the restBound[] array from a text file.  This file has the same format
// as that written out to the "restBound.txt" file, and that file itself can
// be used (the result will be written out to this file again afterwards!)
void LoadRestBound(struct TreeData *td) {
	FILE *f;
	static char buf[256];
	unsigned i, bound, line = 2;
	
	if (!(f = fopen(td->restBoundInputFname, "rb"))) {
		fprintf(stderr, "Unable to open file '%s' to retrieve restBound[] array, aborting.\n", td->restBoundInputFname);
		exit(1);
	}
	
	fgets(buf, sizeof buf, f);		// Skip header line
	while (fgets(buf, sizeof buf, f)) {
		if (sscanf(buf, "%u %u", &i, &bound) < 2) {
			fprintf(stderr, "Unable to parse line %u of '%s', aborting.\n", line, td->restBoundInputFname);
			fclose(f);
			exit(1);
		}
		
		td->restBound[i] = bound;		//HACK: Bounds checking?  What's that?
		++line;
	}
	
	fclose(f);
	fprintf(stderr, "Read %u lines from '%s'.\n", line, td->restBoundInputFname);
}

// This is the answer to EVERYTHING.  Uses aspects of PartitionColumns() and
// InitIncompatBoundRest(): Find the set of sites that are non-PI on the first
// n taxa; these sites have known cost x on ANY tree on these n taxa.  Then, use
// whatever lower-bounding technique you want on those sites to tell you what
// the minimum possible score of those sites on the full tree is (call this y).
// Now, you can always add (y-x) to the weight of any tree on those n taxa!
// We use a 1-connected partition bound approach for the lower bounding
// operation, as it seems to do much better than incompatibility-based (i.e.
// matchings or cliques).
//void InitKickAssBoundRest_BEFORE_REEXPANDING_HEAYV_COLUMNS(struct TreeData *td) {
void InitKickAssBoundRest(struct TreeData *td) {
	unsigned i, j, k, nonPiWeight, numSites, nonPiWeightTotal, finalWeight;
	unsigned partitionWeight, mstWeight;
	unsigned char *lowerBoundBuf;
	unsigned char *buf;
	unsigned *addScore;
	unsigned *weightsBuf;
	char fnameBuf[80];		//DEBUG
	FILE *nonPiFile;
	
	buf = (unsigned char *) malloc(td->numTaxa * 2);
	lowerBoundBuf = (unsigned char *) malloc(td->numTaxa * td->seqLen);
	addScore = (unsigned *) malloc(td->numTaxa * sizeof (unsigned));
	weightsBuf = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
	
	for (i = 3; i < td->numTaxa; ++i) {
		numSites = 0;
		nonPiWeightTotal = 0;
		for (j = 0; j < td->seqLen; ++j) {
			for (k = 0; k < td->numTaxa; ++k) {
//				buf[k] = GetBaseAt(k, j, td);
				buf[k] = td->charMat[k * td->seqLen + j];
			}
			
			// Is site j non-PI on a tree containing the first i taxa?
			if ((nonPiWeight = ParsimoniouslyUninformativeWeight_b1t(buf, i)) != -1) {
				// Yes: Add this site to the list
				for (k = 0; k < td->numTaxa; ++k) {
					lowerBoundBuf[k * td->seqLen + numSites] = buf[k];
				}
				weightsBuf[numSites] = td->weights[j];
				
				++numSites;
				nonPiWeightTotal += nonPiWeight;
			}
		}
		
		// Find a lower bound for the weight contributed by the sites in
		// the list when all taxa are on the tree.
		
		//HACK: ScorePartition() does a terrible job when there are many
		// sites, as it makes no attempt to subdivide.
		if (numSites > 1) {		//HACK: should be "> 0" but this seems to fail...
			partitionWeight = PartitionAndScoreColumns(lowerBoundBuf, numSites, td->seqLen, td->numTaxa, weightsBuf);		//HACK: partMask and weights params are dummies
			mstWeight = CalcMstWeight_b1(lowerBoundBuf, numSites, td->seqLen, td->numTaxa, weightsBuf, 1) / 2;
			if (partitionWeight > mstWeight) {
				finalWeight = partitionWeight;
				DBGPRINT4("Choosing partition bound of %u over MST bound of %u for %u taxa.\n", partitionWeight, mstWeight, i);
			} else {
				finalWeight = mstWeight;		//HACK
				DBGPRINT4("Choosing MST bound of %u over partition bound of %u for %u taxa.\n", mstWeight, partitionWeight, i);
			}
			
			if (finalWeight > nonPiWeightTotal) {
				addScore[i] = finalWeight - nonPiWeightTotal;
			} else {
				addScore[i] = 0;
			}
		} else {
			addScore[i] = 0;
		}
		
		//DEBUG: Write out a file containing the non-PI-down-to-taxon-(i+1) sites
		sprintf(fnameBuf, "nonpi%u.phy", i);
		nonPiFile = fopen(fnameBuf, "w");
		fprintf(nonPiFile, "%u %u\n", td->numTaxa, numSites);
		for (j = 0; j < td->numTaxa; ++j) {
			fprintf(nonPiFile, "T%02.2u       ", j + 1);
			PrintMaskSeq_b1(lowerBoundBuf + j * td->seqLen, numSites, nonPiFile);
			fputc('\n', nonPiFile);
		}
		fclose(nonPiFile);
		
		DBGPRINT6("InitKickAssBoundRest(): nonPiWeightTotal = %u, finalWeight = %u: Bound for %u taxa is %u (used %u sites).\n", nonPiWeightTotal, finalWeight, i, addScore[i], numSites);
	}
	
	// If it happens that a bound computed for n > m taxa is greater than that
	// computed for m taxa, adjust the bound for m taxa.
	for (i = td->numTaxa - 1; i >= 3; --i) {
		for (j = i - 1; j >= 3; --j) {
			if (addScore[i] > addScore[j]) {
				DBGPRINT5("InitKickAssBoundRest(): Bound for %u taxa (%u) is better than that for %u taxa (%u)!  Updating latter bound.\n", i, addScore[i], j, addScore[j]);
				addScore[j] = addScore[i];
			}
		}
	}
	
	for (i = 3; i < td->numTaxa; ++i) {
		if (addScore[i] > td->restBound[i - 1]) {
			DBGPRINT4("InitKickAssBoundRest(): Bound of %u for %u taxa is now %u better than before!\n", addScore[i], i, addScore[i] - td->restBound[i - 1]);
			td->restBound[i - 1] = addScore[i];
		}
	}
	
	free(buf);
	free(lowerBoundBuf);
	free(addScore);
	free(weightsBuf);
}

// This is the answer to EVERYTHING.  Uses aspects of PartitionColumns() and
// InitIncompatBoundRest(): Find the set of sites that are non-PI on the first
// n taxa; these sites have known cost x on ANY tree on these n taxa.  Then, use
// whatever lower-bounding technique you want on those sites to tell you what
// the minimum possible score of those sites on the full tree is (call this y).
// Now, you can always add (y-x) to the weight of any tree on those n taxa!
// We use a 1-connected partition bound approach for the lower bounding
// operation, as it seems to do much better than incompatibility-based (i.e.
// matchings or cliques).
void NEWInitKickAssBoundRest(struct TreeData *td) {
	unsigned i, j, k, nonPiWeight, numSites, nonPiWeightTotal, finalWeight;
	unsigned m, bufWidth;
	unsigned partitionWeight, mstWeight;
	unsigned char *lowerBoundBuf;
	unsigned char *buf;
	unsigned *addScore;
	unsigned *weightsBuf;
	char fnameBuf[80];		//DEBUG
	FILE *nonPiFile;
	
	buf = (unsigned char *) malloc(td->numTaxa * 2);
	bufWidth = 0;
	for (i = 0; i < td->seqLen; ++i) {
		bufWidth += td->weights[i];
	}
	lowerBoundBuf = (unsigned char *) malloc(td->numTaxa * bufWidth);
	addScore = (unsigned *) malloc(td->numTaxa * sizeof (unsigned));
	weightsBuf = (unsigned *) malloc(bufWidth * sizeof (unsigned));
	
	for (i = 3; i < td->numTaxa; ++i) {
		numSites = 0;
		nonPiWeightTotal = 0;
		for (j = 0; j < td->seqLen; ++j) {
			for (k = 0; k < td->numTaxa; ++k) {
//				buf[k] = GetBaseAt(k, j, td);
				buf[k] = td->charMat[k * td->seqLen + j];
			}
			
			// Is site j non-PI on a tree containing the first i taxa?
			if ((nonPiWeight = ParsimoniouslyUninformativeWeight_b1t(buf, i)) != -1) {
				// Now x weight-1 columns, where x is the weight of the site
				// Yes: Add this site to the list
				for (m = 0; m < td->weights[j]; ++m) {
					for (k = 0; k < td->numTaxa; ++k) {
						lowerBoundBuf[k * bufWidth + numSites] = buf[k];
					}
					
					weightsBuf[numSites] = 1;
					++numSites;
				}
				nonPiWeightTotal += nonPiWeight * td->weights[j];
			}
		}
		
		// Find a lower bound for the weight contributed by the sites in
		// the list when all taxa are on the tree.
		
		//HACK: ScorePartition() does a terrible job when there are many
		// sites, as it makes no attempt to subdivide.
		if (numSites > 1) {		//HACK: should be "> 0" but this seems to fail...
			partitionWeight = PartitionAndScoreColumns(lowerBoundBuf, numSites, bufWidth, td->numTaxa, weightsBuf);		//HACK: partMask and weights params are dummies
			mstWeight = CalcMstWeight_b1(lowerBoundBuf, numSites, bufWidth, td->numTaxa, weightsBuf, 1) / 2;
			if (partitionWeight > mstWeight) {
				finalWeight = partitionWeight;
				DBGPRINT4("Choosing partition bound of %u over MST bound of %u for %u taxa.\n", partitionWeight, mstWeight, i);
			} else {
				finalWeight = mstWeight;		//HACK
				DBGPRINT4("Choosing MST bound of %u over partition bound of %u for %u taxa.\n", mstWeight, partitionWeight, i);
			}
			
			if (finalWeight > nonPiWeightTotal) {
				addScore[i] = finalWeight - nonPiWeightTotal;
			} else {
				addScore[i] = 0;
			}
		} else {
			addScore[i] = 0;
		}
		
		//DEBUG: Write out a file containing the non-PI-down-to-taxon-(i+1) sites
		sprintf(fnameBuf, "nonpi%u.phy", i);
		nonPiFile = fopen(fnameBuf, "w");
		fprintf(nonPiFile, "%u %u\n", td->numTaxa, numSites);
		for (j = 0; j < td->numTaxa; ++j) {
			fprintf(nonPiFile, "T%02.2u       ", j + 1);
			PrintMaskSeq_b1(lowerBoundBuf + j * td->seqLen, numSites, nonPiFile);
			fputc('\n', nonPiFile);
		}
		fclose(nonPiFile);
		
		DBGPRINT6("InitKickAssBoundRest(): nonPiWeightTotal = %u, finalWeight = %u: Bound for %u taxa is %u (used %u sites).\n", nonPiWeightTotal, finalWeight, i, addScore[i], numSites);
	}
	
	// If it happens that a bound computed for n > m taxa is greater than that
	// computed for m taxa, adjust the bound for m taxa.
	for (i = td->numTaxa - 1; i >= 3; --i) {
		for (j = i - 1; j >= 3; --j) {
			if (addScore[i] > addScore[j]) {
				DBGPRINT5("InitKickAssBoundRest(): Bound for %u taxa (%u) is better than that for %u taxa (%u)!  Updating latter bound.\n", i, addScore[i], j, addScore[j]);
				addScore[j] = addScore[i];
			}
		}
	}
	
	for (i = 3; i < td->numTaxa; ++i) {
		if (addScore[i] > td->restBound[i - 1]) {
			DBGPRINT4("InitKickAssBoundRest(): Bound of %u for %u taxa is now %u better than before!\n", addScore[i], i, addScore[i] - td->restBound[i - 1]);
			td->restBound[i - 1] = addScore[i];
		}
	}
	
	free(buf);
	free(lowerBoundBuf);
	free(addScore);
	free(weightsBuf);
}

void InitKA2BoundRest(struct TreeData *td) {
	unsigned i, j, k, m, numSites, numSitesUsed, pos;
	unsigned *upperBoundWeights, *weightsBuf;
	int partitionWeight, mstWeight, finalWeight;
	unsigned char *lowerBoundBuf;
	unsigned char *buf;
	unsigned *addScore;
	char fnameBuf[80];		//DEBUG
	FILE *nonPiFile;
	
	buf = (unsigned char *) malloc(td->numTaxa * 2);
	addScore = (unsigned *) malloc(td->numTaxa * sizeof (unsigned));
	
	// We can include each column j weight[j] times!  If columns have
	// very high weights then the later random column choosing will
	// not be very efficient, however.
	numSites = 0;
	for (i = 0; i < td->seqLen; ++i) {
		numSites += td->weights[i];
	}
	
	lowerBoundBuf = (unsigned char *) malloc(td->numTaxa * numSites);		// Rectangular buffer: each row is a taxon, each site occupies 1 byte
	upperBoundWeights = (unsigned *) malloc(numSites * sizeof (unsigned));
	weightsBuf = (unsigned *) malloc(numSites * sizeof (unsigned));		//HACK: nothing is ever copied in here!
	
	for (i = 3; i < td->numTaxa; ++i) {
		pos = 0;
		for (j = 0; j < td->seqLen; ++j) {
			//HACK: following can probably be changed to "for (k = 0; k < j; ...)"
			for (k = 0; k < td->numTaxa; ++k) {
//				buf[k] = GetBaseAt(k, j, td);
				buf[k] = td->charMat[k * td->seqLen + j];
//				GetBaseFromMask(buf[k]);		//DEBUG
			}
			
			// Get an upper bound on site j on any MP tree containing the first i taxa
			{	//DEBUG
				unsigned nonPiWeight = ParsimoniouslyUninformativeWeight_b1t(buf, i);
				unsigned ubWeight = SiteWeightUpperBound_b1t(buf, i);
//				DBGPRINT5("upperBoundWeights[%u]=%u, nonPiWeight=%u (for tree with first %u taxa on it)\n", numSites, ubWeight, nonPiWeight, i);
				assert(nonPiWeight == -1 || nonPiWeight == ubWeight);
			}
			for (m = 0; m < td->weights[j]; ++m) {
				upperBoundWeights[pos] = SiteWeightUpperBound_b1t(buf, i);		//HACK: compute this once
				for (k = 0; k < td->numTaxa; ++k) {
					lowerBoundBuf[k * numSites + pos] = buf[k];
				}
				
				++pos;
			}
		}
		
		// Find a lower bound for the weight contributed by the sites in
		// the list when all taxa are on the tree.
		
		//HACK: ScorePartition() does a terrible job when there are many
		// sites, as it makes no attempt to subdivide.
		if (numSites > 1) {		//HACK: should be "> 0" but this seems to fail...
			// With KA2, not all sites may be used.  The returned weight has had the
			// upper bound weights subtracted off it already; we need to adjust
			// the mstWeight in the same way to make them comparable.
			// KA2PartitionAndScoreColumns() changes the upperBoundWeights[] array.
			// IMPORTANT!  Make a copy of the weights array and use that.  Before,
			// we were allowing KA2PartitionAndScoreColumns() to adjust td->weights[]
			// without allowing it to adjust td->charMat[]!
			for (k = 0; k < numSites; ++k) {			// Now that we have reexpanded heavy columns
				weightsBuf[k] = 1;
			}
			partitionWeight = KA2PartitionAndScoreColumns(lowerBoundBuf, numSites, numSites, td->numTaxa, weightsBuf, upperBoundWeights, &numSitesUsed);		//HACK: partMask and weights params are dummies
			//DEBUG: see whether MST weights are screwing us up
			//CHECK2009
//			mstWeight = CalcMstWeight_b1(lowerBoundBuf, numSitesUsed, td->seqLen, td->numTaxa, td->weights, 1) / 2;		//HACK: partMask and weights params are dummies
//			for (j = 0; j < numSitesUsed; ++j) {
//				mstWeight -= upperBoundWeights[j];
//			}
//			if (partitionWeight > mstWeight) {
				finalWeight = partitionWeight;
//				DBGPRINT4("Choosing partition bound of %u over MST bound of %u for %u taxa.\n", partitionWeight, mstWeight, i);
//			} else {
//				finalWeight = mstWeight;		//HACK
//				DBGPRINT4("Choosing MST bound of %u over partition bound of %u for %u taxa.\n", mstWeight, partitionWeight, i);
//			}
			
			if (finalWeight > 0) {
				addScore[i] = finalWeight;
			} else {
				addScore[i] = 0;
			}
		} else {
			addScore[i] = 0;
		}
		
		//DEBUG: Write out a file containing the non-PI-down-to-taxon-(i+1) sites
		sprintf(fnameBuf, "nonpi%u.phy", i);
		nonPiFile = fopen(fnameBuf, "w");
		fprintf(nonPiFile, "%u %u\n", td->numTaxa, numSites);
		for (j = 0; j < td->numTaxa; ++j) {
			fprintf(nonPiFile, "T%02.2u       ", j + 1);
			PrintMaskSeq_b1(lowerBoundBuf + j * td->seqLen, numSites, nonPiFile);		// Remember: only the first numSitesUsed of these sites are included
			fputc('\n', nonPiFile);
		}
		fclose(nonPiFile);
		
		DBGPRINT5("InitKA2BoundRest(): finalWeight = %u: Bound for %u taxa is %u (used %u sites).\n", finalWeight, i, addScore[i], numSitesUsed);
	}
	
	// If it happens that a bound computed for n > m taxa is greater than that
	// computed for m taxa, adjust the bound for m taxa.
	for (i = td->numTaxa - 1; i >= 3; --i) {
		for (j = i - 1; j >= 3; --j) {
			if (addScore[i] > addScore[j]) {
				DBGPRINT5("InitKickAssBoundRest(): Bound for %u taxa (%u) is better than that for %u taxa (%u)!  Updating latter bound.\n", i, addScore[i], j, addScore[j]);
				addScore[j] = addScore[i];
			}
		}
	}
	
	for (i = 3; i < td->numTaxa; ++i) {
		if (addScore[i] > td->restBound[i - 1]) {
			DBGPRINT4("InitKickAssBoundRest(): Bound of %u for %u taxa is now %u better than before!\n", addScore[i], i, addScore[i] - td->restBound[i - 1]);
			td->restBound[i - 1] = addScore[i];
		}
	}
	
	free(buf);
	free(lowerBoundBuf);
	free(addScore);
	free(upperBoundWeights);
	free(weightsBuf);
}

// Must be called pre-SQUEEZEBASES.
//#define PROB_KEEP_WORSE RAND_MAX / 100		// Not yet implemented properly
#define NREPS 500
//#define PARTITION_SIZE 8						// Fixed size of partition
#define PARTITION_SIZE 4						// Fixed size of partition
//HACK
//#ifdef SQUEEZEBASES
//#if PARTITION_SIZE != 8
//#error PARTITION_SIZE must be 8 if using SQUEEZEBASES.
//#endif
//#else	// not SQUEEZEBASES
//#if PARTITION_SIZE != 4
//#error PARTITION_SIZE must be 4 if not using SQUEEZEBASES.
//#endif
//#endif	// SQUEEZEBASES
#define PROB_KEEP_WORSE 0						// Not yet implemented properly
//#define PROB_KEEP_WORSE RAND_MAX / 100			// Not yet implemented properly
//#define DELTA_PROB_KEEP_WORSE RAND_MAX / 300
#define DELTA_PROB_KEEP_WORSE 0
#define NUM_MUTATIONS 2							// Number of site swaps to perform in between scoring operations
#define REALLY_GREEDY 1							// Only considers score of current partition, ignores score of remaining "site pool"
#define GREEDINESS 4							// Multiplier for weight of current partition (only relevant when REALLY_GREEDY not #defined)
#define FINAL_SWAPS 500						// Number of site swaps to try after initial processing
unsigned KA2PartitionAndScoreColumns(unsigned char *data, unsigned dataWidth, unsigned width, unsigned height, unsigned *weights, unsigned *upperBoundWeights, unsigned *pNumSitesUsed) {
	unsigned numParts, partWidth = PARTITION_SIZE;
	unsigned char *partBuf, *bestPartBuf;
	int score, newScore, bestScore;
	unsigned *weightsBuf, *bestWeightsBuf, *upperBoundWeightsBuf, *bestUpperBoundWeightsBuf, *partWidths;
	unsigned i, j, nReps = NREPS;
	int total = 0;
	unsigned halvesScore, bestPartMask;
	int restScore, newRestScore, bestRestScore;
	unsigned offsetA, offsetB;
	unsigned partA, partB;
	unsigned partAWidth, partBWidth;
	int *partScores;
	int newPartAScore, newPartBScore;
	int checkTotal1 = 0, checkTotal2 = 0;
	unsigned probKeepWorse = PROB_KEEP_WORSE;
	unsigned pos = 0;
	//DEBUG
	char partFileName[80];
	FILE *partFile;
	
	partBuf = (unsigned char *) malloc(width * height);
	weightsBuf = (unsigned *) malloc(width * sizeof (unsigned));
	bestWeightsBuf = (unsigned *) malloc(width * sizeof (unsigned));
	partScores = (int *) malloc(width * sizeof (int));
	bestPartBuf = (unsigned char *) malloc(width * height);		//HACK
	upperBoundWeightsBuf = (unsigned *) malloc(width * sizeof (unsigned));
	bestUpperBoundWeightsBuf = (unsigned *) malloc(width * sizeof (unsigned));
	partWidths = (unsigned *) malloc(width * sizeof (unsigned));
	
#ifdef DEBUG
	for (i = 0; i < height; ++i) {
		PrintMaskSeq_b1(data + i * width, dataWidth, dbgfile);
		DBGPRINT1("\n");
	}
#endif	//DEBUG

	bestScore = 1;		// Necessary to force "if (bestScore <= 0)" not to run on first iteration
	pos = 0;
	numParts = 0;
	while (1) {
		// Adjust the partition width.  This is initialised to PARTITION_SIZE
		// and will stay at that unless either (a) we are at the end, or (b) the
		// previous score was negative.
		// (a): Set the width to the remaining width.
		// (b): Decrement the width by one and try again.
		if (bestScore <= 0) {
			if (!partWidth) {
				break;					// Can't shrink partitions any further
			}
			pos -= partWidth;			// We will reprocess the last partition
			total -= bestScore;			// Undo adding of negative score
			--partWidth;
			--numParts;
		}
		if (pos + partWidth > dataWidth) {
			partWidth = dataWidth - pos;
		}
		if (!partWidth) break;				// Loop termination condition
		partWidths[numParts] = partWidth;	// Important!  Now that partition widths can vary, we need to record each partition's width
		
		DBGPRINT3("Starting partition processing loop. pos=%u, partWidth=%u.\n", pos, partWidth);
		
		// Evaluate lower-bound score of this partition
		bestScore = score = 0;
#ifdef REALLY_GREEDY
		bestRestScore = 0;
#else	// not REALLY_GREEDY
#error Not implemented!	//HACK
#endif	// not REALLY_GREEDY
		if (pos + partWidth < dataWidth) {
			for (i = 0; i < nReps; ++i) {
				// Create a copy of the current state
				memcpy(partBuf, data, width * height);		//HACK: some unnecessary copying
				memcpy(weightsBuf, weights, dataWidth * sizeof (unsigned));
				memcpy(upperBoundWeightsBuf, upperBoundWeights, dataWidth * sizeof (unsigned));
				
				// Randomly alter this new version
				for (j = 0; j < NUM_MUTATIONS; ++j) {
					KA2MutatePartition(partBuf + pos, partWidth, dataWidth - pos, width, height, weightsBuf + pos, NULL, NULL, upperBoundWeightsBuf + pos);
				}
				
				// Evaluate new partition.  This actually uses the same old ExhaustiveScorePartition() as for KICKASSBOUND,
				// since the site score upper bounds are independent of the way these scores are calculated.
				newScore = ExhaustiveScorePartition(partBuf + pos, partWidth, (1 << partWidth) - 1, width, height, weightsBuf + pos, &bestPartMask, 0);
				// Subtract out site upper bound scores
				for (j = 0; j < partWidth; ++j) {
					newScore -= upperBoundWeightsBuf[pos + j];
				}
#ifdef REALLY_GREEDY
				newRestScore = 0;
#else	// not REALLY_GREEDY
				// Can't afford to do exhaustive search on potentially very wide partition remaining!
#error Not implemented!	//HACK
#endif	// REALLY_GREEDY
				
				if (GREEDINESS * newScore + newRestScore > GREEDINESS * bestScore + bestRestScore) {
					DBGPRINT6("New score of %u (restScore=%u, combined=%u) and partMask of %x beats old score of %u!\n", newScore, restScore, newScore + restScore, bestPartMask, score);
					memcpy(bestPartBuf, partBuf, width * height);		//HACK: some unnecessary copying
					memcpy(bestWeightsBuf, weightsBuf, dataWidth * sizeof (unsigned));		//HACK: some unnecessary copying
					memcpy(bestUpperBoundWeightsBuf, upperBoundWeightsBuf, dataWidth * sizeof (unsigned));
					memcpy(data, partBuf, width * height);		//HACK: some unnecessary copying
					memcpy(weights, weightsBuf, dataWidth * sizeof (unsigned));
					memcpy(upperBoundWeights, upperBoundWeightsBuf, dataWidth * sizeof (unsigned));
					bestScore = score = newScore;
					bestRestScore = restScore = newRestScore;
					probKeepWorse = PROB_KEEP_WORSE;
				} else if (GREEDINESS * newScore + newRestScore > GREEDINESS * score + restScore || rand() < probKeepWorse) {
					DBGPRINT4("The new score of %u and partMask of %x is WORSE THAN (or equal to) the old score of %u, but let's keep it anyway.\n", newScore, bestPartMask, score);
					memcpy(data, partBuf, width * height);		//HACK: some unnecessary copying
					memcpy(weights, weightsBuf, dataWidth * sizeof (unsigned));
					memcpy(upperBoundWeights, upperBoundWeightsBuf, dataWidth * sizeof (unsigned));
					score = newScore;
					restScore = newRestScore;
					if (probKeepWorse > DELTA_PROB_KEEP_WORSE) {
						probKeepWorse -= DELTA_PROB_KEEP_WORSE;		// So the probability of making a non-optimal choice decreases the more such choices we make
					} else {
						probKeepWorse = 0;		// So the probability of making a non-optimal choice decreases the more such choices we make
					}
				}
			}
		} else {
			// This is the final partition.  We don't (cannot!) perform any site swaps.
			
			// Evaluate new partition.  This actually uses the same old ExhaustiveScorePartition() as for KICKASSBOUND,
			// since the site score upper bounds are independent of the way these scores are calculated.
			newScore = ExhaustiveScorePartition(data + pos, partWidth, (1 << partWidth) - 1, width, height, weights + pos, &bestPartMask, 0);		// IMPORTANT to use weights instead of weightsBuf!
			// Subtract out site upper bound scores
			for (j = 0; j < partWidth; ++j) {
				newScore -= upperBoundWeights[pos + j];
			}
			bestScore = score = newScore;
			// Need to populate the best...[] arrays so that their contents can be copied back to the "main" ones.
			if (bestScore > 0) {
				memcpy(bestPartBuf, data, width * height);		//HACK: some unnecessary copying
				memcpy(bestWeightsBuf, weights, dataWidth * sizeof (unsigned));		//HACK: some unnecessary copying
				memcpy(bestUpperBoundWeightsBuf, upperBoundWeights, dataWidth * sizeof (unsigned));
			}
		}
		
		DBGPRINT4("The score of the best partition for partition #%u is %d.  This has width %u.\n", numParts, bestScore, partWidth);
		if (bestScore > 0) {		// Otherwise bestPartBuf etc. will not have been filled out
			memcpy(data, bestPartBuf, width * height);		//HACK: some unnecessary copying
			memcpy(weights, bestWeightsBuf, dataWidth * sizeof (unsigned));
			memcpy(upperBoundWeights, bestUpperBoundWeightsBuf, dataWidth * sizeof (unsigned));
		}
		
//		//DEBUG
//#ifdef DEBUG
//		for (i = 0; i < height; ++i) {
//			PrintMaskSeq_b1(data + i * width, dataWidth, dbgfile);
//			DBGPRINT1("\n");
//		}
//#endif	//DEBUG
		
		//DEBUG
		sprintf(partFileName, "NEWpartition%u_partscore=%u.phy", numParts, bestScore);
		partFile = fopen(partFileName, "w");
		fprintf(partFile, "%u %u\n", height, partWidth);
		for (i = 0; i < height; ++i) {
			fprintf(partFile, "T%02u       ", i + 1);
			for (j = 0; j < partWidth; ++j) {
				fputc(GetBaseFromMask(data[i * width + pos + j]), partFile);
			}
			fputc('\n', partFile);
		}
		fclose(partFile);
		
		partScores[numParts] = bestScore;
		total += bestScore;
		pos += partWidth;
		++numParts;
	}
	
	pos = 0;
	for (i = 0; i < numParts; ++i) {
//		checkTotal1 += ExhaustiveScorePartition(data + pos, partWidths[i], (1 << partWidths[i]) - 1, width, height, weights + pos, NULL, 0);
		score = ExhaustiveScorePartition(data + pos, partWidths[i], (1 << partWidths[i]) - 1, width, height, weights + pos, NULL, 0);
		DBGPRINT5("Starting calculation of checkTotal1 on partition #%u.  ExhaustiveScorePartition(%u, %u, ...)=%u.\n", i, pos, partWidths[i], score);
		checkTotal1 += score;
		for (j = 0; j < partWidths[i]; ++j) {
//			checkTotal1 -= upperBoundWeightsBuf[pos + j];
			DBGPRINT3("Subtracting upperBoundWeight[%u]=%u from checkTotal1...\n", pos + j, upperBoundWeights[pos + j]);
			checkTotal1 -= upperBoundWeights[pos + j];
		}
		DBGPRINT3("Running total for checkTotal1 = %u after processing partition #%u.\n", checkTotal1, i);
		checkTotal2 += partScores[i];
		pos += partWidths[i];
	}
	DBGPRINT3("checkTotal1 = %u, checkTotal2 = %u.\n", checkTotal1, checkTotal2);
	assert(checkTotal1 == checkTotal2);		//DEBUG
	
	DBGPRINT4("Completed KA2PartitionColumns() with total bound of %u.  %u sites were used (from the total of %u sites).\n", total, pos, dataWidth);
	
	free(partBuf);
	free(bestPartBuf);
	free(weightsBuf);
	free(bestWeightsBuf);
	free(partScores);
	free(upperBoundWeightsBuf);
	free(bestUpperBoundWeightsBuf);
	free(partWidths);
	
	*pNumSitesUsed = pos;
	return total;
}






unsigned PartitionAndScoreColumns(unsigned char *data, unsigned dataWidth, unsigned width, unsigned height, unsigned *weights) {
	unsigned part, numParts, thisPartWidth, partWidth = PARTITION_SIZE;
	unsigned char *partBuf, *bestPartBuf;
	unsigned score, newScore, bestScore;
	unsigned *weightsBuf, *bestWeightsBuf;
	unsigned i, j, nReps = NREPS;
	unsigned total = 0;
	unsigned halvesScore, bestPartMask, restScore, newRestScore, bestRestScore;
	unsigned offsetA, offsetB;
	unsigned partA, partB;
	unsigned partAWidth, partBWidth;
	unsigned *partScores;
	unsigned newPartAScore, newPartBScore;
	unsigned checkTotal1 = 0, checkTotal2 = 0;
	unsigned probKeepWorse = PROB_KEEP_WORSE;
	//DEBUG
	char partFileName[80];
	FILE *partFile;
	
	numParts = (dataWidth - 1) / partWidth + 1;		//TODO: If numParts can be 0, this needs extra handling.
	partBuf = (unsigned char *) malloc(width * height);
	weightsBuf = (unsigned *) malloc(width * sizeof (unsigned));
	bestWeightsBuf = (unsigned *) malloc(width * sizeof (unsigned));
	partScores = (unsigned *) malloc(width * sizeof (unsigned));
	bestPartBuf = (unsigned char *) malloc(width * height);		//HACK
	
#ifdef DEBUG
	for (i = 0; i < height; ++i) {
		PrintMaskSeq_b1(data + i * width, dataWidth, dbgfile);
		DBGPRINT1("\n");
	}
#endif	//DEBUG
	
	for (part = 0; part < numParts - 1; ++part) {
		thisPartWidth = partWidth;
		
		// Evaluate lower-bound score of this partition
		bestScore = score = 0;
#ifdef REALLY_GREEDY
		bestRestScore = 0;
#else	// not REALLY_GREEDY
		bestRestScore = ScorePartition(data + (part + 1) * partWidth, dataWidth - (part + 1) * partWidth, (1 << PARTITION_SIZE) - 1, width, height, weights + part * partWidth);
#endif	// not REALLY_GREEDY
		for (i = 0; i < nReps; ++i) {
			// Create a copy of the current state
			memcpy(partBuf, data, width * height);		//HACK: some unnecessary copying
			memcpy(weightsBuf, weights, dataWidth * sizeof (unsigned));
			
			// Randomly alter this new version
			for (j = 0; j < NUM_MUTATIONS; ++j) {
				MutatePartition(partBuf + part * partWidth, thisPartWidth, dataWidth - part * partWidth, width, height, weightsBuf + part * partWidth, NULL, NULL);		// WTJW 11/6/2004
			}
			
			// Evaluate new partition
			newScore = ExhaustiveScorePartition(partBuf + part * partWidth, thisPartWidth, (1 << PARTITION_SIZE) - 1, width, height, weights + part * partWidth, &bestPartMask, 0);
#ifdef REALLY_GREEDY
			newRestScore = 0;
#else	// not REALLY_GREEDY
			// Can't afford to do exhaustive search on potentially very wide partition remaining!
			newRestScore = ScorePartition(data + (part + 1) * partWidth, dataWidth - (part + 1) * partWidth, (1 << PARTITION_SIZE) - 1, width, height, weights + part * partWidth);
#endif	// not REALLY_GREEDY
			
			if (GREEDINESS * newScore + newRestScore > GREEDINESS * bestScore + bestRestScore) {
				DBGPRINT6("New score of %u (restScore=%u, combined=%u) and partMask of %x beats old score of %u!\n", newScore, restScore, newScore + restScore, bestPartMask, score);
				memcpy(bestPartBuf, partBuf, width * height);		//HACK: some unnecessary copying
				memcpy(bestWeightsBuf, weightsBuf, dataWidth * sizeof (unsigned));		//HACK: some unnecessary copying
				memcpy(data, partBuf, width * height);		//HACK: some unnecessary copying
				memcpy(weights, weightsBuf, dataWidth * sizeof (unsigned));
				bestScore = score = newScore;
				bestRestScore = restScore = newRestScore;
				probKeepWorse = PROB_KEEP_WORSE;
			} else if (GREEDINESS * newScore + newRestScore > GREEDINESS * score + restScore || rand() < probKeepWorse) {
				DBGPRINT4("The new score of %u and partMask of %x is WORSE THAN (or equal to) the old score of %u, but let's keep it anyway.\n", newScore, bestPartMask, score);
				memcpy(data, partBuf, width * height);		//HACK: some unnecessary copying
				memcpy(weights, weightsBuf, dataWidth * sizeof (unsigned));
				score = newScore;
				restScore = newRestScore;
				if (probKeepWorse > DELTA_PROB_KEEP_WORSE) {
					probKeepWorse -= DELTA_PROB_KEEP_WORSE;		// So the probability of making a non-optimal choice decreases the more such choices we make
				} else {
					probKeepWorse = 0;		// So the probability of making a non-optimal choice decreases the more such choices we make
				}
			}
		}
		
		DBGPRINT3("The score of the best partition for partition #%u is %u.\n", part, bestScore);
		memcpy(data, bestPartBuf, width * height);		//HACK: some unnecessary copying
		memcpy(weights, bestWeightsBuf, dataWidth * sizeof (unsigned));
		
		//DEBUG
		sprintf(partFileName, "NEWpartition%u_partscore=%u.phy", part, bestScore);
		partFile = fopen(partFileName, "w");
		fprintf(partFile, "%u %u\n", height, thisPartWidth);
		for (i = 0; i < height; ++i) {
			fprintf(partFile, "T%02u       ", i + 1);
			for (j = 0; j < thisPartWidth; ++j) {
				fputc(GetBaseFromMask(data[i * width + part * partWidth + j]), partFile);
			}
			fputc('\n', partFile);
		}
		fclose(partFile);
		
		partScores[part] = bestScore;
		total += bestScore;
	}
	
	// Evaluate the final (possibly narrower) partition
	bestScore = score = ExhaustiveScorePartition(data + part * partWidth, dataWidth - (numParts - 1) * partWidth, (1 << PARTITION_SIZE) - 1, width, height, weights + part * partWidth, NULL, 0);
	DBGPRINT3("The score of the (only) partition for partition #%u is %u.\n", part, score);
	partScores[part] = bestScore;
	total += bestScore;
	
	for (part = 0; part < numParts; ++part) {
		checkTotal1 += ExhaustiveScorePartition(data + part * partWidth, (part == numParts - 1) ? (dataWidth - (numParts - 1) * partWidth) : partWidth, (1 << PARTITION_SIZE) - 1, width, height, weights + part * partWidth, NULL, 0);
		checkTotal2 += partScores[part];
	}
	DBGPRINT3("checkTotal1 = %u, checkTotal2 = %u.\n", checkTotal1, checkTotal2);
	assert(checkTotal1 == checkTotal2);		//DEBUG
	
	DBGPRINT2("Completed PartitionColumns() with total bound of %u.\n", total);
	
	free(partBuf);
	free(bestPartBuf);
	free(weightsBuf);
	free(bestWeightsBuf);
	free(partScores);
	
	return total;
}

//CHECK2009
//#define EXHAUSTIVESCOREPARTITION_RECURSIVE 1
unsigned ExhaustiveScorePartition(unsigned char *part, unsigned partWidth, unsigned partMask, unsigned width, unsigned height, unsigned *weights, unsigned *pBestPartMask, unsigned depth) {
	unsigned j, score, bestScore = 0, bestPartMask;
	unsigned onBits = 0;
	unsigned indices[32];
	unsigned mask, notMask;
	unsigned subdivideScore;
	
	// Figure out which bits are on in the partMask and how many there are in total
	for (j = 0; j < partWidth; ++j) {
		if ((partMask >> j) & 1) {
			indices[onBits++] = j;
		}
	}
	
//	DBGPRINT6("%*sExhaustiveScorePartition(): depth = %u, partMask = %x, onBits = %u\n", depth * 2, "", depth, partMask, onBits);
	
	// Consider all 2^n possible bipartitions of this partition
	for (j = 0; j < (1 << onBits); ++j) {
		mask = TurnOnBits(j, onBits, indices);
		notMask = TurnOnBits(~j, onBits, indices);
		
		score = ScorePartition(part, partWidth, mask, width, height, weights) +
			ScorePartition(part, partWidth, notMask, width, height, weights);
		if (score > bestScore) {
			bestScore = score;
			bestPartMask = mask;
		}
#ifdef EXHAUSTIVESCOREPARTITION_RECURSIVE
		if (onBits > 1 && j != 0 && j != (1 << onBits) - 1) {
//			DBGPRINT5("%*sExhaustiveScorePartition() recursing.  mask=%x, notMask=%x.\n", depth * 2, "", mask, notMask);
			subdivideScore =
				ExhaustiveScorePartition(part, partWidth, mask, width, height, weights, NULL, depth + 1) +
				ExhaustiveScorePartition(part, partWidth, notMask, width, height, weights, NULL, depth + 1);
			
			if (subdivideScore > bestScore) {
				bestScore = subdivideScore;
			}
		}
#endif	// EXHAUSTIVESCOREPARTITION_RECURSIVE
	}
	
	if (pBestPartMask) {
//		DBGPRINT5("pBestPartMask (at addr %p, current val=%u) is not NULL, so I will assign %u to it.  bestScore=%u, by the way.\n", pBestPartMask, *pBestPartMask, bestPartMask, bestScore);
		*pBestPartMask = bestPartMask;
	}
	
//	DBGPRINT4("%*sbestScore = %u.\n", depth * 2, "", bestScore);
	return bestScore;
}

// Use a nifty log2(nBits) algorithm to count the number of on bits in a 32-bit integer.
unsigned CountOnBits(unsigned x) {
	x = (0x55555555 & x) + (0x55555555 & (x >> 1));
	x = (0x33333333 & x) + (0x33333333 & (x >> 2));
	x = (0x0f0f0f0f & x) + (0x0f0f0f0f & (x >> 4));
	x = (0x00ff00ff & x) + (0x00ff00ff & (x >> 8));
	x = (0x0000ffff & x) + (0x0000ffff & (x >> 16));
	return x;
}

// Turn on bit indices[i] if the ith bit in bits is on.  Affects n bits altogether.
unsigned TurnOnBits(unsigned bits, unsigned n, unsigned *indices) {
	unsigned i, x = 0;
	
	for (i = 0; i < n; ++i) {
		x |= ((bits >> i) & 1) << indices[i];
	}
	
	return x;
}

// WTJW 24/3/2005: Revamped this.  Now use a Kruskal-style algo to figure out
// 1-connectedness.
//HACK: what the hell to do with them site weights?  Currently assumes every site has weight 1 (i.e. ALLWEIGHT1 is #defined)
unsigned ScorePartition(unsigned char *part, unsigned partWidth, unsigned partMask, unsigned width, unsigned height, unsigned *weights) {
	unsigned *counted;
	unsigned nDistinct = 1;
	unsigned i, j, k, differences;
	unsigned **edgeList;		// 2D array of edges, indexed by distance
	unsigned *edgeCount;
	unsigned curLen, nEdges, maxDifferences = 0;
	unsigned bI, bJ, lastEdgeWeight, temp;
	unsigned *baseOf;			// Records lowest-numbered node in component that each node belongs to.
	
//	DBGPRINT5("Entering ScorePartition(partWidth=%d, partMask=%x, width=%d, height=%d).\n", partWidth, partMask, width, height);
	
	counted = (unsigned *) malloc(height * sizeof (unsigned));
	memset(counted, 0, height * sizeof (unsigned));
	
	// There can be at most partWidth differences between any two sequences.
	edgeList = (unsigned **) malloc((partWidth + 1) * sizeof (unsigned *));
	edgeCount = (unsigned *) malloc((partWidth + 1) * sizeof (unsigned));
	for (i = 0; i < partWidth + 1; ++i) {
		// There could be as many as (height * (height - 1) / 2) edges in any
		// given list, and each edge requires 2 unsigned integers.
		edgeList[i] = (unsigned *) malloc(height * (height - 1) * sizeof (unsigned));
		edgeCount[i] = 0;
	}
	
	for (i = 1; i < height; ++i) {
		if (counted[i]) continue;
		
		for (j = 0; j < i; ++j) {
			if (counted[j]) continue;
			differences = 0;
			for (k = 0; k < partWidth; ++k) {
				if ((partMask >> k) & 1) {
					// Use a liberal definition of "different", since we want a lower bound
					if (!(part[i * width + k] & part[j * width + k])) {
						++differences;
					}
				}
			}
			
			if (!differences) {
				counted[i] = 1;				// This seq is identical to the one being examined
				break;				// Would be nice if we could say "Continue outer loop" like in Perl
			}
			
			// Add this edge to the "sorted list"
			edgeList[differences][edgeCount[differences] * 2] = i;
			edgeList[differences][edgeCount[differences] * 2 + 1] = j;
			++edgeCount[differences];
			
			if (differences > maxDifferences) {
				maxDifferences = differences;
			}
		}
		
		if (!counted[i]) {
			++nDistinct;
		}
	}
	
	if (nDistinct == 1) {
		for (i = 0; i < partWidth + 1; ++i) {
			free(edgeList[i]);
		}
		free(edgeCount);
		free(edgeList);
		free(counted);
		return 0;
	}
	
	// Now perform Kruskal's algorithm, pulling out one edge at a time in
	// order of increasing length.  If the entire graph becomes connected while
	// the edge length is still 1, then it is 1-connected, otherwise it is not.
	// In the latter case, we can look through the minimum distances between
	// each component and some other component, choose the largest, and add this to the MST
	// length.  The reason being that at least that many "metres of road"
	// must be laid down to connect that component to some other component on
	// any SMT, and in the worst case (for us calculating lower bounds) such a
	// road may enable all components to be connected.  This turns out to be
	// the weight of the last edge added to the MST.
	
	baseOf = (unsigned *) malloc(height * sizeof (unsigned));
	for (i = 0; i < height; ++i) {
		baseOf[i] = i;
	}
	
	curLen = 1;
	nEdges = 0;
	while (nEdges < nDistinct - 1) {
		// Choose next edge
		if (!edgeCount[curLen]) {
			++curLen;
			continue;
		}
		
		--edgeCount[curLen];
		i = edgeList[curLen][edgeCount[curLen] * 2];
		j = edgeList[curLen][edgeCount[curLen] * 2 + 1];
		
		// Make sure this edge would not create a cycle, by finding the lowest-
		// numbered node contained in this component.
		bI = i;
		while (baseOf[bI] != bI) {
			bI = baseOf[bI];
		}
		
		bJ = j;
		while (baseOf[bJ] != bJ) {
			bJ = baseOf[bJ];
		}
		
		if (bI == bJ) {
			continue;
		}
		
		// We know that the edge joins separate components, so join them now.
		if (bI < bJ) {
			// Compress paths to speed up future connectivity tests
			bJ = j;
			while (baseOf[bJ] != bJ) {
				temp = baseOf[bJ];
				baseOf[bJ] = bI;
				bJ = temp;
			}
			baseOf[bJ] = bI;
		} else {
			// Compress paths to speed up future connectivity tests
			bI = i;
			while (baseOf[bI] != bI) {
				temp = baseOf[bI];
				baseOf[bI] = bJ;
				bI = temp;
			}
			baseOf[bI] = bJ;
		}
		
		++nEdges;
		lastEdgeWeight = curLen;
	}
	
	free(baseOf);
	
	for (i = 0; i < partWidth + 1; ++i) {
		free(edgeList[i]);
	}
	
	free(edgeCount);
	free(edgeList);
	free(counted);
	
	return nDistinct - 2 + lastEdgeWeight;
}

void MutatePartition(unsigned char *part, unsigned partWidth, unsigned remainingWidth, unsigned width, unsigned height, unsigned *weights, unsigned *pRIn, unsigned *pROut) {
	unsigned rIn, rOut, tempWeight, i;
	unsigned char tempBase;
	
	// Pick a site inside the partition and a site outside of the partition at
	// random, and swap them
	rIn = random(partWidth);
	rOut = random(remainingWidth - partWidth) + partWidth;
	
	for (i = 0; i < height; ++i) {
		tempBase = part[width * i + rIn];
		part[width * i + rIn] = part[width * i + rOut];
		part[width * i + rOut] = tempBase;
	}
	
	tempWeight = weights[rIn];
	weights[rIn] = weights[rOut];
	weights[rOut] = tempWeight;
	
	// Pass back the positions of the sites we swapped
	if (pRIn) *pRIn = rIn;
	if (pROut) *pROut = rOut;
}

void KA2MutatePartition(unsigned char *part, unsigned partWidth, unsigned remainingWidth, unsigned width, unsigned height, unsigned *weights, unsigned *pRIn, unsigned *pROut, unsigned *upperBoundWeights) {
	unsigned rIn, rOut, tempWeight, tempUpperBoundWeight, i;
	unsigned char tempBase;
	
	// Pick a site inside the partition and a site outside of the partition at
	// random, and swap them
	rIn = random(partWidth);
	rOut = random(remainingWidth - partWidth) + partWidth;
	
	for (i = 0; i < height; ++i) {
		tempBase = part[width * i + rIn];
		part[width * i + rIn] = part[width * i + rOut];
		part[width * i + rOut] = tempBase;
	}
	
	tempWeight = weights[rIn];
	weights[rIn] = weights[rOut];
	weights[rOut] = tempWeight;
	tempUpperBoundWeight = upperBoundWeights[rIn];
	upperBoundWeights[rIn] = upperBoundWeights[rOut];
	upperBoundWeights[rOut] = tempUpperBoundWeight;
	
	// Pass back the positions of the sites we swapped
	if (pRIn) *pRIn = rIn;
	if (pROut) *pROut = rOut;
}

void SwapSites(unsigned char *data, unsigned a, unsigned b, unsigned width, unsigned height, unsigned *weights) {
	unsigned tempWeight, i;
	unsigned char tempBase;
	
	for (i = 0; i < height; ++i) {
		tempBase = data[width * i + a];
		data[width * i + a] = data[width * i + b];
		data[width * i + b] = tempBase;
		tempWeight = weights[a];
		weights[a] = weights[b];
		weights[b] = tempWeight;
	}
}

// The MP tree cannot be worse than an MST.
//HACK: Actually not true sometimes if there are ambiguous bases, since the
// MST algo lets an ambiguous base simultaneously agree with two distinct
// bases (e.g. it would allow A-N-T with a cost of zero, whereas in fact one
// of these edges must have a mutation).  The brute force solution would be
// to try assigning a fixed base to each ambiguous base, and work out the MST
// on each such dataset, and take the minimum of all these, but this would
// result in 4**n different MSTs being calculated if there were n 'N' bases!
// This isn't feasible.  It may be possible to extend the MST algo so that it
// will still efficiently produce the exact MST even with the existence of
// ambiguous bases, but this seems *hard*.  Since we are looking for an upper
// bound, we can safely compute a value that is too high, so we can use this
// approximation instead: (crappy partially thought-out algo deleted)
//
// Another, possibly easier way: when building the MST, as soon as a sequence
// containing an ambiguous base is added to the tree, *shrink* it down to the
// intersection of {what it was} and {what it just joined with}.  (Pretty
// similar to how the Fitch algo works.)  Only do this if the ambiguous base
// just added overlaps with the base it just joined with (same as for Fitch).
// This is a greedy approach, but it should work, and it *may* even give the
// actual optimal tree (although I would need to think about that pretty hard).
// TODO!  (This is still broken at the moment!)
// WTJW 23/3/2005: Now ask CalcTreeMstWeight() to overestimate the MST whenever
// ambiguous bases are involved.  Should now be safe to use.
void ImproveInitialUpperBound(struct TreeData *td) {
	unsigned mstWeight = CalcTreeMstWeight_b1(td, 0) + td->constantWeight;
	
	if (mstWeight < td->bound) {
		if (td->options & OPT_DISPLAYNEWBOUNDS) fprintf(stderr, "Improving initial upper bound to %u, based on minimum spanning tree weight.\n", mstWeight);
		td->bound = mstWeight;
	}
}

//TODO: Perform TBR on the tree produced by BuildGreedyTree().  Also try several different taxon orderings.
//HACK: Unlike ImproveInitialUpperBound(), which must be called before BlockifyTaxa(), this must actually be called after.
void ImproveInitialUpperBoundViaTbr(struct TreeData *td) {
	unsigned *taxonOrder = malloc(td->numTaxa * sizeof (unsigned));
	unsigned i;
	struct tree *root;
	unsigned score;
	void *arenaSavePoint;

	// Start by trying all taxa in the order in which they currently exist (which is most likely not
	// the order they appeared in the input file, since OrderTaxa() has already been called)
	for (i = 0; i < td->numTaxa; ++i) {
		taxonOrder[i] = i;
	}

	arenaSavePoint = ArenaGetPos(&td->alignedMem);		// Enable "one big free()" at the end
	root = BuildGreedyTree(td, taxonOrder, &score);
	assert(VerifyTree(root->nodes[LEFT], root));

	TreeBisectionReconnection(root, td);
	assert(VerifyTree(root->nodes[LEFT], root));

	DestroyTree(root);			// We're actually not really interested in the tree
	ArenaSetPos(&td->alignedMem, arenaSavePoint);
	free(taxonOrder);
}

// Only works on pre-SQUEEZEBASES data.
// Returns a lower bound on the number of distinct bases.
// If there are no ambiguous bases, the exact number is of course returned.
// Otherwise, a score of zero could be returned (if all sites are ambiguous).
// This could possibly be improved in future by looking for bases that, although
// ambiguous, do not intersect -- but that seems expensive.  (Consider the 3
// ambiguous bases [AC], [CG] and [GT]: [AC] and [GT] don't overlap, so there are
// at least 2 distinct bases here, but I can't see how to discover that without
// an all-vs-all comparison which is maybe overkill.)
unsigned NumDistinctBasesLowerBound_b1t(unsigned char *site, unsigned numTaxa) {
	unsigned nBase[16] = { 0 };
	unsigned i;
	
	for (i = 0; i < numTaxa; ++i) {
		assert(site[i] <= 15);			//DEBUG
		++nBase[site[i]];
	}
	
	return !!nBase[1] + !!nBase[2] + !!nBase[4] + !!nBase[8];
}

// This routine expects the data in TRANSPOSED FORM.
// effNumTaxa is the number of taxa to consider, starting from taxon zero.  This
// allows you to consider whether a site is PI or not with respect to just an
// initial subset of taxa.  Returns the constant weight that would be added
// to any tree by this site (if it is non-PI), or -1 if it is PI.  If the
// site is constant, the weight returned will be 0.
int ParsimoniouslyUninformativeWeight_b1t(unsigned char *buf, unsigned effNumTaxa) {
	unsigned i, j;
	unsigned containsTheRest, disjoint;
	unsigned char theRest;
	
	// Is this site parsimoniously uninformative?
	// When only dealing with a single base for each taxon, this
	// question translates to: count the frequencies of each base
	// at this site; if all but 1 of the frequencies are 0 or 1,
	// the site is not parsimoniously informative.
	// When ambiguous base sets are allowed for each taxon, it's a
	// bit harder.  See below.
	for (theRest = 1; theRest <= 16; theRest <<= 1) {
		// Posit base theRest as "the rest", and see whether all base sets
		// not containing this base are pairwise disjoint from each
		// other.  If so, then this column is parsimoniously
		// uninformative.
		disjoint = 1;
		containsTheRest = 0;
		for (i = 0; i < effNumTaxa && disjoint; ++i) {
			if (buf[i] & theRest) {
				++containsTheRest;
				continue;		// Consider only base sets not containing theRest
			}
			
			for (j = i + 1; j < effNumTaxa; ++j) {
				if (buf[j] & theRest) continue;		// Consider only base sets not containing theRest
				if (buf[i] & buf[j]) {
					// These two base sets are not pairwise disjoint
					disjoint = 0;
					break;
				}
			}
		}
		
		if (disjoint) {
			return effNumTaxa - containsTheRest;			// Parsimoniously uninformative, having this constant weight
		}
	}
	
	return -1;		// Parsimoniously informative
}

// This routine is similar to ParsimoniouslyUninformativeWeight_b1t() above.  It
// expects the data in TRANSPOSED FORM.  The difference is that instead of
// returning -1 for PI sites, it returns an upper bound for this site's weight
// on any MP tree.  It so happens that a non-PI site's weight is equal to this
// upper bound on ANY MP tree.
unsigned SiteWeightUpperBound_b1t(unsigned char *site, unsigned height) {
	unsigned i;
	unsigned biggest, biggestOverall = 0;
	unsigned char theRest;
	
	// If only non-ambiguous bases were concerned, the easiest way would be
	// to find the frequencies of each base and set m to the largest.
	// We could do that here, to provide a conservative upper bound, however
	// we can do better by allowing any ambiguous base set to be cast to
	// the most popular base it contains.
	for (theRest = 1; theRest <= 8; theRest <<= 1) {
		biggest = 0;
		for (i = 0; i < height; ++i) {
			if (site[i] & theRest) {
				++biggest;
			}
		}
		
		if (biggest > biggestOverall) {
			biggestOverall = biggest;
		}
	}
	
	return height - biggestOverall;
}

// Functionally identical to SiteWeightUpperBound_b1t() except that it expects the data in non-transposed form.
unsigned SiteWeightUpperBound_b1(unsigned char *data, unsigned iSite, unsigned height, unsigned memWidth) {
	unsigned i;
	unsigned biggest, biggestOverall = 0;
	unsigned char theRest;
	
	// If only non-ambiguous bases were concerned, the easiest way would be
	// to find the frequencies of each base and set m to the largest.
	// We could do that here, to provide a conservative upper bound, however
	// we can do better by allowing any ambiguous base set to be cast to
	// the most popular base it contains.
	for (theRest = 1; theRest <= 8; theRest <<= 1) {
		biggest = 0;
		for (i = 0; i < height; ++i) {
			if (data[i * memWidth + iSite] & theRest) {
				++biggest;
			}
		}
		
		if (biggest > biggestOverall) {
			biggestOverall = biggest;
		}
	}
	
	return height - biggestOverall;
}

// Computes the number of distinct bases minus 1.  Currently handles ambiguous
// bases safely but not optimally.
//HACK: Probably don't need both this function and NumDistinctBases(), but for
// now I can't remember whether the latter **needs** to ignore ambiguous
// bases -- if it doesn't then I should amalgamate the functions, otherwise I
// should keep them separate and optimise this function to handle ambiguous
// bases better.
unsigned SiteWeightLowerBound_b1t(unsigned char *site, unsigned height) {
	unsigned x = NumDistinctBasesLowerBound_b1t(site, height);
	if (x) {		// NumDistinctBasesLowerBound() can return 0 if all bases are ambiguous.
		--x;
	}

	return x;
}

// dest is allowed to point to s1 or s2, or to neither.
void UnionOfSeqs(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned len) {
	unsigned i;
	
	for (i = 0; i < len; ++i) {
		dest[i] = s1[i] | s2[i];
	}
}

// dest is allowed to point to s1 or s2, or to neither.  Will compute Hamming
// distances between sequences where sets of bases are allowed at each site.
unsigned HammingDist(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len) {
	unsigned i;
	unsigned dist = 0;
	
	for (i = 0; i < len; ++i) {
#ifdef SQUEEZEBASES
		if (!((s1[i] & s2[i]) & 0x0F)) {
			dist += weights[i * 2];
		}
		if (!((s1[i] & s2[i]) & 0xF0)) {
			dist += weights[i * 2 + 1];
		}
#else	// not SQUEEZEBASES
		if (!(s1[i] & s2[i])) {
			dist += weights[i];
		}
#endif	// not SQUEEZEBASES
	}
	
	return dist;
}


// Used as a comparator in calls to SortSitesAndWeights().
// Skips over the 1st sizeof (unsigned) bytes, which contain the site's weight.
int ColumnCompareIncludingWeights(const void *a, const void *b) {
	return memcmp((unsigned char *) a + sizeof (unsigned), (unsigned char *) b + sizeof (unsigned), _ColumnCompareWidth);		//HACK: uses global: _ColumnCompareWidth.  **Slightly** better than before.
}

// Used as a comparator in calls to SortSitesAndWeights().
// The site's weight is stored in the first (sizeof (int)) bytes of each column.
// Sort primarily by that, then lexicographically by character states to ensure stability
// even if sequences are given in a different order (all other things being equal).
int ColumnCompareWeights(const void *a, const void *b) {
	int x;
	x = *((int *) a) - *((int *) b);		// Technically they're unsigned, but C guarantees these types are the same size.
	if (!x) {
		x = -memcmp((unsigned char *) a + sizeof (unsigned), (unsigned char *) b + sizeof (unsigned), _ColumnCompareWidth);	//HACK: Marginally better than using gTD.numTaxa...
	}
#ifdef FLIPCOLS
	return x;
#else	// not FLIPCOLS
	return -x;
#endif	// not FLIPCOLS
}

// NOTE: transposes bytes.
void Transpose_b1(unsigned char *data, unsigned w, unsigned h) {
	unsigned char *buf = (unsigned char *) malloc(w * h);
	unsigned i, j;
	
	for (i = 0; i < h; ++i) {
		for (j = 0; j < w; ++j) {
			buf[j * h + i] = data[i * w + j];
		}
	}
	
	memcpy(data, buf, w * h);
	free(buf);
}

//WTJW: Now reads interleaved format files, and gives much better error reporting.
void ReadPhylipAlignment(FILE *input_file, struct CharData *cd) {
	char buf[80];		// Only needs to be big enough to hold the first line
	unsigned i, col, charsRead = 0, line = 2;
	int blankLine;
	int pos, blockWidth;
	static char *allowedChars = "ACGTURYSWKMBDHV-N?";
	char c;
	
	// Read first line
	fgets(buf, sizeof(buf), input_file);
	if (feof(input_file)) {
		fprintf(stderr, "Error reading input file: no first line!  Aborting.\n");
		exit(1);
	}
	
	sscanf(buf, "%d %d", &(cd->numTaxa), &(cd->seqLen));
	
	cd->charMat = (unsigned char *) malloc(cd->numTaxa * cd->seqLen);
	cd->labels = (char **) malloc(cd->numTaxa * sizeof (char *));
	
	while (charsRead < cd->seqLen) {
		for (i = 0; i < cd->numTaxa; ++i, ++line) {
			if (feof(input_file)) {
				if (i) {
					fprintf(stderr, "Error reading input file: %u sites specified in top line, but input only contains %u sites for first %u taxa and %u sites for remaining %u taxa, aborting.\n", cd->seqLen, charsRead + blockWidth, i, charsRead, cd->numTaxa - i);
				} else {
					fprintf(stderr, "Error reading input file: %u sites specified in top line, but input only contains %u sites, aborting.\n", cd->seqLen, charsRead);
				}
				exit(1);
			}
			
			// Read sequence data
			pos = 0;
			blankLine = 1;
			col = 0;
			while (1) {
				c = fgetc(input_file);
				
				if (feof(input_file) || c == '\n') {
					break;
				}
				
				if (!isspace(c)) {
					blankLine = 0;
				}

				if (!charsRead && col < 10) {
					// Label
					buf[col++] = c;
				} else {
					// Sequence data
					if (strchr(allowedChars, toupper(c))) {
						if (charsRead + pos == cd->seqLen) {
							fprintf(stderr, "Error reading input file: too many characters on line %u, aborting.\n", line);
							exit(1);
						} else {
							cd->charMat[i * cd->seqLen + charsRead + pos++] = toupper(c);
						}
					} else if (isspace(c) || isdigit(c)) {		// PHYLIP format explicitly allows digits, which are ignored.
						// Do nothing.
					} else {
						fprintf(stderr, "Error reading input file: invalid character '%c' on line %u, aborting.\n", c, line);
						exit(1);
					}
				}
			}
			
			// Ignore blank lines
			if (blankLine) {
				--i;
				continue;
			}

			if (!i) {
				blockWidth = pos;
			}

			if (!charsRead) {
				// This way allows spaces *within* the taxon label, and also handles short labels, trailing \r and/or \n etc.
				buf[col--] = 0;
				for (; !buf[col] || isspace(buf[col]); --col) {
					buf[col] = 0;
				}

				cd->labels[i] = malloc(strlen(buf) + 1);
				strcpy(cd->labels[i], buf);
			}

			if (blockWidth != pos) {
				fprintf(stderr, "Error reading input file: line %u contains %u %s sites than the first line in this group, aborting.\n", line, pos > blockWidth ? pos - blockWidth : blockWidth - pos, pos > blockWidth ? "more" : "fewer");
				exit(1);
			}
		}
		
		charsRead += pos;
	}
	
	return;
}

//WTJW: now updates parent pointers #ifdef BACKPOINTERS.  The root node
// will have a parent of NULL.
// 7/2/2003: Now knows to deep-copy internal nodes and shallow-copy leaves.
// blockified == 0 until BlockifyTaxa() is called.
// Assumes the data has been blockified
struct tree* CopyTree(struct tree *root, struct TreeData *td, int saveSeq) {
	struct tree *copy;
	
	assert(!saveSeq);		// We don't handle this yet.  We probably don't want to save the entire history of each seq.

	INCMEASURE(copyTreeCount);
	
	copy = malloc(sizeof(struct tree));
	copy->contents = root->contents;
	copy->seqsByTreeSize = NULL;			// Doing this can save a lot of memory when there are thousands of trees!
#ifdef BACKPOINTERS
	copy->nodes[PARENT] = NULL;	// Will be overwritten unless this is the root
#endif	// BACKPOINTERS
#ifdef SMARTFITCH
	copy->cumScore = root->cumScore;
#endif	// SMARTFITCH
	if (root->nodes[LEFT] == NULL) {
		copy->nodes[LEFT] = NULL;
	} else {
		copy->nodes[LEFT] = CopyTree(root->nodes[LEFT], td, saveSeq);
#ifdef BACKPOINTERS
		copy->nodes[LEFT]->nodes[PARENT] = copy;		// WTJW 4/4/2003: Fixed bug
#endif	// BACKPOINTERS
	}
	if (root->nodes[RIGHT] == NULL) {
		copy->nodes[RIGHT] = NULL;
	} else {
		copy->nodes[RIGHT] = CopyTree(root->nodes[RIGHT], td, saveSeq);
#ifdef BACKPOINTERS
		copy->nodes[RIGHT]->nodes[PARENT] = copy;		// WTJW 4/4/2003: Fixed bug
#endif	// BACKPOINTERS
	}
	
	return copy;
}

//WTJW: Recursively deallocate tree nodes
void DestroyTree(struct tree *root) {
	INCMEASURE(destroyTreeCount);
	
	if (!root) return;			// We've gone one past the end
	
	// 7/2/2003: Free sequence memory allocated to internal nodes only
	if (root->contents == -1) {
		// WTJW 8/2/2005: Only free memory if it was already allocated!
//		free(root->sequence);
		// Sequences are allocated by the arena, and freed all at once.
		//if (root->sequence) {
		//	free(root->sequence);
		//}
	}
	DestroyTree(root->nodes[LEFT]);
	DestroyTree(root->nodes[RIGHT]);
	free(root);
}

void DestroyCharData(struct CharData *cd) {
	unsigned i;

	free(cd->charMat);
	
	for (i = 0; i < cd->numTaxa; ++i) {
		free(cd->labels[i]);
	}
	free(cd->labels);
}

void DestroyTreeData(struct TreeData *td) {
#ifdef SITESCORES
	unsigned i;
#endif	// SITESCORES
	free(td->restBound);
#ifdef FASTLOWWEIGHTFITCH
	free(td->byteWeights);
#endif	// FASTLOWWEIGHTFITCH
	free(td->taxonMap);
	TreeList_Destroy(&td->bestList);
	// td->charMat and td->weights are actually (re-)allocated in BlockifyTaxa() using ArenaAllocate()!
	//ALIGNED_FREE(td->weights);
	//free(td->charMat);
	ArenaDestroy(&td->alignedMem);
#ifdef SITESCORES
	free(td->totalSiteScores);
	for (i = 0; i < td->numTaxa; ++i) {
		free(td->boundSiteScoresByTaxon[i]);
	}
	free(td->boundSiteScoresByTaxon);
#endif	// SITESCORES
#ifdef LOGTREES
	if (td->treeLogFile) {
		fclose(td->treeLogFile);
		td->treeLogFile = NULL;
	}
#endif	// LOGTREES
	if (td->workStack) {
		free(td->workStack);
	}
#ifdef TARGETMULTI
	assert(!td->mpiRequests);		// This should be cleaned up by the routine that allocated it.
#endif	// TARGETMULTI
	if (td->tempTreeFile) {
		fclose(td->tempTreeFile);
		td->tempTreeFile = NULL;
		assert(td->tempTreeFName);
		if (!(td->options & OPT_KEEPTEMPFILES)) {
			remove(td->tempTreeFName);
			DBGPRINT2("Removed temporary tree file %s.\n", td->tempTreeFName);
		}
	}
	if (td->tempTreeFName) {
		free(td->tempTreeFName);
		td->tempTreeFName = NULL;
	}
}

//HACK: very inefficient.  I should create a reverse taxon map once and store it in an array, just like taxonMap[].
unsigned ReverseMapTaxon(unsigned label, struct TreeData *td) {
	unsigned i;
	
	for (i = 0; i < td->numTaxa; ++i) {
		if (td->taxonMap[i] == label) return i;
	}
	
	return INT_MAX;				// Not found
}

void TreeList_Init(TreeList *ptl) {
	*ptl = NULL;
}

// Uses simple, probably inefficient recursive algorithm -- but efficiency is
// not too important here.
void TreeList_Destroy(TreeList *ptl) {
	if (!*ptl) return;
	
	INCMEASURE(treeList_DestroyCount);
	
	TreeList_Destroy(&(*ptl)->next);
	DestroyTree((*ptl)->t);			// This tree was a copy
	free(*ptl);
	*ptl = NULL;
}

void TreeList_AddTree(TreeList *ptl, struct tree *t, struct TreeData *td) {
	struct TreeListNode *tl = (struct TreeListNode *) malloc(sizeof (struct TreeListNode));
	
	INCMEASURE(treeList_AddTreeCount);
	
	// WTJW 8/2/2005: Don't save the sequence data of trees.
//	tl->t = CopyTree(t, td);		// Important to make a copy of the tree!
	tl->t = CopyTree(t, td, 0);		// Important to make a copy of the tree!
	tl->next = *ptl;
	*ptl = tl;
}

TreeListIterator TreeList_GetIterator(TreeList tl) {
	return tl;		// Actually, we do nothing fancy here for the time being.
}

// Passing in the TreeList pointer allows for implementation of lightweight
// iterators in the future.
#include "nowarnunusedparam.h"
struct tree *TreeList_GetTree(PARAMNOTUSED(TreeList tl), TreeListIterator pos) {
	return pos->t;
}
#include "warnpop.h"

// See comment for TreeList_GetTree().  Not really HasNext(), more like HasThis().
#include "nowarnunusedparam.h"
int TreeList_HasNext(PARAMNOTUSED(TreeList tl), TreeListIterator pos) {
	return pos != NULL;
}
#include "warnpop.h"

// See comment for TreeList_GetTree().
#include "nowarnunusedparam.h"
TreeListIterator TreeList_GetNext(PARAMNOTUSED(TreeList tl), TreeListIterator pos) {
	return pos->next;
}
#include "warnpop.h"

void SaveTreeList(struct TreeData *td, FILE *f) {
	unsigned i;
	TreeListIterator tli;
	
	if (td->options & OPT_NEXUSFMT) {
		fprintf(f,
			"#NEXUS\n"
			"[ Trees produced by " PROGLONGNAME "\n"
			"> %u tree(s), each having score %u. ]\n"
			"BEGIN TREES;\n"
			"\tTRANSLATE\n", td->numTrees, td->bound);
		
		for (i = 0; i < td->numTaxa; ++i) {
			fprintf(f, "\t\t%u\t%s%c\n", i + 1, td->cd->labels[i], i < td->numTaxa - 1 ? ',' : ';');
		}
		
		i = 1;
	}
	
	for (tli = TreeList_GetIterator(td->bestList); TreeList_HasNext(td->bestList, tli); tli = TreeList_GetNext(td->bestList, tli)) {
		if (td->options & OPT_NEXUSFMT) {
			fprintf(f, "\tTREE FASTDNAMP_%u = [&U] ", i++);
		}
		
		PrintTree(TreeList_GetTree(td->bestList, tli), td, f);
		fprintf(f, ";\n");
	}
	
	if (td->options & OPT_NEXUSFMT) {
		fprintf(f, "END;\n");
	}
}

#ifdef DEBUG
void DebugPrintSwitches(void) {
	DBGPRINT1("--- Switches enabled at compile-time ---\n");
#ifdef TARGETMULTI
	DBGPRINT1("TARGETMULTI\n");
#endif	// TARGETMULTI
	DBGPRINT2("SITESPERBYTE=%u\n", SITESPERBYTE);
	DBGPRINT2("BYTESPERBLOCK=%u\n", BYTESPERBLOCK);
	DBGPRINT1("FITCHWHOLESEQ\n");
#ifdef EXHAUSTIVE
	DBGPRINT1("EXHAUSTIVE\n");
#endif	// EXHAUSTIVE
#ifdef FASTCFITCH
	DBGPRINT1("FASTCFITCH\n");
#endif	// FASTCFITCH
#ifdef FASTC64FITCH
	DBGPRINT1("FASTC64FITCH\n");
#endif	// FASTC64FITCH
#ifdef X86ASMFITCH
	DBGPRINT1("X86ASMFITCH\n");
#endif	// X86ASMFITCH
#ifdef MMXASMFITCH
	DBGPRINT1("MMXASMFITCH\n");
#endif	// MMXASMFITCH
#ifdef SSE2ASMFITCH
	DBGPRINT1("SSE2ASMFITCH\n");
#endif	// SSE2ASMFITCH
#ifdef FASTWEIGHT1FITCH
	DBGPRINT1("FASTWEIGHT1FITCH\n");
#endif	// FASTWEIGHT1FITCH
#ifdef FASTLOWWEIGHTFITCH
	DBGPRINT1("FASTLOWWEIGHTFITCH\n");
#endif	// FASTLOWWEIGHTFITCH
#ifdef SMARTFITCH
	DBGPRINT1("SMARTFITCH\n");
#endif	// SMARTFITCH
#ifdef EMPTYFITCH
	DBGPRINT1("EMPTYFITCH\n");
#endif	// EMPTYFITCH
#ifdef REVERSETAXA
	DBGPRINT1("REVERSETAXA\n");
#endif	// REVERSETAXA
#ifdef MAXMINIORDER
	DBGPRINT1("MAXMINIORDER\n");
#endif	// MAXMINIORDER
#ifdef MAXMINITREEORDER
	DBGPRINT1("MAXMINITREEORDER\n");
#endif	// MAXMINITREEORDER
#ifdef CACHEOPT
	DBGPRINT3("CACHEOPT (CACHESIZE=%d, CACHEADDRMASK=%x)\n", CACHESIZE, CACHEADDRMASK);
#endif	// CACHEOPT
#ifdef SQUEEZEBASES
	DBGPRINT1("SQUEEZEBASES\n");
#endif	// SQUEEZEBASES
#ifdef BIGENDIAN
	DBGPRINT1("BIGENDIAN\n");
#endif	// BIGENDIAN
#ifdef LITTLEENDIAN
	DBGPRINT1("LITTLEENDIAN\n");
#endif	// LITTLEENDIAN
#ifdef DYNAMICORDER
	DBGPRINT1("DYNAMICORDER\n");
#endif	// DYNAMICORDER
#ifdef DYNORDER_SUMWEIGHTS
	DBGPRINT1("DYNORDER_SUMWEIGHTS\n");
#endif	// DYNORDER_SUMWEIGHTS
#ifdef DYNORDER_MAXWEIGHT
	DBGPRINT1("DYNORDER_MAXWEIGHT\n");
#endif	// DYNORDER_MAXWEIGHT
#ifdef DYNORDER_SUMTHENMAX
	DBGPRINT1("DYNORDER_SUMTHENMAX\n");
#endif	// DYNORDER_SUMTHENMAX
#ifdef DYNORDER_MAXTHENSUM
	DBGPRINT1("DYNORDER_MAXTHENSUM\n");
#endif	// DYNORDER_MAXTHENSUM
#ifdef DYNORDER_QUICK
	DBGPRINT1("DYNORDER_QUICK\n");
#endif	// DYNORDER_QUICK
#ifdef DETECTREDUNDANTFITCH
	DBGPRINT1("DETECTREDUNDANTFITCH\n");
#endif	// DETECTREDUNDANTFITCH
	DBGPRINT1("YANBADER\n");
#ifdef STOPEARLY
	DBGPRINT1("STOPEARLY\n");
#endif	// STOPEARLY
#ifdef WEIGHT1MULTIPLY
	DBGPRINT1("WEIGHT1MULTIPLY\n");
#endif	// WEIGHT1MULTIPLY
#ifdef MEASURE
	DBGPRINT1("MEASURE\n");
#endif	// MEASURE
	DBGPRINT1("--- End of switches enabled at compile-time ---\n");
}
#endif	// DEBUG

// The unsigned long datatype affords calculations for trees with up to 19 taxa.
unsigned long CountSubTrees(unsigned numTaxaOnTree, unsigned numTaxa) {
	unsigned long n = 1;
	unsigned i;
	
	if (numTaxa <= 2) return 0;					// No trees on <= 2 taxa
	if (numTaxaOnTree <= 2) numTaxaOnTree = 3;	// Choose any three taxa for the base tree
	for (i = numTaxaOnTree; i < numTaxa; ++i) {
		n *= (i - 2) * 2 + 1;
	}
	
	return n;
}
