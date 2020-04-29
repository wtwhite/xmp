#ifndef __COMMON_H
#define __COMMON_H
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "switches.h"
#include "timer.h"
#include "elapsedtime.h"
#include "sitepaircosts.h"
#include "partbound.h"		// Need this for enum partbound_strategy_type
#include "arena.h"
#ifdef TARGETMULTI
#include "mpi.h"
#endif	// TARGETMULTI

#define PROGNAME "fastDNAmp"
#define PROGLONGNAME "fastDNAmp: Exact Maximum Parsimony Phylogenetic Tree Search v1.1"

// EMMS matters only for x86 versions that use MMX registers.  That means MMXASMFITCH, and possibly SSE2ASMFITCH
// (although it probably only uses XMM registers, in which case it doesn't need EMMS protection either).
#if defined(MMXASMFITCH) || defined(SSE2ASMFITCH)
#ifndef __linux__
#define EMMS __asm emms
#else	// not __linux__
#define EMMS asm("emms")
#endif	// not __linux__
#else	// NOT defined(MMXASMFITCH) || defined(SSE2ASMFITCH)
#define EMMS
#endif	// NOT defined(MMXASMFITCH) || defined(SSE2ASMFITCH)

//WTJW: MAXTAXA and MAXSEQLEN exist only so that certain functions can use local
// arrays and avoid calling malloc() on every iteration.  I'll remove them
// eventually.
#define LEFT 0
#define RIGHT 1
#ifdef BACKPOINTERS
// Incorporate backpointer as a 3rd entry in this array
#define MAXDEGREE 3
#define PARENT 2
#else	// not BACKPOINTERS
#define MAXDEGREE 2
#endif	// not BACKPOINTERS
#define MAXTAXA 100	//HACK

// A 2Mb buffer seems reasonable.  If TARGETMULTI is defined, a buffer this size will be used to read chunks
// from each worker; otherwise, a buffer this size will be used to copy our own temporary tree file
// to the final destination.
#define TREEBUFSIZE (2 * 1024 * 1024)

// Bitmask #defines for TreeData.options.
#define OPT_NEXUSFMT 1
#define OPT_DISCARDMISSAMBIG 2
#define OPT_DISPLAYNEWBOUNDS 4
#define OPT_DISPLAYNEWTREES 8
#define OPT_QUIET 16
#define OPT_USELOADRESTBOUND 32
#define OPT_USEDISTBASESBOUNDREST 64
#define OPT_USEINCOMPATBOUNDREST 128
#define OPT_USEMSTBOUNDREST 256
#define OPT_USEKICKASSBOUNDREST 512
#define OPT_USEKA2BOUNDREST 1024
#define OPT_USEGREEDYHITTINGSETBOUND 2048
#define OPT_IMPROVEINITUPPERBOUND 4096
#define OPT_USEPARTBOUND 8192
#define OPT_ONLYCOMPUTELOWERBOUND 16384
#define OPT_KEEPTEMPFILES 32768
#define OPT_SAVEREORDEREDDATASET 65536
#define OPT_TBR 131072

// Bitmask #defines for use with WritePhylipAlignment().
#define WPA_PRINTWEIGHTS 1
#define WPA_ALIGN 2

//WTJW: #defines for base bitmasks
#define BM_A 1
#define BM_C 2
#define BM_G 4
#define BM_T 8

// Every tree containing this many taxa will be counted as one "unit" of search
// space by proportionComplete and actualTreesExamined in the TreeData structure.
//#define COUNT_UNIT 8
#define COUNT_UNIT 9
//#define COUNT_UNIT 10
// This is how many trees of size COUNT_UNIT there are.  Should be (2*COUNT_UNIT-5)!!.
//#define COMPLETE_UNITS 10395
#define COMPLETE_UNITS 135135
//#define COMPLETE_UNITS 2027025

// Each edge in a tree has 3 sequences associated with it for the purposes of the
// Yan-Bader Fitch algorithm.  These are stored in an array of struct edge_seqs,
// and indexed by the constants below.  Each node stores the 3 sequences for the
// edge that is vertically above it.
#define MIDPOINT 0		// Put this at offset 0 because it's used by both FitchBases() and FitchScore().
#define UP 1
#define DOWN 2

struct edge_seqs {
	unsigned char *seqs[3];		// An array of 3 pointers to unsigned char
};

struct tree {
	int contents;
	struct tree *nodes[MAXDEGREE];
	//HACK: yes, I know this change (now, 8/10/2009, *back* to "unsigned char *") is likely to quietly break a heap of things...
	//unsigned char *sequence;
	// seqsByTreeSize[3] contains pointers to the 3 seqs for this edge with 3 taxa on the tree.
	// Although this is slightly wasteful because edges introduced later in the tree will never access
	// the lower-indexed elements, it simplifies access.  (And perhaps we can subtract some amount off
	// the pointer to the allocated space to pack things more tightly if we use our own arena-based
	// memory allocation system.)
	struct edge_seqs *seqsByTreeSize;
#ifdef FITCHDETECTSAME
	unsigned score;			//WTJW: score increment incurred by Fitching the two immediate children of this node
#endif
#ifdef SMARTFITCH
	unsigned cumScore;		// What is the cost, up to and including this node?
#endif
};

struct edge {
	struct tree *parent;	
	int child;	//LEFT or RIGHT
};

//WTJW: produced from a FILE* by readPhylipAlignment()
struct CharData {
	unsigned char *charMat;
	char **labels;
	unsigned numTaxa;
	unsigned seqLen;
};

// Tree Lists
// ==========
// The idea is that TreeList and TreeListIterator are opaque types, and that
// all accesses to trees in the list should be done through the "member"
// functions starting with "TreeList_".  This is so that the underlying
// implementation of the tree list can be changed later without affecting
// clients of the list.  The current implementation, which is just a basic
// linked list, is functional but will be pretty slow (and incur non-negligible
// overhead) if many trees will be put in the list.
struct TreeListNode {
	struct tree *t;
	struct TreeListNode *next;
};

typedef struct TreeListNode *TreeList;
typedef TreeList TreeListIterator;

//WTJW: produced from a CharData by initTreeData()
struct TreeData {
	unsigned char *charMat;		// Contains enough space for all internal nodes
	unsigned numTaxa;
	unsigned seqLen;			// Length of sequence in bases.  Not necessarily the same as the length in memory in bytes!
	//unsigned memSeqLen;			// Length in memory.  Can be multiplied by i and added to charMat to get a pointer to sequence #(i+1).
	unsigned seqLenBlocks;		// Length in BYTESPERBLOCK-byte blocks.  Only valid after BlockifyTaxa() has been called.
	unsigned *weights;			// Array of seqLen weights
	unsigned constantWeight;	// Non-PI sites are responsible for this much weight in any tree
	unsigned *taxonMap;			// taxonMap[42] is the original taxon number for taxon #43 after taxa are reordered by InitTreeData()
#ifdef FASTWEIGHT1FITCH
	unsigned heavySeqLen;
#endif	// FASTWEIGHT1FITCH
#ifdef FASTLOWWEIGHTFITCH
	unsigned char *byteWeights;	// Each weight is 1 byte; each group of 4 weights must sum to <= 255
	unsigned heavySeqLen;
#endif	// FASTLOWWEIGHTFITCH
	unsigned *restBound;		// restBound[42] is the weight that must be added by taxa #44, #45, #46, ..., #numTaxa.
	char *restBoundInputFname;	// Name of file to load restBound[] array from
	enum partbound_strategy_type partBoundStrategy;		// Only used by -Bp
	int partBoundMaxItersNoImprovement;					// Only used by -Bp
	struct CharData *cd;		// Whence this TreeData was created
//	struct TreeList *bestList;
	TreeList bestList;
	unsigned options;
	unsigned maxTrees;
	unsigned numTrees;
	unsigned optimalTreesMissed;	// Has the maxTrees limit caused us to miss other optimal trees?
	unsigned startCol;			// Not used after charMat matrix created
	unsigned endCol;			// Not used after charMat matrix created
	unsigned proportionComplete;		// Actually counts a tree containing COUNT_UNIT taxa as 1
	unsigned treesActuallyExamined;		// Actually counts a tree containing COUNT_UNIT taxa as 1
	unsigned updateProgressNow;	// Set to 1 periodically by timer-actioned routine
	unsigned updateIntervalSecs;
	TimerHandle timerHandle;
	TimePoint startTime;				// Program start time.  Durations can be found by subtracting from the next stage's start time.
	TimePoint orderTaxaStartTime;		// Time at which OrderTaxa() is called
	TimePoint upperBoundStartTime;		// Time at which upper-bound (heuristic tree) computations start
	TimePoint lowerBoundStartTime;		// Time at which lower-bound computations start
	TimePoint branchAndBoundStartTime;	// Time at which B&B search starts
	TimePoint outputResultsStartTime;	// Time at which output file is generated
	TimePoint endTime;					// Just before program exits
	unsigned bound;
#ifdef DYNAMICORDER
	unsigned char *treeStateCounts[4];
	unsigned char *restStateCounts[4];
#endif	// DYNAMICORDER
#ifdef SITESCORES
	unsigned *totalSiteScores;		// Each entry contains the score for BYTESPERBLOCK sites (or BYTESPERBLOCK * 2 #ifdef SQUEEZEBASES)
	unsigned *boundSiteScores;		// Each entry contains the score for BYTESPERBLOCK sites (or BYTESPERBLOCK * 2 #ifdef SQUEEZEBASES)
	unsigned **boundSiteScoresByTaxon;	// x[3][5] contains the score for the 6th group of BYTESPERBLOCK (* 2) sites for trees containing 3 taxa.  restBound[3] can be safely added.
#endif	// SITESCORES
	//unsigned *lengths;				// Number of blocks of sites, indexed by number of taxa on the tree.  Multiply by BYTESPERBLOCK to get memory sizes, and additionally by SITESPERBYTE (1 or 2) for site counts.
	//unsigned *offsets;				// Cumulative total of above.  Use to access into weights[], or multiply by three and use to access into sequence data
	unsigned score;					// The score of the current tree is maintained.
	//int *constantWeights;			// This value will be negative for all indices above 3!  It is an **adjustment** value.
	//unsigned char **taxonSequence;	// taxonSequence[5] points to the sequence data for taxon #6, and has room for the three fields ("-->", "<--", "potential root").  Only the first field is filled out in InitTreeData().
//	unsigned *memBuffer;			// Marker for an "arena" used for allocating sequence buffers for internal nodes.  Should be faster than malloc(), and is valid since we always deallocate in reverse order of allocation.
#ifdef LOGTREES
	FILE *treeLogFile;
#endif	// LOGTREES
	unsigned rngSeed;				// Doesn't really belong here, but needs to be stored so ReportSettings() can report it.
	int rand1;						// See ProcessCmdLineArgs() for an explanation.
	int rand2;
	int rand3;
	// Array of pairs (i, j) used for job stealing.  workStack[k] == (i, j) means
	// "all subproblems with taxon k at edges numbered < i on the k-taxon tree have either been
	// processed or are being processed now; and there are j edges in total." (Taxa
	// and edges considered to number up from 0.)
	unsigned (*workStack)[2];
	unsigned jobTreeSize;			// Number of taxa on the tree for this job
#ifdef TARGETMULTI
	int rank;
	int nWorkers;
	MPI_Request *mpiRequests;		// An array of outstanding non-blocking MPI requests.  Some may be MPI_REQUEST_NULL.
	int nMpiRequests;				// The size of the above array.
	int nStealsFromUs;				// Number of steals that we have responded to but not yet completed (waited on).
	int pendingBoundBroadcast;		// 1 if we should send the boss a MSG_NEWBOUND at the nearest opportunity.
	int pendingAskForJob;			// 1 if we need to ask for a new job at the nearest opportunity.
	//int receivedUB;					// A broadcast bound that we have received
	unsigned *receiveBuf;			// The buffer used for all worker receives
	unsigned (**stealBuf)[2];		// Each element is a pointer to a buffer used to send a job
#endif	// TARGETMULTI
	char *tempTreeFName;			// The path to tempTreeFile.  We must store it because it may depend on info sent by the master node for TARGETMULTI.
	FILE *tempTreeFile;				// The file where current-best trees found by this node are temporarily written.
	struct arena alignedMem;		// Holds an arena for quickly allocating BYTESPERBLOCK-aligned memory in LIFO order
};

#ifdef TARGETMULTI
// Message tags used for both Boss -> Worker and Worker -> Boss MPI communication.
// For messages sent by workers to the boss, we have now combined MSG_ASKFORJOB and MSG_NEWBOUND into a single message
// type, MSG_WORKERREQUEST.  The payload will be 1 unsigned int in length
// if it indicates NEWBOUND (in which case that byte will contain the bound) or
// 0 bytes long if it indicates ASKFORJOB.  Combining these messages is needed
// to guarantee that they arrive in the order sent (otherwise
// deadlocks can occur).  It's just as important that MSG_NEWJOB remain a
// separate message type so that the boss can wait on it only while stealing.
enum message_type {
	MSG_WORKERREQUEST,		// W -> B
	MSG_NEWJOB,				// B -> W (response to MSG_ASKFORJOB), W -> B (response to MSG_STEALJOB)
	MSG_STEALJOB,			// B -> W
	MSG_NEWBOUND,			// B -> W (B tells everyone the new UB).  Note Ws tell B about new UBs they find using MSG_WORKERREQUEST.
	// Actually we should not have a MSG_FINISHED -- instead we should send a special MSG_NEWJOB to the worker.
	// Then on the worker side, every MSG_ASKFORJOB sent produces exactly one MSG_NEWJOB sent back.
//	MSG_FINISHED,			// B -> W
	MSG_TREES,				// W -> B (raw text block of newline-delimited Newick-format trees, sent at the very end)
	MSG_NTREESLEFT,			// W -> W (number of trees remaining; only sent if td->maxTrees != 0)
	MAXMSGTYPE
};

// One item for each type of nonblocking MPI communication request that a worker can simultaneously have active.
// These are stored in the array mpiRequests of TreeData, followed by up to nWorkers more requests for handling
// the completion of sends of stolen jobs.
enum worker_comm_requests {
	WCR_RCV_ALL,
	WCR_SND_ASKFORJOB,
	WCR_SND_NEWBOUND,
	MAX_WCR
};
#endif	// TARGETMULTI

struct tree *NewNode(int label, int nTaxaOnTree, struct TreeData *td);
void AddToPendantNode(struct tree *parentNode, int child, struct tree *childNode);
struct tree * FASTCALL InsertNode(struct tree *from, int child, int label, int nTaxaOnTree, struct TreeData *td);
void FASTCALL DeleteNode(struct tree *from, int child);

void SetupDefaults(struct TreeData *td);
int ProcessCmdLineArgs(int argc, char **argv, struct TreeData *td, char **pInFileName, FILE **pInFile, char **pOutFileName, FILE **pOutFile);
void ReportSettings(struct TreeData *td, char *inFileName, FILE *inFile, char *outFileName, FILE *outFile);
void OutputResults(struct TreeData *td, FILE *outFile);
void ReportTiming(struct TreeData *td);
void EstimateTimeRemaining(struct TreeData *td);
char *DurationToString(double secs, int precise, char *buf);
struct tree *ConstructBaseTree(struct TreeData *td, int t1, int t2, int t3);
void PrintUsage(void);
void ReadPhylipAlignment(FILE *input_file, struct CharData *cd);
void InitTreeData(struct CharData *cd, struct TreeData *td);
void WritePhylipAlignment(unsigned char *buf, unsigned seqLen, unsigned numTaxa, unsigned *weights, unsigned *taxonMap, FILE *outFile, unsigned options);
int RecodeMatrix_b1t(unsigned char *fromBuf, unsigned char *toBuf, unsigned seqLen, unsigned numTaxa, unsigned startCol, unsigned endCol, unsigned **pWeights, unsigned discardMissAmbig);
int RemoveUninfSites_b1t(unsigned char *buf, unsigned *weights, unsigned seqLen, unsigned numTaxa, unsigned effNumTaxa, unsigned *pConstantWeight, unsigned *pNonPiSites, unsigned *pNonPiWeight);
int CombineIdenticalSites_b1t(unsigned char *buf, unsigned *weights, unsigned seqLen, unsigned numTaxa, unsigned effNumTaxa);
void SortSitesAndWeights_b1t(unsigned char *fromBuf, unsigned char *toBuf, unsigned *fromWeights, unsigned *toWeights, unsigned seqLen, unsigned numTaxa, int (*comparator)(const void *, const void *));
unsigned *ReweightToPowersOfTwo(unsigned char **pData, unsigned *pSeqLen, unsigned numTaxa, unsigned **pWeights, unsigned blockWidth);
void LoadRestBound(struct TreeData *td);
void ImproveInitialUpperBound(struct TreeData *td);
void ImproveInitialUpperBoundViaTbr(struct TreeData *td);
void InitKickAssBoundRest(struct TreeData *td);
void InitKA2BoundRest(struct TreeData *td);		// WTJW 11/6/2004
unsigned PartitionAndScoreColumns(unsigned char *data, unsigned dataWidth, unsigned width, unsigned height, unsigned *weights);
unsigned KA2PartitionAndScoreColumns(unsigned char *data, unsigned dataWidth, unsigned width, unsigned height, unsigned *weights, unsigned *upperBoundWeights, unsigned *pNumSitesUsed);
unsigned NumDistinctBasesLowerBound_b1t(unsigned char *site, unsigned numTaxa);
int ParsimoniouslyUninformativeWeight_b1t(unsigned char *buf, unsigned effNumTaxa);
unsigned SiteWeightUpperBound_b1t(unsigned char *site, unsigned height);		// WTJW 11/6/2004
unsigned SiteWeightUpperBound_b1(unsigned char *data, unsigned iSite, unsigned height, unsigned memWidth);
unsigned SiteWeightLowerBound_b1t(unsigned char *site, unsigned height);		// WTJW 9/8/2004
void InitColBound(struct TreeData *td);
int ColumnCompareIncludingWeights(const void *a, const void *b);
int ColumnCompareWeights(const void *a, const void *b);
void DestroyCharData(struct CharData *cd);
void DestroyTreeData(struct TreeData *td);

void TreeList_Init(TreeList *ptl);
void TreeList_Destroy(TreeList *ptl);
void TreeList_AddTree(TreeList *ptl, struct tree *t, struct TreeData *td);
TreeListIterator TreeList_GetIterator(TreeList tl);
struct tree *TreeList_GetTree(TreeList tl, TreeListIterator pos);
int TreeList_HasNext(TreeList tl, TreeListIterator pos);
TreeListIterator TreeList_GetNext(TreeList tl, TreeListIterator pos);

void SaveTreeList(struct TreeData *td, FILE *f);

void UnionOfSeqs(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned len);
unsigned HammingDist(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len);

int VerifyTree(struct tree *root, struct tree *parent);
void PrintTree(struct tree *root, struct TreeData *td, FILE *f);
void PrintMaskSeqH1(unsigned char *seq, unsigned len, FILE *f);
char *TreeToString(struct tree *root, char *buf, struct TreeData *td);
char *StringToTree(char *cp, struct tree *root, struct TreeData *td);
struct tree* CopyTree(struct tree *root, struct TreeData *td, int saveSeq);
void DestroyTree(struct tree *root);
void OrderColumns(struct TreeData *td);
void Transpose_b1(unsigned char *data, unsigned w, unsigned h);
unsigned ReverseMapTaxon(unsigned label, struct TreeData *td);
unsigned long CountSubTrees(unsigned numTaxaOnTree, unsigned numTaxa);
void PackVertical(struct TreeData *td);
unsigned PackVerticalSeq(unsigned char *in, unsigned *out, unsigned numChars);
void PrintSeqVS(unsigned *seq, unsigned len, FILE *f);

void PartitionColumns(struct TreeData *td);
unsigned ExhaustiveScorePartition(unsigned char *part, unsigned partWidth, unsigned partMask, unsigned width, unsigned height, unsigned *weights, unsigned *pBestPartMask, unsigned depth);
unsigned ScorePartition(unsigned char *part, unsigned partWidth, unsigned partMask, unsigned width, unsigned height, unsigned *weights);
//void MutatePartition(unsigned char *part, unsigned partWidth, unsigned remainingWidth, unsigned width, unsigned height, unsigned *weights);
void MutatePartition(unsigned char *part, unsigned partWidth, unsigned remainingWidth, unsigned width, unsigned height, unsigned *weights, unsigned *pRIn, unsigned *pROut);
void KA2MutatePartition(unsigned char *part, unsigned partWidth, unsigned remainingWidth, unsigned width, unsigned height, unsigned *weights, unsigned *pRIn, unsigned *pROut, unsigned *upperBoundWeights);		// WTJW 11/6/2004
void SwapSites(unsigned char *data, unsigned a, unsigned b, unsigned width, unsigned height, unsigned *weights);
unsigned TurnOnBits(unsigned bits, unsigned n, unsigned *indices);
unsigned CountOnBits(unsigned x);

char GetBaseFromMask(unsigned char mask);
unsigned char GetMaskFromBase(char base);
int GetCodeFromBase(char base);

#ifdef DEBUG
void DebugPrintSwitches(void);
#endif	// DEBUG

// Global Variables
extern struct CharData gCD;
extern struct TreeData gTD;
extern FILE *dbgfile;
extern unsigned _ColumnCompareWidth;			// Used by qsort() callbacks in at least 2 modules.

// Static (and hopefully inlined) functions
// Safely divide two numbers, rounding up the answer.  Handles the x == 0 case.
static INLINE unsigned UDIVROUNDUP(unsigned x, unsigned div) {
	//return (x - !!x) / div + !!x;				// Works, but my hunch is the code below is faster.
	return (!x - 1) & ((x - 1) / div + 1);
}

static INLINE unsigned UROUNDUP(unsigned x, unsigned toNearest) {
	return UDIVROUNDUP(x, toNearest) * toNearest;
}

static INLINE unsigned GetNumBlocksFromSeqLen(unsigned seqLen) {
	//return (!seqLen - 1) & ((seqLen - 1) / (SITESPERBYTE * BYTESPERBLOCK) + 1);		// Safely round up
	return UDIVROUNDUP(seqLen, SITESPERBYTE * BYTESPERBLOCK);
}

#ifdef DEBUG
// Used for reporting the smallest and largest tree encountered during each updateIntervalSecs period.
extern unsigned branchAndBoundMinNumTaxa;
extern unsigned branchAndBoundMaxNumTaxa;
#endif	// DEBUG
#endif	// #include guard
