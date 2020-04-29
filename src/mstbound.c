#include "common.h"
#include <dbgprint.h>
#include "mstbound.h"
#include "seq.h"

// Used only by this module.  Since the function name is fairly general sounding but it
// could mean a variety of different things (e.g. is the data copied?), make it a static
// function of this module only.
// Specifically, it is called before BlockifyTaxa() is run, so the size of the charMat
// block is src->numTaxa * src->seqLen bytes.
static void CopyTreeData_b1_w1(struct TreeData *dest, struct TreeData *src) {
	// Populate new TreeData with relevant data
	dest->cd = src->cd;					// cd is not copied, as cd is not owned by the TreeData
	dest->bound = src->bound;
	dest->numTaxa = src->numTaxa;
	dest->seqLen = src->seqLen;
	dest->seqLenBlocks = src->seqLenBlocks;
	dest->charMat = (unsigned char *) malloc(dest->numTaxa * dest->seqLen);
	memcpy(dest->charMat, src->charMat, dest->numTaxa * dest->seqLen);
	dest->weights = (unsigned *) malloc(dest->seqLen * sizeof (unsigned));
	memcpy(dest->weights, src->weights, dest->seqLen * sizeof (unsigned));
	dest->taxonMap = (unsigned *) malloc(dest->numTaxa * sizeof (unsigned));
	memcpy(dest->taxonMap, src->taxonMap, dest->numTaxa * sizeof (unsigned));
	dest->restBound = (unsigned *) malloc(dest->numTaxa * sizeof (unsigned));
	memcpy(dest->restBound, src->restBound, dest->numTaxa * sizeof (unsigned));
#ifdef FASTLOWWEIGHTFITCH
	dest->byteWeights = (unsigned *) malloc(dest->seqLen);
	memcpy(dest->byteWeights, src->byteWeights, dest->seqLen);
#endif	// FASTLOWWEIGHTFITCH
	dest->bestList = NULL;				// NOTE: we don't copy the list of best trees
	dest->numTrees = 0;
	dest->optimalTreesMissed = 0;
	dest->maxTrees = src->maxTrees;
#ifdef LOGTREES
	dest->treeLogFile = NULL;
#endif	// LOGTREES
	if (src->workStack) {
		dest->workStack = (unsigned (*)[2]) malloc((src->numTaxa + 1) * 2 * sizeof (unsigned));
		memcpy(dest->workStack, src->workStack, (src->numTaxa + 1) * 2 * sizeof (unsigned));
	}
	dest->jobTreeSize = src->jobTreeSize;
	dest->rngSeed = src->rngSeed;
	dest->rand1 = src->rand1;			// Pushing the boundaries of pointlessness...  :)
	dest->rand2 = src->rand2;
	dest->rand3 = src->rand3;
#ifdef TARGETMULTI
	// It's almost certainly wrong to copy these across.
	dest->rank = -1;			// Indicate it's uninitialised.
	dest->nWorkers = 0;
	dest->mpiRequests = NULL;
	dest->nMpiRequests = 0;
	dest->pendingBoundBroadcast = 0;
	dest->receiveBuf = NULL;
	dest->stealBuf = NULL;
#endif	// TARGETMULTI
	dest->tempTreeFName = NULL;
	dest->tempTreeFile = NULL;
}

//HACK: I believe it is possible to violate the inequality w(SMT) >= w(MST)/2
// when some of the bases are ambiguous.  I'm sure I actually ran into this
// behaviour at some point and that I did some thinking about it and decided
// that this (i.e. ambiguous bases) was the reason.  Therefore this routine
// needs to be fixed by ignoring all sites containing ambiguous bases.
// WTJW 23/3/2005: Now ask CalcTreeMstWeight() to underestimate the MST weight
// whenever there are ambiguous bases, so this should be safe to use.
void InitMstBoundRest_b1(struct TreeData *td) {
	unsigned i;
	unsigned theMstBound;
	unsigned char *origCharMat;
	struct TreeData newTd;
	
	CopyTreeData_b1_w1(&newTd, td);
	origCharMat = newTd.charMat;
	
	for (i = 0; i < td->numTaxa - 1; ++i) {
		theMstBound = CalcTreeMstWeight_b1(&newTd, 1) / 2;		// Integer divide rounds down, as desired
		DBGPRINT3("theMstBound[%u] = <%u>\n", i, theMstBound);		//DEBUG
		if (theMstBound > td->restBound[i]) {
			DBGPRINT3("mstBound beats current lower bound for trees of first %u taxa by %u\n", i + 1, theMstBound - td->restBound[i]);
			td->restBound[i] = theMstBound;
		}
		
		// Form the new union of all taxa on the tree so far, and update the
		// character matrix pointer to point to this new "reduced" taxon set
		UnionOfSeqs(newTd.charMat, newTd.charMat + newTd.seqLen, newTd.charMat + newTd.seqLen, newTd.seqLen);
		newTd.charMat += newTd.seqLen;
		--newTd.numTaxa;
	}
	
	newTd.charMat = origCharMat;
	DestroyTreeData(&newTd);
}
