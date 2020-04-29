#include "common.h"
#include "greedytree.h"
#include "seq.h"
#include "yanbader.h"
#include <dbgprint.h>

struct greedy_context {
	unsigned bestScoreIncrease;
	struct tree *bestEdge;
};

//HACK: Most of this is heavily duplicated from BandBScanEdge().
// ctx->bestScoreIncrease should start very large; it will be decreased as positions are tried.
// Note that iTaxon is 0-based.
void FASTCALL GreedyFindPlaceForTaxon(struct tree *t, unsigned iTaxon, unsigned numTaxa, struct TreeData *td, struct greedy_context *ctx) {
	unsigned scoreIncrease;
	
	//TODO: Make separate measures for these
	//INCMEASURE(branchAndBoundScanEdgeCount);
	//INCMEASURE(branchAndBoundScanEdgeCountByDepth[numTaxa]);
	
	// Is there an edge?
	if (!t) return;

	// Compute the score that adding the next taxon to this edge would incur,
	// before actually adding the node to the tree.
	scoreIncrease = CHOSEN_FitchScore(
		t->seqsByTreeSize[numTaxa].seqs[MIDPOINT],
		td->charMat + iTaxon * td->seqLenBlocks * BYTESPERBLOCK,
		td->weights,
		td->seqLenBlocks,
		ctx->bestScoreIncrease		// Possibly unused.  Will be interpreted as unsigned, but that only a makes a rare case slightly slower
	);
	
	if (scoreIncrease < ctx->bestScoreIncrease) {
		ctx->bestEdge = t;
		ctx->bestScoreIncrease = scoreIncrease;
	}

	// Recurse
	//TODO: Transform one of these recursions into a loop (i.e. explicit tail-recursion)
	GreedyFindPlaceForTaxon(t->nodes[LEFT], iTaxon, numTaxa, td, ctx);
	GreedyFindPlaceForTaxon(t->nodes[RIGHT], iTaxon, numTaxa, td, ctx);
}

// Quickly build a tree by adding taxa in the order given.  (I.e. the 1st taxon added is taxon taxonOrder[0].)
// Each taxon is added in the place where it adds the least weight.
// pScore is where the score of the tree will be saved to.
// The taxon numbers in taxonOrder are 0-based.
struct tree *BuildGreedyTree(struct TreeData *td, int *taxonOrder, unsigned *pScore) {
	unsigned oldScore = td->score;			// Need to save this, since ConstructBaseTree() changes it
	struct tree *root;
	struct greedy_context ctx;
	int i;

	// ConstructBaseTree() adds to td->score.
	td->score = 0;
	root = ConstructBaseTree(td, taxonOrder[0] + 1, taxonOrder[1] + 1, taxonOrder[2] + 1);		//HACK: Silly ConstructBaseTree() wants 1-based taxon numbers...

	for (i = 3; i < td->numTaxa; ++i) {
		DBGPRINT3("BuildGreedyTree(): Weight after adding %u taxa is %u.\n", i, td->score);
		ctx.bestEdge = NULL;
		ctx.bestScoreIncrease = INT_MAX;
		GreedyFindPlaceForTaxon(root->nodes[LEFT], taxonOrder[i], i, td, &ctx);

		td->score += ctx.bestScoreIncrease;

		// It fits, so add it to the tree
		//HACK: Parameters needed by InsertNode() force clumsiness...
		InsertNode(ctx.bestEdge->nodes[PARENT], ctx.bestEdge == ctx.bestEdge->nodes[PARENT]->nodes[LEFT] ? LEFT : RIGHT, taxonOrder[i] + 1, i + 1, td);
		
		// Update the sequences on the tree
		UpdateTreeAfterInsertion(ctx.bestEdge->nodes[PARENT], i + 1, td);		// The new internal node is the new parent of the bestEdge node

#ifdef LOGTREES
		fprintf(td->treeLogFile, "%u\t%u\t%u\t%u\t%u\t", numTaxa + 1, td->score, td->restBound[numTaxa], td->score + td->restBound[numTaxa], td->bound);
#ifdef LOGTREESALL
		PrintTree(root, td, td->treeLogFile);
#else	// not LOGTREESALL
		fputc('-', td->treeLogFile);
#endif	// LOGTREESALL
		fputc('\n', td->treeLogFile);
#endif	// LOGTREES
	}

#ifdef DEBUG
	DBGPRINT2("BuildGreedyTree(): Found tree of weight %u: ", td->score);
	PrintTree(root, td, dbgfile);
	fputc(';', dbgfile);
	fputc('\n', dbgfile);
#endif	// DEBUG
	if (td->score < td->bound) {
		// Update the bound and clear out the current list of best trees
		td->bound = td->score;
		if (td->options & OPT_DISPLAYNEWBOUNDS) fprintf(stderr, "BuildGreedyTree(): New bound = %u\n", td->bound);
#ifdef TARGETMULTI
		DBGPRINT3("W%04d: BuildGreedyTree(): Found new UB of %u.\n", td->rank, td->bound);
#endif	// TARGETMULTI
	}

	*pScore = td->score;
	td->score = oldScore;
	return root;
}
