#include "common.h"
#include "seq.h"
#include "yanbader.h"
#include <dbgprint.h>
#include "arena.h"
#include <limits.h>

// Although we only need a single sequence triple per node rather than an array of them, rather than define a new and marginally different struct, we just use struct tree, with
// seqsByTreeSize[] containing a single element.

// NOTE: Although ImproveInitialUpperBound() must be called BEFORE BlockifyTaxa(), TreeBisectionReconnection() must be called AFTER this as it depends
// on sequence data being in the chosen high-performance format.

#define ICHILD(t) (LEFT + (RIGHT - LEFT) * ((t) == (t)->nodes[PARENT]->nodes[RIGHT]))
#define ISIBLING(t) (LEFT + (RIGHT - LEFT) * ((t) == (t)->nodes[PARENT]->nodes[LEFT]))
#define SIBLING(t) ((t)->nodes[PARENT]->nodes[ISIBLING(t)])

struct tbr_context {
	struct tree *root;
	struct tree *cut;			// The lower node on the edge that is to be deleted
	struct tree *curBelow;		// The thing currently chosen from below the cut edge
	unsigned char *curBelowSeq;	// The sequence for the thing currently chosen from below the cut edge, to be Fitched with the corresponding sequence from the thing above
	struct tree *bestCut;		// The cut edge resulting in the best combination so far
	struct tree *bestBelow;		// The thing chosen from below the cut edge in the best combination so far
	struct tree *bestAbove;		// The thing chosen from above the cut edge in the best combination so far
	unsigned nTaxa;				// Take from td->numTaxa
	int bestScoreImprovement;			// The score of the original cut edge minus the score of the best combination so far.  Signed because initially -1.
	unsigned oldEdgeScore;		// The score of the original edge.  Can be used with STOPEARLY.
	unsigned *weights;			// Taken from td->weights
	unsigned seqLenBlocks;		// Taken from td->seqLenBlocks
#ifdef DEBUG
	struct TreeData *td;		// You might think we need this for td->weights and td->seqLenBlocks, but actually we copied those -- now its only needed for PrintTree().
#endif	// DEBUG
	struct arena tbrMem;		// Need to make deep copies of all sequences to prevent changes to one changing another
	struct tree *elbow;			// Set by CutAt().
	struct tree *sibling;		// Set by CutAt();
	unsigned char *oldSibDown;		// Set by CutAt().  Points to the original sibling->seqsByTreeSize[ctx->nTaxa].seqs[DOWN].
	unsigned char *oldSibMidpoint;	// Set by CutAt().  Points to the original sibling->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT].
	unsigned char *oldCutLeftDown;	// Set by CutAt().  Iff the lower subtree is not a leaf, points to the original cut->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN].
	unsigned char *oldCutRightDown;	// Set by CutAt().  Iff the lower subtree is not a leaf, points to the original cut->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN].

	// After CutAt(), the following all point to distinct buffers that can be reused:
	// cut->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];
	// elbow->seqsByTreeSize[ctx->nTaxa].seqs[UP];
	// elbow->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];
	// cut->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];		// Iff the lower tree is not a leaf
	// cut->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];		// Iff the lower tree is not a leaf

	// The following are only set before performing a final chosen cut.  That's because if RotateTree() has been called, they're not reachable via cut anymore.
	unsigned char *oldCutMidPoint;
	unsigned char *oldCutLeftMidpoint;		// Iff the lower tree is not a leaf
	unsigned char *oldCutRighttMidpoint;		// Iff the lower tree is not a leaf

	// And the following are still available even after RotateTree(), so we don't need to store them separately here:
	//unsigned char *elbow->seqsByTreeSize[ctx->nTaxa].seqs[UP];
	//unsigned char *elbow->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];
};

//DEBUG
int VerifyTreeMIDPOINTSeqs(struct tree *t, struct tbr_context *ctx, int *piEdge) {
	int ok = 1;

	if (!t) {
		return 1;
	}

	if (t->nodes[PARENT]) {
		void *oldPos = ArenaGetPos(&ctx->tbrMem);
		unsigned char *tmp = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);

		//DBGPRINT2("VerifyTreeMIDPOINTSeqs(): Testing edge %d...\n", *piEdge);

		CHOSEN_FitchBases(
			t->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			tmp,
			ctx->seqLenBlocks
		);

		if (memcmp(tmp, t->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT], ctx->seqLenBlocks * BYTESPERBLOCK)) {
			DBGPRINT2("VerifyTreeMIDPOINTSeqs(): Edge %d has a bad MIDPOINT seq!\n", *piEdge);
			ok = 0;
			assert(0);
		}

		ArenaSetPos(&ctx->tbrMem, oldPos);

		++*piEdge;
	}

	ok &= VerifyTreeMIDPOINTSeqs(t->nodes[LEFT], ctx, piEdge);
	ok &= VerifyTreeMIDPOINTSeqs(t->nodes[RIGHT], ctx, piEdge);

	return ok;
}

// DAMMIT, we need to do postorder and that messes up the edge order...
int VerifyTreeUPSeqs(struct tree *t, struct tbr_context *ctx, int *piEdge) {
	int ok = 1;

	if (!t) {
		return 1;
	}

	ok &= VerifyTreeUPSeqs(t->nodes[LEFT], ctx, piEdge);
	ok &= VerifyTreeUPSeqs(t->nodes[RIGHT], ctx, piEdge);

	if (t->nodes[LEFT]) {
		// Count the edge leading from the root but don't test it -- we can't, as the result isn't stored anywhere!
		if (t->nodes[RIGHT]) {
			void *oldPos = ArenaGetPos(&ctx->tbrMem);
			unsigned char *tmp = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);

			//DBGPRINT2("VerifyTreeUPSeqs(): Testing postorder node %d...\n", *piEdge);

			CHOSEN_FitchBases(
				t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
				t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
				tmp,
				ctx->seqLenBlocks
			);

			if (memcmp(tmp, t->seqsByTreeSize[ctx->nTaxa].seqs[UP], ctx->seqLenBlocks * BYTESPERBLOCK)) {
				DBGPRINT2("VerifyTreeUPSeqs(): Postorder node %d has a bad UP seq!\n", *piEdge);
				ok = 0;
				assert(0);
			}

			ArenaSetPos(&ctx->tbrMem, oldPos);
		}

		++*piEdge;
	}

	return ok;
}

// DAMMIT, we need to do postorder and that messes up the edge order...
int VerifyTreeDOWNSeqs(struct tree *t, struct tbr_context *ctx, int *piEdge) {
	int ok = 1;

	if (!t) {
		return 1;
	}

	if (t->nodes[PARENT]) {
		//DBGPRINT2("VerifyTreeDOWNSeqs(): Testing edge %d...\n", *piEdge);

		if (!t->nodes[PARENT]->nodes[RIGHT]) {
			// Root's only child.  Should be identical to root's UP sequence.
			assert(t->nodes[PARENT] == ctx->root);
			if (memcmp(t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN], ctx->root->seqsByTreeSize[ctx->nTaxa].seqs[UP], ctx->seqLenBlocks * BYTESPERBLOCK)) {
				DBGPRINT2("VerifyTreeDOWNSeqs(): Edge %d (from root to its child) has a bad DOWN seq!\n", *piEdge);
				ok = 0;
				assert(0);
			}
		} else {
			// Regular edge
			void *oldPos = ArenaGetPos(&ctx->tbrMem);
			unsigned char *tmp = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);

			CHOSEN_FitchBases(
				t->nodes[PARENT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
				SIBLING(t)->seqsByTreeSize[ctx->nTaxa].seqs[UP],
				tmp,
				ctx->seqLenBlocks
			);

			if (memcmp(tmp, t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN], ctx->seqLenBlocks * BYTESPERBLOCK)) {
				DBGPRINT2("VerifyTreeDOWNSeqs(): Edge %d has a bad DOWN seq!\n", *piEdge);
				ok = 0;
				assert(0);
			}

			ArenaSetPos(&ctx->tbrMem, oldPos);
		}

		++*piEdge;
	}

	ok &= VerifyTreeDOWNSeqs(t->nodes[LEFT], ctx, piEdge);
	ok &= VerifyTreeDOWNSeqs(t->nodes[RIGHT], ctx, piEdge);

	return ok;
}

// All piEdge does is number the edges.  It needs to be point at 0 on the 1st call.
// Actually piEdge is totally ignored!
int VerifyTreeSeqs(struct tree *t, struct tbr_context *ctx, int *piEdge) {
	int actualEdgeCount = 0;
	int ok = 1;

	ok = VerifyTreeMIDPOINTSeqs(t, ctx, &actualEdgeCount);
	actualEdgeCount = 0;
	ok &= VerifyTreeUPSeqs(t, ctx, &actualEdgeCount);
	actualEdgeCount = 0;
	//DEBUG: VerifyTreeDOWNSeqs() doesn't work for bi-rooted trees yet so temporarily turn it off so RerootTreeAbove() can work...
	//ok &= VerifyTreeDOWNSeqs(t, ctx, &actualEdgeCount);
	return ok;
}

// I'm kinda surprised that I can't calculate this en passant while doing the incremental updates...  But I'm pretty sure I can't.
unsigned FitchScoreTree(struct tree *t, struct tbr_context *ctx) {
	//static unsigned char testBuf[5000];			//DEBUG
	unsigned score;		//DEBUG

	assert(t);

	if (t->nodes[LEFT]) {
		if (t->nodes[RIGHT]) {
			// Regular non-root node
			return (score = FitchScoreTree(t->nodes[LEFT], ctx) +
				FitchScoreTree(t->nodes[RIGHT], ctx) +
				CHOSEN_FitchScore(
					t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
					t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
					ctx->weights,
					ctx->seqLenBlocks,
					INT_MAX
				), printf("<%u>", score), score);
		} else {
			// Root
			return (score = FitchScoreTree(t->nodes[LEFT], ctx) +
				CHOSEN_FitchScore(
					t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
					t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],		// The root's fixed sequence is actually stored here
					ctx->weights,
					ctx->seqLenBlocks,
					INT_MAX
				), printf("R<%u>", score), score);
		}
	} else {
		// Leaf
		printf("L<%u>", 0);
		return 0;
	}
}

// Rotate a bi-rooted tree so that t has t->nodes[iFrom] as a parent.  In the usual case, this means swapping its current
// parent with that child and leaving the other child alone; in the case of the root, the root node is actually
// removed from the tree and should be freed by the caller.  Returns the node that should be linked as a child to the 
// original t->nodes[iFrom].  This is a tricky dance that requires a lot of care in the order that updates are made.
// Also note that because after rotation every sequence is still pointed to by exactly one node, we are allowed
// to shallow-copy them -- there's no risk of creating pointer aliasing.
struct tree *RotateTree(struct tree *t, int iFrom, struct tbr_context *ctx) {
	int iOther = LEFT + (RIGHT - LEFT) * (iFrom == LEFT);
	struct tree *from = t->nodes[iFrom];
	struct tree *other = t->nodes[iOther];

	if (!t->nodes[PARENT]) {
		struct tree *oldRoot = t;

		// The root is handled specially -- it actually disappears.
		assert(other);					// We must be a bi-rooted tree.
		other->nodes[PARENT] = from;

		other->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = from->seqsByTreeSize[ctx->nTaxa].seqs[UP];		// UP doesn't change
		other->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = oldRoot->seqsByTreeSize[ctx->nTaxa].seqs[UP];			// Think about it
		//free(oldRoot);			// Caller should free it
		return other;
	} else {
		// A regular node
		struct tree *rest = RotateTree(t->nodes[PARENT], ICHILD(t), ctx);		// Must pass the original "version" of this node to the recursive instance
		t->nodes[PARENT] = from;
		t->nodes[iFrom] = rest;			// Other doesn't change

		t->seqsByTreeSize[ctx->nTaxa].seqs[UP] = from->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];
		t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = from->seqsByTreeSize[ctx->nTaxa].seqs[UP];
		t->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = from->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];
		return t;
	}
}

// Slightly simpler than UpdateAboveTbr() because we know that we will always be proceeding downwards through the tree,
// updating only DOWN and MIDPOINT sequences.
unsigned UpdateBelowTbr(struct tree *t, struct tbr_context *ctx) {
	unsigned score = 0;

	assert(t);			// We should never call ourselves with a NULL node

	if (t->nodes[LEFT]) {
		// Non-leaf whose DOWN sequence was just updated.
		assert(t->nodes[RIGHT]);		// If a non-root node has at least one child, it must have 2.

		// Update LEFT's DOWN sequence.
		CHOSEN_FitchBases(
			t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			ctx->seqLenBlocks
		);

		//HACK: It would be nicer and faster if this could be combined with the previous into a single "FitchBasesAndScore()" function...
		score += CHOSEN_FitchScore(
			t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			ctx->weights,
			ctx->seqLenBlocks,
			INT_MAX
		);

		// Update LEFT's MIDPOINT sequence.
		CHOSEN_FitchBases(
			t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT],
			ctx->seqLenBlocks
		);

		score += UpdateBelowTbr(t->nodes[LEFT], ctx);			// Recurse down & left

		// Update RIGHT's DOWN sequence.
		CHOSEN_FitchBases(
			t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			ctx->seqLenBlocks
		);

		//HACK: It would be nicer and faster if this could be combined with the previous into a single "FitchBasesAndScore()" function...
		score += CHOSEN_FitchScore(
			t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			ctx->weights,
			ctx->seqLenBlocks,
			INT_MAX
		);

		// Update RIGHT's MIDPOINT sequence.
		CHOSEN_FitchBases(
			t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT],
			ctx->seqLenBlocks
		);

		score += UpdateBelowTbr(t->nodes[RIGHT], ctx);			// Recurse down & right
	}

	return score;
}

// This is designed to be called from any internal node on a tree whose root consists of a single
// taxon with a LEFT child but no RIGHT child (such as trees produced by ConstructBaseTree(), BuildGreedyTree() or BranchAndBound()).
// It can also be called from a leaf by specifying either LEFT or RIGHT as the edge just updated.
// It won't work properly if called on a tree with a "binary" root (2 children).
// iFrom is one of LEFT or RIGHT, indicating the edge directed towards t that has just been updated.
unsigned UpdateAboveTbr(struct tree *t, int iFrom, struct tbr_context *ctx) {
	unsigned score = 0;
	int iSibling;

	assert(t);		// We should never call ourselves recursively with a NULL node.
	assert(iFrom != PARENT);		// Should now be handled by UpdateBelowTbr();

	// Nothing to do at the root
	if (!t->nodes[PARENT]) {
		return 0;
	}

	iSibling = LEFT + (RIGHT - LEFT) * (iFrom == LEFT);

	// If we've been called from a non-leaf node, we need to go down the other side.
	if (t->nodes[iSibling]) {
		// Update sibling's DOWN sequence.
		CHOSEN_FitchBases(
			t->nodes[iFrom]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			t->nodes[iSibling]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			ctx->seqLenBlocks
		);

		//HACK: It would be nicer and faster if this could be combined with the previous into a single "FitchBasesAndScore()" function...
		score += CHOSEN_FitchScore(
			t->nodes[iFrom]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			ctx->weights,
			ctx->seqLenBlocks,
			INT_MAX
		);

		// Update sibling's MIDPOINT sequence.
		CHOSEN_FitchBases(
			t->nodes[iSibling]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
			t->nodes[iSibling]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
			t->nodes[iSibling]->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT],
			ctx->seqLenBlocks
		);

		score += UpdateBelowTbr(t->nodes[iSibling], ctx);			// Recurse down
	}

	// Update our UP sequence.
	CHOSEN_FitchBases(
		t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		t->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		ctx->seqLenBlocks
	);

	//HACK: It would be nicer and faster if this could be combined with the previous into a single "FitchBasesAndScore()" function...
	score += CHOSEN_FitchScore(
		t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		ctx->weights,
		ctx->seqLenBlocks,
		INT_MAX
	);

	// Update our MIDPOINT sequence.
	CHOSEN_FitchBases(
		t->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
		t->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT],
		ctx->seqLenBlocks
	);

	score += UpdateAboveTbr(t->nodes[PARENT], LEFT + (RIGHT - LEFT) * (t == t->nodes[PARENT]->nodes[RIGHT]), ctx);				// Recurse up

	return score;
}

void EnumThingsAbove(struct tree *t, struct tbr_context *ctx) {
	unsigned score;

	if (!t) {
		return;
	}

	// Now that we actually physically cut the edge, there's no special case for the elbow on the upper subtree.
	// Just an ordinary edge.  Use the MIDPOINT sequence.
	score = CHOSEN_FitchScore(
		ctx->curBelowSeq,
		t->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT],
		ctx->weights,
		ctx->seqLenBlocks,
		//ctx->oldEdgeScore		// Possibly unused.  Will be interpreted as unsigned, but that only a makes a rare case slightly slower
		INT_MAX
	);

	DBGPRINT4("EnumThingsAbove(): Found new edge from lower %s to upper edge with weight %u (improvement: %d).\n", ctx->curBelow == ctx->cut ? "root" : "edge", score, (int) ctx->oldEdgeScore - (int) score);
	if ((int) ctx->oldEdgeScore - (int) score > ctx->bestScoreImprovement) {
		DBGPRINT2("EnumThingsAbove(): Found new best-improvement edge with improvement %d!\n", (int) ctx->oldEdgeScore - (int) score);
		ctx->bestCut = ctx->cut;
		ctx->bestBelow = ctx->curBelow;
		ctx->bestAbove = t;
		ctx->bestScoreImprovement = (int) ctx->oldEdgeScore - (int) score;
	}

	EnumThingsAbove(t->nodes[LEFT], ctx);
	EnumThingsAbove(t->nodes[RIGHT], ctx);
}

// Enumerate every "thing" in the subtree below the cut edge.  A "thing" is either the root of this subtree (which may be a single node)
// or an edge.  Every thing enumerated is passed to EnumThingsBelow(), which compares it with every "thing" in the subtree above the cut edge.
void EnumThingsBelow(struct tree *t, struct tbr_context *ctx) {
	if (!t) {
		return;
	}

	if (t == ctx->cut) {
		// Special case: we are the binary root of the lower subtree.  (We may or may not be a leaf.)
		// Either way, we want to use the UP sequence.  Also, we want to skip regular processing of our
		// 2 edges and proceed straight to having our grandchildren each examine their parental edge.
		ctx->curBelow = t;
		ctx->curBelowSeq = ctx->curBelow->seqsByTreeSize[ctx->nTaxa].seqs[UP];
		EnumThingsAbove(ctx->root->nodes[LEFT], ctx);			// Compare this with every "thing" in the subtree above the cut.

		if (t->nodes[LEFT]) {
			EnumThingsBelow(t->nodes[LEFT]->nodes[LEFT], ctx);
			EnumThingsBelow(t->nodes[LEFT]->nodes[RIGHT], ctx);
			assert(t->nodes[RIGHT]);		// If one child is present, both should be
			EnumThingsBelow(t->nodes[RIGHT]->nodes[LEFT], ctx);
			EnumThingsBelow(t->nodes[RIGHT]->nodes[RIGHT], ctx);
		}
	} else {
		// Just an ordinary edge.  Use the MID sequence.
		ctx->curBelow = t;
		ctx->curBelowSeq = ctx->curBelow->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];
		EnumThingsAbove(ctx->root->nodes[LEFT], ctx);			// Compare this with every "thing" in the subtree above the cut.

		EnumThingsBelow(t->nodes[LEFT], ctx);
		EnumThingsBelow(t->nodes[RIGHT], ctx);
	}
}

//DEBUG
static int DEBUG_skipEdgesCount = 0;
static unsigned char CHECKcutUP[5000];
static unsigned char CHECKcutDOWN[5000];
static unsigned char CHECKcutMIDPOINT[5000];
static unsigned char CHECKelbowUP[5000];
static unsigned char CHECKelbowDOWN[5000];
static unsigned char CHECKelbowMIDPOINT[5000];

//void CutAt(struct tree *t, struct tbr_context *ctx) {
//	unsigned belowScore;
//	unsigned lowAboveScore;
//	unsigned highAboveScore;
//	unsigned fusedEdgeScore;
//	struct tree *sibling;
//	struct tree *elbow;
//#ifdef DEBUG
//	int edgeCount = 0;
//#endif	// DEBUG
//	// Don't consider deleting the edge leading to the root taxon, since that would lead to exploring trees
//	// that we already explore in other ways.  (Every tree can be "picked up" by the root taxon.)
//	//HACK: I suppose we could store these scores for each edge, rather than recomputing them, but
//	// I doubt that will give much speedup since for every edge cut, we computing O(nTaxa^2) Fitch scores.
//	ctx->oldEdgeScore = CHOSEN_FitchScore(
//		t->seqsByTreeSize[ctx->nTaxa].seqs[UP],
//		t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
//		ctx->weights,
//		ctx->seqLenBlocks,
//		INT_MAX				//HACK: Should probably be UINT_MAX, but let's do this just in case ..._FitchScore() has a bug...
//	);
//
//	DBGPRINT2("EnumCuts(): Cutting an edge with weight %u.\n", ctx->oldEdgeScore);
//#ifdef DEBUG
//	DBGPRINT1("Full tree: ");
//	PrintTree(ctx->root, ctx->td, dbgfile);
//	DBGPRINT1(";\n");
//
//	DBGPRINT1("Tree below this edge: ");
//	PrintTree(t, ctx->td, dbgfile);
//	DBGPRINT1(";\n");
//#endif	//DEBUG
//
//	DBGPRINT2("Original score of entire tree before cut: %u\n", FitchScoreTree(ctx->root, ctx));
//	ctx->cut = t;
//
//	assert(VerifyTreeSeqs(ctx->root, ctx, &edgeCount));
//
//	//DEBUG: Try skipping over the first few edge
//	if (DEBUG_skipEdgesCount) {
//		--DEBUG_skipEdgesCount;
//	} else {
//
//	elbow = t->nodes[PARENT];
//	sibling = SIBLING(t);
//
//#ifdef DEBUG
//	//DEBUG: Test an assumption...
//	memcpy(CHECKcutUP, t->seqsByTreeSize[ctx->nTaxa].seqs[UP], ctx->seqLenBlocks * BYTESPERBLOCK);
//	memcpy(CHECKcutDOWN, t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN], ctx->seqLenBlocks * BYTESPERBLOCK);
//	memcpy(CHECKcutMIDPOINT, t->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT], ctx->seqLenBlocks * BYTESPERBLOCK);
//	memcpy(CHECKelbowUP, elbow->seqsByTreeSize[ctx->nTaxa].seqs[UP], ctx->seqLenBlocks * BYTESPERBLOCK);
//	memcpy(CHECKelbowDOWN, elbow->seqsByTreeSize[ctx->nTaxa].seqs[DOWN], ctx->seqLenBlocks * BYTESPERBLOCK);
//	memcpy(CHECKelbowMIDPOINT, elbow->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT], ctx->seqLenBlocks * BYTESPERBLOCK);
//#endif	// DEBUG
//	// Join the sibling to its grandparent
//	sibling->nodes[PARENT] = t->nodes[PARENT]->nodes[PARENT];		// OK since our parent must have a parent
//	t->nodes[PARENT]->nodes[PARENT]->nodes[ICHILD(t->nodes[PARENT])] = sibling;
//
//	// Sibling inherits parent's DOWN seq; its UP seq is unchanged.  Its MIDPOINT is also known already.
//	//DEBUG: Let's try deep-copying these sequences so that UpdateBelowTbr() and UpdateAboveTbr() avoid scribbling over other nodes' sequences...
//	//sibling->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = t->nodes[PARENT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];
//	//sibling->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];		// Think about it...
//	memcpy(sibling->seqsByTreeSize[ctx->nTaxa].seqs[DOWN], t->nodes[PARENT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN], ctx->seqLenBlocks * BYTESPERBLOCK);
//	memcpy(sibling->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT], t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN], ctx->seqLenBlocks * BYTESPERBLOCK);
//
//	// Upper subtree = "fused" edge joining sibling to grandparent + subtree above the edge + subtree below the edge.
//	fusedEdgeScore = CHOSEN_FitchScore(
//		sibling->seqsByTreeSize[ctx->nTaxa].seqs[UP],
//		sibling->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
//		ctx->weights,
//		ctx->seqLenBlocks,
//		INT_MAX
//	);
//
//	DBGPRINT1("Tree above this edge: ");
//	PrintTree(ctx->root, ctx->td, dbgfile);
//	DBGPRINT1(";\n");
//
//	lowAboveScore = UpdateBelowTbr(sibling, ctx);
//	highAboveScore = UpdateAboveTbr(sibling->nodes[PARENT], ICHILD(sibling), ctx);
//
//	edgeCount = 0;
//	assert(VerifyTreeSeqs(ctx->root, ctx, &edgeCount));		// Verify the upper tree
//
//	if (t->nodes[LEFT]) {
//		assert(t->nodes[RIGHT]);
//		// We don't actually need the UP or MIDPOINT sequences for our LEFT and RIGHT children (and they don't have any
//		// sensible meaning anyway) so we don't bother making any changes here.  We need a special case anyway to skip
//		// over the immediate children of the root, so we might as well just leave the root "thing"'s relevant sequence
//		// in t->seqsByTreeSize[ctx->nTaxa].seqs[UP].
//		//t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP];
//		//t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP];
//		//t->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = t->seqsByTreeSize[ctx->nTaxa].seqs[UP];		// This just allows us to always compare with the MIDPOINT seq, avoiding a special case.
//
//		// Actually I think we do need this for UpdateBelowTbr() to work!  Deep-copy them for safety.  But don't worry about updating the MIDPOINT as I'm
//		// pretty sure we don't need that.
//		memcpy(t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN], t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP], ctx->seqLenBlocks * BYTESPERBLOCK);
//		memcpy(t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN], t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP], ctx->seqLenBlocks * BYTESPERBLOCK);
//		belowScore = UpdateBelowTbr(t->nodes[LEFT], ctx);
//		belowScore += UpdateBelowTbr(t->nodes[RIGHT], ctx);
//	} else {
//		belowScore = 0;
//	}
//
//	// Sever the lower subtree
//	t->nodes[PARENT] = NULL;
//
//	//belowScore = UpdateBelowTbr(t, ctx);
//
//	DBGPRINT5("Weight of upper subtree = %u + %u + %u = %u.\n", fusedEdgeScore, lowAboveScore, highAboveScore, fusedEdgeScore + lowAboveScore + highAboveScore);
//	DBGPRINT2("Weight of bottom subtree = %u.\n", belowScore);
//
//	// Verify the lower tree, in 2 halves (use a single edge count tho...)
//	//edgeCount = 0;
//	//assert(VerifyTreeSeqs(ctx->cut->nodes[LEFT], ctx, &edgeCount));
//	//assert(VerifyTreeSeqs(ctx->cut->nodes[RIGHT], ctx, &edgeCount));
//
//	//HACK: Now recalculate it... correctly... ;)
//	lowAboveScore = FitchScoreTree(ctx->root, ctx);		// Don't need highAboveScore anymore
//	belowScore = FitchScoreTree(t, ctx);
//	//DBGPRINT4("CORRECT weight of upper subtree = %u + %u = %u.\n", fusedEdgeScore, lowAboveScore, fusedEdgeScore + lowAboveScore);
//	DBGPRINT2("CORRECT weight of upper subtree = %u.\n", lowAboveScore);
//	DBGPRINT2("CORRECT weight of bottom subtree = %u.\n", belowScore);
//	}
//}

void CutAt(struct tree *t, struct tbr_context *ctx) {
	int iGrandchild;
#ifdef DEBUG
	int edgeCount = 0;		// Only used for VerifyTreeSeqs()
#endif	// DEBUG

	ctx->oldEdgeScore = CHOSEN_FitchScore(
		t->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
		ctx->weights,
		ctx->seqLenBlocks,
		INT_MAX				//HACK: Should probably be UINT_MAX, but let's do this just in case ..._FitchScore() has a bug...
	);

	DBGPRINT2("EnumCuts(): Cutting an edge with weight %u.\n", ctx->oldEdgeScore);
#ifdef DEBUG
	DBGPRINT1("Full tree: ");
	PrintTree(ctx->root, ctx->td, dbgfile);
	DBGPRINT1(";\n");

	DBGPRINT1("Tree below this edge: ");
	PrintTree(t, ctx->td, dbgfile);
	DBGPRINT1(";\n");
#endif	//DEBUG

	DBGPRINT2("Original score of entire tree before cut: %u\n", FitchScoreTree(ctx->root, ctx));

	assert(VerifyTreeSeqs(ctx->root, ctx, &edgeCount));

	ctx->cut = t;
	ctx->elbow = t->nodes[PARENT];
	ctx->sibling = SIBLING(t);
	iGrandchild = ICHILD(ctx->elbow);

	// Save pointers to buffers that will not be referenced by the new subtrees
	ctx->oldSibDown = ctx->sibling->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];
	ctx->oldSibMidpoint = ctx->sibling->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];
	
	if (t->nodes[LEFT]) {
		assert(t->nodes[RIGHT]);
		// The DOWN and MIDPOINT sequences of the cut's LEFT and RIGHT children are now also undefined.  We need to change
		// the pointers for the DOWN sequences, so we have to save the original buffer pointers explicitly, but we can
		// just leave the MIDPOINT pointers unchanged to remember those buffers.
		ctx->oldCutLeftDown = t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];
		ctx->oldCutRightDown = t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];

		// Update sequence buffers
		// Each DOWN sequence here must refer to its sibling's UP sequence, but note it does not OWN that buffer.
		t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[UP];
		t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[UP];

		UpdateBelowTbr(t->nodes[LEFT], ctx);
		UpdateBelowTbr(t->nodes[RIGHT], ctx);
	}

	// Rewire the tree
	ctx->sibling->nodes[PARENT] = ctx->elbow->nodes[PARENT];
	ctx->elbow->nodes[PARENT]->nodes[iGrandchild] = ctx->sibling;
	t->nodes[PARENT] = NULL;

	// Update sequence buffers
	// Sibling's UP doesn't change.
	ctx->sibling->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = ctx->elbow->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];
	ctx->sibling->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];		// Think about it

	// Elbow's DOWN sequence is now owned by sibling's DOWN.  Leave the other 2 pointers unchanged to track these buffers.
	// Cut's UP doesn't change.  Its DOWN and MIDPOINT are now undefined (but leave the pointers unchanged to keep track of the buffers)

	UpdateAboveTbr(ctx->sibling->nodes[PARENT], iGrandchild, ctx);
	UpdateBelowTbr(ctx->sibling, ctx);
}

// Glue the original edge back in.  All the necessary info has been saved in ctx.
void UnCut(struct tbr_context *ctx) {
	int iSibling = ICHILD(ctx->sibling);

	// Wire the tree back together.  The node pointers in elbow are unchanged so don't need setting.
	ctx->cut->nodes[PARENT] = ctx->elbow;
	ctx->sibling->nodes[PARENT]->nodes[iSibling] = ctx->elbow;
	ctx->sibling->nodes[PARENT] = ctx->elbow;

	// Set up the sequence buffers as they were before
	// Elbow's UP and MIDPOINT point to buffers which have not been touched, so leave them unchanged.
	ctx->elbow->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = ctx->sibling->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];

	// Cut's UP and MIDPOINT haven't changed (even though its MIDPOINT was technically undefined after the cut).
	ctx->cut->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = ctx->sibling->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];

	// Sibling's UP hasn't changed.
	ctx->sibling->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = ctx->oldSibDown;
	ctx->sibling->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = ctx->oldSibMidpoint;

	if (ctx->cut->nodes[LEFT]) {
		assert(ctx->cut->nodes[RIGHT]);

		// UP and MIDPOINT sequences haven't changed (even though MIDPOINT sequences were technically undefined after the cut).
		ctx->cut->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = ctx->oldCutLeftDown;
		ctx->cut->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = ctx->oldCutRightDown;

		UpdateBelowTbr(ctx->cut->nodes[LEFT], ctx);
		UpdateBelowTbr(ctx->cut->nodes[RIGHT], ctx);
	}

	UpdateBelowTbr(ctx->sibling, ctx);
	UpdateAboveTbr(ctx->elbow->nodes[PARENT], iSibling, ctx);

	assert(VerifyTree(ctx->root, NULL));
}

void EnumCuts(struct tree *t, struct tbr_context *ctx) {
	if (!t) {
		return;
	}

	if (t->nodes[PARENT] != ctx->root) {
		CutAt(t, ctx);

		// Enumerate all pairs of "things" from the two subtrees and find the best.
		EnumThingsBelow(t, ctx);

		UnCut(ctx);
	}

	EnumCuts(t->nodes[LEFT], ctx);
	EnumCuts(t->nodes[RIGHT], ctx);
}

// We need to be able to alter sequences without worrying that they are also pointed to by other nodes.
// We only actually do this for the sequence at position ctx->nTaxa, since that's the only sequence we care about.
// The tricky part is that not all sequences at all nodes are defined...
// Should only be run on mono-rooted trees.
void DeepCopyTreeSeqs(struct tree *t, struct tbr_context *ctx) {
	unsigned char *newSeq;
	int i;

	if (!t) {
		return;
	}
	
	// We don't copy any sequences at the root as they are all undefined.
	if (t->nodes[PARENT]) {
		// Copy UP, DOWN, and MIDPOINT sequences
		for (i = 0; i < 3; ++i) {
			// Don't copy anything at the root node -- all of its sequences are undefined.
			// At the root's child, deep-copy UP and MIDPOINT, and leave DOWN shallow-copied (this points to the taxon sequence for the root).
			// At other internal nodes, deep-copy all sequences.
			// At leaves, deep-copy DOWN and MIDPOINT sequences; UP sequences point to taxon sequences and can remain as shallow copies.
			if (t->nodes[PARENT]->nodes[PARENT] ? (t->nodes[LEFT] && t->nodes[RIGHT]) || (t->nodes[PARENT] && i != UP) : i != DOWN) {
				//DBGPRINT3("DeepCopyTreeSeqs(): Allocating copy buffer for %s for node %p.\n", i == UP ? "UP" : i == DOWN ? "DOWN" : "MIDPOINT", t);		//DEBUG
				newSeq = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);
				memcpy(newSeq, t->seqsByTreeSize[ctx->nTaxa].seqs[i], ctx->seqLenBlocks * BYTESPERBLOCK);
				t->seqsByTreeSize[ctx->nTaxa].seqs[i] = newSeq;
			}
		}
	}
	
	DeepCopyTreeSeqs(t->nodes[LEFT], ctx);
	DeepCopyTreeSeqs(t->nodes[RIGHT], ctx);
}

// Given a bi-rooted tree, reroot it by inserting a new bi-root above the given node.
// The root must be specified as ctx->root will normally refer to the root of the upper subtree, and we will
// be calling this function on the lower subtree.
struct tree *RerootTreeAbove(struct tree *t, struct tree *root, struct tbr_context *ctx) {
	struct tree *newRoot = malloc(sizeof (struct tree));
	int edgeCount = 0;
	int iChild = ICHILD(t);
	int iSibling = ISIBLING(t);

	newRoot->contents = -1;
	//TODO: Not only is the following being a massive leak, bear in mind that we don't allocate ANY room in the arena for storing struct edge_seq objects.
	newRoot->seqsByTreeSize = ArenaAllocate(&ctx->tbrMem, UROUNDUP(sizeof (struct edge_seqs) * (ctx->nTaxa + 1), BYTESPERBLOCK));		//TODO: Massive leak!
	newRoot->nodes[ICHILD(t)] = t;
	// Tricky!  The following line doesn't work because the compiler can evaluate the LHS of the assignment AFTER calling
	// RotateTree(), which results in ISIBLING(t) then equalling ICHILD(t) in the previous statement!
	//newRoot->nodes[ISIBLING(t)] = RotateTree(t->nodes[PARENT], ICHILD(t), ctx);
	newRoot->nodes[iSibling] = RotateTree(t->nodes[PARENT], ICHILD(t), ctx);
	newRoot->nodes[PARENT] = NULL;

	// One of these has been set by RotateTree() and needs correcting; the other just needs to be set.
	newRoot->nodes[LEFT]->nodes[PARENT] = newRoot;
	newRoot->nodes[RIGHT]->nodes[PARENT] = newRoot;

	// The only sequence we need to populate for the new root node is the UP sequence.  Coincidentally, that's
	// the only sequence allocated in the current root, so steal its seqsByTreeSize array and sequence memory.
	// Basically, the DOWN and MIDPOINT seqs of the original t are "spread" to the UP seqs of its new sibling and root, respectively.
	//TODO: Use the line below to get rid of the leak
	//newRoot->seqsByTreeSize = root->seqsByTreeSize;
	newRoot->nodes[iSibling]->seqsByTreeSize[ctx->nTaxa].seqs[UP] = t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];
	newRoot->seqsByTreeSize[ctx->nTaxa].seqs[UP] = t->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];
	free(root);

	assert(VerifyTree(newRoot, NULL));
	assert(VerifyTreeSeqs(newRoot, ctx, &edgeCount));
	return newRoot;
}

// Splice the bi-rooted tree rooted at t into the single-rooted host tree above the node host, by adding
// another node above host who gets t as a left child.  Similar to InsertNode() from common.c.
//TODO: Acquire the memory for seqsByTreeSize from the node in the upper subtree that was deleted.  Must do this to avoid a massive memory leak!
//TODO: Although CutAt() will never cut off the root node, it's still possible to join the lower subtree immediately below the root,
// which necessitates some special cases that we aren't handling yet.
void SpliceIn(struct tree *t, struct tree *host, struct tbr_context *ctx) {
	struct tree *newNode = malloc(sizeof (struct tree));
	int iSibling = ICHILD(host);
	int edgeCount = 0;

	newNode->contents = -1;
	newNode->seqsByTreeSize = ArenaAllocate(&ctx->tbrMem, UROUNDUP(sizeof (struct edge_seqs) * (ctx->nTaxa + 1), BYTESPERBLOCK));		//TODO: Scavenge existing mem
	newNode->nodes[LEFT] = t;
	newNode->nodes[RIGHT] = host;
	newNode->nodes[PARENT] = host->nodes[PARENT];
	host->nodes[PARENT]->nodes[iSibling] = newNode;		// Tricky: Need to calculate iSibling earlier to avoid a "race condition" where the index in the LHS is evaluated AFTER the RHS!  Damn sequence points...
	host->nodes[PARENT] = newNode;
	t->nodes[PARENT] = newNode;

	newNode->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = host->seqsByTreeSize[ctx->nTaxa].seqs[DOWN];
	newNode->seqsByTreeSize[ctx->nTaxa].seqs[UP] = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);		//TODO: Scavenge existing mem
	CHOSEN_FitchBases(
		t->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		host->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		newNode->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		ctx->seqLenBlocks
	);

	newNode->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);		//TODO: Scavenge existing mem
	CHOSEN_FitchBases(
		newNode->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		newNode->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
		newNode->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT],
		ctx->seqLenBlocks
	);

	// t currently only has an UP sequence allocated.  Need to borrow or allocate memory for the other seqs.
	t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = host->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT];		// Think about it
	t->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);		//TODO: Scavenge existing mem
	CHOSEN_FitchBases(
		t->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		t->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
		t->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT],
		ctx->seqLenBlocks
	);

	if (t->nodes[LEFT]) {
		// The children of a bi-root also only have an UP sequence allocated.  Need to allocate memory for the other seqs.
		assert(t->nodes[RIGHT]);
		t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);		//TODO: Scavenge existing mem
		t->nodes[LEFT]->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);		//TODO: Scavenge existing mem
		t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);		//TODO: Scavenge existing mem
		t->nodes[RIGHT]->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);		//TODO: Scavenge existing mem
		UpdateBelowTbr(t, ctx);
	}

	UpdateAboveTbr(newNode->nodes[PARENT], ISIBLING(newNode->nodes[PARENT]), ctx);

	// We have pointed other nodes at these sequences so we need to get fresh memory for them to avoid aliasing.
	host->seqsByTreeSize[ctx->nTaxa].seqs[DOWN] = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);		//TODO: Scavenge existing mem
	host->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT] = ArenaAllocate(&ctx->tbrMem, ctx->seqLenBlocks * BYTESPERBLOCK);		//TODO: Scavenge existing mem
	CHOSEN_FitchBases(
		t->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		newNode->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
		host->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
		ctx->seqLenBlocks
	);
	CHOSEN_FitchBases(
		host->seqsByTreeSize[ctx->nTaxa].seqs[UP],
		host->seqsByTreeSize[ctx->nTaxa].seqs[DOWN],
		host->seqsByTreeSize[ctx->nTaxa].seqs[MIDPOINT],
		ctx->seqLenBlocks
	);

	UpdateBelowTbr(host, ctx);

	assert(VerifyTree(ctx->root, NULL));
	assert(VerifyTreeSeqs(ctx->root, ctx, &edgeCount));
}

//TODO: The Enum...() functions all expect seqsByTreeSize[] in every tree node to contain the relevant sequence
// in the first position, whereas in fact this will be at position td->numTaxa.  Could either look this up inside every
// call, or alternatively (for a possibly very modest speed gain) temporarily "slide" all of these arrays forward
// by td->numTaxa positions before calling them, then "slide" them all back.
void TreeBisectionReconnection(struct tree *root, struct TreeData *td) {
	struct tbr_context ctx;
	int edgeCount = 0;		// Only used for VerifyTreeSeqs()
	struct tree *elbow;
	struct tree *sibling;
	struct tree *lowerRoot;
	const int slack = 1;	// The number of additional sequences our memory arena should accommodate.  At least 1 is needed for the Verify...Seqs() routines.

	ctx.bestCut = NULL;
	ctx.bestAbove = NULL;
	ctx.bestBelow = NULL;
	ctx.bestScoreImprovement = -1;		// Guarantees we will find the original edge placement at least
	ctx.root = root;
	ctx.nTaxa = td->numTaxa;
	ctx.weights = td->weights;
	ctx.seqLenBlocks = td->seqLenBlocks;
#ifdef DEBUG
	ctx.td = td;
#endif	// DEBUG

	// Because we need to modify sequences in an order that can't be predicted ahead of time, and because some sequences
	// may be aliased between nodes (e.g. ConstructBaseTree() does this for the root's UP sequence and its child's DOWN
	// sequence) we need to use our own arena allocator and deep-copy all the non-fixed sequences.
	//HACK: The precise number of sequences we allocate of course depends on implementation details of all the code we run.
	// But in general, DeepCopyTreeSeqs() sets the baseline.  It requires room for 3 sequences per internal node and 2 per leaf node (including the root).
	// We should only need a small fixed number of sequences on top of that for shuffling around or local operations.
	ArenaCreate(&ctx.tbrMem,
		//(ctx.nTaxa * 2 - 3) * ((ctx.nTaxa + 1) * sizeof (struct edge_seqs)) +			// Don't need to allocate this as we actually use what has already been allocated
		((ctx.nTaxa - 3) * 3 + ctx.nTaxa * 2 + slack) * (ctx.seqLenBlocks * BYTESPERBLOCK)		// For sequence data.
	);
	DeepCopyTreeSeqs(root, &ctx);

	DBGPRINT1("Finished DeepCopyTreeSeqs().\n");

	assert(VerifyTreeSeqs(root, &ctx, &edgeCount));

	EnumCuts(root->nodes[LEFT], &ctx);

	DBGPRINT3("TBR(): Old edge had weight %u; best TBR rearrangement found has score improvement %d.\n", ctx.oldEdgeScore, ctx.bestScoreImprovement);

	elbow = ctx.bestCut->nodes[PARENT];
	sibling = SIBLING(ctx.bestCut);
	CutAt(ctx.bestCut, &ctx);

#ifdef DEBUG
	DBGPRINT1("Lower subtree before rerooting: ");
	PrintTree(ctx.bestCut, ctx.td, dbgfile);
	fputc(';', dbgfile);
	fputc('\n', dbgfile);
#endif	// DEBUG

	// Rearrange the lower subtree if necessary
	if (ctx.bestCut == ctx.bestBelow) {
		DBGPRINT1("Best lower subtree is already rooted correctly...\n");
		lowerRoot = ctx.bestCut;
	} else {
		DBGPRINT1("Rerooting lower subtree...\n");
		lowerRoot = RerootTreeAbove(ctx.bestBelow, ctx.bestCut, &ctx);
	}

#ifdef DEBUG
	DBGPRINT1("Lower subtree after rerooting: ");
	PrintTree(lowerRoot, ctx.td, dbgfile);
	fputc(';', dbgfile);
	fputc('\n', dbgfile);
#endif	// DEBUG

	// Finish the job
	SpliceIn(lowerRoot, ctx.bestAbove, &ctx);

#ifdef DEBUG
	DBGPRINT1("After splicing lower subtree back in, entire tree is: ");
	PrintTree(ctx.root, ctx.td, dbgfile);
	fputc(';', dbgfile);
	fputc('\n', dbgfile);
#endif	// DEBUG

	DBGPRINT2("Score of new tree: %u\n", FitchScoreTree(ctx.root, &ctx));

	// This of course means that all sequence data in the resulting tree is no longer accessible...  But we don't actually care about that anyway.
	ArenaDestroy(&ctx.tbrMem);
}
