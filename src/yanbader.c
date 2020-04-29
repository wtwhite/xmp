#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "common.h"
#include <dbgprint.h>
#include "yanbader.h"
#include "seq.h"
#include "measure.h"
#ifdef TARGETMULTI
#include "mpiworker.h"
#include "mpi.h"
#endif	// TARGETMULTI

static void FASTCALL RecordNewTree(struct TreeData *td, struct tree *root);
static void FASTCALL UpdateTree(struct tree *node, struct tree *newNode, struct tree *oldNode, int seqChanged, int nTaxaOnTree, struct TreeData *td);

//TODO: InitTreeData() stuff:
// Set up weights[], lengths[], offsets[], constantWeights[], taxonSequences[][].

void FASTCALL BranchAndBound(struct tree *root, unsigned numTaxa, struct TreeData *td) {
	unsigned fastForwardIdx = INT_MAX;			// Assume we won't be fast-forwarding.
#ifdef TARGETMULTI
	int idx, flag;
	MPI_Status mpiStatus;
#endif	// TARGETMULTI
	
	if (td->updateProgressNow) {
		td->updateProgressNow = 0;
		EMMS;		// Allow crossing over from MMX mode to floating-point mode if MMXASMFITCH is in effect
#ifdef DEBUG
		fprintf(stderr, "%.2f%% complete.  %u trees having length %u found so far.  numTaxa=%u (%u <= numTaxa <= %u).  ", ((double) td->proportionComplete) / COMPLETE_UNITS * 100, td->numTrees, td->bound, numTaxa, branchAndBoundMinNumTaxa, branchAndBoundMaxNumTaxa);
		branchAndBoundMinNumTaxa = branchAndBoundMaxNumTaxa = numTaxa;
#else	// not DEBUG
		fprintf(stderr, "%.2f%% complete.  %u trees having length %u found so far.  ", ((double) td->proportionComplete) / COMPLETE_UNITS * 100, td->numTrees, td->bound);
#endif	// not DEBUG
		
		// Display an estimate of how much longer this will take
		EstimateTimeRemaining(td);
		
#ifdef MEASUREUPDATE
		ShowDepthsForPeriod(td);
#endif	// MEASUREUPDATE
	}
	//DBGPRINT1("Call to BranchAndBound\n");	//DEBUG
	
	if (numTaxa >= td->numTaxa) {
		return;
	}
	
	INCMEASURE(branchAndBoundCount);

#ifdef DEBUG
	if (numTaxa < branchAndBoundMinNumTaxa) {
		branchAndBoundMinNumTaxa = numTaxa;
	}
	if (numTaxa > branchAndBoundMaxNumTaxa) {
		branchAndBoundMaxNumTaxa = numTaxa;
	}
#endif	// DEBUG

	// Start scanning the tree
	// 3/3/2005: We now add in the td->constantWeights[numTaxa] value inside
	// BandBScanEdge(), at the same time as ScoreFitch() is computed, to
	// enable allowableScoreIncrease to be computed correctly.
	// (Alternatively, we could leave things as they were, and take off
	// (and then re-add) the td->constantWeights[numTaxa] value from td->bound.)
//	td->score += td->constantWeights[numTaxa];
////	if (numTaxa == 3) {
////		td->score += td->constantWeights[numTaxa];
////	}
////	td->score += td->constantWeights[numTaxa + 1];

	if (numTaxa > td->jobTreeSize) {
		// Keep track of which edges have yet to be processed
		td->workStack[numTaxa][0] = 0;
		td->workStack[numTaxa][1] = 2 * numTaxa - 3;
#ifdef TARGETMULTI
		// Process any steal requests or new bounds that have come up.
		flag = 0;
		// See the explanation in Worker() for this hack.
		if (td->nStealsFromUs < td->nWorkers) {
			MPI_Testany(td->nMpiRequests, td->mpiRequests, &idx, &flag, &mpiStatus);
		} else {
			// Hold your horses!  We need to wait for a send to complete before doing anything else.
			MPI_Testany(td->nMpiRequests - MAX_WCR, td->mpiRequests + MAX_WCR, &idx, &flag, &mpiStatus);
			idx += MAX_WCR;
		}
		if (flag) {
			DBGPRINT4("W%04d: Received message while working.  idx = %d, tag = %d.\n", td->rank, idx, mpiStatus.MPI_TAG);
			assert(!(idx == WCR_RCV_ALL && mpiStatus.MPI_TAG == MSG_NEWJOB));		// It should be impossible to receive a new job at this point.
			WorkerMain(td, root, idx, &mpiStatus, 0);
		}
#endif	// TARGETMULTI
	} else {
		// We may be fast-forwarding through some of these edges.
		// td->workStack[] already contains the pairs set up as they should be.
		fastForwardIdx = 0;
	}

	assert(td->workStack[numTaxa][0] <= 2 * numTaxa - 3);
	assert(td->workStack[numTaxa][0] <= 2 * numTaxa - 3);
	BandBScanEdge(root, root, LEFT, numTaxa, td, &fastForwardIdx);

	// Work-stealing may cause the 2nd value to have changed.
	assert(td->workStack[numTaxa][0] == td->workStack[numTaxa][1]);
//	td->score -= td->constantWeights[numTaxa];
////	td->score -= td->constantWeights[numTaxa + 1];
}

// td->score should already contain the non-PI weight of the first numTaxa + 1
// taxa.
void FASTCALL BandBScanEdge(struct tree *root, struct tree *parent, int iChild, unsigned numTaxa, struct TreeData *td, unsigned *pFastForwardIdx) {
	struct tree *child = parent->nodes[iChild];
	int scoreIncrease;
	int allowableScoreIncrease;			// Must be signed so that "if (scoreIncrease <= allowableScoreIncrease)" handles the case where the starting tree is above the initial bound
	unsigned i, proportion;
	
	INCMEASURE(branchAndBoundScanEdgeCount);
	INCMEASURE(branchAndBoundScanEdgeCountByDepth[numTaxa]);
	
	// Is there an edge?
	if (!child) return;

//	//DEBUG
//	{
//		unsigned trueScore;
//		
//		trueScore = Fitch(root, td);
//		if (trueScore + td->constantWeight != td->score) {
//			DBGPRINT4("ZOINKS!  td->score = %d, but trueScore + td->constantWeight = %d!  (%d taxa on tree)\n", td->score, trueScore + td->constantWeight, numTaxa);
//		}
//	}
	
	// Should we "fast-forward" through this edge in order to build an initial tree to work on from a stolen job?
	if (*pFastForwardIdx < td->workStack[numTaxa][0]) {
		++*pFastForwardIdx;

		// Recurse
		BandBScanEdge(root, child, LEFT, numTaxa, td, pFastForwardIdx);
		BandBScanEdge(root, child, RIGHT, numTaxa, td, pFastForwardIdx);
		return;
	}

	// Terminate early if others have stolen all the remaining subproblems for this tree.
	if (td->workStack[numTaxa][0] == td->workStack[numTaxa][1]) {
		return;
	}

	// Compute the score that adding the next taxon to this edge would incur,
	// before actually adding the node to the tree.
	allowableScoreIncrease = td->bound - td->score - td->restBound[numTaxa];			// Used by early-termination versions of ScoreFitch() function
	//scoreIncrease = td->constantWeights[numTaxa] + CHOSEN_FitchScore(
	//	td->taxonSequence[numTaxa] + dataOffset * BYTESPERBLOCK,
	//	parent->nodes[iChild]->sequence + (dataOffset + td->lengths[numTaxa] * 2) * BYTESPERBLOCK,
	//	td->weights + td->offsets[numTaxa] * BYTESPERBLOCK * SITESPERBYTE,
	//	td->lengths[numTaxa],
	//	allowableScoreIncrease		// Possibly unused.  Will be interpreted as unsigned, but that only a makes a rare case slightly slower
	//);
	scoreIncrease = CHOSEN_FitchScore(
		parent->nodes[iChild]->seqsByTreeSize[numTaxa].seqs[MIDPOINT],
		td->charMat + numTaxa * td->seqLenBlocks * BYTESPERBLOCK,
		td->weights,
		td->seqLenBlocks,
		allowableScoreIncrease		// Possibly unused.  Will be interpreted as unsigned, but that only a makes a rare case slightly slower
	);
	
	td->score += scoreIncrease;
	
#ifdef MEASURE
	if (numTaxa + 1 == td->numTaxa) {
		INCMEASURE(fullTreesConsideredCount);
	}
#endif	// MEASURE

	if (scoreIncrease <= allowableScoreIncrease) {
		void *arenaSavePoint = ArenaGetPos(&td->alignedMem);		// Enable "one big free()" at the end

		// It fits, so add it to the tree
		InsertNode(parent, iChild, numTaxa + 1, numTaxa + 1, td);
		
		// Update the sequences on the tree
		UpdateTreeAfterInsertion(parent->nodes[iChild], numTaxa + 1, td);

#ifdef LOGTREES
		fprintf(td->treeLogFile, "%u\t%u\t%u\t%u\t%u\t", numTaxa + 1, td->score, td->restBound[numTaxa], td->score + td->restBound[numTaxa], td->bound);
#ifdef LOGTREESALL
		PrintTree(root, td, td->treeLogFile);
#else	// not LOGTREESALL
		fputc('-', td->treeLogFile);
#endif	// LOGTREESALL
		fputc('\n', td->treeLogFile);
#endif	// LOGTREES
		
		if (numTaxa + 1 == td->numTaxa) {
			assert(td->score <= td->bound);			// Could only fail if we have calculated illegal restBound[] values
			if (td->score < td->bound) {
				// Update the bound and clear out the current list of best trees
				td->bound = td->score;
				if (td->options & OPT_DISPLAYNEWBOUNDS) fprintf(stderr, "New bound = %u\n", td->bound);
				DeleteOldTrees(td);
#ifdef TARGETMULTI
				DBGPRINT3("W%04d: Found new UB of %u.\n", td->rank, td->bound);
				td->pendingBoundBroadcast = 1;
				if (td->mpiRequests[WCR_SND_NEWBOUND] == MPI_REQUEST_NULL) {
					// It's safe to tell everyone now.
					SendNewBound(td);
				}
#endif	// TARGETMULTI
			}

			RecordNewTree(td, root);
		} else {
			// Tree is not yet complete: recurse
			BranchAndBound(root, numTaxa + 1, td);
			if (numTaxa + 1 == COUNT_UNIT) {
				++td->proportionComplete;
				++td->treesActuallyExamined;
			}
		}
		
		DeleteNode(parent, iChild);
		ArenaSetPos(&td->alignedMem, arenaSavePoint);
	} else {
		// How many full trees would have been generated by this subtree?
		// (Remember that we count any tree containing COUNT_UNIT taxa
		// as 1.)
		if (numTaxa + 1 <= COUNT_UNIT) {
			proportion = 1;
			for (i = numTaxa + 1; i < COUNT_UNIT; ++i) {
				proportion *= i * 2 - 3;
			}
			td->proportionComplete += proportion;
		}
	}
	
	td->score -= scoreIncrease;
	++td->workStack[numTaxa][0];		// The subproblem in which taxon numTaxa+1 was added to this edge has been solved.

	// Recurse
	//TODO: Transform one of these recursions into a loop (i.e. explicit tail-recursion)
	BandBScanEdge(root, child, LEFT, numTaxa, td, pFastForwardIdx);
	BandBScanEdge(root, child, RIGHT, numTaxa, td, pFastForwardIdx);
}

// Should be called after any operation that lowers td->bound -- either a tree we discovered ourselves,
// or a bound reported by another worker.
void FASTCALL DeleteOldTrees(struct TreeData *td) {
	//TreeList_Destroy(&td->bestList);
	if (td->tempTreeFile) {
		fclose(td->tempTreeFile);
		td->tempTreeFile = NULL;
		remove(td->tempTreeFName);
	}

	td->numTrees = 0;
	td->optimalTreesMissed = 0;		// We start afresh
}

void FASTCALL RecordNewTree(struct TreeData *td, struct tree *root) {
	if (td->options & OPT_DISPLAYNEWTREES) {
		fprintf(stderr, "Found tree ");
		PrintTree(root, td, stderr);
#ifdef MEASURE
		// Only reports sensible numbers of trees that fit into an unsigned int; always works for up to 12 taxa.
		fprintf(stderr, " with score %u.  %u full trees considered so far.\n", td->score, fullTreesConsideredCount);
#else	// not MEASURE
		fprintf(stderr, " with score %u\n", td->score);
#endif	// not MEASURE
	}

	if (td->numTrees < td->maxTrees) {
		// Add this tree to the list
		//TreeList_AddTree(&td->bestList, root, td);
		if (!td->tempTreeFile) {
			if (!(td->tempTreeFile = fopen(td->tempTreeFName, "w+b"))) {		// We need read access too, and non-Windows systems should ignore "b".
				perror(td->tempTreeFName);
				exit(1);
			}
		}

		// The trees we write out are already in Newick format; if we need NEXUS format, some rejigging
		// will be needed.  We don't bother writing out NEXUS format yet, since that requires overhead for
		// writing out the header, and #ifdef TARGETMULTI the final tree numbers will depend on the number
		// of trees produced by other workers anyway.
		PrintTree(root, td, td->tempTreeFile);
		fputc(';', td->tempTreeFile);
		fputc('\n', td->tempTreeFile);
	} else {
		td->optimalTreesMissed = 1;		// Only true if the current bound is really the lowest bound at the end.
	}

	++td->numTrees;
}

// For convenience.
#if LEFT + 1 == RIGHT
#define OTHERCHILD(node, x) (node)->nodes[LEFT + ((x) == (node)->nodes[LEFT])]
#define OTHERSIBLING(x) OTHERCHILD((x)->nodes[PARENT], (x))
#else
#error The definition of OTHERCHILD requires LEFT + 1 == RIGHT.
#endif

// This is the heart of the new "half-Yan-Bader" algorithm for updating the tree's sequences
// when a new taxon is added.  It should be called twice from the location of the newly
// inserted node.
// newNode is the node that has just been updated.  oldNode is the other node that we'll need.
// If we are entering from the parent, (newNode, oldNode) will be (PARENT, sibling) in some order;
// otherwise we are entering from a child, and (newNode, oldNode) will be (LEFT, RIGHT) in some order.
// Correctly handles the cases where node == NULL (parent was a leaf) or node == root.  Both just stop recursion.
// NOTE: We do not update any sequences of the root!  While this is fine for B&B because we never access any of
// these sequences, other code may need to propagate the UP sequence populated by ConstructBaseTree().  (That's
// the only sequence it populates.)
static void FASTCALL UpdateTree(struct tree *node, struct tree *newNode, struct tree *oldNode, int seqChanged, int nTaxaOnTree, struct TreeData *td) {
	// I have thought about the cases and it's correct to just stop recursing if we are the root, since the
	// only way we can be called is from our LEFT child, and all the seqs for the edge below us have just
	// been updated by the caller.
	if (node && node->nodes[PARENT]) {
		// Did we enter this node from the PARENT?
		int enterFromParent = (newNode == node->nodes[PARENT] || oldNode == node->nodes[PARENT]);
		if (enterFromParent) {
			// Entering from the parent.  We need to update our DOWN seq.
			// If the node visited prior to our parent was our parent's parent, we need to use our parent's new DOWN
			// and our sibling's old UP; otherwise, the node visited prior to our parent was our sibling, in which
			// case we need to use our parent's old DOWN and our sibling's new UP.
			// To accomplish this, just decide UP or DOWN based on whether the given node is the PARENT.
			if (seqChanged) {
				node->seqsByTreeSize[nTaxaOnTree].seqs[DOWN] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
				INCMEASURE(preorderFitchCallsFitchBasesCount);
				//CHOSEN_FitchBases(
				//	newNode->seqsByTreeSize[nTaxaOnTree].seqs[UP + (DOWN - UP) * (newNode == node->nodes[PARENT])],
				//	oldNode->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP + (DOWN - UP) * (oldNode == node->nodes[PARENT])],
				//	node->seqsByTreeSize[nTaxaOnTree].seqs[DOWN],
				//	td->seqLenBlocks
				//);
				seqChanged = CHOSEN_FitchBasesAndCompare(
					newNode->seqsByTreeSize[nTaxaOnTree].seqs[UP + (DOWN - UP) * (newNode == node->nodes[PARENT])],
					oldNode->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP + (DOWN - UP) * (oldNode == node->nodes[PARENT])],
					node->seqsByTreeSize[nTaxaOnTree].seqs[DOWN],
					node->seqsByTreeSize[nTaxaOnTree - 1].seqs[DOWN],
					td->seqLenBlocks
				);

				// Is the new DOWN sequence we just calculated identical to the previous one?
				//if (!memcmp(node->seqsByTreeSize[nTaxaOnTree].seqs[DOWN], node->seqsByTreeSize[nTaxaOnTree - 1].seqs[DOWN], td->seqLenBlocks * BYTESPERBLOCK)) {
				//if (!CompareSeqs_b2_w16(node->seqsByTreeSize[nTaxaOnTree].seqs[DOWN], node->seqsByTreeSize[nTaxaOnTree - 1].seqs[DOWN], td->seqLenBlocks)) {
				//	seqChanged = 0;
				//}
			}

			if (seqChanged) {		// Remember, it might have changed just now.
				// Our new UP is the same as before.
				node->seqsByTreeSize[nTaxaOnTree].seqs[UP] = node->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP];

				// To increase cache locality, update the MIDPOINT before recursing (even though this duplicates code...)
				node->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
				INCMEASURE(preorderFitchCallsFitchBasesCount);
				CHOSEN_FitchBases(
					node->seqsByTreeSize[nTaxaOnTree].seqs[UP],
					node->seqsByTreeSize[nTaxaOnTree].seqs[DOWN],
					node->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT],
					td->seqLenBlocks
				);
			} else {
				// The "new" sequence from the previous node visited hasn't changed, so neither will ours.
				// Just pointer-copy.
				INCMEASURE(updateTreePointerCopiesCount);
				node->seqsByTreeSize[nTaxaOnTree] = node->seqsByTreeSize[nTaxaOnTree - 1];	// Copy all 3 pointers
			}

			// Recurse on LEFT child.
			UpdateTree(node->nodes[LEFT], node, node->nodes[RIGHT], seqChanged, nTaxaOnTree, td);

			// Recurse on RIGHT child.
			UpdateTree(node->nodes[RIGHT], node, node->nodes[LEFT], seqChanged, nTaxaOnTree, td);
		} else {
			int origSeqChanged = seqChanged;		// Needed for pushing to child.

			// Entering from one of the children.  We need to update our UP seq.
			// We expect that newNode is the child we are entering from, and oldNode is the other child
			// (note the semantics of oldNode is different from in the case where we enter from PARENT).
			// Update our UP seq from the entering child's new UP and the non-entering child's old UP.
			if (seqChanged) {
				node->seqsByTreeSize[nTaxaOnTree].seqs[UP] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
				INCMEASURE(preorderFitchCallsFitchBasesCount);
				//CHOSEN_FitchBases(
				//	newNode->seqsByTreeSize[nTaxaOnTree].seqs[UP],
				//	oldNode->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP],
				//	node->seqsByTreeSize[nTaxaOnTree].seqs[UP],
				//	td->seqLenBlocks
				//);
				seqChanged = CHOSEN_FitchBasesAndCompare(
					newNode->seqsByTreeSize[nTaxaOnTree].seqs[UP],
					oldNode->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP],
					node->seqsByTreeSize[nTaxaOnTree].seqs[UP],
					node->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP],
					td->seqLenBlocks
				);

				// Is the new UP sequence we just calculated identical to the previous one?
				//if (!memcmp(node->seqsByTreeSize[nTaxaOnTree].seqs[UP], node->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP], td->seqLenBlocks * BYTESPERBLOCK)) {
				//if (!CompareSeqs_b2_w16(node->seqsByTreeSize[nTaxaOnTree].seqs[UP], node->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP], td->seqLenBlocks)) {
				//	seqChanged = 0;
				//}
			}

			if (seqChanged) {		// Remember, it might have changed just now.
				// Our new DOWN is the same as before.
				node->seqsByTreeSize[nTaxaOnTree].seqs[DOWN] = node->seqsByTreeSize[nTaxaOnTree - 1].seqs[DOWN];

				// To increase cache locality, update the MIDPOINT before recursing (even though this duplicates code...)
				node->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
				INCMEASURE(preorderFitchCallsFitchBasesCount);
				CHOSEN_FitchBases(
					node->seqsByTreeSize[nTaxaOnTree].seqs[UP],
					node->seqsByTreeSize[nTaxaOnTree].seqs[DOWN],
					node->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT],
					td->seqLenBlocks
				);
			} else {
				// The "new" sequence from the previous node visited hasn't changed, so neither will ours.
				// Just pointer-copy.
				INCMEASURE(updateTreePointerCopiesCount);
				node->seqsByTreeSize[nTaxaOnTree] = node->seqsByTreeSize[nTaxaOnTree - 1];	// Copy all 3 pointers
			}

			// Recurse on PARENT.  If PARENT is the root, OTHERSIBLING(node) will be NULL, but that
			// doesn't matter since UpdateTree() will bail out.
			UpdateTree(node->nodes[PARENT], node, OTHERSIBLING(node), seqChanged, nTaxaOnTree, td);

			// Recurse on other child.
			// Tricky: Although it feels like we should be the "new" node for the following call to UpdateTree(),
			// actually we are the "old" node, since that UpdateTree() instance will only look at our DOWN seq
			// (which hasn't changed), and needs the fresh new UP seq from the child we entered from.
			// Just as tricky: the child cares about whether the **other child's UP** has changed, not
			// about whether our UP has changed.
			UpdateTree(oldNode, newNode, node, origSeqChanged, nTaxaOnTree, td);
		}

		// For the following assert()s to be useful, freshly allocated memory must be set to 0xcdcdcdcd.
		assert(node->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] != 0xcdcdcdcd);
		assert(node->seqsByTreeSize[nTaxaOnTree].seqs[UP] != 0xcdcdcdcd);
		assert(node->seqsByTreeSize[nTaxaOnTree].seqs[DOWN] != 0xcdcdcdcd);
	}
}

//HACK: This relies on InsertNode() adding the new leaf on the left-hand side.
// It should be called immediately after InsertNode().
// Now actually needed for greedy tree construction in tbr.c as well as by branch & bound here.
void FASTCALL UpdateTreeAfterInsertion(struct tree *newInternal, int nTaxaOnTree, struct TreeData *td) {
	// We inherit our DOWN sequence from our right child.
	newInternal->seqsByTreeSize[nTaxaOnTree].seqs[DOWN] = newInternal->nodes[RIGHT]->seqsByTreeSize[nTaxaOnTree - 1].seqs[DOWN];

	// Our right child will need to know our "old" DOWN sequence (which is the same)
	newInternal->seqsByTreeSize[nTaxaOnTree - 1].seqs[DOWN] = newInternal->seqsByTreeSize[nTaxaOnTree].seqs[DOWN];

	// Compute our own UP sequence
	newInternal->seqsByTreeSize[nTaxaOnTree].seqs[UP] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
	INCMEASURE(preorderFitchCallsFitchBasesCount);
	CHOSEN_FitchBases(
		newInternal->nodes[LEFT]->seqsByTreeSize[nTaxaOnTree].seqs[UP],			// The new sequence from the new node on the left
		newInternal->nodes[RIGHT]->seqsByTreeSize[nTaxaOnTree - 1].seqs[UP],	// The old sequence from the right
		newInternal->seqsByTreeSize[nTaxaOnTree].seqs[UP],
		td->seqLenBlocks
	);

	// Compute our own MIDPOINT sequence
	newInternal->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
	INCMEASURE(preorderFitchCallsFitchBasesCount);
	CHOSEN_FitchBases(
		newInternal->seqsByTreeSize[nTaxaOnTree].seqs[UP],
		newInternal->seqsByTreeSize[nTaxaOnTree].seqs[DOWN],
		newInternal->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT],
		td->seqLenBlocks
	);

	// The new leaf's DOWN sequence will be our right child's old MIDPOINT sequence.  (Think about it...)
	newInternal->nodes[LEFT]->seqsByTreeSize[nTaxaOnTree].seqs[DOWN] = newInternal->nodes[RIGHT]->seqsByTreeSize[nTaxaOnTree - 1].seqs[MIDPOINT];

	// Compute the new leaf's MIDPOINT sequence
	newInternal->nodes[LEFT]->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT] = ArenaAllocate(&td->alignedMem, td->seqLenBlocks * BYTESPERBLOCK);
	INCMEASURE(preorderFitchCallsFitchBasesCount);
	CHOSEN_FitchBases(
		newInternal->nodes[LEFT]->seqsByTreeSize[nTaxaOnTree].seqs[UP],
		newInternal->nodes[LEFT]->seqsByTreeSize[nTaxaOnTree].seqs[DOWN],
		newInternal->nodes[LEFT]->seqsByTreeSize[nTaxaOnTree].seqs[MIDPOINT],
		td->seqLenBlocks
	);

	// Push changes below and to the right of the new internal node (there may be a large subtree here)
	UpdateTree(newInternal->nodes[RIGHT], newInternal->nodes[LEFT], newInternal, 1, nTaxaOnTree, td);

	// Push changes above the new internal node.  UpdateTree() handles the node == root case correctly.
	UpdateTree(newInternal->nodes[PARENT], newInternal, OTHERSIBLING(newInternal), 1, nTaxaOnTree, td);
}

#undef OTHERCHILD
#undef OTHERSIBLING
