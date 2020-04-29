#include <string.h>
#include <limits.h>		// For INT_MAX
#include "common.h"
#include <dbgprint.h>
#include "ordertaxa.h"
#include "seq.h"

static int CompareUnsigned(const void *a, const void *b);

#ifdef MAXMINITREEORDER
#error 21/2/2005: The MAXMINITREEORDER version of OrderTaxa() needs to be updated to use the 1-base-per-byte functions, and the ...ScanEdge...() functions need to stop using InsertNode(), which will produce a node of the wrong size.
#error 2/2/2010: Also, we need to figure out what the calls to CopyTree() should actually do.
void OrderTaxa(struct TreeData *td) {
	unsigned char *seqbuf = (unsigned char *) malloc(td->seqLen);
	unsigned i, j, n, oldBound;
	unsigned furthesti = 0, furthestj = 1, furthestdist = 0;
	unsigned totaldist, dummy, temp;
	struct tree *root, *furthesttree;
//	struct TreeList *oldBestList;
	TreeList oldBestList;
	unsigned nBestEqual;
	unsigned *distBuf, *furthestDistBuf;
	unsigned distBufPos;
	
//	DBGPRINT1("Before OrderTaxa():\n");
//	for (i = 0; i < td->numTaxa; ++i)
//	{
//		PrintMaskSeq(td->charMat + i * td->seqLen, td->seqLen);
//		DBGPRINT1("\n");
//	}
	
	td->taxonMap[td->numTaxa - 1] = td->numTaxa - 1;		// This is true initially
	nBestEqual = 0;
	for (i = 0; i < td->numTaxa - 1; ++i) {
		td->taxonMap[i] = i;						// This is true initially
		
		for (j = i + 1; j < td->numTaxa; ++j) {
			totaldist = Fitch2SeqsCost(td->charMat + i * td->seqLen, td->charMat + j * td->seqLen, td->weights, td->seqLen);
			if (totaldist > furthestdist) {
				furthesti = i;
				furthestj = j;
				furthestdist = totaldist;
				nBestEqual = 1;
			} else if (totaldist == furthestdist) {
				++nBestEqual;
			}
		}
	}
	
	DBGPRINT4("Most-different pair are %d and %d with a distance of %d.\n", furthesti, furthestj, furthestdist);
	DBGPRINT3("There were %u best-equal pairs of taxa (each having weight %u).\n", nBestEqual, furthestdist);
	
	// Swap most-different pair with taxa at front
	if (furthesti != 0) {
		memcpy(seqbuf, td->charMat, td->seqLen);
		memcpy(td->charMat, td->charMat + furthesti * td->seqLen, td->seqLen);
		memcpy(td->charMat + furthesti * td->seqLen, seqbuf, td->seqLen);
		td->taxonMap[0] = furthesti;
		td->taxonMap[furthesti] = 0;
	}
	if (furthestj != 1) {
		memcpy(seqbuf, td->charMat + td->seqLen, td->seqLen);
		memcpy(td->charMat + td->seqLen, td->charMat + furthestj * td->seqLen, td->seqLen);
		memcpy(td->charMat + furthestj * td->seqLen, seqbuf, td->seqLen);
		temp = td->taxonMap[furthestj];
		td->taxonMap[furthestj] = td->taxonMap[1];
		td->taxonMap[1] = temp;
	}
	
	// Construct a rooted tree with these two taxa as children
	root = NewNode(1, &gTD);
	AddToPendantNode(root, LEFT, NewNode(-1, &gTD));
	AddToPendantNode(root->nodes[LEFT], LEFT, NewNode(2, &gTD));
	AddToPendantNode(root->nodes[LEFT], RIGHT, NewNode(3, &gTD));
#ifdef SMARTFITCH
	root->nodes[LEFT]->cumScore = Fitch2Seqs(root->nodes[LEFT]->nodes[LEFT]->sequence, root->nodes[LEFT]->nodes[RIGHT]->sequence, root->nodes[LEFT]->sequence, td->weights, td->seqLen);
	root->cumScore = Fitch2SeqsCost(root->nodes[LEFT]->sequence, root->sequence, td->weights, td->seqLen) + root->nodes[LEFT]->cumScore;
DBGPRINT2("root->nodes[LEFT]->cumScore = <%d>\n", root->nodes[LEFT]->cumScore);		//DEBUG
DBGPRINT2("root->cumScore = <%d>\n", root->cumScore);		//DEBUG
#endif	// SMARTFITCH

	// Now to find the third taxon, form a 3-taxon tree using each candidate
	// taxon in turn, and pick the one that gives the heaviest tree.
	furthesti = 2;
	furthestdist = 0;
	nBestEqual = 0;
	for (n = 2; n < td->numTaxa; ++n) {
		root->nodes[LEFT]->nodes[RIGHT]->contents = n + 1;
		root->nodes[LEFT]->nodes[RIGHT]->sequence = td->charMat + n * td->seqLen;
		totaldist = Fitch(root, td);
		
		if (totaldist > furthestdist) {
			furthesti = n;
			furthestdist = totaldist;
			nBestEqual = 1;
		} else if (totaldist == furthestdist) {
			++nBestEqual;
		}
	}
	
	if (furthesti != 2) {
		memcpy(seqbuf, td->charMat + td->seqLen * 2, td->seqLen);
		memcpy(td->charMat + td->seqLen * 2, td->charMat + furthesti * td->seqLen, td->seqLen);
		memcpy(td->charMat + furthesti * td->seqLen, seqbuf, td->seqLen);
		temp = td->taxonMap[furthesti];
		td->taxonMap[furthesti] = td->taxonMap[2];
		td->taxonMap[2] = temp;
	}
	
	root->nodes[LEFT]->nodes[RIGHT]->contents = 3;
	root->nodes[LEFT]->nodes[RIGHT]->sequence = td->charMat + td->seqLen * 2;
	
#ifdef DEBUG
	DBGPRINT1("The maxmini tree on 3 taxa is: ");
	PrintTree(root, td, dbgfile);
	DBGPRINT1("\n");
	DBGPRINT3("There were %u best-equal taxa to choose from (each having weight %u).\n", nBestEqual, furthestdist);
#endif	// DEBUG
	
	// Now to find the nth taxon, sum the distances from each candidate taxon
	// to all the n-1 taxa already decided upon.  The taxon with the largest
	// sum becomes the nth taxon.
	oldBestList = td->bestList;			// Will restore later
	oldBound = td->bound;				// Will restore later
//	td->bestList = (struct TreeList *) malloc(sizeof (struct TreeList));
//	td->bestList->t = NULL;
	TreeList_Destroy(&td->bestList);
	TreeList_Init(&td->bestList);
	distBuf = (unsigned *) malloc(td->numTaxa * 2 * sizeof (unsigned));
	furthestDistBuf = (unsigned *) malloc(td->numTaxa * 2 * sizeof (unsigned));
	for (n = 3; n < td->numTaxa - 1; ++n) {
		furthesti = 3;
//		furthestdist = 0;
		memset(furthestDistBuf, 0, td->numTaxa * 2 * sizeof (unsigned));
		nBestEqual = 0;
		for (i = n; i < td->numTaxa; ++i) {
			// We will misuse the "bound" field for storing the lowest score obtained
			td->bound = INT_MAX;
//			OrderTaxaScanEdge(root, root, LEFT, i, td);
			distBufPos = 0;
			OrderTaxaScanEdgeThorough(root, root, LEFT, i, td, distBuf, &distBufPos);
			qsort(distBuf, distBufPos, sizeof (unsigned), CompareUnsigned);
//			DBGPRINT2("VerifyTree(td->bestList->t) = %d\n", VerifyTree(td->bestList->t, NULL));
			
//			if (td->bound > furthestdist) {
//				furthesti = i;
//				furthestdist = td->bound;
////				furthesttree = td->bestList->t;
//				furthesttree = CopyTree(TreeList_GetTree(td->bestList, TreeList_GetIterator(td->bestList)), td);
//				nBestEqual = 1;
//			} else if (td->bound == furthestdist) {
//				++nBestEqual;
//			}
			for (j = 0; j < distBufPos; ++j) {
				if (distBuf[j] > furthestDistBuf[j]) {
					// We have a new maximum
					furthesti = i;
					memcpy(furthestDistBuf, distBuf, distBufPos * sizeof (unsigned));
					furthesttree = CopyTree(TreeList_GetTree(td->bestList, TreeList_GetIterator(td->bestList)), td, 1);
					nBestEqual = 1;
					break;
				} else if (distBuf[j] < furthestDistBuf[j]) {
					// This is worse than the current best; ignore it
					break;
				}
			}
			
			if (j == distBufPos) {
				++nBestEqual;			// Shouldn't happen too often!
			}
			
//			td->bestList->t = NULL;
			TreeList_Destroy(&td->bestList);
			TreeList_Init(&td->bestList);
		}
		
#ifdef DEBUG
//		DBGPRINT4("The furthest taxon is #%d with distance %d.  The maxmini tree on %d taxa is: ", furthesti + 1, furthestdist, n + 1);
		fprintf(dbgfile, "The furthest taxon is #%d with distances %d,%d,%d,%d,%d.  The maxmini tree on %d taxa is: ", furthesti + 1, furthestDistBuf[0], furthestDistBuf[1], furthestDistBuf[2], furthestDistBuf[3], furthestDistBuf[4], n + 1);
		PrintTree(furthesttree, td, dbgfile);
		DBGPRINT1("\n");
		DBGPRINT2("There were %u best-equal taxa to choose from.\n", nBestEqual);
#endif	// DEBUG
		
		// Swap most-distant taxon with nth taxon
//		DBGPRINT3("Most-different taxon is %d with a total distance of %d.\n", furthesti, furthestdist);
		if (furthesti != n) {
			memcpy(seqbuf, td->charMat + n * td->seqLen, td->seqLen);
			memcpy(td->charMat + n * td->seqLen, td->charMat + furthesti * td->seqLen, td->seqLen);
			memcpy(td->charMat + furthesti * td->seqLen, seqbuf, td->seqLen);
			temp = td->taxonMap[furthesti];
			td->taxonMap[furthesti] = td->taxonMap[n];
			td->taxonMap[n] = temp;
			RenameTaxon(furthesttree, furthesti + 1, n + 1);	// After reordering taxa, this taxon is always taxon #n
		}
		
		DestroyTree(root);
		root = furthesttree;
	}
	
	free(distBuf);
	free(furthestDistBuf);
	
//	DBGPRINT1("After OrderTaxa():\n");
//	for (i = 0; i < td->numTaxa; ++i)
//	{
//		PrintMaskSeq(td->charMat + i * td->seqLen, td->seqLen);
//		DBGPRINT1("\n");
//	}
	
//	free(td->bestList);
	TreeList_Destroy(&td->bestList);
	free(seqbuf);
	
	td->bestList = oldBestList;
	td->bound = oldBound;
}

// Change the first occurrence of taxon #from to #to on a tree.  Needed by MAXMINITREEORDER.
void RenameTaxon(struct tree *root, int from, int to) {
	if (!root) return;
	
	if (root->contents == from) {
//		DBGPRINT3("RenameTaxon(): Renamed taxon #%d to #%d.\n", from, to);
		root->contents = to;
		return;
	}
	
	RenameTaxon(root->nodes[LEFT], from, to);
	RenameTaxon(root->nodes[RIGHT], from, to);
}

void OrderTaxaScanEdge(struct tree *root, struct tree *parent, int iChild, unsigned numTaxa, struct TreeData *td) {
	struct tree *newParent;
	struct tree *child = parent->nodes[iChild];		//HACK
	unsigned score;
#ifdef SMARTFITCH
	struct tree *p, *tempTree;
#endif	// SMARTFITCH
	
	// Is there an edge?
	if (!child) return;
	
	// Recurse
	OrderTaxaScanEdge(root, child, LEFT, numTaxa, td);
	OrderTaxaScanEdge(root, child, RIGHT, numTaxa, td);
	
	// Now we do the actual work
	newParent = InsertNode(parent, iChild, numTaxa+1, td);		//HACK: not using newParent yet, but I will...
	
#ifdef SMARTFITCH
	//WTJW 31/3/2003: Recompute Fitch only along path from new node to root.
	for (p = newParent; p != root; p = p->nodes[PARENT]) {
		p->cumScore = Fitch2Seqs(p->nodes[LEFT]->sequence, p->nodes[RIGHT]->sequence, p->sequence, td->weights, td->seqLen) + p->nodes[LEFT]->cumScore + p->nodes[RIGHT]->cumScore;
	}
	root->cumScore = Fitch2SeqsCost(root->nodes[LEFT]->sequence, root->sequence, td->weights, td->seqLen) + root->nodes[LEFT]->cumScore;
	score = root->cumScore;
#else	// not SMARTFITCH
	score = Fitch(root, td);
#endif	// not SMARTFITCH
	score += td->constantWeight;	// Non-PI sites contribute this much weight to any tree
	if (score < td->bound) {
		td->bound = score;
//		if (td->bestList->t)
//			DestroyTree(td->bestList->t);			//WTJW: Remove previous best tree
//		td->bestList->t = CopyTree(root, td);
		TreeList_Destroy(&td->bestList);			//WTJW: Remove previous best tree
		TreeList_AddTree(&td->bestList, root, td);
	}	
	
	DeleteNode(parent, iChild);
	
	//WTJW 31/3/2003: Recompute Fitch only along path from new node to root.
#ifdef SMARTFITCH
	for (p = parent; p != root; p = p->nodes[PARENT]) {
		p->cumScore = Fitch2Seqs(p->nodes[LEFT]->sequence, p->nodes[RIGHT]->sequence, p->sequence, td->weights, td->seqLen) + p->nodes[LEFT]->cumScore + p->nodes[RIGHT]->cumScore;
	}
	root->cumScore = Fitch2SeqsCost(root->nodes[LEFT]->sequence, root->sequence, td->weights, td->seqLen) + root->nodes[LEFT]->cumScore;
#endif	// SMARTFITCH
}

void OrderTaxaScanEdgeThorough(struct tree *root, struct tree *parent, int iChild, unsigned numTaxa, struct TreeData *td, unsigned *distBuf, unsigned *distBufPosPtr) {
	struct tree *newParent;
	struct tree *child = parent->nodes[iChild];		//HACK
	unsigned score;
#ifdef SMARTFITCH
	struct tree *p, *tempTree;
#endif	// SMARTFITCH
	
	// Is there an edge?
	if (!child) return;
	
//	DBGPRINT4("OrderTaxaScanEdgeThorough(numTaxa=%u, distBuf=%u, *distBufPosPtr=%u) called.\n", numTaxa, distBuf, *distBufPosPtr);
	
	// Recurse
	OrderTaxaScanEdgeThorough(root, child, LEFT, numTaxa, td, distBuf, distBufPosPtr);
	OrderTaxaScanEdgeThorough(root, child, RIGHT, numTaxa, td, distBuf, distBufPosPtr);
	
	// Now we do the actual work
	newParent = InsertNode(parent, iChild, numTaxa+1, td);		//HACK: not using newParent yet, but I will...
	
#ifdef SMARTFITCH
	//WTJW 31/3/2003: Recompute Fitch only along path from new node to root.
	for (p = newParent; p != root; p = p->nodes[PARENT]) {
		p->cumScore = Fitch2Seqs(p->nodes[LEFT]->sequence, p->nodes[RIGHT]->sequence, p->sequence, td->weights, td->seqLen) + p->nodes[LEFT]->cumScore + p->nodes[RIGHT]->cumScore;
	}
	root->cumScore = Fitch2SeqsCost(root->nodes[LEFT]->sequence, root->sequence, td->weights, td->seqLen) + root->nodes[LEFT]->cumScore;
	score = root->cumScore;
#else	// not SMARTFITCH
	score = Fitch(root, td);
#endif	// not SMARTFITCH
	score += td->constantWeight;	// Non-PI sites contribute this much weight to any tree
	// The following is still necessary, as we still need to keep track of the best tree for this taxon.
	if (score < td->bound) {
		td->bound = score;
//		if (td->bestList->t)
//			DestroyTree(td->bestList->t);			//WTJW: Remove previous best tree
//		td->bestList->t = CopyTree(root, td);
		TreeList_Destroy(&td->bestList);			//WTJW: Remove previous best tree
		TreeList_AddTree(&td->bestList, root, td);
	}	
	
	distBuf[*distBufPosPtr] = score;
	++(*distBufPosPtr);
	
	DeleteNode(parent, iChild);
	
	//WTJW 31/3/2003: Recompute Fitch only along path from new node to root.
#ifdef SMARTFITCH
	for (p = parent; p != root; p = p->nodes[PARENT]) {
		p->cumScore = Fitch2Seqs(p->nodes[LEFT]->sequence, p->nodes[RIGHT]->sequence, p->sequence, td->weights, td->seqLen) + p->nodes[LEFT]->cumScore + p->nodes[RIGHT]->cumScore;
	}
	root->cumScore = Fitch2SeqsCost(root->nodes[LEFT]->sequence, root->sequence, td->weights, td->seqLen) + root->nodes[LEFT]->cumScore;
#endif	// SMARTFITCH
}
#else	// not MAXMINITREEORDER
void OrderTaxa(struct TreeData *td) {
	unsigned char *seqbuf = (unsigned char *) malloc(td->seqLen);
	unsigned i, j, n, temp;
	unsigned furthesti = 0, furthestj = 1, furthestdist = 0;
	unsigned totaldist;
	unsigned nBestEqual;
#ifdef MAXMINIORDER
	unsigned *distBuf, *furthestDistBuf;
#endif	// MAXMINIORDER
	
	td->taxonMap = (unsigned *) malloc(td->numTaxa * sizeof (unsigned));		// 28/2/2005: Moved allocation into this function
	td->taxonMap[td->numTaxa - 1] = td->numTaxa - 1;		// This is true initially
	nBestEqual = 0;
	for (i = 0; i < td->numTaxa - 1; ++i) {
		td->taxonMap[i] = i;						// This is true initially
		
		for (j = i + 1; j < td->numTaxa; ++j) {
			totaldist = FitchScore_b1_w1_basic(td->charMat + i * td->seqLen, td->charMat + j * td->seqLen, td->weights, td->seqLen);
			if (totaldist > furthestdist) {
				furthesti = i;
				furthestj = j;
				furthestdist = totaldist;
				nBestEqual = 1;
			} else if (totaldist == furthestdist) {
				++nBestEqual;
			}
		}
	}
	
	DBGPRINT4("Most-different pair are %d and %d with a distance of %d.\n", furthesti + 1, furthestj + 1, furthestdist);
	DBGPRINT3("There were %u best-equal pairs of taxa (each having weight %u).\n", nBestEqual, furthestdist);
	
	// Swap most-different pair with taxa at front
	if (furthesti != 0) {
		memcpy(seqbuf, td->charMat, td->seqLen);
		memcpy(td->charMat, td->charMat + furthesti * td->seqLen, td->seqLen);
		memcpy(td->charMat + furthesti * td->seqLen, seqbuf, td->seqLen);
		td->taxonMap[0] = furthesti;		// Don't need a temporary for the very first swap.
		td->taxonMap[furthesti] = 0;
	}
	if (furthestj != 1) {
		memcpy(seqbuf, td->charMat + td->seqLen, td->seqLen);
		memcpy(td->charMat + td->seqLen, td->charMat + furthestj * td->seqLen, td->seqLen);
		memcpy(td->charMat + furthestj * td->seqLen, seqbuf, td->seqLen);
		temp = td->taxonMap[furthestj];
		td->taxonMap[furthestj] = td->taxonMap[1];
		td->taxonMap[1] = temp;
	}
	
	// When MAXMINIORDER is #defined (which is recommended), we find the nth taxon by measuring distances from each
	// remaining taxon to all the n-1 taxa already decided upon and taking the taxon whose
	// nearest distance is greatest.  If two taxa have equal nearest distances, we look at
	// the next nearest distance, and so on -- this means we need to keep track of all distances.
	// If MAXMINIORDER is not #defined, to find the nth taxon, we sum the distances from each candidate taxon
	// to all the n-1 taxa already decided upon.  The taxon with the largest
	// sum becomes the nth taxon.
#ifdef MAXMINIORDER
	// We won't actually use totaldist or furthestdist beyond this point.
	distBuf = (unsigned *) malloc(td->numTaxa * sizeof (unsigned));		// We don't actually record the other taxon here,
	furthestDistBuf = (unsigned *) malloc(td->numTaxa * sizeof (unsigned));	// although we could for better diagnostics.
#endif	// MAXMINIORDER
	for (n = 2; n < td->numTaxa - 1; ++n) {
#ifdef MAXMINIORDER
		memset(furthestDistBuf, 0, td->numTaxa * sizeof (unsigned));
#else	// not MAXMINIORDER
		furthestdist = 0;			// Not used by MAXMINIORDER; it uses furthestDistBuf instead.
#endif	// not MAXMINIORDER
		furthesti = n;
		nBestEqual = 0;
		for (i = n; i < td->numTaxa; ++i) {
			// Let's consider making taxon i the nth taxon.
#ifndef MAXMINIORDER
			totaldist = 0;
#endif	// not MAXMINIORDER
			for (j = 0; j < n; ++j) {
#ifdef MAXMINIORDER
				distBuf[j] =
#else	// not MAXMINIORDER
				totaldist +=
#endif	// not MAXMINIORDER
					FitchScore_b1_w1_basic(td->charMat + i * td->seqLen, td->charMat + j * td->seqLen, td->weights, td->seqLen);
			}

#ifdef MAXMINIORDER
			// Look at all distances starting from the nearest, until we find a difference.
			qsort(distBuf, n, sizeof (unsigned), CompareUnsigned);
			for (j = 0; j < n; ++j) {
				if (distBuf[j] > furthestDistBuf[j]) {
					// We have a new maximum
					furthesti = i;
					memcpy(furthestDistBuf, distBuf, n * sizeof (unsigned));
					nBestEqual = 1;
					break;
				} else if (distBuf[j] < furthestDistBuf[j]) {
					// This is worse than the current best; ignore it
					break;
				}
			}
			
			if (j == n) {
				++nBestEqual;			// Shouldn't happen too often!
			}
#else	// not MAXMINIORDER
			// Just compare the distance sums.
			if (totaldist > furthestdist) {
				// We have a new maximum
				furthesti = i;
				furthestdist = totaldist;
				nBestEqual = 1;
			} else if (totaldist == furthestdist) {
				++nBestEqual;
			}
#endif	// not MAXMINIORDER
		}
		
		// By this time, furthesti identifies the most-distant taxon.
		// Swap most-distant taxon with nth taxon
#ifdef MAXMINIORDER
		// We don't actually record which other taxon furthesti is furthest from, so don't try to report it!
		DBGPRINT3("Most-different taxon is %u with a nearest distance of %u from a taxon already considered.\n", td->taxonMap[furthesti] + 1, furthestDistBuf[0]);
		DBGPRINT3("There were %u best-equal pairs of taxa (each having nearest distance %u).\n", nBestEqual, furthestDistBuf[0]);
#else	// not MAXMINIORDER
		DBGPRINT4("Most-different taxon is %u with a total distance of %u from taxon %u.\n", td->taxonMap[furthesti] + 1, furthestdist, td->taxonMap[furthestj] + 1);
		DBGPRINT3("There were %u best-equal pairs of taxa (each having weight %u).\n", nBestEqual, furthestdist);
#endif	// not MAXMINIORDER
		if (furthesti != n) {
			memcpy(seqbuf, td->charMat + n * td->seqLen, td->seqLen);
			memcpy(td->charMat + n * td->seqLen, td->charMat + furthesti * td->seqLen, td->seqLen);
			memcpy(td->charMat + furthesti * td->seqLen, seqbuf, td->seqLen);
			temp = td->taxonMap[furthesti];
			td->taxonMap[furthesti] = td->taxonMap[n];
			td->taxonMap[n] = temp;
		}
	}
	
#ifdef MAXMINIORDER
	free(furthestDistBuf);
	free(distBuf);
#endif	// MAXMINIORDER
	free(seqbuf);
}
#endif	// not MAXMINITREEORDER

int CompareUnsigned(const void *a, const void *b) {
	return *((unsigned *) a) - *((unsigned *) b);
//	return *((unsigned *) b) - *((unsigned *) a);
}
