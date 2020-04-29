#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
//#include <mem.h>
#include <string.h>
#include <assert.h>
#include "common.h"
#include <dbgprint.h>
#include "seq.h"
#include "measure.h"

// Not too efficient but makes things easier...
unsigned char FASTCALL GetBaseAt_b2_w1(unsigned taxon, unsigned site, struct TreeData *td) {
	if (site & 1) {
		return td->charMat[taxon * UDIVROUNDUP(td->seqLen, 2) + site / 2] >> 4;
	} else {
		return td->charMat[taxon * UDIVROUNDUP(td->seqLen, 2) + site / 2] & 0x0F;
	}
}

void FASTCALL PrintMaskSeq_b2_w1(unsigned char *seq, unsigned memSeqLen, FILE *f) {
	unsigned i;
	unsigned char b;
	
	for (i = 0; i < memSeqLen * 2; ++i) {
		b = (i & 1) ? (seq[i / 2] >> 4) : (seq[i / 2] & 0x0F);
		switch (b) {
		case BM_A:
		case BM_C:
		case BM_G:
		case BM_T:
			fputc(GetBaseFromMask(b), f);
			break;
		
		default:
#ifdef SHOWSHORTBASESETS
			fputc(GetBaseFromMask(b), f);
#else	// not SHOWSHORTBASESETS
			fputc('[', f);
			if (b & BM_A) fputc('A', f);
			if (b & BM_C) fputc('C', f);
			if (b & BM_G) fputc('G', f);
			if (b & BM_T) fputc('T', f);
			fputc(']', f);
#endif	// SHOWSHORTBASESETS
		}
	}
}

// Form a complete graph, with sequences as vertices and edges weighted by
// pairwise Hamming distances.  Compute a minimum bound on the parsimony weight
// by taking half the weight of a MST on this graph.
unsigned FASTCALL CalcMstWeight_b2_w1(unsigned char *data, unsigned dataWidth, unsigned width, unsigned height, unsigned *weights) {
	struct tree root, left, right;		// Actual structures (not pointers)
	unsigned *hDist, *onTree;
	unsigned dist;
	unsigned minDist, minI, minJ, mstWeight = 0;
	unsigned i, j, k;
	
	hDist  = (unsigned *) malloc(height * height * sizeof (unsigned));
	onTree = (unsigned *) malloc(height * sizeof (unsigned));
	memset(onTree, 0, height * sizeof (unsigned));	//HACK: lazy
	memset(hDist, 0, height * height * sizeof (unsigned));	//DEBUG
	
	// Use Fitch() to compute Hamming distances between pairs of sequences
	for (i = 0; i < height - 1; ++i) {
		for (j = i + 1; j < height; ++j) {
			// Find Hamming distance between i-th and j-th taxa
//			dist = Fitch2Seqs_BASIC(td->charMat + i * td->memSeqLen, td->charMat + j * td->memSeqLen, NULL, td->weights, td->memSeqLen);
//			dist = Fitch2SeqsCost_b1_w1_basic(td->charMat + i * td->memSeqLen, td->charMat + j * td->memSeqLen, td->weights, td->memSeqLen);
			dist = FitchScore_b2_w1(data + i * width, data + j * width, weights, dataWidth);
			
			//HACK: lazy, don't need to put it in both places
			hDist[i * height + j] = hDist[j * height + i] = dist;
		}
	}
	
	// Find the first edge.
	minDist = INT_MAX;
	for (i = 0; i < height - 1; ++i) {
		for (j = i + 1; j < height; ++j) {
			if (hDist[i * height + j] < minDist) {
				minDist = hDist[i * height + j];
				minI = i;
				minJ = j;
			}
		}
	}
	
	mstWeight += minDist;
	onTree[minI] = 1;
	onTree[minJ] = 1;
	
	DBGPRINT4("First edge is (%u,%u) with weight %u.\n", minI, minJ, minDist);
	
	// Find all remaining edges.  One end must be on the tree, the other off.
	for (k = 0; k < height - 2; ++k) {		// One edge is already present
		// Find the next edge.
		minDist = INT_MAX;
		for (i = 0; i < height; ++i) {
			if (onTree[i]) {
				for (j = 0; j < height; ++j) {
					if (!onTree[j]) {
						// The edge from i to j is a candidate
						if (hDist[i * height + j] < minDist) {
							minDist = hDist[i * height + j];
							minI = i;
							minJ = j;
						}
					}
				}
			}
		}
		
		mstWeight += minDist;
		onTree[minI] = 1;
		onTree[minJ] = 1;
		DBGPRINT4("Added edge (%u,%u) with weight %u.\n", minI, minJ, minDist);
	}
	
	//DEBUG
	for (i = 0; i < height; ++i) {
		assert(onTree[i]);
	}
	DBGPRINT1("Hamming distances:\n   ");
	for (i = 0; i < height; ++i) {
		DBGPRINT2(" %2d", i);
	}
	DBGPRINT1("\n");
	for (i = 0; i < height; ++i) {
		DBGPRINT2("%2d:", i);
		for (j = 0; j < height; ++j) {
			DBGPRINT2(" %2d", hDist[i * height + j]);
		}
		DBGPRINT1("\n");
	}
	
	free(hDist);
	free(onTree);
	
	return mstWeight;
}

unsigned FASTCALL CalcTreeMstWeight_b2_w1(struct TreeData *td) {
	return CalcMstWeight_b2(td->charMat, UDIVROUNDUP(td->seqLen, 2), UDIVROUNDUP(td->seqLen, 2), td->numTaxa, td->weights);
}


///////////////////////////////////////////////////////////////////////////////
// Functions for all versions of BASIC Fitch computations.
///////////////////////////////////////////////////////////////////////////////

unsigned FASTCALL FitchBoth_b2_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned *weights, unsigned nBlocks) {
	unsigned thisScore = 0, i;
	unsigned char newVal;
	unsigned char newVal2;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < nBlocks; ++i) {
		if (!(newVal = s1[i] & s2[i] & 0x0F)) {
			// No intersection of base sets: form union
			newVal = (s1[i] | s2[i]) & 0x0F;
			thisScore += weights[i * 2];
		}
		if (!(newVal2 = s1[i] & s2[i] & 0xF0)) {
			// No intersection of base sets: form union
			newVal2 = (s1[i] | s2[i]) & 0xF0;
			thisScore += weights[i * 2 + 1];
		}
		newVal |= newVal2;
		
		dest[i] = newVal;
	}
#endif	// EMPTYFITCH
	return thisScore;
}

void FASTCALL FitchBases_b2_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned nBlocks) {
	unsigned i;
	unsigned char newVal;
	unsigned char newVal2;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < nBlocks; ++i) {
		if (!(newVal = s1[i] & s2[i] & 0x0F)) {
			// No intersection of base sets: form union
			newVal = (s1[i] | s2[i]) & 0x0F;
		}
		if (!(newVal2 = s1[i] & s2[i] & 0xF0)) {
			// No intersection of base sets: form union
			newVal2 = (s1[i] | s2[i]) & 0xF0;
		}
		newVal |= newVal2;
		
		dest[i] = newVal;
	}
#endif	// EMPTYFITCH
}

unsigned FASTCALL FitchScoreWeight1_b2_w1_basic(unsigned char *s1, unsigned char *s2, unsigned nBlocks) {
	unsigned thisScore = 0, i;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < nBlocks; ++i) {
		if (!(s1[i] & s2[i] & 0x0F)) {
			// No intersection of base sets: form union
			++thisScore;
		}
		if (!(s1[i] & s2[i] & 0xF0)) {
			// No intersection of base sets: form union
			++thisScore;
		}
	}
#endif	// EMPTYFITCH
	return thisScore;
}

unsigned FASTCALL FitchScore_b2_w1_basic(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned nBlocks) {
	unsigned thisScore = 0, i;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < nBlocks; ++i) {
		if (!(s1[i] & s2[i] & 0x0F)) {
			// No intersection of base sets: form union
			thisScore += weights[i * 2];
		}
		if (!(s1[i] & s2[i] & 0xF0)) {
			// No intersection of base sets: form union
			thisScore += weights[i * 2 + 1];
		}
	}
#endif	// EMPTYFITCH
	return thisScore;
}
