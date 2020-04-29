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
unsigned char FASTCALL GetBaseAt_b1_w1(unsigned taxon, unsigned site, struct TreeData *td) {
	return td->charMat[taxon * td->seqLen + site];
}

void FASTCALL PrintMaskSeq_b1_w1(unsigned char *seq, unsigned memSeqLen, FILE *f) {
	unsigned i;
	unsigned char b;
	
	for (i = 0; i < memSeqLen; ++i) {
		b = seq[i];
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
// WTJW 23/3/2005: If underestimate is 1, the standard Fitch function will be
// used to calculate Hamming distances, meaning that for example the path
// A-N-T will be zero-cost.  This is necessary when computing lower bounds.
// If underestimate is 0, a modified Hamming distance function is used that
// incurs a cost of 1 any time two base sets are not exactly equal, so that for
// example the path A-N will be have cost 1.  This is necessary when computing
// an initial upper bound.
unsigned FASTCALL CalcMstWeight_b1_w1(unsigned char *data, unsigned dataWidth, unsigned width, unsigned height, unsigned *weights, int underestimate) {
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
//			dist = Fitch2SeqsCost_b1(data + i * width, data + j * width, weights, dataWidth);
			if (underestimate) {
				dist = FitchScore_b1_w1(data + i * width, data + j * width, weights, dataWidth);
			} else {
				dist = StrictlyEqualFitch2SeqsCost_b1(data + i * width, data + j * width, weights, dataWidth);
			}
			
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

unsigned FASTCALL CalcTreeMstWeight_b1_w1(struct TreeData *td, int underestimate) {
//	return CalcMstWeight_b1(td->charMat, td->memSeqLen, td->memSeqLen, td->numTaxa, td->weights, underestimate);
	DBGPRINT6("Calling CalcMstWeight_b1(%d, %d, %d, %d, %d)\n", td->charMat, td->seqLen, td->seqLen, td->numTaxa, td->weights);
	DBGPRINT2("underestimate = %d\n", underestimate);
	return CalcMstWeight_b1(td->charMat, td->seqLen, td->seqLen, td->numTaxa, td->weights, underestimate);
}


unsigned FASTCALL FitchBoth_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned *weights, unsigned nBlocks) {
	unsigned thisScore = 0, i;
	unsigned char newVal;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < nBlocks; ++i) {
		if (!(newVal = s1[i] & s2[i])) {
			// No intersection of base sets: form union
			newVal = s1[i] | s2[i];
			thisScore += weights[i * 2];
		}
		
		dest[i] = newVal;
	}
#endif	// EMPTYFITCH
	return thisScore;
}

void FASTCALL FitchBases_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned nBlocks) {
	unsigned i;
	unsigned char newVal;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < nBlocks; ++i) {
		if (!(newVal = s1[i] & s2[i])) {
			// No intersection of base sets: form union
			newVal = s1[i] | s2[i];
		}
		
		dest[i] = newVal;
	}
#endif	// EMPTYFITCH
}

unsigned FASTCALL FitchScoreWeight1_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned nBlocks) {
	unsigned thisScore = 0, i;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < nBlocks; ++i) {
		if (!(s1[i] & s2[i])) {
			// No intersection of base sets: form union
			++thisScore;
		}
	}
#endif	// EMPTYFITCH
	return thisScore;
}

unsigned FASTCALL FitchScore_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len) {
	unsigned thisScore = 0, i;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < len; ++i) {
		if (!(s1[i] & s2[i])) {
			// No intersection of base sets: form union
			thisScore += weights[i];
		}
	}
#endif	// EMPTYFITCH
	return thisScore;
}




//TODO: Rename/fix the functions below.

// The non-SQUEEZEBASES half of the function formerly known as Fitch2Seqs_BASIC().
unsigned FASTCALL Fitch2Seqs_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned *weights, unsigned len) {
	unsigned thisScore = 0, i;
	char newVal;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < len; ++i) {
		if (!(newVal = s1[i] & s2[i])) {
			// No intersection of base sets: form union
			newVal = s1[i] | s2[i];
			thisScore += weights[i];
		}
		
		dest[i] = newVal;
	}
#endif	// EMPTYFITCH
	return thisScore;
}

// This operates like Fitch2SeqsCost_b1_w1_basic(), but incurs cost whenever
// corresponding base sets are not exactly equal.
unsigned FASTCALL StrictlyEqualFitch2SeqsCost_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len) {
	unsigned thisScore = 0, i;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	for (i = 0; i < len; ++i) {
		if (s1[i] != s2[i]) {
			// Base sets are not equal: therefore there *could* be a unit cost
			// incurred along this edge.
			thisScore += weights[i];
		}
	}
#endif	// EMPTYFITCH
	return thisScore;
}

// (Hopefully) slightly faster version in C.
unsigned FASTCALL Fitch2Seqs_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned *weights, unsigned len) {
	unsigned thisScore = 0, i;
	unsigned *dws1, *dws2, *dwdest;
	unsigned andVal, orVal, maskVal, newVal;
	unsigned hiBitsOn;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	dws1 = (unsigned *) s1;
	dws2 = (unsigned *) s2;
	dwdest = (unsigned *) dest;
	len >>= 2;
	
	for (i = 0; i < len; ++i) {
		andVal = dws1[i] & dws2[i];
		hiBitsOn = ((andVal + 0x0F0F0F0F) >> 4) & 0x01010101;
		orVal = dws1[i] | dws2[i];
		maskVal = (hiBitsOn + 0x7F7F7F7F) ^ 0x80808080;
		newVal = andVal | (orVal & maskVal);
//			thisScore += weights[i * 4]     & -(maskVal & 1);
//			thisScore += weights[i * 4 + 1] & -((maskVal >> 8) & 1);
//			thisScore += weights[i * 4 + 2] & -((maskVal >> 16) & 1);
//			thisScore += weights[i * 4 + 3] & -((maskVal >> 24) & 1);
#ifdef BIGENDIAN
		thisScore += ((unsigned) weights[i * 4])     & -((int) ((maskVal >> 24) & 1));
		thisScore += ((unsigned) weights[i * 4 + 1]) & -((int) ((maskVal >> 16) & 1));
		thisScore += ((unsigned) weights[i * 4 + 2]) & -((int) ((maskVal >> 8) & 1));
		thisScore += ((unsigned) weights[i * 4 + 3]) & -((int) (maskVal & 1));
#else	// not BIGENDIAN
		thisScore += ((unsigned) weights[i * 4])     & -((int) (maskVal & 1));
		thisScore += ((unsigned) weights[i * 4 + 1]) & -((int) ((maskVal >> 8) & 1));
		thisScore += ((unsigned) weights[i * 4 + 2]) & -((int) ((maskVal >> 16) & 1));
		thisScore += ((unsigned) weights[i * 4 + 3]) & -((int) ((maskVal >> 24) & 1));
#endif	// BIGENDIAN
		dwdest[i] = newVal;
	}
#endif	// EMPTYFITCH
	return thisScore;
}

// (Hopefully) slightly faster version in C.
unsigned FASTCALL Fitch2SeqsCost_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len) {
	unsigned thisScore = 0, i;
	unsigned *dws1, *dws2;
	unsigned andVal, maskVal;
	unsigned hiBitsOn;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	dws1 = (unsigned *) s1;
	dws2 = (unsigned *) s2;
	len >>= 2;
	
	for (i = 0; i < len; ++i) {
		andVal = dws1[i] & dws2[i];
		hiBitsOn = ((andVal + 0x0F0F0F0F) >> 4) & 0x01010101;
		maskVal = (hiBitsOn + 0x7F7F7F7F) ^ 0x80808080;
//			thisScore += weights[i * 4]     & -(maskVal & 1);
//			thisScore += weights[i * 4 + 1] & -((maskVal >> 8) & 1);
//			thisScore += weights[i * 4 + 2] & -((maskVal >> 16) & 1);
//			thisScore += weights[i * 4 + 3] & -((maskVal >> 24) & 1);
#ifdef BIGENDIAN
		thisScore += ((unsigned) weights[i * 4])     & -((int) ((maskVal >> 24) & 1));
		thisScore += ((unsigned) weights[i * 4 + 1]) & -((int) ((maskVal >> 16) & 1));
		thisScore += ((unsigned) weights[i * 4 + 2]) & -((int) ((maskVal >> 8) & 1));
		thisScore += ((unsigned) weights[i * 4 + 3]) & -((int) (maskVal & 1));
#else	// not BIGENDIAN
		thisScore += ((unsigned) weights[i * 4])     & -((int) (maskVal & 1));
		thisScore += ((unsigned) weights[i * 4 + 1]) & -((int) ((maskVal >> 8) & 1));
		thisScore += ((unsigned) weights[i * 4 + 2]) & -((int) ((maskVal >> 16) & 1));
		thisScore += ((unsigned) weights[i * 4 + 3]) & -((int) ((maskVal >> 24) & 1));
#endif	// BIGENDIAN
#ifdef OHGODFASTCFITCHISSCREWED
		if (thisScore != Fitch2Seqs_BASIC(s1, s2, dest, weights, (i + 1) * 4)) {
			DBGPRINT2("OH MY GOD: FASTCFITCH has failed while processing DWORD #%u.\n", i + 1);
			DBGPRINT5("s1: %c%c%c%c\n", GetBaseFromMask(s1[0]), GetBaseFromMask(s1[1]), GetBaseFromMask(s1[2]), GetBaseFromMask(s1[3]));
			DBGPRINT5("s2: %c%c%c%c\n", GetBaseFromMask(s2[0]), GetBaseFromMask(s2[1]), GetBaseFromMask(s2[2]), GetBaseFromMask(s2[3]));
			DBGPRINT5("Weights: %u %u %u %u\n", weights[0], weights[1], weights[2], weights[3]);
			DBGPRINT2("s1 in hex: %X\n", *((unsigned *) s1));
			DBGPRINT2("s2 in hex: %X\n", *((unsigned *) s2));
			DBGPRINT2("Score according to FASTCFITCH: %u.\n", thisScore);
			DBGPRINT2("Score according to SLOWCFITCH: %u.\n", Fitch2Seqs_BASIC(s1, s2, dest, weights, (i + 1) * 4));
			exit(42);
		}
#endif	// OHGODFASTCFITCHISSCREWED
	}
#endif	// EMPTYFITCH
	return thisScore;
}

// There is no basic implementation of the ...Weight1() functions.
unsigned FASTCALL Fitch2SeqsWeight1_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned len) {
	unsigned thisScore = 0, i;
	unsigned *dws1, *dws2, *dwdest;
	unsigned andVal, orVal, maskVal, newVal;
	unsigned hiBitsOn;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	dws1 = (unsigned *) s1;
	dws2 = (unsigned *) s2;
	dwdest = (unsigned *) dest;
	len >>= 2;
	
	for (i = 0; i < len; ++i) {
		andVal = dws1[i] & dws2[i];
		hiBitsOn = ((andVal + 0x0F0F0F0F) >> 4) & 0x01010101;
		orVal = dws1[i] | dws2[i];
		maskVal = (hiBitsOn + 0x7F7F7F7F) ^ 0x80808080;
		newVal = andVal | (orVal & maskVal);
#ifndef HACKDONTCOMPUTESCORE
		thisScore += ((maskVal & 0x01010101) * 0x01010101) >> 24;
#endif	// HACKDONTCOMPUTESCORE
		dwdest[i] = newVal;
	}
#endif	// EMPTYFITCH
	return thisScore;
}

unsigned FASTCALL Fitch2SeqsWeight1Cost_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned len) {
	unsigned thisScore = 0, i;
	unsigned *dws1, *dws2;
	unsigned andVal;
	unsigned hiBitsOn;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	dws1 = (unsigned *) s1;
	dws2 = (unsigned *) s2;
	len >>= 2;
	
	for (i = 0; i < len; ++i) {
		andVal = dws1[i] & dws2[i];
		hiBitsOn = (((andVal + 0x0F0F0F0F) >> 4) & 0x01010101) ^ 0x01010101;	// We reverse the sense of the LSB in each byte
#ifndef HACKDONTCOMPUTESCORE
		thisScore += (hiBitsOn * 0x01010101) >> 24;
#endif	// HACKDONTCOMPUTESCORE
	}
#endif	// EMPTYFITCH
	return thisScore;
}

// This can only be used when SQUEEZEBASES is not in use.  There is no "basic" implementation.
unsigned FASTCALL Fitch2SeqsLowWeight_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned char *byteWeights, unsigned len) {
	unsigned thisScore = 0, i;
	unsigned *dws1, *dws2, *dwdest;
	unsigned andVal, orVal, maskVal, newVal;
	unsigned hiBitsOn;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	dws1 = (unsigned *) s1;
	dws2 = (unsigned *) s2;
	dwdest = (unsigned *) dest;
	len >>= 2;
	
	for (i = 0; i < len; ++i) {
		andVal = dws1[i] & dws2[i];
		hiBitsOn = ((andVal + 0x0F0F0F0F) >> 4) & 0x01010101;
		orVal = dws1[i] | dws2[i];
		maskVal = (hiBitsOn + 0x7F7F7F7F) ^ 0x80808080;
		newVal = andVal | (orVal & maskVal);
#ifndef HACKDONTCOMPUTESCORE
//			thisScore += ((maskVal & 0x01010101) * 0x01010101) >> 24;
		thisScore += ((maskVal & ((unsigned *) byteWeights)[i]) * 0x01010101) >> 24;
#endif	// HACKDONTCOMPUTESCORE
		dwdest[i] = newVal;
	}
#endif	// EMPTYFITCH
	return thisScore;
}

// This can only be used when SQUEEZEBASES is not in use.
unsigned FASTCALL Fitch2SeqsLowWeightCost_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned char *byteWeights, unsigned len) {
	unsigned thisScore = 0, i;
	unsigned *dws1, *dws2;
	unsigned andVal, maskVal;
	unsigned hiBitsOn;
	
	INCMEASURE(fitch2SeqsCount);
	
#ifndef EMPTYFITCH
	dws1 = (unsigned *) s1;
	dws2 = (unsigned *) s2;
	len >>= 2;
	
	for (i = 0; i < len; ++i) {
		andVal = dws1[i] & dws2[i];
		hiBitsOn = ((andVal + 0x0F0F0F0F) >> 4) & 0x01010101;
		maskVal = (hiBitsOn + 0x7F7F7F7F) ^ 0x80808080;
#ifndef HACKDONTCOMPUTESCORE
//			thisScore += ((maskVal & 0x01010101) * 0x01010101) >> 24;
		thisScore += ((maskVal & ((unsigned *) byteWeights)[i]) * 0x01010101) >> 24;
#endif	// HACKDONTCOMPUTESCORE
	}
#endif	// EMPTYFITCH
	return thisScore;
}
