#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "common.h"
#include <dbgprint.h>
#include "distbasesbound.h"
#include "incompatbound.h"

// Since we always add the taxa to the tree in the same order, and where we add
// them makes no difference to the number of distinct bases at any site, there
// are only numTaxa possible distinct values for the minimum-additional-weight:
// the m-a-w when no taxa are on the tree and numTaxa are not; the m-a-w when
// the first taxon is on the tree and the remaining (numTaxa-1) are not; the
// m-a-w when the first 2 taxa are on the tree and the remaining (numTaxa-2) are
// not; and so on.  Thus we can compute these m-a-w values at the start and just
// store them in an array of size numTaxa!
//
// td->restBound[x] = the minimum additional weight that must be added to the
// tree by taxa x+1, x+2, ..., numTaxa.  td->restBound[numTaxa] = 0, since
// at this point, all taxa have been added.
void InitDistBasesPerSiteBoundRest_b1_w1(struct TreeData *td) {
	unsigned i, j, k;
	unsigned *treeStateCounts[4];
	unsigned *restStateCounts[4];
	unsigned *ignore;
	unsigned distBasesPerSiteMinusOne, totalDistBasesPerSiteMinusOne;
	unsigned char c;
#ifdef DEBUG
	unsigned nBases;
	unsigned nBaseSiteFreqs[4] = { 0 }, nBaseWeightFreqs[4] = { 0 };
#endif	// DEBUG
#ifdef SITESCORES
	unsigned delta;
	unsigned *subtractFromBoundSiteScores;
	unsigned checkTotal, checkTotalBoundSiteScores, slack;		//DEBUG
#endif	// SITESCORES
	
	// Allocate and initialise state count arrays
	ignore = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
	memset(ignore, 0, td->seqLen * sizeof (unsigned));
	for (j = 0; j < 4; ++j) {
		treeStateCounts[j] = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
		restStateCounts[j] = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
		for (i = 0; i < td->seqLen; ++i) {
			treeStateCounts[j][i] = 0;
			restStateCounts[j][i] = 0;
		}
	}
	
	// Count the distinct bases at every site
	for (j = 0; j < td->numTaxa; ++j) {
		for (i = 0; i < td->seqLen; ++i) {
			c = td->charMat[j * td->seqLen + i];
			if      (c == BM_A) ++restStateCounts[0][i];
			else if (c == BM_C) ++restStateCounts[1][i];
			else if (c == BM_G) ++restStateCounts[2][i];
			else if (c == BM_T) ++restStateCounts[3][i];
			else ignore[i] = 1;			// Too hard to compute this bound for sites containing ambiguous bases, so just ignore them
		}
	}
	
#ifdef DEBUG
	// Just for fun: count the number of sites having {2, 3, 4} distinct bases
	for (i = 0; i < td->seqLen; ++i) {
		if (ignore[i]) continue;
		
		nBases = 0;
		for (j = 0; j < 4; ++j) {
			if (restStateCounts[j][i] > 0)
				++nBases;
		}
		
		++nBaseSiteFreqs[nBases - 1];
		nBaseWeightFreqs[nBases - 1] += td->weights[i];
	}
	
	for (j = 0; j < 4; ++j) {
		DBGPRINT4("DISTBASESBOUNDREST: there are %u sites (having total weight %u) %u distinct bases.\n", nBaseSiteFreqs[j], nBaseWeightFreqs[j], j + 1);
	}
#endif	// DEBUG
		
#ifdef SITESCORES
	subtractFromBoundSiteScores = (unsigned *) malloc(td->seqLen);
	memset(subtractFromBoundSiteScores, 0, td->seqLen);
#endif	// SITESCORES
	
	// Effectively add each taxon in turn to the tree
	for (j = 0; j < td->numTaxa; ++j) {
		// Count the distinct bases at every site
		distBasesPerSiteMinusOne = 0;
		totalDistBasesPerSiteMinusOne = 0;
		for (i = 0; i < td->seqLen; ++i) {
			c = td->charMat[j * td->seqLen + i];
			if (ignore[i]) continue;
			
			if (c == BM_A) { ++treeStateCounts[0][i]; --restStateCounts[0][i]; }
			if (c == BM_C) { ++treeStateCounts[1][i]; --restStateCounts[1][i]; }
			if (c == BM_G) { ++treeStateCounts[2][i]; --restStateCounts[2][i]; }
			if (c == BM_T) { ++treeStateCounts[3][i]; --restStateCounts[3][i]; }
			
			// For each distinct base at this site that is not on the tree but
			// is in the remaining taxa, count the weight of this site into
			// the bound for this number of taxa.
			distBasesPerSiteMinusOne = 0;
			for (k = 0; k < 4; ++k) {
				if (treeStateCounts[k][i] == 0 && restStateCounts[k][i] > 0) {
					distBasesPerSiteMinusOne += td->weights[i];
				}
			}
			
			totalDistBasesPerSiteMinusOne += distBasesPerSiteMinusOne;
#ifdef SITESCORES
			// Need to counterbalance the bound for this group of sites
			subtractFromBoundSiteScores[i / (4 * BASESPERBYTE)] += distBasesPerSiteMinusOne;
#endif	// SITESCORES
		}
#ifdef SITESCORES
		checkTotalBoundSiteScores = checkTotal = slack = 0;
//		for (i = 0; i < td->seqLen / (4 * BASESPERBYTE); ++i) {
		for (i = 0; i < td->seqLen / 4; ++i) {
			if (td->boundSiteScores[i] >= subtractFromBoundSiteScores[i]) {
				td->boundSiteScoresByTaxon[j][i] = td->boundSiteScores[i] - subtractFromBoundSiteScores[i];
			} else {
				td->boundSiteScoresByTaxon[j][i] = 0;			// i.e. saturate at zero
				DBGPRINT5("Whoa!  With %u taxa on the tree, and for sites %u-%u, the distinct-bases-per-site bound is %u better than the partition bound!\n", j, (i * 4 * BASESPERBYTE) + 1, (i * 4 * BASESPERBYTE) + 4 * BASESPERBYTE, subtractFromBoundSiteScores[i] - td->boundSiteScores[i]);
				slack += subtractFromBoundSiteScores[i] - td->boundSiteScores[i];
			}
			checkTotal += td->boundSiteScoresByTaxon[j][i];
			checkTotalBoundSiteScores += td->boundSiteScores[i];
			subtractFromBoundSiteScores[i] = 0;
		}
		
DBGPRINT2("checkTotal = <%u>\n", checkTotal);		//DEBUG

DBGPRINT2("slack = <%u>\n", slack);		//DEBUG

DBGPRINT2("checkTotalBoundSiteScores = <%u>\n", checkTotalBoundSiteScores);		//DEBUG

DBGPRINT2("totalDistBasesPerSiteMinusOne = <%u>\n", totalDistBasesPerSiteMinusOne);		//DEBUG

		assert(checkTotal - slack == checkTotalBoundSiteScores - totalDistBasesPerSiteMinusOne);		//DEBUG
#endif	// SITESCORES
		
		DBGPRINT3("totalDistBasesPerSiteMinusOne[%u] = <%u>\n", j, totalDistBasesPerSiteMinusOne);		//DEBUG
		if (totalDistBasesPerSiteMinusOne > td->restBound[j]) {
			DBGPRINT3("totalDistBasesPerSiteMinusOne beats current lower bound for trees of first %u taxa by %u\n", j + 1, totalDistBasesPerSiteMinusOne - td->restBound[j]);	//DEBUG
			td->restBound[j] = totalDistBasesPerSiteMinusOne;
		}
	}
	
	for (j = 0; j < 4; ++j) {
		free(treeStateCounts[j]);
		free(restStateCounts[j]);
	}
	free(ignore);
	
#ifdef SITESCORES
	free(subtractFromBoundSiteScores);
#endif	// SITESCORES
}
