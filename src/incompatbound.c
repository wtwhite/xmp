#include <stdio.h>
#include <stdlib.h>
//#include <limits.h>
//#include <time.h>
//#include <ctype.h>
//#include <mem.h>
#include <string.h>
#include <assert.h>
#include "common.h"
#include "seq.h"
#include <procargs.h>		// Part of my own library
#include <dbgprint.h>		// My DBGPRINTn() macros
#include "incompatbound.h"

int GreedyMaxMatchingSiteCompare(const void *a, const void *b);
int CompareIncompatPairs(const void *a, const void *b);
void SaveDimacsGraph(char *fname, unsigned *incompatArray, unsigned numInc, unsigned seqLen, unsigned weighted);

//HACK: god this is awful
// These scores are only additive with respect to the number-of-distinct-bases-
// minus-one bound, not the MST bound.  So if it is run, the only restBound[]-populating
// routine that can be run before it is InitDistBasesPerSiteBoundRest().
// The idea is that if two sites are incompatible with each other, but compatible
// with the tree built so far, then one of them must be duplicated by the time
// the entire tree is built.  If the characters "displayed" by the two sites,
// when considering just the first n taxa, are compatible with each other, then
// it is possible that a tree on n taxa displays both of them without duplication;
// more strongly, if these displayed characters are both parsimoniously uninformative,
// then *any* tree on n taxa must display both of them without duplication (since
// a non-PI site can never be duplicated) and so we can safely say that adding
// the remaining taxa must duplicate one of the sites.  This means that the
// lower bound on the score that must be added by the remaining taxa can safely
// be increased by 1).
void InitIncompatBoundRest_b1(struct TreeData *td) {
	unsigned i, j, k, m;
	unsigned *incompatArray;
	unsigned *addScore;
	unsigned numInc = 0;
	unsigned char *buf;
	int incPos;
	unsigned startPos;
	char fnameBuf[80];		//DEBUG
	int maxSitePairs = 0;

	// This can happen if the entire dataset is non-PI.  Necessary since td->seqLen is unsigned,
	// and later computations rely on td->seqLen - 1 being sensible.
	if (!td->seqLen) {
		return;
	}

	for (i = 0; i < td->seqLen; ++i) {
		maxSitePairs += td->weights[i];
	}
	incompatArray = malloc(maxSitePairs * (maxSitePairs - 1) / 2 * (sizeof (unsigned) * 4));

	buf = (unsigned char *) malloc(td->numTaxa * 2);
	
	for (i = 0; i < td->seqLen - 1; ++i) {
		for (j = i + 1; j < td->seqLen; ++j) {
			for (k = 0; k < td->numTaxa; ++k) {
				buf[k] = GetBaseAt_b1(k, i, td);
			}
			
			// Are sites i and j incompatible?
			if ((incPos = IsIncompatible_b1(i, j, td)) != -1) {
				// A site must be compatible with any tree on up to the first N
				// taxa, where N is the the maximum number of for which
				// the "displayed" character on the first N taxa is non-PI.
				
				// First need to transpose the characters.  NOTE:
				// ParsimoniouslyUninformativeWeight_b1t() is accustomed to
				// dealing with characters prior to SqueezeBases() being called,
				// so each character always occupies one byte.
				for (k = 0; k < td->numTaxa; ++k) {
					buf[k + td->numTaxa] = GetBaseAt_b1(k, j, td);
				}
				
				for (k = 2; k < td->numTaxa; ++k) {			// Must have more than 3 taxa to be P.I.
					// We have put the sites in a transposed buffer, so they
					// are now effectively sites 0 and 1.
					if (ParsimoniouslyUninformativeWeight_b1t(buf, k) == -1 || ParsimoniouslyUninformativeWeight_b1t(buf + td->numTaxa, k) == -1) break;
				}
				
				// Don't bother putting it in if k is too small
				if (k > 3) {
					//HACK: This is a pretty inefficient way to do this, but it should at least work...  And depending
					// on the way that GreedyMaxMatchingScore() scores site pairs, it will possibly provide better
					// (or at least greedier) choices of site pairs as it is possible to avoid using all of a site up
					// in one go.
					for (m = 0; m < (td->weights[i] < td->weights[j] ? td->weights[i] : td->weights[j]); ++m) {
	//					DBGPRINT4("Sites %u and %u are incompatible and compatible with any tree on up to %u taxa!\n", i + 1, j + 1, k - 1);
						incompatArray[numInc++] = i;
						incompatArray[numInc++] = j;
						incompatArray[numInc++] = k;
						incompatArray[numInc++] = CalcSitePairIncompatCost_b1t(buf, td->numTaxa);
					}
				}
			}
		}
	}
	
	free(buf);
	
	//HACK: It would be much cleaner to use an array-of-arrays -- then yucky hacks like the following line
	// wouldn't be necessary.  But things are mostly working now so there's not much point changing things.
	numInc /= 4;
	
	// Try to determine which incompatibilities to use.
	qsort(incompatArray, numInc, sizeof (unsigned) * 4, CompareIncompatPairs);
	
	addScore = malloc(td->numTaxa * sizeof (unsigned));
	for (i = 0; i < 4; ++i) {
		addScore[i] = 0;
	}
	
	startPos = 0;
	for (i = 4; i < td->numTaxa; ++i) {
		// For trees on up to n taxa, we can use all incompatible sites whose
		// characters are non-PI when "displayed" on the first n taxa.
		
		// Advance to start of characters which are non-PI for up to at least n taxa
		while (startPos < numInc && incompatArray[startPos * 4 + 2] <= i) {
			++startPos;
		}
		
//		DBGPRINT4("About to compute greedy maximum matching for trees of size %u (startPos=%u, numInc=%u).\n", i, startPos, numInc);
		//DEBUG
		sprintf(fnameBuf, "incompat%d.dimacs", i);
		SaveDimacsGraph(fnameBuf, incompatArray + startPos * 4, numInc - startPos, td->seqLen, 0);		// Use unweighted edges for now
		addScore[i] = GreedyMaxMatchingScore(incompatArray + startPos * 4, numInc - startPos, td->seqLen);
		DBGPRINT5("Score for greedy maximum matching for trees of size %u (startPos=%u, numInc=%u) is %u.\n", i, startPos, numInc, addScore[i]);
	}
	
	// 26/10/2009: Remember: restBound[3] is effectively added to the score of the initial tree PLUS the cost from wherever
	// taxon #4 (i.e. the 4th taxon) is added -- so it must consist of the minimum weight added by taxa
	// #5, #6, ..., #numTaxa.  If it included taxon #4 then we would be counting that taxon twice!  (Well, once for
	// real and once as an LB.)  It's confusing that it works this way, but the comment alongside its definition in
	// common.h is totally accurate.
	for (i = 1; i < td->numTaxa; ++i) {
		DBGPRINT3("Adding a total of %u to the bound for all trees with %u taxa.\n", addScore[i], i);
		td->restBound[i - 1] += addScore[i];
	}
	
	free(addScore);
	free(incompatArray);
}

//CHECK2009: A tricky one.  I think the termination condition of n_pairs <= 2
// is really just an optimisation, as every scenario I can think of would delete
// every sequence pair for a compatible site pair.
// According to the reticulate.c program of Ingrid Jakobsen, the way to
// calculate incompatibility when more than two character states are allowed
// is as follows:
// 1.  Determine all distinct pairs of character states from the two sites.
// 2.  Examine the current list of pairs, and eliminate any pairs for which
//     the state for either (or both) of the sites is unique within the list of
//     states at that site.
// 3.  If one or more pairs of states was removed, go back to 2.
// 4.  If no pairs of states were removed, examine the number of pairs remaining
//     in the list:
//     - If two or fewer pairs remain, the sites are compatible
//     - Otherwise, the sites are incompatible.
// I will augment this by considering what happens when ambiguous sets of
// character states are present -- in this case, if the state pair you are
// currently checking for uniqueness could POSSIBLY be unique, eliminate it.
// E.g. if you have this data:
//
//   1    2
//   A    A
//   A    C
//   C    A
//   [AC] C
//
// Then when checking the third state pair, C in site 1 could be unique if
// A is chosen for site 1 of the fourth pair, so it can be removed.
int IsIncompatible_b1(unsigned site1, unsigned site2, struct TreeData *td) {
	unsigned i, j, k;
	unsigned char fullIb1, fullIb2, ib1, ib2, jb1, jb2;
	unsigned *live;
	unsigned changes, numLive, unique1, unique2;
#ifdef DEBUG
	unsigned ambiguous = 0;
#endif	// DEBUG
	
	live = (unsigned *) malloc(td->numTaxa * sizeof (unsigned));
	numLive = 0;
	for (i = 0; i < td->numTaxa; ++i) {
		// Is this the first incidence of this base-pair?
		for (j = 0; j < i; ++j) {
			if (GetBaseAt_b1(i, site1, td) == GetBaseAt_b1(j, site1, td) && GetBaseAt_b1(i, site2, td) == GetBaseAt_b1(j, site2, td)) break;
		}
		
		// This way we only make live one of every set of identical pairs of bases.
		if (j == i) {
			live[i] = 1;
			++numLive;
		} else {
			live[i] = 0;
		}
	}
	
	do {
		changes = 0;
		for (i = 0; i < td->numTaxa; ++i) {
			if (!live[i]) continue;
			
			fullIb1 = GetBaseAt_b1(i, site1, td);
			fullIb2 = GetBaseAt_b1(i, site2, td);
			
#ifdef DEBUG
			if ((fullIb1 != BM_A &&
				fullIb1 != BM_C &&
				fullIb1 != BM_G &&
				fullIb1 != BM_T) ||
				(fullIb2 != BM_A &&
				fullIb2 != BM_C &&
				fullIb2 != BM_G &&
				fullIb2 != BM_T)) ambiguous = 1;
#endif	// DEBUG
			
			// In the case of ambiguous base sets, try every possible
			// "assignment" of each base
			for (k = 1; k <= 16; k <<= 1) {
				ib1 = fullIb1 & k;
				ib2 = fullIb2 & k;
				if (ib1 == 0 && ib2 == 0) continue;		// Just saves time
				
				// Is either base in this pair unique for its site?
				unique1 = unique2 = 1;
				for (j = 0; j < td->numTaxa; ++j) {
					if (!live[j] || j == i) continue;
					
					if (ib1) {
						jb1 = GetBaseAt_b1(j, site1, td);
						
						// ib1 is now always a single base, so if jb1 is an
						// ambiguous base set, then a different base can always
						// be chosen from it to preserve the uniqueness of ib1.
						if (ib1 == jb1) {
							unique1 = 0;
						}
					}
					if (ib2) {
						jb2 = GetBaseAt_b1(j, site2, td);
						if (ib2 == jb2) {
							unique2 = 0;
						}
					}
				}
				
				// Was a unique pair found?
				if ((ib1 && unique1) || (ib2 && unique2)) {
					// Remove this pair from the list
					live[i] = 0;
					++changes;
					--numLive;
					break;
				}
			}
		}
	} while (changes && numLive > 2);
	
	free(live);
	
	if (numLive <= 2) {
#ifdef DEBUG
		if (ambiguous) {
//			DBGPRINT3("Sites %u and %u contain ambiguous base sets, and I've found them to be compatible.\n", site1 + 1, site2 + 1);
		}
#endif	// DEBUG
		return -1;			// Compatible
	} else {
#ifdef DEBUG
		if (ambiguous) {
//			DBGPRINT4("Sites %u and %u contain ambiguous base sets, but I've found them to be incompatible anyway (numLive=%u)!\n", site1 + 1, site2 + 1, numLive);
		}
#endif	// DEBUG
		return 0;			// Incompatible
	}
}

//CHECK2009: I'm certain that we can in fact simply remove pairs containing ambiguous bases and
// still obtain a valid lower bound on the cost of the site pair -- you can never decrease the
// score of a partition by adding sequences!  So this will yield better LBs.  Should also
// rename the function to explicitly show that we compute an LB on the cost, not necessarily the exact cost.
//HACK: If the site contains ambiguous base sets, return a fixed cost of 1
// (which is always a conservative estimate, since we know that the two sites
// ARE incompatible).
unsigned CalcSitePairIncompatCost_b1t(unsigned char *site, unsigned numTaxa) {
	unsigned distBound1 = SiteWeightLowerBound_b1t(site, numTaxa);
	unsigned distBound2 = SiteWeightLowerBound_b1t(site + numTaxa, numTaxa);
	unsigned index = 0, ambig = 0, rawCost, cost;
	int base1, base2;
	unsigned i;
	static int map[16] = {
		-1,	// 0
		0,	// 1
		1,	// 2
		-1,	// 3
		2,	// 4
		-1,	// 5
		-1,	// 6
		-1,	// 7
		3,	// 8
		-1,	// 9
		-1,	// 10
		-1,	// 11
		-1,	// 12
		-1,	// 13
		-1,	// 14
		-1,	// 15
	};
	
	// Compute the 16-bit index into the sitePairCost[] array.
	for (i = 0; i < numTaxa; ++i) {
		base1 = map[site[i]];
		base2 = map[site[i + numTaxa]];
		
		if (base1 >= 0 && base2 >= 0) {
			index |= 1 << (base1 * 4 + base2);
		} else {
			//HACK: can't deal with this yet
			// WTJW 22/3/2005: I have a feeling that we could just ignore the
			// ambiguous base pair, effectively removing this taxon from
			// consideration, but I haven't thought this through 100%.
			ambig = 1;
			break;
		}
	}
	
	if (ambig) {
		return 1;		// A conservative estimate of the cost, since we know that the two sites ARE incompatible
	} else {
		rawCost = sitePairCost[index - 1];
		cost = rawCost - distBound1 - distBound2;
		
		return cost;
	}
}

unsigned GreedyMaxMatchingScore(unsigned *incompatArray, unsigned numInc, unsigned numSites) {
	unsigned *sortedArray = (unsigned *) malloc(numInc * 4 * sizeof (unsigned));
	unsigned *siteArray = (unsigned *) malloc(numSites * sizeof (unsigned));
	unsigned i, edges, writeTo, site;
	
	DBGPRINT2("GreedyMaxMatchingScore() called with %u incompatible pairs of sites.\n", numInc);
	
	memset(siteArray, 0, numSites * sizeof (unsigned));
	
	// Each node in the "graph" corresponds to a site.  Determine the degree of each.
	for (i = 0; i < numInc; ++i) {
		assert(incompatArray[i * 4] < numSites);		//DEBUG
		assert(incompatArray[i * 4 + 1] < numSites);		//DEBUG
		
		++siteArray[incompatArray[i * 4]];
		++siteArray[incompatArray[i * 4 + 1]];
		sortedArray[i * 4] = incompatArray[i * 4];
		sortedArray[i * 4 + 1] = incompatArray[i * 4 + 1];
		sortedArray[i * 4 + 3] = incompatArray[i * 4 + 3];		// The weight of this incompatibility (>= 1)
	}
	
	// Now edges at the top of the list are those that are least likely to
	// interfere with other edges (I hope...)
	edges = 0;
	while (numInc > edges) {
		// Determine a "score" for each edge based on the degree of the nodes it
		// connects together.  You could use e.g. the product or the sum of the
		// degrees.
		for (i = 0; i < numInc; ++i) {
			sortedArray[i * 4 + 2] = siteArray[sortedArray[i * 4]] + siteArray[sortedArray[i * 4 + 1]];
		}
		
		// So that sortedArray will live up to its name
		qsort(sortedArray, numInc, 4 * sizeof (unsigned), GreedyMaxMatchingSiteCompare);
		
		// Take the edge off the top and add it to the matching (we don't
		// actually record the edges in the matching, only their scores).
		++edges;
		
		// Now remove all other edges that used either endpoint of this edge
		writeTo = edges;
		for (i = edges; i < numInc; ++i) {
			//HACK: since the lower-numbered site is always stored first, we
			// may be able to get away with fewer comparisons than this
			if (sortedArray[i * 4] != sortedArray[(edges - 1) * 4] &&
				sortedArray[i * 4] != sortedArray[(edges - 1) * 4 + 1] &&
				sortedArray[i * 4 + 1] != sortedArray[(edges - 1) * 4] &&
				sortedArray[i * 4 + 1] != sortedArray[(edges - 1) * 4 + 1]) {
				// This edge does not use either endpoint of the edge just added,
				// so it can be retained
				
				// Only bother copying if we actually have to
				if (i != writeTo) {
					memcpy(sortedArray + writeTo * 4, sortedArray + i * 4, 4 * sizeof (unsigned));
				}
				
				++writeTo;
			} else {
				// This edge had as an endpoint one of the endpoints of the edge
				// just added to the matching, so it cannot be added.  Reduce
				// the degree of the other vertex in this edge accordingly.
				if (sortedArray[i * 4] == sortedArray[(edges - 1) * 4] ||
					sortedArray[i * 4] == sortedArray[(edges - 1) * 4 + 1]) site = sortedArray[i * 4];
				if (sortedArray[i * 4 + 1] == sortedArray[(edges - 1) * 4] ||
					sortedArray[i * 4 + 1] == sortedArray[(edges - 1) * 4 + 1]) site = sortedArray[i * 4 + 1];
				
				--siteArray[site];
			}
		}
		
		numInc = writeTo;
	}
	
	free(sortedArray);
	free(siteArray);
	return edges;
}

// Determine a "score" for each edge based on the degree of the nodes it
// connects together.  You could use e.g. the product or the sum of the
// degrees.  Puts the list in ascending order.
int GreedyMaxMatchingSiteCompare(const void *a, const void *b) {
	if (*((unsigned *) a + 3) == *((unsigned *) b + 3)) {
		// Weights are equal: compare on "interference" value
		return *((unsigned *) a + 2) - *((unsigned *) b + 2);
	} else {
		// Prefer the edge with highest weight
		//CHECK2009: I'm not sure whether this is necessarily a good idea...
		return *((unsigned *) a + 3) - *((unsigned *) b + 3);
	}
}

// Order sites so that those sites that become parsimoniously informative the earliest come first.
int CompareIncompatPairs(const void *a, const void *b) {
//	return ((unsigned *) b)[2] - ((unsigned *) a)[2];
	return ((unsigned *) a)[2] - ((unsigned *) b)[2];
}

//HACK: weighted is currently ignored (don't know what to do with it)
void SaveDimacsGraph(char *fname, unsigned *incompatArray, unsigned numInc, unsigned seqLen, unsigned weighted) {
	FILE *f;
	unsigned i;
	
	if (!(f = fopen(fname, "w"))) {
		fprintf(stderr, "Unable to open file '%s' for output, aborting.\n", fname);
		exit(1);
	}
	
	// Write header line
	fprintf(f, "p clq %u %u\n", numInc, seqLen);
	
	for (i = 0; i < numInc; ++i) {
		fprintf(f, "e %u %u\n", incompatArray[i * 4], incompatArray[i * 4 + 1]);
	}
	
	fclose(f);
}
