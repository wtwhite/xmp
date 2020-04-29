#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
//#include <mem.h>
#include <string.h>
#include <assert.h>
#include "common.h"
#include "seq.h"
#include "yanbader.h"
#include <procargs.h>		// Part of my own library
#include <dbgprint.h>		// My DBGPRINTn() macros
#include "hittingsetbound.h"

int ColumnCompare(const void *a, const void *b);

// This is called when the -Bh option (OPT_USEGREEDYHITTINGSETBOUND) is given.
// We create n(n-1)/2 sets, one for each pair of taxa.  Each set contains a
// list of the sites that distinguish (separate) those two taxa.  We then try
// to greedily cover these sets with sites.  If we manage to cover every set
// with k sites, then we have a k-site part of the dataset with all taxa
// distinct, so we can apply the 1-connected bound.
// Starting with an empty set of sites, we pick as the next site to add the
// site which hits the most sets, and which.
// In particular, a site must hit more sets than the upper bound on its cost;
// since otherwise, it is no more useful than when it is considered as its own
// part of size 1 (in the case of equality) and could be worse (in the case
// where its cost upper bound exceeds the number of sets it hits).
//
// Handling duplicate columns: The current scheme will work when there are
// multiple identical sites -- as soon as one of a set of identical sites has
// been added to the cover, all of the other sites in that set will hit zero
// pairs in future rounds and consequently won't be selected.  We keep a copy
// of the td->weights[] array and decrement a site weight whenever a site
// is added to a cover; when a site's weight is zero, it is not eligible for
// being selected.  This is potentially a slow way to do things (when only a
// few sites remain, we still accumulate totals into hits[] for all sites) but
// it does have advantages -- notably, site positions do not change so the
// lists of sites in sets[] does not have to be updated at any stage.  (An
// alternative would be to have each element of sets[] be a pointer, and
// maintain a linked list of sites from which sites are deleted as they are
// added to the cover.)
void InitGreedyHittingSetBoundRest(struct TreeData *td) {
	unsigned i, j, k, n, numSets;
	unsigned *sets;
	unsigned *setSizes;
	unsigned *covered;
	int *hits;
	unsigned *inCover;
	unsigned *weights;
	unsigned setNo;				//HACK: For intermediate computations.
	unsigned sitesRemaining;
	unsigned pairsRemaining;
	unsigned bestSite;
	int bestHits;
	unsigned *bestSiteByUB;
	int *bestHitsByUB;
	unsigned char *colBuf;
	unsigned char *buf;
	unsigned *allOnes;
	int oneConnBound, oneConnFat, totalFat;
	unsigned ka2Weight, totalUpperBound, sitesInCover, numSitesUsed;
	unsigned *upperBounds;
	unsigned HACK;		//DEBUG
	unsigned temp;
	unsigned *upperBoundWeights, *lowerBoundWeights;
	
	DBGPRINT1("Entering InitGreedyHittingSetBoundRest().  Initial dataset:\n");
	for (i = 0; i < td->numTaxa; ++i) {
		for (j = 0; j < td->seqLen; ++j) {
			putc(GetBaseFromMask(GetBaseAt(i, j, td)), dbgfile);
		}
//		PrintMaskSeq_b1(data + i * width, dataWidth, stdout);
		DBGPRINT1("\n");
	}
	
	// Determine lists of sites that distinguish each taxon pair
	//HACK: Could be using only n(n-1)/2 elements, but that makes indexing harder...
	numSets = td->numTaxa * td->numTaxa;
	sets = (unsigned *) malloc(numSets * td->seqLen * sizeof (unsigned));			//HACK: this could potentially burn some serious memory, but for 50 taxa and 100 sites it's only 1Mb
	setSizes = (unsigned *) malloc(numSets * sizeof (unsigned));
	
	for (i = 0; i < numSets; ++i) {
		setSizes[i] = 0;
	}
	
	//HACK: possibly more efficient to run through sites in the outermost loop
	// but this gets the job done too.
	for (i = 0; i < td->numTaxa - 1; ++i) {
		for (j = i + 1; j < td->numTaxa; ++j) {
			setNo = i * td->numTaxa + j;
			for (k = 0; k < td->seqLen; ++k) {
				// Does site k distinguish taxa i and j?
				if (!(GetBaseAt(i, k, td) & GetBaseAt(j, k, td))) {			//HACK: we should probably exclude sites containing ambiguous bases altogether -- will need to think about this.
					sets[setNo * td->seqLen + setSizes[setNo]++] = k;
				}
			}
		}
	}
	
	hits = (int *) malloc(td->seqLen * sizeof (unsigned));			// Needs to be int so that max-finding loop can use -1 as a sentinel
	inCover = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
	weights = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
	upperBounds = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
	upperBoundWeights = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
	lowerBoundWeights = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
	covered = (unsigned *) malloc(numSets * sizeof (unsigned));
	
	buf = (unsigned char *) malloc(td->numTaxa * td->seqLen);
	colBuf = (unsigned char *) malloc(td->numTaxa);
	
	// The upper bound weight of a site can be at most the number of taxa minus 1.
	bestSiteByUB = (unsigned *) malloc(td->numTaxa * sizeof (unsigned));
	bestHitsByUB = (int *) malloc(td->numTaxa * sizeof (unsigned));
	
	// Further down, KA2PartitionAndScoreColumns() requires a weight vector.  This will be all ones.
	allOnes = (unsigned *) malloc(td->seqLen * sizeof (unsigned));
	for (i = 0; i < td->seqLen; ++i) {
		allOnes[i] = 1;
	}
	
	for (n = 3; n < td->numTaxa; ++n) {
		DBGPRINT2("Starting construction of site partition for %d taxa.\n", n);
		
		// Find a greedy cover for the first n taxa.
		//HACK: we are currently creating the hitting set in the exact same way
		// on each iteration of this loop!  Need to adjust this so that the
		// upper bound of the cost of each site on the first n taxa is factored
		// into choosing which site to add next.
		sitesRemaining = 0;
		for (i = 0; i < td->seqLen; ++i) {
			weights[i] = td->weights[i];		//HACK: at this point we should set the weights of sites containing ambiguous bases to zero.  Maybe.
			sitesRemaining += weights[i];
		}
		
		DBGPRINT1("Populating ubw[] and lbw[]...\n");
		for (i = 0; i < td->seqLen; ++i) {
			for (j = 0; j < td->numTaxa; ++j) {
				colBuf[j] = GetBaseAt(j, i, td);						// For SiteWeightUpperBound_b1t()
			}
			
			upperBoundWeights[i] = SiteWeightUpperBound_b1t(colBuf, n);
			lowerBoundWeights[i] = SiteWeightLowerBound_b1t(colBuf, n);		// Number of distinct bases minus 1
//DBGPRINT2("upperBoundWeights[i] = <%u>\n", upperBoundWeights[i]);		//DEBUG
//DBGPRINT2("lowerBoundWeights[i] = <%u>\n", lowerBoundWeights[i]);		//DEBUG
		}
		
//DBGPRINT2("0 upperBoundWeights[0] = <%u>\n", upperBoundWeights[0]);		//DEBUG
		// Each iteration of this while loop constructs a subset of sites
		// (i.e. a part in the partition).
		while (sitesRemaining) {
			DBGPRINT2("Starting construction of site subset with %d sites in hand.\n", sitesRemaining);
			
			for (i = 0; i < td->seqLen; ++i) {
				inCover[i] = 0;
			}
			
			for (i = 0; i < numSets; ++i) {
				covered[i] = 0;
			}
			
			sitesInCover = 0;
			pairsRemaining = td->numTaxa * (td->numTaxa - 1) / 2;
//DBGPRINT2("1 upperBoundWeights[0] = <%u>\n", upperBoundWeights[0]);		//DEBUG
			while (sitesRemaining && pairsRemaining && sitesInCover < 6) {			//HACK
				// Find a site which separates the maximum number of pairs.
				for (i = 0; i < td->seqLen; ++i) {
					hits[i] = 0;
				}
				
				for (i = 0; i < td->numTaxa - 1; ++i) {
					for (j = i + 1; j < td->numTaxa; ++j) {
						setNo = i * td->numTaxa + j;
						
						// Following loop is a no-op for sets already hit
						if (!covered[setNo]) {
							for (k = 0; k < setSizes[setNo]; ++k) {
								++hits[sets[setNo * td->seqLen + k]];
							}
						}
					}
				}
				
//DBGPRINT2("2 upperBoundWeights[0] = <%u>\n", upperBoundWeights[0]);		//DEBUG
				// Find a site which separates the maximum number of pairs.
				//HACK: not yet enforcing constraint that this site have minimum possible
				// upper bound weight.
				// Now we find the best site for each upper bound weight.
				for (i = 0; i < td->numTaxa; ++i) {
					bestHitsByUB[i] = -1;
				}
				
//		DBGPRINT1("Finding best site in each UB class...\n");
				for (i = 0; i < td->seqLen; ++i) {
					// Only consider sites not already in the cover/hitting set, and
					// which have not been "used up" in other parts yet.
					// Notice that we are comparing site i only to other sites in its UB weight class.
//DBGPRINT2("upperBoundWeights[i] = <%u>\n", upperBoundWeights[i]);		//DEBUG
					if (weights[i] && !inCover[i] && hits[i] > bestHitsByUB[upperBoundWeights[i]]) {
						bestSiteByUB[upperBoundWeights[i]] = i;
						bestHitsByUB[upperBoundWeights[i]] = hits[i];
					}
				}
				
				// The best site overall is the one with the lowest upper bound
				// weight and a number of hits greater than its lower bound
				// weight (that being the number of distinct bases minus 1) and
				// also greater than 1.
				bestSite = 0;
				bestHits = -1;
				for (i = 0; i < td->numTaxa; ++i) {
					if (bestHitsByUB[i] > 1) {
						if (bestHitsByUB[i] > lowerBoundWeights[bestSiteByUB[i]]) {
							bestSite = bestSiteByUB[i];
							bestHits = bestHitsByUB[i];
							break;
						}
					}
				}
				
				if (bestHits <= 0) {			// Tricky!
					// No more improvements can be made!
					DBGPRINT1("No more improvements can be made.\n");
					break;
				}
				
				DBGPRINT5("The best site is site #%d, which covers %d pairs and has lower bound %d and upper bound %d.\n", bestSite, bestHits, lowerBoundWeights[bestSite], upperBoundWeights[bestSite]);
				
				// Cover the relevant sets
				//HACK: doing this the slowest possible way at the mo...
				HACK = 0;		//DEBUG
				for (i = 0; i < td->numTaxa - 1; ++i) {
					for (j = i + 1; j < td->numTaxa; ++j) {
						setNo = i * td->numTaxa + j;
						if (!covered[setNo]) {
							// This pair is not already covered.
							// Does site bestSite distinguish taxa i and j?
							if (!(GetBaseAt(i, bestSite, td) & GetBaseAt(j, bestSite, td))) {			//HACK: we should probably exclude sites containing ambiguous bases altogether -- will need to think about this.
								covered[setNo] = 1;
								++HACK;
							}
						}
					}
				}
				
				assert(HACK == bestHits);		//DEBUG
				
				inCover[bestSite] = 1;
				--weights[bestSite];			// When this goes to zero, the site has been "used up."
				--sitesRemaining;
				++sitesInCover;
				pairsRemaining -= bestHits;
			}
			
			if (!sitesInCover) {
				// We weren't able to find any useful sites.
				DBGPRINT2("The remaining %d sites could not be used to produce a hitting set.  This is the end for this partition.\n", sitesRemaining);
				break;
			}
			
			// We have a hitting set.  Now find a LB on the parsimony score of this
			// set of sites, and take off the sum of the UBs of each participating
			// site on the 1st n taxa.  What remains is the "fat" that can be added
			// to the score of any tree on the 1st n taxa.
			
			totalUpperBound = 0;
			for (i = 0; i < td->seqLen; ++i) {
				upperBounds[i] = 0;
			}
			
			k = 0;
			for (i = 0; i < td->seqLen; ++i) {
				if (inCover[i]) {
					for (j = 0; j < td->numTaxa; ++j) {
						colBuf[j] = GetBaseAt(j, i, td);						// For SiteWeightUpperBound_b1t()
						buf[j * sitesInCover + k] = GetBaseAt(j, i, td);		// For KA2PartitionAndScoreColumns()
					}
					
					upperBounds[k] = SiteWeightUpperBound_b1t(colBuf, n);
					totalUpperBound += upperBounds[k];
					++k;
				}
			}
			
			DBGPRINT4("%d sites in the cover.  Upper bound on site cost up to %d taxa = %d.\n", sitesInCover, n, totalUpperBound);
			// Determine an LB for this site set on all taxa.
			if (!pairsRemaining) {
				// This makes it easy.  All taxa are distinct, so in the worst case
				// the distance from any taxon to its nearest neighbour is 1.  Thus
				// a 1-connected bound can be used (although it may be possible to
				// get a tighter LB by working a bit harder).
				DBGPRINT1("All pairs of taxa have been separated, so using the easy 1-conn bound.\n");
				oneConnBound = td->numTaxa - 1;
			} else {
				// Other possibilities: sort the taxa and count the number of
				// distinct rows.
				DBGPRINT2("%d pairs of taxa remain unseparated, so can't do the easy 1-conn bound!\n", pairsRemaining);
				
				// Sort the taxa.  The call to KA2PartitionAndScoreColumns() later on does not mind that the rows have been sorted.
				//HACK: Using the global variable _ColumnCompareWidth.
				_ColumnCompareWidth = sitesInCover;
				qsort(buf, td->numTaxa, sitesInCover, ColumnCompare);
				
				oneConnBound = 0;
				for (i = 1; i < td->numTaxa; ++i) {
					if (memcmp(buf + (i - 1) * sitesInCover, buf + i * sitesInCover, sitesInCover)) {
						++oneConnBound;
					}
				}
				
//				//DEBUG: Actually list the unseparated pairs of taxa
//				DBGPRINT1("The following pairs of taxa remain unseparated:\n");
//				for (i = 0; i < td->numTaxa - 1; ++i) {
//					for (j = i + 1; j < td->numTaxa; ++j) {
//						setNo = i * td->numTaxa + j;
//						
//						// Following loop is a no-op for sets already hit
//						if (!covered[setNo]) {
//							DBGPRINT3("(%u,%u)\n", i, j);
//						}
//					}
//				}
			}
			
			oneConnFat = oneConnBound - totalUpperBound;
			DBGPRINT3("1-conn bound for all taxa = %d.  After taking off UB on sites, this leaves: %d.\n", oneConnBound, oneConnFat);
			
			//HACK: We will ignore numSitesUsed for now.
			// NOTE: KA2PartitionAndScoreColumns() swaps columns around!  Luckily we don't need the data later on.
			ka2Weight = (int) KA2PartitionAndScoreColumns(buf, sitesInCover, sitesInCover, td->numTaxa, allOnes, upperBounds, &numSitesUsed);
			DBGPRINT2("KA2 bound for all taxa = %d.  (This is after taking off UB on sites.)\n", ka2Weight);
			
			if (ka2Weight > oneConnFat) {
				totalFat = ka2Weight;
			} else {
				totalFat = oneConnFat;
			}
			
			if (totalFat > (int) td->restBound[n - 1]) {
				DBGPRINT4("InitGreedyHittingSetBoundRest(): Bound of %u for %u taxa is now %u better than before!\n", totalFat, n, totalFat - (int) td->restBound[n - 1]);
				td->restBound[n - 1] = (unsigned) totalFat;
			}
		}
	}
	
	free(sets);
	free(setSizes);
	free(hits);
	free(inCover);
	free(weights);
	free(covered);
	free(colBuf);
	free(buf);
	free(allOnes);
	free(upperBounds);
	free(upperBoundWeights);
	free(lowerBoundWeights);
	free(bestSiteByUB);
	free(bestHitsByUB);
}

int ColumnCompare(const void *a, const void *b) {
	return memcmp(a, b, _ColumnCompareWidth);		//HACK: uses global _ColumnCompareWidth from common.c.  **Slightly** better than before.
}
