#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "common.h"
#include "partbound.h"
#include "seq.h"
#include "maxmin.h"
#include <dbgprint.h>		// Important: must #include switches.h (or common.h, which #includes it) before dbgprint.h!

#define PARTBLOCKSIZE 8		// Purely for efficiency reasons; not tied to BYTESPERBLOCK.

// Uncomment the following #define to turn on logging of every move considered by OptimisePartitionOnePass().  This can produce
// hundreds of Mb of logging output in files with names of the form "OptimisePartitionOnePass.%dparts.%dtaxa.%dtaxaOnTree.dump.txt".
// Note: Unfortunately this naming scheme is not sufficient to distinguish multiple passes over a partition that do
// not reduce the number of parts, so in this case earlier files will be overwritten.
//#define PARTBOUND_LOG_EVERY_MOVE_CONSIDERED

// Uncomment the following #define to turn on more verbose logging to stderr showing each improved move and each
// move taken.  Note that you need to have DEBUG #defined as well to get the DBGPRINTnn() macros to do anything.
//#define PARTBOUND_VERBOSE_LOGGING

// Uncomment the following #define to dump out partitions of the data, and interesting statistics about them,
// before and after calling OptimisePartition().
// The files have names "partDump.%dnToT.chunk%d.before.part" and "partDump.%dnToT.chunk%d.after.part".  The naming system
// is sufficient to distinguish all partitions produced during a run, so overwriting will not occur, although much
// output may be produced.
//#define PARTBOUND_DUMP_PARTITIONS

// The maximum number of sites that will be added to any part in a partition in a single move.  We need to reserve
// this many extra columns of space in each part to enable efficient score calculation.
#define PARTBOUND_MAX_SITES_MOVED 1

struct kruskal_edge {
	int w;
	int v[2];
};

struct kruskal_component {
	int p;			// Index of parent.  If comps[i] == i, then this is a root.
	int w;			// If a root, weight of this component; meaningless otherwise
};

struct compute_part_lower_bound_buffers {
	struct kruskal_component *comps;
	struct kruskal_edge *edges;
	int *uniqueSeqs;
	int *furthestPairSeen;
	int *nEdgesWithWeight;
	unsigned char *dataBuf;
};

//HACK: This global variable records a bitmask showing the bound types that were maximal in the last call
// to ComputePartLowerBound_preallocated_b1_w1().
int gBestBoundMask;

//HACK: This global variable is needed by _CompareMem() to know how wide the memory block is.
int _CompareMemWidth;

static int _CompareEdgesByWeight(const struct kruskal_edge *e1, const struct kruskal_edge *e2) {
	return e1->w - e2->w;
}

static int find_component(struct kruskal_component *comps, int v) {
	int u = v;
	int temp;

	while (comps[v].p != v) {
		v = comps[v].p;
	}

	// Path compression: make every node on the path to the root point directly to the root.
	// Guarantees O(nlog n) time for any sequence of n find_component() and union_components() operations.
	while (comps[u].p != v) {
		temp = comps[u].p;
		comps[u].p = v;
		u = temp;
	}

	return v;
}

// v1 and v2 must be component roots.
// w is the weight of the edge connecting them.
// Returns the name of the new component (although currently this isn't used).
//TODO: Use ranking to speed things up.
static int union_components(struct kruskal_component *comps, int v1, int v2, int w) {
	comps[v1].p = v2;
	comps[v2].w += comps[v1].w + w;
	return v2;
}

// This will work if the sequences in pt are not sorted, but it will then take O(nTaxa^2) time --
// if they are sorted, it will take only O(nUniq^2) time.
// The BT_NDISTINCT and BT_FURTHESTPAIR bounds are no longer used, since they are dominated by other bounds.
enum bound_types {
	BT_ONECONN,
	BT_3COMPS,
	BT_3COMPSSTEINER,
	BT_MST,
	MAX_BT			// Used for dimensioning arrays etc.
};

// Computes the desired "final" score according to the specified strategy.
static INLINE int selectFinalScore(enum partbound_strategy_type strategy, int lb, int ub, int additive) {
	if (strategy == PBS_CLIPPEDANDADDITIVE) {
		return iMax(additive, lb - ub);		// We assume additive >= 0 to avoid an extra test
	} else if (strategy == PBS_CLIPPED) {
		return iMax(0, lb - ub);
	} else {
		assert(strategy == PBS_LB);
		return lb;
	}
}

// Useful for diagnostics now that we have removed the calculated value from struct scored_partition.
static int getTotalClippedScore(struct scored_partition *ptn) {
	int i, total = 0;

	for (i = 0; i < ptn->nParts; ++i) {
		total += iMax(0, ptn->parts[i].scoreLB - ptn->parts[i].scoreUB);
	}

	return total;
}

void SetupBuffersForComputePartLowerBound(struct compute_part_lower_bound_buffers *buffers, int nTaxa, int maxDist) {
	// Set up fast union/find data structure
	buffers->comps = malloc(nTaxa * sizeof (struct kruskal_component));
	// Bucket sorting of edges means we need space for twice as many edges as there could possibly be.
	buffers->edges = malloc(nTaxa * (nTaxa - 1) * sizeof (struct kruskal_edge));
	buffers->uniqueSeqs = malloc(nTaxa * sizeof (int));
	buffers->furthestPairSeen = malloc(nTaxa * nTaxa * sizeof (int));
	buffers->nEdgesWithWeight = malloc((maxDist + 1) * sizeof (int));
	buffers->dataBuf = malloc(nTaxa * maxDist);
}

void TearDownBuffersForComputePartLowerBound(struct compute_part_lower_bound_buffers *buffers) {
	free(buffers->nEdgesWithWeight);
	free(buffers->furthestPairSeen);
	free(buffers->uniqueSeqs);
	free(buffers->edges);
	free(buffers->comps);
	free(buffers->dataBuf);
}

int _CompareMem(const void *a, const void *b) {
	return memcmp(a, b, _CompareMemWidth);
}

// If you want to call this as a one-off, call ComputePartLowerBound_b1_w1().  This version exists to speed
// things up when run inside a tight loop by avoiding repetitive calls to malloc()/free().
// comps[i] is the parent of seq i.
int ComputePartLowerBound_preallocated_b1_w1(struct part pt, int nTaxa, struct kruskal_component *comps, struct kruskal_edge *edges, int *uniqueSeqs, int *furthestPairSeen, int *nEdgesWithWeight, unsigned char *dataBuf) {
	int i, j;
	int nComps = 0;			// NOT the size of comps, but the number of components in it.
	int nEdges = 0;
	int edgeComp[2];
	int firstLongEdge = 0;
	int bounds[MAX_BT] = { 0 };
	int totalForestWeight = 0;
	int testEdgeComp[2];
	int temp;
	int furthestPairBest = 0;
	int nFurthestPairsSeen = 0;
	int best = 0;
	int stateCount[16] = { 0 };
	int maxEdgeWeight = 0;		// At all times, the first maxEdgeWeight+1 entries in nEdgesWithWeight[] must be correct.
	int nOrigComps;				// Unlike nComps, this doesn't change.

	gBestBoundMask = 0;			// Need to clear this up here in case we return early via a special case

	// Special-case some simple, common cases for speed
	if (!pt.nSites) {
		return 0;
	} else if (pt.nSites == 1) {
		for (i = 0; i < nTaxa; ++i) {
			++stateCount[pt.data[i * pt.memWidth]];
		}

		if (stateCount[1] + stateCount[2] + stateCount[4] + stateCount[8] == nTaxa) {
			// This single site contains no ambiguous states, so a tight bound is given by the number of
			// distinct states minus 1.
			return !!stateCount[1] + !!stateCount[2] + !!stateCount[4] + !!stateCount[8] - 1;
		}
	} else if (pt.nSites == 2) {
		// If this site pair contains no ambiguous states, we can look up the answer in a precomputed table.
		//HACK: This duplicates a lot of code in CalcSitePairIncompatCost_b1t().
		unsigned index = 0, ambig = 0;
		int base1, base2;
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
		for (i = 0; i < nTaxa; ++i) {
			base1 = map[pt.data[i * pt.memWidth]];
			base2 = map[pt.data[i * pt.memWidth + 1]];
			
			if (base1 >= 0 && base2 >= 0) {
				index |= 1 << (base1 * 4 + base2);
			} else {
				ambig = 1;
				break;
			}
		}

		if (!ambig) {
			return sitePairCost[index - 1];
		}
	}

	nEdgesWithWeight[0] = 0;

	// Build an array of unique sequence numbers as we go.
	// First sort the sequences into a separate buffer first to improve things for high taxon counts.
	//TODO: Use a bucket sort instead of qsort().
	for (i = 0; i < nTaxa; ++i) {
		memcpy(dataBuf + i * pt.nSites, pt.data + i * pt.memWidth, pt.nSites);
	}

	_CompareMemWidth = pt.nSites;
	qsort(dataBuf, nTaxa, pt.nSites, _CompareMem);

	// Helper macro to reduce typing
//#define SEQ(i) (pt.data + (i) * pt.memWidth)
#define SEQ(i) (dataBuf + (i) * pt.nSites)
	for (i = 0; i < nTaxa; ++i) {
		// Is this sequence definitely distinct?
		if (!i || memcmp(SEQ(i), SEQ(i - 1), pt.nSites)) {
			for (j = 0; j < nComps; ++j) {
				// Measure the minimum-possible distance between seqs i and uniqueSeqs[j].
				// Even though we excluded exact duplicates, this can be 0, due to seqs with ambiguous bases.
				// After thinking about this carefully, we do need to add 0-length edges to the list for Kruskal
				// to process -- imagine if there are 2 sequences having mindist 1 (e.g. "A" and "C") and then a
				// 3rd sequence having mindist 0 to both is found (e.g. "[AC]").  We need to union all of those
				// components in this case, which we might as well do in Kruskal.
				// This won't slow things down much as we don't expect there to be many 0-length edges anyway --
				// if the seqs have been sorted, they will only be produced by seqs having ambiguous bases.
				edges[nEdges].w = FitchScoreWeight1_b1_w1(SEQ(uniqueSeqs[j]), SEQ(i), pt.nSites);
				edges[nEdges].v[0] = j;			// Label vertices by unique seq numbers so we can efficiently
				edges[nEdges].v[1] = nComps;	// scan through all components later.

				// "Incrementally initialise" nEdgesWithWeight, to avoid expensive zeroing for small parts.
				if (edges[nEdges].w > maxEdgeWeight) {
					memset(nEdgesWithWeight + maxEdgeWeight + 1, 0, (edges[nEdges].w - maxEdgeWeight) * sizeof (int));
					maxEdgeWeight = edges[nEdges].w;
				}

				++nEdgesWithWeight[edges[nEdges].w];		// Used for bucket sorting later
				++nEdges;
			}

			uniqueSeqs[nComps] = i;
			comps[nComps].p = nComps;			// Also label components by unique seq numbers.
			comps[nComps].w = 0;
			++nComps;
		}
	}
#undef SEQ

	nOrigComps = nComps;

	// Bucket-sort edges by weight.
	// First convert edge counts to offsets-to-last-entry-of-that-weight-plus-one.
	for (i = 1; i <= maxEdgeWeight; ++i) {
		nEdgesWithWeight[i] += nEdgesWithWeight[i - 1];
	}

	// Backfill edges into the 2nd half of the edges array.
	//HACK: For now, just rely on edges having enough space to hold 2 copies of the edge array, and HACKfully advance
	// the pointer to point to the 2nd copy.  This would break things if we still deallocated edges at the end of
	// this function.
	for (i = 0; i < nEdges; ++i) {
		edges[nEdges + --nEdgesWithWeight[edges[i].w]] = edges[i];
	}
	edges += nEdges;		//HACK: see comment above.

	// Run the Kruskal algorithm to build an MST.
	// We can stop running when the graph is connected.
	//for (i = 0; i < nEdges; ++i) {
	for (i = 0; nComps > 1; ++i) {
		// Would this edge create a cycle?
		edgeComp[0] = find_component(comps, edges[i].v[0]);
		edgeComp[1] = find_component(comps, edges[i].v[1]);

		if (edgeComp[0] != edgeComp[1]) {
			// No: we can combine these components.

			// First check if we are 1-connected so far
			if (edges[i].w > 1 && !firstLongEdge) {
				firstLongEdge = edges[i].w;
				bounds[BT_ONECONN] = totalForestWeight + firstLongEdge;

				// The old BT_FURTHESTPAIR bound is equivalent to BT_ONECONN if exactly 2 components, so don't bother testing for it.
				// Also, by the triangle inequality, BT_FURTHESTPAIR must be dominated by BT_3COMPS on 3 components.

				// BT_3COMPS:
				// If there are 3 or more 1-conn components, they must either be connected by the 2 shortest of the 3
				// shortest intercomponent edges, or by a Steiner point.  Choose the best set of 3 components.
				//TODO: Index the components so that they have no "holes" -- then we can run them through much faster
				// during the cubic step.
				if (nComps >= 3) {
					int v[3];

					memset(furthestPairSeen, 0, nTaxa * nTaxa * sizeof (int));
					for (j = i; nFurthestPairsSeen < nComps * (nComps - 1) / 2; ++j) {
						assert(j < nEdges);
						testEdgeComp[0] = find_component(comps, edges[j].v[0]);
						testEdgeComp[1] = find_component(comps, edges[j].v[1]);

						if (testEdgeComp[0] != testEdgeComp[1]) {
							// A separate pair of components: measure them.
							if (testEdgeComp[0] > testEdgeComp[1]) {
								temp = testEdgeComp[0];
								testEdgeComp[0] = testEdgeComp[1];
								testEdgeComp[1] = temp;
							}

							if (!furthestPairSeen[testEdgeComp[0] * nTaxa + testEdgeComp[1]]) {
								furthestPairSeen[testEdgeComp[0] * nTaxa + testEdgeComp[1]] = edges[j].w;
								++nFurthestPairsSeen;
							}
						}
					}

					// Consider every triple of components.
					for (v[0] = 0; v[0] < nOrigComps - 2; ++v[0]) {
						if (comps[v[0]].p != v[0]) {
							continue;
						}

						for (v[1] = v[0] + 1; v[1] < nOrigComps - 1; ++v[1]) {
							if (comps[v[1]].p != v[1]) {
								continue;
							}

							for (v[2] = v[1] + 1; v[2] < nOrigComps; ++v[2]) {
								int best3CompNoSteiner;
								int best3CompSteiner;

								if (comps[v[2]].p != v[2]) {
									continue;
								}

								// Find the shortest way to connect these 3 components directly (i.e. shortest edge pair)
								if (furthestPairSeen[v[0] * nTaxa + v[1]] < furthestPairSeen[v[0] * nTaxa + v[2]]) {
									if (furthestPairSeen[v[0] * nTaxa + v[2]] < furthestPairSeen[v[1] * nTaxa + v[2]]) {
										best3CompNoSteiner = furthestPairSeen[v[0] * nTaxa + v[1]] + furthestPairSeen[v[0] * nTaxa + v[2]];
									} else {
										best3CompNoSteiner = furthestPairSeen[v[0] * nTaxa + v[1]] + furthestPairSeen[v[1] * nTaxa + v[2]];
									}
								} else if (furthestPairSeen[v[0] * nTaxa + v[1]] < furthestPairSeen[v[1] * nTaxa + v[2]]) {
									best3CompNoSteiner = furthestPairSeen[v[0] * nTaxa + v[1]] + furthestPairSeen[v[0] * nTaxa + v[2]];
								} else {
									best3CompNoSteiner = furthestPairSeen[v[1] * nTaxa + v[2]] + furthestPairSeen[v[0] * nTaxa + v[2]];
								}

								// Regardless of whether we use a Steiner point, we can always add the weight
								// of the components.
								best3CompNoSteiner += comps[v[0]].w + comps[v[1]].w + comps[v[2]].w;
								
								// What if we use a Steiner point in the middle?
								best3CompSteiner = (furthestPairSeen[v[0] * nTaxa + v[1]] +
									furthestPairSeen[v[1] * nTaxa + v[2]] +
									furthestPairSeen[v[0] * nTaxa + v[2]] + 1) / 2 +
									comps[v[0]].w + comps[v[1]].w + comps[v[2]].w;

								if (best3CompNoSteiner <= best3CompSteiner && best3CompNoSteiner > bounds[BT_3COMPS]) {
									// An optimal solution can be reached without using a Steiner point.
									bounds[BT_3COMPS] = best3CompNoSteiner;
								}

								if (best3CompSteiner <= best3CompNoSteiner && best3CompSteiner > bounds[BT_3COMPSSTEINER]) {
									// An optimal solution can be reached using a Steiner point.
									bounds[BT_3COMPSSTEINER] = best3CompSteiner;
								}
							}
						}
					}
				}
			}

			union_components(comps, edgeComp[0], edgeComp[1], edges[i].w);
			--nComps;
			totalForestWeight += edges[i].w;
		}
	}

	if (!firstLongEdge) {
		// The entire taxon set is 1-connected.
		bounds[BT_ONECONN] = totalForestWeight;
	} else {
		// Can use half the MST as a bound.  If nSites <= 2, David Bryant proved in "A Subdivision approach to
		// maximum parsimony" (Annals of Combinatorics, in press as of 2009) that we can actually use the entire MST.
		bounds[BT_MST] = totalForestWeight;
		if (pt.nSites > 2) {
			bounds[BT_MST] /= 2;
		}
	}

	// Find the maximum bound
	for (i = 0; i < MAX_BT; ++i) {
		if (bounds[i] > best) {
			best = bounds[i];
			gBestBoundMask = 0;
		}

		// Interesting to know which bounds do the best.
		if (bounds[i] == best) {
			gBestBoundMask |= 1 << i;
		}
	}

	return best;
}

// For convenience.  Currently, all existing code calls ComputePartLowerBound_preallocated_b1_w1() directly instead.
int ComputePartLowerBound_b1_w1(struct part pt, int nTaxa) {
	int result;
	struct compute_part_lower_bound_buffers buffers;

	SetupBuffersForComputePartLowerBound(&buffers, nTaxa, pt.nSites);
	result = ComputePartLowerBound_preallocated_b1_w1(pt, nTaxa, buffers.comps, buffers.edges, buffers.uniqueSeqs, buffers.furthestPairSeen, buffers.nEdgesWithWeight, buffers.dataBuf);
	TearDownBuffersForComputePartLowerBound(&buffers);

	return result;
}

// The only safe upper bound on the weight of a part is given by summing the UBs for each column.
// If we knew that the entire problem was this part then better UBs could be found, but they are invalid
// if the part is part of a bigger partition since other sites can force more sites in this part to break.
int ComputePartUpperBound_b1_w1(struct part pt, int nTaxa) {
	int i, score = 0;
	for (i = 0; i < pt.nSites; ++i) {
		score += SiteWeightUpperBound_b1(pt.data, i, nTaxa, pt.memWidth);
	}

	return score;
}

// Computes the minimum score that must be added by the remaining (nTaxa - nTaxaOnTree) taxa.
// We use the number of distinct bases in each site that have not yet been seen.  This can result
// in bounds better than LB(full) - UB(partial), because we do not make the costly assumption
// that the UB is always reached on the partial tree.  E.g. consider the site CCAAAG: the UB on
// the 1st 5 taxa is 2, since 2 separate mutations *may* be needed for the Cs on some tree; the
// LB on the entire dataset is also 2, so LB - UB = 0.  But the DISTBASESBOUND tells us that
// **1 more** mutation will always be required for the final taxon, regardless of whether 1
// or 2 mutations are needed for the 1st 5.
// Incorporating this into the PARTBOUND means that (unlike previously!) it will always dominate
// DISTBASESBOUND, since in the worst case we will choose either all-sites-additive (equivalent
// to DISTBASESBOUND) or all-sites-nonadditive (equivalent to the earlier PARTBOUND).
//HACK: This function reimplements functionality similar to that already in InitDistBasesPerSiteBoundRest_b1()
// and SiteWeightLowerBound_b1t().  But it's not identical, and I'm too lazy to figure out exactly
// how to best refactor things.
int ComputePartAdditiveScore_b1_w1(struct part pt, int nTaxa, int nTaxaOnTree) {
	int i, j, nDistinct = 0;

	// This approach doesn't maximise the score we get from ambiguous sites -- e.g. if a transposed
	// site looks like A[AC]C then for nTaxaOnTree == 1 we would return just 0, even though clearly
	// we could return 1 since the 2nd taxon cannot be both A and C at the same time.
	for (i = 0; i < pt.nSites; ++i) {
		char seen = 0;			// A bitmask containing all bases already seen at the current site
		for (j = 0; j < nTaxa; ++j) {
			if (j >= nTaxaOnTree && pt.data[pt.memWidth * j + i] & ~seen) {
				++nDistinct;
			}

			seen |= pt.data[pt.memWidth * j + i];
		}
	}

	return nDistinct;
}

// Resizes and reshapes a part so that it has memWidth **at least** minWidth (i.e. it may be more, for efficiency reasons).
// This may shrink or grow the partition as necessary.
// If pt->memWidth == 0, this is treated as a first-time allocation.
// If minWidth == 0, this is treated as a deallocation.
void PartResize(struct part *pt, int nTaxa, int minWidth) {
	int i;
	minWidth = (!minWidth - 1) & ((minWidth - 1) / PARTBLOCKSIZE + 1) * PARTBLOCKSIZE;		// Round up to a multiple of PARTBLOCKSIZE bytes.
	if (pt->memWidth) {
		// Possibly resizing (and reshaping) an existing part.
		if (minWidth > pt->memWidth) {
			// Growing.  Can upsize in-place after realloc().
			pt->data = realloc(pt->data, nTaxa * minWidth);
			for (i = nTaxa - 1; i > 0; --i) {
				// Need memmove() not memcpy() since in general the regions will overlap.
				// In theory we only need to copy pt->nSites bytes, but copying full blocks may actually be faster.
				memmove(pt->data + i * minWidth, pt->data + i * pt->memWidth, pt->memWidth);
			}
			pt->memWidth = minWidth;
		} else if (minWidth < pt->memWidth) {
			if (minWidth > PARTBLOCKSIZE && minWidth + PARTBLOCKSIZE < pt->memWidth) {
				// Shrinking.  Use hysteresis to avoid many resizings across a resize boundary.
				// Can downsize in-place before realloc().
				for (i = 1; i < nTaxa; ++i) {
					// Need memmove() not memcpy() since in general the regions will overlap.
					// In theory we only need to copy pt->nSites bytes, but copying full blocks may actually be faster.
					memmove(pt->data + i * minWidth, pt->data + i * pt->memWidth, pt->memWidth);
				}
				pt->data = realloc(pt->data, nTaxa * minWidth);
				pt->memWidth = minWidth;
			} else if (!minWidth) {
				// Interpret this as a deallocation request.
				free(pt->data);
				pt->data = NULL;
				pt->memWidth = 0;
			}
		}
	} else {
		// Interpret this as a request to allocate space for a new part.
		pt->memWidth = minWidth;		// Round up to a multiple of PARTBLOCKSIZE bytes.
		pt->data = malloc(nTaxa * pt->memWidth);
	}
}

// Copy a rectangular block of sites from one part to another.  Assumes the parts have adequate memory allocated.
// Assumes the blocks do not overlap.  (If they do, you'll want to rewrite this function with memmove() instead of memcpy().)
void PartCopySites(struct part dst, int dSite, struct part src, int sSite, int nTaxa, int nSites) {
	int i;
	assert(sSite + nSites <= src.memWidth);
	assert(dSite + nSites <= dst.memWidth);
	for (i = 0; i < nTaxa; ++i) {
		memcpy(dst.data + i * dst.memWidth + dSite, src.data + i * src.memWidth + sSite, nSites);
	}
}

// Swap a single pair of sites between two parts.
void PartSwapSites(struct part dst, int dSite, struct part src, int sSite, int nTaxa) {
	int i;
	unsigned char temp;

	if (src.data == dst.data && sSite == dSite) {
		return;
	}

	for (i = 0; i < nTaxa; ++i) {
		temp = dst.data[i * dst.memWidth + dSite];
		dst.data[i * dst.memWidth + dSite] = src.data[i * src.memWidth + sSite];
		src.data[i * src.memWidth + sSite] = temp;
	}
}

void PartMoveSiteAndResize(struct scored_part *parts, int *move, int nTaxa) {
	++parts[move[1]].pt.nSites;
	PartResize(&parts[move[1]].pt, nTaxa, parts[move[1]].pt.nSites + PARTBOUND_MAX_SITES_MOVED);

	// Append the new site to the destination partition
	PartCopySites(parts[move[1]].pt, parts[move[1]].pt.nSites - 1, parts[move[0]].pt, move[2], nTaxa, 1);

	// Remove the original site from the source partition by swapping with the last position if necessary
	if (move[2] != parts[move[0]].pt.nSites - 1) {
		PartSwapSites(parts[move[0]].pt, move[2], parts[move[0]].pt, parts[move[0]].pt.nSites - 1, nTaxa);
	}

	if (--parts[move[0]].pt.nSites) {
		PartResize(&parts[move[0]].pt, nTaxa, parts[move[0]].pt.nSites + PARTBOUND_MAX_SITES_MOVED);
	} else {
		// This part has shrunk to empty: deallocate memory for it.
		PartResize(&parts[move[0]].pt, nTaxa, 0);
	}

	// Record the fact that we have changed these parts
	parts[move[0]].touched = 1;
	parts[move[1]].touched = 1;
#ifdef PARTBOUND_VERBOSE_LOGGING
	DBGPRINT4("Moved P%d[%d] -> P%d.\n", move[0], move[2], move[1]);
	DBGPRINT7("P%d now has %d sites (mW=%d); P%d now has %d sites (mW=%d).\n", move[0], parts[move[0]].pt.nSites, parts[move[0]].pt.memWidth, move[1], parts[move[1]].pt.nSites, parts[move[1]].pt.memWidth);
#endif	// PARTBOUND_VERBOSE_LOGGING
}

// Attempt to move sites around between parts so as to maximise the final score of the entire partition.  This
// is found by taking, for each part, the maximum of the clipped score (LB - UB, clipped at 0) and the additive
// score found by counting the number of distinct bases in taxa not yet present on the tree, and summing these
// across all parts in the partition.
// If greediness == 3, immediately perform any move that improves the clipped or raw score.  (TODO)
// If greediness == 2, immediately perform the best move found between any pair of parts.
// If greediness == 1, gather up all best moves for each part pair and then greedily choose a "maximal matching" of
// moves to make on distinct part pairs (since these will be independent).  (TODO)
// Return value: the index of the 1st touched part.  All touched parts will be at or to the right of this position.
// Iff the return value == (the new value of) ptn->nParts, then no parts were touched -- i.e. no moves could be performed.
int OptimisePartitionOnePass(struct scored_partition *ptn, int nTaxa, int nTaxaOnTree, int greediness, enum partbound_strategy_type strategy, int iFirstTouchedPart, int maxDist) {
	int i, j;
	int iPart[2];
	int scoreLB[2];
	int scoreUB[2];
	int scoreAdditive[2];
	unsigned maxBoundTypes[2];
	int maxRaw;
	int final;
	int bestFinalScoreImprovement;
	int bestFinalClippedScoreLB[2];
	int bestFinalClippedScoreUB[2];
	int bestFinalAdditiveScore[2];
	unsigned bestFinalClippedMaxBoundTypes[2];
	int bestFinalMove[3];			// {i, j, k}: Move site k from part i to part j.
	int bestRawScore;
	int bestRawScoreImprovement;
	int bestRawScoreLB[2];
	int bestRawScoreUB[2];
	int bestRawAdditiveScore[2];
	unsigned bestRawMaxBoundTypes[2];
	int bestRawMove[3];				// {i, j, k}: Move site k from part i to part j.
	int k;
	int d;
	int moved;
	int nMoves = 0;
	struct compute_part_lower_bound_buffers buffers;
	int siteUB;
#ifdef PARTBOUND_LOG_EVERY_MOVE_CONSIDERED
	char possibleMoveLogFName[80];
	FILE *possibleMoveLogFile;

	sprintf(possibleMoveLogFName, "OptimisePartitionOnePass.%dparts.%dtaxa.%dtaxaOnTree.dump.txt", ptn->nParts, nTaxa, nTaxaOnTree);
	possibleMoveLogFile = fopen(possibleMoveLogFName, "w");			//DEBUG
#endif	// PARTBOUND_LOG_EVERY_MOVE_CONSIDERED

	DBGPRINT5("   OptimisePartitionOnePass(nParts=%d, nTaxa=%d, nTaxaOnTree=%d, iFirstTouchedPart=%d) called.\n", ptn->nParts, nTaxa, nTaxaOnTree, iFirstTouchedPart);

	assert(greediness == 2);		// Currently the only implemented option

	SetupBuffersForComputePartLowerBound(&buffers, nTaxa, maxDist);

	bestFinalScoreImprovement = 0;
	bestRawScore = INT_MIN;			// This will be negative in general, because we exclude moves that improve the clipped score.
	bestRawScoreImprovement = 0;
	iPart[0] = iFirstTouchedPart;					// Whoops!  Need to remember to initialise this!
	for (i = iFirstTouchedPart; i < ptn->nParts; ++i) {
		if (!ptn->parts[i].pt.nSites) {
			continue;
		}

		// Slide the next part forward if necessary
		if (iPart[0] != i) {
			ptn->parts[iPart[0]] = ptn->parts[i];
		}

		// i: The location of the next part to read from.  Assume: ptn->parts[i].pt.nSites > 0.
		// iPart[0]: The location of the next part, copied from ptn->parts[i] if they are different.
		// Only advanced when the result of an outer loop iteration leaves a nonempty part.

		//DBGPRINT3("   OptimisePartitionOnePass(): Considering moves to/from part %d (which now occupies position %d).\n", i, iPart[0]);
		ptn->parts[iPart[0]].touched = 0;			// This will be set to 1 if this part is changed in any way.
		for (iPart[1] = 0; iPart[1] < iPart[0]; ++iPart[1]) {
#ifdef PARTBOUND_VERBOSE_LOGGING
			DBGPRINT8("   ===== Considering parts P%d (originally P%d; nS=%d, mW=%d) and P%d (nS=%d, mW=%d) =====\n", iPart[0], i, ptn->parts[iPart[0]].pt.nSites, ptn->parts[iPart[0]].pt.memWidth, iPart[1], ptn->parts[iPart[1]].pt.nSites, ptn->parts[iPart[1]].pt.memWidth);
#endif	// PARTBOUND_VERBOSE_LOGGING
			// Try moves in both directions, unless both parts contain exactly one site.
			for (d = 0; d < 2 - (ptn->parts[iPart[0]].pt.nSites == 1 && ptn->parts[iPart[1]].pt.nSites == 1); ++d) {
				// Old raw, clipped and final scores can be calculated outside of the site loop.
				int oldMaxRaw = iMax(
					ptn->parts[iPart[d]].scoreLB - ptn->parts[iPart[d]].scoreUB,
					ptn->parts[iPart[1-d]].scoreLB - ptn->parts[iPart[1-d]].scoreUB
				);

				int oldFinal =
					selectFinalScore(strategy, ptn->parts[iPart[d]].scoreLB,   ptn->parts[iPart[d]].scoreUB,   ptn->parts[iPart[d]].scoreAdditive) +
					selectFinalScore(strategy, ptn->parts[iPart[1-d]].scoreLB, ptn->parts[iPart[1-d]].scoreUB, ptn->parts[iPart[1-d]].scoreAdditive);

				// Try all 1-site moves from part iPart[d] to part iPart[1-d].
				// We assume here that both parts have room for an additional site.
				--ptn->parts[iPart[d]].pt.nSites;
				++ptn->parts[iPart[1-d]].pt.nSites;

				for (k = 0; k < ptn->parts[iPart[d]].pt.nSites + 1; ++k) {		// Remember ptn->parts[iPart[d]].nSites has already been decremented
					// Append the site to part iPart[1-d].
					PartCopySites(ptn->parts[iPart[1-d]].pt, ptn->parts[iPart[1-d]].pt.nSites - 1, ptn->parts[iPart[d]].pt, k, nTaxa, 1);	// Remember ptn->parts[iPart[1-d]].nSites has already been incremented

					// To "remove" site k from part iPart[d], swap it with the last site.
					PartSwapSites(ptn->parts[iPart[d]].pt, k, ptn->parts[iPart[d]].pt, ptn->parts[iPart[d]].pt.nSites, nTaxa);		// Remember ptn->parts[iPart[d]].nSites has already been decremented

					// Compute lower bounds (hard work)
					scoreLB[0] = ComputePartLowerBound_preallocated_b1_w1(ptn->parts[iPart[d]].pt, nTaxa, buffers.comps, buffers.edges, buffers.uniqueSeqs, buffers.furthestPairSeen, buffers.nEdgesWithWeight, buffers.dataBuf);
					maxBoundTypes[0] = gBestBoundMask;
					scoreLB[1] = ComputePartLowerBound_preallocated_b1_w1(ptn->parts[iPart[1-d]].pt, nTaxa, buffers.comps, buffers.edges, buffers.uniqueSeqs, buffers.furthestPairSeen, buffers.nEdgesWithWeight, buffers.dataBuf);
					maxBoundTypes[1] = gBestBoundMask;

					// Compute upper bounds (easy, since a UB is always exactly the sum of per-site UBs)
					siteUB = SiteWeightUpperBound_b1(ptn->parts[iPart[d]].pt.data, ptn->parts[iPart[d]].pt.nSites, nTaxaOnTree, ptn->parts[iPart[d]].pt.memWidth);
					scoreUB[0] = ptn->parts[iPart[d]].scoreUB - siteUB;
					scoreUB[1] = ptn->parts[iPart[1-d]].scoreUB + siteUB;

					// Compute additive bounds using the distinct-bases-remaining approach
					scoreAdditive[0] = ComputePartAdditiveScore_b1_w1(ptn->parts[iPart[d]].pt, nTaxa, nTaxaOnTree);
					scoreAdditive[1] = ComputePartAdditiveScore_b1_w1(ptn->parts[iPart[1-d]].pt, nTaxa, nTaxaOnTree);

					// Raw scores can be negative
					// To avoid moving every site for which UB(site) > LB(site) into one enormous part, which will
					// slow down later passes enormously, don't consider the raw score of empty parts.
					maxRaw = iMax(ptn->parts[iPart[d]].pt.nSites ? scoreLB[0] - scoreUB[0] : INT_MIN, scoreLB[1] - scoreUB[1]);

					// The final score defaults to the maximum of the sum of the clipped scores and
					// the sum of the additive scores, though it can also be just the clipped score,
					// just the additive score (equivalent to DISTBASESBOUND) or just the lower bound
					// (in which case it can be negative).
					final =
						selectFinalScore(strategy, scoreLB[0], scoreUB[0], scoreAdditive[0]) +
						selectFinalScore(strategy, scoreLB[1], scoreUB[1], scoreAdditive[1]);

#ifdef PARTBOUND_LOG_EVERY_MOVE_CONSIDERED
					fprintf(possibleMoveLogFile, "Considering move: P%d[%d] -> P%d.\n", iPart[d], k, iPart[1-d]);
					fprintf(possibleMoveLogFile, "scoreLB[0]=%d\n", scoreLB[0]);
					fprintf(possibleMoveLogFile, "scoreLB[1]=%d\n", scoreLB[1]);
					fprintf(possibleMoveLogFile, "scoreUB[0]=%d\n", scoreUB[0]);
					fprintf(possibleMoveLogFile, "scoreUB[1]=%d\n", scoreUB[1]);
					fprintf(possibleMoveLogFile, "raw[0]=%d\n", raw[0]);
					fprintf(possibleMoveLogFile, "raw[1]=%d\n", raw[1]);
					fprintf(possibleMoveLogFile, "clipped[0]=%d\n", clipped[0]);
					fprintf(possibleMoveLogFile, "clipped[1]=%d\n", clipped[1]);
					fprintf(possibleMoveLogFile, "scoreAdditive[0]=%d\n", scoreAdditive[0]);
					fprintf(possibleMoveLogFile, "scoreAdditive[1]=%d\n", scoreAdditive[1]);
					fprintf(possibleMoveLogFile, "final=%d\n", final);
#endif	// PARTBOUND_LOG_EVERY_MOVE_CONSIDERED

					// Is this the best move so far?
					if (final - oldFinal > bestFinalScoreImprovement) {
						bestFinalScoreImprovement = final - oldFinal;
						bestFinalClippedScoreLB[0] = scoreLB[0];
						bestFinalClippedScoreUB[0] = scoreUB[0];
						bestFinalClippedScoreLB[1] = scoreLB[1];
						bestFinalClippedScoreUB[1] = scoreUB[1];
						bestFinalAdditiveScore[0] = scoreAdditive[0];
						bestFinalAdditiveScore[1] = scoreAdditive[1];
						bestFinalMove[0] = iPart[d];
						bestFinalMove[1] = iPart[1-d];
						bestFinalMove[2] = k;
						bestFinalClippedMaxBoundTypes[0] = maxBoundTypes[0];
						bestFinalClippedMaxBoundTypes[1] = maxBoundTypes[1];
#ifdef PARTBOUND_VERBOSE_LOGGING
						DBGPRINT7("   Found better clipped/additive move (move site %d from part %d to part %d) with final score %d (clipped score %d, additive score %d)!\n", k, iPart[d], iPart[1-d], bestFinalScore, clipped[0] + clipped[1], scoreAdditive[0] + scoreAdditive[1]);
#endif	// PARTBOUND_VERBOSE_LOGGING
					} else if (final >= oldFinal && maxRaw <= 0 && maxRaw > oldMaxRaw && maxRaw > bestRawScore) {
						// Note that we must check final >= oldFinal above now that we allow strategies that
						// could produce a -ve final score (PBS_LB).
						bestRawScore = maxRaw;
						bestRawScoreImprovement = maxRaw - oldMaxRaw;
						bestRawScoreLB[0] = scoreLB[0];
						bestRawScoreUB[0] = scoreUB[0];
						bestRawScoreLB[1] = scoreLB[1];
						bestRawScoreUB[1] = scoreUB[1];
						bestRawAdditiveScore[0] = scoreAdditive[0];
						bestRawAdditiveScore[1] = scoreAdditive[1];
						bestRawMove[0] = iPart[d];
						bestRawMove[1] = iPart[1-d];
						bestRawMove[2] = k;
						bestRawMaxBoundTypes[0] = maxBoundTypes[0];
						bestRawMaxBoundTypes[1] = maxBoundTypes[1];
#ifdef PARTBOUND_VERBOSE_LOGGING
						DBGPRINT6("   Found better raw move (move site %d from part %d to part %d) with raw score %d!  (Old score=%d)\n", k, iPart[d], iPart[1-d], bestRawScore, oldRaw[0]);
#endif	// PARTBOUND_VERBOSE_LOGGING
					}

					// "Unremove" the site in part i.
					PartSwapSites(ptn->parts[iPart[d]].pt, k, ptn->parts[iPart[d]].pt, ptn->parts[iPart[d]].pt.nSites, nTaxa);		// Remember ptn->parts[iPart[d]].nSites has already been decremented
				}

				++ptn->parts[iPart[d]].pt.nSites;
				--ptn->parts[iPart[1-d]].pt.nSites;
			}

			moved = 0;
			if (greediness == 2) {
				// Grab the best move for this part pair right away, if it improves anything.
				if (bestFinalScoreImprovement > 0) {
					PartMoveSiteAndResize(ptn->parts, bestFinalMove, nTaxa);

					// Update totals
					ptn->totalFinalScore += bestFinalScoreImprovement;
					ptn->scoreLB +=
						bestFinalClippedScoreLB[0] - ptn->parts[bestFinalMove[0]].scoreLB +
						bestFinalClippedScoreLB[1] - ptn->parts[bestFinalMove[1]].scoreLB;

					// Save scores
					ptn->parts[bestFinalMove[0]].scoreLB = bestFinalClippedScoreLB[0];
					ptn->parts[bestFinalMove[0]].scoreUB = bestFinalClippedScoreUB[0];
					ptn->parts[bestFinalMove[1]].scoreLB = bestFinalClippedScoreLB[1];
					ptn->parts[bestFinalMove[1]].scoreUB = bestFinalClippedScoreUB[1];
					ptn->parts[bestFinalMove[0]].scoreAdditive = bestFinalAdditiveScore[0];
					ptn->parts[bestFinalMove[1]].scoreAdditive = bestFinalAdditiveScore[1];
					ptn->parts[bestFinalMove[0]].bestBoundsMask = bestFinalClippedMaxBoundTypes[0];
					ptn->parts[bestFinalMove[1]].bestBoundsMask = bestFinalClippedMaxBoundTypes[1];
					
					moved = 1;
				} else if (bestRawScoreImprovement > 0) {
					PartMoveSiteAndResize(ptn->parts, bestRawMove, nTaxa);

					// Update totals (just the LB for raw moves; there's no sum that can be sensibly maintained)
					ptn->scoreLB +=
						bestRawScoreLB[0] - ptn->parts[bestRawMove[0]].scoreLB +
						bestRawScoreLB[1] - ptn->parts[bestRawMove[1]].scoreLB;

					// Save scores
					ptn->parts[bestRawMove[0]].scoreLB = bestRawScoreLB[0];
					ptn->parts[bestRawMove[0]].scoreUB = bestRawScoreUB[0];
					ptn->parts[bestRawMove[1]].scoreLB = bestRawScoreLB[1];
					ptn->parts[bestRawMove[1]].scoreUB = bestRawScoreUB[1];
					ptn->parts[bestRawMove[0]].scoreAdditive = bestRawAdditiveScore[0];
					ptn->parts[bestRawMove[1]].scoreAdditive = bestRawAdditiveScore[1];
					ptn->parts[bestRawMove[0]].bestBoundsMask = bestRawMaxBoundTypes[0];
					ptn->parts[bestRawMove[1]].bestBoundsMask = bestRawMaxBoundTypes[1];

					moved = 1;
				}
			}

			if (moved) {
				// Reset "best" measures
				bestFinalScoreImprovement = 0;
				bestRawScore = INT_MIN;			// This will be negative in general, because we exclude moves that improve the clipped score.
				bestRawScoreImprovement = 0;
				++nMoves;

				// Make sure that if a part was removed entirely by the move, it is part i.
				// This is for efficiency only: parts to the left of part i will be rescanned
				// many times by the iPart[1] inner loop, so having holes there is more expensive.
				if (ptn->parts[iPart[0]].pt.nSites && !ptn->parts[iPart[1]].pt.nSites) {
					struct scored_part tempPart = ptn->parts[iPart[1]];
					ptn->parts[iPart[1]] = ptn->parts[iPart[0]];
					ptn->parts[iPart[0]] = tempPart;
				}

				// Has part[iPart[0]] been removed entirely?  If so, decrement the loop index
				// to enable it to be backfilled on the next pass, and skip to the next outer loop pass.
				if (!ptn->parts[iPart[0]].pt.nSites) {
#ifdef PARTBOUND_VERBOSE_LOGGING
					DBGPRINT2("   Part P%d is empty and will be backfilled.\n", iPart[0]);
#endif	// PARTBOUND_VERBOSE_LOGGING
					--iPart[0];
					break;
				}
			}
		}

		++iPart[0];
	}

	TearDownBuffersForComputePartLowerBound(&buffers);

#ifdef PARTBOUND_LOG_EVERY_MOVE_CONSIDERED
	fclose(possibleMoveLogFile);		//DEBUG
#endif	// PARTBOUND_LOG_EVERY_MOVE_CONSIDERED

	// We may have deleted parts along the way.
	ptn->nParts = iPart[0];
	DBGPRINT3("   Performed %d moves in this pass.  %d parts remain.\n", nMoves, ptn->nParts);

	// If any moves were performed, we need to move the changed parts to the far end so that we can
	// efficiently process every pair of parts containing >= 1 changed part on the next pass.
	if (nMoves) {
		struct scored_part tempPart;

		// I've thought about the loop logic carefully.  We never change i twice, or j twice, or i and j once each
		// without a test in between, so != and == can be used.  They make it clear that the
		// only way the loop exits is when i == j, which simplifies things a bit.
		// Invariants true at all times:
		//    - Everything in positions > j is touched.
		//    - Everything in positions < i is untouched.
		i = 0;								// i: Where we will next look for a touched part
		j = ptn->nParts - 1;				// j: Where we will put the next touched part
		while (1) {
			// Skip over parts that we don't need to move
			while (!ptn->parts[i].touched && i != j) {
				++i;
			}

			if (i == j) {
				break;
			}

			// Skip over parts already in the right place
			while (ptn->parts[j].touched && i != j) {
				--j;
			}

			if (i == j) {
				break;
			}

			// Swap the touched part over to the right
			tempPart = ptn->parts[i];
			ptn->parts[i] = ptn->parts[j];
			ptn->parts[j] = tempPart;
		}

		// At this point, i (== j) is either the last untouched part, or the 1st touched part.
		if (!ptn->parts[i].touched) {
			++i;
		}

		DBGPRINT3("   After placing touched parts at the end, the first touched part is P%d (i.e. %d parts were touched).\n", i, ptn->nParts - i);
		return i;
	} else {
		// Return "one past the end" to indicate no moves occurred, so there should be no more passes.
		return ptn->nParts;
	}
}

// Repeatedly call OptimisePartitionOnePass() until no more positive changes are possible.
// At completion, ptn->totalFinalScore will be the greatest value (that we can find) that can safely
// be added to any tree with the 1st nTaxaOnTree taxa on it and still produce a lower bound on
// the MP tree score on the full dataset (assuming strategy != PBS_LB).
void OptimisePartition(struct scored_partition *ptn, int nTaxa, int nTaxaOnTree, int greediness, enum partbound_strategy_type strategy, int maxItersNoImprovement, int maxDist) {
	int i;
	int iFirstTouchedPart = 0;
	int nItersNoImprovement = 0;

	for (i = 0; iFirstTouchedPart < ptn->nParts; ++i) {
		int origTotalFinalScore = ptn->totalFinalScore;

		DBGPRINT4("  OptimisePartition(): Beginning pass #%d from part #%d of %d...\n", i + 1, iFirstTouchedPart, ptn->nParts);
		iFirstTouchedPart = OptimisePartitionOnePass(ptn, nTaxa, nTaxaOnTree, greediness, strategy, iFirstTouchedPart, maxDist);
		DBGPRINT6("  OptimisePartition(): After %d passes: LB=%d, UB=%d, clipped score = %d, final score = %d\n", i + 1, ptn->scoreLB, ptn->scoreUB, getTotalClippedScore(ptn), ptn->totalFinalScore);

		if (ptn->totalFinalScore == origTotalFinalScore) {
			if (maxItersNoImprovement && ++nItersNoImprovement == maxItersNoImprovement) {
				DBGPRINT2("  OptimisePartition(): Reached %d passes with no improvement in the final score: stopping.\n", maxItersNoImprovement);
				break;
			}
		} else {
			nItersNoImprovement = 0;
		}
	}

	DBGPRINT1("  OptimisePartition(): No more good moves.\n");
}

void PartitionDump(struct scored_partition *ptn, int nTaxa, FILE *f) {
	int i, j, k;

	for (i = 0; i < nTaxa; ++i) {
		for (j = 0; j < ptn->nParts; ++j) {
			if (j) {
				fputc('|', f);
			}

			for (k = 0; k < ptn->parts[j].pt.nSites; ++k) {
				fputc(GetBaseFromMask(ptn->parts[j].pt.data[i * ptn->parts[j].pt.memWidth + k]), f);
			}
		}

		fputc('\n', f);
	}
}

void PartitionDumpStats(struct scored_partition *ptn, enum partbound_strategy_type strategy, FILE *f) {
	int i;
	int totalRawScore = 0;
	int totalClippedScore = 0;
	int totalAdditiveScore = 0;
	int totalFinalScore = 0;
	int totalLB = 0;
	int totalUB = 0;
	int narrowestPart = 0;
	int widestPart = 0;
	int totalWidth = 0;

	fprintf(f, "  part#  nSites memWdth scoreLB scoreUB     raw clipped additve   final  bestBTs\n");
	fprintf(f, "-------+-------+-------+-------+-------+-------+-------+-------+-------+-------+\n");
	for (i = 0; i < ptn->nParts; ++i) {
		int rawScore = ptn->parts[i].scoreLB - ptn->parts[i].scoreUB;
		int clippedScore = rawScore > 0 ? rawScore : 0;
		int finalScore = selectFinalScore(strategy, ptn->parts[i].scoreLB, ptn->parts[i].scoreUB, ptn->parts[i].scoreAdditive);
		fprintf(f, "%7d %7d %7d %7d %7d %7d %7d %7d %7d   0x%04x\n",
			i,
			ptn->parts[i].pt.nSites,
			ptn->parts[i].pt.memWidth,
			ptn->parts[i].scoreLB,
			ptn->parts[i].scoreUB,
			rawScore,
			clippedScore,
			ptn->parts[i].scoreAdditive,
			finalScore,
			ptn->parts[i].bestBoundsMask
		);

		totalWidth += ptn->parts[i].pt.nSites;
		totalRawScore += rawScore;
		totalClippedScore += clippedScore;
		totalAdditiveScore += ptn->parts[i].scoreAdditive;
		totalFinalScore += finalScore;
		totalLB += ptn->parts[i].scoreLB;
		totalUB += ptn->parts[i].scoreUB;
		if (ptn->parts[i].pt.nSites < ptn->parts[narrowestPart].pt.nSites) {
			narrowestPart = i;
		}
		if (ptn->parts[i].pt.nSites > ptn->parts[widestPart].pt.nSites) {
			widestPart = i;
		}
	}

	fprintf(f, "Number of parts:      %6d\n", ptn->nParts);
	fprintf(f, "Total width:          %6d\n", totalWidth);
	fprintf(f, "Total raw score:      %6d\n", totalRawScore);
	// We don't score the total raw score in the scored_partition object, so we can't verify this.
	fprintf(f, "Total clipped score:  %6d\n", totalClippedScore);
	fprintf(f, "Total additive score: %6d\n", totalAdditiveScore);
	fprintf(f, "Total final score:    %6d\n", totalFinalScore);
	if (totalFinalScore != ptn->totalFinalScore) {
		fprintf(f, "ZOINKS!  The total final score stored inside the partition is %d!\n", ptn->totalFinalScore);
		assert(0);
	}
	fprintf(f, "Total LB:            %6d\n", totalLB);
	if (totalLB != ptn->scoreLB) {
		fprintf(f, "ZOINKS!  The total LB stored inside the partition is %d!\n", ptn->scoreLB);
		assert(0);
	}
	fprintf(f, "Total UB:            %6d\n", totalUB);
	if (totalUB != ptn->scoreUB) {
		fprintf(f, "ZOINKS!  The total UB stored inside the partition is %d!\n", ptn->scoreUB);
		assert(0);
	}
	if (ptn->nParts) {
		fprintf(f, "Narrowest part:      %6d with width %6d\n", narrowestPart, ptn->parts[narrowestPart].pt.nSites);
		fprintf(f, "Widest part:         %6d with width %6d\n", widestPart, ptn->parts[widestPart].pt.nSites);
		fprintf(f, "Average width:       %6.2lf\n", (double) totalWidth / ptn->nParts);
	}
}

// And, actually calculate upper bounds on those parts, and update the boundRest array!
// And, run repeatedly until all duplicate sites have been used up!
// Requires that sites are ordered by weight.
//HACK: Split up the TreeData * argument.
int PartBoundSubtree(struct TreeData *td, struct scored_part fullAlignment, int nTaxaOnTree) {
	int i;
	int *siteCounts;
	struct scored_partition ptn;
	int totalScoreLB = 0;
	int totalScoreUB = 0;
	int totalClippedScore = 0;
	int totalFinalScore = 0;
	int rightmostNonzeroSite = td->seqLen - 1;
	int multiplier;
	int totalMultiplier = 0;
	int iChunk;
#ifdef PARTBOUND_DUMP_PARTITIONS
	char fname[80];
	FILE *partdumpfile;
#endif	// PARTBOUND_DUMP_PARTITIONS

	DBGPRINT3(" PartBoundSubtree(nTaxa=%d, nTaxaOnTree=%d) called.\n", td->numTaxa, nTaxaOnTree);

	// Make repeated passes, each time extracting 1 copy of each remaining site into a "trivial partition" where each site
	// is in its own part, until all sites have been used up.  If the smallest site weight is > 1, we can speed up the
	// process up by removing that many sites in 1 go.
	ptn.parts = malloc(td->seqLen * sizeof (struct scored_part));
	siteCounts = malloc(td->seqLen * sizeof (int));
	for (i = 0; i < td->seqLen; ++i) {
		siteCounts[i] = td->weights[i];
		assert(!i || td->weights[i] <= td->weights[i - 1]);		//DEBUG: Ensure sites ordered by weight
	}

	//while (td->weights[0] > totalMultiplier) {
	for (iChunk = 0; td->weights[0] > totalMultiplier; ++iChunk) {
		// Create the new partition.  Need to add sites in reverse order to capitalise on weight ordering of sites.
		DBGPRINT2(" PartBoundSubtree(): Beginning chunk #%d.\n", iChunk + 1);
		ptn.scoreLB = 0;
		ptn.scoreUB = 0;
		ptn.totalFinalScore = 0;
		multiplier = td->weights[rightmostNonzeroSite] - totalMultiplier;
		totalMultiplier += multiplier;
		ptn.nParts = 0;
		for (i = rightmostNonzeroSite; i >= 0; --i) {
			ptn.parts[ptn.nParts].pt.nSites = 1;
			ptn.parts[ptn.nParts].pt.memWidth = 0;
			ptn.parts[ptn.nParts].pt.data = NULL;
			PartResize(&ptn.parts[ptn.nParts].pt, td->numTaxa, 1 + PARTBOUND_MAX_SITES_MOVED);
			PartCopySites(ptn.parts[ptn.nParts].pt, 0, fullAlignment.pt, i, td->numTaxa, 1);

			// Perform an initial measurement
			ptn.parts[ptn.nParts].scoreLB = ComputePartLowerBound_b1_w1(ptn.parts[ptn.nParts].pt, td->numTaxa);
			ptn.parts[ptn.nParts].scoreUB = ComputePartUpperBound_b1_w1(ptn.parts[ptn.nParts].pt, nTaxaOnTree);
			ptn.parts[ptn.nParts].scoreAdditive = ComputePartAdditiveScore_b1_w1(ptn.parts[ptn.nParts].pt, nTaxaOnTree, td->numTaxa);
			ptn.parts[ptn.nParts].touched = 1;		// Currently ignored by OptimisePartition()'s 1st pass, but anyway.
			ptn.parts[ptn.nParts].bestBoundsMask = gBestBoundMask;

			// Update partition totals
			ptn.scoreLB += ptn.parts[ptn.nParts].scoreLB;
			ptn.scoreUB += ptn.parts[ptn.nParts].scoreUB;
			ptn.totalFinalScore += selectFinalScore(td->partBoundStrategy, ptn.parts[ptn.nParts].scoreLB, ptn.parts[ptn.nParts].scoreUB, ptn.parts[ptn.nParts].scoreAdditive);

			++ptn.nParts;

			// Find the rightmost nonzero site for the next pass
			if (i < td->seqLen - 1 && td->weights[i] - totalMultiplier > 0 && td->weights[i + 1] - totalMultiplier == 0) {
				rightmostNonzeroSite = i;
			}
		}

		DBGPRINT5(" Before calling OptimisePartition(): LB=%d, UB=%d, totalClippedScore=%d, totalFinalScore=%d.\n", ptn.scoreLB, ptn.scoreUB, getTotalClippedScore(&ptn), ptn.totalFinalScore);

#ifdef PARTBOUND_DUMP_PARTITIONS
		sprintf(fname, "partDump.%dnToT.chunk%d.before.part", nTaxaOnTree, iChunk + 1);
		partdumpfile = fopen(fname, "w");
		PartitionDump(&ptn, td->numTaxa, partdumpfile);
		PartitionDumpStats(&ptn, td->partBoundStrategy, partdumpfile);
		fclose(partdumpfile);
#endif	// PARTBOUND_DUMP_PARTITIONS

		DBGPRINT3(" Calling OptimisePartition() with a %d-part partition that will be multiplied by %d...\n", ptn.nParts, multiplier);
		OptimisePartition(&ptn, td->numTaxa, nTaxaOnTree, 2, td->partBoundStrategy, td->partBoundMaxItersNoImprovement, ptn.nParts);		// Max dist between any seqs = # of sites = # of parts

		totalScoreLB += ptn.scoreLB * multiplier;
		totalScoreUB += ptn.scoreUB * multiplier;
		totalClippedScore += getTotalClippedScore(&ptn) * multiplier;
		totalFinalScore += ptn.totalFinalScore * multiplier;
		DBGPRINT6(" PartBoundSubtree(): After chunk #%d: Total LB = %d, total UB = %d, total clipped score = %d, total final score = %d.\n", iChunk + 1, totalScoreLB, totalScoreUB, totalClippedScore, totalFinalScore);

#ifdef PARTBOUND_DUMP_PARTITIONS
		sprintf(fname, "partDump.%dnToT.chunk%d.after.part", nTaxaOnTree, iChunk + 1);
		partdumpfile = fopen(fname, "w");
		PartitionDump(&ptn, td->numTaxa, partdumpfile);
		PartitionDumpStats(&ptn, td->partBoundStrategy, partdumpfile);
		fclose(partdumpfile);
#endif	// PARTBOUND_DUMP_PARTITIONS

		for (i = 0; i < ptn.nParts; ++i) {
			PartResize(&ptn.parts[i].pt, td->numTaxa, 0);		// Deallocate
		}
	}

	free(siteCounts);
	free(ptn.parts);
	return totalFinalScore;
}

void InitPartBoundRest_b1_w1(struct TreeData *td) {
	struct scored_part fullAlignment = { { td->seqLen, td->seqLen, td->charMat }, 0, 0, 0, 0, 0 };
	int fullAlignmentLB;
	int i;
	int bound;

	// For kicks, let's get an LB on the entire alignment.  For any realistic-size alignment, the 1/2-MST bound
	// will be used.
	fullAlignmentLB = ComputePartLowerBound_b1_w1(fullAlignment.pt, td->numTaxa);
	DBGPRINT2("InitPartBoundRest_b1_w1(): LB for entire alignment = %d.\n", fullAlignmentLB);

	for (i = 4; i < td->numTaxa; ++i) {
		DBGPRINT2("InitPartBoundRest_b1_w1(): About to calculate restBound for trees with <= %d taxa...\n", i);
		bound = PartBoundSubtree(td, fullAlignment, i);
		DBGPRINT3("InitPartBoundRest_b1_w1(): Our best estimate for trees with <= %d taxa is %d.\n", i, bound);
		if (bound > td->restBound[i - 1]) {
			DBGPRINT4("This beats the existing bound for <= %d taxa of %d by %d.\n", i, td->restBound[i - 1], bound - td->restBound[i - 1]);
			td->restBound[i - 1] = bound;
		}
	}
}

void generateSitePairCostsFile(void) {
	unsigned index, i;
	struct part pt = { 2, 0, NULL };
	int nSeqs, score;
	FILE *outfile;

	DBGPRINT1("Generating a new sitepaircosts.c file...\n");

	outfile = fopen("sitepaircosts.c", "w");
	fprintf(outfile,
		"/* GENERATED using generateSitePairCostsFile().\n"
		" * See sitepaircosts.h for an explanation of the purpose of this array.\n"
		" * $Id: partbound.c 430 2010-03-13 13:24:59Z wtwhite $\n"
		" */\n"
		"\n"
		"unsigned char sitePairCost[] = {\n"
	);

	PartResize(&pt, 16, 2);

	// Loop logic is tricky, since unsigned could in theory be 16 bits...
	index = 1;
	do {
		// Produce a single part containing this site pair.
		nSeqs = 0;
		for (i = 0; i < 16; ++i) {
			if (index & (1 << i)) {
				pt.data[nSeqs * pt.memWidth] = 1 << (i >> 2);
				pt.data[nSeqs * pt.memWidth + 1] = 1 << (i & 3);
				++nSeqs;
			}
		}

		// Score this site pair.
		score = ComputePartLowerBound_b1_w1(pt, nSeqs);
		fprintf(outfile, "\t%d%s\t/* #%u */\n", score, index != 0xFFFF ? "," : "", index);
		++index;
	} while (index & 0xFFFF);

	PartResize(&pt, 0, 0);
	fprintf(outfile, "};\n");
	fclose(outfile);
}
