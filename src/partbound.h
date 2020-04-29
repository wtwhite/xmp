#ifndef __PARTBOUND_H
#define __PARTBOUND_H
struct TreeData;			// Forward declaration

enum partbound_strategy_type {
	PBS_CLIPPEDANDADDITIVE,	// The standard, usually-best strategy
	PBS_CLIPPED,			// The old (pre-rev-376) PARTBOUND strategy
	PBS_LB					// This strategy can be used for comparison with the minmax squeeze program.
};

struct part {
	int nSites;
	int memWidth;		// Must be >= nSites always.
	char *data;
};

// We don't store the raw score or the clipped score here, since they can be trivially computed from
// scoreLB and scoreUB.  We also don't store the "final" score, which is simply the maximum of the
// clipped score and scoreAdditive.
struct scored_part {
	struct part pt;
	int scoreLB;
	int scoreUB;
	int scoreAdditive;	// The minimum score that must be added by all remaining taxa
	int touched;		// Was this part modified in the previous loop iteration?
	unsigned bestBoundsMask;		// Each 1 bit indicates a bound type that produced a best-equal LB.
};

struct scored_partition {
	struct scored_part *parts;
	int nParts;
	int scoreLB;
	int scoreUB;
	int totalFinalScore;	// For each part, usually this takes the maximum of the clipped score and scoreAdditive.
};

void InitPartBoundRest_b1_w1(struct TreeData *td);
int ComputePartLowerBound_b1_w1(struct part pt, int nTaxa);
void generateSitePairCostsFile(void);		//HACK: Should only ever be called once!
#endif	// #include guard
