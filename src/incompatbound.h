#include "common.h"

void InitIncompatBoundRest_b1(struct TreeData *td);
int IsIncompatible_b1(unsigned site1, unsigned site2, struct TreeData *td);
unsigned CalcSitePairIncompatCost_b1t(unsigned char *site, unsigned numTaxa);
unsigned GreedyMaxMatchingScore(unsigned *incompatArray, unsigned numInc, unsigned numSites);
