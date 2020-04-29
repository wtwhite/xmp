#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "common.h"
#include "bound.h"
#include "distbasesbound.h"
#include "incompatbound.h"
#include "mstbound.h"
#include "hittingsetbound.h"
#include "partbound.h"

void ComputeBounds(struct TreeData *td) {
	td->restBound = (unsigned *) malloc(td->numTaxa * sizeof (unsigned));
	memset(td->restBound, 0, td->numTaxa * sizeof (unsigned));
	
	if (td->options & OPT_USEDISTBASESBOUNDREST) {
		InitDistBasesPerSiteBoundRest_b1_w1(td);
	}
	if (td->options & OPT_USEINCOMPATBOUNDREST) {
		// This will *add* to td->restBound[] elements instead of taking maximums, so the *only* method
		// that may be called before it is InitDistBasesPerSiteBoundRest_b1().
		InitIncompatBoundRest_b1(td);
	}
	if (td->options & OPT_USELOADRESTBOUND) {
		LoadRestBound(td);
	}
	if (td->options & OPT_USEMSTBOUNDREST) {
		InitMstBoundRest_b1(td);
	}
	if (td->options & OPT_USEKICKASSBOUNDREST) {
		InitKickAssBoundRest(td);
	}
	if (td->options & OPT_USEKA2BOUNDREST) {
		InitKA2BoundRest(td);
	}
	if (td->options & OPT_USEGREEDYHITTINGSETBOUND) {
		InitGreedyHittingSetBoundRest(td);
	}
	if (td->options & OPT_USEPARTBOUND) {
		InitPartBoundRest_b1_w1(td);
	}

#ifdef TARGETMULTI
	// Only the master creates the restBound.txt file.
	if (!td->rank) {
#endif	// TARGETMULTI
	{	//DEBUG
		FILE *brFile;
		unsigned i;
		brFile = fopen("restBound.txt", "wt");
		fprintf(brFile, "i\tMin weight added by taxa i+1, i+2, ..., n-1.\n");
		for (i = 2; i < td->numTaxa; ++i) {
			fprintf(brFile, "%u\t%u\n", i, td->restBound[i]);
		}
		fclose(brFile);
	}
#ifdef TARGETMULTI
	}
#endif	// TARGETMULTI
}
