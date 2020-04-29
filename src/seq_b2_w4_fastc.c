#include <assert.h>
#include "common.h"
#include "seq.h"
#include "measure.h"

///////////////////////////////////////////////////////////////////////////////
// Functions for all _w4 versions of FASTCFITCH Fitch computations.
///////////////////////////////////////////////////////////////////////////////

// This function does not compute weights, only preliminary base sets.
// len should be a 32-bit word count.
void FitchBases_b2_w4_fastc(uint32_t *s1, uint32_t *s2, uint32_t *dest, unsigned nBlocks) {
	unsigned i;
	uint32_t needUnion;
	
	INCMEASURE(fitchBasesCount);
	
	for (i = 0; i < nBlocks; ++i) {
		// Use complex expressions: let the compiler figure out the best order
		// of computation and the best common subexpressions to reuse (mostly).
		needUnion = ((((s1[i] & s2[i] & 0x77777777) + 0x77777777) | (s1[i] & s2[i])) & 0x88888888) >> 3;
		dest[i] = (s1[i] & s2[i]) | ((s1[i] | s2[i]) & ((needUnion + 0x77777777) ^ 0x88888888));
	}
}

// This function does not compute weights, only preliminary base sets.
// len should be a 32-bit word count.
unsigned FitchBasesAndCompare_b2_w4_fastc(uint32_t *s1, uint32_t *s2, uint32_t *dest, uint32_t *prev, unsigned nBlocks) {
	unsigned i;
	uint32_t needUnion;
	uint32_t changed = 0;
	
	INCMEASURE(fitchBasesCount);
	
	for (i = 0; i < nBlocks; ++i) {
		// Use complex expressions: let the compiler figure out the best order
		// of computation and the best common subexpressions to reuse (mostly).
		needUnion = ((((s1[i] & s2[i] & 0x77777777) + 0x77777777) | (s1[i] & s2[i])) & 0x88888888) >> 3;
		dest[i] = (s1[i] & s2[i]) | ((s1[i] | s2[i]) & ((needUnion + 0x77777777) ^ 0x88888888));
		changed |= dest[i] ^ prev[i];
	}

	assert(sizeof (unsigned) >= sizeof (uint32_t));		// Necessary for the statement below to avoid truncation!
	return changed;		// Will be 0 iff newly-created sequence is identical to prev
}

// Generate all combinations of Fitch algorithms using the X-macro technique.
// I *believe* this will ultimately be more maintainable than maintaining separate blobs of code
// for each different combination of switches.
#include "generate_all_fitch_algos_b2_w4_fastc.inc"
