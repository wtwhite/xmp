#include "common.h"
#include "seq.h"
#include "measure.h"

///////////////////////////////////////////////////////////////////////////////
// Functions for all _w8 versions of FASTCFITCH Fitch computations.
///////////////////////////////////////////////////////////////////////////////

// This function does not compute weights, only preliminary base sets.
// len should be a 64-bit word count.
void FitchBases_b2_w8_fastc64(uint64_t *s1, uint64_t *s2, uint64_t *dest, unsigned nBlocks) {
	unsigned i;
	uint64_t needUnion;
	
	INCMEASURE(fitchBasesCount);
	
	for (i = 0; i < nBlocks; ++i) {
		// Use complex expressions: let the compiler figure out the best order
		// of computation and the best common subexpressions to reuse (mostly).
		needUnion = ((((s1[i] & s2[i] & 0x7777777777777777ULL) + 0x7777777777777777ULL) | (s1[i] & s2[i])) & 0x8888888888888888ULL) >> 3;
		dest[i] = (s1[i] & s2[i]) | ((s1[i] | s2[i]) & ((needUnion + 0x7777777777777777ULL) ^ 0x8888888888888888ULL));
	}
}

// This function does not compute weights, only preliminary base sets.
// len should be a 64-bit word count.
unsigned FitchBasesAndCompare_b2_w8_fastc64(uint64_t *s1, uint64_t *s2, uint64_t *dest, uint64_t *prev, unsigned nBlocks) {
	unsigned i;
	uint64_t needUnion;
	uint64_t changed = 0;
	
	INCMEASURE(fitchBasesCount);
	for (i = 0; i < nBlocks; ++i) {
		// Use complex expressions: let the compiler figure out the best order
		// of computation and the best common subexpressions to reuse (mostly).
		needUnion = ((((s1[i] & s2[i] & 0x7777777777777777ULL) + 0x7777777777777777ULL) | (s1[i] & s2[i])) & 0x8888888888888888ULL) >> 3;
		dest[i] = (s1[i] & s2[i]) | ((s1[i] | s2[i]) & ((needUnion + 0x7777777777777777ULL) ^ 0x8888888888888888ULL));
		changed |= dest[i] ^ prev[i];
	}

	// Yes, we need the "!= 0" below, since unsigned may well be shorter than uint64_t =>
	// possible truncation => declaring seqs identical when they aren't!
	return changed != 0;		// Will be 0 iff newly-created sequence is identical to prev
}

// Generate all combinations of Fitch algorithms using the X-macro technique.
// I *believe* this will ultimately be more maintainable than maintaining separate blobs of code
// for each different combination of switches.
#include "generate_all_fitch_algos_b2_w8_fastc64.inc"
