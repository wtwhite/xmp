// Used as an X-macro to build various Fitch algorithm variants.
// This file should not be #included directly, rather you should #include "generate_all_fitch_b2_w8_fastc64.inc",
// which will build appropriately-named declarations (if PARAM_HEADER is #defined) or definitions
// (if not) for all combinations of STOPEARLY and FASTWEIGHT1FITCH.

///////////////////////////////////////////////////////////////////////////////
// The functions past this point do the "hard work."
///////////////////////////////////////////////////////////////////////////////

#define PASTE2(x, y) x ## y
#define PASTE(x, y) PASTE2(x, y)

//#pragma message("-----Beginning X-macro!-----")

// Make sure we're not stamping on anyone's toes
#if defined(TEMP_SUFFIX1) || defined(TEMP_SUFFIX2) || defined(TEMP_SUFFIX3)
#error "At least one of the macros TEMP_SUFFIX1, TEMP_SUFFIX2 and TEMP_SUFFIX3 is already in use!"
#endif	// defined(TEMP_SUFFIX1) || defined(TEMP_SUFFIX2) || defined(TEMP_SUFFIX3)

#define TEMP_SUFFIX1 _b2_w8_fastc64

#ifdef PARAM_STOPEARLY
#define TEMP_SUFFIX2 PASTE(_stopearly, TEMP_SUFFIX1)
//#pragma message("PARAM_STOPEARLY")
#else	// not PARAM_STOPEARLY
#define TEMP_SUFFIX2 TEMP_SUFFIX1
//#pragma message("no PARAM_STOPEARLY")
#endif	// not PARAM_STOPEARLY

#ifdef PARAM_FASTWEIGHT1FITCH
#define TEMP_SUFFIX3 PASTE(_fw1, TEMP_SUFFIX2)
//#pragma message("PARAM_FASTWEIGHT1FITCH")
#else	// not PARAM_FASTWEIGHT1FITCH
#define TEMP_SUFFIX3 TEMP_SUFFIX2
//#pragma message("no PARAM_FASTWEIGHT1FITCH")
#endif	// not PARAM_FASTWEIGHT1FITCH

// We don't produce a special version for FASTWEIGHT1FITCH.
#ifndef PARAM_FASTWEIGHT1FITCH
// This function does not compute preliminary base sets, only weights.
// If the score exceeds the limit passed in and STOPEARLY is #defined, the function returns immediately.
// len should be a 32-bit word count.
unsigned PASTE(FitchScoreWeight1, TEMP_SUFFIX2)(uint64_t *s1, uint64_t *s2, unsigned nBlocks
#ifdef PARAM_STOPEARLY
	, unsigned maxScore
#endif	// PARAM_STOPEARLY
)
#ifdef PARAM_HEADER
;
#else	// not PARAM_HEADER
{
	unsigned i;
	uint64_t needUnion;
	unsigned score = 0;

	// Drop into FASTWEIGHT1FITCH-mode
	for (i = 0; i < nBlocks; ++i) {
		needUnion = ~((((s1[i] & s2[i] & 0x7777777777777777ULL) + 0x7777777777777777ULL) | (s1[i] & s2[i])) >> 3);
#ifdef PARAM_WEIGHT1MULTIPLY
		//score += ((needUnion & 0x1111111111111111ULL) * 0x1111111111111111ULL) >> 60;
		score += (((needUnion & 0x1111111111111110ULL) * 0x1111111111111111ULL) >> 60) + (needUnion & 0x0F);		// Need separate additions to avoid overflow
#else	// not PARAM_WEIGHT1MULTIPLY
		needUnion &= 0x1111111111111111ULL;
		needUnion += (needUnion >> 4);
		needUnion += (needUnion >> 8);
		needUnion += (needUnion >> 16);
		//needUnion += (needUnion >> 32);
		//score += needUnion & 0xF;
		score += (needUnion & 0x0F) + ((needUnion >> 32) & 0x0F);		// Need separate additions to avoid overflow
#endif	// not PARAM_WEIGHT1MULTIPLY
		
#ifdef PARAM_STOPEARLY
		if (score > maxScore) {
			return score;
		}
#endif	// PARAM_STOPEARLY
	}

	return score;
}
#endif	// not PARAM_HEADER
#endif	// not PARAM_FASTWEIGHT1FITCH

// This function does not compute preliminary base sets, only weights.
// If the score exceeds the limit passed in and STOPEARLY is #defined, the function returns immediately.
// len should be a 32-bit word count.
// WTJW 14/2/2005: This function correctly implements all 16 possible
// permutations of the flags STOPEARLY, FASTWEIGHT1FITCH, WEIGHT1MULTIPLY and BIGENDIAN.
unsigned PASTE(FitchScore, TEMP_SUFFIX3)(uint64_t *s1, uint64_t *s2, unsigned *weights, unsigned nBlocks
#ifdef PARAM_STOPEARLY
	, unsigned maxScore
#endif	// PARAM_STOPEARLY
)
#ifdef PARAM_HEADER
;
#else	// not PARAM_HEADER
{
	unsigned i;
	uint64_t needUnion;
	unsigned score = 0;
	
	INCMEASURE(scoreFitchCount);
	
	for (i = 0; i < nBlocks; ++i) {
#ifdef PARAM_FASTWEIGHT1FITCH
		if (weights[i * 16] == 1) {
			return score + PASTE(FitchScoreWeight1, TEMP_SUFFIX2)(s1 + i, s2 + i, nBlocks - i
#ifdef PARAM_STOPEARLY
				, maxScore
#endif	// PARAM_STOPEARLY
			);
		}
#endif	// PARAM_FASTWEIGHT1FITCH
		
		// Turn on the LSB of each nibble iff all bits in that nibble are off.
		needUnion = (((s1[i] & s2[i] & 0x7777777777777777ULL) + 0x7777777777777777ULL) | (s1[i] & s2[i])) >> 3;
		
		score +=
#ifdef BIGENDIAN
			((((needUnion >> 56) & 1) - 1) & weights[i * 16]) +
			((((needUnion >> 60) & 1) - 1) & weights[i * 16 + 1]) +
			((((needUnion >> 48) & 1) - 1) & weights[i * 16 + 2]) +
			((((needUnion >> 52) & 1) - 1) & weights[i * 16 + 3]) +
			((((needUnion >> 40) & 1) - 1) & weights[i * 16 + 4]) +
			((((needUnion >> 44) & 1) - 1) & weights[i * 16 + 5]) +
			((((needUnion >> 32) & 1) - 1) & weights[i * 16 + 6]) +
			((((needUnion >> 36) & 1) - 1) & weights[i * 16 + 7]) +
			((((needUnion >> 24) & 1) - 1) & weights[i * 16 + 8]) +
			((((needUnion >> 28) & 1) - 1) & weights[i * 16 + 9]) +
			((((needUnion >> 16) & 1) - 1) & weights[i * 16 + 10]) +
			((((needUnion >> 20) & 1) - 1) & weights[i * 16 + 11]) +
			((((needUnion >> 8) & 1) - 1)  & weights[i * 16 + 12]) +
			((((needUnion >> 12) & 1) - 1) & weights[i * 16 + 13]) +
			(((needUnion & 1) - 1)         & weights[i * 16 + 14]) +
			((((needUnion >> 4) & 1) - 1)  & weights[i * 16 + 15]);
#else	// not BIGENDIAN
			(((needUnion & 1) - 1)         & weights[i * 16]) +
			((((needUnion >> 4) & 1) - 1)  & weights[i * 16 + 1]) +
			((((needUnion >> 8) & 1) - 1)  & weights[i * 16 + 2]) +
			((((needUnion >> 12) & 1) - 1) & weights[i * 16 + 3]) +
			((((needUnion >> 16) & 1) - 1) & weights[i * 16 + 4]) +
			((((needUnion >> 20) & 1) - 1) & weights[i * 16 + 5]) +
			((((needUnion >> 24) & 1) - 1) & weights[i * 16 + 6]) +
			((((needUnion >> 28) & 1) - 1) & weights[i * 16 + 7]) +
			((((needUnion >> 32) & 1) - 1) & weights[i * 16 + 8]) +
			((((needUnion >> 36) & 1) - 1) & weights[i * 16 + 9]) +
			((((needUnion >> 40) & 1) - 1) & weights[i * 16 + 10]) +
			((((needUnion >> 44) & 1) - 1) & weights[i * 16 + 11]) +
			((((needUnion >> 48) & 1) - 1) & weights[i * 16 + 12]) +
			((((needUnion >> 52) & 1) - 1) & weights[i * 16 + 13]) +
			((((needUnion >> 56) & 1) - 1) & weights[i * 16 + 14]) +
			((((needUnion >> 60) & 1) - 1) & weights[i * 16 + 15]);
#endif	// not BIGENDIAN
		
#ifdef PARAM_STOPEARLY
		if (score > maxScore) {
			return score;
		}
#endif	// PARAM_STOPEARLY
	}
	
	return score;
}
#endif	// not PARAM_HEADER

// Clean up
#undef TEMP_SUFFIX1
#undef TEMP_SUFFIX2
#undef TEMP_SUFFIX3
