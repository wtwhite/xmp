#ifndef __SEQ_B2_W4_FASTC
#define __SEQ_B2_W4_FASTC

#include <stdint.h>

// Declarations for FASTCFITCH functions.
void FitchBases_b2_w4_fastc(uint32_t *s1, uint32_t *s2, uint32_t *dest, unsigned nBlocks);
unsigned FitchBasesAndCompare_b2_w4_fastc(uint32_t *s1, uint32_t *s2, uint32_t *dest, uint32_t *prev, unsigned nBlocks);

// Generate declarations for all combinations of Fitch algorithms using the X-macro technique.
#define PARAM_HEADER
#include "generate_all_fitch_algos_b2_w4_fastc.inc"
#undef PARAM_HEADER

// If we have been selected as the implementation, export macros for briefer names.
#ifdef FASTCFITCH
#define FitchBases_b2_w4 FitchBases_b2_w4_fastc
#define FitchBasesAndCompare_b2_w4 FitchBasesAndCompare_b2_w4_fastc
#define FitchScore_b2_w4 FitchScore_b2_w4_fastc
#define FitchScoreWeight1_b2_w4 FitchScoreWeight1_b2_w4_fastc
#endif	// FASTCFITCH
#endif	// #include guard
