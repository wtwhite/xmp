#ifndef __SEQ_B2_W8_FASTC64
#define __SEQ_B2_W8_FASTC64

#include <stdint.h>

// Declarations for FASTCFITCH64 functions.
void FitchBases_b2_w8_fastc64(uint64_t *s1, uint64_t *s2, uint64_t *dest, unsigned nBlocks);
unsigned FitchBasesAndCompare_b2_w8_fastc64(uint64_t *s1, uint64_t *s2, uint64_t *dest, uint64_t *prev, unsigned nBlocks);

// Generate declarations for all combinations of Fitch algorithms using the X-macro technique.
#define PARAM_HEADER
#include "generate_all_fitch_algos_b2_w8_fastc64.inc"
#undef PARAM_HEADER

// If we have been selected as the implementation, export macros for briefer names.
#ifdef FASTC64FITCH
#define FitchBases_b2_w8 FitchBases_b2_w8_fastc64
#define FitchBasesAndCompare_b2_w8 FitchBasesAndCompare_b2_w8_fastc64
#define FitchScore_b2_w8 FitchScore_b2_w8_fastc64
#define FitchScoreWeight1_b2_w8 FitchScoreWeight1_b2_w8_fastc64
#endif	// FASTC64FITCH
#endif	// #include guard
