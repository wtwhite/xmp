//#define SSE2ALIGN
#ifdef _MSC_VER
#define SSE2ALIGN __declspec(align(16))
#else	// not _MSC_VER
#define SSE2ALIGN __attribute__((aligned(16)))
//#define SSE2ALIGN "FOONLY"
#endif	// not _MSC_VER

#include "switches.h"
#include <stdint.h>

// There probably is no uint128_t type, so use void * to avoid misinterpretation.
void FASTCALL FitchBases_b2_w16_sse2asm(void *s1, void *s2, void *dest, unsigned nBlocks);
unsigned FASTCALL FitchBasesAndCompare_b2_w16_sse2asm(void *s1, void *s2, void *dest, void *prev, unsigned nBlocks);
unsigned CompareSeqs_b2_w16_sse2asm(void *s1, void *s2, unsigned nBlocks);

// Generate declarations for all combinations of Fitch algorithms using the X-macro technique.
#define PARAM_HEADER
#include "generate_all_fitch_algos_b2_w16_sse2asm.inc"
#undef PARAM_HEADER

// If we have been selected as the implementation, export macros for briefer names.
#ifdef SSE2ASMFITCH
#define FitchBases_b2_w16 FitchBases_b2_w16_sse2asm
#define FitchScore_b2_w16 FitchScore_b2_w16_sse2asm
#define FitchScore_stopearly_b2_w16 FitchScore_stopearly_b2_w16_sse2asm
#define FitchScore_fw1_b2_w16 FitchScore_fw1_b2_w16_sse2asm
#define FitchScore_fw1_stopearly_b2_w16 FitchScore_fw1_stopearly_b2_w16_sse2asm
#define FitchBasesAndCompare_b2_w16 FitchBasesAndCompare_b2_w16_sse2asm
#define CompareSeqs_b2_w16 CompareSeqs_b2_w16_sse2asm
//#define FitchScoreWeight1_b2_w16 FitchScoreWeight1_b2_w16_sse2asm
#endif	// SSE2ASMFITCH
