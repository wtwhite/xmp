#ifndef __SEQ_H
#define __SEQ_H
#include "switches.h"

#if defined(BASICFITCH) + defined(FASTCFITCH) + defined(FASTC64FITCH) + defined(X86ASMFITCH) + defined(MMXASMFITCH) + defined(SSE2ASMFITCH) != 1
#error Exactly one of BASICFITCH, FASTCFITCH, FASTC64FITCH, X86ASMFITCH, MMXASMFITCH and SSE2ASMFITCH must be specified.
#endif

// Include headers for all implementations.  It doesn't matter if a particular implementation is not
// supported on this platform -- its source file just won't be built.
// Each implementation header file is responsible for checking whether it is the selected implementation
// (e.g. seq_b2_sse2asm.h should check for SSE2ASMFITCH) and if so, "nominating" any Fitch functions that
// it supplies by creating macros (e.g. a macro FitchBases_b2_w16() that forwards to FitchBases_b2_w16_sse2asm()).
// This way, by omitting the final implementation suffix, the caller can choose the "selected" implementation.
// Defaults (from either BASICFITCH or FASTCFITCH) will be provided below for any macro that does not exist
// in the selected implementation.
#include "seq_b2_w1_basic.h"
#include "seq_b2_w4_fastc.h"
#include "seq_b2_w8_fastc64.h"
#include "seq_b2_sse2asm.h"

// Produce CHOSEN_...() macros depending on what config options were chosen.
// The CHOSEN_...() names are for calling during the main loop.  They can be #defined
// to incorporate extra functional constraints (currently SITESPERBYTE, BYTESPERBLOCK, FASTWEIGHT1FITCH
// and STOPEARLY) which must be obeyed by the input data.  For convenience, CHOSEN_...() macros
// for all scoring algorithm variants take a maxScore parameter which is needed for STOPEARLY but
// is otherwise discarded -- just pass INT_MAX for this if you want the full score.  Doing it this
// way avoids the need to write #ifdefs in calling code.
#define PASTE2(x, y) x ## y
#define PASTE(x, y) PASTE2(x, y)

#define FITCHSUFFIX1 PASTE(_b, PASTE(SITESPERBYTE, PASTE(_w, BYTESPERBLOCK)))

#ifdef MMXASMFITCH
#error "This hasn't been updated since big changes were made in rev 195."
#ifdef STOPEARLY
//#define ScoreFitch(a, b, c, d, e) ScoreFitch_b2_w8_mmxasm(a, b, c, d, e)
#define ScoreFitch(a, b, c, d, e) ScoreFitch_stopearly_b2_w8_mmxasm(a, b, c, d, e)
#else	// not STOPEARLY
#define ScoreFitch(a, b, c, d, e) ScoreFitch_b2_w8_mmxasm(a, b, c, d)
#endif	// not STOPEARLY
#endif	// MMXASMFITCH

#ifdef STOPEARLY
#define FITCHSUFFIX2 PASTE(_stopearly, FITCHSUFFIX1)
#else	// not STOPEARLY
#define FITCHSUFFIX2 FITCHSUFFIX1
#endif	// not STOPEARLY

#ifdef FASTWEIGHT1FITCH
#define FITCHSUFFIX3 PASTE(_fw1, FITCHSUFFIX2)
#else	// not FASTWEIGHT1FITCH
#define FITCHSUFFIX3 FITCHSUFFIX2
#endif	// not FASTWEIGHT1FITCH

// We now have everything we need to define the macros.

#define CHOSEN_FitchBases PASTE(FitchBases, FITCHSUFFIX1)
#define CHOSEN_FitchBasesAndCompare PASTE(FitchBasesAndCompare, FITCHSUFFIX1)
#ifdef STOPEARLY
#define CHOSEN_FitchScoreWeight1(a, b, c, d) PASTE(FitchScoreWeight1, FITCHSUFFIX2)(a, b, c, d)
#define CHOSEN_FitchScore(a, b, c, d, e) PASTE(FitchScore, FITCHSUFFIX3)(a, b, c, d, e)
#else	// not STOPEARLY
#define CHOSEN_FitchScoreWeight1(a, b, c, d) PASTE(FitchScoreWeight1, FITCHSUFFIX2)(a, b, c)
#define CHOSEN_FitchScore(a, b, c, d, e) PASTE(FitchScore, FITCHSUFFIX3)(a, b, c, d)
#endif	// not STOPEARLY


// For any functions not claimed, plug in a default implementation.  Usually that's the FASTCFITCH implementation.

// _b1_w1: Default is BASICFITCH.
#ifndef FitchBases_b1_w1
#define FitchBases_b1_w1 FitchBases_b1_w1_basic
#endif

#ifndef FitchBasesAndCompare_b1_w1
#define FitchBasesAndCompare_b1_w1 FitchBasesAndCompare_b1_w1_basic
#endif

#ifndef FitchScore_b1_w1
#define FitchScore_b1_w1 FitchScore_b1_w1_basic
#endif

#ifndef FitchScoreWeight1_b1_w1
#define FitchScoreWeight1_b1_w1 FitchScoreWeight1_b1_w1_basic
#endif

#ifndef FitchScore_stopearly_b1_w1
#define FitchScore_stopearly_b1_w1(a, b, w, n, m) FitchScore_b1_w1(a, b, w, n)
#endif

#ifndef FitchScoreWeight1_stopearly_b1_w1
#define FitchScoreWeight1_stopearly_b1_w1(a, b, n, m) FitchScoreWeight1_b1_w1(a, b, n)
#endif

#ifndef FitchScore_fw1_b1_w1
#define FitchScore_fw1_b1_w1 FitchScore_b1_w1
#endif

#ifndef FitchScoreWeight1_fw1_b1_w1
#define FitchScoreWeight1_fw1_b1_w1 FitchScoreWeight1_b1_w1
#endif

#ifndef FitchScore_fw1_stopearly_b1_w1
#define FitchScore_fw1_stopearly_b1_w1 FitchScore_stopearly_b1_w1
#endif

#ifndef FitchScoreWeight1_fw1_stopearly_b1_w1
#define FitchScoreWeight1_fw1_stopearly_b1_w1 FitchScoreWeight1_stopearly_b1_w1
#endif


// _b1_w4: Default is whatever _b1_w1 is (since we know that will work).
#ifndef FitchBases_b1_w4
#define FitchBases_b1_w4(a, b, c, n) FitchBases_b1_w1(a, b, c, (n) * 4)
#endif

#ifndef FitchBasesAndCompare_b1_w4
#define FitchBasesAndCompare_b1_w4(a, b, c, d, n) FitchBasesAndCompare_b1_w1(a, b, c, d, (n) * 4)
#endif

#ifndef FitchScore_b1_w4
#define FitchScore_b1_w4(a, b, w, n) FitchScore_b1_w1(a, b, w, (n) * 4)
#endif

#ifndef FitchScoreWeight1_b1_w4
#define FitchScoreWeight1_b1_w4(a, b, n) FitchScoreWeight1_b1_w1(a, b, (n) * 4)
#endif

#ifndef FitchScore_stopearly_b1_w4
#define FitchScore_stopearly_b1_w4(a, b, w, n, m) FitchScore_stopearly_b1_w1(a, b, w, (n) * 4, m)
#endif

#ifndef FitchScoreWeight1_stopearly_b1_w4
#define FitchScoreWeight1_stopearly_b1_w4(a, b, n, m) FitchScoreWeight1_stopearly_b1_w1(a, b, (n) * 4, m)
#endif

#ifndef FitchScore_fw1_b1_w4
#define FitchScore_fw1_b1_w4(a, b, w, n) FitchScore_fw1_b1_w1(a, b, w, (n) * 4)
#endif

#ifndef FitchScoreWeight1_fw1_b1_w4
#define FitchScoreWeight1_fw1_b1_w4(a, b, n) FitchScoreWeight1_fw1_b1_w1(a, b, (n) * 4)
#endif

#ifndef FitchScore_fw1_stopearly_b1_w4
#define FitchScore_fw1_stopearly_b1_w4(a, b, w, n, m) FitchScore_fw1_stopearly_b1_w1(a, b, w, (n) * 4, m)
#endif

#ifndef FitchScoreWeight1_fw1_stopearly_b1_w4
#define FitchScoreWeight1_fw1_stopearly_b1_w4(a, b, n, m) FitchScoreWeight1_fw1_stopearly_b1_w1(a, b, (n) * 4, m)
#endif


// _b1_w8: Default is whatever _b1_w4 is (since we know that will work).
#ifndef FitchBases_b1_w8
#define FitchBases_b1_w8(a, b, c, n) FitchBases_b1_w4(a, b, c, (n) * 2)
#endif

#ifndef FitchBasesAndCompare_b1_w8
#define FitchBasesAndCompare_b1_w8(a, b, c, d, n) FitchBasesAndCompare_b1_w4(a, b, c, d, (n) * 2)
#endif

#ifndef FitchScore_b1_w8
#define FitchScore_b1_w8(a, b, w, n) FitchScore_b1_w4(a, b, w, (n) * 2)
#endif

#ifndef FitchScoreWeight1_b1_w8
#define FitchScoreWeight1_b1_w8(a, b, n) FitchScoreWeight1_b1_w4(a, b, (n) * 2)
#endif

#ifndef FitchScore_stopearly_b1_w8
#define FitchScore_stopearly_b1_w8(a, b, w, n, m) FitchScore_stopearly_b1_w4(a, b, w, (n) * 2, m)
#endif

#ifndef FitchScoreWeight1_stopearly_b1_w8
#define FitchScoreWeight1_stopearly_b1_w8(a, b, n, m) FitchScoreWeight1_stopearly_b1_w4(a, b, (n) * 2, m)
#endif

#ifndef FitchScore_fw1_b1_w8
#define FitchScore_fw1_b1_w8(a, b, w, n) FitchScore_fw1_b1_w4(a, b, w, (n) * 2)
#endif

#ifndef FitchScoreWeight1_fw1_b1_w8
#define FitchScoreWeight1_fw1_b1_w8(a, b, n) FitchScoreWeight1_fw1_b1_w4(a, b, (n) * 2)
#endif

#ifndef FitchScore_fw1_stopearly_b1_w8
#define FitchScore_fw1_stopearly_b1_w8(a, b, w, n, m) FitchScore_fw1_stopearly_b1_w4(a, b, w, (n) * 2, m)
#endif

#ifndef FitchScoreWeight1_fw1_stopearly_b1_w8
#define FitchScoreWeight1_fw1_stopearly_b1_w8(a, b, n, m) FitchScoreWeight1_fw1_stopearly_b1_w4(a, b, (n) * 2, m)
#endif


// _b1_w16: Default is whatever _b1_w8 is (since we know that will work).
#ifndef FitchBases_b1_w16
#define FitchBases_b1_w16(a, b, c, n) FitchBases_b1_w8(a, b, c, (n) * 2)
#endif

#ifndef FitchBasesAndCompare_b1_w16
#define FitchBasesAndCompare_b1_w16(a, b, c, d, n) FitchBasesAndCompare_b1_w8(a, b, c, d, (n) * 2)
#endif

#ifndef FitchScore_b1_w16
#define FitchScore_b1_w16(a, b, w, n) FitchScore_b1_w8(a, b, w, (n) * 2)
#endif

#ifndef FitchScoreWeight1_b1_w16
#define FitchScoreWeight1_b1_w16(a, b, n) FitchScoreWeight1_b1_w8(a, b, (n) * 2)
#endif

#ifndef FitchScore_stopearly_b1_w16
#define FitchScore_stopearly_b1_w16(a, b, w, n, m) FitchScore_stopearly_b1_w8(a, b, w, (n) * 2, m)
#endif

#ifndef FitchScoreWeight1_stopearly_b1_w16
#define FitchScoreWeight1_stopearly_b1_w16(a, b, n, m) FitchScoreWeight1_stopearly_b1_w8(a, b, (n) * 2, m)
#endif

#ifndef FitchScore_fw1_b1_w16
#define FitchScore_fw1_b1_w16(a, b, w, n) FitchScore_fw1_b1_w8(a, b, w, (n) * 2)
#endif

#ifndef FitchScoreWeight1_fw1_b1_w16
#define FitchScoreWeight1_fw1_b1_w16(a, b, n) FitchScoreWeight1_fw1_b1_w8(a, b, (n) * 2)
#endif

#ifndef FitchScore_fw1_stopearly_b1_w16
#define FitchScore_fw1_stopearly_b1_w16(a, b, w, n, m) FitchScore_fw1_stopearly_b1_w8(a, b, w, (n) * 2, m)
#endif

#ifndef FitchScoreWeight1_fw1_stopearly_b1_w16
#define FitchScoreWeight1_fw1_stopearly_b1_w16(a, b, n, m) FitchScoreWeight1_fw1_stopearly_b1_w8(a, b, (n) * 2, m)
#endif


// _b2_w1: Default is BASICFITCH.
#ifndef FitchBases_b2_w1
#define FitchBases_b2_w1 FitchBases_b2_w1_basic
#endif

#ifndef FitchBasesAndCompare_b2_w1
#define FitchBasesAndCompare_b2_w1 FitchBasesAndCompare_b2_w1_basic
#endif

#ifndef FitchScore_b2_w1
#define FitchScore_b2_w1 FitchScore_b2_w1_basic
#endif

#ifndef FitchScoreWeight1_b2_w1
#define FitchScoreWeight1_b2_w1 FitchScoreWeight1_b2_w1_basic
#endif

#ifndef FitchScore_stopearly_b2_w1
#define FitchScore_stopearly_b2_w1(a, b, w, n, m) FitchScore_b2_w1(a, b, w, n)
#endif

#ifndef FitchScoreWeight1_stopearly_b2_w1
#define FitchScoreWeight1_stopearly_b2_w1(a, b, n, m) FitchScoreWeight1_b2_w1(a, b, n)
#endif

#ifndef FitchScore_fw1_b2_w1
#define FitchScore_fw1_b2_w1 FitchScore_b2_w1
#endif

#ifndef FitchScoreWeight1_fw1_b2_w1
#define FitchScoreWeight1_fw1_b2_w1 FitchScoreWeight1_b2_w1
#endif

#ifndef FitchScore_fw1_stopearly_b2_w1
#define FitchScore_fw1_stopearly_b2_w1 FitchScore_stopearly_b2_w1
#endif

#ifndef FitchScoreWeight1_fw1_stopearly_b2_w1
#define FitchScoreWeight1_fw1_stopearly_b2_w1 FitchScoreWeight1_stopearly_b2_w1
#endif


// _b2_w4: Default is FASTCFITCH.
#ifndef FitchBases_b2_w4
#define FitchBases_b2_w4 FitchBases_b2_w4_fastc
#endif

#ifndef FitchBasesAndCompare_b2_w4
#define FitchBasesAndCompare_b2_w4 FitchBasesAndCompare_b2_w4_fastc
#endif

#ifndef FitchScore_b2_w4
#define FitchScore_b2_w4 FitchScore_b2_w4_fastc
#endif

#ifndef FitchScoreWeight1_b2_w4
#define FitchScoreWeight1_b2_w4 FitchScoreWeight1_b2_w4_fastc
#endif

#ifndef FitchScore_stopearly_b2_w4
#define FitchScore_stopearly_b2_w4 FitchScore_stopearly_b2_w4_fastc
#endif

#ifndef FitchScoreWeight1_stopearly_b2_w4
#define FitchScoreWeight1_stopearly_b2_w4 FitchScoreWeight1_stopearly_b2_w4_fastc
#endif

#ifndef FitchScore_fw1_b2_w4
#define FitchScore_fw1_b2_w4 FitchScore_fw1_b2_w4_fastc
#endif

#ifndef FitchScoreWeight1_fw1_b2_w4
#define FitchScoreWeight1_fw1_b2_w4 FitchScoreWeight1_b2_w4
#endif

#ifndef FitchScore_fw1_stopearly_b2_w4
#define FitchScore_fw1_stopearly_b2_w4 FitchScore_fw1_stopearly_b2_w4_fastc
#endif

#ifndef FitchScoreWeight1_fw1_stopearly_b2_w4
#define FitchScoreWeight1_fw1_stopearly_b2_w4 FitchScoreWeight1_stopearly_b2_w4
#endif


// _b2_w8: Default is whatever _b2_w4 is (since we know that will work).
#ifndef FitchBases_b2_w8
#define FitchBases_b2_w8(a, b, c, n) FitchBases_b2_w4(a, b, c, (n) * 2)
#endif

#ifndef FitchBasesAndCompare_b2_w8
#define FitchBasesAndCompare_b2_w8(a, b, c, d, n) FitchBasesAndCompare_b2_w4(a, b, c, d, (n) * 2)
#endif

#ifndef FitchScore_b2_w8
#define FitchScore_b2_w8(a, b, w, n) FitchScore_b2_w4(a, b, w, (n) * 2)
#endif

#ifndef FitchScoreWeight1_b2_w8
#define FitchScoreWeight1_b2_w8(a, b, n) FitchScoreWeight1_b2_w4(a, b, (n) * 2)
#endif

#ifndef FitchScore_stopearly_b2_w8
#define FitchScore_stopearly_b2_w8(a, b, w, n, m) FitchScore_stopearly_b2_w4(a, b, w, (n) * 2, m)
#endif

#ifndef FitchScoreWeight1_stopearly_b2_w8
#define FitchScoreWeight1_stopearly_b2_w8(a, b, n, m) FitchScoreWeight1_stopearly_b2_w4(a, b, (n) * 2, m)
#endif

#ifndef FitchScore_fw1_b2_w8
#define FitchScore_fw1_b2_w8(a, b, w, n) FitchScore_fw1_b2_w4(a, b, w, (n) * 2)
#endif

#ifndef FitchScoreWeight1_fw1_b2_w8
#define FitchScoreWeight1_fw1_b2_w8(a, b, n) FitchScoreWeight1_fw1_b2_w4(a, b, (n) * 2)
#endif

#ifndef FitchScore_fw1_stopearly_b2_w8
#define FitchScore_fw1_stopearly_b2_w8(a, b, w, n, m) FitchScore_fw1_stopearly_b2_w4(a, b, w, (n) * 2, m)
#endif

#ifndef FitchScoreWeight1_fw1_stopearly_b2_w8
#define FitchScoreWeight1_fw1_stopearly_b2_w8(a, b, n, m) FitchScoreWeight1_fw1_stopearly_b2_w4(a, b, (n) * 2, m)
#endif

// _b2_w16: Default is whatever _b2_w8 is (since we know that will work).
#ifndef FitchBases_b2_w16
#define FitchBases_b2_w16(a, b, c, n) FitchBases_b2_w8(a, b, c, (n) * 2)
#endif

#ifndef FitchBasesAndCompare_b2_w16
#define FitchBasesAndCompare_b2_w16(a, b, c, d, n) FitchBasesAndCompare_b2_w8(a, b, c, d, (n) * 2)
#endif

#ifndef FitchScore_b2_w16
#define FitchScore_b2_w16(a, b, w, n) FitchScore_b2_w8(a, b, w, (n) * 2)
#endif

#ifndef FitchScoreWeight1_b2_w16
#define FitchScoreWeight1_b2_w16(a, b, n) FitchScoreWeight1_b2_w8(a, b, (n) * 2)
#endif

#ifndef FitchScore_stopearly_b2_w16
#define FitchScore_stopearly_b2_w16(a, b, w, n, m) FitchScore_stopearly_b2_w8(a, b, w, (n) * 2, m)
#endif

#ifndef FitchScoreWeight1_stopearly_b2_w16
#define FitchScoreWeight1_stopearly_b2_w16(a, b, n, m) FitchScoreWeight1_stopearly_b2_w8(a, b, (n) * 2, m)
#endif

#ifndef FitchScore_fw1_b2_w16
#define FitchScore_fw1_b2_w16(a, b, w, n) FitchScore_fw1_b2_w8(a, b, w, (n) * 2)
#endif

#ifndef FitchScoreWeight1_fw1_b2_w16
#define FitchScoreWeight1_fw1_b2_w16(a, b, n) FitchScoreWeight1_fw1_b2_w8(a, b, (n) * 2)
#endif

#ifndef FitchScore_fw1_stopearly_b2_w16
#define FitchScore_fw1_stopearly_b2_w16(a, b, w, n, m) FitchScore_fw1_stopearly_b2_w8(a, b, w, (n) * 2, m)
#endif

#ifndef FitchScoreWeight1_fw1_stopearly_b2_w16
#define FitchScoreWeight1_fw1_stopearly_b2_w16(a, b, n, m) FitchScoreWeight1_fw1_stopearly_b2_w8(a, b, (n) * 2, m)
#endif

// Because C preprocessor macros are expanded "dynamically" (i.e. they always expand to their
// defined expansion text, which may be further expanded if it contains macro calls),
// we must leave FITCHSUFFIX1 etc. defined.
//#undef FITCHSUFFIX1
//#undef FITCHSUFFIX2
//#undef FITCHSUFFIX3






// WTJW 23/3/2005
#define StrictlyEqualFitch2SeqsCost_b1 StrictlyEqualFitch2SeqsCost_b1_w1_basic

// ===============
// IMPLEMENTATIONS
// ===============
// There are currently 3 different "implementations" of the Fitch routines.  All
// provide a pair of "generic" methods, having names starting with "Fitch2Seqs"
// and "Fitch2SeqsCost", as well as possibly more methods for dealing with
// certain common situations (such as all sites having a weight of 1).

// The BASICFITCH implementation is slow but easy-to-understand (and presumably
// reliable).  It does not provide any special functions for handling sequences
// with low weights or all weights equal to 1.
// All IncFitch() variations are currently unimplemented.
#ifdef BASICFITCH
#define Fitch2Seqs_b1_w1                Fitch2Seqs_b1_w1_basic
#define Fitch2SeqsCost_b1_w1            Fitch2SeqsCost_b1_w1_basic
#define Fitch2Seqs_b1_w4                Fitch2Seqs_b1_w1_basic								// Safe to reuse w1 version for w4
#define Fitch2SeqsCost_b1_w4            Fitch2SeqsCost_b1_w1_basic							// Safe to reuse w1 version for w4
#define Fitch2SeqsWeight1_b1_w1         Fitch2SeqsWeight1_b1_w1_basic_IS_NOT_ALLOWED		// None of the accelerated versions are allowed for BASICFITCH
#define Fitch2SeqsWeight1Cost_b1_w1     Fitch2SeqsWeight1Cost_b1_w1_basic_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1_b1_w4         Fitch2SeqsWeight1_b1_w4_basic_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1Cost_b1_w4     Fitch2SeqsWeight1Cost_b1_w4_basic_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b1_w1       Fitch2SeqsLowWeight_b1_w1_basic_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b1_w1   Fitch2SeqsLowWeightCost_b1_w1_basic_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b1_w4       Fitch2SeqsLowWeight_b1_w4_basic_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b1_w4   Fitch2SeqsLowWeightCost_b1_w4_basic_IS_NOT_ALLOWED
#define IncFitch_b1_w1                  IncFitch_b1_w1_basic_IS_NOT_ALLOWED
#define IncFitchRoot_b1_w1              IncFitchRoot_b1_w1_basic_IS_NOT_ALLOWED
#define IncFitch_b1_w4                  IncFitch_b1_w4_basic_IS_NOT_YET_IMPLEMENTED			//HACK! (use FASTCFITCH version)
#define IncFitchRoot_b1_w4              IncFitchRoot_b1_w4_basic_IS_NOT_YET_IMPLEMENTED		//HACK! (use FASTCFITCH version)
#define Fitch2Seqs_b2_w1                Fitch2Seqs_b2_w1_basic
#define Fitch2SeqsCost_b2_w1            Fitch2SeqsCost_b2_w1_basic
#define Fitch2Seqs_b2_w4                Fitch2Seqs_b2_w1_basic								// Reuse w1 version for w4
#define Fitch2SeqsCost_b2_w4            Fitch2SeqsCost_b2_w1_basic							// Reuse w1 version for w4
#define Fitch2SeqsWeight1_b2_w1         Fitch2SeqsWeight1_b2_w1_basic_IS_NOT_ALLOWED		// None of the accelerated versions are allowed for BASICFITCH
#define Fitch2SeqsWeight1Cost_b2_w1     Fitch2SeqsWeight1Cost_b2_w1_basic_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1_b2_w4         Fitch2SeqsWeight1_b2_w4_basic_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1Cost_b2_w4     Fitch2SeqsWeight1Cost_b2_w4_basic_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b2_w1       Fitch2SeqsLowWeight_b2_w1_basic_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b2_w1   Fitch2SeqsLowWeightCost_b2_w1_basic_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b2_w4       Fitch2SeqsLowWeight_b2_w4_basic_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b2_w4   Fitch2SeqsLowWeightCost_b2_w4_basic_IS_NOT_ALLOWED
#define IncFitch_b2_w1                  IncFitch_b2_w1_basic_IS_NOT_ALLOWED
#define IncFitchRoot_b2_w1              IncFitchRoot_b2_w1_basic_IS_NOT_ALLOWED
#define IncFitch_b2_w4                  IncFitch_b2_w4_basic_IS_NOT_YET_IMPLEMENTED			//HACK! (use FASTCFITCH version)
#define IncFitchRoot_b2_w4              IncFitchRoot_b2_w4_basic_IS_NOT_YET_IMPLEMENTED		//HACK! (use FASTCFITCH version)
// The "fast" version of the standard operations, taking two extra parameters, "heavySeqLen"
// and "byteWeights", is just the original version.  Both parameters are ignored,
// and are eliminated from the token stream passed to the compiler (so they need
// not refer to actual symbols).
#define FastFitch2Seqs(s1, s2, dest, weights, len, heavySeqLen, byteWeights) Fitch2Seqs((s1), (s2), (dest), (weights), (len))
#define FastFitch2SeqsCost(s1, s2, weights, len, heavySeqLen, byteWeights) Fitch2SeqsCost((s1), (s2), (weights), (len))
#endif	// BASICFITCH

// The FASTCFITCH implementation is fast and portable, using a variety of tricks
// that significantly reduce performance losses due to branch mispredictions in
// typical modern processors.  Functions are provided for speedier handling of
// sequences where all weights are 1, or sequences where the weights are all
// less than 16 (the latter requires that there be only 1 base per byte, i.e.
// SQUEEZEBASES is NOT #defined.)  This implementation requires a word size that
// is a multiple of 4.
#ifdef FASTCFITCH
#if BYTESPERBLOCK % 4 != 0
#error FASTCFITCH requires that the BYTESPERBLOCK be set to a multiple of 4.
#endif	// BYTESPERBLOCK % 4 != 0
//#define Fitch2Seqs_b2_w8				Fitch2Seqs_b2_w8_fastc
//#define Fitch2SeqsCost_b2_w8			Fitch2SeqsCost_b2_w8_fastc
//#define Fitch2SeqsWeight1_b2_w8         Fitch2SeqsWeight1_b2_w8_fastc
//#define Fitch2SeqsWeight1Cost_b2_w8     Fitch2SeqsWeight1Cost_b2_w8_fastc
#define Fitch2Seqs_b2_w8				Fitch2Seqs_b2_w4_fastc
#define Fitch2SeqsCost_b2_w8			Fitch2SeqsCost_b2_w4_fastc
#define Fitch2SeqsWeight1_b2_w8         Fitch2SeqsWeight1_b2_w4_fastc
#define Fitch2SeqsWeight1Cost_b2_w8     Fitch2SeqsWeight1Cost_b2_w4_fastc
#define Fitch2Seqs_b1_w1                Fitch2Seqs_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsCost_b1_w1            Fitch2SeqsCost_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2Seqs_b1_w4                Fitch2Seqs_b1_w4_fastc
#define Fitch2SeqsCost_b1_w4            Fitch2SeqsCost_b1_w4_fastc
#define Fitch2SeqsWeight1_b1_w1         Fitch2SeqsWeight1_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1Cost_b1_w1     Fitch2SeqsWeight1Cost_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1_b1_w4         Fitch2SeqsWeight1_b1_w4_fastc
#define Fitch2SeqsWeight1Cost_b1_w4     Fitch2SeqsWeight1Cost_b1_w4_fastc
#define Fitch2SeqsLowWeight_b1_w1       Fitch2SeqsLowWeight_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b1_w1   Fitch2SeqsLowWeightCost_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b1_w4       Fitch2SeqsLowWeight_b1_w4_fastc
#define Fitch2SeqsLowWeightCost_b1_w4   Fitch2SeqsLowWeightCost_b1_w4_fastc
#define IncFitch_b1_w1                  IncFitch_b1_w1_fastc_IS_NOT_ALLOWED
#define IncFitchRoot_b1_w1              IncFitchRoot_b1_w1_fastc_IS_NOT_ALLOWED
#define IncFitch_b1_w4                  IncFitch_b1_w4_fastc
#define IncFitchRoot_b1_w4              IncFitchRoot_b1_w4_fastc
#define Fitch2Seqs_b2_w1                Fitch2Seqs_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsCost_b2_w1            Fitch2SeqsCost_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2Seqs_b2_w4                Fitch2Seqs_b2_w4_fastc
#define Fitch2SeqsCost_b2_w4            Fitch2SeqsCost_b2_w4_fastc
#define Fitch2SeqsWeight1_b2_w1         Fitch2SeqsWeight1_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1Cost_b2_w1     Fitch2SeqsWeight1Cost_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1_b2_w4         Fitch2SeqsWeight1_b2_w4_fastc
#define Fitch2SeqsWeight1Cost_b2_w4     Fitch2SeqsWeight1Cost_b2_w4_fastc
#define Fitch2SeqsLowWeight_b2_w1       Fitch2SeqsLowWeight_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b2_w1   Fitch2SeqsLowWeightCost_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b2_w4       Fitch2SeqsLowWeight_b2_w4_fastc_IS_NOT_ALLOWED		// The ...LowWeight() functions can't work with SQUEEZEBASES
#define Fitch2SeqsLowWeightCost_b2_w4   Fitch2SeqsLowWeightCost_b2_w4_fastc_IS_NOT_ALLOWED
#define IncFitch_b2_w1                  IncFitch_b2_w1_fastc_IS_NOT_ALLOWED
#define IncFitchRoot_b2_w1              IncFitchRoot_b2_w1_fastc_IS_NOT_ALLOWED
#define IncFitch_b2_w4                  IncFitch_b2_w4_fastc
#define IncFitchRoot_b2_w4              IncFitchRoot_b2_w4_fastc
#endif	// FASTCFITCH

// The X86ASMFITCH implementation will only run on 80386-compatible processors,
// but is potentially the fastest (although at the time of writing, most routines
// are about the same speed as their FASTCFITCH counterparts on my Pentium IV,
// and some are actually slower!)  Functions are provided for speedier handling
// of sequences where all weights are 1, or sequences where the weights are all
// less than 16 (the latter requires that there be only 1 base per byte, i.e.
// SQUEEZEBASES is NOT #defined.)  This implementation requires a word size that
// is a multiple of 4.  You must also compile and link in one or both of the
// separate file(s) seq_b<n>_x86asm.c, where n is either 1 or 2.
// All IncFitch() variations are currently unimplemented.
#ifdef X86ASMFITCH
#if BYTESPERBLOCK % 4 != 0
#error X86ASMFITCH requires that the BYTESPERBLOCK be set to a multiple of 4.
#endif	// BYTESPERBLOCK % 4 != 0
#define Fitch2Seqs_b1_w1                Fitch2Seqs_b1_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsCost_b1_w1            Fitch2SeqsCost_b1_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2Seqs_b1_w4                Fitch2Seqs_b1_w4_x86asm
#define Fitch2SeqsCost_b1_w4            Fitch2SeqsCost_b1_w4_x86asm
#define Fitch2SeqsWeight1_b1_w1         Fitch2SeqsWeight1_b1_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1Cost_b1_w1     Fitch2SeqsWeight1Cost_b1_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1_b1_w4         Fitch2SeqsWeight1_b1_w4_x86asm
#define Fitch2SeqsWeight1Cost_b1_w4     Fitch2SeqsWeight1Cost_b1_w4_x86asm
#define Fitch2SeqsLowWeight_b1_w1       Fitch2SeqsLowWeight_b1_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b1_w1   Fitch2SeqsLowWeightCost_b1_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b1_w4       Fitch2SeqsLowWeight_b1_w4_x86asm
#define Fitch2SeqsLowWeightCost_b1_w4   Fitch2SeqsLowWeightCost_b1_w4_x86asm
#define IncFitch_b1_w1                  IncFitch_b1_w1_x86asm_IS_NOT_YET_IMPLEMENTED
#define IncFitchRoot_b1_w1              IncFitchRoot_b1_w1_x86asm_IS_NOT_YET_IMPLEMENTED
#define IncFitch_b1_w4                  IncFitch_b1_w4_x86asm_IS_NOT_YET_IMPLEMENTED
#define IncFitchRoot_b1_w4              IncFitchRoot_b1_w4_x86asm_IS_NOT_YET_IMPLEMENTED
#define Fitch2Seqs_b2_w1                Fitch2Seqs_b2_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsCost_b2_w1            Fitch2SeqsCost_b2_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2Seqs_b2_w4                Fitch2Seqs_b2_w4_x86asm
#define Fitch2SeqsCost_b2_w4            Fitch2SeqsCost_b2_w4_x86asm
#define Fitch2SeqsWeight1_b2_w1         Fitch2SeqsWeight1_b2_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1Cost_b2_w1     Fitch2SeqsWeight1Cost_b2_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1_b2_w4         Fitch2SeqsWeight1_b2_w4_x86asm
#define Fitch2SeqsWeight1Cost_b2_w4     Fitch2SeqsWeight1Cost_b2_w4_x86asm
#define Fitch2SeqsLowWeight_b2_w1       Fitch2SeqsLowWeight_b2_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b2_w1   Fitch2SeqsLowWeightCost_b2_w1_x86asm_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b2_w4       Fitch2SeqsLowWeight_b2_w4_x86asm_IS_NOT_ALLOWED		// The ...LowWeight() functions can't work with SQUEEZEBASES
#define Fitch2SeqsLowWeightCost_b2_w4   Fitch2SeqsLowWeightCost_b2_w4_x86asm_IS_NOT_ALLOWED
#define IncFitch_b2_w1                  IncFitch_b2_w1_x86asm_IS_NOT_YET_IMPLEMENTED
#define IncFitchRoot_b2_w1              IncFitchRoot_b2_w1_x86asm_IS_NOT_YET_IMPLEMENTED
#define IncFitch_b2_w4                  IncFitch_b2_w4_x86asm_IS_NOT_YET_IMPLEMENTED
#define IncFitchRoot_b2_w4              IncFitchRoot_b2_w4_x86asm_IS_NOT_YET_IMPLEMENTED
#endif	// X86ASMFITCH

//HACK
#ifdef MMXASMFITCH
#if BYTESPERBLOCK % 8 != 0
#error MMXASMFITCH requires that the BYTESPERBLOCK be set to a multiple of 8.
#endif	// BYTESPERBLOCK % 8 != 0
#define Fitch2Seqs_b2_w8				Fitch2Seqs_b2_w8_mmxasm
#define Fitch2SeqsCost_b2_w8			Fitch2SeqsCost_b2_w8_mmxasm
#define Fitch2SeqsWeight1_b2_w8         Fitch2SeqsWeight1_b2_w8_mmxasm
#define Fitch2SeqsWeight1Cost_b2_w8     Fitch2SeqsWeight1Cost_b2_w8_mmxasm
//HACK: The remainder will just fall through to FASTCFITCH.
#define Fitch2Seqs_b1_w1                Fitch2Seqs_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsCost_b1_w1            Fitch2SeqsCost_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2Seqs_b1_w4                Fitch2Seqs_b1_w4_fastc
#define Fitch2SeqsCost_b1_w4            Fitch2SeqsCost_b1_w4_fastc
#define Fitch2SeqsWeight1_b1_w1         Fitch2SeqsWeight1_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1Cost_b1_w1     Fitch2SeqsWeight1Cost_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1_b1_w4         Fitch2SeqsWeight1_b1_w4_fastc
#define Fitch2SeqsWeight1Cost_b1_w4     Fitch2SeqsWeight1Cost_b1_w4_fastc
#define Fitch2SeqsLowWeight_b1_w1       Fitch2SeqsLowWeight_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b1_w1   Fitch2SeqsLowWeightCost_b1_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b1_w4       Fitch2SeqsLowWeight_b1_w4_fastc
#define Fitch2SeqsLowWeightCost_b1_w4   Fitch2SeqsLowWeightCost_b1_w4_fastc
#define IncFitch_b1_w1                  IncFitch_b1_w1_fastc_IS_NOT_ALLOWED
#define IncFitchRoot_b1_w1              IncFitchRoot_b1_w1_fastc_IS_NOT_ALLOWED
#define IncFitch_b1_w4                  IncFitch_b1_w4_fastc
#define IncFitchRoot_b1_w4              IncFitchRoot_b1_w4_fastc
#define Fitch2Seqs_b2_w1                Fitch2Seqs_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsCost_b2_w1            Fitch2SeqsCost_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2Seqs_b2_w4                Fitch2Seqs_b2_w4_fastc
#define Fitch2SeqsCost_b2_w4            Fitch2SeqsCost_b2_w4_fastc
#define Fitch2SeqsWeight1_b2_w1         Fitch2SeqsWeight1_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1Cost_b2_w1     Fitch2SeqsWeight1Cost_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsWeight1_b2_w4         Fitch2SeqsWeight1_b2_w4_fastc
#define Fitch2SeqsWeight1Cost_b2_w4     Fitch2SeqsWeight1Cost_b2_w4_fastc
#define Fitch2SeqsLowWeight_b2_w1       Fitch2SeqsLowWeight_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeightCost_b2_w1   Fitch2SeqsLowWeightCost_b2_w1_fastc_IS_NOT_ALLOWED
#define Fitch2SeqsLowWeight_b2_w4       Fitch2SeqsLowWeight_b2_w4_fastc_IS_NOT_ALLOWED		// The ...LowWeight() functions can't work with SQUEEZEBASES
#define Fitch2SeqsLowWeightCost_b2_w4   Fitch2SeqsLowWeightCost_b2_w4_fastc_IS_NOT_ALLOWED
#define IncFitch_b2_w1                  IncFitch_b2_w1_fastc_IS_NOT_ALLOWED
#define IncFitchRoot_b2_w1              IncFitchRoot_b2_w1_fastc_IS_NOT_ALLOWED
#define IncFitch_b2_w4                  IncFitch_b2_w4_fastc
#define IncFitchRoot_b2_w4              IncFitchRoot_b2_w4_fastc
#endif	// MMXASMFITCH


// ==================
// WORD SIZE DEFAULTS
// ==================
// If no "_w<n>" suffix is specified in a function's name, assume a default
// based on the word size (SEPACKSIZE) setting.
// Note: Currently, word sizes other than 1 or 4 bytes are permitted, but
// no defaults are set up for them so you will need to specify the "w_<n>" suffix.

#if BYTESPERBLOCK == 1
#define GetBaseAt_b1                    GetBaseAt_b1_w1
#define PrintMaskSeq_b1                 PrintMaskSeq_b1_w1
#define CalcMstWeight_b1                CalcMstWeight_b1_w1
#define CalcTreeMstWeight_b1            CalcTreeMstWeight_b1_w1
#define Fitch2Seqs_b1                   Fitch2Seqs_b1_w1
#define Fitch2SeqsCost_b1               Fitch2SeqsCost_b1_w1
#define Fitch2SeqsWeight1_b1            Fitch2SeqsWeight1_b1_w1							// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsWeight1Cost_b1        Fitch2SeqsWeight1Cost_b1_w1						// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsLowWeight_b1          Fitch2SeqsLowWeight_b1_w1						// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsLowWeightCost_b1      Fitch2SeqsLowWeightCost_b1_w1					// In fact this is not allowed, but a lower level handles this fact
#define GetBaseAt_b2                    GetBaseAt_b2_w1
#define PrintMaskSeq_b2                 PrintMaskSeq_b2_w1
#define CalcMstWeight_b2                CalcMstWeight_b2_w1
#define CalcTreeMstWeight_b2            CalcTreeMstWeight_b2_w1
#define Fitch2Seqs_b2                   Fitch2Seqs_b2_w1
#define Fitch2SeqsCost_b2               Fitch2SeqsCost_b2_w1
#define Fitch2SeqsWeight1_b2            Fitch2SeqsWeight1_b2_w1							// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsWeight1Cost_b2        Fitch2SeqsWeight1Cost_b2_w1						// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsLowWeight_b2          Fitch2SeqsLowWeight_b2_w1						// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsLowWeightCost_b2      Fitch2SeqsLowWeightCost_b2_w1					// In fact this is not allowed, but a lower level handles this fact
#define IncFitch_b1                     IncFitch_b1_w1									// In fact this is not yet implemented, but a lower level handles this fact
#define IncFitchRoot_b1                 IncFitchRoot_b1_w1								// In fact this is not yet implemented, but a lower level handles this fact
#define IncFitch_b2                     IncFitch_b2_w1									// In fact this is not yet implemented, but a lower level handles this fact
#define IncFitchRoot_b2                 IncFitchRoot_b2_w1								// In fact this is not yet implemented, but a lower level handles this fact
#endif	// BYTESPERBLOCK == 1

#if BYTESPERBLOCK == 4
#define GetBaseAt_b1                    GetBaseAt_b1_w4
#define PrintMaskSeq_b1                 PrintMaskSeq_b1_w4
#define CalcMstWeight_b1                CalcMstWeight_b1_w4
#define CalcTreeMstWeight_b1            CalcTreeMstWeight_b1_w4
#define Fitch2Seqs_b1                   Fitch2Seqs_b1_w4
#define Fitch2SeqsCost_b1               Fitch2SeqsCost_b1_w4
#define Fitch2SeqsWeight1_b1            Fitch2SeqsWeight1_b1_w4
#define Fitch2SeqsWeight1Cost_b1        Fitch2SeqsWeight1Cost_b1_w4
#define Fitch2SeqsLowWeight_b1          Fitch2SeqsLowWeight_b1_w4
#define Fitch2SeqsLowWeightCost_b1      Fitch2SeqsLowWeightCost_b1_w4
#define GetBaseAt_b2                    GetBaseAt_b2_w4
#define PrintMaskSeq_b2                 PrintMaskSeq_b2_w4
#define CalcMstWeight_b2                CalcMstWeight_b2_w4
#define CalcTreeMstWeight_b2            CalcTreeMstWeight_b2_w4
#define Fitch2Seqs_b2                   Fitch2Seqs_b2_w4
#define Fitch2SeqsCost_b2               Fitch2SeqsCost_b2_w4
#define Fitch2SeqsWeight1_b2            Fitch2SeqsWeight1_b2_w4
#define Fitch2SeqsWeight1Cost_b2        Fitch2SeqsWeight1Cost_b2_w4
#define Fitch2SeqsLowWeight_b2          Fitch2SeqsLowWeight_b2_w4						// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsLowWeightCost_b2      Fitch2SeqsLowWeightCost_b2_w4					// In fact this is not allowed, but a lower level handles this fact
#define IncFitch_b1                     IncFitch_b1_w4
#define IncFitchRoot_b1                 IncFitchRoot_b1_w4
#define IncFitch_b2                     IncFitch_b2_w4
#define IncFitchRoot_b2                 IncFitchRoot_b2_w4
#endif	// BYTESPERBLOCK == 4

//HACK: this is very hack
#if BYTESPERBLOCK == 8
#define GetBaseAt_b1                    GetBaseAt_b1_w4
#define PrintMaskSeq_b1                 PrintMaskSeq_b1_w4
#define CalcMstWeight_b1                CalcMstWeight_b1_w4
#define CalcTreeMstWeight_b1            CalcTreeMstWeight_b1_w4
#define Fitch2Seqs_b1                   Fitch2Seqs_b1_w4
#define Fitch2SeqsCost_b1               Fitch2SeqsCost_b1_w4
#define Fitch2SeqsWeight1_b1            Fitch2SeqsWeight1_b1_w4
#define Fitch2SeqsWeight1Cost_b1        Fitch2SeqsWeight1Cost_b1_w4
#define Fitch2SeqsLowWeight_b1          Fitch2SeqsLowWeight_b1_w4
#define Fitch2SeqsLowWeightCost_b1      Fitch2SeqsLowWeightCost_b1_w4
#define GetBaseAt_b2                    GetBaseAt_b2_w4
#define PrintMaskSeq_b2                 PrintMaskSeq_b2_w4
#define CalcMstWeight_b2                CalcMstWeight_b2_w4
#define CalcTreeMstWeight_b2            CalcTreeMstWeight_b2_w4
#define Fitch2Seqs_b2                   Fitch2Seqs_b2_w4
#define Fitch2SeqsCost_b2               Fitch2SeqsCost_b2_w4
#define Fitch2SeqsWeight1_b2            Fitch2SeqsWeight1_b2_w8				// This one is actually correct!
#define Fitch2SeqsWeight1Cost_b2        Fitch2SeqsWeight1Cost_b2_w8			// So is this one!
#define Fitch2SeqsLowWeight_b2          Fitch2SeqsLowWeight_b2_w4						// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsLowWeightCost_b2      Fitch2SeqsLowWeightCost_b2_w4					// In fact this is not allowed, but a lower level handles this fact
#define IncFitch_b1                     IncFitch_b1_w4
#define IncFitchRoot_b1                 IncFitchRoot_b1_w4
#define IncFitch_b2                     IncFitch_b2_w4
#define IncFitchRoot_b2                 IncFitchRoot_b2_w4
#endif	// BYTESPERBLOCK == 8

#if BYTESPERBLOCK == 16
#define GetBaseAt_b1                    GetBaseAt_b1_w4
#define PrintMaskSeq_b1                 PrintMaskSeq_b1_w4
#define CalcMstWeight_b1                CalcMstWeight_b1_w4
#define CalcTreeMstWeight_b1            CalcTreeMstWeight_b1_w4
#define GetBaseAt_b2                    GetBaseAt_b2_w4
#define PrintMaskSeq_b2                 PrintMaskSeq_b2_w4
#define CalcMstWeight_b2                CalcMstWeight_b2_w4
#define CalcTreeMstWeight_b2            CalcTreeMstWeight_b2_w4
#endif	// BYTESPERBLOCK == 16

// =====================
// BYTE PACKING DEFAULTS
// =====================
// If no "_b<n>" suffix is specified in a function's name, assume a default
// based on the byte packing (SQUEEZEBASES) setting.
#ifdef SQUEEZEBASES
#define BASESPERBYTE 2
#define GetBaseAt                       GetBaseAt_b2
#define PrintMaskSeq                    PrintMaskSeq_b2
#define CalcMstWeight                   CalcMstWeight_b2
#define CalcTreeMstWeight               CalcTreeMstWeight_b2
#define Fitch2Seqs                      Fitch2Seqs_b2
#define Fitch2SeqsCost                  Fitch2SeqsCost_b2
#define Fitch2SeqsWeight1               Fitch2SeqsWeight1_b2
#define Fitch2SeqsWeight1Cost           Fitch2SeqsWeight1Cost_b2
#define Fitch2SeqsLowWeight             Fitch2SeqsLowWeight_b2							// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsLowWeightCost         Fitch2SeqsLowWeightCost_b2						// In fact this is not allowed, but a lower level handles this fact
#define IncFitch                        IncFitch_b2
#define IncFitchRoot                    IncFitchRoot_b2
#else	// not SQUEEZEBASES
#define BASESPERBYTE 1
#define GetBaseAt                       GetBaseAt_b1
#define PrintMaskSeq                    PrintMaskSeq_b1
#define CalcMstWeight                   CalcMstWeight_b1
#define CalcTreeMstWeight               CalcTreeMstWeight_b1
#define Fitch2Seqs                      Fitch2Seqs_b1
#define Fitch2SeqsCost                  Fitch2SeqsCost_b1
#define Fitch2SeqsWeight1               Fitch2SeqsWeight1_b1
#define Fitch2SeqsWeight1Cost           Fitch2SeqsWeight1Cost_b1
#define Fitch2SeqsLowWeight             Fitch2SeqsLowWeight_b1							// In fact this is not allowed, but a lower level handles this fact
#define Fitch2SeqsLowWeightCost         Fitch2SeqsLowWeightCost_b1						// In fact this is not allowed, but a lower level handles this fact
#define IncFitch                        IncFitch_b1
#define IncFitchRoot                    IncFitchRoot_b1
#endif	// SQUEEZEBASES

// The fast version of the standard operations calls the ...Weight1 or ...LowWeight
// version on the "light" part of the sequence, if possible.  If we are using
// the BASICFITCH implementation, or if neither of the FAST... macros have been
// #defined, we just use the standard function, which ignores both parameters.
#if !defined(BASICFITCH) && (defined(FASTWEIGHT1FITCH) || defined(FASTLOWWEIGHTFITCH))
#ifdef FASTWEIGHT1FITCH
// Ignore byteWeights
#ifdef ALLWEIGHT1
#define FastFitch2Seqs(s1, s2, dest, weights, len, heavySeqLen, byteWeights) \
	Fitch2SeqsWeight1((s1), (s2), (dest), (len))
#define FastFitch2SeqsCost(s1, s2, weights, len, heavySeqLen, byteWeights) \
	Fitch2SeqsWeight1Cost((s1), (s2), (len))
#else	// not ALLWEIGHT1
#define FastFitch2Seqs(s1, s2, dest, weights, len, heavySeqLen, byteWeights) \
	Fitch2Seqs((s1), (s2), (dest), (weights), (heavySeqLen)) + \
	Fitch2SeqsWeight1((s1) + (heavySeqLen), (s2) + (heavySeqLen), (dest) + (heavySeqLen), (len) - (heavySeqLen))
#define FastFitch2SeqsCost(s1, s2, weights, len, heavySeqLen, byteWeights) \
	Fitch2SeqsCost((s1), (s2), (weights), (heavySeqLen)) + \
	Fitch2SeqsWeight1Cost((s1) + (heavySeqLen), (s2) + (heavySeqLen), (len) - (heavySeqLen))
#endif	// ALLWEIGHT1
#endif	// FASTWEIGHT1FITCH
#ifdef FASTLOWWEIGHTFITCH
#ifdef ALLWEIGHT1
#define FastFitch2Seqs(s1, s2, dest, weights, len, heavySeqLen) \
	Fitch2SeqsLowWeight((s1), (s2), (dest), (byteWeights), (len))
#define FastFitch2SeqsCost(s1, s2, weights, len, heavySeqLen) \
	Fitch2SeqsLowWeightCost((s1), (s2), (byteWeights), (len))
#else	// not ALLWEIGHT1
#define FastFitch2Seqs(s1, s2, dest, weights, len, heavySeqLen) \
	Fitch2Seqs((s1), (s2), (dest), (weights), (heavySeqLen)) + \
	Fitch2SeqsLowWeight((s1) + (heavySeqLen), (s2) + (heavySeqLen), (dest) + (heavySeqLen), (byteWeights), (len) - (heavySeqLen))
#define FastFitch2SeqsCost(s1, s2, weights, len, heavySeqLen) \
	Fitch2SeqsCost((s1), (s2), (weights), (heavySeqLen)) + \
	Fitch2SeqsLowWeightCost((s1) + (heavySeqLen), (s2) + (heavySeqLen), (byteWeights), (len) - (heavySeqLen))
#endif	// ALLWEIGHT1
#endif	// FASTLOWWEIGHTFITCH
#else	// not !defined(BASICFITCH) && (defined(FASTWEIGHT1FITCH) || defined(FASTLOWWEIGHTFITCH))
// Either we are using BASICFITCH, or neither of the FAST... constants is #defined.
// Ignore both parameters and just call the standard function
#define FastFitch2Seqs(s1, s2, dest, weights, len, heavySeqLen, byteWeights) Fitch2Seqs((s1), (s2), (dest), (weights), (len))
#define FastFitch2SeqsCost(s1, s2, weights, len, heavySeqLen, byteWeights) Fitch2SeqsCost((s1), (s2), (weights), (len))
#endif	// !defined(BASICFITCH) && (defined(FASTWEIGHT1FITCH) || defined(FASTLOWWEIGHTFITCH))

// Include all files, as we want to be able to access any function by specifying
// its full name.
#include "seq_b1.h"
//#include "seq_b2.h"
#include "seq_b2_sse2asm.h"			// Harmless to include these even on a non-x86 processor

#endif	// #include guard
