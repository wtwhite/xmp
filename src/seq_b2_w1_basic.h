#ifndef __SEQ_B2_W1_BASIC_H
#define __SEQ_B2_W1_BASIC_H
#include <stdio.h>
#include "switches.h"
unsigned Fitch(struct tree *root, struct TreeData *td);

// Some functions are identical for both word sizes: just define one function,
// and then use a #defined macro for each word size.
#define GetBaseAt_b2_w4 GetBaseAt_b2_w1
unsigned char FASTCALL GetBaseAt_b2_w1(unsigned taxon, unsigned site, struct TreeData *td);
#define PrintMaskSeq_b2_w4 PrintMaskSeq_b2_w1
void FASTCALL PrintMaskSeq_b2_w1(unsigned char *seq, unsigned memSeqLen, FILE *f);
#define CalcMstWeight_b2_w4 CalcMstWeight_b2_w1
unsigned FASTCALL CalcMstWeight_b2_w1(unsigned char *data, unsigned dataWidth, unsigned width, unsigned height, unsigned *weights);
#define CalcTreeMstWeight_b2_w4 CalcTreeMstWeight_b2_w1
unsigned FASTCALL CalcTreeMstWeight_b2_w1(struct TreeData *td);


#include <stdint.h>

// Declarations for BASICFITCH functions.
unsigned FASTCALL FitchBoth_b2_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned *weights, unsigned nBlocks);
void FASTCALL FitchBases_b2_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned nBlocks);
unsigned FASTCALL FitchScoreWeight1_b2_w1_basic(unsigned char *s1, unsigned char *s2, unsigned nBlocks);
unsigned FASTCALL FitchScore_b2_w1_basic(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned nBlocks);

// If we have been selected as the implementation, export macros for briefer names.
#ifdef BASICFITCH
#define FitchBases_b2_w1 FitchBases_b2_w1_basic
#define FitchScore_b2_w1 FitchScore_b2_w1_basic
#define FitchScoreWeight1_b2_w1 FitchScoreWeight1_b2_w1_basic
#endif	// BASICFITCH
#endif	// #include guard
