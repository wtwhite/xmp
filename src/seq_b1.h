unsigned Fitch(struct tree *root, struct TreeData *td);

// Some functions are identical for both word sizes: just define one function,
// and then use a #defined macro for each word size.
#define GetBaseAt_b1_w4 GetBaseAt_b1_w1
unsigned char FASTCALL GetBaseAt_b1_w1(unsigned taxon, unsigned site, struct TreeData *td);
#define PrintMaskSeq_b1_w4 PrintMaskSeq_b1_w1
void FASTCALL PrintMaskSeq_b1_w1(unsigned char *seq, unsigned memSeqLen, FILE *f);
#define CalcMstWeight_b1_w4 CalcMstWeight_b1_w1
unsigned FASTCALL CalcMstWeight_b1_w1(unsigned char *data, unsigned dataWidth, unsigned width, unsigned height, unsigned *weights, int underestimate);
#define CalcTreeMstWeight_b1_w4 CalcTreeMstWeight_b1_w1
unsigned FASTCALL CalcTreeMstWeight_b1_w1(struct TreeData *td, int underestimate);

unsigned FASTCALL StrictlyEqualFitch2SeqsCost_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len);

// These functions implement various modes of Fitch processing.
unsigned FASTCALL FitchBoth_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned *weights, unsigned nBlocks);
void FASTCALL FitchBases_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned nBlocks);
unsigned FASTCALL FitchScoreWeight1_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned nBlocks);
unsigned FASTCALL FitchScore_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len);
//unsigned FASTCALL Fitch2Seqs_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned *weights, unsigned len);
//unsigned FASTCALL Fitch2SeqsCost_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len);
//unsigned FASTCALL StrictlyEqualFitch2SeqsCost_b1_w1_basic(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len);		// WTJW 23/3/2005
//unsigned FASTCALL Fitch2Seqs_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned *weights, unsigned len);
//unsigned FASTCALL Fitch2SeqsCost_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned *weights, unsigned len);
//unsigned FASTCALL Fitch2SeqsWeight1_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned len);
//unsigned FASTCALL Fitch2SeqsWeight1Cost_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned len);
//unsigned FASTCALL Fitch2SeqsLowWeight_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned char *dest, unsigned char *byteWeights, unsigned len);
//unsigned FASTCALL Fitch2SeqsLowWeightCost_b1_w4_fastc(unsigned char *s1, unsigned char *s2, unsigned char *byteWeights, unsigned len);




// These functions are concerned with incremental Fitch processing.
