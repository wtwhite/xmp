// This file is an X-macro that generates a family of Fitch algorithm functions.
//
// ScoreFitch_XYZ(), where XYZ is one of:
//    b2_w4_fastc
//    stopearly_b2_w4_fastc
//    fw1_b2_w4_fastc
//    fw1_stopearly_b2_w4_fastc
//
// ScoreFitchWeight1_XYZ(), where XYZ is one of:
//    b2_w4_fastc
//    stopearly_b2_w4_fastc
//
// #define PARAM_HEADER before #include-ing to produce declarations (for inclusion in a header);
// leave it undefined to produce actual definitions (for inclusion in a source code module).

#include "fitchtemplate_b2_w16_sse2asm.inc"
#define PARAM_STOPEARLY
#include "fitchtemplate_b2_w16_sse2asm.inc"
#define PARAM_FASTWEIGHT1FITCH
#include "fitchtemplate_b2_w16_sse2asm.inc"
#undef PARAM_STOPEARLY
#include "fitchtemplate_b2_w16_sse2asm.inc"
#undef PARAM_FASTWEIGHT1FITCH
