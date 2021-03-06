#ifndef __SWITCHES_H
#define __SWITCHES_H
//WTJW: The following series of #defines are switches.

// Create MPI parallel version
/* #undef TARGETMULTI */

// This segment is configured by cmake.

// Comment out next #define to print loads of debugging information.
#define DEBUG
// The number of sites to store in a byte.  Either 1 or 2.  2 is faster of course.
#define SITESPERBYTE 2
// The number of bytes to store in a block.  All routines will expect sequences to
// be a multiple of this block size.
#define BYTESPERBLOCK 4
// A basic but slow version of the Fitch algorithm variants.
/* #undef BASICFITCH */
// Use (hopefully) faster C version of Fitch...() functions.
// Requires BYTESPERBLOCK to be a multiple of 4.
#define FASTCFITCH
// Use (hopefully) faster C version of Fitch...() functions that use 64-bit arithmetic.  While
// this will most likely compile and run correctly on a 32-bit machine (with compiled code
// that emulates the 64-bit arithmetic), it's likely to be slower than simply using FASTCFITCH --
// so we make this a separate implementation, to avoid slowing down e.g. the unimplemented
// routines in MMXASMFITCH or SSE2ASMFITCH which would otherwise default to it.
// Requires BYTESPERBLOCK to be a multiple of 8.
/* #undef FASTC64FITCH */
// Should we use a souped-up assembly-language SSE2 version of the Fitch algorithm?
// Requres BYTESPERBLOCK == 16.
/* #undef SSE2ASMFITCH */
// Handle C compilers (like MSVC++ 2008!) that don't support the "inline" keyword gracefully.
// Usually they will inline small functions anyway on high optimisation levels.
#define INLINE 
// MSVC++ (and possibly other compilers) support a __fastcall modifier on functions that
// will usually be faster.
#define FASTCALL __fastcall

// timer.h needs to know which header to include.
/* #undef HAS_POSIX_TIMER */
// elapsedtime.h needs to know which header to include.
/* #undef HAS_GETTIMEOFDAY */

// If #defined, ScoreFitch...() functions return as soon as the score goes over
// a prespecified maximum passed in as an argument.
#define STOPEARLY
// Optimise Fitch computations for columns having weight = 1?  This switch
// has no effect for SITESPERBYTE == 1 or BYTESPERBLOCK < 4.
#define FASTWEIGHT1FITCH
// If #defined, accumulate weight-1 columns using a single multiply operation;
// if not, use repeated shift-and-adds (3 for 32-bit words)
#define WEIGHT1MULTIPLY
// Comment out next #define to avoid counting invocations of BranchAndBound(), AddNode() and Fitch().
#define MEASURE

// For the MPI version, we have only proved (using MPI-Spin) that a version in which most MPI non-blocking sends
// use the synchronous form (MPI_Issend()) is correct.  If SYNCHRONOUSSENDS is not configured, we will use those,
// otherwise we use the possibly higher-performance standard-mode sends, which are allowed to perform buffering.
// I believe this should still result in a correct implementation.
/* #undef SYNCHRONOUSSENDS */
#ifdef SYNCHRONOUSSENDS
#define MPI_CHOSEN_Isend MPI_Issend
#else	// not SYNCHRONOUSSENDS
#define MPI_CHOSEN_Isend MPI_Isend
#endif	// not SYNCHRONOUSSENDS

// If CMDLINEOVERRIDE is defined, no other switches are defined.  This allows
// the user to specify all switches on the command line using the -D directive.
#ifndef CMDLINEOVERRIDE
//WTJW: Comment out next #define to avoid backpointers.  Needed by IncFitch(), however.
#define BACKPOINTERS 1
//WTJW: Attempt to get the OS to pump up our priority as much as possible, to
// make timing more accurate (Windows only so far).  TODO
//#define HIGHPRIORITY 1
//WTJW: Comment out next #define to avoid giving a breakdown of the BranchAndBound() invocations at every tree size on every progress update.
//#define MEASUREUPDATE 1
//#define DETECTREDUNDANTFITCH 1
//WTJW: Comment out next #define to avoid counting unnecessary invocations of Fitch().
//#define FITCHEXTMEASURE 1
//WTJW: Comment out next #define to turn off sorting of spectrum columns by
// weight, and use of a minimum additional-weight bound for rarely-occurring
// columns (TODO)
//#define SPECTRUMSORT 1
//WTJW: Should we sort columns in ascending order of weight?
//#define FLIPCOLS 1
//WTJW: Should we squeeze 8 bases into 32 bits instead of the usual 4?
#define SQUEEZEBASES 1
//WTJW: Should we record scores per 4/8 bases?
//#define SITESCORES 1
//WTJW: Should we exhaustively enumerate all trees?
//#define EXHAUSTIVE 1
//WTJW: Use smart Fitch?  (Requires backpointers)
//#define SMARTFITCH 1
//WTJW: Turns off actual processing of Fitch, only counts runs.
//#define EMPTYFITCH 1
//WTJW: Try to arrange partitions of columns in such a way as to give each
// partition a high lower bound.
//#define PARTITIONCOLUMNS 1
//WTJW: Attempt to reorder columns for faster incremental Fitch processing?
//#define ORDERCOLUMNS 1
//WTJW: Use slightly different taxon ordering where the score for each taxon
// is the minimum Hamming distance to any taxon already in the list, rather
// than the sum of Hamming distances to all taxa in the list.
#define MAXMINIORDER 1
//WTJW: Use MEGA's (and possibly PAUP*'s) max-mini taxon ordering scheme.
//#define MAXMINITREEORDER 1
//WTJW: What criteria should be used to decide the dynamic ordering of taxa?
#define DYNORDER_SUMWEIGHTS 1
//#define DYNORDER_MAXWEIGHT 1
//#define DYNORDER_SUMTHENMAX 1
//#define DYNORDER_MAXTHENSUM 1
//#define DYNORDER_QUICK 1
//WTJW: Try to optimise things for the cache?
//#define CACHEOPT 1		//BUG: not working yet -- need to store original address passed back by malloc() so we can pass it to free()
#define CACHESIZE 64
//#define CACHEADDRMASK 0x3F
#define CACHEADDRMASK ((CACHESIZE) - 1)
//#define HACKDONTCOMPUTESCORE 1
//WTJW: Look for long stretches of sites having the same weight.  Then use
// Fitch2SeqsWeight1() on these stretches, and multiply by the common weight.
// The #defined value represents the number of sites which must occur having
// the same weight, and probably should be a multiple of BYTESPERBLOCK.
//HACK: not yet implemented
//#define FASTCOMMONWEIGHTFITCH 64
//#define ALLWEIGHT1 1
//WTJW: Columns having up to and including this weight will be re-expanded
// into several weight-1 columns, in the hope that FASTWEIGHT1FITCH will give
// an overall speedup.
//#define REEXPANDUPTO 3
//WTJW: The following switch should usually be on, but must be turned off
// if POWEROFTWOWEIGHTS is in force.
#define SORTCOLSBYWEIGHT 1
//WTJW: Optimise Fitch computations for groups of 4 columns having total weight
// <= 255?  This switch works ONLY with SQUEEZEBASES turned off.
//#define FASTLOWWEIGHTFITCH 1
//WTJW: If turned on, PrintMaskSeq() will display sets of bases using their
// one-letter IUPAC code, rather than a list of characters inside square brackets.
#define SHOWSHORTBASESETS 1
//WTJW: Use dynamic ordering of taxa?  Should speed up analyses for large
// numbers of taxa significantly (I hope...)
//#define DYNAMICORDER 1
// One of the following two #defines must be turned on.  This is important
// for the FASTCFITCH optimisation.  Below, I attempt to determine this if
// it is not explicitly specified here.
//#define LITTLEENDIAN 1
//#define BIGENDIAN 1
//WTJW: For debugging.  If defined, causes PrintTree() and TreeToString() to
// produce taxon labels that have not been translated back via taxonMap[].
//#define DEBUG_SHOWUNRENAMEDTAXONLABELS
//#define OHGODFASTCFITCHISSCREWED 1
#endif	// CMDLINEOVERRIDE
//#define LOGTREES 1
//#define LOGTREESALL 1

// Disallow invalid combinations of switches.  A bit of a hack; certainly not
// exhaustive.

#if !defined(SITESPERBYTE)
#error Must define SITESPERBYTE as 1 or 2.
#endif

#if !defined(BYTESPERBLOCK)
#error Must define BYTESPERBLOCK as an integer.  Typically 1, 4, 8 or 16.
#endif

#if defined(SMARTFITCH)
#if !defined(BACKPOINTERS)
#error SMARTFITCH requires BACKPOINTERS to be specified.
#endif
#endif

#if defined(FASTWEIGHT1FITCH) && defined(FASTLOWWEIGHTFITCH)
#error Only one of FASTWEIGHT1FITCH and FASTLOWWEIGHTFITCH can be used at the same time.
#endif

#if defined(SQUEEZEBASES) && defined(FASTLOWWEIGHTFITCH)
#error You need to turn off SQUEEZEBASES to use FASTLOWWEIGHTFITCH.
#endif

#if defined(X86ASMFITCH) && defined(FASTCFITCH)
#error You can only turn on one of X86ASMFITCH and FASTCFITCH.
#endif

// Quite a few things require the BYTESPERBLOCK to be a multiple of 4.
#if defined(X86ASMFITCH) && (BYTESPERBLOCK % 4 != 0)
#error X86ASMFITCH currently requires that BYTESPERBLOCK be a multiple of 4.
#endif

#if defined(FASTWEIGHT1FITCH) && (BYTESPERBLOCK % 4 != 0)
#error FASTWEIGHT1FITCH currently requires that BYTESPERBLOCK be a multiple of 4.
#endif

#if defined(FASTLOWWEIGHTFITCH) && (BYTESPERBLOCK % 4 != 0)
#error FASTLOWWEIGHTFITCH currently requires that BYTESPERBLOCK be a multiple of 4.
#endif

#if defined(FASTCFITCH) && (BYTESPERBLOCK % 4 != 0)
#error FASTCFITCH currently requires that BYTESPERBLOCK be a multiple of 4.
#endif

#if defined(FASTCOMMONWEIGHTFITCH)
#error FASTCOMMONWEIGHTFITCH is not implemented yet.
#endif

#if defined(SPECTRUMSORT)
#error SPECTRUMSORT is not implemented.  Can''t see the point.
#endif

#if defined(MAXMINIORDER) && defined(MAXMINITREEORDER)
#error At most one of MAXMINIORDER and MAXMINITREEORDER can be specified.
#endif

//WTJW: Now I attempt to detect endianness, if it has not already been specified.
#if !defined(BIGENDIAN) && !defined(LITTLEENDIAN)
#if defined(__BYTE_ORDER)
#if __BYTE_ORDER == __BIG_ENDIAN
#define BIGENDIAN 1
#else
#if __BYTE_ORDER == __LITTLE_ENDIAN
#define LITTLEENDIAN 1
#endif
#endif
#endif
#endif

#if defined(LITTLEENDIAN)
#if defined(BIGENDIAN)
#error At most one of LITTLEENDIAN and BIGENDIAN can be specified.
#endif
#else
#if !defined(BIGENDIAN)
#error You must define either LITTLEENDIAN or BIGENDIAN.
#endif
#endif

#if defined(DYNAMICORDER)
#if defined(DYNORDER_SUMWEIGHTS) + defined(DYNORDER_MAXWEIGHT) + defined(DYNORDER_SUMTHENMAX) + defined(DYNORDER_MAXTHENSUM) != 1
#error Exactly one of DYNORDER_SUMWEIGHTS, DYNORDER_MAXWEIGHT, DYNORDER_SUMTHENMAX or DYNORDER_MAXTHENSUM must be #defined.
#endif
#endif

#if defined(PARTITIONCOLUMNS) && defined(ORDERCOLUMNS)
#error Only one of PARTITIONCOLUMNS and ORDERCOLUMNS can be #defined.
#endif

#if defined(PARTITIONCOLUMNS) && !defined(ALLWEIGHT1)
#error For now, you need to set ALLWEIGHT1 to expand everything back to weight-1 columns to use PARTITIONCOLUMNS.
#endif

//HACK
#if defined(SITESCORES) && !defined(INCFITCH)
#error For now, you need to #define INCFITCH if you want to use SITESCORES.  (This will change in future.)
#endif

//HACK
#if defined(HIGHPRIORITY)
#error HIGHPRIORITY is not implemented yet!
#endif

// If dbgprint.h is also being included, it must be included *after* this file.  That's because it depends
// on whether DEBUG is #defined or not, and we do the defining here!
#ifdef DBGPRINT1
#error If you want to include dbgprint.h, you must do so *after* including switches.h (or a file that includes switches.h).
#endif	// DBGPRINT1
#endif	// #include guard
