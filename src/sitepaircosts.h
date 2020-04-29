// This source file contains a single 65535-entry (64Kb) array, recording the
// parsimony scores of every possible pair of sites.  This is possible because
// duplicate pairs of bases can always be elided when just the score is
// required, so only patterns of distinct base pairs need be considered, and the
// maximum possible number of such distinct pairs of bases is 16.
// To compute an index into the array:
// 1. Take every distinct pair of bases (a, b) in your data
// 2. Set i = 0 if a = 'A', 1 if a = 'C', 2 if a = 'G' or 3 if a = 'T'
// 3. Compute j similarly for b
// 4. Turn on bit i * 4 + j in a 16-bit word.  This word becomes the index value.

extern unsigned char sitePairCost[];
