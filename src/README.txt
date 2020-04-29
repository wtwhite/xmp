FASTDNAMP implementation based on the algorithm of Yan and Bader
================================================================

This directory contains files copied across from the ka2 dir, plus extensive
modifications.  The basic idea is to use the algorithm of Yan and Bader,
described in their 2003 paper "Fast character optimization in parsimony
phylogeny reconstruction."

I have written down on paper what the key
differences and expected improvements are; to summarise, we do a full tree
recalculation (two passes, even) every time a new taxon is to be added to
the tree, but finding the cost of adding this taxon at each edge in the tree
only takes one Fitch operation per edge, leading to an amortized cost of
3 Fitch operations per B&B search node.

There is also the fact that since we are not doing any kind of incremental
computation of the tree, we are free to compact together columns that appear
the same on the first m+1 taxa, when only looking at trees containing m
taxa.  This means that the "effective width" of the character data decreases
at lower depths in the search tree!

-- WTJW 19/1/2005

8/2/2005
========

Seems to be working!  Timing on it029712 compiling with /O2 and using INCOMPATBOUND:

mt-10.phy: 17 seconds!  (Used to be 30-something seconds, MINIMUM.)

Profiling with gprof on it027879:
  %   cumulative   self              self     total
 time   seconds   seconds    calls  ms/call  ms/call  name
 59.74     10.43    10.43  2173623     0.00     0.00  ScoreFitch_b2_w8_fastc_nostopearly
 39.06     17.25     6.82  5214159     0.00     0.00  PrelimFitch_b2_w8_fastc

Clearly, it is the summing of weights that is the time-consuming part when
there are many sites.  A FASTWEIGHT1FITCH-type trick should speed it up a LOT
on this example, since about 70% of the sites have weight 1.

When STOPEARLY and FASTWEIGHT1FITCH are used, drops to 9s!

its36.dat
---------
62370 optimal trees were found, each having a score of 233.
All optimal trees were found.
Time taken: 47 minutes, 18 seconds.

The PAUP* time (as reported in its36_paup_bandb_bif_all.tre, possibly run on
an older machine) was 01:10:05.3!  We caned!

After implementing FASTWEIGHT1FITCH and STOPEARLY, and removing BACKPOINTERS:
62370 optimal trees were found, each having a score of 233.
All optimal trees were found.
Time taken: 38 minutes, 17 seconds.

Still to implement:
- Shrink representation for lower numbers of taxa
	- Lotsa work since:
		- All restBound[]-computing routines have to be rewritten
		- Column-collapsing code (and possibly other things) have to be
		  factored out of InitTreeData()
- FASTWEIGHT1FITCH equivalent
- Implement MMX and SSE2 versions in assembly
- Cleaning up things (e.g. checking that CopyTreeData() works as expected, etc.)

3/3/2005
========
Based on looking at treelog.txt for mt-10.phy, I'm now pretty certain that the
score computed for each partial tree is correct (the total weight of all sites
that are non-PI on the entire taxon set is added to this score).

Using bounds from a file sitting around called restBound_KABOUND.txt on IT029712:

62370 optimal trees were found, each having a score of 233.
All optimal trees were found.
Time taken: 31 minutes, 39 seconds.

Running on IT027879:

805.47 PrelimFitch
456.79 ScoreFitch
386.06 BandBScanEdge
  7.30 InsertNode (including malloc())
  3.14 DeleteNode (including free())

On IT027879, before using branch probabilities:
32:45

After using branch probabilities:
32:19

22/3/2005
=========
After implementing the SSE2 assembly version, we get the following.  Note that
BACKPOINTERS is defined.

62370 optimal trees were found, each having a score of 233.
All optimal trees were found.
Time taken: 30 minutes, 20 seconds.

Quite a small improvement unfortunately.
