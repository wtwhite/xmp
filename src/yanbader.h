#ifndef __YANBADER_H
#define __YANBADER_H

// Automatically generated prototype list
void FASTCALL BranchAndBound(struct tree *root, unsigned numTaxa, struct TreeData *td);
void FASTCALL BandBScanEdge(struct tree *root, struct tree *parent, int iChild, unsigned numTaxa, struct TreeData *td, unsigned *pFastForwardIdx);
void FASTCALL DeleteOldTrees(struct TreeData *td);
void FASTCALL UpdateTreeAfterInsertion(struct tree *newInternal, int nTaxaOnTree, struct TreeData *td);
#endif	// #include guard
