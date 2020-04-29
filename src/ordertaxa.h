#ifndef __ORDERTAXA_H
#define __ORDERTAXA_H
void OrderTaxa(struct TreeData *td);
void OrderTaxaScanEdge(struct tree *root, struct tree *parent, int iChild, unsigned numTaxa, struct TreeData *td);
void OrderTaxaScanEdgeThorough(struct tree *root, struct tree *parent, int iChild, unsigned numTaxa, struct TreeData *td, unsigned *distBuf, unsigned *distBufPosPtr);
void ComputeTaxonAdditionWeight(struct tree *root, struct tree *parent, int iChild, unsigned numTaxa, struct TreeData *td, unsigned *taxonScores);
void RenameTaxon(struct tree *root, int from, int to);
#endif	// #include guard
