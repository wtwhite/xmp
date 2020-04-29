#include "interface.h"

/*From the per-order encodeing, generate the true tree, for each node, find its 
parent and brother*/ 
int     encode_to_true_tree(tNode_T encode,int tree_size,T_tree_node *true_tree);


/*Generate tree*/
int gen_tree(int level, tNode_T *tree_stack, int *pos);

/*Get the Neighbor-Joining tree onthe taxa (i)...n*/
int Adjusted_SK_NJ(int *dist_matrix,tNode_T tree,int i);

/*Get the minimum distance between  one taxa from 1...(i-1) and
one from i..n*/
int compute_min_dist(int *dist_matrix,int i);

/*traverse true_tree in in-order and save the tags of visited nodes in
tree started from i, return next avaliable positionin tree*/
int in_order_walk(struct tree *true_tree, tNode_T tree, int i);

/*traverse true_tree in pre-order and save the tags of visited nodes in
tree started from i, return next avaliable position in tree*/
int pre_order_walk(struct tree *true_tree, tNode_T tree, int i);

/*change the rooted binary tree to unrooted tree, suppose a terminal as
root and return the pseudo root*/
struct tree   *rooted_to_unrooted(struct tree *true_tree);


/*insert the partial tree with lowerbound into the heap of
frontier[frontier_id]*/
struct hNode * enqueue(int frontier_id,tNode_T partial_tree,int 
lowerbound);

/*decide wether the local change on tree will exceed val, if yes, 
reurn 1*/
int execede_gap(int pos, int parent,int new_taxa);


/*compute the maximum unresolved discrepandcy*/
int compue_max_undis(tNode_T tree,int tree_size,int level,int gap);

