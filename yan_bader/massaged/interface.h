/*this include file contains the protype of functions for
/*cross-reference, they are excerpted from other inculde file such as
fitch.h*/
#include <assert.h>
#include "structs.h"

#define DEBUG
#undef INC_TREE_LEN  /*define incrementally compute tree lengthh*/

/*Use  Fitch's method to find the possible minimal state assignment for a
given tree and return the  parsimony tree length*/
int fitch_pass(tNode_T tree,int tree_size);

/*Use parsimoy informative site*/
int new_fitch_pass(tNode_T tree,int level);

/*Use parsimoy informative site,not only return cost, but alsoreutnr p,b*/
int new_fitch_pass_2(tNode_T tree,int level,int*p,int *b);

/*not only return cost, but alsoreutnr p,b*/
int fitch_pass_2(tNode_T tree,int tree_size,int*p,int *b);

/*called by static_reorder,process tree use full seq_len,return 
tree_len*/
int get_branch_interval(tNode_T tree,int tree_size);


/*fitch_pass on the previous level tree must be proceedey before it*/
int fitch_incremental_pass(T_tree_node *tree,int pos,int new_taxa_id);

int old_fitch_incremental_pass(T_tree_node *tree,int pos,int 
new_taxa_id);

/*Use Eck and Dayhoff(1966)'s heuristic to find a feasible
solution: add one specis at a time and it is inserted to the best
position, it is called "Wagner tree"?*/
int heuristic_Eck_Day();

/*Decide the order of sequence addition staticly by max-mimi rules, order is arrya of
size num_seq*/
int static_reorder(int *order);

/*Decide the order of sequence addition staticly by max-mimi rules*/
int static_reorder_max_mini(int *order);

/*Decide the order of sequence addition staticly by max-variance rules*/
int static_reorder_max_variance(int *order);

/*Decide the order of sequence addition staticly by max-mini_lb 
rules*/
int static_reorder_max_mini_lb(int *order);

/*Decide the order of sequence addition staticly by max-mini_lb 
rules, use constant Fitch algorithm*/
int static_reorder_max_mini_lb_const(int *order);


/*Decide the order of sequence addition staticly by max-rank rules*/
int static_reorder_max_rank(int *order);

int static_reorder_max_rank1(int *order);

int static_reorder_max_rank1_lb(int *order);

int static_reorder_mini_rank1_lb(int *order);

int static_reorder_rank_tournament(int *order);

/*Compute lower_bpind2[i]: the number of character state present in taxa
{s_i,..., s_n}, but no appear before*/
int lb_part2(int *order);

int lb2_part2(int *order);

/*Reorder column of sequence_matrix, remove constant sites and 
parsimonous-uninformative sites*/
int reorder_cite(int *num_state);

int reorder_cite1();

/*Intialize the frointer*/
int initalize_search(int int_gb);

/*Remove the minimal item heap[0] from the heap*/
int dequeue(int id);

/*Expand par_sol: inserting all of its children to the frontier[id] if
they pass the lowerbound test*/
int expand(struct hNode *par_sol, int id);

/*Expand par_sol: inserting all of its children to the frontier[id] if
they pass the lowerbound test, use the amortized constant algorithm to
compute the cost of new tree*/
int expand_cons_fit(struct hNode *sel_hnode, int id);


/*Compute upper_bound2: at level i, it is (the cost of tree on the
remaining taxa generated by Neighbor_Joining) + (the minimum distance
between one taxa from the partial tree and one taxa from the remaininh
taxa*/
int ub_part2();


/*compute the distance matrix*/
int  compute_dist_matrix(int *dist_matrix);

/*rearrange the seq_matrix according toorder*/
int rearrange_seq(int *order);


/*Decide the order of sequence addition staticly by max-mini_lb rules*/
int static_reorder_max_mini_lb_cross(int *order);

int old_static_reorder_max_mini_lb_cross(int *order);

/*Decide the order of sequence addition staticly by max lowerbound of 
the partial tree*/
int static_reorder_max_lb(int *order);

/*Decide the order of sequence addition staticly by max-mini_lb rules
use constant Fitch algorithm*/
int static_reorder_max_mini_lb_cons(int *order);

/* Greedy algorithm, each time insert the taxa at the best position*/
int greedy_fixed_order(tNode_T tree);

/*compute tree length, save internal node states, only consider for those
parsimonious informative sites*/
int new_fitch_pass_3(tNode_T tree,int level);

