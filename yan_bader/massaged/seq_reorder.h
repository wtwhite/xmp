#include "interface.h"

/*compute lower_bound2[i]: the number of character state
present in taxa s_i but not appear before s_i*/
int state_discrepancy();

/*From the distance matrix, choose the pair with maximum distance*/
int max_distance(int *dist_matrix, int *max_i, int *max_j);

/*Combine any 2 sites and encode the combination, shape of
comb_seq:[num_seq][C(seq_len,2)+seq_len], for each taxon, the sites 
sequence 
is:(1,1),(1,2),(1,seq_len),(2,2),...(2,se_len),...,(seq_len.seq_len)*/
int encode_2_sites(T_char_state *comb_seq);

/*obtain the weight for each pair of sites*/
int obtain_weight(T_char_state *or1, int row_size, int level, int *w);


/*get the maximum matching for graphon n vertices, with weight w*/
int maximum_matching(int *w, int *total_w,int n,int *combine);


int maximum_matching2(int *w, int n,int *combine);


/*compute the lowerbond by maximum matching*/
int compute_lb(T_char_state *or_state,int *w,int *includei,int *score2);


/*Decide the best core tree that has the maximum tree cost,        
order0 ,order1,order2*/
int decide_core_tree(T_char_state *temp_state,int *orde0,int 
*order1,int *order2);

