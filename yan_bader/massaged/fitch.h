
#include "interface.h"

/*Compute the Farris interval for a given tree and return the parsimony  
 tree length*/
int compute_farris_int(tNode_T tree,int tree_size,int *tree_len);


/*compute the intersection of the two character state sets,maintaining the
increasing order of list,Assume s1,s2 are both in increasing order
int intersect(T_char_state_s *s1, T_char_state_s *s2,T_char_state_s *result);

/*Combine the union of the two non-intersected character state sets,
maintaining the increase order,assume s1,s2 are both in increasing order
int union_set(T_char_state_s *s1, T_char_state_s *s2,T_char_state_s *result);
*/
