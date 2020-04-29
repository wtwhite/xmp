#define STAT

#define MAX_COST  200000 /*the maximum possible score of a tree*/

#define NUM_PROCESSOR  1 /*The current version is sequential*/

#define POOL_SIZE 200 /*size of tree pool*/

#define  TRUE  1
#define  FALSE 0


#define RED  0
#define BLACK 1

typedef unsigned char T_char_state;  /*binary ebcoding of 
character states,
each states occupy one bit, 1 menas it has thatstates, in future can the 
type can be developed into multiple bytes*/

typedef  T_char_state * tNode_T;/*maximum possible n is 255,array
of size (2n-2), preorder traverse of a true tree: for terminals
: tag>0; internal node: tag=0*/

struct tree{
	struct tree *parent;
	struct tree *lChild;
	struct tree *rChild;
	unsigned char tag;
};


typedef struct sol_list_t{
	tNode_T  tree; /*array of size (2n-2), it is pre-order traverse of the
	true tree*/
	struct sol_list_t  *next;
}SOL_LIST_T;


struct hNode{
	int lowerbound; /*lowerbound of the trees consistent with the partial
	tree*/
	int upperbound;
	tNode_T par_sol;  /*per-order traaversal of the unrooted tree*/
	T_char_state  *internal_states;/*size :(num_seq-1)*seq_len*/ 
};

typedef struct fNode {
	int tree_size; /*at level i, it is 2*i-2*, i=3..n*/
	int array_size;
	int heap_size;
	struct hNode *heap;
} fNode_t;


extern int *order;
extern T_char_state *sum_state; /* shape:[num_seq*seq_len]*/ 
extern int *infor_len;
extern int *infor_site;
extern int global_lb;
extern T_char_state *seq_matrix; /*shape: seq_matrix[2*num_seq-2][seq_len],
the first "num_seq" rows are states of input, the remaining rows are states of the 
inernal nods whose tag >=num_seq*/
extern int num_seq;
extern int seq_len;
extern int thresh_level; /*From 3 to thresh_level, all the nodes are
generated*/
extern struct fNode *frontier; /*active nodes in frontier*/
extern int global_ub; /*Global upperbound*/
extern SOL_LIST_T  *opt_list; /*The list of best solution, phylogenies
with global_ub*/
extern int *lower_bound2; /*lowerbound of each partial tree =score of the
partial tree + lower_bound2 at this level. size:(n-3+1)*/
extern int *upper_bound2;/*upperbound of each partial tree =score of the
partial tree + upper_bound2 at this level. size:(n-3+1)*/

extern int sum_sinlegton_site; /*score of sinlegton sites*/

extern int *subtree_len; /*size: 2*num_seq-2. for leaves including the pseudo root,
it is 0, for internal node(with tag>=num_seq), it is length of the subtree with root
at the internal_node*/

typedef struct{
  int parent;
  int brother;
} T_tree_node;

typedef struct{
	int *p; /*size :2*num_seq, parent of each node*/
 	int *b;/*size :2*num_seq, brother of each node*/
	int *stack; /*size:n. defined for computing tree lenght*/
	T_char_state *temp1_state; /*size :seq_len*/
        T_char_state *temp2_state; /*size :seq_len*/
	T_tree_node *temp_true_tree; /*size: 2*num_seq-2*/
	tNode_T	tree1;
	tNode_T tree2;
	tNode_T  *tree_pool; /*to save the tree in opt_list; size:POOL_SIZE*/
	int tree_pool_ptr; /*point to next avaliable unit in tree_pool*/
	struct tree *true_tree_pool; /*array of struct tree, size :2*n-1 for rooted binary tree*/
	int *dist_matrix; /*distance matrix*/
} temp_memory_T;

extern temp_memory_T temp_memory;


#ifdef  STAT
typedef struct {
	int num_of_decomposed; /*the number of decomposed nodes in the B&B
	tree*/

}stat_data_T;

extern stat_data_T stat_data; /*the statistic data to analyze the
performance of the B&B*/
#endif
