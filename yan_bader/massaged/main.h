#include "interface.h"

/*process the input options*/
int   process_opt(int argc,char **argv,char **in_file,int *int_gb);


/*get the input information from a plain file, every row of the input file
is a sequence of some species except the first row, the first row of the
input file is format information: number of species, length of sequence*/
int get_input(char *filename,int **num_state);

/*Output the tree in a list*/
int output(SOL_LIST_T *list, int num_seq,int *order);

/*print the subtree with root at tree[i] in NEXUS form and return the next avaliable
position in tree*/
int print_tree(tNode_T tree, int i,int *order);


/*output the trees with minimum score,save them in opt_list*/
int extract_min(struct hNode *heap,int heap_size,int self_id);
