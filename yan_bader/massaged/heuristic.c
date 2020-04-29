#include <stdlib.h>
#include "heuristic.h"

/* Greedy algorithm, each time insert the taxa at the best position*/
int greedy_fixed_order(tNode_T tree)
{ 
	tNode_T temp_tree;
	int i,j,l;
	int cost,mini_cost,best_pos;
	int tree_size,tree_len;
	int internal_id;
	T_char_state temp_state;
	T_char_state *ptr_seq1,*ptr_seq2;

	temp_tree=temp_memory.tree1;
	internal_id=num_seq;
	tree[0]=order[0];
	tree[1]=internal_id;
	tree[2]=order[1];
	tree[3]=order[2];


	internal_id++;
	/*Decide which position to be added at the i_th step*/
	for(i=3;i<num_seq;i++){
		tree_size=2*i-2;
		mini_cost=MAX_COST;

		tree_len=get_branch_interval(tree,tree_size);
		ptr_seq1=seq_matrix+order[i]*seq_len;
		/*try every position j, insert (internal_id,i) right before 
tree[j]*/
		for(j=1;j<tree_size;j++){
		 	cost=0;
			ptr_seq2=seq_matrix+(tree[j]+2*num_seq)*seq_len;
			for(l=0;l<seq_len;l++){
                        temp_state=ptr_seq1[l] & ptr_seq2[l];
                        if (temp_state==0) cost++;
                       }

			if(cost<mini_cost){
				mini_cost=cost;
				best_pos=j;
			}
		}

		mini_cost+=tree_len;

		/*build the partial tree, insert (internal_id,order[i]) at 
the branch
		(best_pos,best_pos->parent)*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=order[i];
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		internal_id++;
	}

	return(mini_cost);
}



/*Use Eck and Dayhoff(1966)'s heuristic to find a feasible
solution: add one specis at a time amd it is inserted to the best
position, it is called "Wagner tree"?;return a tree and its parsimony
tree length.*/
/*mini-mini algorithm*/
int heuristic_Eck_Day()
{ 
	int *order;

	int temp;
	int i,j,k,l;
	T_char_state mask;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*min_pos;  /*mininum cost and minimum position for each taxa*/
	int *ptr_cos,*ptr_pos;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	T_char_state temp_state1;
	int tree_len,delta;
	T_char_state *ptr_seq1,*ptr_seq2,*ptr_seq3;


	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	min_pos=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	
	order=(int *)malloc((2*num_seq-2)*sizeof(int));
for(i=0;i<2*num_seq-2;i++) order[i]=i;


	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	decide_core_tree(temp_state,&order0,&order1,&order2);

	/*reorder the sequence matrix*/
	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;


	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;


	temp_order=order[order2];
	order[order2]=order[2];
	order[2]=temp_order;

	/*initialize the first core tree*/
	tree[0]=order[0];
	tree[1]=internal_id;
	tree[2]=order[1];
	tree[3]=order[2];
	tree_size=4;

	internal_id++;

	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<=num_seq;i++){
		/*Preprocess tree to get branch inerval*/
		tree_len=get_branch_interval(tree,tree_size);

	        /*for each remaning sequences (s_i,....s_i), compute their best
 position (minimal cost) to insert into the current partial tree*/
		for(j=i;j<=num_seq;j++){
		  ptr_seq1=seq_matrix+order[j-1]*seq_len;
		  for(k=1;k<tree_size;k++){
		    /*insert s_j before k*/
		    delta=0;
                    ptr_seq2=seq_matrix+(tree[k]+2*num_seq)*seq_len;
		    for(l=0;l<seq_len;l++){
			temp_state1=ptr_seq1[l] & ptr_seq2[l];
			if (temp_state1==0) delta++;
		    }
		    *(min_cost+(j-1)*(2*num_seq-2)+k)=tree_len+delta;
		  }
		}

		/*for each taxa, compute the position of its minimum cost 
					and the median cost*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);

			ptr_pos=min_pos+j-1;
			*ptr_pos=1;  /*best position for taxa j*/
			for(k=2;k<tree_size;k++)
				if(ptr_cos[*ptr_pos]>ptr_cos[k])
					*ptr_pos=k;
		}

		/*find the minimum mini cost*/
		best_taxa=i;
		max_mini=*(min_cost+(i-1)*(2*num_seq-2)+min_pos[i-1]);
		for(j=i+1;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			if(max_mini > ptr_cos[min_pos[j-1]]) {
				best_taxa=j;
				max_mini=ptr_cos[min_pos[j-1]];
			}
		}

		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
		}

		best_pos=min_pos[best_taxa-1];
/*		
printf("i:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);*/

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;

		/*build the partial tree, insert the best_taxa at the branch
							(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
							tree[best_pos]*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=order[i-1];
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		tree_size+=2;
		internal_id++;
	}


	free(min_cost);
	free(min_pos);
	free(temp_state);
	return(max_mini);
 
}


