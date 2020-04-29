/* This file contains ths fucntions to reorder the addition ofsequences*/
// WTJW 28/4/2004: Following #include line added.
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "seq_reorder.h"

/*rearrange the seq_matrix according toorder*/
int rearrange_seq(int *order)
{
	T_char_state *out_seq,*ptr_seq1,*ptr_seq2;
	int i;

	out_seq=(T_char_state *) malloc(2*num_seq*seq_len*sizeof(T_char_state));

	for(i=0;i<num_seq;i++){
		ptr_seq1=seq_matrix+(order[i]-1)*seq_len;
		ptr_seq2=out_seq+i*seq_len;
		memcpy(ptr_seq2,ptr_seq1,seq_len*sizeof(T_char_state));
	}


	free(seq_matrix);
	seq_matrix=out_seq;

	return(0);
}


/*Given a sequence matrix,  return the pair with maximum distance
"(max_i,max_j)"*/
int compute_dist_matrix(int *dist_matrix)
{
	int i,j,k;
	int *ptr;
	T_char_state *ptr_seq;

	/*Initialize dist_matrix*/
	ptr=dist_matrix;
	for(i=0;i<num_seq;i++)
		for(j=0;j<num_seq;j++,ptr++)
			*ptr=0;

	/*Compute distance matricx, actually half of the distance matrix is
			enough*/
	ptr_seq=seq_matrix;
	for(k=0;k<seq_len;k++){
		for(i=0;i<num_seq;i++)
			for(j=0;j<num_seq;j++)
				if (*(ptr_seq+i*seq_len+k) != *(ptr_seq+j*seq_len+k))
					*(dist_matrix+i*num_seq+j)=*(dist_matrix+i*num_seq+j)+1;
	}


	return(0);
}


/*Decide the order of sequence addition staticly by max-mimi rules*/
int static_reorder(int *order)
{
	int order0,order1,temp_order;
	int i,j,k;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*min_pos;  /*mininum cost and minimum position for each taxa*/
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;

	min_cost=(int*)malloc(sizeof(int)*num_seq);
	min_pos=(int*)malloc(sizeof(int)*num_seq);

	/*Get the pair with the maximum distance*/
	compute_dist_matrix(temp_memory.dist_matrix);

	/*choose the pair with maximum distance**/
	max_distance(temp_memory.dist_matrix,&order0,&order1);

	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	/*Reorder seq_matrix, let the sequnce indexed by order0 to be the
			first sequence of seq_matrix*/
	seq_ptr1=seq_matrix;
	seq_ptr2=seq_matrix+order0*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));

	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;

	/*Let the sequence indexed by order1 to be the second sequence of
			seq_matrix */
	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;
	seq_ptr1=seq_matrix+1*seq_len;
	seq_ptr2=seq_matrix+order1*seq_len;

	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));


	/*Intialize the partial tree by the first 2 sequences*/
	tree=temp_memory.tree1;
	tree[0]=0;
	tree[1]=1;
	tree_size=2;
	temp_tree=temp_memory.tree2;

	internal_id=num_seq;
	/*Decide which sequence to be added at the i_th step*/
	for(i=3;i<num_seq;i++){

		for(j=i;j<=num_seq;j++) min_cost[j-1]=MAX_COST;

		/*for each remaning sequences (s_i,....s_i), compute their best
						position (minimal cost) to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=j-1;
				cur_cost=fitch_pass(temp_tree,2*i-2);
				if(cur_cost<min_cost[j-1]){
					min_cost[j-1]=cur_cost;
					min_pos[j-1]=k;
				}
			}
		}

		best_taxa=i;
		max_mini=0;

		/*Update the best_taxa and best position*/
		for(j=i;j<=num_seq;j++)
			if(min_cost[j-1]>max_mini){
				max_mini=min_cost[j-1];
				best_taxa=j;
			}
		best_pos=min_pos[best_taxa-1];


		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr1=seq_matrix+(i-1)*seq_len;
		seq_ptr2=seq_matrix+(best_taxa-1)*seq_len;

		memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));



		/*build the partial tree, insert the best_taxa at the branch
						(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
						tree[best_pos]*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=i-1;
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		tree_size+=2;
		internal_id++;
	}


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(temp_state);
	return(0);
}



/*Decide the order of sequence addition staticly by max-mimi rules*/
int static_reorder_max_mini(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*min_pos;  /*mininum cost and minimum position for each taxa*/
	int *ptr_cos,*ptr_pos,*med,*var;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;

	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	min_pos=(int*)malloc(sizeof(int)*num_seq);
	med=(int*)malloc(sizeof(int)*num_seq);
	var=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

	/*reorder the sequence matrix*/
	seq_ptr1=seq_matrix;
	seq_ptr2=seq_matrix+order0*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;


	seq_ptr1=seq_matrix+seq_len;
	seq_ptr2=seq_matrix+order1*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;


	seq_ptr1=seq_matrix+2*seq_len;
	seq_ptr2=seq_matrix+order2*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order2];
	order[order2]=order[2];
	order[2]=temp_order;

	/*initialize the first core tree*/
	tree[0]=0;
	tree[1]=internal_id;
	tree[2]=1;
	tree[3]=2;
	tree_size=4;

	internal_id++;
	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){


		/*for each remaning sequences (s_i,....s_i), compute their best
						position (minimal cost) to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=j-1;
				*(min_cost+(j-1)*(2*num_seq-2)+k)
				    =fitch_pass(temp_tree,2*i-2);
			}
		}

		/*for each taxa, compute the position of its minimum cost 
				and the median cost*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);

			/*compute median of the cost*/
			med[j-1]=0;
			for(k=1;k<tree_size;k++)
				med[j-1]+=ptr_cos[k];
			med[j-1]=med[j-1]/(tree_size-1);

			/*compute variance of the cost*/
			var[j-1]=0;
			for(k=1;k<tree_size;k++)
				var[j-1]+=(ptr_cos[k]-med[j-1])^2;

			ptr_pos=min_pos+j-1;
			*ptr_pos=1;  /*best position for taxa j*/
			for(k=2;k<tree_size;k++)
				if(ptr_cos[*ptr_pos]>ptr_cos[k])
					*ptr_pos=k;
		}

		/*find the maximum mini cost*/
		best_taxa=i;
		max_mini=*(min_cost+(i-1)*(2*num_seq-2)+min_pos[i-1]);
		for(j=i+1;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			if(max_mini < ptr_cos[min_pos[j-1]]) {
				best_taxa=j;
				max_mini=ptr_cos[min_pos[j-1]];
			}
		}

		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			if(max_mini == ptr_cos[min_pos[j-1]]) {
				if(var[j-1]>var[best_taxa-1])  best_taxa=j;
			}
		}

		best_pos=min_pos[best_taxa-1];
#ifdef DEBUG1
		printf("i:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);
#endif
		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr1=seq_matrix+(i-1)*seq_len;
		seq_ptr2=seq_matrix+(best_taxa-1)*seq_len;

		memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));



		/*build the partial tree, insert the best_taxa at the branch
						(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
						tree[best_pos]*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=i-1;
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		tree_size+=2;
		internal_id++;
	}


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(med);
	free(var);
	free(min_cost);
	free(min_pos);
	free(temp_state);
	return(0);
}




/*Decide the order of sequence addition staticly by max-rank rules*/
int static_reorder_max_rank(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*min_pos;  /*mininum cost and minimum position for each taxa*/
	int *ptr_cos,*ptr_pos,*med,*var;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	int *rank,*ptr_rank;


	rank=(int *)malloc(sizeof(int)*num_seq);
	ptr_rank=(int *)malloc(sizeof(int)*num_seq);

	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	min_pos=(int*)malloc(sizeof(int)*num_seq);
	med=(int*)malloc(sizeof(int)*num_seq);
	var=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

	/*reorder the sequence matrix*/
	seq_ptr1=seq_matrix;
	seq_ptr2=seq_matrix+order0*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;


	seq_ptr1=seq_matrix+seq_len;
	seq_ptr2=seq_matrix+order1*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;


	seq_ptr1=seq_matrix+2*seq_len;
	seq_ptr2=seq_matrix+order2*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order2];
	order[order2]=order[2];
	order[2]=temp_order;

	/*initialize the first core tree*/
	tree[0]=0;
	tree[1]=internal_id;
	tree[2]=1;
	tree[3]=2;
	tree_size=4;

	internal_id++;
	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){


		/*for each remaning sequences (s_i,....s_i), compute their best
						position (minimal cost) to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=j-1;
				*(min_cost+(j-1)*(2*num_seq-2)+k)
				    =fitch_pass(temp_tree,2*i-2);
			}
		}

		/*for each taxon, sort the cost in nondecresasing order 
				min_cost, and save the bese postion in min_pos*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			/*sort ptr_cos[1..(tree_size-1)] by insertion sorting*/
			ptr_cos[0]=-1;
			min_pos[j-1]=1;
			for(k=2;k<tree_size;k++){
				max_mini=ptr_cos[k];
				for(l=k-1;l>=0;l--){
					if(ptr_cos[l]>max_mini)  ptr_cos[l+1]=ptr_cos[l];
					else break;
				}
				ptr_cos[l+1]=max_mini;
				if(l==0) min_pos[j-1]=k;
			}
		}

		/*compute rank for each taxon*/
		for(j=i;j<=num_seq;j++) {
			rank[j-1]=0;
			ptr_rank[j-1]=1;
		}


		order0=0;
		for(k=0;k<(tree_size-1)*(num_seq-i+1);k++){
			max_mini=MAX_COST;
			for(l=i;l<=num_seq;l++){
				if(ptr_rank[l-1]<tree_size){
					temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
					if(max_mini >temp){
						max_mini=temp;
					}
				}
			}

			for(l=i;l<=num_seq;l++){
				if(ptr_rank[l-1]<tree_size){
					temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
					while(max_mini ==temp){
						rank[l-1]+=order0;
						ptr_rank[l-1]++;
						if(ptr_rank[l-1]<tree_size)
							temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
						else break;
					}
				}
			}

			if(max_mini>global_ub)
			{
				for(l=i;l<=num_seq;l++)
					rank[l-1]+=(tree_size-ptr_rank[l-1])*order0;
				break;
			}
			else order0++;

		}

		/*find the taxa with maximum  rank*/
		best_taxa=i;
		max_mini=rank[i-1];
		for(j=i+1;j<=num_seq;j++){
			if(max_mini < rank[j-1]) {
				best_taxa=j;
				max_mini=rank[j-1];
			}
		}


		best_pos=min_pos[best_taxa-1];
		printf("\ni:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr1=seq_matrix+(i-1)*seq_len;
		seq_ptr2=seq_matrix+(best_taxa-1)*seq_len;

		memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));



		/*build the partial tree, insert the best_taxa at the branch
						(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
						tree[best_pos]*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=i-1;
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		tree_size+=2;
		internal_id++;
	}


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);			// WTJW: was "Printf()", unbelievably

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(med);
	free(var);
	free(min_cost);
	free(min_pos);
	free(temp_state);
	return(0);
}



int static_reorder_max_rank1(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*max_pos;  /*mininum cost and maximum position for 
		each taxa*/
	int *ptr_cos,*ptr_pos,*med,*var;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	int *rank,*ptr_rank;


	rank=(int *)malloc(sizeof(int)*num_seq);
	ptr_rank=(int *)malloc(sizeof(int)*num_seq);

	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	max_pos=(int*)malloc(sizeof(int)*num_seq);
	med=(int*)malloc(sizeof(int)*num_seq);
	var=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

	/*reorder the sequence matrix*/
	seq_ptr1=seq_matrix;
	seq_ptr2=seq_matrix+order0*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;


	seq_ptr1=seq_matrix+seq_len;
	seq_ptr2=seq_matrix+order1*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;


	seq_ptr1=seq_matrix+2*seq_len;
	seq_ptr2=seq_matrix+order2*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order2];
	order[order2]=order[2];
	order[2]=temp_order;

	/*initialize the first core tree*/
	tree[0]=0;
	tree[1]=internal_id;
	tree[2]=1;
	tree[3]=2;
	tree_size=4;

	internal_id++;
	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){


		/*for each remaning sequences (s_i,....s_i), compute their best
						position (minimal cost) to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=j-1;
				*(min_cost+(j-1)*(2*num_seq-2)+k)
				    =fitch_pass(temp_tree,2*i-2);
			}
		}

		/*Save the worst position in max_pos*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			max_pos[j-1]=1;
			temp=ptr_cos[1];
			for(k=2;k<tree_size;k++){
				if(ptr_cos[k]<temp) {
					temp=ptr_cos[k];
					max_pos[j-1]=k;
				}
			}
		}

		/*for each taxon, sort the cost in nondecresasing order 
				min_cost*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			/*sort ptr_cos[1..(tree_size-1)] by insertion sorting*/
			ptr_cos[0]=-1;
			for(k=2;k<tree_size;k++){
				max_mini=ptr_cos[k];
				for(l=k-1;l>=0;l--){
					if(ptr_cos[l]>max_mini)  ptr_cos[l+1]=ptr_cos[l];
					else break;
				}
				ptr_cos[l+1]=max_mini;
			}
		}

		/*compute rank for each taxon*/
		for(j=i;j<=num_seq;j++) {
			rank[j-1]=0;
			ptr_rank[j-1]=1;
		}


		order0=0;
		for(k=0;k<(tree_size-1)*(num_seq-i+1);k++){
			max_mini=MAX_COST;
			for(l=i;l<=num_seq;l++){
				if(ptr_rank[l-1]<tree_size){
					temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
					if(max_mini >temp){
						max_mini=temp;
					}
				}
			}

			if(max_mini > global_ub) break;
			else order0++;

			for(l=i;l<=num_seq;l++){
				if(ptr_rank[l-1]<tree_size){
					temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
					while(max_mini ==temp){
						rank[l-1]+=order0;
						ptr_rank[l-1]++;
						if(ptr_rank[l-1]<tree_size)
							temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
						else break;
					}
				}
			}
		}


		/*find the taxa with minimum  ptr_rank, which means most 
				branches will be pruned*/
		best_taxa=i;
		temp=ptr_rank[i-1];
		for(j=i+1;j<=num_seq;j++){
			if(temp > ptr_rank[j-1]) {
				temp=ptr_rank[j-1];
			}
		}

		/*Choose taxon with maximum rank from those whose 
				ptr_rank==temp*/
		max_mini=-1;
		for(j=i;j<=num_seq;j++)
			if ((ptr_rank[j-1]==temp) && (max_mini<rank[j-1]))
			{
				max_mini=rank[j-1];
				best_taxa=j;
			}


		best_pos=max_pos[best_taxa-1];
		printf("\ni:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr1=seq_matrix+(i-1)*seq_len;
		seq_ptr2=seq_matrix+(best_taxa-1)*seq_len;

		memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));



		/*build the partial tree, insert the best_taxa at the branch
						(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
						tree[best_pos]*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=i-1;
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		tree_size+=2;
		internal_id++;
	}


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(med);
	free(var);
	free(min_cost);
	free(max_pos);
	free(temp_state);
	return(0);
}



int static_reorder_max_rank1_lb(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state mask;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*max_pos;  /*mininum cost and maximum position for 
		each taxa*/
	int *ptr_cos,*ptr_pos,*med,*var;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	int *rank,*ptr_rank;
	T_char_state *all_flag,*cur_flag,*temp_flag;

	all_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	cur_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	temp_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);



	rank=(int *)malloc(sizeof(int)*num_seq);
	ptr_rank=(int *)malloc(sizeof(int)*num_seq);

	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	max_pos=(int*)malloc(sizeof(int)*num_seq);
	med=(int*)malloc(sizeof(int)*num_seq);
	var=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

	/*reorder the sequence matrix*/
	seq_ptr1=seq_matrix;
	seq_ptr2=seq_matrix+order0*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;


	seq_ptr1=seq_matrix+seq_len;
	seq_ptr2=seq_matrix+order1*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;


	seq_ptr1=seq_matrix+2*seq_len;
	seq_ptr2=seq_matrix+order2*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order2];
	order[order2]=order[2];
	order[2]=temp_order;

	/*initialize the first core tree*/
	tree[0]=0;
	tree[1]=internal_id;
	tree[2]=1;
	tree[3]=2;
	tree_size=4;

	internal_id++;

	for(i=0;i<seq_len;i++) cur_flag[i]=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<3;i++)
		for(j=0;j<seq_len;j++,seq_ptr1++)
			cur_flag[j]=cur_flag[j] | (*seq_ptr1);
	for(i=0;i<seq_len;i++) all_flag[i]=cur_flag[i];
	for(i=3;i<num_seq;i++)
		for(j=0;j<seq_len;j++,seq_ptr1++)
			all_flag[j]=all_flag[j] | (*seq_ptr1);

	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){
		/*for each taxon, compute the lb2*/
		for(j=i;j<=num_seq;j++){
			seq_ptr1=seq_matrix+(j-1)*seq_len;
			lower_bound2[j]=0;
			for(k=0;k<seq_len;k++){
				temp_flag[k]=cur_flag[k] | seq_ptr1[k];
				temp_flag[k]=all_flag[k]-temp_flag[k];
				for(l=0,mask=1;l<8*sizeof(T_char_state);l++,mask=mask<<1){
					if(temp_flag[k] & mask) lower_bound2[j]++;
				}
			}
		}


		/*for each remaning sequences (s_i,....s_i), compute their best
							position (minimal cost) to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=j-1;
				*(min_cost+(j-1)*(2*num_seq-2)+k)
				    =fitch_pass(temp_tree,2*i-2)
				+lower_bound2[j];
			}
		}

		/*Save the worst position in max_pos*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			max_pos[j-1]=1;
			temp=ptr_cos[1];
			for(k=2;k<tree_size;k++){
				if(ptr_cos[k]<temp) {
					temp=ptr_cos[k];
					max_pos[j-1]=k;
				}
			}
		}

		/*for each taxon, sort the cost in nondecresasing order 
					min_cost*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			/*sort ptr_cos[1..(tree_size-1)] by insertion sorting*/
			ptr_cos[0]=-1;
			for(k=2;k<tree_size;k++){
				max_mini=ptr_cos[k];
				for(l=k-1;l>=0;l--){
					if(ptr_cos[l]>max_mini)  ptr_cos[l+1]=ptr_cos[l];
					else break;
				}
				ptr_cos[l+1]=max_mini;
			}
		}

		/*compute rank for each taxon*/
		for(j=i;j<=num_seq;j++) {
			rank[j-1]=0;
			ptr_rank[j-1]=1;
		}


		order0=0;
		for(k=0;k<(tree_size-1)*(num_seq-i+1);k++){
			max_mini=MAX_COST;
			for(l=i;l<=num_seq;l++){
				if(ptr_rank[l-1]<tree_size){
					temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
					if(max_mini >temp){
						max_mini=temp;
					}
				}
			}

			if(max_mini > global_ub) break;
			else order0++;

			for(l=i;l<=num_seq;l++){
				if(ptr_rank[l-1]<tree_size){
					temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
					while(max_mini ==temp){
						rank[l-1]+=order0;
						ptr_rank[l-1]++;
						if(ptr_rank[l-1]<tree_size)
							temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
						else break;
					}
				}
			}
		}


		/*find the taxa with minimum  ptr_rank, which means most 
					branches will be pruned*/
		best_taxa=i;
		temp=ptr_rank[i-1];
		for(j=i+1;j<=num_seq;j++){
			if(temp > ptr_rank[j-1]) {
				temp=ptr_rank[j-1];
			}
		}

		/*Choose taxon with maximum rank from those whose 
					ptr_rank==temp*/
		max_mini=-1;
		for(j=i;j<=num_seq;j++)
			if ((ptr_rank[j-1]==temp) && (max_mini<rank[j-1]))
			{
				max_mini=rank[j-1];
				best_taxa=j;
			}


		best_pos=max_pos[best_taxa-1];
		printf("\ni:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr1=seq_matrix+(i-1)*seq_len;
		seq_ptr2=seq_matrix+(best_taxa-1)*seq_len;

		/*update cur_flag*/
		for(k=0;k<seq_len;k++)
			cur_flag[k]=cur_flag[k] | seq_ptr2[k];



		memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));



		/*build the partial tree, insert the best_taxa at the branch
							(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
							tree[best_pos]*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=i-1;
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		tree_size+=2;
		internal_id++;
	}


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(med);
	free(var);
	free(min_cost);
	free(max_pos);
	free(temp_state);
	return(0);
}



/*Decide the order of sequence addition staticly by max-mini_lb rules*/
/*This old version physically reorder the taxon (rows of seq_matrix)), which affect fitch_pass and output etc.*/
int old_static_reorder_max_mini_lb(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state mask;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*min_pos;  /*mininum cost and minimum position for each taxa*/
	int *ptr_cos,*ptr_pos,*med,*var;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	T_char_state *all_flag,*cur_flag,*temp_flag;

	all_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	cur_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	temp_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	min_pos=(int*)malloc(sizeof(int)*num_seq);
	med=(int*)malloc(sizeof(int)*num_seq);
	var=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

	/*reorder the sequence matrix*/
	seq_ptr1=seq_matrix;
	seq_ptr2=seq_matrix+order0*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;


	seq_ptr1=seq_matrix+seq_len;
	seq_ptr2=seq_matrix+order1*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;


	seq_ptr1=seq_matrix+2*seq_len;
	seq_ptr2=seq_matrix+order2*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
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

	for(i=0;i<seq_len;i++) cur_flag[i]=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<3;i++)
		for(j=0;j<seq_len;j++,seq_ptr1++)
			cur_flag[j]=cur_flag[j] | (*seq_ptr1);
	for(i=0;i<seq_len;i++) all_flag[i]=cur_flag[i];
	for(i=3;i<num_seq;i++)
		for(j=0;j<seq_len;j++,seq_ptr1++)
			all_flag[j]=all_flag[j] | (*seq_ptr1);

	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){

		/*for each taxon, compute the lb2*/
		for(j=i;j<=num_seq;j++){
			seq_ptr1=seq_matrix+(j-1)*seq_len;
			lower_bound2[j]=0;
			for(k=0;k<seq_len;k++){
				temp_flag[k]=cur_flag[k] | seq_ptr1[k];
				temp_flag[k]=all_flag[k]-temp_flag[k];
				for(l=0,mask=1;l<8*sizeof(T_char_state);l++,mask=mask<<1){
					if(temp_flag[k] & mask) lower_bound2[j]++;
				}
			}
		}

		/*for each remaning sequences (s_i,....s_i), compute their best
							position (minimal cost) to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,order[j-1]) right before 
tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=order[j-1];
				*(min_cost+(j-1)*(2*num_seq-2)+k)
				    =fitch_pass(temp_tree,2*i-2)+lower_bound2[j];
			}
		}

		/*for each taxa, compute the position of its minimum cost 
					and the median cost*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);

			/*compute median of the cost*/
			med[j-1]=0;
			for(k=1;k<tree_size;k++)
				med[j-1]+=ptr_cos[k];
			med[j-1]=med[j-1]/(tree_size-1);

			/*compute variance of the cost*/
			var[j-1]=0;
			for(k=1;k<tree_size;k++)
				var[j-1]+=(ptr_cos[k]-med[j-1])^2;

			ptr_pos=min_pos+j-1;
			*ptr_pos=1;  /*best position for taxa j*/
			for(k=2;k<tree_size;k++)
				if(ptr_cos[*ptr_pos]>ptr_cos[k])
					*ptr_pos=k;
		}

		/*find the maximum mini cost*/
		best_taxa=i;
		max_mini=*(min_cost+(i-1)*(2*num_seq-2)+min_pos[i-1]);
		for(j=i+1;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			if(max_mini < ptr_cos[min_pos[j-1]]) {
				best_taxa=j;
				max_mini=ptr_cos[min_pos[j-1]];
			}
		}

		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			if(max_mini == ptr_cos[min_pos[j-1]]) {
				if(var[j-1]>var[best_taxa-1])  best_taxa=j;
			}
		}

		best_pos=min_pos[best_taxa-1];
		printf("i:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr1=seq_matrix+(i-1)*seq_len;
		seq_ptr2=seq_matrix+(best_taxa-1)*seq_len;

		/*update cur_flag*/
		for(k=0;k<seq_len;k++)
			cur_flag[k]=cur_flag[k] | seq_ptr2[k];

		memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));



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


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(med);
	free(var);
	free(min_cost);
	free(min_pos);
	free(temp_state);
	return(0);
}




/*Decide the order of sequence addition staticly by max-variance rules*/
int static_reorder_max_variance(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*min_pos;  /*mininum cost and minimum position for each taxa*/
	int *ptr_cos,*ptr_pos,*med,*var;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;

	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	min_pos=(int*)malloc(sizeof(int)*num_seq);
	med=(int*)malloc(sizeof(int)*num_seq);
	var=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

	/*reorder the sequence matrix*/
	seq_ptr1=seq_matrix;
	seq_ptr2=seq_matrix+order0*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;


	seq_ptr1=seq_matrix+seq_len;
	seq_ptr2=seq_matrix+order1*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;


	seq_ptr1=seq_matrix+2*seq_len;
	seq_ptr2=seq_matrix+order2*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order2];
	order[order2]=order[2];
	order[2]=temp_order;

	/*initialize the first core tree*/
	tree[0]=0;
	tree[1]=internal_id;
	tree[2]=1;
	tree[3]=2;
	tree_size=4;

	internal_id++;
	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){


		/*for each remaning sequences (s_i,....s_i), compute their best
							position (minimal cost) to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=j-1;
				*(min_cost+(j-1)*(2*num_seq-2)+k)
				    =fitch_pass(temp_tree,2*i-2);
			}
		}

		/*for each taxa, compute the position of its minimum cost 
					and the median cost*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);

			/*compute median of the cost*/
			med[j-1]=0;
			for(k=1;k<tree_size;k++)
				med[j-1]+=ptr_cos[k];
			med[j-1]=med[j-1]/(tree_size-1);

			/*compute variance of the cost*/
			var[j-1]=0;
			for(k=1;k<tree_size;k++)
				var[j-1]+=(ptr_cos[k]-med[j-1])^2;

			ptr_pos=min_pos+j-1;
			*ptr_pos=1;  /*best position for taxa j*/
			for(k=2;k<tree_size;k++)
				if(ptr_cos[*ptr_pos]>ptr_cos[k])
					*ptr_pos=k;
		}

		/*find the maximum mini cost*/
		best_taxa=i;
		max_mini=var[i-1];
		for(j=i+1;j<=num_seq;j++){
			if(max_mini < var[j-1]) {
				best_taxa=j;
				max_mini=var[j-1];
			}
		}


		best_pos=min_pos[best_taxa-1];
		printf("i:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr1=seq_matrix+(i-1)*seq_len;
		seq_ptr2=seq_matrix+(best_taxa-1)*seq_len;

		memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));



		/*build the partial tree, insert the best_taxa at the branch
							(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
							tree[best_pos]*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=i-1;
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		tree_size+=2;
		internal_id++;
	}


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(med);
	free(var);
	free(min_cost);
	free(min_pos);
	free(temp_state);

	return(0);
}


/*compute lower_bound2[i]: the number of character state
	present in taxa s_i but not appear before s_i; num_state[j]:the number
	of state of the j_th character.*/
int state_discrepancy()
{
	int i,j,k;
	int absent;
	T_char_state *or_state;
	T_char_state *ptr_seq;

	or_state=(T_char_state *)malloc(seq_len*sizeof(T_char_state));
	for(j=0;j<seq_len;j++) or_state[j]=0;

	ptr_seq=seq_matrix;

	for(i=0;i<num_seq;i++){
		lower_bound2[i]=0;
		for(j=0;j<seq_len;j++,ptr_seq++){
			if ((or_state[j] & (*ptr_seq)) == 0) lower_bound2[i]++;
			or_state[j] =or_state[j] | (*ptr_seq);
		}
	}


	free(or_state);
	return(0);
}


/*Compute lower_bound2[i]: the number of character state present in taxa
	{s_i,..., s_n}, but no appear before*/
int lb_part2(int *order)
 {
	int i,j,k,score;
	T_char_state temp_state,mask;
	T_char_state *or_state,*ptr_seq;
	
	or_state=(T_char_state *)malloc(seq_len*sizeof(T_char_state));
	
	/*compute lower_bound2[i]: the # of states from s_1..s_i*/
	ptr_seq=seq_matrix+order[0]*seq_len;
	for(i=0;i<seq_len;i++)
		or_state[i]=ptr_seq[i];
	lower_bound2[1]=seq_len;
	
	for(i=2;i<=num_seq;i++){
	  ptr_seq=seq_matrix+order[i-1]*seq_len;
	  score=0;
	  for(j=0;j<seq_len;j++){
		or_state[j]=or_state[j] |ptr_seq[j];
		temp_state=or_state[j];
		mask=1;
		for(k=0;k<8*sizeof(T_char_state);k++,mask=mask<<1)
		  if ((mask & temp_state) !=0)  score++;
	  }
	  lower_bound2[i]=score;
	}
	

	for(i=1;i<num_seq;i++){
	  lower_bound2[i]=lower_bound2[num_seq]-lower_bound2[i];
	}
	lower_bound2[num_seq]=0;

	free(or_state);

	return(0);
}



/*Reorder column of sequence_matrix, remove constant sites and 
	parsimonous-uninformative sites*/
int reorder_cite(int *num_state)
{
	T_char_state *out_seq_matrix;
	T_char_state *ptr_in_seq,*ptr_out_seq;
	int i,j,max_state,temp,temp_state;
	int new_seq_len;
	int *order;

	order=(int *)malloc(sizeof(int)*seq_len);
	for(i=0;i<seq_len;i++) order[i]=i;

	/*Get the sorted order of num_state: bubble sort*/
	for(i=0;i<seq_len;i++){
		max_state=num_state[i];
		for(j=i+1;j<seq_len;j++) {
			if(max_state<num_state[j]){
				num_state[i]=num_state[j];
				num_state[j]=max_state;
				max_state=num_state[i];

				temp=order[i];
				order[i]=order[j];
				order[j]=temp;
			}
		}
	}


	/*Remove constant sites*/
	new_seq_len=seq_len;
	for(i=seq_len-1;(i>=0) && (num_state[i]==1);i--,new_seq_len--);

	/*Permute the columns of sequence matrix*/
	out_seq_matrix=(T_char_state *)malloc(sizeof(T_char_state)*4*num_seq*seq_len);
	ptr_in_seq=seq_matrix;
	ptr_out_seq=out_seq_matrix;
	for(j=0;j<new_seq_len;j++)
		for(i=0;i<num_seq;i++)
			*(ptr_out_seq+i*new_seq_len+j)=*(ptr_in_seq+i*seq_len+order[j]);

	free(seq_matrix);
	seq_matrix=out_seq_matrix;
	seq_len=new_seq_len;
	printf("seq_len:%d\n",seq_len);
	free(order);

#ifdef DEBUG1
	printf("\n the sequence matrix after reordering cite\n");
	ptr_out_seq=out_seq_matrix;
	for(i=0;i<num_seq;i++) {
		for(j=0;j<seq_len;j++,ptr_out_seq++)
			printf("%d  ",*ptr_out_seq);
		printf("\n");
	}
#endif

	return(0);
}


/*From the distance matrix, choose the pair with maximum distance*/
int max_distance(int *dist_matrix, int *max_i, int *max_j)
{
	int max_dist,dist;
	int i,j;

	max_dist=0;
	*max_i=0;
	*max_j=1;
	for(i=0;i<num_seq;i++)
		for(j=i+1;j<num_seq;j++){
			dist=*(dist_matrix+i*num_seq+j);
			if(max_dist < dist)
			{
				max_dist=dist;
				*max_i=i;
				*max_j=j;
			}
		}

	return(0);
}

int static_reorder_mini_rank1_lb(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state mask;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*max_pos;  /*mininum cost and maximum position for 
		each taxa*/
	int *ptr_cos,*ptr_pos,*med,*var;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	int *rank,*ptr_rank;
	T_char_state *all_flag,*cur_flag,*temp_flag;

	all_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	cur_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	temp_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);



	rank=(int *)malloc(sizeof(int)*num_seq);
	ptr_rank=(int *)malloc(sizeof(int)*num_seq);

	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	max_pos=(int*)malloc(sizeof(int)*num_seq);
	med=(int*)malloc(sizeof(int)*num_seq);
	var=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

	/*reorder the sequence matrix*/
	seq_ptr1=seq_matrix;
	seq_ptr2=seq_matrix+order0*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;


	seq_ptr1=seq_matrix+seq_len;
	seq_ptr2=seq_matrix+order1*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;


	seq_ptr1=seq_matrix+2*seq_len;
	seq_ptr2=seq_matrix+order2*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order2];
	order[order2]=order[2];
	order[2]=temp_order;

	/*initialize the first core tree*/
	tree[0]=0;
	tree[1]=internal_id;
	tree[2]=1;
	tree[3]=2;
	tree_size=4;

	internal_id++;

	for(i=0;i<seq_len;i++) cur_flag[i]=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<3;i++)
		for(j=0;j<seq_len;j++,seq_ptr1++)
			cur_flag[j]=cur_flag[j] | (*seq_ptr1);
	for(i=0;i<seq_len;i++) all_flag[i]=cur_flag[i];
	for(i=3;i<num_seq;i++)
		for(j=0;j<seq_len;j++,seq_ptr1++)
			all_flag[j]=all_flag[j] | (*seq_ptr1);

	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){
		/*for each taxon, compute the lb2*/
		for(j=i;j<=num_seq;j++){
			seq_ptr1=seq_matrix+(j-1)*seq_len;
			lower_bound2[j]=0;
			for(k=0;k<seq_len;k++){
				temp_flag[k]=cur_flag[k] | seq_ptr1[k];
				temp_flag[k]=all_flag[k]-temp_flag[k];
				for(l=0,mask=1;l<sizeof(T_char_state);l++,mask=mask<<1){
					if(temp_flag[k] & mask) lower_bound2[j]++;
				}
			}
		}


		/*for each remaning sequences (s_i,....s_i), compute their best
							position (minimal cost) to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=j-1;
				*(min_cost+(j-1)*(2*num_seq-2)+k)
				    =fitch_pass(temp_tree,2*i-2);
				+lower_bound2[j];
			}
		}

		/*Save the best position in max_pos*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			max_pos[j-1]=1;
			temp=ptr_cos[1];
			for(k=2;k<tree_size;k++){
				if(ptr_cos[k]<temp) {
					temp=ptr_cos[k];
					max_pos[j-1]=k;
				}
			}
		}

		/*for each taxon, sort the cost in nondecresasing order 
					min_cost*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			/*sort ptr_cos[1..(tree_size-1)] by insertion sorting*/
			ptr_cos[0]=-1;
			for(k=2;k<tree_size;k++){
				max_mini=ptr_cos[k];
				for(l=k-1;l>=0;l--){
					if(ptr_cos[l]>max_mini)  ptr_cos[l+1]=ptr_cos[l];
					else break;
				}
				ptr_cos[l+1]=max_mini;
			}
		}

		/*compute rank for each taxon*/
		for(j=i;j<=num_seq;j++) {
			rank[j-1]=0;
			ptr_rank[j-1]=1;
		}


		order0=0;
		for(k=0;k<(tree_size-1)*(num_seq-i+1);k++){
			max_mini=MAX_COST;
			for(l=i;l<=num_seq;l++){
				if(ptr_rank[l-1]<tree_size){
					temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
					if(max_mini >temp){
						max_mini=temp;
					}
				}
			}

			if(max_mini > global_ub) break;
			else order0++;

			for(l=i;l<=num_seq;l++){
				if(ptr_rank[l-1]<tree_size){
					temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
					while(max_mini ==temp){
						rank[l-1]+=order0;
						ptr_rank[l-1]++;
						if(ptr_rank[l-1]<tree_size)
							temp=*(min_cost+(l-1)*(2*num_seq-2)+ptr_rank[l-1]);
						else break;
					}
				}
			}
		}


		/*find the taxa with minimum  ptr_rank, 
		which means most 
					branches will be pruned*/
		best_taxa=i;
		temp=ptr_rank[i-1];
		for(j=i+1;j<=num_seq;j++){
			if(temp > ptr_rank[j-1]) {
				temp=ptr_rank[j-1];
			}
		}

		/*Choose taxon with minimum rank from those whose 
					ptr_rank==temp*/
		max_mini=MAX_COST;
		for(j=i;j<=num_seq;j++)
			if ((ptr_rank[j-1]==temp) && 
			    (max_mini>rank[j-1]))
			{
				max_mini=rank[j-1];
				best_taxa=j;
			}


		best_pos=max_pos[best_taxa-1];
		printf("\ni:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr1=seq_matrix+(i-1)*seq_len;
		seq_ptr2=seq_matrix+(best_taxa-1)*seq_len;

		/*update cur_flag*/
		for(k=0;k<seq_len;k++)
			cur_flag[k]=cur_flag[k] | seq_ptr2[k];



		memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));



		/*build the partial tree, insert the best_taxa at the branch
							(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
							tree[best_pos]*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=i-1;
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		tree_size+=2;
		internal_id++;
	}


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(med);
	free(var);
	free(min_cost);
	free(max_pos);
	free(temp_state);
	return(0);
}


/*Set infor_site[seq_len] and infor_len[seq_len]*/
int reorder_cite1()
{
	T_char_state *out_seq_matrix,*or_state;
	T_char_state *ptr_in_seq,*ptr_out_seq;
	T_char_state  temp_state;
	int i,j,k,temp,sum_len;
	int *nonsingle_level,*cite_order;
	T_char_state *non_singleton_state;

	nonsingle_level=(int *)malloc(sizeof(int)*seq_len);
	for(i=0;i<seq_len;i++)
		nonsingle_level[i]=seq_len;

	cite_order=(int *)malloc(sizeof(int)*seq_len);

	or_state=(T_char_state 
	    *)malloc(sizeof(T_char_state)*seq_len);

	non_singleton_state=(T_char_state
	    *)malloc(sizeof(T_char_state)*seq_len);
	for(i=0;i<seq_len;i++)
		non_singleton_state[i]=0;


	/*get or_state from the first 4 taxon*/
	ptr_in_seq=seq_matrix+order[0]*seq_len;
	for(i=0;i<seq_len;i++,ptr_in_seq++) or_state[i]=*ptr_in_seq;
	for(i=1;i<4;i++){
		 ptr_in_seq=seq_matrix+order[i]*seq_len;
		for(j=0;j<seq_len;j++,ptr_in_seq++){
			temp_state=*ptr_in_seq;
			if (temp_state & or_state[j])
				non_singleton_state[j]|=temp_state;
			or_state[j]=or_state[j]|temp_state;
		}
	}

	/*for each site, check wether it is  parsimonious informative site*/
	temp=0;
	sum_len=0;
	for(i=0;i<seq_len;i++){
		temp_state=non_singleton_state[i];

		if (temp_state>0){
			while((temp_state &0x1)==0)
				temp_state=temp_state>>1;
			if(temp_state !=1) {
				nonsingle_level[i]=3;
				temp++;
			}
		}


		if(nonsingle_level[i]>3){
			/*for non-informative sites, find how many 
			states in temp_state*/
			temp_state=or_state[i];
			sum_len--;
			for(k=0;k<8*sizeof(T_char_state);k++,temp_state=temp_state>>1)
				if(temp_state & 0x1) sum_len++;
		}
	}
	infor_site[3]=temp;
	infor_len[3]=sum_len;

	/*for level i where to add the i_th taxon, find the new 
	 parsimonious informative sites*/
	for(i=4;i<num_seq;i++){
		ptr_in_seq=seq_matrix+order[i]*seq_len;
		for(j=0;j<seq_len;j++,ptr_in_seq++) {
			temp_state=*ptr_in_seq;
			if (temp_state & or_state[j])
				non_singleton_state[j]|=temp_state;
			or_state[j]=or_state[j] |  temp_state;
		}
		temp=0;
		sum_len=0;
		for(j=0;j<seq_len;j++){
			/*check wether it is new informative site*/
					if (nonsingle_level[j]>i){
					  temp_state=non_singleton_state[j];
					  if(temp_state>0){
			                   while((temp_state &0x1)==0)
			                     temp_state=temp_state>>1;   
			                   if(temp_state !=1) {
						nonsingle_level[j]=i;
						temp++;
					    }
					   }
				        }
			
/*for non-informative sites,
 find how many states in temp_state*/
					if(nonsingle_level[j]>i){
			temp_state=or_state[j];
			sum_len--;
			for(k=0;k<8*sizeof(T_char_state);k++,temp_state=temp_state>>1)
				if(temp_state & 0x1) sum_len++;
		}

	}
	infor_site[i]=infor_site[i-1]+temp;
	infor_len[i]=sum_len;
}


/*Get the cite_order of sites by nonsinleton_level*/
temp=0;
for(i=3;i<num_seq;i++)
for(j=0;j<seq_len;j++)
if(nonsingle_level[j]==i){
	cite_order[j]=temp;
	temp++;
}



/*Permute the columns of sequence matrix*/
out_seq_matrix=(T_char_state 
*)malloc(sizeof(T_char_state)*4*num_seq*seq_len);
ptr_in_seq=seq_matrix;
ptr_out_seq=out_seq_matrix;
for(j=0;j<seq_len;j++)
for(i=0;i<num_seq;i++)
*(ptr_out_seq+i*seq_len+cite_order[j])=*(ptr_in_seq+i*seq_len+j);

free(seq_matrix);
seq_matrix=out_seq_matrix;

free(nonsingle_level);
free(cite_order);
free(or_state);

#ifdef DEBUG1
printf("\n the sequence matrix after recite_ordering cite\n");
ptr_out_seq=out_seq_matrix;
for(i=0;i<num_seq;i++) {
	for(j=0;j<seq_len;j++,ptr_out_seq++)
		printf("%d  ",*ptr_out_seq);
	printf("\n");
}
#endif

return(0);
}



int static_reorder_rank_tournament(int *order)
{
	int temp;
	int i,j,k,l,m;
	T_char_state mask;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *cost,*min_pos;  /*cost at eachposition for each taxa*/
	int *ptr_cos1,*ptr_cos2,*ptr_pos;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	int *rank,*beat;
	T_char_state *all_flag,*cur_flag,*temp_flag;

	all_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	cur_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	temp_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);



	rank=(int *)malloc(sizeof(int)*num_seq);
	beat=(int*)malloc(sizeof(int)*num_seq*num_seq);

	cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	min_pos=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

	/*reorder the sequence matrix*/
	seq_ptr1=seq_matrix;
	seq_ptr2=seq_matrix+order0*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order0];
	order[order0]=order[0];
	order[0]=temp_order;


	seq_ptr1=seq_matrix+seq_len;
	seq_ptr2=seq_matrix+order1*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order1];
	order[order1]=order[1];
	order[1]=temp_order;


	seq_ptr1=seq_matrix+2*seq_len;
	seq_ptr2=seq_matrix+order2*seq_len;
	memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
	memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));
	temp_order=order[order2];
	order[order2]=order[2];
	order[2]=temp_order;

	/*initialize the first core tree*/
	tree[0]=0;
	tree[1]=internal_id;
	tree[2]=1;
	tree[3]=2;
	tree_size=4;

	internal_id++;

	for(i=0;i<seq_len;i++) cur_flag[i]=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<3;i++)
		for(j=0;j<seq_len;j++,seq_ptr1++)
			cur_flag[j]=cur_flag[j] | (*seq_ptr1);
	for(i=0;i<seq_len;i++) all_flag[i]=cur_flag[i];
	for(i=3;i<num_seq;i++)
		for(j=0;j<seq_len;j++,seq_ptr1++)
			all_flag[j]=all_flag[j] | (*seq_ptr1);

	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){
		/*for each taxon, compute the lb2*/
		for(j=i;j<=num_seq;j++){
			seq_ptr1=seq_matrix+(j-1)*seq_len;
			lower_bound2[j]=0;
			for(k=0;k<seq_len;k++){
				temp_flag[k]=cur_flag[k] | seq_ptr1[k];
				temp_flag[k]=all_flag[k]-temp_flag[k];
				for(l=0,mask=1;l<8*sizeof(T_char_state);l++,mask=mask<<1){
					if(temp_flag[k] & mask) lower_bound2[j]++;
				}
			}
		}


		/*for each remaning sequences (s_i,....s_i), compute costs to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=j-1;
				*(cost+(j-1)*(2*num_seq-2)+k)
				    =fitch_pass(temp_tree,2*i-2)
				+lower_bound2[j];
			}
		}

		/*Save the best position in min_pos*/
		for(j=i;j<=num_seq;j++){
			ptr_cos1=cost+(j-1)*(2*num_seq-2);
			min_pos[j-1]=1;
			temp=ptr_cos1[1];
			for(k=2;k<tree_size;k++){
				if(ptr_cos1[k]<temp) {
					temp=ptr_cos1[k];
					min_pos[j-1]=k;
				}
			}
		}

		/*compute beat[i,j], which is the # of times s_i beat s_j when comparing 
their lowerbounds at each position, beat[i,j]+beat[j,i]=2*(tree_size-1)*(tree_size-1)*/
		for(j=i;j<=num_seq;j++){
		   *(beat+(j-1)*num_seq+(j-1))=0; 
		   for(k=j+1;k<=num_seq;k++){
			temp=0; /*beat[j][k]=0*/
			ptr_cos1=cost+(j-1)*(2*num_seq-2);
			ptr_cos2=cost+(k-1)*(2*num_seq-2);

			/*compare taxons s_j and s_k for each pair of positions*/
			for(l=1;l<tree_size;l++)
			  for(m=1;m<tree_size;m++)
			    if (ptr_cos1[l]>ptr_cos2[m]) temp+=2;
			    else if (ptr_cos1[l]==ptr_cos2[m]) temp++;

			*(beat+(j-1)*num_seq+(k-1))=temp;
			*(beat+(k-1)*num_seq+(j-1))=2*(tree_size-1)*(tree_size-1)-temp;  
		   }
		}


		/*compute rank[i], how many taxons s_i beat*/
		temp=(tree_size-1)*(tree_size-1);
		for(j=i;j<=num_seq;j++){
		  rank[j-1]=0;
		  ptr_cos1=beat+(j-1)*num_seq;
		  for(k=i;k<=num_seq;k++)
		     rank[j-1]+=ptr_cos1[k-1];
/*
		    if(ptr_cos1[k-1]>temp) rank[j-1]++;
*/
		}


		/*Find the taxon with maximum rank*/
		best_taxa=i;
		temp=rank[i-1];
		for(j=i+1;j<=num_seq;j++){
		  if(temp<rank[j-1]){
			temp=rank[j-1];
			best_taxa=j;
		  }
		}

/*		best_pos=min_pos[best_taxa-1];
*/
best_pos=1;		printf("\ni:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr1=seq_matrix+(i-1)*seq_len;
		seq_ptr2=seq_matrix+(best_taxa-1)*seq_len;

		/*update cur_flag*/
		for(k=0;k<seq_len;k++)
			cur_flag[k]=cur_flag[k] | seq_ptr2[k];



		memcpy(temp_state,seq_ptr1,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr1,seq_ptr2,seq_len*sizeof(T_char_state));
		memcpy(seq_ptr2,temp_state,seq_len*sizeof(T_char_state));



		/*build the partial tree, insert the best_taxa at the branch
							(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
							tree[best_pos]*/
		memcpy(temp_tree,tree+best_pos,
		    (tree_size-best_pos));
		tree[best_pos]=internal_id;
		tree[best_pos+1]=i-1;
		memcpy(tree+best_pos+2,temp_tree,
		    (tree_size-best_pos));
		tree_size+=2;
		internal_id++;
	}


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(cost);
	free(temp_state);
	free(all_flag);
	free(cur_flag);
	free(temp_flag);
	free(rank);
	free(beat);
	free(min_pos);

	return(0);
}



/*Decide the order of sequence addition staticly by max-mini_lb rules*/
int static_reorder_max_mini_lb(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state mask;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*min_pos;  /*mininum cost and minimum position for each taxa*/
	int *ptr_cos,*ptr_pos,*med,*var;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	T_char_state *all_flag,*cur_flag,*temp_flag;

	all_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	cur_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	temp_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	min_pos=(int*)malloc(sizeof(int)*num_seq);
	med=(int*)malloc(sizeof(int)*num_seq);
	var=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

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

	for(i=0;i<seq_len;i++) cur_flag[i]=0;
	for(i=0;i<3;i++){
		seq_ptr1=seq_matrix+order[i]*seq_len;
		for(j=0;j<seq_len;j++,seq_ptr1++)
			cur_flag[j]=cur_flag[j] | (*seq_ptr1);
        }

	for(i=0;i<seq_len;i++) all_flag[i]=cur_flag[i];
	for(i=3;i<num_seq;i++){
		seq_ptr1=seq_matrix+order[i]*seq_len;
		for(j=0;j<seq_len;j++,seq_ptr1++)
			all_flag[j]=all_flag[j] | (*seq_ptr1);
	}

	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){

		/*for each taxon, compute the lb2*/
		for(j=i;j<=num_seq;j++){
			seq_ptr1=seq_matrix+order[j-1]*seq_len;
			lower_bound2[j]=0;
			for(k=0;k<seq_len;k++){
				temp_flag[k]=cur_flag[k] | seq_ptr1[k];
				temp_flag[k]=all_flag[k]-temp_flag[k];
				for(l=0,mask=1;l<8*sizeof(T_char_state);l++,mask=mask<<1){
					if(temp_flag[k] & mask) lower_bound2[j]++;
				}
			}
		}

		/*for each remaning sequences (s_i,....s_i), compute their best
							position (minimal cost) to insert into the current partial tree*/
		for(k=1;k<tree_size;k++){
			/*insert (internal_id,order[j]) right before 
tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));

			for(j=i;j<=num_seq;j++){
				temp_tree[k+1]=order[j-1];
				*(min_cost+(j-1)*(2*num_seq-2)+k)
				    =fitch_pass(temp_tree,2*i-2)+lower_bound2[j];
			}
		}

		/*for each taxa, compute the position of its minimum cost 
					and the median cost*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);

			/*compute median of the cost*/
			med[j-1]=0;
			for(k=1;k<tree_size;k++)
				med[j-1]+=ptr_cos[k];
			med[j-1]=med[j-1]/(tree_size-1);

			/*compute variance of the cost*/
			var[j-1]=0;
			for(k=1;k<tree_size;k++)
				var[j-1]+=(ptr_cos[k]-med[j-1])^2;

			ptr_pos=min_pos+j-1;
			*ptr_pos=1;  /*best position for taxa j*/
			for(k=2;k<tree_size;k++)
				if(ptr_cos[*ptr_pos]>ptr_cos[k])
					*ptr_pos=k;
		}

		/*find the maximum mini cost*/
		best_taxa=i;
		max_mini=*(min_cost+(i-1)*(2*num_seq-2)+min_pos[i-1]);
		for(j=i+1;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			if(max_mini < ptr_cos[min_pos[j-1]]) {
				best_taxa=j;
				max_mini=ptr_cos[min_pos[j-1]];
			}
		}

		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			if(max_mini == ptr_cos[min_pos[j-1]]) {
				if(var[j-1]>var[best_taxa-1])  best_taxa=j;
			}
		}

		best_pos=min_pos[best_taxa-1];
		lower_bound2[i]=lower_bound2[best_taxa];
		printf("i:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr2=seq_matrix+order[i-1]*seq_len;

		/*update cur_flag*/
		for(k=0;k<seq_len;k++)
			cur_flag[k]=cur_flag[k] | seq_ptr2[k];


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

	lower_bound2[num_seq]=0;

#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(med);
	free(var);
	free(min_cost);
	free(min_pos);
	free(temp_state);
	return(0);
}



/*Decide the order of sequence addition staticly by max-mini_lb rules*/
int old_static_reorder_max_mini_lb_cross(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state mask;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,min_pos;
	int cost,flag,cur_cost;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	T_char_state *all_flag,*cur_flag,*temp_flag;
	tNode_T  *best_tree; /* best partial tree */
	int *min_cost,*prune;

	all_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	cur_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	temp_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	best_tree=(tNode_T *)malloc(num_seq*sizeof(tNode_T));
	for(i=0;i<num_seq;i++)
	  best_tree[i]=(T_char_state 
*)malloc((2*num_seq-2)*sizeof(T_char_state));


	min_cost=(int *)malloc(sizeof(int)*num_seq);
	prune=(int *)malloc(sizeof(int)*num_seq);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

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

	for(i=0;i<seq_len;i++) cur_flag[i]=0;
	for(i=0;i<3;i++){
		seq_ptr1=seq_matrix+order[i]*seq_len;
		for(j=0;j<seq_len;j++,seq_ptr1++)
			cur_flag[j]=cur_flag[j] | (*seq_ptr1);
        }

	for(i=0;i<seq_len;i++) all_flag[i]=cur_flag[i];
	for(i=3;i<num_seq;i++){
		seq_ptr1=seq_matrix+order[i]*seq_len;
		for(j=0;j<seq_len;j++,seq_ptr1++)
			all_flag[j]=all_flag[j] | (*seq_ptr1);
	}

	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){

		/*for each taxon, compute the lb2*/
		for(j=i;j<=num_seq;j++){
			seq_ptr1=seq_matrix+order[j-1]*seq_len;
			lower_bound2[j]=0;
			for(k=0;k<seq_len;k++){
				temp_flag[k]=cur_flag[k] | seq_ptr1[k];
				temp_flag[k]=all_flag[k]-temp_flag[k];
				for(l=0,mask=1;l<8*sizeof(T_char_state);l++,mask=mask<<1){
					if(temp_flag[k] & mask) lower_bound2[j]++;
				}
			}
		}

		/*for each remaning sequences (s_i,....s_i), compute their best
							position (minimal cost) to insert into the current partial tree*/
		for(j=i;j<=num_seq;j++){
			/*find the best_tree for (0..i-1) and (j-1)*/

		   min_cost[j-1]=MAX_COST;
		   prune[j-1]=0;
		   /*insert (j-1) into tree:best_tree[order[i-2]]*/
		   for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before 
tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));
			temp_tree[k+1]=order[j-1];
			cost=fitch_pass(temp_tree,2*i-2)+lower_bound2[j];
			if(cost>global_ub) prune[j-1]++;
			if(cost<min_cost[j-1]){	
			   flag=1;
			   min_cost[j-1]=cost;
			   min_pos=k;
			}
		   }	
			/*insert (i-2) into best_tree[order[j-1]]*/
		   if(i>4) 
		     for(k=1;k<tree_size;k++){
                        /*insert (internal_id,i-2) right before
tree[k]*/
                        memcpy(temp_tree,best_tree[order[j-1]],k);
                        temp_tree[k]=internal_id;
                        memcpy(temp_tree+k+2,best_tree[order[j-1]]+k,
(tree_size-k));
                        temp_tree[k+1]=order[i-2];
                        cost=fitch_pass(temp_tree,2*i-2)+lower_bound2[j];
			if(cost>global_ub) prune[j-1]++;
                        if(cost<=min_cost[j-1]){
			   flag=2;
                           min_cost[j-1]=cost;
                           min_pos=k;
                        }
		     }
			

			/*save the best tree*/
			if (flag==1){
			  memcpy(temp_tree,tree,min_pos);
			 temp_tree[min_pos]=internal_id;
			 temp_tree[min_pos+1]=order[j-1];
			 memcpy(temp_tree+min_pos+2,tree+min_pos,
tree_size-min_pos);
			}
			else{
			 memcpy(temp_tree,best_tree[order[j-1]],min_pos);
                         temp_tree[min_pos]=internal_id;
                         temp_tree[min_pos+1]=order[i-2];
                         memcpy(temp_tree+min_pos+2,
best_tree[order[j-1]]+min_pos,tree_size-min_pos);

			}
		
			memcpy(best_tree[order[j-1]],temp_tree,
tree_size+2);

		}



		/*find the maximum mini cost*/
		best_taxa=i;
		max_mini=min_cost[i-1];
	        flag=prune[i-1];
		for(j=i+1;j<=num_seq;j++){
/*
			if((max_mini==min_cost[j-1]) && (flag<prune[j-1]))
			    { best_taxa=j; flag=prune[j-1];}
			if(max_mini < min_cost[j-1]) {
				best_taxa=j;
				max_mini=min_cost[j-1];
				flag=prune[j-1];
			}
*/
		if(max_mini<min_cost[j-1]){
		    max_mini=min_cost[j-1];
			best_taxa=j;
		}
/*		if (flag<prune[j-1]){
		     best_taxa=j;
                                max_mini=min_cost[j-1];
                                flag=prune[j-1];
		  }*/
		}


		lower_bound2[i]=lower_bound2[best_taxa];
		printf("i:%d,max_min1:%d,best_taxa:%d\n",i,max_mini,best_taxa);



		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr2=seq_matrix+order[i-1]*seq_len;

		/*update cur_flag*/
		for(k=0;k<seq_len;k++)
			cur_flag[k]=cur_flag[k] | seq_ptr2[k];


		/*build the partial tree, insert the best_taxa at the branch
							(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
							tree[best_pos]*/
		tree_size+=2;
		memcpy(tree,best_tree[order[i-1]],tree_size);
		internal_id++;
	}

	lower_bound2[num_seq]=0;

#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	for(i=0;i<num_seq;i++)
	   free(best_tree[i]);
	free(best_tree);

	free(min_cost);
	free(temp_state);
	free(all_flag);
	free(cur_flag);
	free(temp_flag);

	return(0);
}


/*Decide the order of sequence addition staticly by max-mini_lb rules*/
int static_reorder_max_mini_lb_cross(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state mask;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,min_pos;
	int cost,flag,cur_cost;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id,missed;
	int order0,order1,order2,temp_order;
	T_char_state *cur_flag,*temp_flag;
	tNode_T  *best_tree; /* best partial tree */
	int *min_cost,*prune;

	T_char_state *comb_seq;
        int row_size;
        int *include,*w,*score2;

       include=(int *)malloc(seq_len*sizeof(int));
       score2=(int *)malloc(num_seq*sizeof(int));
       row_size=seq_len*(seq_len+1)/2;
        w=(int *)malloc(seq_len*seq_len*sizeof(int));
          comb_seq=(T_char_state *)malloc(sizeof(T_char_state
)*num_seq*row_size);
	encode_2_sites(comb_seq);

	cur_flag=(T_char_state *)malloc(sizeof(T_char_state)*row_size);
	temp_flag=(T_char_state *)malloc(sizeof(T_char_state)*row_size);

	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	best_tree=(tNode_T *)malloc(num_seq*sizeof(tNode_T));
	for(i=0;i<num_seq;i++)
	  best_tree[i]=(T_char_state 
*)malloc((2*num_seq-2)*sizeof(T_char_state));


	min_cost=(int *)malloc(sizeof(int)*num_seq);
	prune=(int *)malloc(sizeof(int)*num_seq);

	internal_id=num_seq;
	tree=temp_memory.tree1;
	temp_tree=temp_memory.tree2;

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					order0=i;
					order1=j;
					order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

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

	for(i=0;i<row_size;i++) cur_flag[i]=0;
	for(i=0;i<3;i++){
		seq_ptr1=comb_seq+order[i]*row_size;
		for(j=0;j<row_size;j++,seq_ptr1++)
			cur_flag[j]=cur_flag[j] | (*seq_ptr1);
        }


	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){

		/*for each taxon, compute the lb2*/
		for(j=i;j<=num_seq;j++){
			seq_ptr1=comb_seq+order[j-1]*row_size;
			lower_bound2[j]=0;
			for(k=0;k<row_size;k++)
				temp_flag[k]=cur_flag[k] | seq_ptr1[k];
				lower_bound2[j]=compute_lb(temp_flag,
w,include,score2+j-1);
		}

		/*for each remaning sequences (s_i,....s_i), compute their best
							position (minimal cost) to insert into the current partial tree*/
		for(j=i;j<=num_seq;j++){
			/*find the best_tree for (0..i-1) and (j-1)*/

		   min_cost[j-1]=MAX_COST;
		   prune[j-1]=0;
		   /*insert (j-1) into tree:best_tree[order[i-2]]*/
		   for(k=1;k<tree_size;k++){
			/*insert (internal_id,j-1) right before 
tree[k]*/
			memcpy(temp_tree,tree,k);
			temp_tree[k]=internal_id;
			memcpy(temp_tree+k+2,tree+k,(tree_size-k));
			temp_tree[k+1]=order[j-1];

			cost=fitch_pass(temp_tree,2*i-2);
/*-lower_bound2[j];
*/			if(cost>global_ub) prune[j-1]++;
			if(cost<min_cost[j-1]){	
			   flag=1;
			   min_cost[j-1]=cost;
			   min_pos=k;
			}
		   }	
			/*insert missed into best_tree[order[j-1]]*/
		   if(i>4) 
		     for(k=1;k<tree_size;k++){
                        /*insert (internal_id,i-2) right before
tree[k]*/
                        memcpy(temp_tree,best_tree[order[j-1]],k);
                        temp_tree[k]=internal_id;
                        memcpy(temp_tree+k+2,best_tree[order[j-1]]+k,
(tree_size-k));
                        temp_tree[k+1]=order[missed];
                        cost=fitch_pass(temp_tree,2*i-2);
/*-lower_bound2[j];
*/			if(cost>global_ub) prune[j-1]++;
                        if(cost<=min_cost[j-1]){
			   flag=2;
                           min_cost[j-1]=cost;
                           min_pos=k;
                        }
		     }
			
			/*save the best tree*/
			if (flag==1){
			  memcpy(temp_tree,tree,min_pos);
			 temp_tree[min_pos]=internal_id;
			 temp_tree[min_pos+1]=order[j-1];
			 memcpy(temp_tree+min_pos+2,tree+min_pos,
tree_size-min_pos);
			}
			else{
			 memcpy(temp_tree,best_tree[order[j-1]],min_pos);
                         temp_tree[min_pos]=internal_id;
                         temp_tree[min_pos+1]=order[missed];
                         memcpy(temp_tree+min_pos+2,
best_tree[order[j-1]]+min_pos,tree_size-min_pos);

			}
		
			memcpy(best_tree[order[j-1]],temp_tree,
tree_size+2);

		}



		/*find the maximum mini cost*/
		best_taxa=i;
		max_mini=min_cost[i-1];
		for(j=i+1;j<=num_seq;j++){
		if(max_mini<min_cost[j-1]){
		    max_mini=min_cost[j-1];
			best_taxa=j;
		}
/*
		if (flag<prune[j-1]){
		     best_taxa=j;
                                max_mini=min_cost[j-1];
                                flag=prune[j-1];
		  }
*/
		}


		lower_bound2[i]=lower_bound2[best_taxa];
		printf("i:%d,max_min1:%d,best_taxa:%d\n",i,max_mini,best_taxa);



		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		seq_ptr2=seq_matrix+temp_order*seq_len;

                /*update cur_flag*/   
                for(k=0;k<seq_len;k++)
                        cur_flag[k]=cur_flag[k] | seq_ptr2[k];

		order[best_taxa-1]=order[i-1];
/*		if (prune[best_taxa-1]=1){*/
		  order[i-1]=temp_order;
		  missed=i-1;
/*		}else {
		   order[i-1]=order[i-2];
		   order[i-2]=temp_order;
		   missed=i-2;
		}
*/

		/*build the partial tree, insert the best_taxa at the branch
							(best_pos,best_pos->parent), i.e. insert(internal_id,i-1) right before
							tree[best_pos]*/
		tree_size+=2;
		memcpy(tree,best_tree[order[i-1]],tree_size);
		internal_id++;
	}

	lower_bound2[num_seq]=0;

#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	for(i=0;i<num_seq;i++)
	   free(best_tree[i]);
	free(best_tree);

	free(min_cost);
	free(temp_state);
	free(cur_flag);
	free(temp_flag);


	free(comb_seq);
	free(w);
	free(include);
	free(score2);

	return(0);
}


/*Decide the order of sequence addition staticly by max lowerbound of 
the partial tree*/
int static_reorder_max_lb(int *order)
{
	int i,j,k,l;
	T_char_state  *temp_state,*or_state;
	T_char_state *comb_seq;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int row_size,temp;
	int *cost,*include,*w;
	int best_taxa,max_lb,cur_cost;
	int internal_id;
	int order0,order1,order2,temp_order;
	int *score2;


	cost=(int *)malloc(num_seq*sizeof(int));
	include=(int *)malloc(seq_len*sizeof(int));
	score2=(int *)malloc(num_seq*sizeof(int));

	row_size=seq_len*(seq_len+1)/2;
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*row_size);
	or_state=(T_char_state *)malloc(sizeof(T_char_state)*row_size);
	w=(int *)malloc(seq_len*seq_len*sizeof(int));
	  comb_seq=(T_char_state *)malloc(sizeof(T_char_state 
)*num_seq*row_size);
	encode_2_sites(comb_seq);


	/*Decide the best core tree*/
	max_lb=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){

			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if(temp_state[l]==0){
				   temp++;
				   temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);		
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;
				for(l=0;l<seq_len;l++){
				     or_state[l]=temp_state[l] & seq_ptr3[l]; 
				    if(or_state[l]==0){
					cur_cost++;
					or_state[l]=temp_state[l] | seq_ptr3[l];
				     }
				}

				if(cur_cost>max_lb){
				  max_lb=cur_cost;
				  order0=i;
				  order1=j;
				  order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}



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

	for(i=0;i<row_size;i++) temp_state[i]=0;
	for(i=0;i<3;i++){
		seq_ptr1=comb_seq+order[i]*row_size;
		for(j=0;j<row_size;j++,seq_ptr1++)
			temp_state[j]=temp_state[j] | (*seq_ptr1);
        }


	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<num_seq;i++){

		/*for each taxon, compute the lowerbound*/
		for(j=i;j<=num_seq;j++){
			seq_ptr1=comb_seq+order[j-1]*row_size;
			for(k=0;k<row_size;k++)
			  or_state[k]=temp_state[k] | seq_ptr1[k];
			cost[j-1]=compute_lb(or_state,w,include,score2+j-1);	
		}

		/*find the best taxa which has maximum cost
		best_taxa=i;
		max_lb=score2[i-1];;
		for(j=i+1;j<=num_seq;j++){
			if(max_lb < score2[j-1]) {
				best_taxa=j;
				max_lb=score2[j-1];
			}
		}

		for(j=i;j<=num_seq;j++)
		  if(score2[j-1]==max_lb)
		{  if(cost[best_taxa-1]<cost[j-1])
			best_taxa=j;
		}
	*/

		/*find the taxa with maximum cost*/
               best_taxa=i;
                max_lb=cost[i-1];;
                for(j=i+1;j<=num_seq;j++){
                        if(max_lb < cost[j-1]) {
                                best_taxa=j;
                                max_lb=cost[j-1];
                        }
                }
		
		
printf("i:%d,max_lb:%d,best_taxa:%d\n",i,max_lb,best_taxa);

		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;

		/*update temp_state*/
		seq_ptr2=comb_seq+order[i-1]*row_size;
		for(k=0;k<row_size;k++)
			temp_state[k]=temp_state[k] | seq_ptr2[k];


	}


#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

#endif

	free(cost);
	free(or_state);
	free(temp_state);
	free(comb_seq);
	free(w);
	free(include);
	free(score2);

	return(0);
}




/*compute the lowerbond by maximum matching*/
int compute_lb(T_char_state *or_state,int *w,int *include,int *score2)
{ 	int row_size,score,weight;
	T_char_state mask,*ptr_or;
	int i,j,k,diff;
	int w_ii,w_ij,w_jj;

	row_size=seq_len*(seq_len+1)/2;

	/*compute the weight*/
	ptr_or=or_state;
	for(i=0;i<seq_len;i++)
    	    for(j=i;j<seq_len;j++,ptr_or++){
        	/*compute #of differtn states when combining (i,j)*/
        	weight=0;
        	mask=1;
        	for(k=0;k<8*sizeof(T_char_state);k++,mask=mask<<1)
          	    if(((*ptr_or) & mask) !=0) weight++;

        	/*w[i][j]=w[j][i]=weight*/
       		*(w+i*seq_len+j)=weight;
        	*(w+j*seq_len+i)=weight;
	}

	/*find the nearly maximum matching*/
	for(i=0;i<seq_len;i++) include[i]=-1;
	score=0;
	*score2=0;
	for(i=0;i<seq_len;i++)
	  if (include[i]<0){
	    w_ii=*(w+i*seq_len+i);
	    for(j=i+1;j<seq_len;j++)
	        if (include[j]<0){
          		w_ij=*(w+i*seq_len+j);
          		w_jj=*(w+j*seq_len+j);
			diff=w_ij-w_ii-w_jj+1;
			if(diff>0){
			               include[i]=j;   
     		        include[j]=i;
            		 score=score+w_ij-1;
			 *score2=(*score2)+diff;
            		 break;
			}
		}
	  }



	for(i=0;i<seq_len;i++)
    	   if(include[i]<0){
        	w_ii=*(w+i*seq_len+i);
        	score+=w_ii-1;
	   }



	return(score);	
}


/*Decide the order of sequence addition staticly by max-mini_lb rules
use constant Fitch algorithm*/
int static_reorder_max_mini_lb_cons(int *order)
{
	int temp;
	int i,j,k,l;
	T_char_state mask;
	T_char_state  *temp_state;
	T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;
	int best_taxa,best_pos,max_mini,cur_cost;
	int *min_cost,*min_pos;  /*mininum cost and minimum position for each taxa*/
	int *ptr_cos,*ptr_pos,*med,*var;
	tNode_T  tree,temp_tree;
	int tree_size;
	int internal_id;
	int order0,order1,order2,temp_order;
	T_char_state *all_flag,*cur_flag,*temp_flag;
	T_char_state temp_state1;
	int tree_len,delta;
	T_char_state *ptr_seq1,*ptr_seq2,*ptr_seq3;

	all_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	cur_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);
	temp_flag=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

	min_cost=(int*)malloc(sizeof(int)*num_seq*(2*num_seq-2));
	min_pos=(int*)malloc(sizeof(int)*num_seq);
	med=(int*)malloc(sizeof(int)*num_seq);
	var=(int*)malloc(sizeof(int)*num_seq);
	temp_state=(T_char_state *)malloc(sizeof(T_char_state)*seq_len);

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

	for(i=0;i<seq_len;i++) cur_flag[i]=0;
	for(i=0;i<3;i++){
		seq_ptr1=seq_matrix+order[i]*seq_len;
		for(j=0;j<seq_len;j++,seq_ptr1++)
			cur_flag[j]=cur_flag[j] | (*seq_ptr1);
        }

	for(i=0;i<seq_len;i++) all_flag[i]=cur_flag[i];
	for(i=3;i<num_seq;i++){
		seq_ptr1=seq_matrix+order[i]*seq_len;
		for(j=0;j<seq_len;j++,seq_ptr1++)
			all_flag[j]=all_flag[j] | (*seq_ptr1);
	}


	/*Decide which sequence to be added at the i_th step*/
	for(i=4;i<=num_seq;i++){

		/*for each taxon, compute the lb2*/
		for(j=i;j<=num_seq;j++){
			seq_ptr1=seq_matrix+order[j-1]*seq_len;
			lower_bound2[j]=0;
			for(k=0;k<seq_len;k++){
				temp_flag[k]=cur_flag[k] | seq_ptr1[k];
				temp_flag[k]=all_flag[k]-temp_flag[k];
				for(l=0,mask=1;l<8*sizeof(T_char_state);l++,mask=mask<<1){
					if(temp_flag[k] & mask) lower_bound2[j]++;
				}
			}
		}

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
		    *(min_cost+(j-1)*(2*num_seq-2)+k)=tree_len+delta+
lower_bound2[j];
		  }
		}

		/*for each taxa, compute the position of its minimum cost 
					and the median cost*/
		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);

			/*compute median of the cost*/
			med[j-1]=0;
			for(k=1;k<tree_size;k++)
				med[j-1]+=ptr_cos[k];
			med[j-1]=med[j-1]/(tree_size-1);

			/*compute variance of the cost*/
			var[j-1]=0;
			for(k=1;k<tree_size;k++)
				var[j-1]+=(ptr_cos[k]-med[j-1])^2;

			ptr_pos=min_pos+j-1;
			*ptr_pos=1;  /*best position for taxa j*/
			for(k=2;k<tree_size;k++)
				if(ptr_cos[*ptr_pos]>ptr_cos[k])
					*ptr_pos=k;
		}

		/*find the maximum mini cost*/
		best_taxa=i;
		max_mini=*(min_cost+(i-1)*(2*num_seq-2)+min_pos[i-1]);
		for(j=i+1;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			if(max_mini < ptr_cos[min_pos[j-1]]) {
				best_taxa=j;
				max_mini=ptr_cos[min_pos[j-1]];
			}
		}

		for(j=i;j<=num_seq;j++){
			ptr_cos=min_cost+(j-1)*(2*num_seq-2);
			if(max_mini == ptr_cos[min_pos[j-1]]) {
				if(var[j-1]>var[best_taxa-1])  best_taxa=j;
			}
		}

		best_pos=min_pos[best_taxa-1];
		lower_bound2[i]=lower_bound2[best_taxa];

#ifdef DEBUG1
		printf("i:%d,max_min:%d,best_taxa:%d\n",i,max_mini,best_taxa);
#endif
		/*Reorder the sequence, let the i-th taxa to be best_taxa*/
		temp_order=order[best_taxa-1];
		order[best_taxa-1]=order[i-1];
		order[i-1]=temp_order;
		seq_ptr2=seq_matrix+order[i-1]*seq_len;

		/*update cur_flag*/
		for(k=0;k<seq_len;k++)
			cur_flag[k]=cur_flag[k] | seq_ptr2[k];


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

	lower_bound2[num_seq]=0;

#ifdef DEBUG1
	printf("\n The order of addition:\n");
	for(j=0;j<num_seq;j++)
		printf("%d ",order[j]);

	seq_ptr1=seq_matrix;
	printf("\n\nThe sequence matrix after reordering:\n");
	for(i=0;i<seq_len;i++){
		for(j=0;j<num_seq;j++,seq_ptr1++)
			printf("%d,",*seq_ptr1);
		printf("\n");
	}
#endif

	free(med);
	free(var);
	free(min_cost);
	free(min_pos);
	free(temp_state);
	return(max_mini);
}


/*called by static_reorder,process tree use full seq_len,return 
tree_len*/
int get_branch_interval(tNode_T tree,int tree_size)
{       int i,j;
	 int *p,*b;      
        int tree_len,delta;
        T_char_state *ptr_seq1,*ptr_seq2,*ptr_seq3;
	T_char_state temp_state;	
        
        p=temp_memory.p;
        b=temp_memory.b; 

        /*Pass par_sol in post-order*/
        tree_len=fitch_pass_2(tree,tree_size,p,b);
                                        
        /*Pass par_sol in pre-order*/
        ptr_seq1=seq_matrix+tree[0]*seq_len;
        ptr_seq2=seq_matrix+(tree[1]+2*num_seq)*seq_len;
        for(j=0;j<seq_len;j++){ ptr_seq2[j]=ptr_seq1[j];
        }       
	
	for(i=2;i<tree_size;i++){
          ptr_seq1=seq_matrix+b[tree[i]]*seq_len;
          ptr_seq2=seq_matrix+(p[tree[i]]+2*num_seq)*seq_len;
          ptr_seq3=seq_matrix+(tree[i]+2*num_seq)*seq_len;

          for(j=0;j<seq_len;j++)
          {
                ptr_seq3[j]=ptr_seq1[j] & ptr_seq2[j];
                if (ptr_seq3[j]==0) ptr_seq3[j]=ptr_seq1[j] | 
ptr_seq2[j];
          }
        }

	/*compute branch interval*/
        for(i=1;i<tree_size;i++){
          ptr_seq1=seq_matrix+tree[i]*seq_len;
          ptr_seq2=seq_matrix+(tree[i]+2*num_seq)*seq_len;
          for(j=0;j<seq_len;j++){
                temp_state=ptr_seq1[j] & ptr_seq2[j];
                if (temp_state==0) temp_state=ptr_seq1[j] | ptr_seq2[j];
                ptr_seq2[j]=temp_state;
          }
        }

	return(tree_len);
}




/*Decide the best core tree that has the maximum tree cost, 
order0 ,order1,order2*/
int decide_core_tree(T_char_state *temp_state,int *order0,int *order1,int *order2)
{
	int temp,max_mini,cur_cost;
	int i,j,k,l;
        T_char_state *seq_ptr1,*seq_ptr2,*seq_ptr3;	

	/*Decide the best core tree*/
	max_mini=0;
	seq_ptr1=seq_matrix;
	for(i=0;i<num_seq;i++){
		seq_ptr2=seq_ptr1+seq_len;
		for(j=i+1;j<num_seq;j++){
			/*computer Farris interval between i and j*/
			temp=0;
			for(l=0;l<seq_len;l++){
				temp_state[l]=(seq_ptr1[l] & seq_ptr2[l]);
				if (temp_state[l] == 0){
					temp++;
					temp_state[l]=(seq_ptr1[l] | seq_ptr2[l]);
				}
			}

			/*add the third taxon*/
			seq_ptr3=seq_ptr2+seq_len;
			for(k=j+1;k<num_seq;k++){
				cur_cost=temp;

				for(l=0;l<seq_len;l++)
					if ((temp_state[l] & seq_ptr3[l]) == 0)
						cur_cost++;

				if(cur_cost>max_mini){
					max_mini=cur_cost;
					*order0=i;
					*order1=j;
					*order2=k;
				}
				seq_ptr3+=seq_len;
			}
			seq_ptr2+=seq_len;
		}
		seq_ptr1+=seq_len;
	}

  return(0);
}
