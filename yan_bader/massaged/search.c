/*This file contains functions for  different search strategies*/
#include <string.h>
#include <stdio.h>
#include "search.h"

/* The first choice is to combine depth-first and best-first;
choose the node with maximum depth,when there is tie, choos the node with
minimum lowerbound. i.e. greedy + back track*/

/*Intialize the frointer*/
int initalize_search(int int_gb)
{ 
	tNode_T *tree; /*Used as a stack to save partial tree at each level*/
	int *pos; /*the insert postion in the tree, organized as a stack*/
	int lowerbound,temp_gub;
	int i,j,ival,tree_size,internal_id;
	tNode_T temp_tree;
	struct hNode *hnode;

	opt_list=NULL;
	temp_memory.tree_pool_ptr=0;

	/*Intialize frontier:Each lever has a pointer to a priority queue with
	the key is lowerbound*/
	for(i=thresh_level;i<=num_seq;i++)
		frontier[i-thresh_level].heap_size=0;

	/*Generate all the nodes from 3 to thresh_level*/
	/*malloc memory space for tree:the stack*/
	pos=(int *)malloc(sizeof(int)*(thresh_level-2));
	tree=(T_char_state 
**)malloc(sizeof(T_char_state *)*(thresh_level-2));
	for(i=3;i<=thresh_level;i++)
		tree[i-3]=(T_char_state 
*)malloc((2*i-2)*sizeof(T_char_state));

	tree_size=2*thresh_level-2;
	for(i=3;i<=thresh_level;i++) pos[i-3]=1;
	/*generate the first tree*/
	internal_id=num_seq;
	tree[0][0]=order[0];
	tree[0][1]=internal_id;
	tree[0][2]=order[1];
	tree[0][3]=order[2];
	internal_id++;

	for(i=3;i<thresh_level;i++){
		memcpy(&tree[i-2][0],&tree[i-3][0],1);
		tree[i-2][1]=internal_id++;
		tree[i-2][2]=order[i];
		memcpy(&tree[i-2][3],&tree[i-3][1],2*i-3);
	}

	if(thresh_level==3){
		temp_tree=frontier[0].heap[frontier[0].heap_size].par_sol;
		memcpy(temp_tree,tree[0],tree_size);
/*		lowerbound=
		    fitch_pass(temp_tree,tree_size)+lower_bound2[3];*/
		
lowerbound=new_fitch_pass_3(temp_tree,thresh_level)+infor_len[thresh_level]+
lower_bound2[thresh_level];
		hnode=enqueue(0,temp_tree,lowerbound);
		memcpy(hnode->internal_states,seq_matrix+num_seq*seq_len,
(thresh_level-2)*seq_len);
	}
	else {
		ival=0;
		for(i=0;i<frontier[0].array_size;i++){
			if(ival<0) break;

			temp_tree=frontier[0].heap[frontier[0].heap_size].par_sol;
			memcpy(temp_tree,tree[thresh_level-3],tree_size);
/*
			lowerbound=
			    fitch_pass(temp_tree,tree_size)+lower_bound2[thresh_level];
*/
lowerbound=new_fitch_pass_3(temp_tree,thresh_level)+infor_len[thresh_level]+
lower_bound2[thresh_level];

			if(lowerbound <= global_ub){
			   hnode=enqueue(0,temp_tree,lowerbound);
			   memcpy(hnode->internal_states,
			seq_matrix+num_seq*seq_len,(thresh_level-2)*seq_len);
			}
			ival=gen_tree(thresh_level,tree,pos);
		}
	}




	return(0);
}

/*Remove the minimal item heap[0] from the heap*/
int dequeue(int i)
{ 
	int smallest,s_id,p_id,l_id,r_id;
	struct hNode *heap,*temp,*min;
	int heap_size;
	tNode_T temp_par_sol;

	heap_size=frontier[i].heap_size;
	frontier[i].heap_size--;
	heap=frontier[i].heap;

	if(frontier[i].heap_size==0)  return(0);

	temp_par_sol=heap[0].par_sol; /* don't let the memory space pointer by heap[0].par_sol
	leak*/

	/*maintain the property of heap*/
	temp=&heap[heap_size-1];
	p_id=1;
	l_id=2*p_id;
	r_id=2*p_id+1;

	s_id=l_id;
	if(r_id <heap_size)
		if( heap[s_id-1].lowerbound > heap[r_id-1].lowerbound)
			s_id=r_id;

	while (heap[s_id-1].lowerbound < temp->lowerbound){
		/*copy smallest to heap[p_id-1]*/
		heap[p_id-1].lowerbound=heap[s_id-1].lowerbound;
		heap[p_id-1].upperbound=heap[s_id-1].upperbound;
		heap[p_id-1].par_sol=heap[s_id-1].par_sol;

		p_id=s_id;

		/*Get the smaller from the its children*/
		l_id=2*p_id;
		if(l_id >=heap_size) break;
		r_id=2*p_id+1;
		s_id=l_id;
		if (r_id<heap_size)
			if( heap[s_id-1].lowerbound > heap[r_id-1].lowerbound)
				s_id=r_id;


	}

	/*copy temp to heap[p_id-1]*/
	heap[p_id-1].lowerbound=temp->lowerbound;
	heap[p_id-1].upperbound=temp->upperbound;
	heap[p_id-1].par_sol=temp->par_sol;

	temp->par_sol=temp_par_sol;

	return(0);
}


/*Expand par_sol: inserting all of its children to the frontier[id] if
they pass the lowerbound test*/
int expand(struct hNode *sel_hnode, int id)
{ 
	tNode_T par_sol,temp_tree;
	int pos,tree_size;
	int tree_len,lowerbound,upperbound,taxa;
	struct hNode *heap;
	int h_id,p_id,par_score,internal_id;
	SOL_LIST_T *temp_list;
	T_tree_node *temp_true_tree;
	int gap,ival;
	int temp_tree_len,old_tree_len;

/*	if(sel_hnode->lowerbound > global_ub) return (0);  not apply 
for loose lowerbound
*/
	par_sol=sel_hnode->par_sol;
	taxa=id+thresh_level;

	tree_size=frontier[id-1].tree_size;
	heap=frontier[id].heap;
	internal_id=num_seq+taxa-3;
/*	tree_len=fitch_pass(par_sol,tree_size);*/
	old_tree_len=new_fitch_pass(par_sol,taxa-1);

	gap=global_ub-lower_bound2[taxa];

#if 0
   printf("At level :%d, gap1:%d,gap2:%d\n",taxa,global_ub-tree_len,
global_ub-sel_hnode->lowerbound);
#endif

	/*change the per-order encoding to true tree*/
	temp_true_tree=temp_memory.temp_true_tree;
	encode_to_true_tree(par_sol,tree_size,temp_true_tree);

	for(pos=1;pos<tree_size;pos++){
		/*insert (internal_id,taxa-1) before pos*/

/*		
ival=execede_gap(par_sol[pos],temp_true_tree[par_sol[pos]].parent,taxa-1);

		tree_len=ival+old_tree_len+infor_len[taxa-1];
		if (tree_len>=gap)
*/


tree_len=new_fitch_incremental_pass(temp_true_tree,par_sol[pos],
taxa-1);
lowerbound=tree_len+lower_bound2[taxa];


/*
                /*insert (internal_id,taxa-1) before pos
                temp_tree=heap[frontier[id].heap_size].par_sol;
                memcpy(temp_tree,par_sol,pos);
                temp_tree[pos]=internal_id;
                temp_tree[pos+1]=order[taxa-1];
                
memcpy(temp_tree+pos+2,par_sol+pos,tree_size-pos);

		tree_len=fitch_pass(temp_tree,tree_size+2);
		lowerbound=tree_len+lower_bound2[taxa];

/*
if(( lowerbound<global_ub-2) && (1)&& (taxa<=num_seq-3))
{ /*compute the maximum unresolved discrepandcy
 tree_len=compue_max_undis(temp_tree,tree_size+2,taxa,global_ub-lower_bound2[taxa]);
lowerbound=tree_len+lower_bound2[taxa];
}
*/

		if(lowerbound <= global_ub) {
                /*insert (internal_id,taxa-1) before pos*/
                temp_tree=heap[frontier[id].heap_size].par_sol;   
                memcpy(temp_tree,par_sol,pos);
                temp_tree[pos]=internal_id;
                temp_tree[pos+1]=order[taxa-1];
                memcpy(temp_tree+pos+2,par_sol+pos,tree_size-pos);

			if(taxa==num_seq){
				if(lowerbound ==global_ub){
					/*insert par_sol into the list of optimal solution:opt_list*/
					temp_list=(SOL_LIST_T *)malloc(sizeof(SOL_LIST_T));
					assert(temp_memory.tree_pool_ptr < POOL_SIZE);
					/*		if(temp_memory.tree_pool_ptr >= POOL_SIZE){
							   p_error(__FILE__,__LINE__,"Tree Pool Overflow");
							}
					*/
					temp_list->tree=temp_memory.tree_pool[temp_memory.tree_pool_ptr];
					temp_memory.tree_pool_ptr++;
					memcpy(temp_list->tree,temp_tree,tree_size+2);

					temp_list->next=opt_list;
					opt_list=temp_list;
				}
				else {
					global_ub=lowerbound;
					/*Delete the old  opt_list*/
					while(opt_list){
						temp_list=opt_list->next;
						free(opt_list);
						opt_list=temp_list;
					}
					/*insert temp_tree into opt_list*/
					opt_list=(SOL_LIST_T *)malloc(sizeof(SOL_LIST_T));
					opt_list->tree=temp_memory.tree_pool[0];
					memcpy(opt_list->tree,temp_tree,tree_size+2);
					opt_list->next=NULL;
					temp_memory.tree_pool_ptr=1;
				}
			} else{
				/*insert par_sol into the priority queue at frontier[id]*/
				enqueue(id,temp_tree,lowerbound);
			}
		}
	}


	return(0);
}


/*Compute upper_bound2: at level i, it is (the cost of tree on the remaining taxa
generated by Neighbor_Joining) + (the minimum distance between one taxa from the
partial tree and one taxa from the remaininh taxa*/
int ub_part2()
{ 
	int i,tree_len,min_dist;
	tNode_T tree;
	int *dist_matrix;

	dist_matrix=temp_memory.dist_matrix;
	tree=temp_memory.tree1;

	/*compute the distance matrix*/
	compute_dist_matrix(dist_matrix);

	for(i=3;i<num_seq-1;i++){
		/*Get the Neighbor-Joining tree onthe taxa (i)...n*/
		Adjusted_SK_NJ(dist_matrix,tree,i);
		tree_len=fitch_pass(tree,2*(num_seq-i));

		/*Get the minimum distance between  one taxa from 1...(i-1) and
		one from i..n*/
		min_dist=compute_min_dist(dist_matrix,i);
		upper_bound2[i-1]=tree_len+min_dist; /*pper_bond2[0] is the NJ
		tree on all taxa*/

	}

	/*computet upper_bound2[num_seq-3]
	  tree_len=*(dist_matrix+(num_seq-2)*num_seq+num_seq-1);
	  min_dist=compute_min_dist(dist_matrix,num_seq-2);
	  upper_bound2[num_seq-3]=tree_len+min_dist;*/

	/*computet upper_bound2[num_seq-2]*/
	min_dist=compute_min_dist(dist_matrix,num_seq-1);
	tree_len=*(dist_matrix+(num_seq-2)*num_seq+(num_seq-1));
	upper_bound2[num_seq-2]=min_dist+tree_len;

	min_dist=compute_min_dist(dist_matrix,num_seq);
	upper_bound2[num_seq-1]=min_dist;


#ifdef DEBUG1
	printf("\nThe upper_bound2 from taxa 3:\n");
	for(i=2;i<num_seq;i++)
		printf("%d  ",upper_bound2[i]);
	printf("\n");
	fflush(stdout);
#endif

	return(0);
}


/*Generate tree*/
int gen_tree(int level, tNode_T *tree, int *pos)
{ 
	int i,j;


	/*find the next insert position*/
	i=level-1;
	pos[i-3]++;
	while(pos[i-3] > (2*i-3)){
		pos[i-3]=1;
		i--;
		if(i<3) return (-1);
		else pos[i-3]++;
	}


	for(j=i;j<=level-1;j++){
		/*insert (num_seq+j-2,j) in tree[j-3] right before
		pos[j-3] and get tree[j+1-3]*/
		memcpy(&tree[j-2][0],&tree[j-3][0],pos[j-3]);
		tree[j-2][pos[j-3]]=num_seq+j-2,
		    tree[j-2][pos[j-3]+1]=order[j];
		memcpy(tree[j-2]+pos[j-3]+2,tree[j-3]+pos[j-3],
		    2*j-2-pos[j-3]);

	}

	return(0);
}


/*Get the Neighbor-Joining tree onthe taxa start_taxa...n*/
int Adjusted_SK_NJ(int *dist_matrix,tNode_T tree,int start_taxa)
{ 
	int *temp_matrix; /*distance matrix on taxa i..n*/
	int *index; /*point to the index in tree_pool,size (n-start_taxa+1),
	when a taxa is discard, it is -1*/
	struct tree *new_node,*tree_pool,*true_tree,*root;
	int taxa,i,j,k,tree_ptr,num_taxa;
	int *R;
	int min_i,min_j,min_s,s;
	int d_iu,d_ju,d_ij,d_ik,d_jk,d_uk;

	num_taxa=num_seq-start_taxa+1;
	temp_matrix=(int *)malloc(sizeof(int)*num_taxa*num_taxa);
	/*get the submatrix from dist_matrix*/
	for(i=0;i<num_taxa;i++)
		for(j=0;j<num_taxa;j++)
			*(temp_matrix+i*num_taxa+j)=
			    *(dist_matrix+(i+start_taxa-1)*num_seq+(j+start_taxa-1));

#ifdef DEBUG1
	printf("\n temp_matrix:\n");
	R=temp_matrix;
	for(i=0;i<num_taxa;i++){
		for(j=0;j<num_taxa;j++,R++)
			printf("%d ",*R);
		printf("\n");
	}
#endif

	index=(int *)malloc(sizeof(int)*num_taxa);
	R=(int *)malloc(sizeof(int)*num_taxa);
	tree_pool=temp_memory.true_tree_pool;

	for(i=0;i<num_taxa;i++){
		index[i]=i;
		tree_pool[i].parent=NULL;
		tree_pool[i].lChild=NULL;
		tree_pool[i].rChild=NULL;
		tree_pool[i].tag=i+start_taxa;
	}
	tree_ptr=num_taxa;

	/*Generated a rooted binary tree by Neighbor-Joining*/
	for(taxa=0;taxa<num_taxa-2;taxa++){
		/*compute R_i=sum of D_ik*/
		for(i=0;i<num_taxa;i++)
			if(index[i]>=0)
			{
				R[i]=0;
				for(k=0;k<num_taxa;k++)
					if(index[k]>=0) R[i]+=*(temp_matrix+i*num_taxa+k);

			}

		/*initialize min_i.min_j,min_s*/
		for(i=0;i<num_taxa;i++)
			if(index[i]>=0) {
				min_i=i;
				break;
			}
		for(j=i+1;j<num_taxa;j++)
			if(index[j]>=0){
				min_j=j;
				break;
			}
		min_s=(num_taxa-2-taxa)*(*(temp_matrix+min_i*num_taxa+min_j)-R[min_i]
		    -R[min_j]);

		/*choose the smallest S_ij*/
		for(i=min_i;i<num_taxa;i++)
			if(index[i]>=0)
				for(j=i+1;j<num_taxa;j++)
					if(index[j]>=0){
						s=(num_taxa-2-taxa)*(*(temp_matrix+i*num_taxa+j)-R[i]
						    -R[j]);
						if(s<min_s){
							min_s=s;
							min_i=i;
							min_j=j;
						}
					}

		/*Cluster min_i and min_j  as U*/
		new_node=&tree_pool[tree_ptr];
		new_node->parent=NULL;
		new_node->lChild=&tree_pool[index[min_i]];
		new_node->rChild=&tree_pool[index[min_j]];
		new_node->tag=0;
		tree_pool[index[min_i]].parent=new_node;
		tree_pool[index[min_j]].parent=new_node;
		index[min_i]=tree_ptr;
		index[min_j]=-1;
		tree_ptr++;


		/*compute the new branch (min_i,u), (min_j,u)*/
		d_ij=*(temp_matrix+min_i*num_taxa+min_j);
		d_iu=((num_taxa-taxa-2)*d_ij+R[min_i]-R[min_j])/(2*(num_taxa-taxa-2));
		d_ju=((num_taxa-taxa-2)*d_ij-R[min_i]+R[min_j])/(2*(num_taxa-taxa-2));

		/*update temp_matrix*/
		for(k=0;k<num_taxa;k++)
			if (index[k]>=0){
				d_ik=*(temp_matrix+min_i*num_taxa+k);
				d_jk=*(temp_matrix+min_j*num_taxa+k);
				d_uk=d_ik-d_iu;
				if(d_uk < d_jk-d_ju)
					d_uk=d_jk-d_ju; /*SK_NJ IS ADJUSTIFIED HERE, in SK_NJ,
						d_uk=(d_ik-d_iu +d_jk-d_ju)/2*/
				R[k]=d_uk;  /*Let R save D_uk temporarily*/
			}

		/*Copy R to the min_i row and min_column*/
		for(k=0;k<num_taxa;k++)
			if (index[k]>=0){
				*(temp_matrix+min_i*num_taxa+k)=R[k];
				*(temp_matrix+k*num_taxa+min_i)=R[k];
			}

	}

	/*Join the last 2 nodes */
	min_i=-1;
	for(i=0;i<num_taxa;i++)
		if(index[i]>=0)
			if(min_i<0) min_i=i;
			else {
				min_j=i;
				break;
			}

	/*Cluster min_i and min_j  as U*/
	new_node=&tree_pool[tree_ptr];
	new_node->parent=NULL;
	new_node->lChild=&tree_pool[index[min_i]];
	new_node->rChild=&tree_pool[index[min_j]];
	new_node->tag=0;
	tree_pool[index[min_i]].parent=new_node;
	tree_pool[index[min_j]].parent=new_node;
	index[min_i]=tree_ptr;

	/*Get the pre-order traver of the true_tree and save it into tree*/
	true_tree=&tree_pool[tree_ptr];

	/*change the rooted binary tree to unrooted tree, suppose a terminal as
	root*/
	root=rooted_to_unrooted(true_tree);

	/*traverse true_tree->lChild in in-order 
	  i=in_order_walk(true_tree->lChild,tree,0); */
	/*traverse true_tree in pre-order*/
	tree[0]=root->tag;
	pre_order_walk(root->lChild,tree,1);

	free(temp_matrix);
	free(index);
	free(R);

	return(0);
}


/*Get the minimum distance between  one taxa from 1...(i-1) and
one from i..n*/
int compute_min_dist(int *dist_matrix,int i)
{ 
	int min_dist,d;
	int j,k;

	if(i<=1) return(0);

	/*Initialize min_dist=dist_matrix[0,i-1]*/
	min_dist=*(dist_matrix+i-1);

	for(j=0;j<i-1;j++)
		for(k=i-1;k<num_seq;k++){
			d=*(dist_matrix+j*num_seq+k);
			if(min_dist>d) min_dist=d;
		}

	return(min_dist);
}


/*traverse true_tree in in-order and save the tags of visited nodes in
tree started from i, return next avaliable positionin tree*/
int in_order_walk(struct tree *true_tree, tNode_T tree, int i)
{ 
	int pos;

	if(true_tree->tag>0) {
		tree[i]=true_tree->tag;
		return(i+1);
	}

	pos=in_order_walk(true_tree->lChild,tree,i);
	tree[pos]=0;
	pos=in_order_walk(true_tree->rChild,tree,pos+1);

	return(pos);
}


/*traverse true_tree in pre-order and save the tags of visited nodes in
tree started from i, return next avaliable position in tree*/
int pre_order_walk(struct tree *true_tree, tNode_T tree, int i)
{ 
	int pos;

	if(true_tree->tag>0) {
		tree[i]=true_tree->tag;
		return(i+1);
	}

	tree[i]=0;
	pos=pre_order_walk(true_tree->lChild,tree,i+1);
	pos=pre_order_walk(true_tree->rChild,tree,pos);

	return(pos);

}


/*change the rooted binary tree to unrooted tree, suppose a terminal as
root*/
struct tree   *rooted_to_unrooted(struct tree *old_root)
{
	struct tree *old_parent,*old_lChild;
	struct tree *new_root;
	struct tree *temp;

	/* step 1:find the leftmore leaf and let it be the new root*/
	temp=old_root;
	while (temp->lChild) temp=temp->lChild;
	new_root=temp;

	/*step 2: for  the nodes on the path from new_root to old_root, exchange
	its parent and leftchild*/
	temp=new_root->parent;
	new_root->parent=NULL;
	while (temp!=old_root){
		old_parent=temp->parent;

		old_lChild=temp->lChild;
		old_lChild->lChild=temp;
		temp->parent=old_lChild;

		temp=old_parent;
	}

	/*step 3: remove old_root from the tree, let its lChild be the parent
	of its rChild*/
	old_lChild=old_root->lChild;
	temp=old_root->rChild;
	old_lChild->lChild=temp;
	temp->parent=old_lChild;

	return (new_root);
}


/*insert the partial tree with lowerbound into the heap of 
frontier[frontier_id]*/
struct hNode * enqueue(int id,tNode_T temp_tree,int lowerbound)
{ 
	int h_id,p_id;
	struct hNode *heap;

	frontier[id].heap_size++;  /*Assume it never overflow*/
	heap=frontier[id].heap;
	h_id=frontier[id].heap_size;
	p_id=h_id/2;
	/*Sift up*/
	if(p_id>0)
		while(heap[p_id-1].lowerbound> lowerbound){
			heap[h_id-1].lowerbound=heap[p_id-1].lowerbound;

			heap[h_id-1].par_sol=heap[p_id-1].par_sol;
			/*memcpy(heap[h_id-1].par_sol,heap[p_id-1].par_sol,tree_size+2);*/

			h_id=p_id;
			p_id=h_id/2;
			if(p_id<1)break;
		}
	/*Copy the new node at heap[h_id-1]*/
	heap[h_id-1].lowerbound=lowerbound;
	/*memcpy(heap[h_id-1].par_sol,temp_tree,tree_size+2);*/
	heap[h_id-1].par_sol=temp_tree;



	return(&heap[h_id-1]);
}


/*From the per-order encodeing, generate the true tree, for each node, 
find its
parent and brother*/
int     encode_to_true_tree(tNode_T tree,int tree_size,T_tree_node 
*true_tree)
{ int *stack;
  int i,stack_ptr,child1,child2;

  stack=temp_memory.stack;
  stack_ptr=0;

  for(i=tree_size-1;i>0;i--){
    if(tree[i]<num_seq){/*push into stack*/
	stack[stack_ptr]=tree[i];
	stack_ptr++;
    }
    else{
	child1=stack[stack_ptr-1];
	child2=stack[stack_ptr-2];
	true_tree[child1].parent=tree[i];
	true_tree[child2].parent=tree[i];
	true_tree[child1].brother=child2;
	true_tree[child2].brother=child1;

	stack[stack_ptr-2]=tree[i];
	stack_ptr--;
    }
  }

  child1=stack[0];
  true_tree[child1].parent=0;
  true_tree[child1].brother=-1; /*no brother*/

  return(0);
}

/*Combine any 2 sites and encode the combination, shape of 
comb_seq:[num_seq][C(seq_len,2)+seq_len], for each taxon, the sites 
sequence 
is:(1,1),(1,2),(1,seq_len),(2,2),...(2,se_len),...,(seq_len.seq_len)*/
int encode_2_sites(T_char_state *comb_seq)
{ int comb_row_size,comb_row_off;
  T_char_state *code_list;
  int *code_index;
  T_char_state state,state1,state2;
  T_char_state *ptr_seq,*ptr_comb,*ptr_code;
  int i,j,k,l;
  int found;
  
  comb_row_size=seq_len*(seq_len+1)/2;
  code_list=(T_char_state *)malloc(sizeof(T_char_state)*num_seq*comb_row_size*2);
  code_index=(int *)malloc(sizeof(int)*comb_row_size);
 
  for(i=0;i<comb_row_size;i++) code_index[i]=0;

  ptr_seq=seq_matrix;
  ptr_comb=comb_seq;

  for(k=0;k<num_seq;k++){
    comb_row_off=0;
    for(i=0;i<seq_len;i++){
      state1=ptr_seq[i];
      for(j=i;j<seq_len;j++)
      {	state2=ptr_seq[j];

	/*check wether state1 and state2 appears before in code_list, 
wether state1=code_list[l][comb_row_off][0] and 
state2=code_list[l][comb_row_off][1]*/
	found=-1;
	for(l=0;l<code_index[comb_row_off];l++){
	  ptr_code=code_list+l*comb_row_size*2+comb_row_off*2;
	  if((state1==ptr_code[0]) && (state2==ptr_code[1])) {
		found=l;
		break;
          }
	}

	if(found<0){
	  /*save state1 and state2 in code_list*/
	  ptr_code=code_list+code_index[comb_row_off]*comb_row_size*2
+comb_row_off*2;
	  ptr_code[0]=state1;
	  ptr_code[1]=state2;

	  found=code_index[comb_row_off];
	  code_index[comb_row_off]++;
	}
	if (found==0) state=1;
         else state=1<<found;

	*ptr_comb=state;
	comb_row_off++;
	ptr_comb++;
       }
    }
    ptr_seq+=seq_len;
  }


#ifdef DEBUG1
printf("\ncomb_seq:\n");
ptr_comb=comb_seq;
  for(i=0;i<num_seq;i++){
    for(j=0;j<comb_row_size;j++,ptr_comb++)
     printf("%d,",*ptr_comb);
   printf("\n");
 }
#endif

  free(code_list);
  free(code_index);

  return(0);
}


/*decide wether the local change on tree will exceed val, if yes, 
reurn 1*/
int execede_gap(int pos, int parent,int new_taxa)
{T_char_state *child1_state,*child2_state,*parent_state;
  int i,len;
  T_char_state temp_state;

        child1_state=seq_matrix+pos*seq_len;
	child2_state=seq_matrix+new_taxa*seq_len;
        parent_state=seq_matrix+parent*seq_len;

	len=0;
	for(i=0;i<infor_site[new_taxa];i++){
		temp_state=child1_state[i] & parent_state[i];
		if (temp_state==0) 
		temp_state=(child1_state[i] | parent_state[i]) & child2_state[i];
		 else
		  temp_state=temp_state & child2_state[i];

		if(temp_state==0) len++;
	}

	return(len);
}



                
/*compute the maximum unresolved discrepandcy*/
int compue_max_undis(tNode_T tree,int tree_size,int level,int gap)
{int old_tree_len,new_tree_len;
       T_tree_node *temp_true_tree;
	tNode_T		temp_tree;
	int i,j;
	T_char_state *ptr1,*ptr2;	
	int new_states,mini,max_mini;

	temp_tree=temp_memory.tree1;
  old_tree_len=fitch_pass(tree,tree_size);


/*	if(level> 2*num_seq/3) return(old_tree_len);
*/
        /*change the per-order encoding to true tree*/
        temp_true_tree=temp_memory.temp_true_tree;
        encode_to_true_tree(tree,tree_size,temp_true_tree);

	ptr2=sum_state+(level-1)*seq_len;  /*states which appears before*/
	max_mini=0;

	for(i=num_seq-1;i>level+2;i--){
	  /*compute the # of new states in the i_th taxon*/
	  ptr1=seq_matrix+i*seq_len;
	  new_states=0;
	  for(j=0;j<seq_len;j++)
	    if((ptr1[j] & ptr2[j])==0) new_states++;

	  mini=MAX_COST;
	  for(j=1;j<tree_size;j++){
		memcpy(temp_tree,tree,j);
		temp_tree[j]=num_seq+i-2;
		temp_tree[j+1]=i;
		memcpy(temp_tree+j+2,tree+j,tree_size-j);
		new_tree_len=fitch_pass(temp_tree,tree_size+2);
		if(mini> new_tree_len) mini=new_tree_len;
	}
	mini=mini-new_states;
	
	if(mini>gap) return(gap+1);
	if(max_mini<mini) max_mini=mini;
	}

  return(max_mini);
}

/*Expand par_sol: inserting all of its children to the frontier[id] if
they pass the lowerbound test, use the amortized constant algorithm to 
compute the cost of new tree*/
int expand_cons_fit(struct hNode *sel_hnode, int id)
{ 
	int i,j,par_seq_len;
	tNode_T par_sol,temp_tree;
	int pos,tree_size;
	int old_low,lowerbound,upperbound,taxa;
	struct hNode *heap;
	int h_id,p_id,par_score,internal_id;
	SOL_LIST_T *temp_list;
	int ival,delta,tree_len,gap;
	int *p,*b;
	T_char_state temp_state,temp_state1;
	T_char_state *ptr_seq1,*ptr_seq2,*ptr_seq3;

	par_sol=sel_hnode->par_sol;
	taxa=id+thresh_level;
	old_low=sel_hnode->lowerbound-lower_bound2[taxa-1]-infor_len[taxa-2]+
lower_bound2[taxa]+infor_len[taxa-1];
	tree_size=frontier[id-1].tree_size;
	heap=frontier[id].heap;
	internal_id=num_seq+taxa-3;
	par_seq_len=infor_site[taxa-1];	

	p=temp_memory.p;
	b=temp_memory.b;

	/*Pass par_sol in post-order*/
 	tree_len=new_fitch_pass_2(par_sol,taxa-1,p,b);
printf("tree_len:%d\n",tree_len);

	/*Pass par_sol in pre-order*/
	ptr_seq1=seq_matrix+par_sol[0]*seq_len;
	ptr_seq2=seq_matrix+(par_sol[1]+2*num_seq)*seq_len;
	for(j=0;j<par_seq_len;j++){ ptr_seq2[j]=ptr_seq1[j];
	}
	
	for(i=2;i<tree_size;i++){
	  ptr_seq1=seq_matrix+b[par_sol[i]]*seq_len;
	  ptr_seq2=seq_matrix+(p[par_sol[i]]+2*num_seq)*seq_len;
	  ptr_seq3=seq_matrix+(par_sol[i]+2*num_seq)*seq_len;
	
	  for(j=0;j<par_seq_len;j++)
	  {
		ptr_seq3[j]=ptr_seq1[j] & ptr_seq2[j];
		if (ptr_seq3[j]==0) ptr_seq3[j]=ptr_seq1[j] | ptr_seq2[j];
	  }
	}
	
	old_low=tree_len+lower_bound2[taxa]+infor_len[taxa-1];
	gap=global_ub-old_low;
	/*for each new tree inserting before "pos", compute delta of 
cost*/
	ptr_seq1=seq_matrix+order[taxa-1]*seq_len;
	for(pos=1;pos<tree_size;pos++){
		ptr_seq2=seq_matrix+par_sol[pos]*seq_len;
          	ptr_seq3=seq_matrix+(par_sol[pos]+2*num_seq)*seq_len;

		/*insert (internal_id,taxa-1) before pos*/
		delta=0;
		for(j=0;(j<par_seq_len);j++){
		  /*compute branch interval*/
		  temp_state=ptr_seq2[j] & ptr_seq3[j];
		  if (temp_state==0) temp_state=ptr_seq2[j] | ptr_seq3[j];
		  temp_state1= temp_state & ptr_seq1[j];
		  if (temp_state1==0) {
			delta++;
			if (delta >gap) break;
		  }

		}
		/*lowerbound=old_low+delta;*/
/*		
lowerbound=tree_len+delta+lower_bound2[taxa]+infor_len[taxa-1];	
*/
	printf("percentage:%d\n",j*100/par_seq_len);
		if(j<par_seq_len) continue;
		lowerbound=old_low+delta;
printf("pos:%d,lowerbound:%d\n",pos,lowerbound);
		if(lowerbound <= global_ub) {
                /*insert (internal_id,taxa-1) before pos*/
                temp_tree=heap[frontier[id].heap_size].par_sol;   
                memcpy(temp_tree,par_sol,pos);
                temp_tree[pos]=internal_id;
                temp_tree[pos+1]=order[taxa-1];
                memcpy(temp_tree+pos+2,par_sol+pos,tree_size-pos);

			if(taxa==num_seq){
				if(lowerbound ==global_ub){
					/*insert par_sol into the list of optimal solution:opt_list*/
					temp_list=(SOL_LIST_T *)malloc(sizeof(SOL_LIST_T));
					assert(temp_memory.tree_pool_ptr < POOL_SIZE);
					/*		if(temp_memory.tree_pool_ptr >= POOL_SIZE){
							   p_error(__FILE__,__LINE__,"Tree Pool Overflow");
							}
					*/
					temp_list->tree=temp_memory.tree_pool[temp_memory.tree_pool_ptr];
					temp_memory.tree_pool_ptr++;
					memcpy(temp_list->tree,temp_tree,tree_size+2);

					temp_list->next=opt_list;
					opt_list=temp_list;
				}
				else {
					global_ub=lowerbound;
					/*Delete the old  opt_list*/
					while(opt_list){
						temp_list=opt_list->next;
						free(opt_list);
						opt_list=temp_list;
					}
					/*insert temp_tree into opt_list*/
					opt_list=(SOL_LIST_T *)malloc(sizeof(SOL_LIST_T));
					opt_list->tree=temp_memory.tree_pool[0];
					memcpy(opt_list->tree,temp_tree,tree_size+2);
					opt_list->next=NULL;
					temp_memory.tree_pool_ptr=1;
				}
			} else{
				/*insert par_sol into the priority queue at frontier[id]*/
				enqueue(id,temp_tree,lowerbound);
			}
		}
	}


	return(0);
}

