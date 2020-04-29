/*This file contains all the functions based on Fitch's optimization
criterion, i.e. multiple unordered character state*/

#include <stdio.h>

#include "fitch.h"


/*Use  Fitch's method to find the possible minimal state assignment for a
given tree and return the  parsimony tree length, adapted from rooted
binary tree to unrooted binary tree*/
int fitch_pass(tNode_T tree,int tree_size)
{
	int tree_len;
	int i,j,stack_ptr,child1,child2;
	int *stack;
	T_char_state  *ptr_seq;
	T_char_state  *child1_state,*child2_state,*new_node;

	stack=temp_memory.stack;
	ptr_seq=seq_matrix;

	stack_ptr=0;
	for(j=tree_size-1;j>0;j--){
		if(tree[j]<num_seq){
			/*push state into the stack*/
			stack[stack_ptr]=tree[j];
			stack_ptr++;
		} else {
			/*pop the first 2 states and compute their Farris interval*/
			child1=stack[stack_ptr-1];
			child2=stack[stack_ptr-2];
			child1_state=ptr_seq+child1*seq_len;
			child2_state=ptr_seq+child2*seq_len;
			new_node=ptr_seq+tree[j]*seq_len;
			stack[stack_ptr-2]=tree[j];
			stack_ptr--;

			tree_len=0;
			for(i=0;i<seq_len;i++){
				new_node[i]=child1_state[i] & child2_state[i];
				if(new_node[i] == 0){
					tree_len++;
					new_node[i]=child1_state[i] | 
child2_state[i];
				}
			}

			/*save the subtre lenght*/
			subtree_len[tree[j]]=tree_len+subtree_len[child1]+
subtree_len[child2];
		}
	}


	/*Decide wether the state of root is in the state set of
				stack[stack_ptr-1]*/
	child1=tree[0];
	child2=stack[0];
	child1_state=ptr_seq+child1*seq_len;
	child2_state=ptr_seq+child2*seq_len;
	tree_len=0;
	for(i=0;i<seq_len;i++){
		if ((child1_state[i] & child2_state[i])==0)
			tree_len++;
	}

	tree_len+=subtree_len[child1]+subtree_len[child2];

	return (tree_len);
}



/*new_pos: the index of the new taxon in the tree*/
int new_fitch_incremental_pass(T_tree_node *tree,int pos,int new_taxa_id)
{
	T_char_state *temp1_state,*temp2_state,*child1_state,*child2_state;
	int temp_tree_len,farris_int;
	int child1,child2;
	int i;
        int par_seq_len,extra_tree_len;

	par_seq_len=infor_site[new_taxa_id];
	extra_tree_len=infor_len[new_taxa_id]; /*tree len from 
the non-informative sites*/

	temp1_state=temp_memory.temp1_state;
        temp2_state=temp_memory.temp2_state;
	

	child1_state=seq_matrix+pos*seq_len;
	child2_state=seq_matrix+order[new_taxa_id]*seq_len;
	/*compute Farris interval*/
           farris_int=0;
           for(i=0;i<par_seq_len;i++){
                temp1_state[i]=child1_state[i] & 
child2_state[i];
                if(temp1_state[i]==0) {
                  temp1_state[i]=child1_state[i] | child2_state[i];
                  farris_int++;
                }
           }



	temp_tree_len=subtree_len[pos]+farris_int;
	child1=pos;
	child2=tree[pos].brother;
	child1=tree[pos].parent;	

	while(child1>0){
	   /*computer farris intervla between child1 and child2*/
	   child2_state=seq_matrix+child2*seq_len;
	   farris_int=0;
	   for(i=0;i<par_seq_len;i++){
		temp2_state[i]=temp1_state[i] & child2_state[i];
		if(temp2_state[i]==0) {
		  temp2_state[i]=temp1_state[i] | child2_state[i];
		  farris_int++;
		}
           }

	     temp_tree_len+=subtree_len[child2]+farris_int;
	     
	     child2=tree[child1].brother;
	     child1=tree[child1].parent;

	    /*exchange temp1_state and temp2_state such that temp2_state is the 
state of the parent*/
	    child1_state=temp1_state;
	    temp1_state=temp2_state;
	    temp2_state=child1_state;
	}

	child1=order[0];
	child1_state=seq_matrix+child1*seq_len; /*state of the pseudo 
root*/
	for(i=0;i<par_seq_len;i++)
	  if((temp1_state[i] & child1_state[i])==0)
	    temp_tree_len++;

	return (temp_tree_len+extra_tree_len);
}


/*compute the intersection of the two character state sets, maintaining 

/*compute the intersection of the two character state sets, maintaining 
the increasing order of list,Assume s1,s2 are both in increasing
order
int intersect(T_char_state_s *s1, T_char_state_s *s2,T_char_state_s
*result)
{  int i1,i2;
  
   i1=0;
   i2=0;
   result->num_of_state=0;

   while(i1<s1->num_of_state){
     while(i2<s2->num_of_state)
	if(s1->list[i1]>s2->list[i2]) i2++;
	else break;
    if(i2>=s2->num_of_state) return (0);
    if(s1->list[i1] < s2->list[i2]) i1++;
    else {
	result->list[result->num_of_state]=s1->list[i1];
	result->num_of_state=result->num_of_state+1;
        i1++;
	i2++;
    }
   }

   return(0);
}



/*Combine the union of the two non-intersected character state sets
int union_set(T_char_state_s *s1, T_char_state_s *s2,T_char_state_s
*result)
{  int i1,i2;

   i1=0;
   i2=0;
   result->num_of_state=0;

   while((i1<s1->num_of_state) && (i2 < s2->num_of_state)){
     if(s1->list[i1]<s2->list[i2]) {
	result->list[result->num_of_state]=s1->list[i1];
        result->num_of_state=result->num_of_state+1;
	i1++;
     } else if(s1->list[i1]> s2->list[i2]){
	result->list[result->num_of_state]=s2->list[i2];
	result->num_of_state=result->num_of_state+1;
	i2++;
     }
     else{
        result->list[result->num_of_state]=s2->list[i1];
        result->num_of_state=result->num_of_state+1;
	i1++;
	i2++;
     }
  }

  if(i1<s1->num_of_state){
	memcpy(&result->list[result->num_of_state],&s1->list[i1],
(s1->num_of_state-i1)*sizeof(T_char_state));
	result->num_of_state+=s1->num_of_state-i1;
  }
   
  if(i2<s2->num_of_state){
        memcpy(&result->list[result->num_of_state],&s2->list[i2],
(s2->num_of_state-i2)*sizeof(T_char_state));
        result->num_of_state+=s2->num_of_state-i2;
  }


   return(0);
}
*/




/*new_pos: the index of the new taxon in the tree*/
int fitch_incremental_pass(T_tree_node *tree,int pos,int new_taxa_id)
{
	T_char_state *temp1_state,*temp2_state,*child1_state,*child2_state;
	int temp_tree_len,farris_int;
	int child1,child2;
	int i;


	temp1_state=temp_memory.temp1_state;
        temp2_state=temp_memory.temp2_state;
	

	child1_state=seq_matrix+pos*seq_len;
	child2_state=seq_matrix+new_taxa_id*seq_len;
	/*compute Farris interval*/
           farris_int=0;
           for(i=0;i<seq_len;i++){
                temp1_state[i]=child1_state[i] & 
child2_state[i];
                if(temp1_state[i]==0) {
                  temp1_state[i]=child1_state[i] | child2_state[i];
                  farris_int++;
                }
           }



	temp_tree_len=subtree_len[pos]+farris_int;
	child1=pos;
	child2=tree[pos].brother;
	child1=tree[pos].parent;	

	while(child1>0){
	   /*computer farris intervla between child1 and child2*/
	   child2_state=seq_matrix+child2*seq_len;
	   farris_int=0;
	   for(i=0;i<seq_len;i++){
		temp2_state[i]=temp1_state[i] & child2_state[i];
		if(temp2_state[i]==0) {
		  temp2_state[i]=temp1_state[i] | child2_state[i];
		  farris_int++;
		}
           }

	     temp_tree_len+=subtree_len[child2]+farris_int;
	     
	     child2=tree[child1].brother;
	     child1=tree[child1].parent;

	    /*exchange temp1_state and temp2_state such that temp2_state is the 
state of the parent*/
	    child1_state=temp1_state;
	    temp1_state=temp2_state;
	    temp2_state=child1_state;
	}

	child1_state=seq_matrix+order[0]*seq_len; /*state of the pseudo 
root*/
	for(i=0;i<seq_len;i++)
	  if((temp1_state[i] & child1_state[i])==0)
	    temp_tree_len++;

	return (temp_tree_len);
}



int new_fitch_pass(tNode_T tree,int level)
{
	int tree_len;
	int i,j,stack_ptr,child1,child2;
	int *stack;
	T_char_state  *ptr_seq;
	T_char_state  *child1_state,*child2_state,*new_node;
	int tree_size,par_seq_len;

        tree_size=2*level-2;
	par_seq_len=infor_site[level];

	stack=temp_memory.stack;
	ptr_seq=seq_matrix;

	stack_ptr=0;
	for(j=tree_size-1;j>0;j--){
		if(tree[j]<num_seq){
			/*push state into the stack*/
			stack[stack_ptr]=tree[j];
			stack_ptr++;
		} else {
			/*pop the first 2 states and compute their Farris interval*/
			child1=stack[stack_ptr-1];
			child2=stack[stack_ptr-2];
			child1_state=ptr_seq+child1*seq_len;
			child2_state=ptr_seq+child2*seq_len;
			new_node=ptr_seq+tree[j]*seq_len;
			stack[stack_ptr-2]=tree[j];
			stack_ptr--;

			tree_len=0;
			for(i=0;i<par_seq_len;i++){
				new_node[i]=child1_state[i] & child2_state[i];
				if(new_node[i] == 0){
					tree_len++;
					new_node[i]=child1_state[i] | 
child2_state[i];
				}
			}

			/*save the subtre lenght*/
			subtree_len[tree[j]]=tree_len+subtree_len[child1]+
subtree_len[child2];
		}
	}


	/*Decide wether the state of root is in the state set of
				stack[stack_ptr-1]*/
	child1=tree[0];
	child2=stack[0];
	child1_state=ptr_seq+child1*seq_len;
	child2_state=ptr_seq+child2*seq_len;
	tree_len=0;
	for(i=0;i<par_seq_len;i++){
		if ((child1_state[i] & child2_state[i])==0)
			tree_len++;
	}

	tree_len+=subtree_len[child1]+subtree_len[child2];

 	return(tree_len);
}


/*Use parsimoy informative site,not return cost, but return
p,b*/
int new_fitch_pass_2(tNode_T tree,int level,int*p,int *b)
{
	int i,j,stack_ptr,child1,child2;
	int *stack;
	T_char_state  *ptr_seq;
	T_char_state  *child1_state,*child2_state,*new_node;
	int tree_size,par_seq_len;
	int tree_len;


        tree_size=2*level-2;
	par_seq_len=infor_site[level];
	stack=temp_memory.stack;
	ptr_seq=seq_matrix;
	tree_len=0;

	stack_ptr=0;
	for(j=tree_size-1;j>0;j--){
		if(tree[j]<num_seq){
			/*push state into the stack*/
			stack[stack_ptr]=tree[j];
			stack_ptr++;
		} else {
			/*pop the first 2 states and compute their Farris interval*/
			child1=stack[stack_ptr-1];
			child2=stack[stack_ptr-2];
	
			p[child1]=tree[j];
			b[child1]=child2;
			p[child2]=tree[j];
			b[child2]=child1;

			child1_state=ptr_seq+child1*seq_len;
			child2_state=ptr_seq+child2*seq_len;
			new_node=ptr_seq+tree[j]*seq_len;
			stack[stack_ptr-2]=tree[j];
			stack_ptr--;

			for(i=0;i<par_seq_len;i++){
				new_node[i]=child1_state[i] & child2_state[i];
				if(new_node[i] == 0){
					new_node[i]=child1_state[i] | 
child2_state[i];
					tree_len++;
				}
			}

		}
	}

	        child1=tree[0];
        child2=stack[0];
        child1_state=ptr_seq+child1*seq_len;
        child2_state=ptr_seq+child2*seq_len;
        for(i=0;i<par_seq_len;i++)
                if ((child1_state[i] & child2_state[i])==0)
                        tree_len++;


 	return(tree_len);
}


/*Use  Fitch's method to find the possible minimal state assignment for a
given tree and return the  parsimony tree length, adapted from rooted
binary tree to unrooted binary tree,not only return cost, but 
alsoreutnr p,b*/
int fitch_pass_2(tNode_T tree,int tree_size,int *p,int *b)
{
	int tree_len;
	int i,j,stack_ptr,child1,child2;
	int *stack;
	T_char_state  *ptr_seq;
	T_char_state  *child1_state,*child2_state,*new_node;

	stack=temp_memory.stack;
	ptr_seq=seq_matrix;

	stack_ptr=0;
	tree_len=0;
	for(j=tree_size-1;j>0;j--){
		if(tree[j]<num_seq){
			/*push state into the stack*/
			stack[stack_ptr]=tree[j];
			stack_ptr++;
		} else {
			/*pop the first 2 states and compute their Farris interval*/
			child1=stack[stack_ptr-1];
			child2=stack[stack_ptr-2];

			p[child1]=tree[j];
                        b[child1]=child2;
                        p[child2]=tree[j];
                        b[child2]=child1;

			child1_state=ptr_seq+child1*seq_len;
			child2_state=ptr_seq+child2*seq_len;
			new_node=ptr_seq+tree[j]*seq_len;
			stack[stack_ptr-2]=tree[j];
			stack_ptr--;

			for(i=0;i<seq_len;i++){
				new_node[i]=child1_state[i] & child2_state[i];
				if(new_node[i] == 0){
					tree_len++;
					new_node[i]=child1_state[i] | 
child2_state[i];
				}
			}

		}
	}


	/*Decide wether the state of root is in the state set of
				stack[stack_ptr-1]*/
	child1=tree[0];
	child2=stack[0];
	child1_state=ptr_seq+child1*seq_len;
	child2_state=ptr_seq+child2*seq_len;
	for(i=0;i<seq_len;i++){
		if ((child1_state[i] & child2_state[i])==0)
			tree_len++;
	}


	return (tree_len);
}


/*compute tree length, save internal node states, only consider for those
parsimonious informative sites*/
int new_fitch_pass_3(tNode_T tree,int level)
{
	int tree_len;
	int i,j,stack_ptr,child1,child2;
	int *stack;
	T_char_state  *ptr_seq;
	T_char_state  *child1_state,*child2_state,*new_node;
	int tree_size,par_seq_len;

        tree_size=2*level-2;
	par_seq_len=infor_site[level];

	stack=temp_memory.stack;
	ptr_seq=seq_matrix;

	stack_ptr=0;
	tree_len=0;
	for(j=tree_size-1;j>0;j--){
		if(tree[j]<num_seq){
			/*push state into the stack*/
			stack[stack_ptr]=tree[j];
			stack_ptr++;
		} else {
			/*pop the first 2 states and compute their Farris interval*/
			child1=stack[stack_ptr-1];
			child2=stack[stack_ptr-2];
			child1_state=ptr_seq+child1*seq_len;
			child2_state=ptr_seq+child2*seq_len;
			new_node=ptr_seq+tree[j]*seq_len;
			stack[stack_ptr-2]=tree[j];
			stack_ptr--;

			for(i=0;i<par_seq_len;i++){
				new_node[i]=child1_state[i] & child2_state[i];
				if(new_node[i] == 0){
					tree_len++;
					new_node[i]=child1_state[i] | 
child2_state[i];
				}
			}

		}
	}


	/*Decide wether the state of root is in the state set of
				stack[stack_ptr-1]*/
	child1=tree[0];
	child2=stack[0];
	child1_state=ptr_seq+child1*seq_len;
	child2_state=ptr_seq+child2*seq_len;
	for(i=0;i<par_seq_len;i++){
		if ((child1_state[i] & child2_state[i])==0)
			tree_len++;
	}


 	return(tree_len);
}
