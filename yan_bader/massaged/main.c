// Changes:
// WTJW 28/4/2004: Changed hardwired gettimeofday()-based timing to use the
// (slightly) platform-independent timing from COIL.

#include <stdio.h>
//#include <sys/time.h>
#include <timing.h>		//WTJW
#include "getopt.h"		//WTJW.  A Windows version I co-opted from www.codeproject.com.
#include "main.h"
/*#undef STAT*/

struct fNode  *frontier; /*active nodes in frontier*/
int global_ub; /*Global upperbound*/
SOL_LIST_T  *opt_list; /*The list of best solution, phylogenies with
global_ub*/


int *order;
int *infor_len;
int *infor_site;
int global_lb;
T_char_state *seq_matrix;
int num_seq;
int seq_len;
temp_memory_T temp_memory;
int thresh_level; /*From 3 to thresh_level, all the nodes are 
generated*/
int *lower_bound2; /*lowerbound of each partial tree =score of the
partial tree + lower_bound2 at this level. size:(n-3+1)*/

int *upper_bound2;/*upperbound of each partial tree =score of the
partial tree + upper_bound2 at this level. size:(n-3+1)*/
int sum_sinlegton_site; /*score of sinlegton sites*/
T_char_state *sum_state; /* shape:[num_seq*seq_len]*/

int *subtree_len;

#ifdef STAT
stat_data_T stat_data; /*the statistic data to analyze the
performance of the B&B*/
#endif

#define SPEC_NAME_LEN 80
char *spec_name;  /*name of speciese, the maxium lenght is SPEC_NAME_LEN*/

main(int argc, char **argv)
{ 
	char *in_file;
	int best_score,id;
	int *cost,*num_state;
	struct hNode *par_sol;
	int i,j,k,m,tree_size,temp_i,int_gb;
	T_char_state *ptr_state,*prev_ptr_state,*ptr_seq;
//	struct timeval start_time,end_time;
	TIMING_TIME start_time, end_time;		//WTJW
	double elapsed_time;


	process_opt(argc,argv,&in_file,&int_gb);

	/*Get sequence information from the input file*/
	get_input(in_file,&num_state);


        /*initialize subtree_len*/
        subtree_len=(int *)malloc(sizeof(int)*(2*num_seq-2));
        for(i=0;i<num_seq;i++) subtree_len[i]=0;

       /*malloc memory space for temp_memory*/
	temp_memory.p=(int *)malloc(sizeof(int)*2*num_seq);
	temp_memory.b=(int *)malloc(sizeof(int)*2*num_seq);
	temp_memory.stack=(int *)malloc(sizeof(int)*num_seq);
	tree_size=2*num_seq-2;
        temp_memory.temp1_state=(T_char_state *)
malloc(sizeof(T_char_state)*seq_len);
        temp_memory.temp2_state=(T_char_state *)
malloc(sizeof(T_char_state)*seq_len);

	temp_memory.temp_true_tree=(T_tree_node *)
malloc(sizeof(T_tree_node)*tree_size);
	temp_memory.tree1=(tNode_T)malloc(sizeof(T_char_state)*tree_size);
	temp_memory.tree2=(tNode_T)malloc(sizeof(T_char_state)*tree_size);
	temp_memory.tree_pool=(tNode_T *)malloc(sizeof(tNode_T)*POOL_SIZE);
	for(i=0;i<POOL_SIZE;i++)
		temp_memory.tree_pool[i]=(unsigned char *)malloc(sizeof(unsigned
		char)*tree_size);
	temp_memory.true_tree_pool=(struct tree *)malloc(sizeof(struct tree)*(2*num_seq-1));
	temp_memory.dist_matrix=(int *)malloc(sizeof(int)*num_seq*num_seq);

	/*compute thresh_level,in order to balance breadth and depth of the B&B
	tree */
	thresh_level=3;
	temp_i=1;

/*	for(i=3;(i<=num_seq) && (temp_i <= num_seq*num_seq);i++){

		temp_i*=2*i-5;
		thresh_level=i;
	}
*/
	/*malloc memory space for frontier*/
	frontier=(fNode_t *)malloc(sizeof(fNode_t)*(num_seq-thresh_level+1));
	frontier[0].tree_size=2*thresh_level-2;
	frontier[0].array_size=temp_i;
	frontier[0].heap_size=temp_i;
	for(i=thresh_level+1;i<=num_seq;i++){
		frontier[i-thresh_level].tree_size=2*i-2;
		frontier[i-thresh_level].array_size=(2*i-5)*NUM_PROCESSOR;
	}
	for(i=thresh_level;i<=num_seq;i++){
		frontier[i-thresh_level].heap=(struct hNode *)malloc(
		    sizeof(struct hNode)*frontier[i-thresh_level].array_size);
		for(j=0;j<frontier[i-thresh_level].array_size;j++){
			frontier[i-thresh_level].heap[j].par_sol=(unsigned char *)malloc
				    (sizeof(unsigned char)*frontier[i-thresh_level].tree_size);
			frontier[i-thresh_level].heap[j].internal_states=(T_char_state *)
malloc(sizeof(T_char_state)*(i-2)*seq_len);
		}
	}

	/*Reorder character (column of "seq_matrix", let the most variable sites
	present first and remove those uninformative positions*/
	reorder_cite(num_state); 

	        order=(int *)malloc((2*num_seq-2)*sizeof(int));
for(i=0;i<2*num_seq-2;i++) order[i]=i;

	/*get heuristic soliution*/
	global_ub=int_gb-sum_sinlegton_site;
	if(global_ub<0){
		/*only compute the upperbound,do not save the partial
		solution,otherwise it will be repeated*/
		global_ub=greedy_fixed_order(temp_memory.tree_pool[0]);
printf("The best score by greedy method is %d,%d\n",
global_ub,sum_sinlegton_site);
		best_score=heuristic_Eck_Day();
printf("The best score by Eck_Day is %d,%d\n",
best_score,sum_sinlegton_site);
		if(best_score<global_ub) global_ub=best_score;
/*		
global_ub=fitch_pass(temp_memory.tree_pool[0],2*num_seq-2);*/

		/*     Adjusted_SK_NJ(temp_memory.dist_matrix,temp_memory.tree_pool[0],1);
		     temp_gub=fitch_pass(temp_memory.tree_pool[0],2*num_seq-2);
		     printf("The best score by Adjusted SK Neighbor Joining is %d\n",
		temp_gub);
		
		     if(global_ub > temp_gub)
		        global_ub= temp_gub;
		*/
	}
	

	/*Decided the addtion order. Reorder sequences (column of "seq_matrix"),
	let the most unlike sequence present first*/
	lower_bound2=(int *)malloc(sizeof(int)*(num_seq+1));
/*static_reorder_rank_tournament(order);
	static_reorder_max_mini(order);*/ 	
/*	static_reorder_max_mini_lb(order); 

printf("Order of addition sequence:%\n");
for(i=0;i<num_seq;i++)          
printf("%d,",order[i]);         
printf("\n");     	

for(i=0;i<2*num_seq-2;i++) order[i]=i;
*/
	best_score=static_reorder_max_mini_lb_cons(order);
	if(best_score<global_ub) global_ub=best_score;
printf("best score by static reorder:%d,%d\n",
best_score,sum_sinlegton_site);
printf("The global_ub is:%d,%d\n",global_ub,sum_sinlegton_site);
/*static_reorder_max_mini_lb_cross(order); 

/*
old_static_reorder_max_mini_lb_cross(order);
  static_reorder_max_lb(order);
*/


  /*
  order[0]=1;
  order[1]=6;
  order[2]=15;
  order[3]=19;
  order[4]=7; 
  order[5]=4;
  order[6]=23; 
  order[7]=0;
  order[8]=14; 
  order[9]=8;
  order[10]=13; 
  order[11]=3;
  order[12]=17; 
  order[13]=21;
  order[14]=22; 
  order[15]=9;
  order[16]=5;  
  order[17]=18;
  order[18]=16; 
  order[19]=2;
  order[20]=20; 
  order[21]=10;
  order[22]=11; 
  order[23]=12;
*/

printf("Order of addition sequence:%\n");
for(i=0;i<num_seq;i++)
printf("%d,",order[i]);
printf("\n");

	infor_site=(int*)malloc(sizeof(int)*seq_len);
	infor_len=(int*)malloc(sizeof(int)*seq_len);
	reorder_cite1();
/*	not reorder cite
	for(i=0;i<seq_len;i++){
	 infor_site[i]=seq_len;
	 infor_len[i]=0;
	}

*/
#ifdef DEBUG
for(i=3;i<=num_seq;i++)
  printf("infor_site[%d]:%d\n",i-1,infor_site[i-1]);
#endif

	/*compte the statex in s_1,...s_(i-1)*/
	sum_state=(T_char_state 
*)malloc(sizeof(T_char_state)*num_seq*seq_len);
	ptr_state=sum_state;
	ptr_seq=seq_matrix+order[0]*seq_len;
	for(j=0;j<seq_len;j++,ptr_state++,ptr_seq++) 
		*ptr_state=*ptr_seq;
	prev_ptr_state=sum_state;
	for(i=1;i<num_seq;i++){
	  ptr_seq=seq_matrix+order[i]*seq_len;
	  for(j=0;j<seq_len;j++,ptr_state++,prev_ptr_state++,ptr_seq++)
		*ptr_state=(*prev_ptr_state) | (*ptr_seq);
	}
		

	/*compute  lower_bound2*/
	lb_part2(order);
printf("\nlowerbound2 from thresh_level:\n");
for(i=1;i<=num_seq;i++)
  printf("%d, ",lower_bound2[i]);


 
/*lb2_part2(order);
for(i=1;i<=num_seq;i++)
  printf("%d, ",lower_bound2[i]);
*/
	/*Compute upper_bound2
	  upper_bound2=(int *)malloc(sizeof(int)*num_seq);
	  ub_part2();
	*/

	/*Branch and Bound*/
	initalize_search(int_gb);

#ifdef STAT1
	stat_data.num_of_decomposed=0;
#endif

#ifdef DEBUG1
	printf("\nThe partial tree at frontier[0] are:\n");
	tree_size=2*thresh_level-2;
	for(i=0;i<frontier[0].heap_size;i++){
		for(j=0;j<tree_size;j++)
			printf("%d",frontier[0].heap[i].par_sol[j]);
		printf("\n");
	}
#endif

//gettimeofday(&start_time,NULL);
get_elapsed_time(&start_time);		//WTJW

	if(thresh_level==num_seq){
		/*All the phylogeny have been generated and saved in frontier[0].heap*/
		global_ub=frontier[0].heap[0].lowerbound;
		extract_min(frontier[0].heap,frontier[0].heap_size,1);
	}
	else{
		/*find the frontier level*/
		for(i=num_seq-thresh_level;i>=0;i--)
			if(frontier[i].heap_size!=0) break;

		while(i>=0){
#ifdef STAT
			stat_data.num_of_decomposed++;
			if(stat_data.num_of_decomposed % 100 ==0){
				printf("\n*******\n");
				printf(" The number of decomposed nodes is %d\n",
stat_data.num_of_decomposed);
				temp_i=0;
				for(j=0;j<=num_seq-thresh_level;j++)
				   temp_i+=frontier[j].heap_size;
				printf("The number of active nodes:%d\n",
 temp_i);
				/*output the frontier
					for(j=0;j<=num_seq-thresh_level;j++){
					  printf("The number of nodes at 
level %d: %d\n",j+thresh_level,frontier[j].heap_size);
					 for(k=0;k<frontier[j].heap_size;k++){
					    par_sol=&frontier[j].heap[k];
					    printf("LB:%d   ",par_sol->lowerbound);
				            for(m=0;m<frontier[j].tree_size;m++)
						printf("%d ",par_sol->par_sol[m]);
					    printf("\n");
					  }
					}
				*/

				/*output the opt_list
				output(opt_list,num_seq,order);
				printf("\n");
				fflush(stdout);*/
			}
#endif
/*			expand(&frontier[i].heap[0],i+1);*/
			expand_cons_fit(&frontier[i].heap[0],i+1);

			dequeue(i);

			/*find the frontier level*/
			for(i=num_seq-thresh_level;i>=0;i--)
				if(frontier[i].heap_size > 0) break;
		}
	}

//gettimeofday(&end_time,NULL);
//elapsed_time=end_time.tv_sec-start_time.tv_sec+
//(end_time.tv_usec-start_time.tv_usec)/1000000.0;
//printf("\nelpased time is %f seconds\n",elapsed_time);
get_elapsed_time(&end_time);		//WTJW
printf("\nelapsed time is %.2f seconds\n", get_time_difference(&start_time, &end_time));
	/*Output the optimal solution*/
	output(opt_list,num_seq,order);

	return 0;
}

/*Need more work to tolerate not so strict format*/
/*Get input information from the input file. Input:[num_seq][seq_len], but in
seq_matrix, it is organized[num_seq][seq_len] for space locality.  This
function also encoded the states. Encode character state and return an array of 
number of states for each site,*/
int get_input(char *filename,int **num_state)
{ 
	FILE  *in_file;
	char err_msg[80];
	int i,j,k,i_rval;
	unsigned char ch;
	T_char_state state,temp_state;
	T_char_state *ptr_seq;
	unsigned char *code_list;
	int *code_index;
	T_char_state *non_singleton_state,*or_state;
        int sum;

	in_file=fopen(filename,"r");
	assert(in_file);
	/*
	  if(in_file==NULL) {
		sprintf(err_msg,"Fail to open file %s",filename);
		p_error(__FILE__,__LINE__,err_msg);
	  }
	*/

	i_rval=fscanf(in_file,"%d%d\n", &num_seq,&seq_len);
	assert(i_rval>0);
	/* if(i_rval<=0) p_error(__FILE__,__LINE__,"Read File Error");
	*/

	*num_state=(int *)malloc(seq_len*sizeof(int));
	assert(*num_state);

	code_index=*num_state;
	for(i=0;i<seq_len;i++) code_index[i]=0;

	non_singleton_state=(T_char_state *)malloc(seq_len*sizeof(T_char_state));
	for(i=0;i<seq_len;i++) non_singleton_state[i]=0;

	or_state=(T_char_state *)malloc(seq_len*sizeof(T_char_state));
	for(i=0;i<seq_len;i++) or_state[i]=0;

	seq_matrix=(T_char_state *)malloc(sizeof(T_char_state)*(2*num_seq-2)*seq_len);
	assert(seq_matrix);
	/*  if (seq_matrix==NULL) p_error(__FILE__,__LINE__,"No Enough Space"); 
	*/

	code_list=(unsigned char *)malloc(sizeof(unsigned char)*num_seq*seq_len);
	assert(code_list);

	spec_name=(char *) malloc(num_seq*SPEC_NAME_LEN);
	ptr_seq=seq_matrix;

	for(i=0;i<num_seq;i++){
		fscanf(in_file,"%s",spec_name+i*SPEC_NAME_LEN);
		for(j=0;j<seq_len;j++){
			i_rval=fscanf(in_file,"%c",&ch);
			while (i_rval && (ch==' '))
				i_rval=fscanf(in_file,"%c",&ch);
			assert(i_rval>0);

			/*check wether"ch" has appeared before, encode the ith character state to 
			be 1<<i*/
			for(k=0;k<code_index[j];k++)
				if(*(code_list+j*num_seq+k)==ch) break;
			if(k<code_index[j]) {
				/*"ch" appeared*/
				state=1 << k;
			}
			else{
				k=code_index[j];
				*(code_list+j*num_seq+k)=ch;
				state= 1 << k;
				code_index[j]++;

			}

			if((or_state[j] & state) !=0)
				non_singleton_state[j]|=state;
			or_state[j]|=state;
			*ptr_seq=state;
			ptr_seq++;
		}
	}


	/*decide wether a cite has at most 1 nonsingleton states*/
	sum_sinlegton_site=0;
	for(i=0;i<seq_len;i++){ 
		temp_state=non_singleton_state[i];
		if(temp_state > 0)
		{ /*shift right to rmeove the first bit with value "1"*/
			while((temp_state & 0x1) == 0)
				temp_state=temp_state>>1;
		         if (temp_state==1){
			   sum_sinlegton_site+=code_index[i]-1;
                           code_index[i]=1;
			}
		}
       }


#if 0   	
	/*Compute the # of non singleton states for each taxon*/
	ptr_seq=seq_matrix;
	for(i=0;i<num_seq;i++){
	  sum=0;
	  for(j=0;j<seq_len;j++,ptr_seq++)
	    if (code_index[j]!=1){
		temp_state=*ptr_seq;
		if((temp_state & non_singleton_state[j]) !=0) sum++;
	    }
	  printf("i:%d,# of non sonleton states:%d\n",sum,i);*/
	}

   for(j=0;j<seq_len;j++)
	if(code_index[j]!=1) 
printf("num_state:%d\n",code_index[j]);

#endif

	free(code_list);
	free(non_singleton_state);
	free(or_state);
	fclose(in_file);

#ifdef DEBUG1
	ptr_seq=seq_matrix;
	printf("\nThe sequence matrix before reordering:\n");
	for(i=0;i<num_seq;i++){
		for(j=0;j<seq_len;j++,ptr_seq++)
			printf("%d,",*ptr_seq);
		printf("\n");
	}
#endif

	return(0);
}


/*Output the tree in a list*/
int output(SOL_LIST_T *list, int num_seq,int *order)
{ 
	tNode_T tree;
	int j;

	printf("\nThe optimal trees with score %d are:\n",global_ub+sum_sinlegton_site);
	while(list){
		/*print tree in NEXUS form*/
		tree=list->tree;
		/*	printf("(%s,",spec_name+order[tree[0]]*SPEC_NAME_LEN);*/
		printf("(%d,",tree[0]);
		j=print_tree(tree,2,order); /*print root->lChild->lChild*/
		printf(",");
		print_tree(tree,j,order); /*print root->lChild->rChild*/
		printf(")\n");
		list=list->next;
	}

	return(0);
}


/*print the subtree with root at tree[i] in NEXUS form and return the next
avaliable position in tree*/
int print_tree(tNode_T tree,int i,int *order)
{ 
	int j,k;

	if(tree[i]<num_seq) {
		/*     printf("%s",spec_name+order[tree[i]]*SPEC_NAME_LEN);*/
		printf("%d",tree[i]);
		return(i+1);
	}

	printf("(");
	j=print_tree(tree,i+1,order);
	printf(",");
	k=print_tree(tree,j,order);
	printf(")");

	return(k);
}


/*output the trees with minimum score,save them in opt_list*/
int extract_min(struct hNode *heap,int heap_size,int self_id)
{ 
	SOL_LIST_T *temp_list;
	int l_id,r_id;
	int tree_size=2*num_seq-2;

	/*insert the tree at heap[self_id-1] into opt_list*/
	temp_list=(SOL_LIST_T *)malloc(sizeof(SOL_LIST_T));
	assert(temp_memory.tree_pool_ptr < POOL_SIZE);
	/*                if(temp_memory.tree_pool_ptr >= POOL_SIZE){
	                   p_error(__FILE__,__LINE__,"Tree Pool Overflow");
	                }
	*/

	temp_list->tree=temp_memory.tree_pool[temp_memory.tree_pool_ptr];
	temp_memory.tree_pool_ptr++;
	memcpy(temp_list->tree,heap[self_id-1].par_sol,tree_size+2);

	temp_list->next=opt_list;
	opt_list=temp_list;


	l_id=2*self_id;
	if(l_id<=heap_size)
		if(heap[l_id-1].lowerbound==global_ub)
			extract_min(heap,heap_size,l_id);

	r_id=l_id+1;
	if(r_id<=heap_size)
		if(heap[r_id-1].lowerbound==global_ub)
			extract_min(heap,heap_size,r_id);


	return(0);
}


int   process_opt(int argc,char **argv,char **in_file_name,int *int_gb)
{
	int c;
	extern char *optarg;
	extern int optind;

	*in_file_name=NULL;
	*int_gb=-1;

	while((c=getopt(argc,argv,"hf:u:"))!=EOF)
		switch (c) {
		case 'f':
			*in_file_name=optarg;
			break;
		case 'u':
			*int_gb=atoi(optarg);
			break;
		case 'h':
		default:
			/*print usage*/
			printf("Usage: mp -f input_file -u global_upperbound\n");
			exit(0);
			break;
		}

	if(*in_file_name == NULL){
		/*print usage*/
		printf("Usage: mp -f input_file -u global_upperbound\n");
		exit(0);
	}

	return(0);
}



