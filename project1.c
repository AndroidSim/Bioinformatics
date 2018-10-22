/*	CS 626: Project 1
	program to find optimal global and local alignments of pairs of sequences and three
	sequences

	written by Andrew Smith
	started on 3/6/01
*/

#include "project1.h"

/* main */

main(int argc, char *argv[])
{
	int i,l_qseq,l_dbseq,n_seqs;
	int *q_seq,*db_seq;
	char *sequence;

	/* tell the structure pointers to point to right structure */
	seq=&seq_data;
	sfxn=&sfxn_data;
	DP=&DP_data;
	zs=&zs_data;
	file=&file_data;
	acids=&aas_data;
	op=&op_data;
	MA=&MA_data;

	/* read command line and set options */
	/* default */
	op->local=0;				/* global alignement */
	op->n_shuffles=100;			/* number random sequences = 100 (the more the better, but
								   the more the computation cost) */
	op->mult_align=0;			/* multiple alignment turned off */
	command_line(argc,argv,op); 

	/* open project files for reading and writing */
	open_files(file);

	/* read into memory scoring matrix (BLOSUM 50) */
	if(!read_sm(sfxn,file->f_sm))
	{
		printf("unable to read scoring matrix\n");
		exit(1);
	}

	/* if op->mult_align =1, then the muliple alignment of three sequences is to be performed */
	if (op->mult_align)
	{
		/* find optimal global alignment for three sequences using a fxn of the same scoring scheme 
	as for pair alignment 
	*/
		op->local=0; /* only global multiple alignment done */
		three_seq_alignment(MA);
	}
	else /* or do pairwise alignments */
	{

	/* do an initial reading of sequence file to get the number of sequences so as to 
	   allocate space for matrix T(i,j), which stores the optimal scores for aligning
	   sequence i into sequence j, and space for z_score matrix of the same size */
	seq->num_seq=0;
	while (read_seq(seq,file->f_seq,2))
	{		
	}
	/* allocate space */
	n_seqs=seq->num_seq;
	DP->T_i_j=malloc((n_seqs*n_seqs)*sizeof(signed int));
	if (!DP->T_i_j)
	{
		printf("can not allocate memory for T(i,j) matrix\n");
		exit(1);
	}
	zs->z_score=malloc((n_seqs*n_seqs)*sizeof(float));
	if (!zs->z_score)
	{
		printf("can not allocate memory for z_score matrix\n");
		exit(1);
	}

	/* loop over sequences twice (once for using sequences as queries, once for using sequences
	as databases) to get optimal scores for all pairs of sequences i,j = T(i,j) */
	seq->length_seq=0;
	/* loop over sequence file so each sequence serves as a query sequence */
	seq->end_qseqs=0; /* the end of query sequence file has not been reached */
/*	printf("reading query sequences: \n"); */
	while (read_seq(seq,file->f_qseq,0))
	{
		l_qseq=seq->length_seq;
		sequence=seq->sequence;
		/* allocate memory for query sequence */
		q_seq=malloc(l_qseq*sizeof(int));
		if (!q_seq)
		{
			printf("can not allocate memory for integer representation of query sequence\n");
			exit(1);
		}
		/* convert query sequence from an array of chars to an array of ints for ease in
		   accessing values from scoring matrix */
		for (i=0; i<l_qseq; ++i)
		{
			q_seq[i]=aachar2integer(sequence[i],acids);
			/* debugging */
/*			printf("%d ",q_seq[i]); */
		} 
/*		printf("\n");
		printf("length of seq = %d\n\n",l_qseq); */

		seq->end_dbseqs=0; /* the end of database sequence file has not been reached */
/*		printf("reading database sequences: \n"); */
		while (read_seq(seq,file->f_dbseq,1))
		{
			l_dbseq=seq->length_seq;
			sequence=seq->sequence;
			/* allocate memory for database sequence */
			db_seq=malloc(l_dbseq*sizeof(int));
			if (!q_seq)
			{
				printf("can not allocate memory for integer representation of database sequence\n");
				exit(1);
			}
			/* convert database sequence for an array of chars to an array of ints for ease
			   int accessing values from scoring matrix */
			for (i=0; i<l_dbseq; ++i)
			{
				db_seq[i]=aachar2integer(sequence[i],acids);
				/* debugging */
/*				printf("%d ",db_seq[i]); */
			} 
/*			printf("\n");
			printf("length of seq = %d\n\n",l_dbseq); */

			/* pass sequences and scoring matrix to dynamic programming function (for both global and
			   local alignments) and calculate T(i,j) */
			dynamic_program(q_seq,db_seq,l_qseq,l_dbseq,0,0);
			
			/* randomly shuffle query sequence and calculate z-score */
			get_z_score(zs,DP,op,q_seq,db_seq,l_qseq,l_dbseq);

			/* check and see if end of database file; if so, break while read sequences */
			if (seq->end_dbseqs==1) break;
/*			printf("reading next database seq:\n"); */
		}

		/* check and see if end of query file; if so, break while read sequences */
		if (seq->end_qseqs==1) break;
/*		printf("reading next query seq:\n"); */
	}

	/* print T(i,j) matrix and z_scores to the file "scores.txt" */
	prnt_T(DP,seq,file->f_scores,op);
	fprintf(file->f_scores,"\n\n");
	prnt_Z(zs,DP,seq,file->f_scores,op);
	} /* end of else part of if multiple alignment */
	
	printf("\nthis is the end...my only friend, the end\n");
}
/* end of main */

int read_sm(s_fxn *sfxn, FILE *sm_file) 
{
	/* put scoring fxn (Blosum50 matrix) into sfxn->scoring matrix read from sm_file*/
	char *aacharacters; 
	char row_aa;
	int i,j,n,naas,eof;
	int *aaintegers; 
	signed int *scoring_matrix; 

	scoring_matrix=sfxn->scoring_matrix;
		
	/* variable initialization */
	naas=0;

	eof=fscanf(sm_file,"%d\n",&naas);
	if(eof==EOF)
	{
		rewind(sm_file);
		eof=fscanf(sm_file,"%d\n",&naas);
		if (eof==EOF)
		{
			printf("there is nothing in scoring matrix file\n");
			return(FALSE);
		}
		else
		{
			rewind(sm_file);
		}
	}

	/* add 1 to number of amino acids to account for gap "aa" */
	naas++;
	acids->naas=naas;

	/* debugging */
/*	printf("naas = %d\n",naas); */

	/* allocate memory for storing the amino acid types or identities and their integer
	   representation, size = naas + 1 ( the plus 1 for gap "aa") */
	acids->aacharacters=malloc(naas*sizeof(char));
	if (!acids->aacharacters)
	{
		printf("can not allocate memory for storing amino acid characters\n");
		return(FALSE);
	}
	acids->aaintegers=malloc(naas*sizeof(int));
	if (!acids->aaintegers)
	{
		printf("can not allocate memory for storing integer representatives of amino acid characters\n");
		return(FALSE);
	}

	/* read into memory the amino acid types from scoring matrix file */
	for (i=0; i<naas; ++i)
	{
		if (i<naas-1)
		{
			fscanf(sm_file,"%c",acids->aacharacters+i);
		}
		if (i==naas-1)
		{
			/* set last character to gap character == '-' */
			*(acids->aacharacters+i)='-';
		}
		acids->aaintegers[i]=i;
	}
	fscanf(sm_file,"\n");
	
	/* debugging */
	aacharacters=acids->aacharacters;
	aaintegers=acids->aaintegers;
/*	printf("types of aa used in characters:\n");
	for (i=0; i<naas; ++i)
	{
		printf("%c ",aacharacters[i]);
	}
	printf("\n");
	printf("types of aa used in integers:\n");
	for (i=0; i<naas; ++i)
	{
		printf("%d ",aaintegers[i]);
	}
	printf("\n"); */

	/* allocate memory for entire scoring matrix */
	sfxn->scoring_matrix=malloc(((naas-1)*(naas-1))*sizeof(signed int));
	if (!sfxn->scoring_matrix)
	{
		printf("can not allocate memory for scoring matrix\n");
		return(FALSE);
	}

	/* read in actual scoring matrix */
	for (i=0; i<(naas-1); ++i)
	{
		for (j=0; j<(naas-1); ++j)
		{
			n=i*(naas-1)+j;
			fscanf(sm_file,"%d",sfxn->scoring_matrix+n);
		}
		fscanf(sm_file," :%c\n",&row_aa);
	}

	/* debugging */
	scoring_matrix=sfxn->scoring_matrix;
/*	printf("scoring matrix:\n");
	for (i=0; i<(naas-1); ++i)
	{
		for (j=0; j<(naas-1); ++j)
		{
			n=i*(naas-1)+j;
			printf("%d\t",*(scoring_matrix+n));
		}
		printf("\n");
	} */

	/* set constant gap penalty equal to -10 */
	sfxn->gap_penalty=-10;
	
	/* debugging */
/*	printf("gap_penalty = %d\n",sfxn->gap_penalty); */

	/* close scoring matrix file because no longer need it */
	fclose(sm_file);

	return(TRUE);
}

int read_seq(sqs *seq, FILE *seq_file, int qordb)
{
	char ch;
	int i,j,n_seq,i_seq,l_seq,eof;
	
	n_seq=seq->num_seq;
	l_seq=seq->length_seq;
	
	eof=fscanf(seq_file,"%d.\n",&i_seq);
	if(eof==EOF)
	{
		rewind(seq_file);
		eof=fscanf(seq_file,"%d.\n",&i_seq);
/*		if (eof==EOF)
		{
			printf("there is nothing in sequence file\n");
			return(FALSE);
		}
		else
		{
			rewind(seq_file);
		}*/
	}
	if (qordb==0) seq->i_qseq=i_seq;
	if (qordb==1) seq->j_dbseq=i_seq;
	if (qordb==2) seq->num_seq++;
	if (qordb==3)
	{
		/* this one for multiple alignment */
		seq->num_seq++;
		seq->i_seq=i_seq;
	}
	/* debugging */
/*	printf("the number of sequence = %d\n",i_seq); */

	/* start by allocating 1 char size to sequence array */
	seq->sequence=malloc(1*sizeof(char));
	if (!seq->sequence)
	{
		printf("can not allocate memory for sequence\n");
		return(FALSE);
	}
	
	i=1;
	while((ch=getc(seq_file))!=EOF)
	{
		if (ch=='\n')
		{
			ch=getc(seq_file);
			if (isalnum(ch))
			{
				ungetc(ch,seq_file);
				continue;
			}
			else
			{
				seq->sequence=realloc(seq->sequence,--i);
				l_seq=i;
				seq->length_seq=l_seq;
				continue; 
/*				break; */
			}
		}
		if (isdigit(ch))
		{
			ungetc(ch,seq_file);
/*			fscanf(seq_file,"%d.\n",&i_seq); */
			seq->sequence=realloc(seq->sequence,--i);
			l_seq=i;
			seq->length_seq=l_seq;
			++n_seq; 
			break; 
		}
		if (isalpha(ch))
		{
			seq->sequence[i-1]=ch;
			++i;
			seq->sequence=realloc(seq->sequence,i);
		}
	}

	/* debugging */
/*	for (j=0; j<l_seq; ++j)
	{
		printf("%c",seq->sequence[j]);
	}
	printf("\n");
/*	printf("l_seq = %d\n",l_seq); */

	if (ch==EOF)
	{
		/* if qordb==0, then reading query file; if qordb==1, then reading database file */
		/* in this case, query file == database file */
		if (qordb==0)
		{
			seq->end_qseqs=1;
		}
		if (qordb==1)
		{
			seq->end_dbseqs=1;
		}
		/* if qordb==2, then just reading sequence file to get number of seqs and such */
		if (qordb==2)
		{
			fclose(seq_file);
			return(FALSE);
		}
		/* if qordb==3, then reading sequences for multiple alignment (need different 
		combination of data than 2,1, or 0  */
		if (qordb==3)
		{
			seq->end_seqs=1;
		}
	} 

	return(TRUE);
}

int aachar2integer(char aa, amino *acids)
{
	int i,naas;
	int *aaintegers;
	char *aacharacters;

	aacharacters=acids->aacharacters;
	aaintegers=acids->aaintegers;
	naas=acids->naas;

	for (i=0; i<naas; ++i)
	{
		if (aa==aacharacters[i]) break;
	}
	return(aaintegers[i]);
}

char integer2aachar(int integer, amino *acids)
{
	int i,naas;
	int *aaintegers;
	char *aacharacters;

	aacharacters=acids->aacharacters;
	aaintegers=acids->aaintegers;
	naas=acids->naas;

	for (i=0; i<naas; ++i)
	{
		if (integer==aaintegers[i]) break;
	}
	return(aacharacters[i]);
}

void open_files(p_files *file)
{
	strcpy(file->f_sm_name,"Blosum50.txt");
	file->f_sm=fopen(file->f_sm_name,"r");
	if (file->f_sm==NULL)
	{
		printf("can not open file %s\n",file->f_sm_name);
		exit(1);
	}
	strcpy(file->f_seq_name,"sequences.txt");
	file->f_seq=fopen(file->f_seq_name,"r");
	if (file->f_seq==NULL)
	{
		printf("can not open file %s\n",file->f_seq_name);
		exit(1);
	}
	strcpy(file->f_qseq_name,"sequences.txt");
	file->f_qseq=fopen(file->f_qseq_name,"r");
	if (file->f_qseq==NULL)
	{
		printf("can not open file %s\n",file->f_qseq_name);
		exit(1);
	}
	strcpy(file->f_dbseq_name,"sequences.txt");
	file->f_dbseq=fopen(file->f_dbseq_name,"r");
	if (file->f_dbseq==NULL)
	{
		printf("can not open file %s\n",file->f_dbseq_name);
		exit(1);
	}
	strcpy(file->f_scores_name,"scores.txt");
	file->f_scores=fopen(file->f_scores_name,"w");
	if (file->f_scores==NULL)
	{
		printf("can not open file %s\n",file->f_scores_name);
		exit(1);
	}
}

void dynamic_program(int *q_seq,int *db_seq,int l_qseq,int l_dbseq,int rs_or_zs,int twod_or_3d)
{
	/* if rs_or_zs==0, then dynamic programming is for just for alignment of two sequences
	   and regular scoring (rs)
	   if rs_or_zs==1, then dynamic programming is for generation of z_score (zs)
	   get_optimal_alignment does not need to be performed and get_optimal_score needs to 
	   be modified */

	/* allcate memory for and initialize dynamic matrix */
	if (!alloc_init_matrix(DP,sfxn,op,l_qseq,l_dbseq))
	{
		printf("can not allocate and initialize dynamic matrix\n");
		exit(1);
	}

	/* create dynamic matrix according to recursion formula involving the scoring matrix, 
	   which assigns a score for matching an amino acid to an amino acid, and constant gap 
	   penalty for two possible combinations, which are a gap in query or gap in database
	   sequence, either for global or local alignment */
	create_matrix(DP,sfxn,acids,op,seq,q_seq,db_seq,l_qseq,l_dbseq);

	/* after creating dynamic matrix, get optimal score either for global or local 
	   alignments */
	if (twod_or_3d==0) get_optimal_score(DP,op,seq,rs_or_zs);

	/* trace back optimal path through dynamic matrix to get optimal alignment */
	if (rs_or_zs==0 && twod_or_3d==0)
	{
		get_optimal_alignment(DP,acids,q_seq,db_seq,l_qseq,l_dbseq);
		prnt_alignment(DP,acids,seq,file->f_scores);
/*		free(DP->alignment); */
	}
}

int alloc_init_matrix(Dynamic *DP,s_fxn *sfxn,opshuns *op,int l_qseq,int l_dbseq)
{
	int i,j,n,n_rows,n_cols,gap_penalty;
	signed int *dynamic_matrix;
	int *trace_back;

	/* allocate memory for dynamic matrix (l_qseq+1)*(l_dbseq+1), added one to both the
	   length of th query sequence (l_qseq) and the length of the database sequence (l_dbseq)
	   to account for matching gaps at beginning of both seqs and throughout both seqs, also 
	   for nucleation center for start site for creation at top right corner */
	/* for long sequences this may be time consuming as i am allocating memory for the 
	   dynamic matrix for every pair of sequences */
	n_rows=l_qseq+1;
	DP->n_rows=n_rows;
	n_cols=l_dbseq+1;
	DP->n_cols=n_cols;

	DP->dynamic_matrix=malloc((n_rows*n_cols)*sizeof(signed int));
	if (!DP->dynamic_matrix)
	{
		printf("can not allocate memory for dynamic matrix\n");
		return(FALSE);
	}

	/* also allocate memory for trace back matrix (size = same a dynamic matrix) */
	DP->trace_back=malloc((n_rows*n_cols)*sizeof(signed int));
	if (!DP->trace_back)
	{
		printf("can not allocate memory for trace_back matrix\n");
		return(FALSE);
	}

	dynamic_matrix=DP->dynamic_matrix;
	gap_penalty=sfxn->gap_penalty;
	trace_back=DP->trace_back;

	/* initialize elements of dynamic matrix to zero because will be adding values to
	   previous elements. 
	   the elements of the first column and row are initialized by multiples of the gap 
	   penalty depending of the row or column index of that element. in other words, the 
	   value of that element in the first row or column is the sum of the value for the 
	   previous element and the gap_penalty */
	/* initialize elements of trace back matrix to 3, which essentially means null or end of
	   alignment (especially for local alignments), the other path options are the following:
	   0 means the score for the corresponding element in dynamic matrix came for the element
	     to its left
	   1 means the score for the element came from the diagonal element
	   2 means the score came from the element above */
	   
	for (i=0; i<n_rows; ++i)
	{
		for (j=0; j<n_cols; ++j)
		{
			n=i*(n_cols)+j;
			*(dynamic_matrix+n)=0;
			*(trace_back+n)=3;
		}
	}
	if (!op->local)
	{
		for (j=1; j<n_cols; ++j)
		{
	/*		*(dynamic_matrix+j)=*(dynamic_matrix+(j-1))+gap_penalty; */
			*(dynamic_matrix+j)=j*gap_penalty;
		}
		for (i=1; i<n_rows; ++i)
		{
	/*		*(dynamic_matrix+i)=*(dynamic_matrix+(i-1))+gap_penalty; */
			*(dynamic_matrix+i*(n_cols))=i*gap_penalty;
		}
	}
	for (j=1; j<n_cols; ++j)
	{
		*(trace_back+j)=0;
	}
	for (i=1; i<n_rows; ++i)
	{
		*(trace_back+i*(n_cols))=2;
	}

	/* debugging */
/*	printf("initialized dynamic matrix: \n");
	for (i=0; i<n_rows; ++i)
	{
		for (j=0; j<n_cols; ++j)
		{
			n=i*(n_cols)+j;
			printf("%d ",*(dynamic_matrix+n));
		}
		printf("\n");
	} */

	return(TRUE);
}

void create_matrix(Dynamic *DP,s_fxn *sfxn,amino *acids,opshuns *op,sqs *seq,int *q_seq,
				   int *db_seq,int l_qseq,int l_dbseq)
{
	int i,j,n,n_rows,n_cols,alpha,beta,current,diag,left,top,naas,trace;
	signed int *scoring_matrix,*dynamic_matrix;
	signed int gap_penalty,current_score,diag_score,left_score,top_score;
	int *trace_back;

	/* create the dynamic matrix according to the following recursion formula:
	   
	   dm(i,j) = the element of the dynamic matrix of row i and column j
       sm(alpha,beta) = the element of the scoring matrix for matcing amino acid alpha and
	                    amino acid beta
       alpha(i) = amino acid alpha at site i of query sequence
	   beta(j) = amino acid beta at site j of database sequence
       gp = gap penalty

	   if global,
		dm(i,j) = max ( dm(i-1,j-1) + sm(alpha(i),beta(j)),  path along diagonal
						dm(i-1,j) + gp,                      path from top
						dm(i,j-1) + gp,                      path from left
						)
	   if local,
		dm(i,j) = max ( same as above
						same as above
						same as above
						0
					  )

       if there are degeneracies in the three possibilities for the maximization with the 
	   diagonal component then the diagonal component has precedence over the the two gap
	   option.  if the two gap options are the same then the precedence is for gap in the
	   shortest of the two sequences.
	*/
	/* trace_back contructed as follows:
	   if the score for the element came from the element to the left plug gap_penalty, then
	   trace back matrix element = 0
	   if score came from diagonal, then trace back matrix element = 1
	   if score came from above, then trace back matrix element = 2
	*/

	n_rows=DP->n_rows;
	n_cols=DP->n_cols;
	dynamic_matrix=DP->dynamic_matrix;
	trace_back=DP->trace_back;
	scoring_matrix=sfxn->scoring_matrix;
	naas=acids->naas;
	gap_penalty=sfxn->gap_penalty;

	for (i=1; i<n_rows; ++i)
	{
		for (j=1; j<n_cols; ++j)
		{
			/* get amino acids at site i of query sequence and site j of database sequence */
			alpha=*(q_seq+(i-1));
			beta=*(db_seq+(j-1));

			current=i*n_cols+j;
			diag=(i-1)*n_cols+(j-1);
			left=i*n_cols+(j-1);
			top=(i-1)*n_cols+j;

			/* assign current (default) score for path along diagonal  and trace = 1 */
			diag_score=*(dynamic_matrix+diag)+*(scoring_matrix+alpha*(naas-1)+beta);
			current_score=diag_score;
			trace =1;

			left_score=*(dynamic_matrix+left)+gap_penalty;
			top_score=*(dynamic_matrix+top)+gap_penalty;
			if (left_score>diag_score && left_score>top_score)
			{
				current_score=left_score;
				trace=0;
			}
			if (top_score>diag_score && top_score>left_score)
			{
				current_score=top_score;
				trace=2;
			}
			/* the following may be unnecessary or incorrect */
			if (left_score==top_score && left_score>diag_score && top_score>diag_score)
			{
				current_score=left_score;
				trace=0;
				if (l_dbseq<l_qseq)
				{
					current_score=top_score;
					trace=2;
				}
			}

			if (op->local==1 && current_score<0)
			{
				*(dynamic_matrix+current)=0;
				*(trace_back+current)=3;
			}
			else
			{
				*(dynamic_matrix+current)=current_score;
				*(trace_back+current)=trace;
			}
		}
	}

	/* debugging */
/*	if (seq->i_qseq==3 && seq->j_dbseq==5)
	{
	fprintf(file->f_scores,"dynamic matrix:\n");
	for (i=0; i<n_rows; ++i)
	{
		for (j=0; j<n_cols; ++j)
		{
			n=i*n_cols+j;
			fprintf(file->f_scores,"%d ",*(dynamic_matrix+n));
		}
		fprintf(file->f_scores,"\n");
	}
	fprintf(file->f_scores,"trace matrix:\n");
	for (i=0; i<n_rows; ++i)
	{
		for (j=0; j<n_cols; ++j)
		{
			n=i*n_cols+j;
			fprintf(file->f_scores,"%d ",*(trace_back+n));
		}
		fprintf(file->f_scores,"\n");
	}
	} */
}

void get_optimal_score(Dynamic *DP,opshuns *op,sqs *seq,int rs_or_zs)
{
	signed int *T_i_j,*dynamic_matrix;
	signed int max,ith_pos,jth_pos;
	int i,j,n,n_rows,n_cols,i_qseq,j_dbseq,n_seq,max_position,current;

	T_i_j=DP->T_i_j;
	dynamic_matrix=DP->dynamic_matrix;
	n_rows=DP->n_rows;
	n_cols=DP->n_cols;
	i_qseq=seq->i_qseq;
	j_dbseq=seq->j_dbseq;
	n_seq=seq->num_seq;

	current=(i_qseq-1)*n_seq+(j_dbseq-1);
	/* for global alignment (op->lorg==0) optimal score is bottom-right corner of dynamic
	   matrix */
	if (op->local==0)
	{
		/* store optimal alignment in matrix T(i,j), where i indexes query sequences and j
		   indexes database sequences */
		if (rs_or_zs==0)
		{
			*(T_i_j+current)=*(dynamic_matrix+(n_rows*n_cols)-1);
			DP->optimal_score=*(dynamic_matrix+(n_rows*n_cols)-1);
			DP->optimal_i=(n_rows-1);
			DP->optimal_j=(n_cols-1);
		}
		if (rs_or_zs==1)
		{
			DP->opt_rand_score=*(dynamic_matrix+(n_rows*n_cols)-1);
		}
		/* or optimal score is the maximum value in both the bottom row and last column */
	}
	/* for local alignment (op->lorg==1) optimal score is the maximum score throughout the
	   entire dynamic matrix */
	if (op->local==1)
	{
		for (i=0; i<n_rows; ++i)
		{
			for (j=0; j<n_cols; ++j)
			{
				n=i*n_cols+j;
				max=*dynamic_matrix; /* initially, max equals top-left corner */
				if (*(dynamic_matrix+n)>max)
				{
					ith_pos=i;
					jth_pos=j;
					max_position=n;
					max=*(dynamic_matrix+n);
				}
			}
		}
		if (rs_or_zs==0) /* if for regular scoring (original query) */
		{
			*(T_i_j+current)=max;
			DP->optimal_score=max;
			DP->optimal_i=ith_pos;
			DP->optimal_j=jth_pos;
		}
		if (rs_or_zs==1) /* if for random sequences for calculating z-score */
		{
			DP->opt_rand_score=max;
		}
	}
}

void get_optimal_alignment(Dynamic *DP,amino *acids,int *q_seq,int *db_seq,int l_qseq,int l_dbseq)
{
	int i,j,max_l_alignment,end,n,m,path,n_rows,n_cols;
	int l_alignment,alpha,beta,temp_q,temp_db,n_switches;
	int *trace_back;
	int **alignment;
	/* allocate alignment matrix: each column is an alignment pair match of the following types:
	   trace_back = 0 corresponds to gap in query => pair = (- aa)
	   trace_back = 1 corresponds to aa matched to aa => pair = (aa aa)
	   trace_back = 2 corresponds to gap in database => pair = (aa -)
	   max size of alignment matrix = (length of query) + (length of database) */
	max_l_alignment=(l_qseq+l_dbseq);
	DP->alignment=malloc(2*sizeof(int*));
	if (!DP->alignment)
	{
		printf("can not allocate memory for alignment matrix\n");
		exit(1);
	}
	for (i=0; i<2; ++i)
	{
		*(DP->alignment+i)=malloc(max_l_alignment*sizeof(int));
		if (!*(DP->alignment+i))
		{
			printf("can not allocate 2nd dimension of alignment matrix\n");
			exit(1);
		}
	}
	alignment=DP->alignment;
	trace_back=DP->trace_back;

	n_rows=DP->n_rows;
	n_cols=DP->n_cols;
	n=DP->optimal_i;
	m=DP->optimal_j;
	/* trace back alignment path from position of optimal score */
	l_alignment=0;
	end=FALSE;
	while (end==FALSE)
	{
		path=*(trace_back+n*(n_cols)+m);
		if (path==0) /* horizontal path */
		{
			beta=m-1;
			/* the last matched pair of alignment will be at beginning of alignment matrix */
/* query */	alignment[0][l_alignment]=aachar2integer('-',acids); /* put gap */
/* db */	alignment[1][l_alignment]=*(db_seq+beta); /* put amino acid from query */
			--m;
		}
		if (path==1)
		{
			alpha=n-1;
			beta=m-1;
			alignment[0][l_alignment]=*(q_seq+alpha);
			alignment[1][l_alignment]=*(db_seq+beta);
			--m;
			--n;
		}
		if (path==2)
		{
			alpha=n-1;
			alignment[0][l_alignment]=*(q_seq+alpha);
			alignment[1][l_alignment]=aachar2integer('-',acids);
			--n;
		}
		if (path==3)
		{
			/* end of alignment has been reached */
			alpha=n;
			beta=m;
			alignment[0][l_alignment]=*(q_seq+alpha);
			alignment[1][l_alignment]=*(db_seq+beta);
			end=TRUE;
			continue;
		}
		++l_alignment;
	}
	DP->l_alignment=l_alignment;

	/* since end of alignment is beginning of alignment matrix, need to reverse alignment
	   matrix (not needed for mult alignment */
	if ((l_alignment%2)==0)
	{
		n_switches=(l_alignment/2);
	}
	else
	{
		n_switches=(l_alignment-1)/2;
	}
	j=l_alignment-1;
	for (i=0; i<n_switches; ++i)
	{
		temp_q=alignment[0][i];
		temp_db=alignment[1][i];
		alignment[0][i]=alignment[0][j];
		alignment[1][i]=alignment[1][j];
		alignment[0][j]=temp_q;
		alignment[1][j]=temp_db;
		--j;
	}
}

void prnt_alignment(Dynamic *DP,amino *acids,sqs *seq,FILE *a_file)
{
	int i,j,k,l_alignment,n_qprnt,n_dbprnt,n_left,line,position;
	int **alignment;
	char q_name[15],db_name[15];

	l_alignment=DP->l_alignment;
	alignment=DP->alignment;

	line=80;

/*	if (op->threading)
	{
		if (seq->i_qseq==1)
		{
			strcpy(q_name,"leghemoglobin");
		}
		if (seq->j_dbseq==1)
		{
			strcpy(db_name,"leghemoglobin");
		}
		if (seq->i_qseq==2)
		{
			strcpy(q_name,"myoglobin");
		}
		if (seq->j_dbseq==2)
		{
			strcpy(db_name,"myoglobin");
		}
		fprintf(a_file,"alignment of %s sequence into %s structure:\n",q_name,db_name);
		fprintf(a_file,"query sequence is top sequence and database structure is bottom sequence:\n\n");
	}
	else
	{*/
		fprintf(a_file,"alignment for query sequence %d and database structure %d:\n",seq->i_qseq,seq->j_dbseq);
/*	} */

	if (l_alignment>line)
	{
		n_left=line;
		n_qprnt=0;
		n_dbprnt=0;
		while (n_left>0)
		{
			if ((l_alignment-n_qprnt)<=line && (l_alignment-n_dbprnt)<=line)
			{
				i=0;
				for (k=n_qprnt; k<l_alignment; ++k)
				{
					fprintf(a_file,"%c",integer2aachar(alignment[i][k],acids));
					++n_qprnt;
				}
				fprintf(a_file,"\n");
				i=1;
				for (k=n_dbprnt; k<l_alignment; ++k)
				{
					fprintf(a_file,"%c",integer2aachar(alignment[i][k],acids));
					++n_dbprnt;
				}
				fprintf(a_file,"\n\n");
			}
			else
			{
				i=0;
				position=n_qprnt;
				for (k=n_qprnt; k<position+line; ++k)
				{
					fprintf(a_file,"%c",integer2aachar(alignment[i][k],acids));
					++n_qprnt;
				}
				fprintf(a_file,"\n");
				i=1;
				position=n_dbprnt;
				for (k=n_dbprnt; k<position+line; ++k)
				{
					fprintf(a_file,"%c",integer2aachar(alignment[i][k],acids));
					++n_dbprnt;
				}
				fprintf(a_file,"\n\n");
			}
			if (n_qprnt != n_dbprnt)
			{
				printf("there is a problem with printing alignment\n");
				exit(1);
			}
			n_left=l_alignment-n_qprnt;
		}
	}
	else
	{
		for (i=0; i<2; ++i)
		{
			for (j=0; j<l_alignment; ++j)
			{
				fprintf(a_file,"%c",integer2aachar(alignment[i][j],acids));
			}
			fprintf(a_file,"\n");
		}
		fprintf(a_file,"\n\n");
	}
}

/*void prnt_alignment(Dynamic *DP,amino *acids,sqs *seq,FILE *a_file)
{
	int i,j,l_alignment;
	int **alignment;

	l_alignment=DP->l_alignment;
	alignment=DP->alignment;

	fprintf(a_file,"alignment for query sequence %d and database sequence %d\n",seq->i_qseq,seq->j_dbseq);
	for (i=0; i<2; ++i)
	{
		if (i==0)
		{
			fprintf(a_file,"query sequence:\t");
		}
		if (i==1)
		{
			fprintf(a_file,"datab sequence:\t");
		}
		for (j=0; j<l_alignment; ++j)
		{
			fprintf(a_file,"%c",integer2aachar(alignment[i][j],acids));
		}
		fprintf(a_file,"\n");
	}
	fprintf(a_file,"\n\n");
} */

void prnt_T(Dynamic *DP,sqs *seq,FILE *a_scores,opshuns *op)
{
	signed int *T_i_j;
	int i,j,n,n_seq;

	T_i_j=DP->T_i_j;
	n_seq=seq->num_seq;

	fprintf(a_scores,"T(i,j):\n");
	if (op->local==0)
	{
		fprintf(a_scores,"global alignment\n\n");
	}
	if (op->local==1)
	{
		fprintf(a_scores,"local alignment\n\n");
	}
	fprintf(a_scores,"database sequence:\t");
	for (i=0; i<n_seq; ++i)
	{
		fprintf(a_scores,"%3d\t",i+1);
	}
	fprintf(a_scores,"\n");
	fprintf(a_scores,"query sequence:\n");
	for (i=0; i<n_seq; ++i)
	{
		fprintf(a_scores,"\t%d\t\t\t",i+1);
		for (j=0; j<n_seq; ++j)
		{
			n=i*n_seq+j;
			fprintf(a_scores,"%0d\t",*(T_i_j+n));
		}
		fprintf(a_scores,"\n");
	}
}

void command_line(int argc, char *argv[],opshuns *op)
{
	int i,var;

	for (i=0; i<argc; ++i)
	{
		var=0;
		if (i+1<argc)
		{
			if (strncmp(argv[i+1],"-",1)!=0) var=1;
		}
		if (strcmp(argv[i],"-local")==0)
		{
			op->local=1;
		}
		if (strcmp(argv[i],"-random")==0)
		{
			if (var) op->n_shuffles=atoi(argv[i+1]);
		}
		if (strcmp(argv[i],"-multa")==0)
		{
			op->mult_align=1;
		}
	}
}

void get_z_score(Zscore *zs,Dynamic *DP,opshuns *op,int *q_seq,int *db_seq,int l_qseq,int l_dbseq)
{
	int i;
	int **random_seqs;

	/* generate a number of shuffled seqs (op->n_shuffles) of the same composition and length
	   of the query sequence and store them in the matrix zs->*random_seqs */
	randomize_qseq(zs,op,q_seq,l_qseq);

	/* determine the optimal scores for aligning each of the op->n_shuffles random seqs into
	   the database sequence using dynamic programming and scoring matrix */
	/* allocate memory for random sequence scores array */
	zs->rand_scores=malloc((op->n_shuffles)*sizeof(signed int));
	if (!zs->rand_scores)
	{
		printf("can not allocate memory for random scores array\n");
		exit(1);
	}

	random_seqs=zs->random_seqs;

	for (i=0; i<op->n_shuffles; ++i)
	{
		dynamic_program(random_seqs[i],db_seq,l_qseq,l_dbseq,1,0);
		*(zs->rand_scores+i)=DP->opt_rand_score;		
	}

	/* calculate the z-score (measure of statistical significance of an alignment score) */ 
	
	calc_z_score(zs,DP,op,seq);
}

void randomize_qseq(Zscore *zs,opshuns *op,int *q_seq,int l_qseq)
{
	int i,j,k,stime,rand_pos1,rand_pos2,temp;
	int *local_qseq;
	long ltime;
	float rand_num;

	/* make local copy of query sequence */
	local_qseq=malloc(l_qseq*sizeof(int));
	if (!local_qseq)
	{
		printf("can not allocate memory for local copy of query sequence for randomization\n");
		exit(1);
	}
	for (i=0; i<l_qseq; ++i)
	{
		*(local_qseq+i)=*(q_seq+i);
	}

	/* allocate memory for random sequences generate */
/*	zs->random_seqs=malloc((op->n_shuffles*l_qseq)*sizeof(int)); */
	zs->random_seqs=malloc((op->n_shuffles)*sizeof(int*));
	if (!zs->random_seqs)
	{
		printf("can not allocate memory for random sequences matrix generated from randomization\n");
		exit(1);
	}
	for (i=0; i<op->n_shuffles; ++i)
	{
		*(zs->random_seqs+i)=malloc(l_qseq*sizeof(int));
		if (!*(zs->random_seqs+i))
		{
			printf("can not allocate 2nd dimension of random sequences matrix\n");
			exit(1);
		}
	}
	/* initialize random number generator rand() with system time */
	ltime=time(NULL);
	stime=(unsigned int) ltime/2;
	srand(stime);

	/* loop over twice the length of the query sequence making random pair permutations.  
	   do this op->n_shuffles times to generate random_seqs matrix (size = n_shuffles*l_qseq) */

	for (i=0; i<op->n_shuffles; ++i)
	{
		for (j=0; j<2*l_qseq; ++j)
		{
			rand_num=rand()/((float)RAND_MAX);
			rand_pos1=(int)((l_qseq-1)*rand_num);
			/* debugging */
			if (rand_pos1<0 || rand_pos1>=l_qseq)
			{
				--j;
				continue;
/*				printf("0>rand_pos1>=l_qseq, rand_pos1 = %d, l_qseq = %d\n",rand_pos1,l_qseq);
				exit(1); */
			}
			rand_num=rand()/((float)RAND_MAX);
			rand_pos2=(int)((l_qseq-1)*rand_num);
			/* debugging */
			if (rand_pos2<0 || rand_pos2>=l_qseq)
			{
				--j;
				continue;
/*				printf("0>rand_pos2>=l_qseq, rand_pos2 = %d, l_qseq = %d\n",rand_pos2,l_qseq);
				exit(1); */
			}
			temp=*(local_qseq+rand_pos1);
			*(local_qseq+rand_pos1)=*(local_qseq+rand_pos2);
			*(local_qseq+rand_pos2)=temp;
		}
		/* copy shuffled local sequence to random sequences matrix */
		for (k=0; k<l_qseq; ++k)
		{
			zs->random_seqs[i][k]=*(local_qseq+k); 
/*			**(zs->random_seqs+i*(l_qseq)+k)=*(local_qseq+k); */ 
		}
	}

	/* debugging */
/*	printf("random seqs matrix:\n");
	for (i=0; i<op->n_shuffles; ++i)
	{
		for (j=0; j<l_qseq; ++j)
		{
			printf("%d  ",zs->random_seqs[i][j]);
		}
		printf("\n");
	} */
}

void calc_z_score(Zscore *zs,Dynamic *DP,opshuns *op,sqs *seq)
{
	int i,current,i_qseq,j_dbseq,n_seq;
	signed int Tab;
	signed int *Tr;
	float avg_Tr;
	float avg_Tr2;

	Tr=zs->rand_scores;
	i_qseq=seq->i_qseq;
	n_seq=seq->num_seq;
	j_dbseq=seq->j_dbseq;

	/* z_score calculated according to the following formula:
	   Tab = score of aligning sequence a into sequence b
	   <Tr> = the average of the scores of aligning the n_shuffled random sequences into
	          database sequence
	   <Tr2> = the average of the squares of Tr
	   <Tr>2 = the average of Tr squared

       z_score = (Tab - <Tr>)/(sqrt(<Tr2> - <Tr>2)), denominator = zs->std_dev
	*/

	/* Tab is the optimal score from regular alignment (DP->optimal_score) */
	Tab=DP->optimal_score;

	/* calculate <Tr>, zs->avg_rand_scores = <Tr> */
	avg_Tr=0.0;
	for (i=0; i<op->n_shuffles; ++i)
	{
		avg_Tr+=*(Tr+i);
	}
	avg_Tr=(avg_Tr)/(op->n_shuffles);
	zs->avg_rand_scores=avg_Tr;

	/* calculate <Tr2> */
	avg_Tr2=0.0;
	for (i=0; i<op->n_shuffles; ++i)
	{
		avg_Tr2+=((*(Tr+i))*(*(Tr+i)));
	}
	avg_Tr2=(avg_Tr2)/(op->n_shuffles);
	zs->avg_rand_scores_sqrd=avg_Tr2;

	/* calculate standard deviation (denominator of z_score) */
	zs->std_dev=(float)sqrt((double)(avg_Tr2-(avg_Tr*avg_Tr)));

	/* calculate actual z_score */
	current=(i_qseq-1)*n_seq+(j_dbseq-1);
	*(zs->z_score+current)=(Tab-avg_Tr)/(zs->std_dev);
}

void prnt_Z(Zscore *zs,Dynamic *DP,sqs *seq,FILE *a_scores,opshuns *op)
{
	float *z_score;
	int i,j,n,n_seq;

	z_score=zs->z_score;
	n_seq=seq->num_seq;

	fprintf(a_scores,"z_score matrix:\n");
	if (op->local==0)
	{
		fprintf(a_scores,"global alignment\n\n");
	}
	if (op->local==1)
	{
		fprintf(a_scores,"local alignment\n\n");
	}
	fprintf(a_scores,"database sequence:\t  ");
	for (i=0; i<n_seq; ++i)
	{
		fprintf(a_scores,"%d      ",i+1);
	}
	fprintf(a_scores,"\n");
	fprintf(a_scores,"query sequence:\n");
	for (i=0; i<n_seq; ++i)
	{
		fprintf(a_scores,"\t%d\t\t\t",i+1);
		for (j=0; j<n_seq; ++j)
		{
			n=i*n_seq+j;
			fprintf(a_scores,"%05.2f  ",*(z_score+n));
		}
		fprintf(a_scores,"\n");
	}
}

void three_seq_alignment(Mult_A *MA)
{
	int l_seq1,l_seq2,l_seq3;
	int *seq1,*seq2,*seq3;
	/* this function performs a multiple sequence alignment of three sequences using 
	   dynamic programming */

	/* get three sequences */
	get_three_sequences(MA,seq);
	seq1=MA->seq1;
	l_seq1=MA->l_seq1;
	seq2=MA->seq2;
	l_seq2=MA->l_seq2;
	seq3=MA->seq3;
	l_seq3=MA->l_seq3;

	/* allocate and initialize 3-dimensional dynamic matrix */
	setup_3d_matrix(MA,DP,seq,sfxn,acids,seq1,l_seq1,seq2,l_seq2,seq3,l_seq3);

	/* create 3-d dynamic matrix */
	create_3d_matrix(MA,sfxn,acids,seq1,l_seq1,seq2,l_seq2,seq3,l_seq3);

	/* get optimal score for multiple alignment */
	get_optimal_3d_score(MA,l_seq1,l_seq2,l_seq3);

	/* trace back 3d dynamic matrix to get multiple alignment */
	get_mult_alignment(MA,acids,seq1,l_seq1,seq2,l_seq2,seq3,l_seq3);

	/* print multiple alignement */
	prnt_mult_alignment(MA,acids,file->f_scores);
}

void get_three_sequences(Mult_A *MA,sqs *seq)
{
	int i,n,which_seq[3];

	/* sequences to be chosen are 6,7,9.  this can be generalized later */
	which_seq[0]=1;
	which_seq[1]=3;
	which_seq[2]=5;
	n=0;
	seq->end_seqs=0;
	while (read_seq(seq,file->f_seq,3))
	{
/*		if (seq->i_seq==which_seq[n])
		{
			switch(seq->i_seq) {
			case 6: */
		if (seq->i_seq==1)
		{
				MA->l_seq1=seq->length_seq;
				/* allocate memory for sequence */
				MA->seq1=malloc(MA->l_seq1*sizeof(int));
				if (!MA->seq1)
				{
					printf("can not allocate memory for integer representation of seq1\n");
					exit(1);
				}
				/* convert sequence from an array of chars to an array of ints for ease in
				   accessing values from scoring matrix */
				for (i=0; i<MA->l_seq1; ++i)
				{
					MA->seq1[i]=aachar2integer(seq->sequence[i],acids);
					/* debugging */
		/*			printf("%d ",MA->seq1[i]); */
				}
		}
/*				break; */
/*			case 7: */
		if (seq->i_seq==3)
		{
				MA->l_seq2=seq->length_seq;
				/* allocate memory for sequence */
				MA->seq2=malloc(MA->l_seq2*sizeof(int));
				if (!MA->seq2)
				{
					printf("can not allocate memory for integer representation of seq1\n");
					exit(1);
				}
				/* convert sequence from an array of chars to an array of ints for ease in
				   accessing values from scoring matrix */
				for (i=0; i<MA->l_seq2; ++i)
				{
					MA->seq2[i]=aachar2integer(seq->sequence[i],acids);
					/* debugging */
		/*			printf("%d ",MA->seq2[i]); */
				}
		}
/*				break; */
/*			case 9: */
		if (seq->i_seq==5)
		{
				MA->l_seq3=seq->length_seq;
				/* allocate memory for sequence */
				MA->seq3=malloc(MA->l_seq3*sizeof(int));
				if (!MA->seq3)
				{
					printf("can not allocate memory for integer representation of seq1\n");
					exit(1);
				}
				/* convert sequence from an array of chars to an array of ints for ease in
				   accessing values from scoring matrix */
				for (i=0; i<MA->l_seq3; ++i)
				{
					MA->seq3[i]=aachar2integer(seq->sequence[i],acids);
					/* debugging */
		/*			printf("%d ",MA->seq3[i]); */
				}
		}
		/*		break;*/
		/*	}  end of switch */
		if (seq->end_seqs==1) break;
		++n;
	}
}

void setup_3d_matrix(Mult_A *MA,Dynamic *DP,sqs *seq,s_fxn *sfxn,amino *acids,int *seq1,int l_seq1,
					 int *seq2,int l_seq2,int *seq3,int l_seq3)
{
	int i,j,l,m,n,l_1dim,l_2dim,l_3dim,n_rows,n_cols;
	int l_3seqs[3];
	int ***trace_back;
	signed int ***dynamic_matrix;
	signed int gap_penalty;
	/* allocate memory for dynamic matrix (l_qseq+1)*(l_dbseq+1), added one to both the
	   length of th query sequence (l_qseq) and the length of the database sequence (l_dbseq)
	   to account for matching seqs to a string a gaps (only one way to do this), also 
	   for nucleation center for start site for creation at top right corner */
	/* for long sequences this may be time consuming as i am allocating memory for the 
	   dynamic matrix for every pair of sequences */
	
	l_1dim=l_seq1+1;
	l_2dim=l_seq2+1;
	l_3dim=l_seq3+1;
	MA->l_1dim=l_1dim;
	MA->l_2dim=l_2dim;
	MA->l_3dim=l_3dim;

	MA->dynamic_matrix=malloc(l_1dim*sizeof(signed int**));
	if (!MA->dynamic_matrix)
	{
		printf("can not allocate memory for mult alignment dynamic matrix\n");
		exit(1);
	}
	for (i=0; i<l_1dim; ++i)
	{
/*		*(MA->dynamic_matrix+i)=malloc(l_2dim*sizeof(signed int*)); */
		MA->dynamic_matrix[i]=malloc(l_2dim*sizeof(signed int*));
		if (!MA->dynamic_matrix[i])
		{
			printf("can not allocate 2nd dimension of mult alignment dynamic matrix\n");
			exit(1);
		}
		for (j=0; j<l_2dim; ++j)
		{
/*			**(MA->dynamic_matrix+i*(l_2dim)+j)=malloc(l_3dim*sizeof(signed int)); */
			MA->dynamic_matrix[i][j]=malloc(l_3dim*sizeof(signed int));
			if (!MA->dynamic_matrix[i][j])
			{
				printf("can not allocate memory for 3rd dimension of dynamic matrix\n");
				exit(1);
			}
		}
	}
	/* also allocate memory for trace back matrix (size = same a dynamic matrix) */
	MA->trace_back=malloc(l_1dim*sizeof(int**));
	if (!MA->trace_back)
	{
		printf("can not allocate memory for 1st dimension of trace_back matrix\n");
		exit(1);
	}
	for (i=0; i<l_1dim; ++i)
	{
/*		*(MA->trace_back+i)=malloc(l_2dim*sizeof(int*)); */
		MA->trace_back[i]=malloc(l_2dim*sizeof(int*));
		if (!MA->trace_back[i])
		{
			printf("can not allocate 2nd dimension of trace_back matrix\n");
			exit(1);
		}
		for (j=0; j<l_2dim; ++j)
		{
/*			**(MA->trace_back+i*(l_2dim)+j)=malloc(l_3dim*sizeof(int)); */
			MA->trace_back[i][j]=malloc(l_3dim*sizeof(int));
			if (!MA->trace_back[i][j])
			{
				printf("can not allocate memory for 3rd dimension of trace_back matrix\n");
				exit(1);
			}
		}
	}

	dynamic_matrix=MA->dynamic_matrix;
	trace_back=MA->trace_back;
	gap_penalty=sfxn->gap_penalty;
	/* initialize entire dynamic matrix with zeros and trace_back matrix with 7s 
	   (7s are like 3s in 2d case) */
	for (l=0; l<l_1dim; ++l)
	{
		for (m=0; m<l_2dim; ++m)
		{
			for (n=0; n<l_3dim; ++n)
			{
				dynamic_matrix[l][m][n]=0;
				trace_back[l][m][n]=7;
			}
		}
	}

	/* initialize three intersecting faces of the mult alignment cube with 2d dynamic
	   programming of the sequences that make up the edges of each face
	   l = index of 1st dimension (seq1), m = index of second dim (seq2), n = index of 3rd (seq3) */
	/* the faces of the trace_back matrix use the convention of 2d alignments
	   the two conventions for path movement (2d and 3d) have overlap, ie 0 = left move in 2d
	   but some other move in 3d.  as long as we know if we are on a face or in the interior
	   of cube, it should be alright */
	dynamic_program(seq1,seq2,l_seq1,l_seq2,0,1);
	for (l=0; l<l_1dim; ++l)
	{
		for (m=0; m<l_2dim; ++m)
		{
			dynamic_matrix[l][m][0]=*(DP->dynamic_matrix+l*l_2dim+m);
			trace_back[l][m][0]=*(DP->trace_back+l*l_2dim+m);
		}
	}
	dynamic_program(seq1,seq3,l_seq1,l_seq3,0,1);
	for (l=0; l<l_1dim; ++l)
	{
		for (n=0; n<l_3dim; ++n)
		{
			dynamic_matrix[l][0][n]=*(DP->dynamic_matrix+l*l_3dim+n);
			trace_back[l][0][n]=*(DP->trace_back+l*l_3dim+n);
		}
	}
	dynamic_program(seq2,seq3,l_seq2,l_seq3,0,1);
	for (m=0; m<l_2dim; ++m)
	{
		for (n=0; n<l_3dim; ++n)
		{
			dynamic_matrix[0][m][n]=*(DP->dynamic_matrix+m*l_3dim+n);
			trace_back[0][m][n]=*(DP->trace_back+m*l_3dim+n);
		}
	}

	/* debugging */
/*	printf("trace back matrix with l=0:\n");
	for (m=0; m<l_2dim; ++m)
	{
		for (n=0; n<l_3dim; ++n)
		{
			printf("%d ",trace_back[0][m][n]);
		}
		printf("\n\n");
	}
	printf("trace back matrix with m=0:\n");
	for (l=0; l<l_1dim; ++l)
	{
		for (n=0; n<l_3dim; ++n)
		{
			printf("%d ",trace_back[l][0][n]);
		}
		printf("\n\n");
	}
	printf("trace back matrix with n=0:\n");
	for (l=0; l<l_1dim; ++l)
	{
		for (m=0; m<l_2dim; ++m)
		{
			printf("%d ",trace_back[l][m][0]);
		}
		printf("\n\n");
	} */
}

void create_3d_matrix(Mult_A *MA,s_fxn *sfxn,amino *acids,int *seq1,int l_seq1,
					 int *seq2,int l_seq2,int *seq3,int l_seq3)
{
	int i,l,m,n,l_diml,l_dimm,l_dimn,naas,alpha,beta,gamma,trace[7],max_pos;
	signed int ***dynamic_matrix;
	int ***trace_back;
	signed int gp,score[7],max;
	signed int *scoring_matrix;

	/* there are 7 different moves that can be made to one element in cube.
	   the following recursive scoring scheme goes as follows:
			T = score of element in cube
			S = score of aa to aa match (Blosum matrix)
			a,b,c = aa
			gp = gap_penalty

			T(l,m,n) = max ( T(l-1,m-1,n-1) + S(a(l),b(m)) + S(a(l),c(n)) + S(b(m),c(n))
							 T(l-1,m-1,n) + S(a(l),b(m)) + gp + gp
							 T(l-1,m,n-1) + S(a(l),c(n)) + gp + gp
							 T(l,m-1,n-1) + S(b(m),c(n)) + gp + gp
							 T(l,m,n-1) + gp + gp
							 T(l,m-1,n) + gp + gp
							 T(l-1,m,n) + gp + gp
						   )
			for trace_back matrix, (same order as for T above)
			path = element of trace_back matrix

			path = 0 => (a,b,c) move (corresponding to T(l-1,m-1,n-1) + etc above being max
			path = 1 => (a,b,-)
			path = 2 => (a,-,c)
			path = 3 => (-,b,c)
			path = 4 => (-,-,c)
			path = 5 => (-,b,-)
			path = 6 => (a,-,-)
	*/

	l_diml=MA->l_1dim;
	l_dimm=MA->l_2dim;
	l_dimn=MA->l_3dim;
	dynamic_matrix=MA->dynamic_matrix;
	trace_back=MA->trace_back;
	scoring_matrix=sfxn->scoring_matrix;
	naas=acids->naas;
	gp=sfxn->gap_penalty;

	for (l=1; l<l_diml; ++l)
	{
		for (m=1; m<l_dimm; ++m)
		{
			for (n=1; n<l_dimn; ++n)
			{
				/* get amino acids at site l,m,n of sequences */
				alpha=*(seq1+(l-1));
				beta=*(seq2+(m-1));
				gamma=*(seq3+(n-1));

				/* calculate all possible scores */
				score[0]=dynamic_matrix[l-1][m-1][n-1]+*(scoring_matrix+alpha*(naas-1)+beta)+*(scoring_matrix+alpha*(naas-1)+gamma)+*(scoring_matrix+beta*(naas-1)+gamma);
				score[1]=dynamic_matrix[l-1][m-1][n]+*(scoring_matrix+alpha*(naas-1)+beta)+gp+gp;
				score[2]=dynamic_matrix[l-1][m][n-1]+*(scoring_matrix+alpha*(naas-1)+gamma)+gp+gp;
				score[3]=dynamic_matrix[l][m-1][n-1]+*(scoring_matrix+beta*(naas-1)+gamma)+gp+gp;
				score[4]=dynamic_matrix[l][m][n-1]+gp+gp;
				score[5]=dynamic_matrix[l][m-1][n]+gp+gp;
				score[6]=dynamic_matrix[l-1][m][n]+gp+gp;
				
				/* find maximum of all possible scores */
				max=score[0];
				max_pos=0;
				for (i=0; i<7; ++i)
				{
					if (score[i]>max)
					{
						max=score[i];
						max_pos=i;
					}
					trace[i]=i;
				}

				dynamic_matrix[l][m][n]=max;
				trace_back[l][m][n]=trace[max_pos];
			}
		}
	}
}

void get_optimal_3d_score(Mult_A *MA,int l_seq1,int l_seq2,int l_seq3)
{
	MA->optimal_3d_score=MA->dynamic_matrix[l_seq1][l_seq2][l_seq3];
	MA->optimal_l=MA->l_seq1;
	MA->optimal_m=MA->l_seq2;
	MA->optimal_n=MA->l_seq3;
}

void get_mult_alignment(Mult_A *MA,amino *acids,int *seq1,int l_seq1,
					 int *seq2,int l_seq2,int *seq3,int l_seq3)
{
	int i,j,k,max_l_alignment,end,last,q,r,s,path,l_1dim,l_2dim,l_3dim;
	int l_alignment,alpha,beta,gamma,temp_1,temp_2,temp_3,n_switches;
	int ***trace_back;
	int **alignment;
	/* allocate alignment matrix: each column is an alignment triplet of the following types:
			path = 0 => (a,b,c) move (corresponding to T(l-1,m-1,n-1) + etc above being max
			path = 1 => (a,b,-)
			path = 2 => (a,-,c)
			path = 3 => (-,b,c)
			path = 4 => (-,-,c)
			path = 5 => (-,b,-)
			path = 6 => (a,-,-)
	   max size of alignment matrix = (length of query) + (length of database) */
	max_l_alignment=(l_seq1+l_seq2+l_seq3);
	MA->alignment=malloc(3*sizeof(int*));
	if (!MA->alignment)
	{
		printf("can not allocate memory for mult alignment matrix\n");
		exit(1);
	}
	for (i=0; i<3; ++i)
	{
		*(MA->alignment+i)=malloc(max_l_alignment*sizeof(int));
		if (!*(MA->alignment+i))
		{
			printf("can not allocate 2nd dimension of mult alignment matrix\n");
			exit(1);
		}
	}
	alignment=MA->alignment;
	trace_back=MA->trace_back;
	l_1dim=MA->l_1dim;
	l_2dim=MA->l_2dim;
	l_3dim=MA->l_3dim;

	q=MA->optimal_l; /* q => alpha */
	r=MA->optimal_m; /* r => beta */
	s=MA->optimal_n; /* s => gamma */
	
	/* trace back alignment path from position of optimal score */
	l_alignment=0;
	end=FALSE;
	while (end==FALSE)
	{
		if (q==0 && r==0 && s==0)
		{
			break;
		} 
		if (q==0 && r!=0 && s!=0)
		{
			/* reinitialize boundary of trace_back matrix for this particular face to path 
			designation used in 2d problem because now we are doing trace back as in pairwise
			alignments,ie to reorient the trace back because we are now on a face */
			for (j=0; j<l_2dim; ++j)
			{
				trace_back[0][j][0]=2;
			}
			for (k=0; k<l_3dim; ++k)
			{
				trace_back[0][0][k]=0;
			}
			trace_back[0][0][0]=3;

			/* reached face comprising seq2 and seq3, this face of the trace_back matrix
			   has the path convention of the 2d case (0=left,1=diag,2=top) */
			/* trace back alignment path from position of optimal score */
			last=FALSE;
			while (last==FALSE)
			{
				path=trace_back[q][r][s];
				if (path==0) /* horizontal path */
				{
					beta=s-1;
					alignment[0][l_alignment]=aachar2integer('-',acids); 
					alignment[1][l_alignment]=aachar2integer('-',acids);
					alignment[2][l_alignment]=*(seq3+beta);
					--s;
				}
				if (path==1)
				{
					alpha=r-1;
					beta=s-1;
					alignment[0][l_alignment]=aachar2integer('-',acids);
					alignment[1][l_alignment]=*(seq2+alpha);
					alignment[2][l_alignment]=*(seq3+beta);
					--r;
					--s;
				}
				if (path==2)
				{
					alpha=r-1;
					alignment[0][l_alignment]=aachar2integer('-',acids);
					alignment[1][l_alignment]=*(seq2+alpha);
					alignment[2][l_alignment]=aachar2integer('-',acids);
					--r;
				}
				if (path==3)
				{
					/* end of alignment has been reached */
					alpha=r;
					beta=s;
					alignment[0][l_alignment]=aachar2integer('-',acids);
					alignment[1][l_alignment]=*(seq2+alpha);
					alignment[2][l_alignment]=*(seq3+beta);
					last=TRUE;
					end=TRUE;
					continue;
				}
			}
			++l_alignment;
			continue;
		}
		if (q!=0 && r==0 && s!=0)
		{
			/* reinitialize boundary of trace_back matrix for this particular face to path 
			designation used in 2d problem because now we are doing trace back as in pairwise
			alignments,ie to reorient the trace back because we are now on a face */
			for (i=0; i<l_1dim; ++i)
			{
				trace_back[i][0][0]=2;
			}
			for (k=0; k<l_3dim; ++k)
			{
				trace_back[0][0][k]=0;
			}
			trace_back[0][0][0]=3;

			/* reached face comprising seq1 and seq3, this face of the trace_back matrix
			   has the path convention of the 2d case (0=left,1=diag,2=top) */
			/* trace back alignment path from position of optimal score */
			last=FALSE;
			while (last==FALSE)
			{
				path=trace_back[q][r][s];
				if (path==0) /* horizontal path */
				{
					beta=s-1;
					alignment[0][l_alignment]=aachar2integer('-',acids); 
					alignment[1][l_alignment]=aachar2integer('-',acids);
					alignment[2][l_alignment]=*(seq3+beta);
					--s;
				}
				if (path==1)
				{
					alpha=q-1;
					beta=s-1;
					alignment[0][l_alignment]=*(seq1+alpha);
					alignment[1][l_alignment]=aachar2integer('-',acids);
					alignment[2][l_alignment]=*(seq3+beta);
					--q;
					--s;
				}
				if (path==2)
				{
					alpha=q-1;
					alignment[0][l_alignment]=*(seq1+alpha);
					alignment[1][l_alignment]=aachar2integer('-',acids);
					alignment[2][l_alignment]=aachar2integer('-',acids);
					--q;
				}
				if (path==3)
				{
					/* end of alignment has been reached */
					alpha=q;
					beta=s;
					alignment[0][l_alignment]=*(seq1+alpha);
					alignment[1][l_alignment]=aachar2integer('-',acids);
					alignment[2][l_alignment]=*(seq3+beta);
					last=TRUE;
					end=TRUE;
					continue;
				}
			}
			++l_alignment;
			continue;
		}
		if (q!=0 && r!=0 && s==0)
		{
			/* reinitialize boundary of trace_back matrix for this particular face to path 
			designation used in 2d problem because now we are doing trace back as in pairwise
			alignments,ie to reorient the trace back because we are now on a face */
			for (i=0; i<l_1dim; ++i)
			{
				trace_back[i][0][0]=2;
			}
			for (j=0; j<l_2dim; ++j)
			{
				trace_back[0][j][0]=0;
			}
			trace_back[0][0][0]=3;

			/* reached face comprising seq1 and seq2, this face of the trace_back matrix
			   has the path convention of the 2d case (0=left,1=diag,2=top) */
			/* trace back alignment path from position of optimal score */
			last=FALSE;
			while (last==FALSE)
			{
				path=trace_back[q][r][s];
				if (path==0) /* horizontal path */
				{
					beta=r-1;
					alignment[0][l_alignment]=aachar2integer('-',acids); 
					alignment[1][l_alignment]=*(seq2+beta);
					alignment[2][l_alignment]=aachar2integer('-',acids);
					--r;
				}
				if (path==1)
				{
					alpha=q-1;
					beta=r-1;
					alignment[0][l_alignment]=*(seq1+alpha);
					alignment[1][l_alignment]=*(seq2+beta);
					alignment[2][l_alignment]=aachar2integer('-',acids);
					--q;
					--r;
				}
				if (path==2)
				{
					alpha=q-1;
					alignment[0][l_alignment]=*(seq1+alpha);
					alignment[1][l_alignment]=aachar2integer('-',acids);
					alignment[2][l_alignment]=aachar2integer('-',acids);
					--q;
				}
				if (path==3)
				{
					/* end of alignment has been reached */
					alpha=q;
					beta=r;
					alignment[0][l_alignment]=*(seq1+beta);
					alignment[1][l_alignment]=*(seq2+alpha);
					alignment[2][l_alignment]=aachar2integer('-',acids);
					last=TRUE;
					end=TRUE;
					continue;
				}
			}
			++l_alignment;
			continue;
		}
		if (q!=0 && r!=0 && s!=0)
		{
			/* we are in the interior of the dynamic cube so there are 7 different possible 
			paths to a certain element in the cube */
			path=trace_back[q][r][s];
			switch (path) {
			case 0:
				alpha=q-1;
				beta=r-1;
				gamma=s-1;
				alignment[0][l_alignment]=*(seq1+alpha);
				alignment[1][l_alignment]=*(seq2+beta);
				alignment[2][l_alignment]=*(seq3+gamma);
				--q;
				--r;
				--s;
				break;
			case 1:
				alpha=q-1;
				beta=r-1;
				alignment[0][l_alignment]=*(seq1+alpha);
				alignment[1][l_alignment]=*(seq2+beta);
				alignment[2][l_alignment]=aachar2integer('-',acids);
				--q;
				--r;
				break;
			case 2:
				alpha=q-1;
				gamma=s-1;
				alignment[0][l_alignment]=*(seq1+alpha);
				alignment[1][l_alignment]=aachar2integer('-',acids);
				alignment[2][l_alignment]=*(seq3+gamma);
				--q;
				--s;
				break;
			case 3:
				beta=r-1;
				gamma=s-1;
				alignment[0][l_alignment]=aachar2integer('-',acids);
				alignment[1][l_alignment]=*(seq2+beta);
				alignment[2][l_alignment]=*(seq3+gamma);
				--r;
				--s;
				break;
			case 4:
				gamma=s-1;
				alignment[0][l_alignment]=aachar2integer('-',acids);
				alignment[1][l_alignment]=aachar2integer('-',acids);
				alignment[2][l_alignment]=*(seq3+gamma);
				--s;
				break;
			case 5:
				beta=r-1;
				alignment[0][l_alignment]=aachar2integer('-',acids);
				alignment[1][l_alignment]=*(seq2+beta);
				alignment[2][l_alignment]=aachar2integer('-',acids);
				--r;
				break;
			case 6:
				alpha=q-1;
				alignment[0][l_alignment]=*(seq1+alpha);
				alignment[1][l_alignment]=aachar2integer('-',acids);
				alignment[2][l_alignment]=aachar2integer('-',acids);
				--q;
				break;
			case 7:
				/* end of alignment has been reached */
				alpha=q;
				beta=r;
				gamma=s;
				alignment[0][l_alignment]=*(seq1+alpha);
				alignment[1][l_alignment]=*(seq2+beta);
				alignment[2][l_alignment]=*(seq3+gamma);
				end=TRUE;
				break; 
			default:
				printf("path not taken in trace_back of mult alignment\n");
				exit(1);
			}
			++l_alignment;
			continue;
		}
/*		++l_alignment; */
	}
	MA->l_alignment=l_alignment;

	/* since end of alignment is beginning of alignment matrix, need to reverse alignment
	   matrix */
	if ((l_alignment%2)==0)
	{
		n_switches=(l_alignment/2);
	}
	else
	{
		n_switches=(l_alignment-1)/2;
	}
	j=l_alignment-1;
	for (i=0; i<n_switches; ++i)
	{
		temp_1=alignment[0][i];
		temp_2=alignment[1][i];
		temp_3=alignment[2][i];
		alignment[0][i]=alignment[0][j];
		alignment[1][i]=alignment[1][j];
		alignment[2][i]=alignment[2][j];
		alignment[0][j]=temp_1;
		alignment[1][j]=temp_2;
		alignment[2][j]=temp_3;
		--j;
	} 
}

void prnt_mult_alignment(Mult_A *MA,amino *acids,FILE *a_file)
{
	int i,j,l_alignment;
	int **alignment;

	l_alignment=MA->l_alignment;
	alignment=MA->alignment;

	fprintf(a_file,"multiple alignment for sequence %d (top), sequence %d (middle), and sequence %d (bottom):\n",1,3,5);
	for (i=0; i<3; ++i)
	{
		if (i==0)
		{
			fprintf(a_file,"sequence %d:\t",1);
		}
		if (i==1)
		{
			fprintf(a_file,"sequence %d:\t",3);
		}
		if (i==2)
		{
			fprintf(a_file,"sequence %d:\t",5);
		}
		for (j=0; j<l_alignment; ++j)
		{
			fprintf(a_file,"%c",integer2aachar(alignment[i][j],acids));
		}
		fprintf(a_file,"\n");
	}
	fprintf(a_file,"\n\n");
}
