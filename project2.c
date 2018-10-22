/*	CS 626: Project 2
	program to find optimal global alignment of a sequence and a structure = threading

    this program is specifically design to thread the leghemoglobin sequence into the
	myoglobin structure and vice versa.  However, this will be generalized later and
	combined with sequence-sequence alignments in a larger program.

	written by Andrew Smith
	started on 4/20/01
*/

#include "project2.h"

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
	struc=&struc_data;

	/* read command line and set options */
	/* default */
	op->local=0;				/* global alignement */
	op->n_shuffles=100;			/* number random sequences = 100 (the more the better, but
								   the more the computation cost) */
	op->contact=0;              /* default threading potential is contact potential */
	op->threading=1;
	command_line(argc,argv,op); 

	/* open project files for reading and writing */
	open_files(file); 

	/* read into memory scoring matrix (BLOSUM 50) */
	if(!read_sm(sfxn,file->f_sm,acids,struc))
	{
		printf("unable to read scoring matrix\n");
		exit(1);
	}

	/* read pdb files and write info to data structures and files */
	read_data_write_files(file,seq,struc,acids);
	
	/* calculate contact map and contact representation for leghemoglobin */
	calc_contact_map(struc,seq,seq->l_leg,seq->leg_seq,struc->leg_sc_geocntr,struc->leg_contact_map);
	printintmatrix2file(file->f_contact,struc->leg_contact_map,seq->l_leg,seq->l_leg);
	calc_contact_rep(struc,seq,seq->l_leg,seq->leg_seq,struc->leg_contact_rep,struc->leg_contact_map);
	printintvector2file(file->f_crep,struc->leg_contact_rep,seq->l_leg);
	calc_contact_types(struc,seq,acids,seq->l_leg,seq->leg_seq,struc->leg_contact_map,struc->leg_contact_types);
	printintmatrix2file(file->f_ctypes,struc->leg_contact_types,seq->l_leg,(acids->naas)-1);

	/* calculate contact map and contact representation for myoglobin */
	calc_contact_map(struc,seq,seq->l_myo,seq->myo_seq,struc->myo_sc_geocntr,struc->myo_contact_map);
	printintmatrix2file(file->f_contact,struc->myo_contact_map,seq->l_myo,seq->l_myo);
	calc_contact_rep(struc,seq,seq->l_myo,seq->myo_seq,struc->myo_contact_rep,struc->myo_contact_map);
	printintvector2file(file->f_crep,struc->myo_contact_rep,seq->l_myo);
	calc_contact_types(struc,seq,acids,seq->l_myo,seq->myo_seq,struc->myo_contact_map,struc->myo_contact_types);
	printintmatrix2file(file->f_ctypes,struc->myo_contact_types,seq->l_myo,(acids->naas)-1);

	/* do an initial reading of sequence file to get the number of sequences so as to 
	   allocate space for matrix T(i,j), which stores the optimal scores for aligning
	   sequence i into sequence j, and space for z_score matrix of the same size */
	if (op->threading)
	{
		/* close query seq file after writing and reopen for reading */
		fclose(file->f_qseq);
		file->f_qseq=fopen(file->f_qseq_name,"r");
		if (file->f_qseq==NULL)
		{
			printf("can not reopen file %s for reading\n",file->f_qseq_name);
			exit(1);
		}
		/* close database seq file after writing and reopen for reading */
		fclose(file->f_dbseq); 
		file->f_dbseq=fopen(file->f_dbseq_name,"r");
		if (file->f_dbseq==NULL)
		{
			printf("can not reopen file %s for reading\n",file->f_dbseq_name);
			exit(1);
		}
	}
	seq->num_seq=0;
	while (read_seq(seq,file->f_dbseq,2))
	{		
	}
	rewind(file->f_dbseq);
	/* allocate space */
	n_seqs=seq->num_seq;
	DP->T_i_j=malloc((n_seqs*n_seqs)*sizeof(signed int));
	if (!DP->T_i_j)
	{
		printf("can not allocate memory for T(i,j) matrix\n");
		exit(1);
	}
/*	zs->z_score=malloc((n_seqs*n_seqs)*sizeof(float));
	if (!zs->z_score)
	{
		printf("can not allocate memory for z_score matrix\n");
		exit(1);
	} */

	if (op->threading)
	{
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
			dynamic_program(q_seq,db_seq,l_qseq,l_dbseq);
			
			/* check and see if end of database file; if so, break while read sequences */
			if (seq->end_dbseqs==1) break;
/*			printf("reading next database seq:\n"); */
		}

		/* check and see if end of query file; if so, break while read sequences */
		if (seq->end_qseqs==1) break;
/*		printf("reading next query seq:\n"); */
	}
	}

	/* print T(i,j) matrix to the file "scores.txt" */
	prnt_T(DP,seq,file->f_scores,op);
	fprintf(file->f_scores,"\n\n"); 
	
	printf("\nthis is the end...my only friend, the end\n");
}
/* end of main */

int read_sm(s_fxn *sfxn,FILE *sm_file,amino *acids,str *struc) 
{
	/* put scoring fxn (Blosum50 matrix) into sfxn->scoring matrix read from sm_file*/
	char *aacharacters,amino_acid[3],aai[3],aaj[3]; 
	char row_aa;
	int i,j,n,naas,eof,aai_num,n_cont,n_cont_cat;
	int *aaintegers,aaj_num[20]; 
	float *scoring_matrix; 

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
	
	if (op->threading)
	{
		/* if using contact potential open file and read in contact energies to scoring matrix */
		if (op->contact)
		{
			strcpy(file->f_cener_name,"contact_energy.txt");
			file->f_cener=fopen(file->f_cener_name,"r");
			if (file->f_cener==NULL)
			{
				printf("can not open file %s\n",file->f_cener_name);
				return(FALSE);
			}

			eof=fscanf(file->f_cener,"%s",amino_acid);
			if(eof==EOF) 
			{
				rewind(file->f_cener);
				eof=fscanf(file->f_cener,"%d\n",amino_acid);
				if (eof==EOF)
				{
					printf("there is nothing in scoring matrix file\n");
					return(FALSE);
				}
				else
				{
					rewind(file->f_cener);
				}
			}
			else
			{
				rewind(file->f_cener);
			}
			/* allocate memory for entire scoring matrix */
			sfxn->scoring_matrix=malloc(((naas-1)*(naas-1))*sizeof(float));
			if (!sfxn->scoring_matrix)
			{
				printf("can not allocate memory for scoring matrix\n");
				return(FALSE);
			}

			/* read in actual scoring matrix */			
			for (i=0; i<naas; ++i)
			{
				if (i!=0)
				{
					fscanf(file->f_cener,"%s",aai);
					aai_num=aachar2integer(match_3letaa(aai),acids);
				}
				for (j=0; j<(naas-1); ++j)
				{
					if (i==0)
					{
						fscanf(file->f_cener,"%s",aaj);
						aaj_num[j]=aachar2integer(match_3letaa(aaj),acids);
					}
					else
					{
						n=(aai_num*(naas-1))+aaj_num[j];
/*						n=(i-1)*(naas-1)+j; */
						fscanf(file->f_cener,"%f",sfxn->scoring_matrix+n);
					}
				}
				fscanf(file->f_cener,"\n");
			}
			/* debugging */
/*			scoring_matrix=sfxn->scoring_matrix;
			printf("scoring matrix:\n");
			for (i=0; i<(naas-1); ++i)
			{
				for (j=0; j<(naas-1); ++j)
				{
					n=i*(naas-1)+j;
					printf("%d\t",*(scoring_matrix+n));
				}
				printf("\n");
			} */

			fclose(file->f_cener);
		}
		if (!op->contact) /* if profile model */
		{
			strcpy(file->f_profile_name,"profile_model.txt");
			file->f_profile=fopen(file->f_profile_name,"r");
			if (file->f_profile==NULL)
			{
				printf("can not open file %s\n",file->f_profile_name);
				return(FALSE);
			}

			eof=fscanf(file->f_profile,"%s",amino_acid);
			if(eof==EOF) 
			{
				rewind(file->f_profile);
				eof=fscanf(file->f_profile,"%d\n",amino_acid);
				if (eof==EOF)
				{
					printf("there is nothing in scoring matrix file\n");
					return(FALSE);
				}
				else
				{
					rewind(file->f_profile);
				}
			}
			else
			{
				rewind(file->f_profile);
			}
			/* allocate memory for entire scoring matrix */
			/* i do not preread file to get number of contacts categories for a site, i just
			set to 11, which is from 1-10 number of contacts plus 0 number of contact 
			category for a site */
			struc->n_cont_cat=11;
			n_cont_cat=struc->n_cont_cat;
			sfxn->scoring_matrix=malloc((n_cont_cat*(naas-1))*sizeof(float));
			if (!sfxn->scoring_matrix)
			{
				printf("can not allocate memory for scoring matrix\n");
				return(FALSE);
			}

			/* read in actual scoring matrix */			
			for (i=0; i<n_cont_cat; ++i)
			{
				if (i!=0)
				{
					fscanf(file->f_profile,"(%d)",&n_cont);
				}
				for (j=0; j<(naas-1); ++j)
				{
					if (i==0)
					{
						fscanf(file->f_profile,"%s",aaj);
						aaj_num[j]=aachar2integer(match_3letaa(aaj),acids);
						n=(i*(naas-1))+aaj_num[j];
						/* energy for zero number of contacts = 0.0 */
						*(sfxn->scoring_matrix+n)=0.0;
					}
					else
					{
						n=(n_cont*(naas-1))+aaj_num[j];
/*						n=(i-1)*(naas-1)+j; */
						fscanf(file->f_profile,"%f",sfxn->scoring_matrix+n);
					}
				}
				fscanf(file->f_profile,"\n");
			}
			/* debugging */
/*			scoring_matrix=sfxn->scoring_matrix;
			printf("scoring matrix:\n");
			for (i=0; i<n_cont_cat; ++i)
			{
				for (j=0; j<(naas-1); ++j)
				{
					n=i*(naas-1)+j;
					printf("%1.2f ",*(scoring_matrix+n));
				}
				printf("\n");
			} */

			fclose(file->f_profile);
		}
	}

	if (!op->threading)
	{
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
	}

	/* debugging */
/*	scoring_matrix=sfxn->scoring_matrix;
	printf("scoring matrix:\n");
	for (i=0; i<(naas-1); ++i)
	{
		for (j=0; j<(naas-1); ++j)
		{
			n=i*(naas-1)+j;
			printf("%d\t",*(scoring_matrix+n));
		}
		printf("\n");
	} */

	/* set constant gap penalty equal to +1 */
	sfxn->gap_penalty=1.0;
	
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
			if (!op->threading)
			{
				fclose(seq_file);
			}
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
/*	strcpy(file->f_seq_name,"sequences.txt");
	file->f_seq=fopen(file->f_seq_name,"r");
	if (file->f_seq==NULL)
	{
		printf("can not open file %s\n",file->f_seq_name);
		exit(1);
	} */
	strcpy(file->f_qseq_name,"q_sequences.txt");
	file->f_qseq=fopen(file->f_qseq_name,"w");
	if (file->f_qseq==NULL)
	{
		printf("can not open file %s\n",file->f_qseq_name);
		exit(1);
	}
	strcpy(file->f_dbseq_name,"db_sequences.txt");
	file->f_dbseq=fopen(file->f_dbseq_name,"w");
	if (file->f_dbseq==NULL)
	{
		printf("can not open file %s\n",file->f_dbseq_name);
		exit(1);
	}
	strcpy(file->f_qcoord_name,"q_coords.txt");
	file->f_qcoord=fopen(file->f_qcoord_name,"w");
	if (file->f_qcoord==NULL)
	{
		printf("can not open file %s\n",file->f_qcoord_name);
		exit(1);
	}
	strcpy(file->f_dbcoord_name,"db_coords.txt");
	file->f_dbcoord=fopen(file->f_dbcoord_name,"w");
	if (file->f_dbcoord==NULL)
	{
		printf("can not open file %s\n",file->f_dbcoord_name);
		exit(1);
	}
	strcpy(file->f_scores_name,"scores.txt");
	file->f_scores=fopen(file->f_scores_name,"w");
	if (file->f_scores==NULL)
	{
		printf("can not open file %s\n",file->f_scores_name);
		exit(1);
	}
	strcpy(file->f_contact_name,"contact_map.txt");
	file->f_contact=fopen(file->f_contact_name,"w");
	if (file->f_contact==NULL)
	{
		printf("can not open file %s\n",file->f_contact_name);
		exit(1);
	}
	strcpy(file->f_crep_name,"contact_rep.txt");
	file->f_crep=fopen(file->f_crep_name,"w");
	if (file->f_crep==NULL)
	{
		printf("can not open file %s\n",file->f_crep_name);
		exit(1);
	}
	strcpy(file->f_ctypes_name,"contact_types.txt");
	file->f_ctypes=fopen(file->f_ctypes_name,"w");
	if (file->f_ctypes==NULL)
	{
		printf("can not open file %s\n",file->f_ctypes_name);
		exit(1);
	}
}

void dynamic_program(int *q_seq,int *db_seq,int l_qseq,int l_dbseq)
{
	/* if rs_or_zs==0, then dynamic programming is for just for alignment of two sequences
	   and regular scoring (rs)
	   if rs_or_zs==1, then dynamic programming is for generation of z_score (zs)
	   get_optimal_alignment does not need to be performed and get_optimal_score needs to 
	   be modified */

	/* allcate memory for and initialize dynamic matrix */
	if (!alloc_init_matrix(DP,sfxn,op,struc,l_qseq,l_dbseq))
	{
		printf("can not allocate and initialize dynamic matrix\n");
		exit(1);
	}

	/* create dynamic matrix according to recursion formula involving the scoring matrix, 
	   which assigns a score for matching an amino acid to an amino acid, and constant gap 
	   penalty for two possible combinations, which are a gap in query or gap in database
	   sequence, either for global or local alignment */
	create_matrix(DP,sfxn,acids,op,seq,struc,q_seq,db_seq,l_qseq,l_dbseq);

	/* after creating dynamic matrix, get optimal score either for global or local 
	   alignments */
	get_optimal_score(DP,op,seq);

	/* trace back optimal path through dynamic matrix to get optimal alignment */
	get_optimal_alignment(DP,acids,q_seq,db_seq,l_qseq,l_dbseq);
	prnt_alignment(DP,acids,seq,file->f_scores);
/*	free(DP->alignment); */
}

int alloc_init_matrix(Dynamic *DP,s_fxn *sfxn,opshuns *op,str *struc,int l_qseq,int l_dbseq)
{
	int i,j,n,n_rows,n_cols;
	float *dynamic_matrix;
	float gap_penalty;
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

	DP->dynamic_matrix=malloc((n_rows*n_cols)*sizeof(float));
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
			*(dynamic_matrix+n)=(float)0.0;
			*(trace_back+n)=3;
		}
	}
	if (!op->local)
	{
		for (j=1; j<n_cols; ++j)
		{
			/* in threading, the gap_penalty is 1 per contact for both contact and profile
			   potentials.  Therefore, the far left column and top row requires not only
			   the structure of database seq but the structure of the query sequence. 
			   The gap penalty is the score of putting a "gap amino acid" into that
			   structural site and these gap scores are added as you move down the far left
			   column or to the right across the top row.
			   To avoid using structural information of the query sequence to define the
			   gap score for the far left column i could have just used the sum of the
			   contant gap penalty (1.0) as you move down the column */
			if (op->threading)
			{
				*(dynamic_matrix+j)=*(dynamic_matrix+(j-1))+calc_gap_penalty(struc,sfxn,j-1,seq->j_dbseq);
			}
			else
			{
				*(dynamic_matrix+j)=j*gap_penalty;
			}
		}
		for (i=1; i<n_rows; ++i)
		{
			if (op->threading)
			{
				*(dynamic_matrix+i*(n_cols))=*(dynamic_matrix+(i-1)*(n_cols))+calc_gap_penalty(struc,sfxn,i-1,seq->i_qseq);
/*				*(dynamic_matrix+i*(n_cols))=i*gap_penalty; */ /* or this one for no structural info */
			}
			else
			{
				*(dynamic_matrix+i*(n_cols))=i*gap_penalty;
			}
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
			printf("%1.2f ",*(dynamic_matrix+n));
		}
		printf("\n");
	} */

	return(TRUE);
}

void create_matrix(Dynamic *DP,s_fxn *sfxn,amino *acids,opshuns *op,sqs *seq,str *struc,
				   int *q_seq,int *db_seq,int l_qseq,int l_dbseq)
{
	int i,j,n,n_rows,n_cols,alpha,beta,current,diag,left,top,naas,trace;
	float *dynamic_matrix;
	float gap_penalty,current_score,diag_score,left_score,top_score;
	int *trace_back;
	float *scoring_matrix;

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
       
	   for threading, the above procedure holds for both the contact potential using the
	   frozen environment approximation and for the profile model.  However, the optimal
	   score for an element of the dynamic matrix is not the maximum of the three possible
	   options but the minimum.
			
			  the calculation of the scores using contact or profile involves a call to the
			  function calc_energy(query aa,db struc site,int to specify contact or profile):
					
					  for contact:  calc_energy(alpha(i),j,1);
					  for profile:  calc_energy(alpha(i),j,0);
			  
			  calc_energy returns a floating point number that is the energy of putting aa
			  alpha into structural site j of the database structure.

			  (note: query sequence is down the rows and database structure is across the
			         columns)

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
			if (op->threading)
			{
				/* seq->j_dbseq says whether the structure is leghemoglobin or myoglobin */
				diag_score=*(dynamic_matrix+diag)+calc_energy(seq,struc,acids,sfxn,alpha,j-1,seq->j_dbseq);
			}
			else
			{
				diag_score=*(dynamic_matrix+diag)+*(scoring_matrix+alpha*(naas-1)+beta);
			}
			current_score=diag_score;
			trace=1;

			if (op->threading)
			{
				/* score from left involves putting a gap in the sequence */
				left_score=*(dynamic_matrix+left)+calc_gap_penalty(struc,sfxn,j-1,seq->j_dbseq);
				/* score from top involves putting a gap in the structure */
				top_score=*(dynamic_matrix+top)+gap_penalty;
/*				top_score=*(dynamic_matrix+top)+calc_gap_penalty(struc,sfxn,i,seq->i_qseq); */
			}
			else
			{
				left_score=*(dynamic_matrix+left)+gap_penalty;
				top_score=*(dynamic_matrix+top)+gap_penalty;
			}

			if (left_score<diag_score && left_score<top_score)
			{
				current_score=left_score;
				trace=0;
			}
			if (top_score<diag_score && top_score<left_score)
			{
				current_score=top_score;
				trace=2;
			}
			/* the following may be unnecessary or incorrect */
			if (left_score==top_score && left_score<diag_score && top_score<diag_score)
			{
				current_score=left_score;
				trace=0;
				if (l_dbseq<l_qseq)
				{
					current_score=top_score;
					trace=2;
				}
			}

			if (op->local==1 && current_score>0)
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
/*	if (seq->i_qseq==1 && seq->j_dbseq==1)
	{
	fprintf(file->f_scores,"dynamic matrix:\n");
	for (i=0; i<n_rows; ++i)
	{
		for (j=0; j<n_cols; ++j)
		{
			n=i*n_cols+j;
			fprintf(file->f_scores,"%1.2f ",*(dynamic_matrix+n));
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

void get_optimal_score(Dynamic *DP,opshuns *op,sqs *seq)
{
	float *T_i_j,*dynamic_matrix;
	float max;
	signed int ith_pos,jth_pos;
	int i,j,n,n_rows,n_cols,i_qseq,j_dbseq,n_seq,max_position,current,rs_or_zs;

	T_i_j=DP->T_i_j;
	dynamic_matrix=DP->dynamic_matrix;
	n_rows=DP->n_rows;
	n_cols=DP->n_cols;
	i_qseq=seq->i_qseq;
	j_dbseq=seq->j_dbseq;
	n_seq=seq->num_seq;

	rs_or_zs=0;/* for now, set just to compute regular score instead of score for z-score */

	current=(i_qseq-1)*n_seq+(j_dbseq-1);
	/* for global alignment (op->lorg==0) optimal score is bottom-right corner of dynamic
	   matrix */
	if (!op->local)
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
	if (op->local)
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

	if (op->threading)
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
	{
		fprintf(a_file,"alignment for query sequence %d and database structure %d:\n",seq->i_qseq,seq->j_dbseq);
	}

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

void prnt_T(Dynamic *DP,sqs *seq,FILE *a_scores,opshuns *op)
{
	float *T_i_j;
	int i,j,n,n_seq;

	T_i_j=DP->T_i_j;
	n_seq=seq->num_seq;

	fprintf(a_scores,"T(i,j):\n");
	if (op->local==0)
	{
		fprintf(a_scores,"global alignment\n");
		if (op->threading)
		{
			fprintf(a_scores,"1 = leghemoglobin and 2 = myoglobin\n\n");
		}
	}
	if (op->local==1)
	{
		fprintf(a_scores,"local alignment\n\n");
	}
	fprintf(a_scores,"database structure:\t");
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
			fprintf(a_scores,"%0f\t",*(T_i_j+n));
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
		if (strcmp(argv[i],"-profile")==0)
		{
			op->contact=0; 
		}
		if (strcmp(argv[i],"-s2s")==0)
		{
			op->threading=0;
		}
	}
}

void read_data_write_files(p_files *file,sqs *seq,str *struc,amino *acids)
{
	FILE *qseq_file;
	FILE *dbseq_file;
	FILE *qcoord_file;
	FILE *dbcoord_file;

	int i,j,l_leg,l_myo;
	char ch;
	float coord;

	qseq_file=file->f_qseq;
	qcoord_file=file->f_qcoord;
	dbseq_file=file->f_dbseq;
	dbcoord_file=file->f_dbcoord;
	/* this is specific for leghemoglobin and myoglobin, but will be generalized later */

	/* read pdb file of leghemoglobin and write its sequence and coordinates of the 
	geometric centers of its side chains to a query sequence file, database sequence file,
	a query coordinate file (for struc-struc alignments later), and database coordinate
	file */

	/* open leghemoglobin PDB file */
	strcpy(file->f_pdb1_name,"1GDJ.txt");
	file->f_pdb1=fopen(file->f_pdb1_name,"r");
	if (file->f_pdb1==NULL)
	{
		printf("can not open file %s\n",file->f_pdb1_name);
		exit(1);
	}
	/* open myoglobin PDB file */
	strcpy(file->f_pdb2_name,"110M.txt");
	file->f_pdb2=fopen(file->f_pdb2_name,"r");
	if (file->f_pdb2==NULL)
	{
		printf("can not open file %s\n",file->f_pdb2_name);
		exit(1);
	}
	/* read leghemoglobin file */
	fprintf(qseq_file,"1.\n");
	fprintf(dbseq_file,"1.\n");
	fprintf(qcoord_file,"leghemoglobin:\n");
	fprintf(dbcoord_file,"leghemoglobin:\n");
	read_PDB_calc_gc(acids,file->f_pdb1,file->f_pdb1_name,seq,struc);

	/* write leghemoglobin sequence and geometric centers to a query and database sequence 
	   and coordinate files */
	l_leg=seq->length_seq;
	seq->l_leg=l_leg;
	/* allocate memory for leghemoglobin sequence */
	seq->leg_seq=malloc(l_leg*(sizeof(int)));
	if (!seq->leg_seq)
	{
		printf("can not allocate memory for seq->leg_seq\n");
		exit(1);
	}
	/* allocate memory for leghemoglobin side chain geometric centers */
	struc->leg_sc_geocntr=malloc((l_leg*3)*sizeof(float));
	if (!struc->leg_sc_geocntr)
	{
		printf("can not allocate memory for struc->leg_sc_geocntr\n");
		exit(1);
	}
	for (i=0; i<l_leg; ++i)
	{
		ch=*(seq->sequence+i);
		fprintf(qseq_file,"%c",ch);
		fprintf(dbseq_file,"%c",ch);
		*(seq->leg_seq+i)=aachar2integer(ch,acids);
		for (j=0; j<3; ++j)
		{
			coord=*(struc->sc_geocenters+(i*3)+j);
			fprintf(qcoord_file,"%f  ",coord);
			fprintf(dbcoord_file,"%f  ",coord);
			*(struc->leg_sc_geocntr+(i*3)+j)=coord;
		}
		fprintf(qcoord_file,"\n");
		fprintf(dbcoord_file,"\n");
	}
	fprintf(qseq_file,"\n");
	fprintf(dbseq_file,"\n");
	
	/* read pdb file of myoglobin and write its sequence and coordinates of the 
	geometric centers of its side chains to a query sequence file, database sequence file,
	a query coordinate file (for struc-struc alignments later), and database coordinate
	file */
	/* read myoglobin file */
	fprintf(qseq_file,"2.\n");
	fprintf(dbseq_file,"2.\n");
	fprintf(qcoord_file,"myoglobin:\n");
	fprintf(dbcoord_file,"myoglobin:\n");
	read_PDB_calc_gc(acids,file->f_pdb2,file->f_pdb2_name,seq,struc);

	/* write leghemoglobin sequence and geometric centers to a query and database sequence 
	   and coordinate files */
	l_myo=seq->length_seq;
	seq->l_myo=l_myo;
	/* allocate memory for myoglobin sequence */
	seq->myo_seq=malloc(l_myo*(sizeof(int)));
	if (!seq->myo_seq)
	{
		printf("can not allocate memory for seq->myo_seq\n");
		exit(1);
	}
	/* allocate memory for myoglobin side chain geometric centers */
	struc->myo_sc_geocntr=malloc((l_myo*3)*sizeof(float));
	if (!struc->myo_sc_geocntr)
	{
		printf("can not allocate memory for struc->myo_sc_geocntr\n");
		exit(1);
	}
	for (i=0; i<l_myo; ++i)
	{
		ch=*(seq->sequence+i);
		fprintf(qseq_file,"%c",ch);
		fprintf(dbseq_file,"%c",ch);
		*(seq->myo_seq+i)=aachar2integer(ch,acids);
		for (j=0; j<3; ++j)
		{
			coord=*(struc->sc_geocenters+(i*3)+j);
			fprintf(qcoord_file,"%f  ",coord);
			fprintf(dbcoord_file,"%f  ",coord);
			*(struc->myo_sc_geocntr+(i*3)+j)=coord;
		}
		fprintf(qcoord_file,"\n");
		fprintf(dbcoord_file,"\n");
	}
	fprintf(qseq_file,"\n");
	fprintf(dbseq_file,"\n");

	/* allocate other structure data such as leghemoglobin and myoglobin contact map and 
	contact representations */
	allocate_struc_data(struc,acids,l_leg,l_myo);

/*	fprintf(seq_file,"1.\n");
	fprintf(coord_file,"1GDJ:\n"); 
	read_PDB_calc_gc(acids,file->f_pdb1,file->f_pdb1_name,seq,struc);
/*	fclose(file->f_pdb1);

/*	fprintf(seq_file,"2.\n");
	fprintf(coord_file,"110M:\n"); 
	read_PDB_calc_gc(acids,file->f_pdb2,file->f_pdb2_name,seq,struc); */

	free(seq->sequence);
	free(struc->sc_geocenters);
	fclose(file->f_pdb1);
	fclose(file->f_pdb2); 
}

void read_PDB_calc_gc(amino *acids,FILE *pdb_file,char *file_name,sqs *seq,str *struc)
{
	char line[81],*token,aa,*atom_num,*atom_type,*aa_num,*aa_type;
	char *x_coord,*y_coord,*z_coord,prev_aa[3],prev_aa_num[1];
	int i,j,k,r,s,once,num_sca,one_atom;
	float *side_chain_atoms,avg_x_coord,avg_y_coord,avg_z_coord;

	/* this function does the following:
		1. reads the sequence and atomic coordinates in the PDB file of a protein
		2. writes the sequence to a specified sequence file 
	    3. calculates the geometric center of the side chains storing the coordinates of the
		   geometric centers in a matrix (1dim = aa, 2dim = x_coord, y_coord, z_coord)
		4. writes the above matrix to a file
	*/

	/* because TRP has largest number of atoms (14, including backbone) i will set the size 
	of the matrix to store the coordinates of side chain atoms as 14x3.  NSCA under 
	definitions equals 14 for this max number of side chain atoms */
	side_chain_atoms=malloc((NSCA*3)*sizeof(float));
	/* initialize side_chain_atoms matrix because will take an average over it later */
	for (r=0; r<NSCA; ++r)
	{
		for (s=0; s<3; ++s)
		{
			*(side_chain_atoms+(r*3)+s)=0.0;
		}
	}

	i=0;
	once=0;
	num_sca=0;
	fgets(line,82,pdb_file);
/*	printf("%s",line); */
	while (line)
	{
/*		token=strtok(line," ");
		if (strcmp(token,"END")==0) break; */
		if (strncmp(line,"END",3)==0) break;
/*		printf("%s\t",token); */
/*		if (strcmp(token,"SEQRES")==0) */
		if (strncmp(line,"SEQRES",6)==0)
		{
			token=strtok(line," ");
			while (token)
			{
				token=strtok(NULL," ");
				if (token==NULL) break;
				if ((strcmp(token,"153")==0 || strcmp(token,"154")==0) && once==0)
				{
					/* quick determination of length, generalize later */
					seq->length_seq=atoi(token);
					once=1;
					/* allocate local storage of sequence */
					seq->sequence=malloc(seq->length_seq*(sizeof(char)));
					if (!seq->sequence)
					{
						printf("can not allocate memory for seq->sequence\n");
						exit(1);
					}
					/* allocate local storage of coordinate matrix of geometric centers of
					side chains */
					struc->sc_geocenters=malloc((seq->length_seq*3)*sizeof(float));
					if (!struc->sc_geocenters)
					{
						printf("can not allocate memory for struc->sc_geocenters\n");
						exit(1);
					}
					/* initialize matrix to zero */
					for (j=0; j<seq->length_seq; ++j)
					{
						for (k=0; k<3; ++k)
						{
							*(struc->sc_geocenters+(i*3)+j)=0.0;
						}
					}
/*					printf("from pdb, l_seq = %d",l_seq); */
				}
				aa=match_3letaa(token);
				if (isalpha(aa))
				{
					*(seq->sequence+i)=aa;
					++i;
/*					printf("%c  ",aa); */
/*					aa_num=aachar2integer(aa,acids); */
				}
			}
		}
/*		token=strtok(line," "); */
/*		if (strcmp(token,"ATOM")==0) */
		if (strncmp(line,"ATOM",4)==0 || strncmp(line,"TER",3)==0)
		{
			one_atom=0;
			token=strtok(line," ");
/*			fscanf(pdb_file,"\n"); 
			fscanf(pdb_file,"%s%d%s%s%d%f%f%f",string,&a,b,c,&d,&e,&f,&g); */
			while (token)
			{
/*				token=strtok(NULL," "); */
				if (token==NULL) break; 
				if (one_atom==0)
				{	
					atom_num=strtok(NULL," ");
					atom_type=strtok(NULL," ");
					aa_type=strtok(NULL," ");
					if (!isalpha(match_3letaa(aa_type)) && strcmp(token,"TER")!=0)
					{
						one_atom=1;
						continue;
					}
					aa_num=strtok(NULL," ");
					if ((strcmp(atom_type,"N")==0 || strcmp(token,"TER")==0) && strcmp(atom_num,"1")!=0 && strncmp(prev_aa,"GLY",3)!=0)
					{	
						/* calculate geometric center of side chain atoms of previous aa if 
						it was not glycine */
						/* the geometric center of the side chain is the average of the 
						coordinates of all the atoms of the side chain */
/*						if (atoi(aa_num)!=(atoi(prev_aa_num)+1))
						{
							printf("aa_num != prev_aa_num+1");
						} */
						avg_x_coord=0;
						avg_y_coord=0;
						avg_z_coord=0;
						for (j=0; j<num_sca; ++j)
						{
							avg_x_coord+=*(side_chain_atoms+(j*3));
							avg_y_coord+=*(side_chain_atoms+(j*3)+1);
							avg_z_coord+=*(side_chain_atoms+(j*3)+2);
						}
						avg_x_coord=(avg_x_coord/num_sca);
						avg_y_coord=(avg_y_coord/num_sca);
						avg_z_coord=(avg_z_coord/num_sca);
						if (strcmp(file_name,"1GDJ.txt")==0)
						{
							*(struc->sc_geocenters+((atoi(prev_aa_num)-1)*3))=avg_x_coord;
							*(struc->sc_geocenters+((atoi(prev_aa_num)-1)*3)+1)=avg_y_coord;
							*(struc->sc_geocenters+((atoi(prev_aa_num)-1)*3)+2)=avg_z_coord;
						}
						if (strcmp(file_name,"110M.txt")==0)
						{
							*(struc->sc_geocenters+((atoi(prev_aa_num))*3))=avg_x_coord;
							*(struc->sc_geocenters+((atoi(prev_aa_num))*3)+1)=avg_y_coord;
							*(struc->sc_geocenters+((atoi(prev_aa_num))*3)+2)=avg_z_coord;
						}
						/* reinitialize number of side chain atoms (num_sca) to zero for 
						next residue */
						num_sca=0;
						/* reinitialize side_chain_atoms matrix */
						for (r=0; r<NSCA; ++r)
						{
							for (s=0; s<3; ++s)
							{
								*(side_chain_atoms+(r*3)+s)=0;
							}
						}
					}
					strcpy(prev_aa,aa_type);
					strcpy(prev_aa_num,aa_num);
					x_coord=strtok(NULL," ");
					y_coord=strtok(NULL," ");				
					z_coord=strtok(NULL," ");
					one_atom=1;
					/* the geometric center of glycine is the coordinates of the alpha carbon */
					if (strncmp(aa_type,"GLY",3)==0 && strncmp(atom_type,"CA",2)==0)
					{
						if (strcmp(file_name,"1GDJ.txt")==0)
						{
							*(struc->sc_geocenters+((atoi(aa_num)-1)*3))=(float)atof(x_coord);
							*(struc->sc_geocenters+((atoi(aa_num)-1)*3)+1)=(float)atof(y_coord);
							*(struc->sc_geocenters+((atoi(aa_num)-1)*3)+2)=(float)atof(z_coord);
						}
						if (strcmp(file_name,"110M.txt")==0)
						{
							*(struc->sc_geocenters+((atoi(aa_num))*3))=(float)atof(x_coord);
							*(struc->sc_geocenters+((atoi(aa_num))*3)+1)=(float)atof(y_coord);
							*(struc->sc_geocenters+((atoi(aa_num))*3)+2)=(float)atof(z_coord);
						}
						one_atom=1;
					}
					/* the following is the end condition "TER" line 
					   look for length of sequence.  condition only good for specified
					   pdb files, generalize later */
					if (/* atoi(aa_type)!=l_seq || */atoi(aa_type)!=153) 
					{
/*						if (isalpha(match_3letaa(aa_type)))
						{
							/* atoms N,C,O, and CA are backbone atoms not side chain atoms */
							if (strcmp(atom_type,"N")!=0)
							{
								if (strcmp(atom_type,"O")!=0)
								{
									if (strcmp(atom_type,"CA")!=0) 
									{
										if (strcmp(atom_type,"C")!=0)
										{
											*(side_chain_atoms+(num_sca*3))=(float)atof(x_coord);
											*(side_chain_atoms+(num_sca*3)+1)=(float)atof(y_coord);
											*(side_chain_atoms+(num_sca*3)+2)=(float)atof(z_coord);
											one_atom=1;
											++num_sca;
										}
									}
								}
							}
/*						} */
					}
				}
				token=strtok(NULL," "); 
			}
		}
		fgets(line,82,pdb_file);
	}
	rewind(pdb_file);
	free(side_chain_atoms);
}

char match_3letaa(char *token)
{
	if(strncmp("GLY",token,3)==0)return('G');
	if(strncmp("ALA",token,3)==0)return('A');
    if(strncmp("ARG",token,3)==0)return('R');
	if(strncmp("ASN",token,3)==0)return('N');
	if(strncmp("ASP",token,3)==0)return('D');
	if(strncmp("CYS",token,3)==0)return('C');
	if(strncmp("GLN",token,3)==0)return('Q');
	if(strncmp("GLU",token,3)==0)return('E');
	if(strncmp("HIS",token,3)==0)return('H');
	if(strncmp("ILE",token,3)==0)return('I');
	if(strncmp("LEU",token,3)==0)return('L');
	if(strncmp("LYS",token,3)==0)return('K');
	if(strncmp("MET",token,3)==0)return('M');
	if(strncmp("PHE",token,3)==0)return('F');
	if(strncmp("PRO",token,3)==0)return('P');
	if(strncmp("SER",token,3)==0)return('S');
	if(strncmp("THR",token,3)==0)return('T');
	if(strncmp("TRP",token,3)==0)return('W');
	if(strncmp("TYR",token,3)==0)return('Y');
	if(strncmp("VAL",token,3)==0)return('V');

	if(strncmp("gly",token,3)==0)return('G');
	if(strncmp("ala",token,3)==0)return('A');
    if(strncmp("arg",token,3)==0)return('R');
	if(strncmp("asn",token,3)==0)return('N');
	if(strncmp("asp",token,3)==0)return('D');
	if(strncmp("cys",token,3)==0)return('C');
	if(strncmp("gln",token,3)==0)return('Q');
	if(strncmp("glu",token,3)==0)return('E');
	if(strncmp("his",token,3)==0)return('H');
	if(strncmp("ile",token,3)==0)return('I');
	if(strncmp("leu",token,3)==0)return('L');
	if(strncmp("lys",token,3)==0)return('K');
	if(strncmp("met",token,3)==0)return('M');
	if(strncmp("phe",token,3)==0)return('F');
	if(strncmp("pro",token,3)==0)return('P');
	if(strncmp("ser",token,3)==0)return('S');
	if(strncmp("thr",token,3)==0)return('T');
	if(strncmp("trp",token,3)==0)return('W');
	if(strncmp("tyr",token,3)==0)return('Y');
	if(strncmp("val",token,3)==0)return('V');
	return (0);
}

void calc_contact_map(str *struc,sqs *seq,int l_stru,int *struc_seq,float *geo_centers,int *contact_map)
{
	int i,j;
	float dx,dy,dz,dx2,dy2,dz2,*distance_map,value;
	double dij;

	/* this function is to calculate the contact map for both the leghemoglobin structure
	   and myoglobin structure from the distances between the geometric centers of their
	   side chains */
	/* the contact map is as follows:  with r_cutoff = 6.4 angstroms
			
			  c(i,j) = contact between amino acid i and amino acid j
			  d(i,j) = distance between aa i and aa j

			  if d(i,j) <= 6.4, then c(i,j) = 1
			  if d(i,j) > 6.4, then c(i,j) = 0

			  where d(i,j) = sqrt((xi-xj)^2+(yi-yj)^2+(zi-zj)^2)
				
				  with xi,yi,zi = cartesian coordinates of geometric center of sidechain
				                  of amino acid i
				   and xj,yj,zj = coordinates of aa j
	*/

	/* allocate memory for distance map */
	distance_map=malloc((l_stru*l_stru)*sizeof(float));
	if (!distance_map)
	{
		printf("can not allocate memory for distance_map\n");
		exit(1);
	}
	/* allocate memory for contact map */
/*	contact_map=malloc((l_stru*l_stru)*sizeof(int));
	if (!contact_map)
	{
		printf("can not allocate memory for distance_map\n");
		exit(1);
	} */

	/* calculate distances and store in distance_map */
	for (i=0; i<l_stru; ++i)
	{
		for (j=0; j<l_stru; ++j)
		{
			dx=*(geo_centers+(i*3))-*(geo_centers+(j*3));
			dx2=dx*dx;
			dy=*(geo_centers+(i*3)+1)-*(geo_centers+(j*3)+1);
			dy2=dy*dy;
			dz=*(geo_centers+(i*3)+2)-*(geo_centers+(j*3)+2);
			dz2=dz*dz;
			dij=sqrt((double)(dx2+dy2+dz2));
			*(distance_map+(i*l_stru)+j)=(float)dij;
		}
	}

	/* loop through contact map converting values <= 6.4 to 1 and values > 6.4 to 0 */
	for (i=0; i<l_stru; ++i)
	{
		for (j=0; j<l_stru; ++j)
		{
			value=*(distance_map+(i*l_stru)+j);
			if (value <= 6.4)
			{
				*(contact_map+(i*l_stru)+j)=1;
			}
			else
			{
				*(contact_map+(i*l_stru)+j)=0;
			}
		}
	}
}

void calc_contact_rep(str *struc,sqs *seq,int l_stru,int *struc_seq,int *contact_rep,int *contact_map)
{
	int i,j,sum_contacts;

	/* this function calculates the contact representation of a structure, i.e. instead of
	   each structural site containing an amino acid it contains the number of contacts to
	   that structural site --> profile of structure */

	/* allocate memory for contact representation */
/*	contact_rep=malloc(l_stru*sizeof(int));
	if (!contact_rep)
	{
		printf("can not allocate memory for contact_rep\n");
		exit(1);
	} */
	
	for (i=0; i<l_stru; ++i)
	{
		sum_contacts=0;
		for (j=0; j<l_stru; ++j)
		{
			sum_contacts+=*(contact_map+(i*l_stru)+j);
		}
		*(contact_rep+i)=sum_contacts;
	}

}

void calc_contact_types(str *struc,sqs *seq,amino *acids,int l_stru,int *struc_seq,int *contact_map,int *contact_type)
{
	int i,j,naas,contact;

	/* this function calculates the number of contacts of a specific type to a structural 
	   site (such as i-ala,i-glu, etc. where i is a structural site) to be used in the 
	   contact potential withthe frozen environment approximation.  this is done by putting 
	   a query aa at site i,then we have contact of a certain type (ala-ala,ala-val,arg-pro,
	   etc) where the contact potential (20x20 matrix) gives a score for.
	   these contact types are stored in the matrix contact_type */

	naas=acids->naas;

	/* allocate memory for contact_type matrix */
/*	contact_type=malloc((l_stru*naas)*sizeof(int)); */

	/* initialize contact_type matrix to zero */
	for (i=0; i<l_stru; ++i)
	{
		for (j=0; j<naas; ++j)
		{
			*(contact_type+(i*naas)+j)=0;
		}
	}

	for (i=0; i<l_stru; ++i)
	{
		for (j=0; j<l_stru; ++j)
		{
			contact=*(contact_map+(i*l_stru)+j);
			if (contact==1)
			{
				*(contact_type+(i*naas)+*(struc_seq+j))+=1;
			}
		}
	}
}

void allocate_struc_data(str *struc,amino *acids,int l_leg_stru,int l_myo_stru)
{
	int naas;

	naas=acids->naas;

	/* allocate memory for leghemoglobin contact map and contact representation */
	struc->leg_contact_map=malloc((l_leg_stru*l_leg_stru)*sizeof(int));
	if (!struc->leg_contact_map)
	{
		printf("can not allocate memory for struc->leg_contact_map\n");
		exit(1);
	}
	struc->leg_contact_rep=malloc(l_leg_stru*sizeof(int));
	if (!struc->leg_contact_rep)
	{
		printf("can not allocate memory for struc->leg_contact_rep\n");
		exit(1);
	}
	struc->leg_contact_types=malloc((l_leg_stru*(naas-1))*sizeof(int));
	if (!struc->leg_contact_types)
	{
		printf("can not allocate memory for struc->leg_contact_types\n");
		exit(1);
	}
	/* allocate memory for myoglobin contact map and contact representation */
	struc->myo_contact_map=malloc((l_myo_stru*l_myo_stru)*sizeof(int));
	if (!struc->myo_contact_map)
	{
		printf("can not allocate memory for struc->myo_contact_map\n");
		exit(1);
	}
	struc->myo_contact_rep=malloc(l_myo_stru*sizeof(int));
	if (!struc->myo_contact_rep)
	{
		printf("can not allocate memory for struc->myo_contact_rep\n");
		exit(1);
	}
	struc->myo_contact_types=malloc((l_myo_stru*(naas-1))*sizeof(int));
	if (!struc->myo_contact_types)
	{
		printf("can not allocate memory for struc->myo_contact_types\n");
		exit(1);
	}
}

void printintmatrix2file(FILE *some_file,int *matrix,int first_dim,int second_dim)
{
	int i,j,element;

	fprintf(some_file,"\n");
	for (i=0; i<first_dim; ++i)
	{
		for (j=0; j<second_dim; ++j)
		{
			element=*(matrix+(i*second_dim)+j);
			fprintf(some_file,"%d",element);
		}
		fprintf(some_file,"\n");
	}
	fprintf(some_file,"\n");
}

void printintvector2file(FILE *some_file,int *vector,int dim)
{
	int i,element;

	fprintf(some_file,"\n");
	for (i=0; i<dim; ++i)
	{
		element=*(vector+i);
		fprintf(some_file,"%d ",element);
	}
	fprintf(some_file,"\n");
}

void close_files(p_files *file)
{
	fclose(file->f_dbcoord);
	fclose(file->f_dbseq);
	fclose(file->f_qcoord);
	fclose(file->f_qseq);
	fclose(file->f_scores);
	fclose(file->f_sm);
	fclose(file->f_contact);
	fclose(file->f_crep);
}

float calc_energy(sqs *seq,str *struc,amino *acids,s_fxn *sfxn,int alpha,int struc_site,int leg_or_myo)
{
	int i,naas,contact_aa,n_conts;
	int *contact_types,*contact_rep;
	float energy;
	float *scoring_matrix;

	naas=acids->naas;
	scoring_matrix=sfxn->scoring_matrix;
	if (leg_or_myo==1) /* structure is leghemoglobin */
	{
		contact_types=struc->leg_contact_types;
		contact_rep=struc->leg_contact_rep;
	}
	if (leg_or_myo==2) /* structure is myoglobin */
	{
		contact_types=struc->myo_contact_types;
		contact_rep=struc->myo_contact_rep;
	}

	/* this function calculates the energy of putting an query amino acid into the 
	structural site of a database structure either by the contact potential using the 
	frozen environment approximation or by the profile model 

		calculation by contact potential using frozen environment approximation:

			energy(a,i) = sum of over contacts of specific type (V-A,K-D,F-Q,etc.)
						= sum over j[e(ai,bj)] 
	
			where  j = structural sites of database structure in contact with structural 
			           site i
				  ai = query aa at structural site i
				  bj = database aa at structural site j (from native sequence of structure)
	        e(ai,bj) = contact energy between ai and bj (given) 
	
		calculation of energy using profile model:
	
	        energy(a,i) = some value saying how well an aa fits into a site with a certain
	                      number of contacts 
	                    = e(a,ni)
	
	        where  a = query amino acid
	              ni = number of contacts to site i
	         e(a,ni) = a number of putting amino acid a into site i with ni
	*/
	
	if (op->contact) /* contact potential */
	{
		energy=0.0;
		for (i=0; i<naas-1; ++i)
		{
			contact_aa=*(contact_types+struc_site*(naas-1)+i); /* this is number of contacts of type bj */
			energy+=*(scoring_matrix+alpha*(naas-1)+contact_aa);
		}
	}
	if (!op->contact) /* profile model */
	{
		n_conts=*(contact_rep+struc_site);
		energy=*(scoring_matrix+n_conts*(naas-1)+alpha);
	}

	return(energy);
}

float calc_gap_penalty(str *struc,s_fxn *sfxn,int struc_site,int leg_or_myo)
{
	float gap_pen;

	/* this function calculate the gap score at a structural site according to the following
	   formula for both contact and profile potential:

				gap_score(i) = gap_penalty * n(i)

				where  gap_score(i) = the gap score of site i
						gap_penalty = 1.0 (in this case)
							   n(i) = number of contacts of site i
	*/

	if (leg_or_myo==1) /* structure is leghemoglobin */
	{
		gap_pen=(sfxn->gap_penalty)*(float)(*(struc->leg_contact_rep+struc_site));
	}
	if (leg_or_myo==2) /* structure is myoglobin */
	{
		gap_pen=(sfxn->gap_penalty)*(float)(*(struc->myo_contact_rep+struc_site));
	}

	return(gap_pen);
}
