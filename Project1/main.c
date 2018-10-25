/*	CS 626: Project 1
	program to find optimal global and local alignments of pairs of sequences with a single
	fixed gap penalty (gp = 10)

	written by Andrew Smith
	started on 3/6/01
*/

#include "project1.h"

/* main */

main(int argc, char *argv[])
{
	seq=&seq_data;
	sfxn=&sfxn_data;
	align=&align_data;
	zscore=&zs_data;
	file=&file_data;

	/* read into memory scoring matrix (BLOSUM 50) and test sequences */
	read_sm(seq,sfxn,file); 

	/* loop over sequences twice (once for using sequences as queries, once for using sequences
	as databases) to get total optimal scores for all pairs of sequences T(i,j) 

	/* pass sequences and scoring matrix to dynamic programming function (for both global and
	local alignments) and calculate T(i,j)

	/* shuffle query sequence and calculate z-score 

	/* find optimal global alignment for three sequences using the same scoring scheme as for
	pair alignment 
	*/
}

void read_sm(sqs *seq, s_fxn *sfxn, p_files *file)
{
/*	char *sequence,*f_seq_name;
	int *l_seq; */
	char *f_sm_name,*aatypes;
	char naas,eof;
	int i,j,n;
	signed int *scoring_matrix;
	signed int gap_penalty;
/*	FILE *f_seq; */
	FILE *f_sm;

	/* variable initialization */
	naas=0;

/*	sequence=seq->q_seq; 
	f_seq=file->f_seq;
	f_seq_name=file->f_seq_name; */
	f_sm=file->f_sm;
	f_sm_name=file->f_sm_name;
	scoring_matrix=sfxn->scoring_matrix;
	gap_penalty=sfxn_data.gap_penalty;

/*	strcpy(f_seq_name,"sequences.txt");
	f_seq=fopen(f_seq_name,"r"); */

	strcpy(f_sm_name,"Blosum50.txt");
	f_sm=fopen(f_sm_name,"r");

	eof=fscanf(f_sm,"%d\n",naas);
	if(eof==EOF)
	{
		rewind(f_sm);
		fscanf(f_sm,"%d\n",naas);
	}

	aatypes=malloc(naas*sizeof(char));

	for (i=0; i<naas; ++i)
	{
		fscanf(f_sm,"%c",aatypes[i]);
	}
	fscanf(f_sm,"\n");
	
	/* debugging */
	for (i=0; i<naas; ++i)
	{
		printf("%c",aatypes[i]);
	}

	/* read in actual scoring matric */
	for (i=0; i<naas; ++i)
	{
		for (j=0; j<naas; ++j)
		{
			n=i*naas+j;
			fscanf(f_sm,"%d",*(scoring_matrix+n));
		}
	}

	/* debugging */
	for (i=0; i<naas; ++i)
	{
		for (j=0; j<naas; ++j)
		{
			n=i*naas+j;
			printf("%d",*(scoring_matrix+n));
		}
	}

	gap_penalty=-10;
	
	/* debugging */
	printf("gap_penatly = %d",gap_penalty);
}