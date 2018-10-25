/*	CS 626: Project 1
	header file for project1.c

	written by Andrew Smith
	3/6/01
*/

/* use some standard libraries */
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fcntl.h"
#include "Time.h"
#include "ctype.h"

/* some definitions */

#define TRUE 1
#define FALSE 0
#define MAXFILENAME 20

/* declare project structures */
struct options {
	int local;					/* option for performing local alignment (default global)*/
	int n_shuffles;				/* number of shuffled sequences for calculating z-score */
	int mult_align;				/* option for performing multiple alignment on sequence 6,7,9*/
} op_data, *op;
typedef struct options opshuns;

struct amino_acids {
	int naas;			/* number of amino acids (includes the gap character '-') */
	char *aacharacters;	/* array to store the actual amino acid characters (1 letter symbol)*/
	int *aaintegers;	/* array to store the integer corresponding to amino acid character */
} aas_data, *acids;
typedef struct amino_acids amino;

struct sequences {
	char *sequence;		/* array for sequence in amino acid character representation */
	int i_seq;			/* the ith sequence read from sequence file */
	int i_qseq;			/* the ith query sequence read from sequence file */
	int j_dbseq;		/* the jth database sequence read from sequence file */
	int length_seq;		/* stores length of sequence read from sequence file */
	int num_seq;		/* stores the number of sequences read */
	int end_seqs;		/* binary controller, if 1, reached end of sequence file */
	int end_qseqs;		/* same as above but for reading query sequences */
	int end_dbseqs;		/* same as above but for reading database sequences */
} seq_data,*seq; 
typedef struct sequences sqs;

struct scoring_fxn {
	signed int *scoring_matrix;		/* stores scoring function (Blosum50 matrix) */
	signed int gap_penalty;			/* stores constant gap penalty (= -10) */
} sfxn_data,*sfxn;					
typedef struct scoring_fxn s_fxn;

struct Dynam_Program {				/* for 2d alignments */
	signed int *dynamic_matrix;		/* stores dynamic matrix or table */
	signed int *T_i_j;				/* stores matrix of optimal scores for all sequences 
										aligned */
	signed int optimal_score;		/* stores optimal score for one alignment */
	signed int opt_rand_score;		/* stores optimal score for alignment with random sequence*/
	int optimal_i;					/* ith position (query, row) of optimal score */
	int optimal_j;					/* jth position (database, column) of optimal score */
	int n_rows;						/* the number of rows of the dynamic matrix */
	int n_cols;						/* the number of columns of the dynamic matrix */
	int l_alignment;				/* length of the alignment (max = length of query + length 
										of database */
	int **alignment;				/* stores actual alignment as sequence of matched pairs 
										(aa2aa, aa2gap) */
	int *trace_back;				/* stores for each element of dynamic matrix from which its 
										value came from as a integer */
} DP_data,*DP;
typedef struct Dynam_Program Dynamic;

struct Z_score {
	float *z_score;				/* stores z_score for an alignment */
	int **random_seqs;			/* stores random seqs generated for z-score as a matrix 
									(size=op->n_shuffles*length of sequence)*/
	signed int *rand_scores;	/* stores optimal scores for random sequences 
									(size=op->n_shuffles)*/
	float avg_rand_scores;		/* stores average of the optimal scores of random sequences*/
	float avg_rand_scores_sqrd;	/* stores average of the squares of the optimal scores of 
									random sequences*/
	float std_dev;				/* standard deviation, denominator for z-score */
} zs_data,*zs;
typedef struct Z_score Zscore;

struct Multiple_Alignment {
	signed int ***dynamic_matrix;	/* dynamic matrix (or cube) for 3 sequence multiple 
										alignment*/
	signed int optimal_3d_score;	/* optimal score for multiple alignment */
	int ***trace_back;				/* trace matrix to get the alignment, similar scheme as for
										2d (integer corresponding to path)*/
	int l_1dim;						/* length of dimension for sequence 1 for dynamic matrix*/
	int l_2dim;						/* length of dimension for sequence 2 for dynamic matrix*/
	int l_3dim;						/* length of dimension for sequence 3 for dynamic matrix*/
	int l_seq1;						/* length of sequence 1 */
	int l_seq2;						/* length of sequence 2 */
	int l_seq3;						/* length of sequence 3 */
	int optimal_l;					/* lth position (dimension 1) of optimal score */
	int optimal_m;					/* mth position (dimension 2) of optimal score */
	int optimal_n;					/* nth position (dimension 3) of optimal score */
	int l_alignment;				/* length of multiple alignment (max = l_seq1+l_seq2+l_seq3)*/
	int **alignment;				/* stores actual multiple alignment as sequence of matched 
										triplets */
	int *seq1;						/* stores sequence 1 */
	int *seq2;						/* stores sequence 2 */
	int *seq3;						/* stores sequence 3 */
} MA_data, *MA;
typedef struct Multiple_Alignment Mult_A;

struct project_files {
	FILE *f_sm;
	char f_sm_name[MAXFILENAME];
	FILE *f_seq;
	char f_seq_name[MAXFILENAME];
	FILE *f_qseq;
	char f_qseq_name[MAXFILENAME];
	FILE *f_dbseq;
	char f_dbseq_name[MAXFILENAME];
	FILE *f_scores;
	char f_scores_name[MAXFILENAME];
} file_data,*file;
typedef struct project_files p_files; 

/* declare function prototypes */

extern int read_sm(s_fxn *sfxn,FILE *sm_file);

extern int read_seq(sqs *seq, FILE *seq_file, int qordb);

extern int aachar2integer(char aa, amino *acids);

extern char integer2aachar(int integer, amino *acids);

extern void open_files(p_files *file);

extern void dynamic_program(/*Dynamic *DP,*/int *q_seq,int *db_seq,int l_qseq,int l_dbseq,int rs_or_zs,int twod_or_3d);

extern int alloc_init_matrix(Dynamic *DP,s_fxn *sfxn,opshuns *op,int l_qseq,int l_dbseq);

extern void create_matrix(Dynamic *DP,s_fxn *sfxn,amino *acids,opshuns *op,sqs *seq,int *q_seq,int *db_seq,int l_qseq,int l_dbseq);

extern void get_optimal_score(Dynamic *DP,opshuns *op,sqs *seq,int rs_or_zs);

extern void get_optimal_alignment(Dynamic *DP,amino *acids,int *q_seq,int *db_seq,int l_qseq,int l_dbseq);

extern void prnt_alignment(Dynamic *DP,amino *acids,sqs *seq,FILE *a_file);

extern void prnt_T(Dynamic *DP,sqs *seq,FILE *a_scores,opshuns *op);

extern void command_line(int argc, char *argv[],opshuns *op);

extern void get_z_score(Zscore *zs,Dynamic *DP,opshuns *op,int *q_seq,int *db_seq,int l_qseq,int l_dbseq);

extern void randomize_qseq(Zscore *zs,opshuns *op,int *q_seq,int l_qseq);

extern void calc_z_score(Zscore *zs,Dynamic *DP,opshuns *op,sqs *seq);

extern void prnt_Z(Zscore *zs,Dynamic *DP,sqs *seq,FILE *a_scores,opshuns *op);

extern void three_seq_alignment(Mult_A *MA);

extern void get_three_sequences(Mult_A *MA,sqs *seq);

extern void setup_3d_matrix(Mult_A *MA,Dynamic *DP,sqs *seq,s_fxn *sfxn,amino *acids,int *seq1,int l_seq1,
					 int *seq2,int l_seq2,int *seq3,int l_seq3);

extern void create_3d_matrix(Mult_A *MA,s_fxn *sfxn,amino *acids,int *seq1,int l_seq1,
					 int *seq2,int l_seq2,int *seq3,int l_seq3);

extern void get_optimal_3d_score(Mult_A *MA,int l_seq1,int l_seq2,int l_seq3);

extern void get_mult_alignment(Mult_A *MA,amino *acids,int *seq1,int l_seq1,
					 int *seq2,int l_seq2,int *seq3,int l_seq3);

extern void prnt_mult_alignment(Mult_A *MA,amino *acids,FILE *a_file);
