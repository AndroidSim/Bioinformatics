/*	CS 626: Project 2
	header file for project2.c

	written by Andrew Smith
	4/20/01
*/

/* use some standard libraries */
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fcntl.h"
#include "time.h"
#include "ctype.h"

/* some definitions */

#define TRUE 1
#define FALSE 0
#define MAXFILENAME 20
#define NSCA 14 

/* declare project structures */
struct options {
	int local;					/* option for performing local alignment (default global)*/
	int n_shuffles;				/* number of shuffled sequences for calculating z-score */
	int mult_align;				/* option for performing multiple alignment on sequence 6,7,9*/
	int contact;                /* if contact=1 use contact potential, if 0 use profile */
	int threading;              /* if doing sequence 2 structure alignments */
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
	int *leg_seq;
	int l_leg;
	int l_myo;
	int *myo_seq;
} seq_data,*seq; 
typedef struct sequences sqs;

struct structures {
	float *sc_geocenters;
	float *leg_sc_geocntr;
	float *myo_sc_geocntr;
	int *leg_contact_map;
	int *myo_contact_map;
	int *leg_contact_rep;
	int *myo_contact_rep;
	int *leg_contact_types;
	int *myo_contact_types;
	int n_cont_cat;
} struc_data, *struc;
typedef struct structures str;

struct scoring_fxn {
	float *scoring_matrix;			/* stores scoring function (profile or contact) */
	float gap_penalty;			/* stores constant gap penalty (= +1) */
} sfxn_data,*sfxn;					
typedef struct scoring_fxn s_fxn;

struct Dynam_Program {				/* for 2d alignments */
	float *dynamic_matrix;		/* stores dynamic matrix or table */
	float *T_i_j;				/* stores matrix of optimal scores for all sequences 
										aligned */
	float optimal_score;		/* stores optimal score for one alignment */
	float opt_rand_score;		/* stores optimal score for alignment with random sequence*/
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
	FILE *f_pdb1;
	char f_pdb1_name[MAXFILENAME];
	FILE *f_pdb2;
	char f_pdb2_name[MAXFILENAME];
	FILE *f_qcoord;
	char f_qcoord_name[MAXFILENAME];
	FILE *f_dbcoord;
	char f_dbcoord_name[MAXFILENAME];
	FILE *f_contact;
	char f_contact_name[MAXFILENAME];
	FILE *f_crep;
	char f_crep_name[MAXFILENAME];
	FILE *f_ctypes;
	char f_ctypes_name[MAXFILENAME];
	FILE *f_cener;
	char f_cener_name[MAXFILENAME];
	FILE *f_profile;
	char f_profile_name[MAXFILENAME];
} file_data,*file;
typedef struct project_files p_files; 

/* declare function prototypes */

extern int read_sm(s_fxn *sfxn,FILE *sm_file,amino *acids,str *struc);

extern int read_seq(sqs *seq, FILE *seq_file, int qordb);

extern int aachar2integer(char aa, amino *acids);

extern char integer2aachar(int integer, amino *acids);

extern void open_files(p_files *file);

extern void dynamic_program(/*Dynamic *DP,*/int *q_seq,int *db_seq,int l_qseq,int l_dbseq);

extern int alloc_init_matrix(Dynamic *DP,s_fxn *sfxn,opshuns *op,str *struc,int l_qseq,int l_dbseq);

extern void create_matrix(Dynamic *DP,s_fxn *sfxn,amino *acids,opshuns *op,sqs *seq,str *struc,
						  int *q_seq,int *db_seq,int l_qseq,int l_dbseq);

extern void get_optimal_score(Dynamic *DP,opshuns *op,sqs *seq);

extern void get_optimal_alignment(Dynamic *DP,amino *acids,int *q_seq,int *db_seq,int l_qseq,
								  int l_dbseq);

extern void prnt_alignment(Dynamic *DP,amino *acids,sqs *seq,FILE *a_file);

extern void prnt_T(Dynamic *DP,sqs *seq,FILE *a_scores,opshuns *op);

extern void command_line(int argc, char *argv[],opshuns *op);
/*
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
*/
extern void read_data_write_files(p_files *file,sqs *seq,str *struc,amino *acids);

extern void read_PDB_calc_gc(amino *acids,FILE *f_pdb,char *file_name,sqs *seq,str *struc);

extern char match_3letaa(char *token);

extern void calc_contact_map(str *struc,sqs *seq,int l_seq,int *struc_seq,float *geo_centers,
							 int *contact_map);

extern void calc_contact_rep(str *struc,sqs *seq,int l_stru,int *struc_seq,int *contact_rep,
							 int *contact_map);

extern void calc_contact_types(str *struc,sqs *seq,amino *acids,int l_stru,int *struc_seq,
							   int *contact_map,int *contact_type);

extern void allocate_struc_data(str *struc,amino *acids,int l_leg_stru,int l_myo_stru);

extern void printintmatrix2file(FILE *file, int *matrix,int first_dim,int second_dim);

extern void printintvector2file(FILE *file,int *vector,int dim);

extern void close_files(p_files *file);

extern float calc_gap_penalty(str *struc,s_fxn *sfxn,int struc_site,int leg_or_myo);

extern float calc_energy(sqs *seq,str *struc,amino *acids,s_fxn *sfxn,int alpha,int struc_site,int leg_or_myo);