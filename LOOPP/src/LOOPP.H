/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000    v. 2.000                                   */
/*                                                                */
/*                   DEFINITIONS                                  */   

/* standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include "Time.h" 

/* prefixed parameters to setup the defaults */
# define VERSION   2.000
# define NTYPE        20 /* default number of letters in alphabet */
# define NTHOM1       11 /* default number of site types for ein */
# define NTHOM2       15 /* default number of contact types for einm */
# define JUMP          1 /* step to move along scaffold while threading */
# define DELTA         4 /* chem dist excluded from contact counting */
# define MAX        3000 /* used for mps */
# define MAXS        256 /* to set length of file names and other strings */
# define MAXL         60 /* line length for printing of alignments */
# define R_low       1.0 /* default cutoff distances */
# define R_high      6.4
# define RMAX      999.9 /* to put away residues of unknown coordinates */
# define RESC_VDW    3.0 /* to scale distances for vdw */
# define epsilon    1e-6 /* Ax > epsilon for LP contraints */
# define GAP_PNLT    8.0 /* default gap penalty for alignment (BLOSUM50) */
# define GAP_PNLT_EIJ  0.2 /* default gap penalty for frozen env. (TE pot) */
# define SHORT        20 /* to skip sequences shorter than SHORT */
# define FALSE         0
# define TRUE          1
# define L_RATIO_0  0.75 /* to narrow search by length (depth=0) */
# define L_RATIO_1   0.5 /* to narrow search by length (depth=1) */
# define L_RATIO_2  0.25 /* to narrow search by length (depth=2) */
# define L_RATIO_3   0.1 /* to narrow search by length (depth=0) */
# define TH_EMPTY   20.0 /* to skip empty structures */

/* TYPE DECLARATIONS */

/* options will be kept in the following:    */
struct keep_options
 {
   int eij,ein,einm,evdw; /* different models of the potential */ 
   int wDc; /* wDc stands for write Difference in contacts vector */
   int mps; /* to write Difference in contacts in MPS format */
   int noth; /* noth means that there is no threshold for DE */
   int up; /* th_ene is a lower bound if(up) */
   int n_par; /* number of parameters in the Model */
   int gapless; /* to define gapless threading for learn and decoy */
   int gaps; /* to allow gaps in threading for exam */
   int nogaps; /* the opposite for gaps - to be used for exam */
   int cdel; /* to use constant penalty for deletions */
   int strgap; /* to use env dependent gap penalties in seq-seq align */
   int r_cont; /* to allow reading of existing f_cont file */
   int info; /* to get general info with potential used, distributions etc */
   int info2; /* to get general even more detailed info */
   int info3; /* to get general even more detailed info */
   int info4; /* to get general even more detailed info */
   int beg,end; /* to define which seq are to be thread */
   int mat; /* to enforce matrix format of potential file */
   int rep_pow, atr_pow; /* powers for vdw-like potentials */
   int distance; /* to use distance map */
   int col,n_col; /* to choose C_alpha, C_sidech, C_beta 3-col blocks in XYZ */
   int def_pot; /* to use internal potentials */
   /* to choose the mode of operation */
   int learn; /* to make standard threading loop and get inequalities */
   int exam; /* to pick up best plausible stru for an independent seq */
   int decoy; /* to get inequalities for explicitly generated decoys */
   int crd2xyz; /* to get reduced representation from CRD file */
   int pdb2xyz; /* to get reduced representation from PDB file */
   int fngrps; /* to use additional structural fingerprints */
   int debug; /* to debug the code - get addt prints */
   /* to allocate arrays */
   int max_length; /* protein size limit (in terms of residues) */
   int max_cont; /* number of contacts (per residue) limit */
   int max_prot; /* protein number limit (size of the set of scaffolds) */
   int max_histo; /* histograms size limit  */
   int n_type; /* to define how many letters are to be used */
   /* for seq-seq, seq-stru and stru-stru alignment */
   int align; /* flag to initiate seq-seq alignment */
   int strucmp; /* flag to initiate stru-stru alignment */
   int loc_a; /* to perform local alignment */
   float t_length; /* length threshold for local alignments */
   float t_seqid; /* sequence identity threshold for redundancies */
   float t_zscore; /* Zscore threshold for printing the alignments */
   float t_rms; /* RMS threshold for printing the alignments */
   float t_rms_db; /* RMS threshold for database building */
   int onelet; /* one letter code for ama */
   int gap_adr; /* address of the GAP res in the current alph */
   float pref_pen; /* gap penalties: prefix/suffix */
   float gap_pen; /*                 and regular   */
   float th_ene; /* threshold for energy differences to be printed */
   float r_cut; /* contact distance cutoff */
   float resc_pot; /* to rescale (e.g. by -1) the potential function */
   int freeze; /* flag to enable freezing of variables */
   int nfrozen; /* to specify number of frozen variables */
   int kbest; /* reshuffling for kbest best alignments */
   int k_prn; /* print only k_prn best matches  */
   int n_shuffle; /* do the reshuffling n_shuffle times */
   int prn_shuffle; /* to print alignments for randomly shuffled seqs */
   int prnalg_lev; /* to define format for printing the alignments */
   int qfmt; /* to choose the query_seq format */
   char qname[MAXS]; /* query sequence name */
   int pick; /* to pick one query seq from the query file (or str from XYZ) */
   int pick_all_seq; /* to pick all the query sequences */
   char pseq_name[MAXS]; /* picked sequence name */
   char pstr_name[MAXS]; /* picked structure name */
   int noredund; /* to avoid redundant domains in a pdb file */
   int rm_redund; /* remove redundancies from the database */
   int build_db; /* prepare list for a new database */
   int list; /* to use query file to specify file names only */
   int search_depth; /* to define depth of the search */
   int mod_search_depth; /* to indicate that defaults were altered */
   int justq; /* to stop after building query files */
   int iprobe; /* probe seq/stru index */
   int ver; /* to print signature */
 } ;
typedef struct keep_options k_opt;

struct all_files
 {
   FILE *coor,*out1,*out2,*cont,*info,*iprot; 
   char f_coor[MAXS],f_out[MAXS],f_dcon[MAXS],f_cont[MAXS],f_info[MAXS],
        f_iprot[MAXS];
   FILE *seq,*type,*pot,*mod,*mps,*exam,*query,*qxyz,*qcont,*qrij;
   char f_seq[MAXS],f_type[MAXS],f_pot[MAXS],f_mod[MAXS],f_mps[MAXS],
        f_exam[MAXS],f_query[MAXS],f_qxyz[MAXS],f_qcont[MAXS],f_qrij[MAXS];
   FILE *vdw_a,*vdw_b,*rij,*fngrps,*rhs,*best,*balg,*list,*newlist;
   char f_vdw_a[MAXS],f_vdw_b[MAXS],f_rij[MAXS],f_fps[MAXS],f_rhs[MAXS],
        f_best[MAXS],f_balg[MAXS],f_list[MAXS],f_newlist[MAXS];  
   char f_crd[MAXS],f_pdb[MAXS]; /* these two guys are obsolate */
   FILE *fngrps2,*fngrps3,*E_fngrps;
   char f_fps2[MAXS],f_fps3[MAXS],E_fps[MAXS];
   FILE *q_fngrps,*q_fngrps2,*q_fngrps3;
   char f_qfps[MAXS],f_qfps2[MAXS],f_qfps3[MAXS];
   FILE *fps_dv,*fps_norm,*Sfile;
   char f_dv[MAXS],f_norm[MAXS],f_entropy[MAXS];
 } ;
typedef struct all_files a_file;

/* description of the residue environmnet in terms of number of neighbours */
struct solvation_shells 
 {
   int *Ncont; /* first solvation shell - counts neighb. of subsequent res i */
   int *Mcont0; /* second solv shell - counts residues  with small */
   int *Mcont1; /* medium */
   int *Mcont2; /* large number of contacts being a neighb. of i */ 
 } ;
typedef struct solvation_shells solv_sh;

/* for tables of strings of fixed length (MAXS) */
struct my_strings
 {
  char mstr[MAXS]; 
 } ;
typedef struct my_strings str_my;

  
struct keep_info 
 {
   unsigned long int *histo,index,n_prn; /* for energy gaps stat */
   unsigned long int ncont,nres,*res_Stat;
   double ncont_decoys,*con_All_decoys;
   unsigned long int *con_all,*con_all_dev; /* for contact statistics */
   int *Point_to_best,ntop;
   int *Best_str,*toBest;
   float *Best_ene,*Zscore1,*Zscore2,*Score_shuffle; /* for best alignments */
   double av_ene,z_score; /* for energy gaps statistics */
   char fmt[MAXS]; /* format for f_out */
   str_my *Names; /* to store names of the proteins */
   int *noRed; /* to indicate redundant structures in the DB */
   char Sout[MAXS]; /* to store things to be printed */
 } ; 
typedef struct keep_info k_info;

struct keep_mps {
  int value;
  int row;
  struct keep_mps *down;
}  ; 
typedef struct keep_mps k_mps;

struct keep_residue 
{
  int num;
  int n_sc;
  float x_a,y_a,z_a,x_sc,y_sc,z_sc,x_b,y_b,z_b;
  char name[MAXS];
  char cnum[MAXS];
  int secstr;
} ;
typedef struct keep_residue k_res;

struct keep_proteins {
  float *X_res,*Y_res,*Z_res; /* coor of residues in red. repr. */
  float *X_q,*Y_q,*Z_q; /* copy for rmsd */
  int *seq; /* sequence */
  int *native_Seq;
  int *query_Seq; /* arrays storing sequences */
  int *FngrPs; /* structural fingerprints */
  float *R_ene; /* native energies */
  int *IJcont; /* ContactMap */
  float *IJr12,*IJr6,*native_IJrij,*query_IJrij; /* continuous ContactMaps */
  int *conR,*conT; /* protein in terms of contact types */
  float *conRvdw,*conTvdw;  /* protein in terms of cont. contact types */
  k_res *Res;
} ; 
typedef struct keep_proteins k_prot;

struct keep_dynprog_Table {
  float  *Table; /* dynamic programming table */
  int *Trace; /* to define genealogy */
  int *aGaps; /* to trace gaps in alignment */
  int *iniCol; /* dpt column to keep seq/stru to be aligned */
  int *iniRow; /* dpt row to keep seq/stru to be aligned */
  float *C_score,*R_score; /* site scores for stru-stru alignment */
} ; 
typedef struct keep_dynprog_Table k_dynpt;

struct keep_Model {
  float  *Potn; /* to keep potential or scoring function */
  int *Imod; 
  int *Alph,*Adr; /* to define alphabet different from generic */
} ; 
typedef struct keep_Model k_Model;

struct fea_mod {
   int *env_eij; /* stores protein cm according to aa types at each stru site */
   int *tn_cont; /* total number of contacts at a stru site */
   int *flag; /* used to see how much of *env_eij is accessed or used */
   float *r12_env; /* stores protein cm according to aa types at each stru site */
   float *r6_env;  /* for continuous potential and contact map */
} fea_data,*fea; 
typedef struct fea_mod f_mod;

struct fingerprnts {
	int fpsn_type; /* number of aa types used for fingerprints or profiles */
	int fps_potn; /* option for either using eij(==0),evdw(==1),ein(==2),einm(==3) potential */
	double norm_d0; /* parameter used in scoring fxn for stru2stru alignments */
	float *x2x_Potn; /* stores potential used to calculate fingerprints */
	double *s_energy; /* stores energy at each stru site for all aa types according to fea */
	double *profile; /* stores fingerprint or profile calculated using *s_energy */
	float *C_profile; /* stores fingerprint for C_name for stru2stru alignments */
	float *R_profile; /* stores fingerprint for R_name for stru2stru alignments */
} fps_data,*fps;
typedef struct fingerprnts fgrprnt;


/* FUNCTION DECLARATIONS */
extern void prn_signature(int argc, char *argv[], a_file *F, k_opt *op, 
			  k_info *Info);
extern void check_opt(a_file *F, k_opt *op, k_info *Info);
extern void interpret_Model(a_file *F, k_opt *op, k_info *Info);
extern void command_line(int argc, char *argv[], a_file *F, k_opt *op, 
			 k_info *Info);
extern void set_default_opt(a_file *F, k_opt *op, k_info *Info);
extern void set_file_names(a_file *F, k_opt *op, k_info *Info);
extern void exit_end(char *Str_to_out, a_file *F, k_prot *prot, k_opt *op,
		     k_info *Info, k_Model *model, solv_sh *ssh, k_dynpt *dpt);
extern void exit_err(char *Str_to_out, a_file *F, k_prot *prot, k_opt *op,
		     k_info *Info, k_Model *model, solv_sh *ssh, k_dynpt *dpt);
extern int generic2adr(int *i, int *n_type, int *Alph, int *Adr);
extern int mat2vec(int *i, int *j, int *n_type);
extern int read_seq(k_opt *op,FILE *Fseq, int *Lseq, char *T_name, int *seq, int *i_prot,
		    int *n_type, int *Alph, int *Adr);
extern int read_seq_generic(FILE *Fseq, int *Lseq, char *T_name, int *seq,
			    int *i_prot, k_opt *op);
extern int get_length_from_seq(FILE *Fseq, int *Lseq, char *T_name, int *Lmax,
			       int *n_res, int *n_prot, k_opt *op);
extern int seq2dig(char *s);
extern int one2dig(char *s);
extern char *dig2amino(int i);
extern char *dig2one(int i);
extern int read_coor(FILE *Fin, k_opt *op, int *Lstru, char *R_name,
		     float *X, float *Y, float *Z);
extern int read_coor_all3(FILE *Fin, k_opt *op, int *Lstru, char *R_name,
			  float *X, float *Y, float *Z, float *Xa, float *Ya, 
			  float *Za, float *Xb, float *Yb, float *Zb);
extern int alloc_mem(a_file *F, k_prot *prot, k_opt *op, k_info *Info,
		     k_Model *model, solv_sh *ssh, k_dynpt *dpt);
extern void clean_up(a_file *F, k_prot *prot, k_opt *op, k_info *Info,
		     k_Model *model, solv_sh *ssh, k_dynpt *dpt);
extern void get_contMap(int *IJ, int *Lstru, char *R_name,
			float *X,float *Y, float *Z, k_opt *op);
extern int read_contMap(a_file *F, k_opt *op, int *IJ, int *Lstru,
			char *R_name, float *r12ij, float *r6ij, float *rij);
extern void write_contMap(a_file *F, int *IJ, int *Lstru, char *R_name,
			  float *X, float *Y, float *Z, k_opt *op,
			  float *r12ij, float *r6ij, float *rij);
extern void write_Query_CM(a_file *F, int *Lstru, char *R_name,
			   float *X, float *Y, float *Z, k_opt *op);
extern void get_contType(k_opt *op, int *Imod, int *contT, int *IJ, int *start,
		  int *seq,int *Lseq,int *start2,int *Icont,
		  solv_sh *ssh,float *contTvdw,float *r12ij,float *r6ij);
extern void get_env(k_opt *op, int *Imod, int *contT, int *IJ, int *start,
		    int *seq, int *Lstru, int *start2, int *Icont,
		    solv_sh *ssh, float *r12ij, float *r6ij, int *Stru_row);
extern float get_score(k_opt *op, int *Imod, int *n, int *seq,
		       int *Lseq, int *m, int *Icont, float *Potn,
		       int *Mcont0, int *Mcont1, int *Mcont2, float *gpen);
extern void write_contType(a_file *F,k_opt *op,int *Imod,int *contR,int *IJ,
		    int *start,int *seq,int *Lseq,char *R_name,int *Icont,
		    solv_sh *ssh,float *contRvdw,float *r12ij,float *r6ij);
extern int read_contType(a_file *F, k_opt *op, int *contR, char *T_name,
			 float *contRvdw);
extern void get_einm(k_opt *op, int *Imod, int *contR, int *IJ_loc, int *start,
		     int *seq, int *Lseq, int *Ncont, int *Mcont0, int *Mcont1,
		     int *Mcont2);
extern void get_env_ein(k_opt *op, int *Imod, int *IJ, int *Lstru, int *Icont);
extern void get_env_einm(k_opt *op, int *Imod, int *contR, int *IJ_loc,
			 int *start, int *seq, int *Lstru, int *Ncont, 
			 int *Mcont0, int *Mcont1, int *Mcont2);
extern float energy(k_opt *op, int *contR, float *Potn, float *contRvdw);
extern void read_pot(k_opt *op, float *Potn, a_file *F, int *Imod);
extern void write_ineq(a_file *F, k_opt *op, float *Potn, int *contR, 
		       int *contT, unsigned long int *index, float *contRvdw,
		       float *contTvdw, k_mps *var_mps[], k_mps *tail_mps[]);
extern void space(FILE *out, int v, int n);
extern void write_mps(a_file *F, k_opt *op, unsigned long int *index,
		      k_mps *var_mps[]);
extern void print_help(char *argv[], k_info *Info);
extern void info(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, int *Imod,
		 int *Alph, int *Adr);
extern void get_stat_poten(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, 
			   int *Imod, int *Alph, int *Adr);
extern void info_pot(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, 
		     int *Imod, int *Alph, int *Adr);
extern void info_cont_nat(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, 
			  int *Imod, int *Alph, int *Adr);
extern void info_cont_dec(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, 
			  int *Imod, int *Alph, int *Adr);
extern void info_ene_tot(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, 
			 int *Imod, int *Alph, int *Adr);
extern void info_effpair_einm(FILE *Finfo, k_info *Info, float *Potn, 
			      k_opt *op, int *Imod, int *Alph, int *Adr);
extern void prn_res(k_opt *op, a_file *F, float *Potn, int *contR, int *contT,
		    unsigned long int *index, float *T_ene, char *R_name,
		    char *T_name, k_info *Info, int *n, int *Ls,
		    float *contRvdw, float *contTvdw,
		    k_mps *var_mps[], k_mps *tail_mps[]);
extern void prepare_files(a_file *F, k_opt *op, k_info *Info, float *Potn); 
extern void set_options(int argc, char *argv[], a_file *F, k_opt *op, 
			k_info *Info); 
extern void set_Model(a_file *F, k_opt *op, k_info *Info, int *Imod, 
		      int *Alph, int *Adr, float *Potn); 
extern void prn_setup(a_file *F, k_opt *op, int *Alph); 
extern void report_progress(k_info *Info, int *wind, char *DB_name, 
			    int *l_stru);
extern int pick_structure(k_opt *op, char *DB_name, int *db_len, int *q_len);
extern void adjust_gap_pen(k_opt *op, int *Ncont, int *j, float *gap_pen_left,
			   float *gap_pen_top, float *gpen_prev, float *gpen);
extern void adjust_prn_level(k_opt *op);
extern void get_alignment(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot, 
			  int *Imod, int *query_Seq, int *l_qseq, k_info *Info,
			  int *Alph, int *Adr, float *Table, int *Trace, 
			  int *aGaps, int *n_row, int *n_col, int *ind, 
			  int *wind, int *db_Seq, char *DB_name, char *Q_name);
extern void get_row_eij(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot, 
			int *Imod, float *Potn, int *db_Seq, int *Ncont, 
			k_info *Info, float *Table, int *Trace, int *i_row, 
			int *n_col, int *alpha, float *gap_pen_left, 
			float *gap_pen_top);
extern int read_Mod(FILE *Fmod, k_opt *op, int *Imod, int *Alph, int *Adr);
extern int read_head_of_Mod(FILE *Fmod, k_opt *op);
extern int imod2vec(int *i, int *n, int *Imod);
extern void read_err(char *name);
extern void open_err(char *name);
extern void alloc_err(char *name);
extern void build_irrep(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
			k_info *Info, k_Model *model);
extern int check_CM(k_opt *op, a_file *F, k_info *Info, int *IJcont, int *Seq,
		    float *r12ij, float *r6ij, float *rij, int *Alph, int *Adr);
extern int switch_rw_CM(k_opt *op, a_file *F, k_info *Info, int *to_read);
extern void prep_info(k_info *Info, float *T_ene, int *max_histo);
extern void add_res(k_info *Info, int *Seq, int *l_seq, int *Adr, int *Alph);
extern void add_cont(int *n_para, int *contR, unsigned long int *n_cont,
		     unsigned long int *cont_all,
		     unsigned long int *cont_all2);
extern void add_cont_decoy(int *n_para, int *contT, double *x_cont,
			   double *cont_All);
extern void prn_title_bar(FILE *Finfo, char *Title, k_opt *op, int *Alph);
extern void init_info(k_opt *op, a_file *F, k_info *Info,
		      k_mps *var[], k_mps *tail[]);
extern void info_prot(FILE *Finfo, char *R_name, float *ene, float *Potn,
		      k_opt *op,int *Imod,int *contR,int *Alph,int *Adr);
extern void skip(a_file *F, k_opt *op, int *seq, int *i_prot, int *contR,
		 float *contRvdw, int *Alph, int *Adr);
extern int pick(char *amino, char *atom);
extern int pick_a(char *amino, char *atom);
extern int pick_b(char *amino, char *atom);
extern void get_ineq(k_opt *op,a_file *F,solv_sh *ssh,int *contR,int *contT,
		     int *IJcont,int *Imod,float *Potn,float *R_ene,int *seq,
		     k_info *Info,float *contRvdw,float *contTvdw,
		     float *r12ij,float *r6ij,float *rij,int *Alph,int *Adr,
		     k_mps *var_mps[],k_mps *tail_mps[]);
extern void get_ineq_dec(k_opt *op, a_file *F, solv_sh *ssh, int *contR,
			 int *contT, int *IJcont, int *Imod, float *Potn,
			 float *R_ene, int *seq, k_info *Info,
			 float *contRvdw, float *contTvdw, int *Alph, int *Adr,
			 k_mps *var_mps[], k_mps *tail_mps[]);
extern void check_seq(k_opt *op,a_file *F,solv_sh *ssh,int *contR,int *contT,
	       int *IJcont,int *Imod,float *Potn,float *R_ene,
	       k_prot *prot,k_info *Info,float *contRvdw,float *contTvdw,
	       float *r12ij,float *r6ij,float *rij,
	       int *Alph,int *Adr,float *Table,
	       int *Trace,int *aGaps,int *Stru_row);
extern void thread_seq(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
		       k_info *Info, k_Model *model, k_dynpt *dpt,
		       k_mps *var_mps[], k_mps *tail_mps[]);
extern int stru2stru(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
		     k_info *Info, k_Model *model, k_dynpt *dpt);
extern int get_best_rms(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
			k_info *Info, k_Model *model, k_dynpt *dpt, int *C_seq,
			int *l_qstru, float *X_q, float *Y_q, float *Z_q);
extern float get_rms(k_opt *op, a_file *F, k_info *Info, float *R_ene,
		     int *i_str, int *n, int *i, int *j, int *aGaps, 
		     float *C_rij, float *R_rij, int *n_row, int *n_col,
		     float *X_q, float *Y_q, float *Z_q, 
		     float *X, float *Y, float *Z);
extern float rmsd(k_opt *op, int *n, int *i, int *j, int *aGaps, int *n_row, 
		  int *n_col, float *X_q, float *Y_q, float *Z_q, 
		  float *X, float *Y, float *Z);
extern void comp_stru(k_opt *op,a_file *F,solv_sh *ssh,int *contR,int *contT,
		      int *IJcont,int *Imod,float *Potn,float *R_ene,int *seq,
		      k_info *Info,float *r12ij,float *r6ij,float *rij,
		      int *Alph,int *Adr,float *Table,
		      int *Trace,int *aGaps,int *Stru_row,int *Stru_col,
		      float *C_score,float *R_score);
extern void align_seq(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
		      k_info *Info, k_Model *model, k_dynpt *dpt);
extern void align_seq_old(k_opt *op,a_file *F,int *Imod,float *Potn,int *seq,
			  k_info *Info,int *Alph,int *Adr,int *Seq_col,
			  float *Table,int *Trace,int *aGaps);
extern float get_cell(k_opt *op, float *Table, int *Trace, int *n_col, int *i,
		      int *j, float *gap_pen_left, float *gap_pen_top,
		      float *Tij);
extern float find_best_successor(k_opt *op,float *Table,int *n_row,int *n_col,
			  int *n,int *m);
extern void prn_dpt(float *Table, int *n_row, int *n_col);
extern int trace_back(k_opt *op, float *Table, int *n_row, int *n_col,
		      int *i, int *j, float *T_ene, int *Trace, int *aGaps);
extern void prn_align_old(k_opt *op,a_file *F,int *seq,int *Alph,int *Seq_col,
	       int *n,int *i,int *j,float *T_ene,int *aGaps,char *T_name,
	       char *R_name,int *n_row,int *n_col);
extern float prn_align(k_opt *op,FILE *Fout,k_info *Info,int *Seq_row,
	       int *Alph,int *Seq_col,int *n,int *i,int *j,float *T_ene,
	       int *aGaps,char *T_name,char *R_name,int *n_row,int *n_col);
extern void init_dpt(k_opt *op, float *Table, int *Trace, int *n_row,
		     int *n_col, int *aGaps);
extern int get_header_crd(int Fcrd);
extern int read_crd(int turn, int Fcrd, FILE *Fxyz, k_opt *op,
	     float *X,float *Y,float *Z,int *seq);
extern void read_line(int inp,char *w);
extern void get_red_repr(a_file *F, k_opt *op, float *X, float *Y, float *Z,
			 int *seq);
extern void prn_seq_alph(int *seq, int *Lseq, char *Tname, int *Alph);
extern void one2three(a_file *F, k_opt *op);
extern void pick_best_aligns(k_opt *op, k_info *Info, int *no_str, 
			     float *T_ene);
extern void pick_best_zscore(k_opt *op,k_info *Info);
extern void shuffle_seq(k_opt *op, a_file *F, solv_sh *ssh, int *contR,
			int *contT, int *IJcont, int *Imod, float *Potn,
			float *R_ene, int *nat_Seq, int *query_Seq,int *l_qseq,
			k_info *Info, float *contRvdw,
			float *contTvdw, float *r12ij, float *r6ij, float *rij,
			int *Alph,
			int *Adr, float *Table, int *Trace, int *aGaps,
			int *Stru_row);
extern void shuffle_seq_algn(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
			     int *Imod, float *Potn, int *query_Seq, 
			     int *l_qseq, k_info *Info, int *Alph, int *Adr,
			     float *Table, int *Trace, int *aGaps);
extern void reshuffle_seq(int *seq, int *Lseq);
extern void shuffle_rand(int *seq, int *Lseq, int *mult);
extern void print_best_matches(k_opt *op, a_file *F, k_info *Info);
extern void to_file(FILE *Fout, char *Str_to_out);
extern void to_stnd_out(char *Str_to_out);
extern void write_down_seq(a_file *F, k_opt *op, k_info *Info, int *Seq,
			   int *lseq);
extern void seq2query(a_file *F,k_opt *op,k_info *Info,int *Seq,int *lseq);
extern int build_query(a_file *F, k_opt *op, k_info *Info, k_prot *prot,k_Model *model);
extern void prepare_query_1let(a_file *F, k_opt *op, k_info *Info);
extern void prepare_query_fasta(a_file *F, k_opt *op, k_info *Info);
extern int prepare_query_swiss(a_file *F, k_opt *op, k_info *Info, int *Seq);
extern int prepare_query_pdb(a_file *F, k_opt *op, k_info *Info, int *Seq);
extern int prepare_xyz_pdb(a_file *F, k_opt *op, k_info *Info, k_res *Res,
			   int *Seq, int *n_domain);
extern float get_zscore(k_opt *op, a_file *F, k_info *Info, float *R_ene,
		       int *i_str, float *e_ave, float *sigma);
extern float get_score_eij(k_opt *op,int *Imod,int *contT,int *IJ,int *start,
			  int *nat_Seq,int *Lseq,int *qseq_i_alpha,int *j_beta,
			  float *contTvdw,float *r12ij,float *r6ij,float *Potn,
			  float *gpen);
extern void skip_query_warn(k_opt *op, k_info *Info);
extern void extend_name(k_opt *op, char *name, char *chain_name);
extern void get_clean_name(char *name, char *chain_name);
extern void get_xyz_col(char *line, int *beg, float *x, float *y, float *z,
			int *n_dig);
extern int get_line_pdb(a_file *F, k_opt *op, k_info *Info, k_res *Res,
			int *Seq, int *n_res, char *line);
extern void restart_counters(int *sum, int *i_let, int *error);
extern int line2seq(int *Seq, int *i_let, char *line, int *sum, 
		    int *read_next);
extern void write_xyz_line(a_file *F, k_opt *op, k_info *Info, k_res *Res,
			   int *i_res);
extern void init_Res(k_res *Res, int *n_res);
extern void check_Res(k_res *Res, int *n_res, int *n_sc_last);
extern int write_down_xyz(a_file *F, k_opt *op, k_info *Info, k_res *Res,
			  int *Seq, int *l_seq, int *n_res);
extern int get_triplet(int *s1, int *s2, int *s3, int *Seq, int *l_seq,
		       int *i_seq);
extern int match_res(int *s1, int *s2, int *s3, int *r1, int *r2, int *r3);
extern void get_default_pot(a_file *F, k_opt *op, float *Potn);
extern void tobi_elber_eij(float *Potn);
extern void meller_elber_tom2(float *Potn);
extern void meller_elber_tom1(float *Potn);
extern void meller_elber_tom2_av(float *Potn);
extern void meller_elber_vdw(float *Potn);
extern void blosum50(float *Potn);
extern int prepare_query_CM(a_file *F, k_opt *op, k_info *Info, k_prot *prot);
extern int xyz2query_CM(a_file *F, k_opt *op, k_info *Info, k_prot *prot, 
			int *i_prot, int *not_empty);
extern int read_Query_CM(a_file *F, k_opt *op, int *IJ, int *Lstru, 
			 char *R_name, float *rij);
extern double *vector(long nl, long nh);
extern void nrerror(char error_text[]);
extern void free_vector(double *v, long nl, long nh);
extern void jacobi(double **a, int n, double d[], double **v, int *nrot);
extern double **convert_matrix(double *a, long nrl, long nrh, long ncl, 
			       long nch);
extern int diag_kabsch(k_opt *op, double kab[3][3], double *eigval);
extern void free_convert_matrix(double **b, long nrl, long nrh, long ncl, 
				long nch);
extern void eigsrt(double d[], double **v, int n);
extern float rm_cent_mass(k_opt *op, int *l_stru, float *X, float *Y, float *Z);
extern void print_best_stru(k_opt *op, a_file *F, k_info *Info, float *rms_ave,
			    float *sigma);
extern void write_query_XYZ(a_file *F, k_opt *op, k_info *Info, int *lstru,
			    float *X, float *Y, float *Z, float *Xa, float *Ya, 
			    float *Za, float *Xb, float *Yb, float *Zb);
extern float get_hydrophob_ind(k_opt *op, k_info *Info, int *Seq, int *l_seq);
extern int check_for_membrane(k_opt *op, a_file *F, k_info *Info, int *Seq, 
			      int *l_seq);
extern int read_entry_list(a_file *F, k_opt *op, k_info *Info, char *name);
extern void query_file_warning(a_file *F, k_opt *op, k_info *Info);
extern void membrane_warning(a_file *F, k_opt *op, k_info *Info);
extern int prepare_query_crd(a_file *F, k_opt *op, k_info *Info,
			     k_prot *prot, int *Seq);
extern void add_name2list(k_opt *op, a_file *F, k_info *Info, char *name);
extern void prn_type_of_run(a_file *F, k_opt *op, k_info *Info);
extern float check_if_empty(int *Lstru, float *X, float *Y, float *Z);
extern float low_complexity(k_opt *op, a_file *F, k_info *Info, int *Seq, 
			    int *l_seq);

extern void get_env_eij(k_opt *op,f_mod *fea,int *Imod,int *IJ,int *nat_Seq,
						int *Lnatseq,int *nop,float *r12ij,float *r6ij,float *s_p_env);
extern float get_energy_eij(k_opt *op,f_mod *fea,int *Imod,int *IJ,int *contT,int *start,int *j_beta,
					 int *qseq_i_alpha,int *nat_Seq,int *Lnatseq,
                     float *Potn,float *contTvdw,float *r12ij,float *r6ij,
                     float *gpen,int *nop,int *nope);
/*extern void write_fngrps (a_file *F, k_opt *op,f_mod *fea,fgrprnt *fps,int *Imod,int *contR,int *IJ,
					int *seq,int *Lseq,int *Lstru,char *R_name,float *contRvdw,float *r12ij,
					float *r6ij,float *Potn); */
extern void write_fngrps(a_file *F, k_opt *op,f_mod *fea,fgrprnt *fps,int *IJ,
					int *seq,int *Lseq,int *Lstru,char *R_name,float *Potn,int quordata);
extern void init_fngrps(a_file *F, k_opt *op,fgrprnt *fps);
extern void free_up_fea(k_opt *op,f_mod *fea);
/*extern void read_fngrps(FILE *fpsf,char *fname,k_opt *op,float *fngrprnt,int *Lstru,char *X_name); */
extern int read_fngrps(FILE *fpsf,char *fname,k_opt *op,fgrprnt *fps,float *fngrprnt,
						float *rsc,int *Lstru,char *X_name,int quordata);
extern void set_fpp(a_file *F,k_opt *op,fgrprnt *fps);
extern void prepare_fps_files(a_file *F,k_opt *op);







