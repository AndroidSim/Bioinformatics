/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000      v. 2.000                                 */
/*                                                                */
/*                   ALIGNMENTS                                   */
 
#include "loopp.h"
   


/* align_seq is supposed to align two sequences using dynamic prog */
void align_seq(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
	       k_info *Info, k_Model *model, k_dynpt *dpt)
{
  int i_prot,n,m,i,j,k,n_row,n_col,alpha,beta,i_rewind,k_flag,
      ovrlap_flag,n_align;
  int i_best,change,i_rm,i_str,include,*db_Seq,*query_Seq,store_kbest;
  int pick_stru,lpick,ind,wind,l_qseq;
  int *Imod, *Alph, *Adr, *Trace, *aGaps, *IJcont, *Ncont;
  float *Potn, *Table, *r12ij, *r6ij, *R_rij;
  float A_ene,gap_pen,gap_pen_left,gap_pen_top,gpen_prev,gpen,Tij;
  char DB_name[MAXS],Q_name[MAXS];

  /* DB - database or reference (row) sequence; 
     Q - query or probe (column) sequence */ 

 /* define local shortcuts */
 db_Seq=prot->native_Seq;
 query_Seq=prot->query_Seq;
 IJcont=prot->IJcont;
 Imod=model->Imod;
 Potn=model->Potn;
 Alph=model->Alph;
 Adr=model->Adr;
 Table=dpt->Table;
 Trace=dpt->Trace;
 aGaps=dpt->aGaps;
 r12ij=prot->IJr12;
 r6ij=prot->IJr6;
 R_rij=prot->native_IJrij;
 Ncont=ssh->Ncont;
 store_kbest=op->kbest;

 /* print title bar */
 prn_type_of_run(F,op,Info);

 /* begin loop over sequences Q_name from f_query */
 while(read_seq(op,F->query,&l_qseq,Q_name,query_Seq,&i_prot,&op->n_type,Alph,Adr))
   {
     i_prot=0;
     ind=1; wind=0; 
     if(op->info2) 
       {
	 sprintf(Info->Sout,
		 "\n Sequence alignments ...  through (name_length): \n");
	 to_stnd_out(Info->Sout);
       }

     /* prepare tables for kth best alignments and shuffling */
     for(i=0;i<op->max_prot;++i) *(Info->Point_to_best+i)=0; 
     Info->ntop=0;
     op->kbest=store_kbest;
     
     /* begin loop over database sequences (DB_name) */
     while(read_seq(op,F->seq,&n_col,DB_name,db_Seq,&i_prot,&op->n_type,Alph,Adr))
       {  
	 ++wind; 
	 n_row=l_qseq;
	 /* to define environment-dependent gap penalties */
	 if(op->strgap) 
	   read_contMap(F,op,IJcont,&n_col,DB_name,r12ij,r6ij,R_rij);
	 
	 if(op->info2) report_progress(Info,&wind,DB_name,&n_col);
	 i_str=wind-1; /* so, in fact wind enumerates structures */ 
     
	 /* if a given structure does not match pick it will be skipped */
	 if(pick_structure(op,DB_name,&n_col,&l_qseq))
	   {
	     /* now get the environmnet (number of neigh.) of each site */
	     if(op->strgap) get_env_ein(op,Imod,IJcont,&n_col,Ncont); 

	     /* init Dpt and build initial row and col using pre/suf gap_pen */
	     init_dpt(op,Table,Trace,&n_row,&n_col,aGaps);

	     /* build the DP table */
	     gap_pen_left=op->gap_pen;
	     gap_pen_top=op->gap_pen;
	     /* seqs start at the second  r/c */
	     for(i=1;i<n_row;++i) 
	       {
		 alpha=*(query_Seq+i-1); 
		 /* build next row (loop over columns inside) */
		 get_row_eij(op,F,ssh,prot,Imod,Potn,db_Seq,Ncont,Info,
			     Table,Trace,&i,&n_col,&alpha,&gap_pen_left,
			     &gap_pen_top);
	       }

	     /* now find and print the optimal alignment */
	     get_alignment(op,F,ssh,prot,Imod,query_Seq,&l_qseq,Info,Alph,Adr,
			   Table,Trace,aGaps,&n_row,&n_col,&ind,&wind,db_Seq,
			   DB_name,Q_name);

	   } /* end if pick */
	 
       } /* end while db_Seq */

     sprintf(Info->Sout,"\n");
     to_stnd_out(Info->Sout);
     Info->index=ind-1;
     if(Info->ntop<op->kbest) op->kbest=Info->ntop;
 
     /* now get Z_scores for shuffled sequence */
     strcpy(op->qname,Q_name);
     if(op->debug) printf("%s\t%d \n",op->qname,l_qseq);
     shuffle_seq_algn(op,F,ssh,prot,Imod,Potn,query_Seq,&l_qseq,Info,Alph,Adr,
		      Table,Trace,aGaps); 
     
     rewind(F->seq);
     rewind(F->cont);

   } /* end while query_Seq */
 
  
}

/* once DP table is ready one may find the alignment */
void get_alignment(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot, int *Imod,
		   int *query_Seq, int *l_qseq, k_info *Info, int *Alph,
		   int *Adr, float *Table, int *Trace, int *aGaps,
		   int *n_row, int *n_col, int *ind, int *wind, int *db_Seq,
		   char *DB_name, char *Q_name)
{
  int include,i,j,n;
  float A_ene;

  /* print DP table */
  if(op->info4) prn_dpt(Table,n_row,n_col);
  
  /* trace back the alignment */
  A_ene=find_best_successor(op,Table,n_row,n_col,&i,&j); 
  /* so, now we are at the cell (i,j) having the value T_ene */
  n=trace_back(op,Table,n_row,n_col,&i,&j,&A_ene,Trace,aGaps); 
  
  /* prepare best alignments and estimate their significance */
  include=1;
  if(op->loc_a && (((float)n/(float)(*l_qseq))<op->t_length)) include=0;
  if(include) pick_best_aligns(op,Info,wind,&A_ene); 

  /* print aligned seqs - aGaps indicates how to move from (i,j)*/
  /* define first format of prinitng */
  op->prnalg_lev=1; /* no details */
  if(op->info) adjust_prn_level(op);

  prn_align(op,F->out1,Info,db_Seq,Alph,query_Seq,&n,&i,&j,
	    &A_ene,aGaps,Q_name,DB_name,n_row,n_col); 
  
  /* prepare additionally data for energy histograms */
  if(op->info) prep_info(Info,&A_ene,&op->max_histo);
  ++(*ind);
  
}

/* change the content of output */
void adjust_prn_level(k_opt *op)
{

  /* call it if any level of info */

  if(op->info3) op->prnalg_lev=3; /* vertical */
  else op->prnalg_lev=2; /* horizontal for -i and -i 2 */

  /* there is no point adding just scores for original */
  
}

/* build one row of the DP table */
void get_row_eij(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot, int *Imod, 
		 float *Potn, int *db_Seq, int *Ncont, k_info *Info, 
		 float *Table, int *Trace, int *i_row, int *n_col, int *alpha, 
		 float *gap_pen_left, float *gap_pen_top)
{
  int j,beta;
  float Tij,gpen_prev,gpen;  
  
  for(j=1;j<*n_col;++j) 
    {
      beta=*(db_Seq+j-1);  /* db_Seq goes to row */

      /* get the proper scoring matrix element */
      if(*alpha>beta) Tij=*(Potn+mat2vec(&beta,alpha,&op->n_type));
      else Tij=*(Potn+mat2vec(alpha,&beta,&op->n_type)); 

      /* define env dependent gap penalties */
      if(op->strgap) adjust_gap_pen(op,Ncont,&j,gap_pen_left,gap_pen_top,
				    &gpen_prev,&gpen);

      /* now decide how to proceed having Tij and gap_pen */
      *(Table+(*i_row)*(*n_col)+j)=get_cell(op,Table,Trace,n_col,i_row,&j,
					    gap_pen_left,gap_pen_top,&Tij);
    }  
}

/* shuffle_seq performs alignment for best matches only including shuffling */
void shuffle_seq_algn(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
		      int *Imod, float *Potn,
		      int *query_Seq, int *l_qseq, k_info *Info, int *Alph,
		      int *Adr, float *Table, int *Trace, int *aGaps)
{
  int i_prot,n,m,Lstru,Lseq,i,j,k,wind,n_row,n_col,alpha,beta;
  int i_best,change,i_rm,i_str,i_shuffle,l,lpick,*loc_Seq,*db_Seq;
  int *IJcont, *Ncont;
  unsigned long int ind;
  float A_ene,gap_pen_left,gap_pen_top,gpen,Tij,e_ave,sigma,ftmp,gpen_prev;
  float *r12ij, *r6ij, *R_rij, *R_ene;
  char DB_name[MAXS];

 /* shuffling loop -  align shuffled sequences  */
 db_Seq=prot->native_Seq;
 IJcont=prot->IJcont;
 Ncont=ssh->Ncont;
 r12ij=prot->IJr12;
 r6ij=prot->IJr6;
 R_rij=prot->native_IJrij;

 ind=1; wind=0; i_str=0; i_prot=0; 
 if(op->info2)
   {
     sprintf(Info->Sout,"\n Shuffling ...  \n");
     to_stnd_out(Info->Sout);
   }
 Lseq=*l_qseq;
 /* allocate array storing copy of the query_Seq */
 loc_Seq=malloc((Lseq+1)*sizeof(int));
 R_ene=malloc(op->max_prot*sizeof(float));
 if(R_ene==NULL) ;
 if(loc_Seq==NULL) ; /* indicate problem */ 
 for(i=0;i<Lseq;++i) *(loc_Seq+i)=*(query_Seq+i);
 
 
 /* begin loop over database sequences (DB_name) */
 rewind(F->seq);
 rewind(F->cont);
 while(read_seq(op,F->seq,&Lstru,DB_name,db_Seq,&i_prot,&op->n_type,Alph,Adr))
   {   
     if(op->strgap) 
       read_contMap(F,op,IJcont,&Lstru,DB_name,r12ij,r6ij,R_rij);
	 
     if(op->pick) /* skip all except picked if desired */
       {
	 lpick=strlen(op->pstr_name);
	 if(lpick>0 && strncmp(DB_name,op->pstr_name,lpick)!=0)
	   *(Info->Point_to_best+i_str)=0;
       }
     
     /* do the shuffling for good templates only */
     if(*(Info->Point_to_best+i_str))
       {

	 if(op->info2) 
	   {
	     sprintf(Info->Sout,"  %s_%d  ",DB_name,Lstru);
	     to_stnd_out(Info->Sout);
	   }
	 /* get the energy of native alignment */
	 *(R_ene+i_str)=0.0;
	 for(i=0;i<Lstru;++i)
	   {
	     alpha=*(db_Seq+i);
	     *(R_ene+i_str)+=(*(Potn+mat2vec(&alpha,&alpha,&op->n_type)));
	   }

	 /* now get the environmnet (number of neigh.) of each site */
	 if(op->strgap) get_env_ein(op,Imod,IJcont,&Lstru,Ncont); 
	 
	 /* shuffle original sequence */
	 for(i=0;i<20;++i) reshuffle_seq(query_Seq,&Lseq); 
	 e_ave=0.0;
	 sigma=0.0;
	 /* now one more time to regenerate original alignment */
	 for(i_shuffle=0;i_shuffle<=op->n_shuffle;++i_shuffle) 
	   {
	     n=0; m=0;
	     n_row=Lseq;
	     n_col=Lstru;
	     /* init Dpt, build initial row and col, fix pre, suf gap_pen */
	     init_dpt(op,Table,Trace,&n_row,&n_col,aGaps);
	     /* now n_row=Lseq+1, n_col=Lstru+1 */   
	     gap_pen_left=op->gap_pen;
	     gap_pen_top=op->gap_pen; 
	     /* shuffle again */
	     if(i_shuffle<op->n_shuffle) reshuffle_seq(query_Seq,&Lseq);
     
	     /* build the DP table */
	     for(i=1;i<n_row;++i) 
	       {
		 /* seqs start at the second  r/c */
		 if(i_shuffle<op->n_shuffle) alpha=*(query_Seq+i-1);
		 else alpha=*(loc_Seq+i-1);
		 /* build next row (loop over columns inside) */
		 get_row_eij(op,F,ssh,prot,Imod,Potn,db_Seq,Ncont,Info,
			     Table,Trace,&i,&n_col,&alpha,&gap_pen_left,
			     &gap_pen_top);
	       }
	    
	     /* print DP table */
	     if(op->debug) prn_dpt(Table,&n_row,&n_col);
		
	     /* trace back the alignment */
	     A_ene=find_best_successor(op,Table,&n_row,&n_col,&i,&j); 
	     /* so, we are at the cell (i,j) having the value T_ene */
	     n=trace_back(op,Table,&n_row,&n_col,&i,&j,&A_ene,Trace,aGaps); 

	     /* print alignment - aGaps defines how to move from (i,j)*/
	     if(op->info) adjust_prn_level(op);
	     
	     if(i_shuffle<op->n_shuffle)
	       {
		 if(op->info) /* if any level of info */
		   {
		     /* for -i -prnshf scores only for shuffled seqs */
		     if(!op->info2) op->prnalg_lev=1; /* scores only */
		     if(op->prn_shuffle)
		       prn_align(op,F->out1,Info,db_Seq,Alph,query_Seq,&n,
				 &i,&j,&A_ene,aGaps,op->qname,DB_name,
				 &n_row,&n_col); 
		   }
		 e_ave+=A_ene; 
		 sigma+=A_ene*A_ene;
	       }
	     else /* print the original alignment to alignments.log */
	       {
		 ftmp=get_zscore(op,F,Info,R_ene,&i_str,&e_ave,&sigma);
		 if(ftmp>op->t_zscore)
		   prn_align(op,F->balg,Info,db_Seq,Alph,loc_Seq,&n,&i,&j,
			     &A_ene,aGaps,op->qname,DB_name,&n_row,&n_col); 
	       }
	     
	     ++ind;
	     
	   } /* end of shuffling loop */	     

       } /* end if structure among the best matches */

     ++i_str;

   } 
 /* end of alignment loop over all the sequences in the DB */

 /* print best scores */
 pick_best_zscore(op,Info); 
 print_best_matches(op,F,Info); 
 
 free(loc_Seq);
 free(R_ene);
   
} /* end of shuffle_seq_algn() */

/* define environment dependent gap and deletion penalties */
void adjust_gap_pen(k_opt *op, int *Ncont, int *j, float *gap_pen_left,
		    float *gap_pen_top, float *gpen_prev, float *gpen)
{
  float gpen_loc;
	
  /* arbitrary heuristics for seq-seq alignments */
  if(op->align) gpen_loc=(op->gap_pen - (5.0 - *(Ncont + *j - 1)));
  /* use external definition of gpen for threading */
  if(op->exam) gpen_loc=*gpen;
  
  if(*j==1) *gpen_prev=gpen_loc;
  *gap_pen_left=gpen_loc;
  if(op->align) *gap_pen_top=(gpen_loc + *gpen_prev)/2.0;
  if(op->exam && !op->cdel) *gap_pen_top=(gpen_loc + *gpen_prev)/2.0;
  *gpen_prev=gpen_loc;

  /* so, unless chosen otherwise (for threading) deletion penalty  
     is averaged over two nearest insertion penalties */
  
}


/* OLD version !!! Can be useful for two external sequences */
/* align_seq is supposed to align two sequences using dynamic prog */
void align_seq_old(k_opt *op,a_file *F,int *Imod,float *Potn,int *seq,
	       k_info *Info,int *Alph,int *Adr,int *Seq_col,float *Table,
	       int *Trace,int *aGaps)
{
  int i_prot,n,m,i,j,k,n_row,n_col,alpha,beta,i_rewind,k_flag,
      ovrlap_flag,n_align;
  float T_ene,gap_pen,gap_pen_left,gap_pen_top,Tij;
  char T_name[MAXS],R_name[MAXS];

  /* R - reference (row) sequence; T - another (column) seq - 
     use op->beg and op->end to pick them up from F-seq file */ 

 /* skip some sequences if beg>1 */
 i_prot=0;
 while((i_prot<op->beg-1) &&
       read_seq(op,F->seq,&n_row,T_name,Seq_col,&i_prot,&op->n_type,Alph,Adr));
  
 /* get now first seq (T_name) and put it into column Seq_col */
 if(!read_seq(op,F->seq,&n_row,T_name,Seq_col,&i_prot,&op->n_type,Alph,Adr)) 
   read_err(F->f_seq);
 

 if(op->beg==op->end) /* align seq against itself */
   {
     strcpy(R_name,T_name);
     n_col=n_row;
     seq=Seq_col; 
     if(op->debug) prn_seq_alph(seq,&n_col,R_name,Alph);
   }
 else
   {
     /* skip again some sequences */
     while((i_prot<op->end-1) &&
       read_seq(op,F->seq,&n_col,R_name,seq,&i_prot,&op->n_type,Alph,Adr));
 
     /* get now first seq (R_name) and put it into row seq */
     if(!read_seq(op,F->seq,&n_col,R_name,seq,&i_prot,&op->n_type,Alph,Adr)) 
       read_err(F->f_seq);
   }
 
 /* init Dpt and build initial row and col using pre- and suf gap_pen */
 init_dpt(op,Table,Trace,&n_row,&n_col,aGaps);

 /* build the DP table */
 gap_pen_left=op->gap_pen;
 gap_pen_top=op->gap_pen;
 for(i=1;i<n_row;++i) 
   for(j=1;j<n_col;++j) 
     {
       /* seqs start at the second  r/c */
       alpha=*(Seq_col+i-1); 
       beta=*(seq+j-1);
       /* get the proper scoring matrix element */
       if(alpha>beta) Tij=*(Potn+mat2vec(&beta,&alpha,&op->n_type));
       else Tij=*(Potn+mat2vec(&alpha,&beta,&op->n_type)); 
       /* now decide how to proceed having Tij and gap_pen_left/right */
       *(Table+i*n_col+j)=get_cell(op,Table,Trace,&n_col,&i,&j,&gap_pen_left,
				   &gap_pen_top,&Tij);
     }

 /* print DP table */
 if(op->debug) prn_dpt(Table,&n_row,&n_col);

 /* trace back the alignment */
 T_ene=find_best_successor(op,Table,&n_row,&n_col,&i,&j); 
 /* so, now we are at the cell (i,j) having the value T_ene */
 n=trace_back(op,Table,&n_row,&n_col,&i,&j,&T_ene,Trace,aGaps); 
 
 /* print the two aligned seqs - aGaps indicates how to move from (i,j)*/
 prn_align_old(op,F,seq,Alph,Seq_col,&n,&i,&j,&T_ene,aGaps,T_name,R_name,
	   &n_row,&n_col);
 
 /* seq alignment is supposed to be performed separately - finish now */
 exit(0); 
 
}

/* init Dpt and build initial row and col using pre- and suf gap_pen */
void init_dpt(k_opt *op,float *Table,int *Trace,int *n_row,int *n_col,
	      int *aGaps)
{
 int i,j,n_align;
  
 if(op->loc_a) op->pref_pen=0;
 n_align=*n_col+*n_row;
 for(i=0;i<n_align;++i) *(aGaps+i)=0;
 ++*n_row;
 ++*n_col; /* add initial row and column */
 for(i=0;i<*n_row;++i) 
   for(j=0;j<*n_col;++j) 
     {
       *(Table+i*(*n_col)+j)=0.0;
       *(Trace+i*(*n_col)+j)=0;
     }
 for(j=1;j<*n_col;++j) 
   {
     *(Table+j)=*(Table+j-1)+op->pref_pen;
     *(Trace+j)=1; /* so you came from the left in the first row */
   }
 for(i=1;i<*n_row;++i) 
   {
     *(Table+i*(*n_col))=*(Table+(i-1)*(*n_col))+op->pref_pen;
     *(Trace+i*(*n_col))=*n_col;  /* you came from the top in the first col */
   }
}


/* print aligned sequences or structures */
float prn_align(k_opt *op,FILE *Fout,k_info *Info,int *Seq_row,int *Alph,
		int *Seq_col,int *n,int *i,int *j,float *T_ene,int *aGaps,
		char *T_name,char *R_name,int *n_row,int *n_col)
{
 int k,n_l,beg_query,beg_target,i_tens,j_tens,i_prev,j_prev,ii,jj,ident,
     ii_beg,jj_beg;
 char fmt[MAXS],line1[MAXS],line2[MAXS],line3[MAXS],line4[MAXS],
      sym_row[MAXS],sym_col[MAXS],sym_row_res[MAXS],sym_col_res[MAXS],
      ctmp[MAXS],gap_sym[MAXS],gap_row[MAXS],gap_1let[MAXS],ctype[MAXS],
      ver_sym_row[MAXS],ver_sym_col[MAXS];
 float p_ident;
 

 /* both the query sequence (or sequence of aligned structure) 
    and the target sequence (or sequence of aligned structure)
    are passed here as Seq_col (query) and Seq_row (target), respectively */

 /* info level determines op->prnalg_lev option:
    info1/2 (info) level of output is used to print alignments horizontally
    info3 is used to print alignments vertically (no horizontal output then)
    if(!op->info1/2/3) then print general information without alignments */

 /* prepare gap symbols and initiate strings */
 strcpy(fmt,"%s\t%s \n");
 strcpy(gap_sym,"GAP");
 strcpy(gap_row,"GAP");
 strcpy(line1,"");
 strcpy(line2,"");
 strcpy(line3,"");
 strcpy(line4,"");
 strcpy(gap_1let,"-");
 strcpy(ctype,"global(seq-seq)");
 if(op->loc_a) strcpy(ctype,"local(seq-seq)");
 /* if not just seq vs seq alignment alter the symbols */
 if(op->exam) /* seq vs stru */
   {
     sprintf(ctmp,"%d",0);
     strcpy(gap_row,ctmp);
     strcpy(ctype,"global(seq-stru)");
     if(op->loc_a) strcpy(ctype,"local(seq-stru)");
   } 
 if(op->strucmp) /* stru vs stru */
   {
     sprintf(ctmp,"%d",0);
     strcpy(gap_row,ctmp);
     strcpy(gap_sym,ctmp);
     strcpy(ctype,"global(stru-stru)");
     if(op->loc_a) strcpy(ctype,"local(stru-stru)");
   } 
 
 /* print the two aligned seqs - aGaps indicates how to move from (i,j)*/

 /* first print what is supposed to be printed anyway */
 sprintf(Info->Sout,
	 "query= %s\tmatch= %s\tenergy= %6.2f\tlength= %d\ttype= %s \n",
	 T_name,R_name,*T_ene,*n,ctype); 
 to_file(Fout,Info->Sout);
 if(op->prnalg_lev>1)
   {
     sprintf(Info->Sout,"\n");
     to_file(Fout,Info->Sout);
   }

 /* setup all the counting variables */
 beg_query=*i+1;
 beg_target=*j+1;
 i_tens=(beg_query/10)%10+1;
 j_tens=(beg_target/10)%10+1;
 /* assuming beg_ always larger than 0 */
 if(beg_query%10 == 0) --i_tens;
 if(beg_target%10 == 0) --j_tens;
 i_prev=*i;
 j_prev=*j;
 /* not to loose the ends */
 ii_beg=*i;
 jj_beg=*j;
 ii=ii_beg;
 jj=jj_beg;
 ident=0;
 n_l=0;
 for(k=*n-1;k>=0;--k) 
   {

     /* check whether you reached the end of seq/stru */
     if(*i>=*n_row-1) *i=*n_row-2;
     if(*j>=*n_col-1) *j=*n_col-2;
     if(op->debug) 
       {
	 sprintf(Info->Sout," gap counter i=%d aGaps=%d \n",k,*(aGaps+k));
	 to_stnd_out(Info->Sout);
       }

     /* prepare symbols for printing */
     strcpy(sym_row,dig2one(*(Alph+*(Seq_row+*j))));
     strcpy(sym_col,dig2one(*(Alph+*(Seq_col+*i))));
     strcpy(ver_sym_row,dig2amino(*(Alph+*(Seq_row+*j))));
     strcpy(ver_sym_col,dig2amino(*(Alph+*(Seq_col+*i))));
     if(op->strucmp) /* use generic alphabet hoping for good ... */
       {
	 strcpy(sym_row,dig2one(*(Seq_row+*j)));
	 strcpy(sym_col,dig2one(*(Seq_col+*i)));
	 strcpy(ver_sym_row,dig2amino(*(Seq_row+*j)));
	 strcpy(ver_sym_col,dig2amino(*(Seq_col+*i)));	
       }

     /* i_tens and j_tens count how many tens of residues - from 1 to (1)0 */
     if(j_tens%10 == 0) j_tens=0;
     if((*j+1)%10 == 0 && (*j+1)!=j_prev) 
       {
	 sprintf(ctmp,"%d",j_tens);
	 ++j_tens;
	 j_prev=*j+1;
       }
     else sprintf(ctmp,"%s",".");
     strcpy(sym_row_res,ctmp);
     
     if(i_tens%10 == 0) i_tens=0;
     if((*i+1)%10 == 0 && (*i+1)!=i_prev) 
       {
	 sprintf(ctmp,"%d",i_tens);
	 ++i_tens;
	 i_prev=*i+1;
       }
     else sprintf(ctmp,"%s",".");
     strcpy(sym_col_res,ctmp);

     /* modify symbols for vertical printing if target or query is a str */
     if(op->exam) 
       {
	 sprintf(ctmp,"%d",*j+1);
	 strcpy(ver_sym_row,ctmp);
       }
     if(op->strucmp) 
       {
	 sprintf(ctmp,"%d",*j+1);
	 strcpy(ver_sym_row,ctmp);
	 sprintf(ctmp,"%d",*i+1);
	 strcpy(ver_sym_col,ctmp);
       }

     if(*(aGaps+k)==0) /* diagonal */
       {
	 if(op->prnalg_lev == 3)
	   {
	     sprintf(Info->Sout,fmt,ver_sym_col,ver_sym_row);
	     to_file(Fout,Info->Sout);
	   }
	 else if((op->prnalg_lev==2) && n_l<MAXL) 
	   {
	     strcpy(line1,strcat(line1,sym_col_res));
	     strcpy(line2,strcat(line2,sym_col));
	     strcpy(line3,strcat(line3,sym_row));
	     strcpy(line4,strcat(line4,sym_row_res));
	   }
	 ++*i; ++*j; ++ii; ++jj;
	 if(strcmp(sym_col,sym_row)==0) ++ident;
       }
     else 
       {
	 if(*(aGaps+k)==1) /* horizontal */
	   {
	     if(op->prnalg_lev == 3)
	       {
		 sprintf(Info->Sout,fmt,gap_sym,ver_sym_row);
		 to_file(Fout,Info->Sout);
	       }
	     else if((op->prnalg_lev==2) && n_l<MAXL) 
	       {
		 /* strcpy(line1,strcat(line1," ")); */
		 strcpy(line1,strcat(line1,"i"));
		 strcpy(line2,strcat(line2,gap_1let));
		 strcpy(line3,strcat(line3,sym_row));
		 strcpy(line4,strcat(line4,sym_row_res));
	       }
	     ++*j; ++jj;
	   }
	 if(*(aGaps+k)==2) /* vertical */
	   {
	     if(op->prnalg_lev == 3)
	       {
		 sprintf(Info->Sout,fmt,ver_sym_col,gap_row);
		 to_file(Fout,Info->Sout);
	       }
	     else if((op->prnalg_lev==2) && n_l<MAXL) 
	       {
		 /* strcpy(line1,strcat(line1,sym_col_res)); */
		 strcpy(line1,strcat(line1,"d"));
		 strcpy(line2,strcat(line2,sym_col));
		 strcpy(line3,strcat(line3,gap_1let));
		 strcpy(line4,strcat(line4," "));
	       }
	     ++*i; ++ii;
	   }	 
       }
     ++n_l;
     if((op->prnalg_lev != 3) && (op->prnalg_lev==2) && n_l==MAXL-1) 
       {
	 /* add the range of residues for each line */
	 sprintf(ctmp,"%s","   ");
	 strcpy(line1,strcat(line1,ctmp));
	 strcpy(line4,strcat(line4,ctmp));
	 sprintf(ctmp,"%s","   query");
	 strcpy(line2,strcat(line2,ctmp));
	 sprintf(ctmp,"%s","   match");
	 strcpy(line3,strcat(line3,ctmp));
	 if(beg_query>*i) beg_query=*i;
	 if(beg_target>*j) beg_query=*j;
	 sprintf(ctmp,"%d",beg_query);
	 strcpy(line1,strcat(line1,ctmp));
	 sprintf(ctmp,"%d",beg_target);
	 strcpy(line4,strcat(line4,ctmp));
	 sprintf(ctmp,"%s"," - ");
	 strcpy(line1,strcat(line1,ctmp));
	 strcpy(line4,strcat(line4,ctmp));
	 sprintf(ctmp,"%d",*i);
	 strcpy(line1,strcat(line1,ctmp));
	 sprintf(ctmp,"%d",*j);
	 strcpy(line4,strcat(line4,ctmp));
	 beg_query=*i+1;
	 beg_target=*j+1;
	 /* beware of MAXS and MAXL - MAXL > 2 MAXS */
	 sprintf(Info->Sout,"%s\n%s\n",line1,line2);
	 to_file(Fout,Info->Sout);
	 sprintf(Info->Sout,"%s\n%s\n\n",line3,line4);
	 to_file(Fout,Info->Sout);
	 strcpy(line1,"");
	 strcpy(line2,"");
	 strcpy(line3,"");
	 strcpy(line4,"");
	 n_l=0;
       }
   }

 /* calculate the percentage of identity in aligned pairs */
 /* p_ident=ident*100.0/((float)(*n_row-1)); */
 p_ident=ident*100.0/((float)(ii-ii_beg)); /* take only aligned part */
 
 /* print the last piece using one-lett. code */
 if((op->prnalg_lev==2) && (op->prnalg_lev!=3)) 
   {
     /* add the range of residues for each line and identity percentage*/
     sprintf(ctmp,"%s","   ");
     strcpy(line1,strcat(line1,ctmp));
     strcpy(line4,strcat(line4,ctmp));
     sprintf(ctmp,"%s  (%s%5.1f  %s%5.1f)","   query","%align=",
	     (ii-ii_beg)*100.0/((float)(*n_row-1)),"%ident=",p_ident);
     strcpy(line2,strcat(line2,ctmp));
     sprintf(ctmp,"%s  (%s%5.1f)","   match","%align=",
	     (jj-jj_beg)*100.0/((float)(*n_col-1)));
     strcpy(line3,strcat(line3,ctmp));
     if(beg_query>*i) beg_query=*i;
     if(beg_target>*j) beg_query=*j;
     sprintf(ctmp,"%d",beg_query);
     strcpy(line1,strcat(line1,ctmp));
     sprintf(ctmp,"%d",beg_target);
     strcpy(line4,strcat(line4,ctmp));
     sprintf(ctmp,"%s"," - ");
     strcpy(line1,strcat(line1,ctmp));
     strcpy(line4,strcat(line4,ctmp));
     sprintf(ctmp,"%d",ii);
     strcpy(line1,strcat(line1,ctmp));
     sprintf(ctmp,"%d",jj);
     strcpy(line4,strcat(line4,ctmp));
     /* now print the rest */
     sprintf(Info->Sout,"%s\n%s\n",line1,line2);
     to_file(Fout,Info->Sout);
     sprintf(Info->Sout,"%s\n%s\n\n",line3,line4);
     to_file(Fout,Info->Sout);
     sprintf(Info->Sout,"\n%s\n"," ");
     to_file(Fout,Info->Sout);
   }

 return(p_ident);

}


/* print aligned sequences or structures */
void prn_align_old(k_opt *op,a_file *F,int *seq,int *Alph,int *Seq_col,
	       int *n,int *i,int *j,float *T_ene,int *aGaps,char *T_name,
	       char *R_name,int *n_row,int *n_col)
{
 int k,n_l;
 char fmt[MAXS],line1[MAXS],line2[MAXS],sym_row[MAXS],ctmp[MAXS],
      gap_sym[MAXS],gap_row[MAXS],gap_1let[MAXS],ctype[MAXS],sym_col[MAXS];
 
 strcpy(fmt,"%s\t%s \n");
 strcpy(gap_sym,"GAP");
 strcpy(gap_row,"GAP");
 if(op->align) /* seq vs seq, if(!op->info) then print 1let codes */
   {
     strcpy(line1,"");
     strcpy(line2,"");
     strcpy(gap_1let,"-");
     strcpy(ctype,"global(seq-seq)");
     if(op->loc_a) strcpy(ctype,"local(seq-seq)");
   }
 if(op->exam) /* seq vs stru */
   {
     sprintf(ctmp,"%d",0);
     strcpy(gap_row,ctmp);
     strcpy(ctype,"global(seq-stru)");
     if(op->loc_a) strcpy(ctype,"local(seq-stru)");
   } 
 if(op->strucmp) /* stru vs stru */
   {
     sprintf(ctmp,"%d",0);
     strcpy(gap_row,ctmp);
     strcpy(gap_sym,ctmp);
     strcpy(ctype,"global(stru-stru)");
     if(op->loc_a) strcpy(ctype,"local(stru-stru)");
   } 
 
 /* print the two aligned seqs - aGaps indicates how to move from (i,j)*/
 if(op->align)
   {
     printf(" Sequence alignment: \t%s vs %s \n",
	    T_name,R_name); 
     printf(" See file %s \n",F->f_out);
   }
 fprintf(F->out1,"%s\t%s\t length= %d\t score= %6.2f\t type= %s \n",
	 T_name,R_name,*n,*T_ene,ctype); 
 if(op->info && op->align) 
   {
     printf("length= %d\t score= %6.2f\t type= %s \n",*n,*T_ene,ctype); 
     printf("%s\t%s \n",T_name,R_name);
   }
 
 n_l=0;
 for(k=*n-1;k>=0;--k) 
   {
     /* check whether you reached the end of seq/stru */
     if(*i>=*n_row-1) *i=*n_row-2;
     if(*j>=*n_col-1) *j=*n_col-2;
     if(op->debug) printf(" gap counter i=%d aGaps=%d \n",k,*(aGaps+k));
     if(op->align) 
       {
	 strcpy(sym_row,dig2amino(*(Alph+*(seq+*j))));
	 strcpy(sym_col,dig2amino(*(Alph+*(Seq_col+*i))));
       }
     if(op->exam) 
       {
	 sprintf(ctmp,"%d",*j+1);
	 strcpy(sym_row,ctmp);
	 strcpy(sym_col,dig2amino(*(Alph+*(Seq_col+*i))));
       }
     if(op->strucmp) 
       {
	 sprintf(ctmp,"%d",*j+1);
	 strcpy(sym_row,ctmp);
	 sprintf(ctmp,"%d",*i+1);
	 strcpy(sym_col,ctmp);
       }

     if(*(aGaps+k)==0) /* diagonal */
       {
	 if(op->info)
	   {
	     if(op->align) printf(fmt,sym_col,sym_row);
	     fprintf(F->out1,fmt,sym_col,sym_row);
	   }
	 else if(op->align && n_l<MAXL) 
	   {
	     strcpy(line1,strcat(line1,dig2one(*(Alph+*(seq+*j)))));
	     strcpy(line2,strcat(line2,dig2one(*(Alph+*(Seq_col+*i)))));
	   }
	 ++*i; ++*j;
       }
     else 
       {
	 if(*(aGaps+k)==1) /* horizontal */
	   {
	     if(op->info)
	       {
		 if(op->align) printf(fmt,gap_sym,sym_row);
		 fprintf(F->out1,fmt,gap_sym,sym_row);
	       }
	     else if(op->align && n_l<MAXL) 
	       {
		 strcpy(line1,strcat(line1,dig2one(*(Alph+*(seq+*j)))));
		 strcpy(line2,strcat(line2,gap_1let));
	       }
	     ++*j;
	   }
	 if(*(aGaps+k)==2) /* vertical */
	   {
	     if(op->info)
	       {
		 if(op->align) printf(fmt,sym_col,gap_row);
		 fprintf(F->out1,fmt,sym_col,gap_row);
	       }
	     else if(op->align && n_l<MAXL) 
	       {
		 strcpy(line1,strcat(line1,gap_1let));
		 strcpy(line2,strcat(line2,dig2one(*(Alph+*(Seq_col+*i)))));
	       }
	     ++*i; 
	   }	 
       }
     ++n_l;
     if(!op->info && op->align && n_l==MAXL-1) 
       {
	 printf("%s\n%s\n\n",line1,line2);
	 strcpy(line1,"");
	 strcpy(line2,"");
	 n_l=0;
       }
   }
 /* print the last piece using one-lett. code */
 if(op->align && !op->info) printf("%s\n%s\n",line1,line2); 

}


/* prn_dpt prints out the dynamic programming table */
void prn_dpt(float *Table,int *n_row,int *n_col)
{
  int i,j;

  printf("\n DP table:\n");
  for(i=0;i<*n_row;++i)
    {
      for(j=0;j<*n_col;++j) printf(" %3.0f",*(Table+i*(*n_col)+j));
      printf("\n");
    }  
}

/* prn_seq_alph prints seq using current alphabet (not generic one) */
void prn_seq_alph(int *seq,int *Lseq,char *Tname,int *Alph)
{
  int j;
  
  printf(" Sequence %s as defined in the current alphabet: \n",Tname);
  for(j=0;j<*Lseq;++j) printf("%s \n",dig2amino(*(Alph+*(seq+j))));
 
}

int trace_back(k_opt *op,float *Table,int *n_row,int *n_col,
	       int *i,int *j,float *T_ene,int *Trace,int *aGaps)
{
 int k,i_rewind,ovrlap_flag,n,k_flag;
 float Tij;
 
 k=(*i)*(*n_col)+*j; /* i=0, ... n_row-1, j=0, ... n_col-1  */

 /* now use i,j as running indices of seq, Seq_col - they start from the
    second r/c so adjust i,j - they become negative while reaching first r/c*/
 --*i;
 --*j;
 k_flag=1;
 ovrlap_flag=0;
 n=0;
 /* if prefix penalty is zero then do overlap matches */
 if(*(Table+*n_col)==0.0) ovrlap_flag=1; 
 if(op->loc_a && *T_ene==0.0) k_flag=0;
 
 /* start tracing loop */ 
 while(k_flag) 
   {
     i_rewind=*(Trace+k);
     k-=i_rewind;
     Tij=*(Table+k);
     if(op->debug) printf(" value  %3.0f trace %d row %d col %d \n",
			  Tij,i_rewind,*i,*j);

     /* get the flag to stop the loop */
     if(op->loc_a && Tij==0.0) k_flag=0;  /* local alignment */
     else if(ovrlap_flag && (*i<=0 || *j<=0)) k_flag=0;
     else if(*i<=0 && *j<=0) k_flag=0; /* global alignment */

     /* make a step back */
     if(i_rewind<*n_col) /* horizontal */
       {
	 *(aGaps+n)=1;
	 --*j;
       }
     else 
       {
	 if(i_rewind==*n_col)  /* vertical */
	   {
	     *(aGaps+n)=2;
	     --*i; 
	   }
	 else /* diagonal */
	   {
	     --*j; --*i;
	   }	 
       }
     ++n; 
   }

 ++*i; /* get back the starting points in seq and Seq_col */
 ++*j;
 
 return(n); /*number of aligned residues */
 
}

float find_best_successor(k_opt *op,float *Table,int *n_row,int *n_col,
			  int *n,int *m)
{
  float cost,ftmp;
  int i,j;
  
  /* for standard global alignment take right-bottom corner */
  *n=*n_row-1;
  *m=*n_col-1;
  cost=*(Table+(*n_row)*(*n_col)-1);  
 
  /* modify it for best overlap match - check the right and bottom c/r */
  if(*(Table+*n_col)==0.0 && !op->loc_a) 
    {
      for(j=1;j<*n_col;++j) 
	{
	  ftmp=*(Table+(*n_row-1)*(*n_col)+j);
	  if(ftmp<cost)
	    {
	      *n=*n_row-1;
	      *m=j;
	      cost=*(Table+(*n_row-1)*(*n_col)+j);  
	    }
	}
      
      for(i=1;i<*n_row;++i) 
	{
	  ftmp=*(Table+(i+1)*(*n_col)-1);
	  if(ftmp<cost)
	    {
	      *n=i;
	      *m=*n_col-1;
	      cost=*(Table+(i+1)*(*n_col)-1);  
	    }	  
	}    
    }

  /* modify it for best local alignment - check the whole table */
  if(op->loc_a) 
    {
      for(i=1;i<*n_row;++i)
	for(j=1;j<*n_col;++j) 
	  {
	    ftmp=*(Table+(i-1)*(*n_col)+j);
	    if(ftmp<cost)
	      {
		*n=i-1;
		*m=j;
		cost=*(Table+(i-1)*(*n_col)+j);  
	      }
	  }
    }
  
  return(cost);
  
}

float get_cell(k_opt *op,float *Table,int *Trace,int *n_col,int *i,int *j,
	       float *gap_pen_left,float *gap_pen_top,float *Tij)
{
  float value,ftmp;
  
  /* Trace can be equal to one (if we go from the left), n_col (if we go
     from above) or n_col+1 (if we go along diagonal) - one has to rewind
     pointer by this much while going back */

  /* get value if you go along diagonal - if degeneracy go along diag */
  value=*(Table+(*i-1)*(*n_col)+*j-1)+*Tij;
  *(Trace+*i*(*n_col)+*j)=*n_col+1;

  /* get value if you go from the left */
  ftmp=*(Table+*i*(*n_col)+*j-1)+*gap_pen_left;
  if(ftmp<value)
    {
      value=ftmp;
      *(Trace+*i*(*n_col)+*j)=1;
    }

  /* get value if you go from the top */
  ftmp=*(Table+(*i-1)*(*n_col)+*j)+*gap_pen_top;
  if(ftmp<value)
    {
      value=ftmp;
      *(Trace+*i*(*n_col)+*j)=*n_col;
    }

  if(op->loc_a && value>0.0) value=0.0;
  return(value);
  
}

