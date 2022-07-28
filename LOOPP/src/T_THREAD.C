/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000     v. 2.000                                  */
/*                                                                */
/*        DIFFERENT THREADING FUNCTIONS                           */
 
#include "loopp.h"
   
/* thread_seq performs threading - both optimal (gapped) and ungapped */
void thread_seq(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
		k_info *Info, k_Model *model, k_dynpt *dpt,
		k_mps *var_mps[], k_mps *tail_mps[])
{
  int i_prot,n,m,Lstru,Lseq,i,j,k,wind,n_row,n_col,alpha,beta,*IJ_loc;
  int i_best,change,i_rm,i_str,include,*nat_Seq,*query_Seq,store_kbest;
  int pick_stru,lpick,pairwise,tom,ind_short;
  int *contR, *contT, *IJcont, *Imod, *Alph, *Adr, *Trace, *aGaps, *Stru_row;
  int noper,neoper;
  unsigned long int ind;
  float T_ene,gap_pen_left,gap_pen_top,gpen,Tij,E_rm,gpen_prev;
  float *Potn, *R_ene, *contRvdw, *contTvdw, *r12ij, *r6ij, *rij, *Table;
  char T_name[MAXS],Q_name[MAXS];

  time_t time1,time2,time3,time4;
  double sum_dtime,ener_time,se_time;
  float s_pct_env;
  int n_thread,n_dpte;

 /* threading loop - for each structure thread the sequences from f_query */
 /* find an optimal seq-stru alignment when op->gaps (gapless thr. otherwise)*/
 /*      Q - query (Q_name) sequence   */
 /*      T - template or database (T_name) structure  */ 

 /* prepare local shortcuts */
 nat_Seq=prot->native_Seq;
 query_Seq=prot->query_Seq;
 IJcont=prot->IJcont;
 contR=prot->conR;
 contT=prot->conT;
 contRvdw=prot->conRvdw;
 contTvdw=prot->conTvdw;
 Imod=model->Imod;
 Potn=model->Potn;
 R_ene=prot->R_ene;
 Alph=model->Alph;
 Adr=model->Adr;
 Table=dpt->Table;
 Trace=dpt->Trace;
 aGaps=dpt->aGaps;
 r12ij=prot->IJr12;
 r6ij=prot->IJr6;
 rij=prot->native_IJrij;
 Stru_row=dpt->iniRow;

 fea=&fea_data;
 noper=0;
 neoper=0;
 n_thread=0;
 sum_dtime=0.0;
 s_pct_env=0.0;
 se_time=0.0;
 ener_time=0.0;
 n_dpte=0;

 store_kbest=op->kbest;
 pairwise=0;
 tom=0;
 if(op->eij || op->evdw) pairwise=1;
 if(op->ein || op->einm) tom=1;
 if(op->gaps)
   {
     /* print title bar */
     prn_type_of_run(F,op,Info);
   }

 /* begin loop over sequences Q_name from f_query */
 while(read_seq(op,F->query,&Lseq,Q_name,query_Seq,&i_prot,&op->n_type,Alph,Adr))
   {
     i_prot=0;
     ind=1; wind=0; ind_short=1;
     if(op->info2) 
       {
	 sprintf(Info->Sout,"\n Threading ... through (name_length): \n");
	 to_stnd_out(Info->Sout);
       }

     /* prepare tables for kth best alignments and shuffling */
     if(op->gaps)
       {
	 for(i=0;i<op->max_prot;++i) *(Info->Point_to_best+i)=0; 
	 Info->ntop=0;
	 op->kbest=store_kbest;

       }
     
     /* begin loop over structures (T_name) */
     while(read_contMap(F,op,IJcont,&Lstru,T_name,r12ij,r6ij,rij))
       {  
	 pick_stru=1;
	 ++wind; 
	 if(op->info2) report_progress(Info,&wind,T_name,&Lstru);
	 i_str=wind-1; /* so, in fact wind enumerates structures */ 

	 /* read native sequence for the sake of prn_align */
	 read_seq(op,F->seq,&Lstru,T_name,nat_Seq,&i_prot,&op->n_type,Alph,Adr);
      
	 /* if a given structure does not match pick it will be skipped */
	 if(pick_structure(op,T_name,&Lstru,&Lseq))
	   {
	    
	     /* set to zero unknown native contacts */
	     /* R stands here for reference i.e. template seq on tmpl stru */
	     for(i=0;i<op->n_par;++i) *(contR+i)=0;
	     n=0; m=0; /* used to initiate a decoy on the scaffolding */

	     /* calculating and printing energies of alignments */
	     if(op->gaps) /* align seq vs stru using DP */
	       {

		 n_row=Lseq;
		 n_col=Lstru;
		 /* init Dpt, build initial row and col, fix pre/suf gap_pen */
		 init_dpt(op,Table,Trace,&n_row,&n_col,aGaps);
		 /* now n_row=Lseq+1, n_col=Lstru+1 */   

		 /* build the DP table */
		 gap_pen_left=op->gap_pen;
		 gap_pen_top=op->gap_pen; 
		 /* here conT stands for contact types given the fact that
		    query seq is threaded through tmpl stru */
		 if(tom) get_env(op,Imod,contT,IJcont,&n,nat_Seq,&Lstru,&m,
				 ssh->Ncont,ssh,r12ij,r6ij,Stru_row);
		 time1=time(NULL);
		 n_thread++;
		 if(pairwise) get_env_eij(op,fea,Imod,IJcont,nat_Seq,
					  &Lstru,&noper,r12ij,r6ij,&s_pct_env); 
		 /* start loop over DP table */
		 for(i=1;i<n_row;++i) 
		   {
		     /* seqs start at the second  r/c */
		     alpha=i-1; /* BEWARE - now seq goes to col */
		     if(pairwise) alpha=*(query_Seq+alpha);
		     for(j=1;j<n_col;++j) 
		       {
			 beta=j-1; /* stru goes to row */
			 /* get the proper scoring matrix element */
			 if(pairwise)
			   /*Tij=get_score_eij(op,Imod,contT,IJcont,&n,nat_Seq,
					     &Lstru,&alpha,&beta,contTvdw,
					     r12ij,r6ij,Potn,&gpen,&noper);*/ 
			   Tij=get_energy_eij(op,fea,Imod,IJcont,contT,&n,&beta,&alpha,
						 nat_Seq,&Lstru,Potn,contTvdw,r12ij,r6ij,
						 &gpen,&noper,&neoper,&ener_time);    
			 else
			   Tij=get_score(op,Imod,&alpha,query_Seq,&Lseq,&beta,
				       ssh->Ncont,Potn,ssh->Mcont0,ssh->Mcont1,
				       ssh->Mcont2,&gpen,&ener_time);
			 /* redefine gap penalties using Potn(env) */
			 if(op->gap_adr>=0 || pairwise) 
			   adjust_gap_pen(op,ssh->Ncont,&j,&gap_pen_left,
					  &gap_pen_top,&gpen_prev,&gpen);
			 /* decide how to proceed */
			 *(Table+i*n_col+j)=get_cell(op,Table,Trace,&n_col,
						     &i,&j,&gap_pen_left,
						     &gap_pen_top,&Tij);
			 se_time+=ener_time;
			 n_dpte++;
		       }
		   } /* end of loop over DP table */
		 if (pairwise)
		 {
			 free_up_fea(op,fea);
		 }
		 time2=time(NULL);
/*		 printf("\n the time to thread %s into %s took %f sec\n",Q_name,T_name,difftime(time2,time1));
		 sum_dtime+=difftime(time2,time1);
/*		 printf("\n noper= %d",noper);
		 printf("\n nope= %d",neoper);
/*		 if (strcmp(Q_name,"1sis")==0) 
		 {
			 if (strcmp(T_name,"1sis")==0)
			 {
				 for (x=0; x<op->n_par; ++x)
				 {
					 printf("%d",contT[x]);
				 }
			 }
		 } */
	    
		 /* now find and print the optimal alignment */
		 get_alignment(op,F,ssh,prot,Imod,query_Seq,&Lseq,Info,
			       Alph,Adr,Table,Trace,aGaps,&n_row,&n_col,
			       &ind_short,&wind,nat_Seq,T_name,Q_name);
  
	       }
	     else /* gapless alignments */
	       {
		 if(op->debug && Lseq<=Lstru)
		   {
		     sprintf(Info->Sout,"\n seq %s of %d res",Q_name,Lseq);
		     to_stnd_out(Info->Sout);
		   }
		 while(Lseq<=Lstru-n) 
		   {	     
		     get_contType(op,Imod,contT,IJcont,&n,query_Seq,&Lseq,&m,
				  ssh->Ncont,ssh,contTvdw,r12ij,r6ij);
		     T_ene=energy(op,contT,Potn,contTvdw);
		     prn_res(op,F,Potn,contR,contT,&ind,&T_ene,Q_name,T_name,
			     Info,&n,&Lseq,contRvdw,contTvdw,&var_mps[0],
			     &tail_mps[0]);
		     if(op->info) prep_info(Info,&T_ene,&op->max_histo); 
		     ++ind;
		     n+=JUMP;
		   }
	       } /* end gaps/gapless */

	   } /* end if(pick_stru) */
	
	
       } /* end of loop over structures T_name for sequence Q_name  */
         /* end of threading loop over all the structures */

	 printf("\n\n the average time for all threads for seq %s = %f sec\n",Q_name,(sum_dtime/n_thread));
	 printf("\n the average percent of env_eij used for seq %s = %f \n",Q_name,(s_pct_env/n_thread));
/*	 printf("\n the average time for energy calculation for seq %s = %f sec\n",Q_name,(se_time/n_dpte));*/
	 printf("\n n_thread = %d \n",n_thread);

     sprintf(Info->Sout,"\n");
     to_stnd_out(Info->Sout);
     if(op->gaps) Info->index=ind_short-1;
     else Info->index=ind-1;
     if(Info->ntop<op->kbest) op->kbest=Info->ntop;

	 time3=time(NULL);
 
     if(op->gaps) /* now get Z_scores for shuffled sequence */
       {
	 strcpy(op->qname,Q_name);
	 shuffle_seq(op,F,ssh,contR,contT,IJcont,Imod,Potn,R_ene,
		     nat_Seq,query_Seq,&Lseq,Info,contRvdw,contTvdw,
		     r12ij,r6ij,rij,Alph,Adr,Table,Trace,aGaps,Stru_row);  
       }

	 time4=time(NULL);
	 printf("\n the time it took to shuffle = %f sec \n",difftime(time4,time3));

     /* prepare files to read them again */
     rewind(F->cont); 
     rewind(F->seq);

   }  /* end of loop over all the sequences in the f_query file */

}
/* end of check_seq() */


/* shuffle_seq performs threading for best matches only including shuffling */
void shuffle_seq(k_opt *op, a_file *F, solv_sh *ssh, int *contR, int *contT,
		 int *IJcont, int *Imod, float *Potn, float *R_ene, 
		 int *nat_Seq, int *query_Seq, int *l_qseq, k_info *Info, 
		 float *contRvdw, float *contTvdw, float *r12ij, float *r6ij, 
		 float *rij, int *Alph, int *Adr, float *Table, int *Trace, 
		 int *aGaps, int *Stru_row)
{
  int i_prot,n,m,Lstru,Lseq,i,j,k,wind,n_row,n_col,alpha,beta;
  int i_best,change,i_rm,i_str,i_shuffle,l,lpick,*loc_Seq,pairwise,tom;
  unsigned long int ind;
  float T_ene,gap_pen_left,gap_pen_top,gpen,Tij,E_rm,e_ave,sigma,
        ftmp,gpen_prev;
  char Tmpl_name[MAXS];
  int noper=0,neoper=0;
  float s_pct_env=0.0;
 /* shuffling loop -  thread the sequence through plausible structures */
 
 pairwise=0;
 tom=0;
 if(op->eij || op->evdw) pairwise=1;
 if(op->ein || op->einm) tom=1;

 ind=1; wind=0; i_str=0; i_prot=0; 
 if(op->info2)
   {
     sprintf(Info->Sout,"\n Shuffling ...  \n");
     to_stnd_out(Info->Sout);
   }
 Lseq=*l_qseq;
 /* allocate array storing copy of the query_Seq */
 loc_Seq=malloc((Lseq+1)*sizeof(int));
 if(loc_Seq==NULL) alloc_err("local_copy_of_theSequence");  
 for(i=0;i<Lseq;++i) *(loc_Seq+i)=*(query_Seq+i);
 
 
 /* begin loop over structures (Tmpl_name) */
 rewind(F->cont); rewind(F->seq); 
 while(read_contMap(F,op,IJcont,&Lstru,Tmpl_name,r12ij,r6ij,rij))
   {   
     /* read the native sequence of the template Tmpl_name */
     read_seq(op,F->seq,&Lstru,Tmpl_name,nat_Seq,&i_prot,&op->n_type,Alph,Adr);
     /* assuming that 1-1 correspond. between SEQ and XYZ already checked */ 

     if(op->pick) /* skip all except picked if desired */
       {
	 lpick=strlen(op->pstr_name);
	 if(lpick>0 && strncmp(Tmpl_name,op->pstr_name,lpick)!=0)
	   *(Info->Point_to_best+i_str)=0;
       }
     
     /* do the shuffling for good templates only */
     if(*(Info->Point_to_best+i_str))
       {
	 if(op->info2) 
	   {
	     sprintf(Info->Sout,"  %s_%d  ",Tmpl_name,Lstru);
	     to_stnd_out(Info->Sout);
	   }
	 /* set to zero unknown native contacts */
	 for(i=0;i<op->n_par;++i) *(contR+i)=0;
	 /* calculate env by fea for pairwise potentials */
	 if(pairwise) get_env_eij(op,fea,Imod,IJcont,nat_Seq,
					  &Lstru,&noper,r12ij,r6ij,&s_pct_env);
	 /* shuffle original sequence */
	 for(i=0;i<20;++i) reshuffle_seq(query_Seq,&Lseq); 
	 e_ave=0.0;
	 sigma=0.0;
	 /* now one more time to recreate the original alignment */
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
	     if(tom) get_env(op,Imod,contT,IJcont,&n,nat_Seq,&Lstru,&m,
			     ssh->Ncont,ssh,r12ij,r6ij,Stru_row);
	     for(i=1;i<n_row;++i) 
	       {
		 /* seqs start at the second  r/c */
		 alpha=i-1; /* BEWARE - now seq goes to col */
		 if(pairwise) 
		   {
		     if(i_shuffle<op->n_shuffle) alpha=*(query_Seq+alpha);
		     else alpha=*(loc_Seq+alpha);
		   }
		 for(j=1;j<n_col;++j) 
		   {
		     beta=j-1; /* stru goes to row */
		     /* get the proper scoring matrix element */
		     if(tom)
		       {
			 if(i_shuffle<op->n_shuffle)
			   Tij=get_score(op,Imod,&alpha,query_Seq,&Lseq,&beta,
				     ssh->Ncont,Potn,ssh->Mcont0,ssh->Mcont1,
				     ssh->Mcont2,&gpen);
			 else /* original seq to print original alignment */
			   Tij=get_score(op,Imod,&alpha,loc_Seq,&Lseq,&beta,
				     ssh->Ncont,Potn,ssh->Mcont0,ssh->Mcont1,
				     ssh->Mcont2,&gpen);
		       }
		     if(pairwise)
		       /*Tij=get_score_eij(op,Imod,contT,IJcont,&n,nat_Seq,
					 &Lstru,&alpha,&beta,contTvdw,
					 r12ij,r6ij,Potn,&gpen,&noper);*/
			   Tij=get_energy_eij(op,fea,Imod,IJcont,contT,&n,&beta,&alpha,
						 nat_Seq,&Lstru,Potn,contTvdw,r12ij,r6ij,
						 &gpen,&noper,&neoper);
		     /* redefine gap penalties using Potn(env) */
		     if(op->gap_adr>=0 || pairwise) 
			   adjust_gap_pen(op,ssh->Ncont,&j,&gap_pen_left,
					  &gap_pen_top,&gpen_prev,&gpen);
		     /* decide how to proceed */
		     *(Table+i*n_col+j)=get_cell(op,Table,Trace,&n_col,
				      &i,&j,&gap_pen_left,&gap_pen_top,&Tij);
		   }
	       }
	    
	     /* print DP table */
	     if(op->debug) prn_dpt(Table,&n_row,&n_col);
	     /* trace back the alignment */
	     T_ene=find_best_successor(op,Table,&n_row,&n_col,&i,&j); 
	     /* so, we are at the cell (i,j) having the value T_ene */
	     n=trace_back(op,Table,&n_row,&n_col,&i,&j,&T_ene,Trace,aGaps); 

	     /* print alignment - aGaps defines how to move from (i,j)*/
	     if(op->info) adjust_prn_level(op);
	     if(i_shuffle<op->n_shuffle)
	       {
		 if(op->info) /* if any level of info */
		   {
		     /* for -i -prnshf scores only for shuffled seqs */
		     if(!op->info2) op->prnalg_lev=1; /* scores only */
		     if(op->prn_shuffle)
		       prn_align(op,F->out1,Info,nat_Seq,Alph,query_Seq,&n,
				 &i,&j,&T_ene,aGaps,op->qname,Tmpl_name,
				 &n_row,&n_col);
		   }
		 e_ave+=T_ene; 
		 sigma+=T_ene*T_ene;
	       }
	     else /* print the original alignment to alignments.log */
	       {
		 ftmp=get_zscore(op,F,Info,R_ene,&i_str,&e_ave,&sigma);
		 if(ftmp>op->t_zscore)
		   prn_align(op,F->balg,Info,nat_Seq,Alph,loc_Seq,&n,&i,&j,
			     &T_ene,aGaps,op->qname,Tmpl_name,&n_row,&n_col); 
	       }
		 
	     /* if(op->info) prep_info(Info,&T_ene,&op->max_histo); */
	     ++ind;
	     
	   } /* end of shuffling loop */
	 if (pairwise)
	 {
		 free_up_fea(op,fea);
	 }

       } /* end if structure among the best matches */

     ++i_str;

   } 
 /* end of threading loop over all the structures */

 /* print best scores */
 if(op->gaps)
   {
     pick_best_zscore(op,Info); 
     print_best_matches(op,F,Info); 
   }
 
 free(loc_Seq);
 
   
} /* end of shuffle_seq() */

/* print names of templates as the loop proceeds */
void report_progress(k_info *Info, int *wind, char *DB_name, int *l_stru)
{
  sprintf(Info->Sout,"  %s_%d",DB_name,*l_stru);
  to_stnd_out(Info->Sout);
  if((*wind%10)== 0) 
    {
      sprintf(Info->Sout,"\n");
      to_stnd_out(Info->Sout);
    }
}

/* if a given structure does not match pick it will be skipped */
int pick_structure(k_opt *op, char *DB_name, int *db_len, int *q_len)
{
  int lpick,pick_stru;
  float len_ratio,low_ratio,up_ratio;
  
  
  pick_stru=1; /* accept by default */
  
  if(op->pick) /* skip all except picked if desired */
    {
      lpick=strlen(op->pstr_name);
      if(lpick>0 && strncmp(DB_name,op->pstr_name,lpick)!=0)
	pick_stru=0;
      /* here one may extend the pick conditions */
    }

  /* use lengths to narrow the search: -sd 3 means no restrictions here */
  /* matching of subdomains best with no restrictions */
  len_ratio=(float)(*q_len)/(float)(*db_len);

  /* search depth level 2: to be added */

  /* level 1: relatively mild restrictions  */
  if(op->search_depth==1) 
    {
      if(op->loc_a) 
	{
	  low_ratio=L_RATIO_2; /* query may be up to 4 times shorter */
	  up_ratio=1.0/L_RATIO_1; /* and up to 2 times longer than db */
	}
      else /* now global alignemnts */
	{
	  low_ratio=L_RATIO_1;
	  up_ratio=1.0/L_RATIO_1;;
	}
      if((len_ratio<low_ratio) || (len_ratio>up_ratio)) 
	{
	  pick_stru=0;
	  /*
	  printf(" Skipping %s of len %d   len_rat=%f low_rat=%f up_rat=%f \n",
		 DB_name,*db_len,len_ratio,low_ratio,up_ratio);
		 */
	}
      
    }

  /* level 0: narrower search - local essentially only as a filter */
  if(op->search_depth==0) 
    {
      if(op->loc_a) 
	{
	  low_ratio=L_RATIO_0; /* q more than 0.75 of db */
	  up_ratio=1.0+2*L_RATIO_3; /* q up to 1.2 longer than db */
	}
      else /* now for global */
	{
	  low_ratio=L_RATIO_0;
	  up_ratio=1.0+L_RATIO_3; /* q up to 1.1 longer than db */
	}
      if((len_ratio<low_ratio) || (len_ratio>up_ratio)) 
	{
	  pick_stru=0;
	  /* 
	  printf(" Skipping %s of len %d   len_rat=%f low_rat=%f up_rat=%f \n",
		 DB_name,*db_len,len_ratio,low_ratio,up_ratio);
		 */
	}
      
    }
  
     
  return(pick_stru);
  
}

float get_zscore(k_opt *op,a_file *F,k_info *Info,float *R_ene,
		int *i_str,float *e_ave,float *sigma)
{
  int i,j,shift_z2;
  float T_ene,z_score,ftmp;

  /* calculate now average and sigma for shuffled sequnce */
  *e_ave=*e_ave/(float)op->n_shuffle;
  z_score=*sigma/(float)op->n_shuffle;
  ftmp=(*e_ave)*(*e_ave);
  ftmp=z_score-ftmp;
  if(ftmp<epsilon) ftmp=epsilon;
  z_score=1.0/sqrt(ftmp);
  if(op->info2) 
    {
      sprintf(Info->Sout," e_ave= %5.2f   st_dev= %5.2f  ",
	      *e_ave,(1.0/z_score));
      to_stnd_out(Info->Sout);
    }
  /* calculate now the actual z_score and send them via Info */
  for(i=0;i<op->kbest;++i)
    {
      if(*(Info->Best_str+i)==*i_str) 
	{
	  T_ene=*(Info->Best_ene+i);
	  ftmp=-(*(R_ene+*i_str)-T_ene)*z_score;
	  z_score=-(T_ene-*e_ave)*z_score;
	  *(Info->Zscore1+i)=z_score;	    
	  *(Info->Zscore2+i)=ftmp;
	  if(op->info2) 
	    {
	      sprintf(Info->Sout,
		      " e_query= %5.2f   e_nat= %5.2f \n",
		      T_ene,*(R_ene+*i_str));
	      to_stnd_out(Info->Sout);
	    }
	}
    }   
  return(z_score);
  
}

/* to print the list of best structures according to Z_score */
void print_best_matches(k_opt *op, a_file *F, k_info *Info)
{
  int i,j,shift_z2,hits;
  float z1,z2,th_low,th_high_low,th_high,th_very_high;
  char sig[MAXS],match[MAXS],query[MAXS];
  
  
  strcpy(match,"structure");
  strcpy(query,"sequence");
  if(op->align) strcpy(match,"sequence");
  if(op->strucmp) strcpy(query,"structure");
  sprintf(Info->Sout,
	  "\n The best matching %ss (for query %s %s) are: \n \n",
	  match,query,op->qname);
  to_stnd_out(Info->Sout);
  to_file(F->best,Info->Sout);
  hits=0;
  
  /* differentiate meaning of the second Z-score */
  if(op->loc_a) shift_z2=5;
  else shift_z2=1;

  /* define significance thresholds */
  th_low=2.0;
  th_high_low=3.0;
  th_high=4.0;
  th_very_high=5.0;
  if(op->align)
    {
      th_high=5.0;
      th_very_high=8.0;      
    }
  
  for(i=0;i<op->kbest;++i)
    {
      j=*(Info->toBest+i);
      z1=*(Info->Zscore1+j);
      z2=*(Info->Zscore2+j);
      if(z1>op->t_zscore)
	{
	  ++hits;
	  strcpy(sig,"very_low");
	  if(z1>th_low) strcpy(sig,"low");
	  if(z1>th_high_low) 
	    {
	      strcpy(sig,"HIGH/low");
	      /* scale down the significance on the edge using zscore2 */
	      if(z2<-0.01 && fabs(z2)>(shift_z2*z1)) strcpy(sig,"low");
	    }
	  if(z1>th_high) strcpy(sig,"HIGH");
	  if(z1>th_very_high) strcpy(sig,"VERY_HIGH");
	  /* scale down the significance if energy positive */
	  if(*(Info->Best_ene+j)>0.0) strcpy(sig,"very_low");
	  sprintf(Info->Sout,
	      " %s\t ene= %7.1f\t z_score= %5.1f  (%5.1f)  significance= %s\n",
	     (Info->Names+*(Info->Best_str+j))->mstr,
	     *(Info->Best_ene+j),z1,z2,sig);
	  to_stnd_out(Info->Sout);
	  to_file(F->best,Info->Sout);
	  /* print only k_prn top ones */
	  if(hits>=op->k_prn) break;
	}
    } 

  if(hits==0)
    {
      sprintf(Info->Sout,"\t No significant matches detected \n");
      to_stnd_out(Info->Sout);
      to_file(F->best,Info->Sout);
    }
  
  
}
  

void reshuffle_seq(int *seq,int *Lseq)
{
  int i,k,l,itmp1,itmp2,itmp3,mult;

  k=0; l=*Lseq-1; 
  mult=1;
  
  /* perform some random pair permutations */
  shuffle_rand(seq,Lseq,&mult); 
  
  /* shuffle (some) tenths */
  k=l/10;
  for(i=0;i<k;++i)
    {
      itmp1=*(seq+3*k+i);
      *(seq+3*k+i)=*(seq+6*k-i);
      itmp2=*(seq+4*k-i);
      *(seq+4*k-i)=*(seq+8*k-i);
      *(seq+6*k-i)=itmp1;
      *(seq+8*k-i)=itmp2;
    }

  /* perform some random pair permutations */
  shuffle_rand(seq,Lseq,&mult); 

  /* shuffle halves */
  for(i=0;i<l/2;++i)
    {
      itmp1=*(seq+i);
      *(seq+i)=*(seq+l-i);
      *(seq+l-i)=itmp1;
    }

  /* perform some random pair permutations */
  shuffle_rand(seq,Lseq,&mult); 

  /* shuffle (some) tenths */
  k=l/10;
  for(i=0;i<k;++i)
    {
      itmp1=*(seq+2*k+i);
      *(seq+2*k+i)=*(seq+6*k-i);
      itmp2=*(seq+4*k-i);
      *(seq+4*k-i)=*(seq+8*k-i);
      *(seq+6*k-i)=itmp1;
      *(seq+8*k-i)=itmp2;
    }

  /* perform some random pair permutations */
  shuffle_rand(seq,Lseq,&mult);

  /* shuffle thirds */
  k=l/3;
  for(i=0;i<k;++i)
    {
      itmp1=*(seq+i);
      *(seq+i)=*(seq+2*k-i);
      *(seq+2*k-i)=*(seq+3*k-i);
      *(seq+3*k-i)=itmp1;
    }

  /* perform some random pair permutations */
  shuffle_rand(seq,Lseq,&mult);

  /* shuffle (some) tenths */
  k=l/10;
  for(i=0;i<k;++i)
    {
      itmp1=*(seq+k+i);
      *(seq+k+i)=*(seq+9*k-i);
      itmp2=*(seq+5*k-i);
      *(seq+5*k-i)=*(seq+7*k-i);
      *(seq+7*k-i)=itmp2;
      *(seq+9*k-i)=itmp1;
    }

  /* perform some random pair permutations */
  shuffle_rand(seq,Lseq,&mult);

  /* shuffle fourths 1-2, 3-4 */
  k=l/4;
  for(i=0;i<k;++i)
    {
      itmp1=*(seq+i);
      *(seq+i)=*(seq+2*k-i);
      itmp2=*(seq+3*k-i);
      *(seq+3*k-i)=*(seq+4*k-i);
      *(seq+2*k-i)=itmp1;
      *(seq+4*k-i)=itmp2;
    } 

  /* perform some random pair permutations */
  shuffle_rand(seq,Lseq,&mult);

  /* shuffle (some) tenths */
  k=l/10;
  for(i=0;i<k;++i)
    {
      itmp1=*(seq+k+i);
      *(seq+k+i)=*(seq+7*k-i);
      itmp2=*(seq+3*k-i);
      *(seq+3*k-i)=*(seq+5*k-i);
      *(seq+5*k-i)=itmp2;
      *(seq+7*k-i)=itmp1;
    }
  
}

void shuffle_rand(int *seq,int *Lseq,int *mult)
{
  int i,l,i1,i2,itmp1;
  float r;
  
  
  l=*Lseq-1;
  /* shuffle randomly */
  for(i=0;i<*mult*l/2;++i)
    {  
      r=rand()/(RAND_MAX+1.0);
      i1=(int)(l*r);
      r=rand()/(RAND_MAX+1.0);
      i2=(int)(l*r);
      itmp1=*(seq+i1);
      *(seq+i1)=*(seq+i2);
      *(seq+i2)=itmp1;
      /* printf(" %d - %d ",i1,i2); */
    }
 /* printf(" RMAX= %d \n ",RAND_MAX); */ 
}

void pick_best_aligns(k_opt *op,k_info *Info,int *no_str,float *T_ene)
{
  int i_rm,i_str,i;
  float E_rm;
  
  i_str=*no_str-1;
  
  if(Info->ntop<op->kbest)
    {
      *(Info->Best_ene+Info->ntop)=*T_ene;
      *(Info->Best_str+Info->ntop)=i_str;
      *(Info->Point_to_best+i_str)=1;
      ++Info->ntop;
    }
  else
    {
      /* find the worst among the current best */
      i_rm=0;
      E_rm=*(Info->Best_ene+i_rm);
      for(i=1;i<op->kbest;++i)
	{
	  if(*(Info->Best_ene+i)>E_rm)
	    {
	      i_rm=i;
	      E_rm=*(Info->Best_ene+i);
	    }
	}
      /* E_rm could be saved and used when calling again */
      
      /* check whether current score should replace the worst */
      if(*T_ene<E_rm)
	{

	  *(Info->Point_to_best+*(Info->Best_str+i_rm))=0;
	  *(Info->Best_ene+i_rm)=*T_ene;
	  *(Info->Best_str+i_rm)=i_str;
	  *(Info->Point_to_best+i_str)=1; 

	}
    }
	    
}

/* prepare additional pointer to order Info.Best according to Zscore */
void pick_best_zscore(k_opt *op,k_info *Info)
{
  int i_rm,i_str,i,i_rank,j;
  float E_rm;
  
  /* for stru2stru alignments Zscore1 contains RMSDs */
  for(i=0;i<op->kbest;++i)
    {
      i_rank=0;
      for(j=0;j<op->kbest;++j)
	{
	  if(op->strucmp)
	    {
	      if(*(Info->Zscore1+i)>*(Info->Zscore1+j)) ++i_rank;
	    }
	  else
	    {
	      if(*(Info->Zscore1+i)<*(Info->Zscore1+j)) ++i_rank;
	    }
	}
      *(Info->toBest+i_rank)=i; /* point to i which is i_rank-th */
      
    }
	    
}

/* OLD VERSION - obsolate now */
/* comp_stru performs structure vs structure alignments */
void comp_stru(k_opt *op,a_file *F,solv_sh *ssh,int *contR,int *contT,
	       int *IJcont,int *Imod,float *Potn,float *R_ene,int *seq,
	       k_info *Info,
	       float *r12ij,float *r6ij,float *rij,
	       int *Alph,int *Adr,float *Table,
	       int *Trace,int *aGaps,int *Stru_row,int *Stru_col,
	       float *C_score,float *R_score)
{
  int i_prot,n,m,i,j,k,n_row,n_col,alpha,beta;
  unsigned long int ind;
  float T_ene,gap_pen_left,gap_pen_top,gpen,Tij,x;
  char C_name[MAXS],R_name[MAXS];

 /* perform stru-stru alignment */
      /* C - column structure
         R - row structure */ 

 if(op->eij || op->evdw)
   {
     printf(" This option makes sense only with TOM-like pot ... \n");
     exit(0);
   }
 
 printf(" Performing structure vs structure alignment ... \n");
 n=0; m=0; /* used to initiate a decoy on the scaffolding */
 alpha=0;
 
 /* skip some structures if beg>1 */
 i_prot=0; 
 while((i_prot<op->beg-1) && 
       read_contMap(F,op,IJcont,&n_row,C_name,r12ij,r6ij,rij))
   {
     if(op->debug) printf("  %s_%d \n ",C_name,n_row);
     ++i_prot;
   }
  
 /* get now first stru (C_name) and put it into column Stru_col */
 if(!read_contMap(F,op,IJcont,&n_row,C_name,r12ij,r6ij,rij)) 
   read_err(F->f_cont);
 ++i_prot;

 /* now get the environmnet and scores of sites in C_score */
 get_env(op,Imod,contT,IJcont,&n,seq,&n_row,&m,ssh->Ncont,ssh,
	 r12ij,r6ij,Stru_col); 
 
 for(i=1;i<n_row;++i) 
   {
     beta=i-1;
     *(C_score+beta)=get_score(op,Imod,&alpha,seq,&n_row,&beta,
			       ssh->Ncont,Potn,ssh->Mcont0,ssh->Mcont1,
			       ssh->Mcont2,&gpen);
   }

 if(op->beg==op->end) /* align stru against itself */
   {
     strcpy(R_name,C_name);
     n_col=n_row;
     Stru_row=Stru_col; 
     R_score=C_score;
   }
 else
   {
     /* skip again some structures */
     while((i_prot<op->end-1) &&
	   read_contMap(F,op,IJcont,&n_col,R_name,r12ij,r6ij,rij))  
       {
	 if(op->debug) printf("  %s_%d \n ",R_name,n_col);
	 ++i_prot;
       }    
        
     /* get now second stru (R_name) and put it into row  */
     if(!read_contMap(F,op,IJcont,&n_col,R_name,r12ij,r6ij,rij)) 
       read_err(F->f_cont);

     /* now get the environmnet and scores of sites in R_score */
     get_env(op,Imod,contT,IJcont,&n,seq,&n_col,&m,ssh->Ncont,ssh,
	     r12ij,r6ij,Stru_row);
      
     for(j=1;j<n_col;++j) 
       {
	 beta=j-1;
	 *(R_score+beta)=get_score(op,Imod,&alpha,seq,&n_col,&beta,
				   ssh->Ncont,Potn,ssh->Mcont0,ssh->Mcont1,
				   ssh->Mcont2,&gpen);
       }

   }
 
 /* print some information */
 printf("  %s_%d vs %s_%d \n",C_name,n_row,R_name,n_col);
 printf(" See file %s \n",F->f_out);
 
 /* init Dpt, build initial row and col, fix pre, suf gap_pen */
 init_dpt(op,Table,Trace,&n_row,&n_col,aGaps);
 /* now n_row=Lstru1+1, n_col=Lstru2+1 */   

 /* build the DP table */
 gap_pen_left=op->gap_pen;
 gap_pen_top=op->gap_pen; 
 for(i=1;i<n_row;++i) /* stru start at the second  r/c */
   {  
     alpha=i-1; 		
     for(j=1;j<n_col;++j) 
       {
	 beta=j-1; 
	 /* get the proper scoring matrix element */
	 x=*(C_score+alpha)-*(R_score+beta);
	 /* arbitrary heuristics - top of the gaussian rewarded */
	 Tij=10.0*(0.8-exp(-1.0*x*x));
	  /* redefine gap penalty using Potn(env) */
	 if(op->gap_adr>=0) gap_pen_left=gpen;
	 /* decide how to proceed */
	 *(Table+i*n_col+j)=get_cell(op,Table,Trace,&n_col,&i,&j,
				     &gap_pen_left,&gap_pen_top,&Tij);
       }
   }
	    
 /* print DP table */
 if(op->debug) prn_dpt(Table,&n_row,&n_col);
 
 /* trace back the alignment */
 T_ene=find_best_successor(op,Table,&n_row,&n_col,&i,&j); 
 /* so, now we are at the cell (i,j) having the value T_ene */
 n=trace_back(op,Table,&n_row,&n_col,&i,&j,&T_ene,Trace,aGaps); 
 
 /* print alignment - aGaps indicates how to move from (i,j)*/
 prn_align_old(op,F,Stru_row,Alph,Stru_col,&n,&i,&j,&T_ene,aGaps,
	   C_name,R_name,&n_row,&n_col); 
   
}

/* END OF OLD */

/* get_ineq performs main threading loop and generates ineq for LP training */
void get_ineq(k_opt *op,a_file *F,solv_sh *ssh,int *contR,int *contT,
	      int *IJcont,int *Imod,float *Potn,float *R_ene,int *seq,
	      k_info *Info,float *contRvdw,float *contTvdw,
	      float *r12ij,float *r6ij,float *rij,int *Alph,int *Adr,
	      k_mps *var_mps[],k_mps *tail_mps[])
{
  int i_prot,n,m,Lstru,Lseq,i,j,k,wind,i_database;
  unsigned long int ind;
  float T_ene;
  char T_name[MAXS],R_name[MAXS];

 /* threading loop - for each structure thread all the sequences */
      /* R - reference (native) structure
         T - thread (decoy) structure or the one used to generate decoy */ 

 /* skip some sequences if beg>1 */
 i_prot=0;
 while(i_prot<op->beg-1) skip(F,op,seq,&i_prot,contR,contRvdw,Alph,Adr);

 ind=1; wind=0; i_database=0;
 printf("\n Gapless threading ... through (name_length): \n");
 /* begin loop over structures (T_name) */
 while(read_contMap(F,op,IJcont,&Lstru,T_name,r12ij,r6ij,rij))
   {    
    ++i_database;
    if(op->info) printf("  %d:%s_%d",i_database,T_name,Lstru);
    else printf("  %s_%d",T_name,Lstru);
    ++wind; 
    if(op->info) 
      {
	if((wind%8)== 0) printf("\n");
      }
    else if((wind%10)== 0) printf("\n");
    /* begin loop over sequences for structure T_name */
    while(read_seq(op,F->seq,&Lseq,R_name,seq,&i_prot,&op->n_type,Alph,Adr) 
	  && (i_prot<=op->end))
      {
       /* get structure of R_name as represented in F->type */
       if(!read_contType(F,op,contR,R_name,contRvdw)) read_err(R_name);
       if(strcmp(T_name,R_name)!=0)
         {
	  if(op->debug && Lseq<=Lstru) 
	    printf("\n seq %s of %d res",R_name,Lseq); 
          n=0; m=0; /* used to initiate a decoy on the scaffolding */
          /* calculating and printing energies of decoys */
          while(Lseq<=Lstru-n)
	    {	     
	      get_contType(op,Imod,contT,IJcont,&n,seq,&Lseq,&m,
			   ssh->Ncont,ssh,contTvdw,r12ij,r6ij);
	      T_ene=energy(op,contT,Potn,contTvdw)-*(R_ene+i_prot-1);
	      prn_res(op,F,Potn,contR,contT,&ind,&T_ene,R_name,T_name,Info,
		      &n,&Lseq,contRvdw,contTvdw,&var_mps[0],&tail_mps[0]);
	      if(op->info) 
		{
		  prep_info(Info,&T_ene,&op->max_histo); 
		  add_cont_decoy(&op->n_par,contT,&Info->ncont_decoys,
				 Info->con_All_decoys);
		}
	      ++ind;
	      n+=JUMP;
	    }
	 }
      } /* end of loop over sequences (R_name) for a structure T_name */

    rewind(F->seq); rewind(F->type); i_prot=0;
    /* skip again some sequences if begin>1 */
    while(i_prot<op->beg-1) skip(F,op,seq,&i_prot,contR,contRvdw,Alph,Adr);

   } 
 /* end of threading loop over all the structures */
 printf("\n");  
 Info->index=ind-1;
 
}


/* get_ineq_dec generates ineq for LP training for a set of explicit decoys */
void get_ineq_dec(k_opt *op,a_file *F,solv_sh *ssh,int *contR,int *contT,
		  int *IJcont,int *Imod,float *Potn,float *R_ene,int *seq,
		  k_info *Info,float *contRvdw,float *contTvdw,
		  int *Alph,int *Adr,k_mps *var_mps[],k_mps *tail_mps[])
{
  int i_prot,n,m,Lstru,Lseq,i,j,k,wind;
  unsigned long int ind;
  float T_ene;
  char T_name[MAXS],R_name[MAXS];

 /* this function generates differences in contacts for explicitly 
    given decoys */

  /* R - reference (native) structure - should be first one in XYZ file
     T - decoy structures - all the remaining stru from XYZ  */ 

 /* fix native */
 i_prot=0;  
 /* get structure of R_name as represented in F->type */
 /* to get Lseq */
 read_seq(op,F->seq,&Lseq,R_name,seq,&i_prot,&op->n_type,Alph,Adr); 
 if(!read_contType(F,op,contR,R_name,contRvdw)) read_err(R_name);
 printf("\n Assuming that %s is the native structure \n",R_name);

 ind=1; wind=0; n=0;
 printf("\n Differences in contacts ... for decoy : \n");
 /* begin loop over structures (T_name) - as represented in F->type */
 while(read_seq(op,F->seq,&Lseq,T_name,seq,&i_prot,&op->n_type,Alph,Adr))
   {   
    read_contType(F,op,contT,T_name,contTvdw);
    printf("  %s ",T_name);
    ++wind; if((wind%10)== 0) printf("\n");
    if(strcmp(T_name,R_name)!=0)
      {
	T_ene=energy(op,contT,Potn,contTvdw)-*R_ene; /* - R_ene[0] */
	prn_res(op,F,Potn,contR,contT,&ind,&T_ene,R_name,T_name,Info,
		&n,&Lseq,contRvdw,contTvdw,&var_mps[0],&tail_mps[0]);
	if(op->info) prep_info(Info,&T_ene,&op->max_histo); 
	++ind;
      }
   } 
 /* end of threading loop over all the structures */
 printf("\n");  
 Info->index=ind-1;
 
}

