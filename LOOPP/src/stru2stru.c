/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000     v. 2.000                                  */
/*                                                                */
/*        STRUCTURE TO STRUCTURE ALIGNMENT                        */
 
#include "loopp.h"
   



/* str2str performs structure to structure alignments */
int stru2stru(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
	      k_info *Info, k_Model *model, k_dynpt *dpt)
{
  int i_prot,n,m,i,j,k,n_row,n_col,alpha,beta,l_qstru,lpick,pick_stru,wind;
  unsigned long int ind;
  float T_ene,gap_pen_left,gap_pen_top,gpen,Tij,x,y,x1,x2,perc_empty;
  char C_name[MAXS],R_name[MAXS];
  /* local pointers for shortcuts */
  int *contR,*contT,*IJcont,*Imod,*Alph,*Adr,*Trace,*aGaps,
      *Stru_row,*Stru_col,*C_seq,*R_seq,include,store_kbest;
  float *Potn,*R_ene,*r12ij,*r6ij,*C_rij,*R_rij,*C_score,*R_score,*Table,
        *X_q,*Y_q,*Z_q;

  int lstru,r,s,l_align; 
  float *C_profile,*R_profile,*disvec_CR,scale;
  double norm_dv,norm_d0;
  char name[MAXS]; 
  
 /* perform stru-stru alignment */
 /* C - column (query) structure,  R - row (database) structure */ 

 /* get local shortcuts */
 contR=prot->conR;
 contT=prot->conT;
 IJcont=prot->IJcont;
 Imod=model->Imod;
 Potn=model->Potn;
 R_ene=prot->R_ene;
 /* seq=prot->seq; */
 C_seq=prot->query_Seq;
 R_seq=prot->native_Seq;
 r12ij=prot->IJr12;
 r6ij=prot->IJr6;
 C_rij=prot->query_IJrij;
 R_rij=prot->native_IJrij;
 Alph=model->Alph;
 Adr=model->Adr;
 Table=dpt->Table;
 Trace=dpt->Trace;
 aGaps=dpt->aGaps;
 Stru_row=dpt->iniRow;
 Stru_col=dpt->iniCol;
 C_score=dpt->C_score;
 R_score=dpt->R_score;
 X_q=prot->X_q;
 Y_q=prot->Y_q;
 Z_q=prot->Z_q;
 
 fea=&fea_data;
 fps=&fps_data;
 norm_d0=fps->norm_d0;
 /* allocate memory for profile distance vector used later for scoring within
	dp table */

 disvec_CR=malloc((op->n_type)*sizeof(float));
 if (!disvec_CR)
 {
	 printf("unable to allocate memory for disvec_CR for stru2stru\n");
	 exit(1);
 }
 F->q_fngrps=fopen(F->f_qfps,"r");
 if (F->q_fngrps==NULL)
 {
	 printf("unable to open query fingrprnt file for x2x alignment\n");
	 return(0);
 }
 F->fngrps=fopen(F->f_fps,"r");
 if (F->fngrps==NULL)
 {
	 printf("unable to open database fingerprnt file for x2x alignment\n");
	 return(0);
 }

 if(op->eij || op->evdw)
   {
     sprintf(Info->Sout," This option works only with TOM-like pot ... \n");
     to_stnd_out(Info->Sout);
     exit(0);
   }

 /* include all the structures at first */
 /* assumption: Query_CM=CM */
 if(op->rm_redund) for(i=0;i<op->max_prot;++i) *(Info->noRed+i)=1;
 
 n=0; m=0; /* used to initiate a decoy on the scaffolding */
 alpha=0;
 ind=1;
 i_prot=0;
 store_kbest=op->kbest;
 op->iprobe=-1;

 /* print first title bar */
 prn_type_of_run(F,op,Info);
 
 /* begin loop over structures C_name (put in the column) from f_qcont */  
 while(read_Query_CM(F,op,IJcont,&l_qstru,C_name,C_rij))
   {
     ++op->iprobe;
     if(op->info2) 
       {
	 sprintf(Info->Sout,
		 "\n Structure alignments ... through (name_length): \n");
	 to_stnd_out(Info->Sout);
       }

     /* read first query sequence and xyz file */
     /* use generic alphabet to avoid translation into ALL */
     read_seq_generic(F->query,&l_qstru,C_name,C_seq,&i_prot,op);
     /* prn_seq_alph(C_seq,&l_qstru,C_name,Alph); */
     read_coor(F->qxyz,op,&l_qstru,C_name,X_q,Y_q,Z_q);
     perc_empty=rm_cent_mass(op,&l_qstru,X_q,Y_q,Z_q);
     if(perc_empty > 10.0)
       {
	 sprintf(Info->Sout," WARNING - please check query struct %s \n",
		 C_name);
	 to_stnd_out(Info->Sout);
	 sprintf(Info->Sout," No coor assigned for %5.2f percent of resid\n",
		 perc_empty);
	 to_stnd_out(Info->Sout);
       }
     n_row=l_qstru;
     i_prot=0;
     wind=0;
     n=0; m=0; /* used to initiate a decoy on the scaffolding */
     alpha=0;  

     /* prepare tables for kth best alignments and shuffling */
     for(i=0;i<op->max_prot;++i) *(Info->Point_to_best+i)=0; 
     Info->ntop=0;
     op->kbest=store_kbest;
       
     /* now get the environmnet and scores of sites in C_score */
/*     get_env(op,Imod,contT,IJcont,&n,C_seq,&n_row,&m,ssh->Ncont,ssh,
	     r12ij,r6ij,Stru_col);

     for(i=0;i<n_row;++i) 
       {
	 beta=i;
	 *(C_score+beta)=get_score(op,Imod,&alpha,C_seq,&n_row,&beta,
				   ssh->Ncont,Potn,ssh->Mcont0,ssh->Mcont1,
				   ssh->Mcont2,&gpen);
       } */

	 /* read fingerprint of query structure from query fingerprint file */
	 if (op->fngrps)
	 {
		 if (!read_fngrps(F->q_fngrps,F->f_qfps,op,fps,C_profile,&scale,&lstru,name,0))
		 {
			 printf("unable to read fingerprints for query stru in x2x alignment\n");
			 exit(1);
		 }
		 /*debugging*/
		 if ((strcmp(name,C_name)!=0) && lstru!=l_qstru)
		 {
			 printf(" name %s != C_name %s and lstru %d != l_qstru %d",name,R_name,lstru,l_qstru);
			 exit(1);
		 }
		 C_profile=fps->C_profile;
	 }

	 	/* debugging */
/*	for (i=0; i<l_qstru; ++i)
	{
		for (j=0; j<op->n_type; ++j)
		{
/*			printf("%4.3f ",fngrprnt);
			++fngrprnt; */
/*			printf("%4.3f ",C_profile[i*op->n_type+j]); 
		}
		printf("\n");
	} 
     
     /* begin loop over structures R_name (put in the row) from f_cont */
     while(read_contMap(F,op,IJcont,&n_col,R_name,r12ij,r6ij,R_rij))  
       {
	 read_seq_generic(F->seq,&n_col,R_name,R_seq,&i_prot,op);
	 pick_stru=1;
	 ++wind;
	 if(op->info2) report_progress(Info,&wind,R_name,&n_col);
	 ++ind;
	 n=0; m=0;
	 alpha=0;

	 /* read fingerprint of database structure from database fingerprint file 
		now outside of if pick loop */
	 if (op->fngrps)
	 {
		 if (!read_fngrps(F->fngrps,F->f_fps,op,fps,R_profile,&scale,&lstru,name,1))
		 {
			 printf("unable to read fingerprints for database stru in x2x alignment\n");
			 exit(1);
		 }
		 /*debugging*/
		 if ((strcmp(name,R_name)!=0) && lstru!=n_col)
		 {
			 printf(" name %s != R_name %s and lstru %d != n_col %d",name,R_name,lstru,n_col);
			 exit(1);
		 }
		 R_profile=fps->R_profile;
	 }

	 /* if a given structure does not match pick it will be skipped */
	 if(pick_structure(op,R_name,&n_col,&l_qstru))
	   {
	     n_row=l_qstru;
	     ++i_prot;
	     /* now get the environmnet and scores of sites in R_score */
/*	     get_env(op,Imod,contT,IJcont,&n,R_seq,&n_col,&m,ssh->Ncont,ssh,
		     r12ij,r6ij,Stru_row);

	     for(j=0;j<n_col;++j) 
	       {
		 beta=j;
		 *(R_score+beta)=get_score(op,Imod,&alpha,R_seq,&n_col,&beta,
				       ssh->Ncont,Potn,ssh->Mcont0,ssh->Mcont1,
				       ssh->Mcont2,&gpen);
	       } */

		 /* read fingerprint of database structure from database fingerprint file */
/*		 if (op->fngrps)
		 { 
			 read_fngrps(F->fngrps,F->f_fps,op,fps,R_profile,&lstru,name,1);
			 /*debugging*/
/*			 if ((strcmp(name,R_name)!=0) && lstru!=n_col)
			 {
				 printf(" name %s != R_name %s and lstru %d != n_col %d",name,R_name,lstru,n_col);
				 exit(1);
			 }
			 R_profile=fps->R_profile;
		 } */
 
	     /* init Dpt, build initial row and col, fix pre, suf gap_pen */
	     init_dpt(op,Table,Trace,&n_row,&n_col,aGaps); 
	     /* now n_row=Lstru1+1, n_col=Lstru2+1 */
		 
		 fprintf(F->fps_norm,"structural alignment: %s -> %s\n",C_name,R_name);
		 fprintf(F->fps_dv,"structural alignment: %s -> %s\n",C_name,R_name);
		 fprintf(F->fps_norm,"\t");
/*		 for (r=0; r<n_col-1; ++r)
		 {
			fprintf(F->fps_norm," %4d ",r);
		 } 
		 fprintf(F->fps_norm,"\n");*/
	     
	     /* build the DP table */
	     gap_pen_left=op->gap_pen;
	     gap_pen_top=op->gap_pen; 
	     for(i=1;i<n_row;++i) /* stru start at the second  r/c */
	       {  
		 alpha=i-1;
		 fprintf(F->fps_norm,"%4d  ",alpha);
		 for(j=1;j<n_col;++j) 
		   {
		     beta=j-1;

			 /* calculate profile distance vector between site alpha and site beta and its norm */
			 norm_dv=(double)0.0;
			 for (k=0; k<op->n_type; ++k)
			 {
				 *(disvec_CR+k)=*(C_profile+alpha*(op->n_type)+k)-*(R_profile+beta*(op->n_type)+k);
				 norm_dv+=pow(*(disvec_CR+k),2);
			 }
			 norm_dv=sqrt(norm_dv);

/*			 fprintf(F->fps_dv,"C_stru pos %4d  ,R_stru pos %4d  ",alpha,beta);
			 for (s=0; s<op->n_type; ++s)
			 {
				 fprintf(F->fps_dv,"% 05.4f ",*(disvec_CR+s)); 		  
			 }
			 fprintf(F->fps_dv,"\n");
			 if(beta<n_col-2)
			 {
				 fprintf(F->fps_norm,"%4.3f ",norm_dv);
			 } */

			 Tij=-1*((1-(pow((norm_dv/norm_d0),2)))/(2*(1+(pow((norm_dv/norm_d0),2)))));

		     /* get the proper scoring matrix element */
/*		     x=*(C_score+alpha)-*(R_score+beta);	     
		     /* arbitrary heuristics - top of the gaussian rewarded */
/*		     Tij=10.0*(0.7-exp(-9.0*x*x));
		     /* if(i==1) printf(" %7.4f %7.4f %7.4f %7.4f\n",
			                x1,x2,x,,Tij); */
		     /* redefine here gap penalty using Potn(env) */
		     /* decide how to proceed */
		     *(Table+i*n_col+j)=get_cell(op,Table,Trace,&n_col,&i,&j,
					     &gap_pen_left,&gap_pen_top,&Tij);
		   }
		 /* printf("\n start row i=%d\n",i+1); */
/*		 fprintf(F->fps_norm,"\n"); */
	       }
	    
	     /* print DP table */
	     if(op->debug) prn_dpt(Table,&n_row,&n_col);
	     /* trace back the alignment */
	     T_ene=find_best_successor(op,Table,&n_row,&n_col,&i,&j);
	     /* so, now we are at the cell (i,j) having the value T_ene */
	     n=trace_back(op,Table,&n_row,&n_col,&i,&j,&T_ene,Trace,aGaps);
		 
		 /* normalize energy of alignment by dividing by alignment length n */
		 T_ene=(T_ene/(float)n);

	     /* prepare best alignments and estimate their significance */
	     include=1;
	     if(op->loc_a && (((float)n/(float)l_qstru)<op->t_length)) 
	       include=0;
	     if(include) pick_best_aligns(op,Info,&wind,&T_ene); 

	     /* print aligned seqs - aGaps indicates how to move from (i,j)*/
	     op->prnalg_lev=1; /* no details */
	     if(op->info) adjust_prn_level(op);
	     /* print alignment - aGaps indicates how to move from (i,j)*/  
	     prn_align(op,F->out1,Info,R_seq,Alph,C_seq,&n,&i,&j,
		       &T_ene,aGaps,C_name,R_name,&n_row,&n_col); 
	     if(op->info) prep_info(Info,&T_ene,&op->max_histo);

/*	(for read_fngrps withing if pick loop)
		 free(R_profile); */

	   } /* end if(pick) */

	   free(R_profile);
	 
       } /* end of while database structures */

     sprintf(Info->Sout,"\n");
     to_stnd_out(Info->Sout);
     Info->index=ind-1;
     if(Info->ntop<op->kbest) op->kbest=Info->ntop;
      
     /* now get Z_scores for shuffled sequence */
     strcpy(op->qname,C_name);
     get_best_rms(op,F,ssh,prot,Info,model,dpt,C_seq,&l_qstru,X_q,Y_q,Z_q);
     rewind(F->cont);
     rewind(F->rij);
     rewind(F->seq);
	 rewind(F->fngrps);
	 free(C_profile);
     
   } /* end of while all the query structures */  

 /* prepare list of non-redundant entries */
 if(op->rm_redund)
   {
     for(i=0;i<op->max_prot;++i)
       {
	 if(*(Info->noRed+i)==1) 
	   add_name2list(op,F,Info,(Info->Names+i)->mstr);
       }
   }

 free(disvec_CR);
 fclose(F->fps_dv);
 fclose(F->fps_norm);
 
 return(1);
 
}

float rm_cent_mass(k_opt *op, int *l_stru, float *X, float *Y, float *Z)
{
  int i,np_align;
  double xcm,ycm,zcm,norm;
  float perc_empty;
  
  xcm=0.0;
  ycm=0.0;
  zcm=0.0;
  np_align=0;
  
  /* discard "infinities" based on X coor - they are coupled */
  for(i=0;i<*l_stru;++i)
    {
      if(*(X+i)<999.0) 
	{
	  ++np_align;
	  xcm+=*(X+i);
	  ycm+=*(Y+i);
	  zcm+=*(Z+i);
	}
    }
  
  /* calculate the percenatge of site without coordinates assigned */
  perc_empty=100.0*np_align/(float)*l_stru; 
  perc_empty=100.0-perc_empty;
  norm=1.0/(double)np_align;
  xcm=xcm*norm;
  ycm=ycm*norm;
  zcm=zcm*norm;

  /* leave "infinities" untouched */
  for(i=0;i<*l_stru;++i)
    {
      if(*(X+i)<999.0)
	{
	  *(X+i)-=xcm;
	  *(Y+i)-=ycm;
	  *(Z+i)-=zcm;
	}
    }  
  
  return(perc_empty);
  
}

/* one should pass C_column and Q_seq here */
int get_best_rms(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
		 k_info *Info, k_Model *model, k_dynpt *dpt, int *C_seq,
		 int *l_qstru, float *X_q, float *Y_q, float *Z_q)
{
  int i_prot,n,m,i,j,k,n_row,n_col,alpha,beta,lpick,pick_stru,wind;
  unsigned long int ind;
  float T_ene,gap_pen_left,gap_pen_top,gpen,Tij,x,y,x1,x2,ftmp,rms_ave,
        sigma,s_id;
  char C_name[MAXS],R_name[MAXS];
  /* local pointers for shortcuts */
  int *contR,*contT,*IJcont,*Imod,*Alph,*Adr,*Trace,*aGaps,
      *Stru_row,*Stru_col,*R_seq,include,i_str,redundant;
  float *Potn,*R_ene,*r12ij,*r6ij,*C_rij,*R_rij,*C_score,*R_score,*Table,
        *X,*Y,*Z;
  
  int lstru,l_align,r; 
  float *C_profile,*R_profile,*disvec_CR,scale;
  double norm_dv,norm_d0,norm;
  char name[MAXS];

 /* perform stru-stru alignment */
 /* C - column (query) structure,  R - row (database) structure */ 

 /* get local shortcuts */
 contR=prot->conR;
 contT=prot->conT;
 IJcont=prot->IJcont;
 Imod=model->Imod;
 Potn=model->Potn;
 R_ene=prot->R_ene;
 R_seq=prot->native_Seq;
 r12ij=prot->IJr12;
 r6ij=prot->IJr6;
 C_rij=prot->query_IJrij;
 R_rij=prot->native_IJrij;
 Alph=model->Alph;
 Adr=model->Adr;
 Table=dpt->Table;
 Trace=dpt->Trace;
 aGaps=dpt->aGaps;
 Stru_row=dpt->iniRow;
 Stru_col=dpt->iniCol;
 C_score=dpt->C_score;
 R_score=dpt->R_score; 
 X=prot->X_res;
 Y=prot->Y_res;
 Z=prot->Z_res;  
 
 fea=&fea_data;
 fps=&fps_data;
 norm_d0=fps->norm_d0;
 /* allocate memory for profile distance vector used later for scoring within
	dp table */

 disvec_CR=malloc((op->n_type)*sizeof(float));
 if (!disvec_CR)
 {
	 printf("unable to allocate memory for disvec_CR for stru2stru\n");
	 exit(1);
 }
 rewind(F->fngrps);
 C_profile=fps->C_profile;

 if(op->info2)
   {
     sprintf(Info->Sout,
	"\n Calculating RMSD for best structure to structure align ... \n");
     to_stnd_out(Info->Sout);
   }
 
 n=0; m=0; /* used to initiate a decoy on the scaffolding */
 alpha=0;
 ind=1;
 i_str=0;
 i_prot=0;
 redundant=0;
 
 /* op->debug=1; */
 rewind(F->cont);
 rewind(F->rij); 
 rewind(F->seq);
 rewind(F->coor); 
 rms_ave=0.0;
 sigma=0.0;
 /* begin loop over structures R_name (put in the row) from f_cont */
 while(read_contMap(F,op,IJcont,&n_col,R_name,r12ij,r6ij,C_rij))  
   {
     read_seq_generic(F->seq,&n_col,R_name,R_seq,&i_prot,op);
     read_coor(F->coor,op,&n_col,R_name,X,Y,Z);
     pick_stru=1;
     ++wind;
     ++ind;
     n=0; m=0;
     alpha=0;
     if(op->pick) /* skip all except picked if desired */
       {
	 lpick=strlen(op->pstr_name);
	 if(lpick>0 && strncmp(R_name,op->pstr_name,lpick)!=0)
	   *(Info->Point_to_best+i_str)=0;
       }

	 /* read fingerprint of database structure from database fingerprint file */
	 if (op->fngrps)
	 { 
		 if (!read_fngrps(F->fngrps,F->f_fps,op,fps,R_profile,&scale,&lstru,name,1))
		 {
			 printf("unable to read fingerprints for database stru for best rmsd\n");
			 exit(1);
		 }
		 /*debugging*/
		 if ((strcmp(name,R_name)!=0) && lstru!=n_col)
		 {
			 printf(" name %s != R_name %s and lstru %d != n_col %d",name,R_name,lstru,n_col);
			 exit(1);
		 }
		 R_profile=fps->R_profile;
	 }
     
     /* if a given structure does not match pick it will be skipped */
     if(*(Info->Point_to_best+i_str))
       {
	 /* now this is done after alignment is found */
	 /* rm_cent_mass(op,&n_col,X,Y,Z); */
	 n_row=*l_qstru;
	 /* ++i_prot; */
	 if(op->info2) 
	   {
	     sprintf(Info->Sout,"  %s_%d  ",R_name,n_col);
	     to_stnd_out(Info->Sout);
	   }
	 /* get the energy of native alignment - coeff manually fixed!!! */
/*	 *(R_ene+i_str)=0.0;
	 for(i=0;i<n_col;++i)
	   {
	     *(R_ene+i_str)+=10.0*(0.7-1.0);
	   }
	 
	 /* now get the environment and scores of sites in R_score */
/*	 get_env(op,Imod,contT,IJcont,&n,R_seq,&n_col,&m,ssh->Ncont,ssh,
		 r12ij,r6ij,Stru_row);

	 for(j=0;j<n_col;++j) 
	   {
	     beta=j;
	     *(R_score+beta)=get_score(op,Imod,&alpha,R_seq,&n_col,&beta,
				       ssh->Ncont,Potn,ssh->Mcont0,ssh->Mcont1,
				       ssh->Mcont2,&gpen);
	   } 

	 /* read fingerprint of database structure from database fingerprint file 
		before being moved outside if pick loop */
/*	 if (op->fngrps)
	 { 
		 read_fngrps(F->fngrps,F->f_fps,op,fps,R_profile,&lstru,name,1);
		 /*debugging*/
/*		 if ((strcmp(name,R_name)!=0) && lstru!=n_col)
		 {
			 printf(" name %s != R_name %s and lstru %d != n_col %d",name,R_name,lstru,n_col);
			 exit(1);
		 }
		 R_profile=fps->R_profile;
	 } */

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
		 
		 /* calculate profile distance vector between site alpha and site beta and its norm */
		 norm_dv=(double)0.0;
		 for (k=0; k<op->n_type; ++k)
		 {
			 *(disvec_CR+k)=*(C_profile+alpha*(op->n_type)+k)-*(R_profile+beta*(op->n_type)+k);
			 norm_dv+=pow(*(disvec_CR+k),2);
		 }
		 norm_dv=sqrt(norm_dv);

		 Tij=-1*((1-(pow((norm_dv/norm_d0),2)))/(2*(1+(pow((norm_dv/norm_d0),2)))));

		 /* get the proper scoring matrix element */
/*		 x=*(C_score+alpha)-*(R_score+beta);	     
		 /* arbitrary heuristics - top of the gaussian rewarded */
/*		 Tij=10.0*(0.7-exp(-9.0*x*x));
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

	 /* normalize energy of alignment by dividing by alignment length n */
	 T_ene=(T_ene/(float)n); 
	 
	 /* get energy of native alignment */
	 *(R_ene+i_str)=0.0;
	 for (r=0; r<lstru; ++r)
	 {
		 *(R_ene+i_str)+=(-1.0*0.5);
	 }
	 *(R_ene+i_str)=(*(R_ene+i_str)/lstru); 

	 /* print aligned seqs - aGaps indicates how to move from (i,j)*/
	 /* define first format of prinitng */
	 op->prnalg_lev=1; /* no details */
	 if(op->info) adjust_prn_level(op);
	 ftmp=get_rms(op,F,Info,R_ene,&i_str,&n,&i,&j,aGaps,C_rij,R_rij,
		       &n_row,&n_col,X_q,Y_q,Z_q,X,Y,Z);
	 /* computes moments of the distribution of rmsd's */
	 rms_ave+=ftmp; 
	 sigma+=ftmp*ftmp;
	 /* here if smaller than RMS threshold */
	 if(ftmp < op->t_rms)
	   {
	     
	     s_id=prn_align(op,F->balg,Info,R_seq,Alph,C_seq,&n,&i,&j,
			    &T_ene,aGaps,op->qname,R_name,&n_row,&n_col);
	     /* high seq identity or small RMSD result in exclusion */
	     if((s_id>100*op->t_seqid) || (ftmp<op->t_rms_db)) 
	       {
		 redundant=1;
		 /* if nored probe structure excludes target structures */
		 if(op->rm_redund) 
		   {
		     if((*(Info->noRed+op->iprobe)==1) && (op->iprobe!=i_str)) 
		       *(Info->noRed+i_str)=0;
		     /* ones excluded probe can no longer be a killer */
		   }
	       }
	   }
	 

       } /* end if(pick) */
	 
     ++i_str;

	 free(R_profile);

   } /* end of while database structures */

   free(disvec_CR);
 
 /* so only those that were picked had a chance to become "redundant" */
 if(op->build_db && !redundant) add_name2list(op,F,Info,op->qname);

 /* print best scores */
 pick_best_zscore(op,Info); 
 print_best_stru(op,F,Info,&rms_ave,&sigma); 
     
 return(0);
 
}

void add_name2list(k_opt *op, a_file *F, k_info *Info, char *name)
{
  
  fprintf(F->newlist,"%s\n",name);
  
}

float get_rms(k_opt *op, a_file *F, k_info *Info, float *R_ene,
	      int *i_str, int *n, int *i, int *j, int *aGaps, 
	      float *C_rij, float *R_rij, int *n_row, int *n_col,
	      float *X_q, float *Y_q, float *Z_q, float *X, float *Y, float *Z)
{
  int k,ii,jj,shift_z2;
  float T_ene,rms;

  rms=RMAX;
  
  /* calculate now the actual z_score and send them via Info */
  /* compute also average RMS for estimate of significance */
  for(k=0;k<op->kbest;++k)
    {
      if(*(Info->Best_str+k)==*i_str) 
	{
	  ii=*i;
	  jj=*j;
	  T_ene=*(Info->Best_ene+k);
	  /* ftmp=-(*(R_ene+*i_str)-T_ene); */
	  /* perhaps use C_rij, R_rij in the future for dRMS */
	  rms=rmsd(op,n,&ii,&jj,aGaps,n_row,n_col,X_q,Y_q,Z_q,X,Y,Z);
	  *(Info->Zscore1+k)=rms;	    
	  /* *(Info->Zscore2+k)=ftmp; */
	  if(op->info2) 
	    {
	      sprintf(Info->Sout," e_query= %5.2f   e_nat= %5.2f \n",
		      T_ene,*(R_ene+*i_str));
	      to_stnd_out(Info->Sout);
	    }
	}
    }   

  return(rms);
  
}

/* compute RMSD betwen aligned structures */
float rmsd(k_opt *op, int *n, int *i, int *j, int *aGaps, int *n_row, int *n_col,
	   float *X_q, float *Y_q, float *Z_q, float *X, float *Y, float *Z)
{
 int k,ii,jj,kk,i_loc,j_loc,np_align;
 double ftmp,f1,f2,f3;
 double kabsch[3][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
 double kab2[3][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
 double b[3][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
 double rot[3][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
 double eigval[3]= {0.0,0.0,0.0};
 double rms,dtmp,norm,xqcm,yqcm,zqcm,xcm,ycm,zcm;
 
 rms=0.0; 
 i_loc=*i;
 j_loc=*j;

 /* remove center of mass for the aligned parts */
 xqcm=0.0;
 yqcm=0.0;
 zqcm=0.0;
 xcm=0.0;
 ycm=0.0;
 zcm=0.0;
 np_align=0;
 for(k=*n-1;k>=0;--k) 
   {
     /* check whether you reached the end of seq/stru */
     if(*i>=*n_row-1) *i=*n_row-2;
     if(*j>=*n_col-1) *j=*n_col-2;
     if(*(aGaps+k)==0) /* diagonal */
       {
	 ii=*i;
	 jj=*j;
	 if((*i<*n_row-1) && (*j<*n_col-1))
	   {
	     if((*(X_q+ii)<999.0) && (*(X+jj)<999.0))
	       {
		 xqcm += *(X_q+ii);
		 yqcm += *(Y_q+ii);
		 zqcm += *(Z_q+ii);
		 xcm += *(X+jj);
		 ycm += *(Y+jj);
		 zcm += *(Z+jj);
		 ++np_align;
	       }
	   }
	 
	 ++*i; ++*j;
       }
     else 
       {
	 if(*(aGaps+k)==1) ++*j; /* horizontal */
	 if(*(aGaps+k)==2) ++*i; /* vertical */ 
       }
   }
  norm=1.0/(double)np_align;
  xcm=xcm*norm;
  ycm=ycm*norm;
  zcm=zcm*norm;
  xqcm=xqcm*norm;
  yqcm=yqcm*norm;
  zqcm=zqcm*norm;

  /* subtract center of mass of database stru, leave "infinities" untouched */
  for(jj=0;jj<*n_col-1;++jj)
    {
      if(*(X+jj)<999.0)
	{
	  *(X+jj)-=xcm;
	  *(Y+jj)-=ycm;
	  *(Z+jj)-=zcm;
	}
    }
  /*
  printf(" xcm=%f ycm=%f zcm=%f   xqcm=%f yqcm=%f zqcm%f \n",xcm,ycm,zcm,
	 xqcm,yqcm,zqcm);
	 */
 
 /* in the first passage construct the Kabsch matrix */
 *i=i_loc;
 *j=j_loc;
 for(k=*n-1;k>=0;--k) 
   {
     /* check whether you reached the end of seq/stru */
     if(*i>=*n_row-1) *i=*n_row-2;
     if(*j>=*n_col-1) *j=*n_col-2;
     if(op->debug) printf(" gap counter i=%d aGaps=%d \n",k,*(aGaps+k));
     if(*(aGaps+k)==0) /* diagonal */
       {
	 ii=*i;
	 jj=*j;
	 /* add (in column order) to lower triangle */
	 if((*i<*n_row-1) && (*j<*n_col-1))
	   {
	     if((*(X_q+ii)<999.0) && (*(X+jj)<999.0))
	       {
		 f1= *(X_q+ii) - xqcm;
		 f2= *(Y_q+ii) - yqcm;
		 f3= *(Z_q+ii) - zqcm;
		 
		 kabsch[0][0] += f1 * (*(X+jj));
		 kabsch[1][0] += f2 * (*(X+jj));
		 kabsch[2][0] += f3 * (*(X+jj));

		 kabsch[0][1] += f1 * (*(Y+jj));
		 kabsch[1][1] += f2 * (*(Y+jj));
		 kabsch[2][1] += f3 * (*(Y+jj));
		 
		 kabsch[0][2] += f1 * (*(Z+jj));
		 kabsch[1][2] += f2 * (*(Z+jj));
		 kabsch[2][2] += f3 * (*(Z+jj));
	       }
	   }
	 
	 ++*i; ++*j;
       }
     else 
       {
	 if(*(aGaps+k)==1) ++*j; /* horizontal */
	 if(*(aGaps+k)==2) ++*i; /* vertical */ 
       }
   }
 /*
 for (ii=0;ii<3;ii++)     
   {
     for (jj=0;jj<3;jj++) 
       printf(" %6.2f ",kabsch[ii][jj]);
     printf(" \n");
   }
   */

 /* multiply kabsch by its transpose */
 for (ii=0;ii<3;ii++) 
   for (jj=0;jj<3;jj++) 
     for (kk=0;kk<3;kk++)
       kab2[ii][jj] += kabsch[kk][ii] * kabsch[kk][jj];
 /*
  for (ii=0;ii<3;ii++)     
   {
     for (jj=0;jj<3;jj++) 
       printf(" %6.2f ",kab2[ii][jj]);
     printf(" \n");
   }
   */

 /* diagonalize kab2 - eigenvectors are returned as columns of kab2 */
 diag_kabsch(op,kab2,eigval);

 /*
  for (ii=0;ii<3;ii++)     
   {
     for (jj=0;jj<3;jj++) 
       printf(" %6.2f ",kab2[ii][jj]);
     printf(" \n");
   }
   */

 /* compute b vectors */
 for (ii=0;ii<3;ii++)     
   for (jj=0;jj<3;jj++) 
     {
       norm=1.0/sqrt(eigval[jj]); 
       for (kk=0;kk<3;kk++)
	 b[ii][jj] += kabsch[ii][kk] * kab2[kk][jj] * norm;
     }
 
 
 /* now get the actual rotation matrix and modify X, Y, Z */
 for (ii=0;ii<3;ii++)     
   for (jj=0;jj<3;jj++) 
     for (kk=0;kk<3;kk++)
       rot[ii][jj] += b[ii][kk] * kab2[jj][kk];
 /*
 for (ii=0;ii<3;ii++)     
   {
     for (jj=0;jj<3;jj++) 
       printf(" %6.2f ",rot[ii][jj]);
     printf(" \n");
   }
   */

 /* transform coordinates of the template (row structure) */
 for (jj=0;jj<*n_col-1;jj++) 
   {
     if(*(X+jj)<999.0)
       {
	 f1 = *(X+jj)*rot[0][0] + *(Y+jj)*rot[0][1] + *(Z+jj)*rot[0][2];
	 f2 = *(X+jj)*rot[1][0] + *(Y+jj)*rot[1][1] + *(Z+jj)*rot[1][2];
	 f3 = *(X+jj)*rot[2][0] + *(Y+jj)*rot[2][1] + *(Z+jj)*rot[2][2];  
	 *(X+jj)=f1;
	 *(Y+jj)=f2;
	 *(Z+jj)=f3;    
       }
   }
 
 /* in the final passage calculate rms for the aligned pairs */
 *i=i_loc;
 *j=j_loc;
 np_align=0;
 for(k=*n-1;k>=0;--k) 
   {
     /* check whether you reached the end of seq/stru */
     if(*i>=*n_row-1) *i=*n_row-2;
     if(*j>=*n_col-1) *j=*n_col-2;
     if(op->debug) printf(" gap counter i=%d aGaps=%d \n",k,*(aGaps+k));
     if(*(aGaps+k)==0) /* diagonal */
       {
	 ii=*i;
	 jj=*j;
	 if((*i<*n_row-1) && (*j<*n_col-1)) /* skip last */
	   {
	     if((*(X_q+ii)<999.0) && (*(X+jj)<999.0))
	       {
		 f1= *(X_q+ii) - xqcm - *(X+jj);
		 f2= *(Y_q+ii) - yqcm - *(Y+jj);
		 f3= *(Z_q+ii) - zqcm - *(Z+jj);
		 rms += f1*f1 + f2*f2 + f3*f3;
		 ++np_align;
	       }
	   }
	 ++*i; ++*j; 
       }
     else 
       {
	 if(*(aGaps+k)==1) ++*j; /* horizontal */
	 if(*(aGaps+k)==2) ++*i; /* vertical */ 
       }
   }
 dtmp=(double)np_align;
 if(dtmp<epsilon) dtmp=epsilon;
 rms=sqrt(rms/dtmp);
 /* printf(" np_alg= %f   RMS= %f\n",dtmp,rms); */
 return((float)rms);
 
}

/* the next 5 (or so) functions were taken from Numerical Recipes */
/* flowers for the NR guys */
/* note small modifications with respect to original code */
int diag_kabsch(k_opt *op, double kab[3][3], double *eigval)
{
  int i,j,k,kk,l,ll,nrot;
  double a[3][3],b[3][3],c[3][3]; /* let c be a buffer for shift */
  double *d,**v,**e;

  for (j=0;j<3;j++) 
    for (k=0;k<3;k++) 
	a[k][j]=kab[k][j];

  /* allocate such that indices start from 1 instead of 0 */
  d=vector(1,3);
  e=convert_matrix(&a[0][0],1,3,1,3);
  v=convert_matrix(&b[0][0],1,3,1,3);
  jacobi(e,3,d,v,&nrot);
  eigsrt(d,v,3);
  
  /* destroy now the original matrix kab */
  for (j=1;j<=3;j++) *(eigval+j-1)=d[j];
  for (j=1;j<=3;j++) 
    for (k=1;k<=3;k++) 
      kab[k-1][j-1]=v[k][j];

  free_convert_matrix(e,1,3,1,3);
  free_convert_matrix(v,1,3,1,3);
  free_vector(d,1,3);
  return(1);
}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

/* diagonalize matrix using Jacobi method */
void jacobi(double **a, int n, double d[], double **v, int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  b=vector(1,n);
  z=vector(1,n);
  for (ip=1;ip<=n;ip++) 
    {
      for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
      v[ip][ip]=1.0;
    }
  for (ip=1;ip<=n;ip++) 
    {
      b[ip]=d[ip]=a[ip][ip];
      z[ip]=0.0;
    }
  *nrot=0;
  for (i=1;i<=50;i++) 
    {
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) 
	{
	  for (iq=ip+1;iq<=n;iq++)
	    sm += fabs(a[ip][iq]);
	}
      if (sm == 0.0) 
	{
	  free_vector(z,1,n);
	  free_vector(b,1,n);
	  return;
	}
      if (i < 4) tresh=0.2*sm/(n*n);
      else tresh=0.0;
      for (ip=1;ip<=n-1;ip++) 
	{
	  for (iq=ip+1;iq<=n;iq++) 
	    {
	      g=100.0*fabs(a[ip][iq]);
	      if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
		  && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
		a[ip][iq]=0.0;
	      else if (fabs(a[ip][iq]) > tresh) 
		{
		  h=d[iq]-d[ip];
		  if ((double)(fabs(h)+g) == (double)fabs(h))
		    t=(a[ip][iq])/h;
		  else 
		    {
		      theta=0.5*h/(a[ip][iq]);
		      t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
		      if (theta < 0.0) t = -t;
		    }
		  c=1.0/sqrt(1+t*t);
		  s=t*c;
		  tau=s/(1.0+c);
		  h=t*a[ip][iq];
		  z[ip] -= h;
		  z[iq] += h;
		  d[ip] -= h;
		  d[iq] += h;
		  a[ip][iq]=0.0;
		  for (j=1;j<=ip-1;j++) 
		    {
		      ROTATE(a,j,ip,j,iq)
		    }
		  for (j=ip+1;j<=iq-1;j++) 
		    {
		      ROTATE(a,ip,j,j,iq)
		    }
		  for (j=iq+1;j<=n;j++) 
		    {
		      ROTATE(a,ip,j,iq,j)
		    }
		  for (j=1;j<=n;j++) 
		    {
		      ROTATE(v,j,ip,j,iq)
		    }
		  ++(*nrot);
		}
	    }
	}
      for (ip=1;ip<=n;ip++) 
	{
	  b[ip] += z[ip];
	  d[ip]=b[ip];
	  z[ip]=0.0;
	}
    }
  nrerror("Too many iterations in routine jacobi");
}

#undef ROTATE

/* sort eigenvectors in descending eigenvalue order */
void eigsrt(double d[], double **v, int n)
{
  int k,j,i;
  double p;

  for (i=1;i<n;i++) 
    {
      p=d[k=i];
      for (j=i+1;j<=n;j++)
	if (d[j] >= p) p=d[k=j];
      if (k != i) 
	{
	  d[k]=d[i];
	  d[i]=p;
	  for (j=1;j<=n;j++) 
	    {
	      p=v[j][i];
	      v[j][i]=v[j][k];
	      v[j][k]=p;
	    }
	}
    }
}

#define NR_END 1

/* Numerical Recipes standard error handler */
void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

/* allocate a float vector with subscript range v[nl..nh] */
double *vector(long nl, long nh)
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}

/* free a float vector allocated with vector() */
void free_vector(double *v, long nl, long nh)
{
  free((char* ) (v+nl-NR_END));
}

double **convert_matrix(double *a, long nrl,long nrh, long ncl, long nch)
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
     declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
     and ncol=nch-ncl+1. The routine should be called with the address
     &a[0][0] as the first argument. */

  /* allocate pointers to rows */
  m=(double **) malloc((unsigned int) ((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure in convert_matrix()");
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl]=a-ncl;
  for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

/* free a matrix allocated by convert_matrix() */
void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch)
{
  free((char* ) (b+nrl-NR_END));
}

#undef NR_END

/* to print the list of best structures according to Z_score */
void print_best_stru(k_opt *op, a_file *F, k_info *Info, float *rms_ave,
		     float *sigma)
{
  int i,j,shift_z2,hits;
  float z1,z2,norm;
  char sig[MAXS],match[MAXS];
  
  strcpy(match,"structure");
  sprintf(Info->Sout,
	  "\n The best matching %ss (for query structure %s) are: \n \n",
	  match,op->qname);
  to_stnd_out(Info->Sout);
  to_file(F->best,Info->Sout);
  hits=0;
  
  norm=(float)op->kbest;
  if(norm==0.0) norm=1.0;

  *rms_ave=*rms_ave/norm; 
  *sigma=sqrt(fabs(*sigma/norm-((*rms_ave)*(*rms_ave))));
  if(*sigma<epsilon) *sigma=epsilon;
  *sigma=1.0/(*sigma);
  /* printf(" norm=%f rms_av=%f sigma=%f \n",norm,*rms_ave,*sigma); */
  

  for(i=0;i<op->kbest;++i)
    {
      j=*(Info->toBest+i);
      z1=*(Info->Zscore1+j); /* so z1 is actually an rms dist */
      z2=(*rms_ave-z1)*(*sigma);
      
      if(z1<op->t_rms)
	{
	  ++hits;
	  /* setup the siginificance level using RMSD (z1) and Z-score (z2) */
	  strcpy(sig,"very_low");
	  if(z1<15.0) 
	    {
	      strcpy(sig,"low");
	      if(z2<1.0) strcpy(sig,"very_low");
	    }
	  if(z1<12.0) 
	    {
	      strcpy(sig,"HIGH/low");
	      if(z2<1.5) strcpy(sig,"low");	      
	    }
	  if(z1<8.0) 
	    {
	      strcpy(sig,"HIGH");
	      if(z2<2.0) strcpy(sig,"HIGH/low");
	    }
	  if(z1<4.0) strcpy(sig,"VERY_HIGH");
	  /* scale down the significance if energy positive */
	  if(*(Info->Best_ene+j)>0.0) strcpy(sig,"very_low");
	  sprintf(Info->Sout,
	      " %s\t ene= %8.4f\t rmsd= %5.1f (z_sc= %5.1f)  similarity= %s \n",
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
  else if(op->info2)
    {
      sprintf(Info->Sout," average RMS= %5.2f \n",*rms_ave);      
      to_stnd_out(Info->Sout);      
    }
  
  
}
  
