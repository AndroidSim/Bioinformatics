/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000   v. 2.000                                    */
/*                                                                */
/*       MEMORY ALLOCATION, OPTION SETUP, FILES PREPARATION       */
 
#include "loopp.h"
   


/* alloc_mem reads f_seq file to define memory requirements and allocates */
int alloc_mem(a_file *F, k_prot *prot, k_opt *op, k_info *Info,
	      k_Model *model, solv_sh *ssh, k_dynpt *dpt)
{
 int n_prot,n_res,m,Lmax,Lseq,i_buffer,length_histo[12],ibin;
 char T_name[MAXS];

 /* check max_length and get some info in terms of lengths distribution */
 F->seq=fopen(F->f_seq,"r");
 if(F->seq==NULL)
   {
     sprintf(Info->Sout," I cannot allocate memory - file %s not found \n",
	     F->f_seq);
     to_stnd_out(Info->Sout);
     exit(1);
   }
 Lmax=0;
 n_res=0;
 n_prot=0;
 for(ibin=0;ibin<12;++ibin) length_histo[ibin]=0; 
 while(get_length_from_seq(F->seq,&Lseq,T_name,&Lmax,&n_res,&n_prot,op))
   {
     /* collect further stat here - rough histogram of seq lengths */
     if(Lseq<=100)
       {
	 if(Lseq<=50) ++length_histo[0];
	 else ++length_histo[1];
       }
     else
       {
	 ibin=1+Lseq/100;
	 if(ibin>11) ibin=11;
	 ++length_histo[ibin];
       }
   }
 fclose(F->seq);
 op->max_length=Lmax;
 op->max_prot=n_prot;
 /* adjust the range of threading loop */
 if(op->end==0) op->end=op->max_prot;
 if(op->beg>op->max_prot) op->beg=op->max_prot;

 /* WARNING!!! - this way of fixing max_cont may fail */
 if(op->r_cut>50.0) op->r_cut=50.0;
 if(op->r_cut<15.0) op->max_cont=(int)(0.5*op->r_cut*op->r_cut);
 else op->max_cont=(int)((op->r_cut/50.0)*op->max_length);

 /* print info about the setup */
 if(op->info2)
   {
     sprintf(Info->Sout,
	     " Allocation of memory (based on files %s and Model): \n",
	     F->f_seq);
     to_stnd_out(Info->Sout);
     sprintf(Info->Sout,"\t n_type (of letters) = %d \n",op->n_type);
     to_stnd_out(Info->Sout);
     sprintf(Info->Sout,"\t n_para = %d \n",op->n_par); 
     to_stnd_out(Info->Sout);
     sprintf(Info->Sout,"\t n_prot = %d \n",op->max_prot); 
     to_stnd_out(Info->Sout);
     sprintf(Info->Sout,"\t max_length = %d \n",op->max_length);
     to_stnd_out(Info->Sout);
     sprintf(Info->Sout,"\t max_cont = %d \n",op->max_cont);
     to_stnd_out(Info->Sout);
     sprintf(Info->Sout,"\t length_dist:   up to 50: %d  51-100: %d \n\t",
	     length_histo[0],length_histo[1]);
     to_stnd_out(Info->Sout);
     for(ibin=2;ibin<10;++ibin) 
       {
	 sprintf(Info->Sout," %d-%d: %d ",1+100*(ibin-1),100*ibin,
		 length_histo[ibin]); 
	 to_stnd_out(Info->Sout);
       }
     sprintf(Info->Sout,"\n\t 901-1000: %d   more than 1001: %d \n",
	     length_histo[10],length_histo[11]);
     to_stnd_out(Info->Sout);
   }

 /* 
 i_buffer=10;
 op->max_length=op->max_length+i_buffer;
 op->max_prot=op->max_prot+i_buffer;
 */
 
 /* allocate arrays storing alphabet when different from generic */
 model->Alph=malloc(op->n_type*sizeof(int));
 if(model->Alph==NULL) return(0); 

 model->Adr=malloc(MAXS*sizeof(int)); 
 /* assumption: MAXS is becoming a limit for the size of generic alphabet */
 if(model->Adr==NULL) return(0); 
 
 model->Potn=malloc((op->n_par)*sizeof(float));
 if(model->Potn==NULL) return(0);

 model->Imod=malloc(op->n_type*sizeof(int));
 if(model->Imod==NULL) return(0);

 /* allocate arrays storing seq and coordinates */
 prot->seq=malloc(op->max_length*sizeof(int));
 if(prot->seq==NULL) return(0); 

 prot->native_Seq=malloc(op->max_length*sizeof(int));
 if(prot->native_Seq==NULL) return(0); 

 prot->query_Seq=malloc(op->max_length*sizeof(int));
 if(prot->query_Seq==NULL) return(0); 

 /* allocate arrays for dynamic programming */
 if(op->align || op->gaps || op->strucmp)
   {
     m=(op->max_length+1)*(op->max_length+1);

     dpt->iniCol=malloc(op->max_length*sizeof(int));
     if(dpt->iniCol==NULL) return(0); 

     dpt->iniRow=malloc(op->max_length*sizeof(int));
     if(dpt->iniRow==NULL) return(0); 

     dpt->Trace=malloc(m*sizeof(int));
     if(dpt->Trace==NULL) return(0);

     dpt->aGaps=malloc(m*sizeof(int));
     if(dpt->aGaps==NULL) return(0);
 
     dpt->Table=malloc(m*sizeof(float));
     if(dpt->Table==NULL) return(0); 
     
     dpt->C_score=malloc(op->max_length*sizeof(float));
     if(dpt->C_score==NULL) return(0); 

     dpt->R_score=malloc(op->max_length*sizeof(float));
     if(dpt->R_score==NULL) return(0);     

     Info->Point_to_best=malloc(op->max_prot*sizeof(int));
     if(Info->Point_to_best==NULL) return(0);

     Info->toBest=malloc(op->kbest*sizeof(int));
     if(Info->toBest==NULL) return(0);
     
     Info->Best_str=malloc(op->kbest*sizeof(int));
     if(Info->Best_str==NULL) return(0);
     
     Info->Best_ene=malloc(op->kbest*sizeof(float));
     if(Info->Best_ene==NULL) return(0);

     Info->Zscore1=malloc(op->kbest*sizeof(float));
     if(Info->Zscore1==NULL) return(0);

     Info->Zscore2=malloc(op->kbest*sizeof(float));
     if(Info->Zscore2==NULL) return(0);

     Info->Score_shuffle=malloc(op->n_shuffle*sizeof(float));
     if(Info->Score_shuffle==NULL) return(0);
     
   }
 
 prot->X_res=malloc(op->max_length*sizeof(float));
 if(prot->X_res==NULL) return(0); 

 prot->Y_res=malloc(op->max_length*sizeof(float));
 if(prot->Y_res==NULL) return(0); 

 prot->Z_res=malloc(op->max_length*sizeof(float));
 if(prot->Z_res==NULL) return(0); 

 if(op->strucmp)
   {
     prot->X_q=malloc(op->max_length*sizeof(float));
     if(prot->X_q==NULL) return(0); 

     prot->Y_q=malloc(op->max_length*sizeof(float));
     if(prot->Y_q==NULL) return(0); 

     prot->Z_q=malloc(op->max_length*sizeof(float));
     if(prot->Z_q==NULL) return(0); 
   }
 
 /* allocate arrays storing FingerPrints */
 if(op->fngrps)
   {
     prot->FngrPs=malloc(op->max_length*sizeof(int));
     if(prot->FngrPs==NULL) return(0); 
   }
 
 /* allocate arrays storing ContMaps */
 prot->IJcont=malloc((op->max_cont*op->max_length)*sizeof(int));
 if(prot->IJcont==NULL) return(0);

 /* allocate arrays storing proteins according to contact types */
 prot->conR=malloc(op->n_par*sizeof(int));
 if(prot->conR==NULL) return(0);

 prot->conT=malloc(op->n_par*sizeof(int));
 if(prot->conT==NULL) return(0);

 /* allocate array storing native energies */
 prot->R_ene=malloc(op->max_prot*sizeof(float));
 if(prot->R_ene==NULL) return(0);

 /* allocate arrays storing proteins in case of continuous models (evdw) */
 if(op->evdw)
   {
     /* allocate arrays storing proteins according to contact types */
     prot->conRvdw=malloc(op->n_par*sizeof(float));
     if(prot->conRvdw==NULL) return(0);

     prot->conTvdw=malloc(op->n_par*sizeof(float));
     if(prot->conTvdw==NULL) return(0);

     /* allocate arrays storing ContMaps */
     prot->IJr12=malloc((op->max_cont*op->max_length)*sizeof(float));
     if(prot->IJr12==NULL) return(0);

     prot->IJr6=malloc((op->max_cont*op->max_length)*sizeof(float));
     if(prot->IJr6==NULL) return(0);

   }
 
 if(op->distance)
   {

     prot->native_IJrij=malloc((op->max_cont*op->max_length)*sizeof(float));
     if(prot->native_IJrij==NULL) return(0);

     prot->query_IJrij=malloc((op->max_cont*op->max_length)*sizeof(float));
     if(prot->query_IJrij==NULL) return(0);

   }
   
 /* allocate arrays storing solvation shells description in terms
    of number of neighbours (Ncont) and number of neigh having 
    themselves specific number of contacts*/

 ssh->Ncont=malloc(op->max_length*sizeof(int)); /* use this in seq-seq */
 if(ssh->Ncont==NULL) return(0);

 if(op->ein || op-> einm)
   {
     ssh->Mcont0=malloc(op->max_length*sizeof(int));
     if(ssh->Mcont0==NULL) return(0); 

     ssh->Mcont1=malloc(op->max_length*sizeof(int));
     if(ssh->Mcont1==NULL) return(0); 

     ssh->Mcont2=malloc(op->max_length*sizeof(int));
     if(ssh->Mcont2==NULL) return(0); 
   }
 
 /* allocate arrays for contact statistics and histograms */
 Info->con_all=malloc(op->n_par*sizeof(unsigned long int));
 if(Info->con_all==NULL) return(0);

 Info->con_all_dev=malloc(op->n_par*sizeof(unsigned long int));
 if(Info->con_all_dev==NULL) return(0);

 Info->con_All_decoys=malloc(op->n_par*sizeof(double));
 if(Info->con_All_decoys==NULL) return(0);

 Info->histo=malloc(op->max_histo*sizeof(unsigned long int));
 if(Info->histo==NULL) return(0);

 Info->res_Stat=malloc(op->n_type*sizeof(unsigned long int));
 if(Info->res_Stat==NULL) return(0);

 Info->Names=malloc(op->max_prot*sizeof(str_my));
 if(Info->Names==NULL) return(0);

 Info->noRed=malloc(op->max_prot*sizeof(int));
 if(Info->noRed==NULL) return(0);

 /* allocate space for residue data */
 prot->Res=malloc(op->max_length*sizeof(k_res));
 if(prot->Res==NULL) return(0); 

 /*
 op->max_length=op->max_length-i_buffer;
 op->max_prot=op->max_prot-i_buffer;
 */
 
 return(1);
 
} /* end of allocate_mem */

void prn_type_of_run(a_file *F, k_opt *op, k_info *Info)
{
  int i;
  char r_type[MAXS],glob_loc[MAXS];
  
  strcpy(r_type,"THREADING - sequence to structure");
  strcpy(glob_loc,"global");
  if(op->loc_a) strcpy(glob_loc,"local");
  if(op->align) strcpy(r_type,"SEQUENCE to sequence");
  if(op->strucmp) strcpy(r_type,"STRUCTURE to structure");
  sprintf(Info->Sout,"\n --- %s (%s) alignments RESULTS --- \n \n",
	  r_type,glob_loc);
  to_stnd_out(Info->Sout);
  to_file(F->best,Info->Sout);

  sprintf(Info->Sout,"\n --- %s (%s) ALIGNMENTS (not ordered)  --- \n \n",
	  r_type,glob_loc);
  if(op->info) to_file(F->balg,Info->Sout);

}


void prn_signature(int argc, char *argv[], a_file *F, k_opt *op, k_info *Info)
{
 int i;


 sprintf(Info->Sout,"\n Ver. %5.3f (as of %s) ",VERSION,__DATE__);
 to_stnd_out(Info->Sout);
 to_file(F->out1,Info->Sout);
 if(op->info) to_file(F->info,Info->Sout);
 if(op->info2) to_file(F->iprot,Info->Sout);
 if(op->exam || op->align || op->strucmp) 
   {
     to_file(F->best,Info->Sout);
     to_file(F->balg,Info->Sout);
   } 

 sprintf(Info->Sout,"\n Options that were used:  ");
 to_stnd_out(Info->Sout);
 to_file(F->out1,Info->Sout);
 if(op->info) to_file(F->info,Info->Sout);
 if(op->info2) to_file(F->iprot,Info->Sout);
 if(op->exam || op->align || op->strucmp) 
   {
     to_file(F->best,Info->Sout);
     to_file(F->balg,Info->Sout);
   } 
  
 for(i=0;i<argc;++i) 
   {
     sprintf(Info->Sout,"%s ",argv[i]);
     to_stnd_out(Info->Sout);
     to_file(F->out1,Info->Sout);
     if(op->info) to_file(F->info,Info->Sout);
     if(op->info2) to_file(F->iprot,Info->Sout);
     if(op->exam || op->align || op->strucmp) 
       {
	 to_file(F->best,Info->Sout);
	 to_file(F->balg,Info->Sout);
       } 
   }

 sprintf(Info->Sout,"\n Database of %d structures was used. ",op->max_prot);
 to_stnd_out(Info->Sout);
 to_file(F->out1,Info->Sout);
 if(op->info) to_file(F->info,Info->Sout);
 if(op->info2) to_file(F->iprot,Info->Sout);
 if(op->exam || op->align || op->strucmp) 
   {
     to_file(F->best,Info->Sout);
     to_file(F->balg,Info->Sout);
   }

 sprintf(Info->Sout,"\n Please, report bugs to meller@cs.cornell.edu \n");
 to_stnd_out(Info->Sout);
 
}


/* free memory and close files */
void clean_up(a_file *F, k_prot *prot, k_opt *op, k_info *Info,
	      k_Model *model, solv_sh *ssh, k_dynpt *dpt)
{
 int n_prot,n_res,Lmax,Lseq;
 char T_name[MAXS];

 sprintf(Info->Sout,"\n Bye, bye! \n");
 to_stnd_out(Info->Sout);

 /* close the files */
 fclose(F->coor);
 fclose(F->seq);
 fclose(F->type);
 fclose(F->cont);
 if(op->evdw) fclose(F->vdw_a);
 if(op->evdw) fclose(F->vdw_b);
 if(!op->def_pot) fclose(F->pot);
 fclose(F->out1);
  
 if(op->wDc) fclose(F->out2);
 if(op->mps) fclose(F->mps);
 if(op->info) fclose(F->info);
 if(op->info2) fclose(F->iprot);
 if(op->list) fclose(F->list);
 if(op->build_db || op->rm_redund) fclose(F->newlist);
 if(op->distance) fclose(F->rij);
 if(op->fngrps) fclose(F->fngrps);
 if(op->freeze && op->wDc) fclose(F->rhs);
 
 if(op->exam || op->align || op->strucmp) 
   {
     fclose(F->exam);
     fclose(F->best);
     fclose(F->balg); 
     fclose(F->query);
     fclose(F->qxyz); 
     fclose(F->qcont);
     fclose(F->qrij); 
   }
 
 /* free arrays storing seq and coordinates */
 free(prot->native_Seq);
 free(prot->query_Seq);
 free(prot->seq);
 free(prot->X_res);
 free(prot->Y_res);
 free(prot->Z_res);
 if(op->strucmp)
   {
     free(prot->X_q);
     free(prot->Y_q);
     free(prot->Z_q);
   }


 
 if(op->fngrps)
   {
     free(prot->FngrPs);
   }
 
 /* free arrays storing ContMaps */
 free(prot->IJcont);
 /* free arrays storing proteins according to contact types */
 free(prot->conR);
 free(prot->conT);
 /* free array storing native energies */
 free(prot->R_ene);
 /* free arrays storing proteins in case of continuous models (evdw) */
 if(op->evdw)
   {
     /* free arrays storing proteins according to contact types */
     free(prot->conRvdw);
     free(prot->conTvdw);
     /* free arrays storing ContMaps */
     free(prot->IJr12);
     free(prot->IJr6);
   }
 /* free arrays storing alphabet when different from generic */
 free(model->Alph);
 free(model->Adr);
 free(model->Potn);
 free(model->Imod);

 if(op->distance)
   {
     free(prot->native_IJrij);
     free(prot->query_IJrij);
   }
 
 /* free arrays storing solvation shells */
 free(ssh->Ncont);
 if(op->ein || op-> einm)
   {
     free(ssh->Mcont0);
     free(ssh->Mcont1);
     free(ssh->Mcont2);
   }

 /* free arrays for contact statistics and histograms */
 free(Info->con_all);
 free(Info->con_all_dev);
 free(Info->con_All_decoys);
 free(Info->histo);
 free(Info->res_Stat);
 free(Info->Names);
 free(Info->noRed);

 free(prot->Res);
 
 /* free arrays for dynamic programming */
 if(op->align || op->gaps || op->strucmp)
   {
     free(dpt->iniCol) ; 
     free(dpt->iniRow) ; 
     free(dpt->Trace) ;
     free(dpt->aGaps) ;
     free(dpt->Table) ; 
     free(dpt->C_score) ; 
     free(dpt->R_score) ;     
     free(Info->Point_to_best) ;
     free(Info->toBest) ;
     free(Info->Best_str) ;
     free(Info->Best_ene) ;
     free(Info->Zscore1) ;
     free(Info->Zscore2) ;
     free(Info->Score_shuffle) ;
   }
  
}

/* just small error function */
void alloc_err(char *name)
{
  printf("Sorry - cannot allocate array %s\n",name);
  exit(1);  
}

/* read seq file to define memory requirements */
int get_length_from_seq(FILE *Fseq,int *Lseq,char *T_name,int *Lmax,
			int *n_res,int *n_prot,k_opt *op)
{
 int i,j,eof;
 char amino[4];


 eof=fscanf(Fseq,"%s\n%d",T_name,Lseq);
 if(eof==EOF) return(0);
 for(i=0;i<*Lseq;++i)
   {
     fscanf(Fseq,"%s",amino);
     j=seq2dig(amino);
     if(j<0) read_err(T_name);     
     ++(*n_res);
   }

 ++(*n_prot);
 if(*Lseq>*Lmax) *Lmax=*Lseq;
 return(1);

}


/* read_Mod interprets model description as given in file Fmod */
int read_Mod(FILE *Fmod, k_opt *op, int *Imod, int *Alph, int *Adr)
{
 int i,j,eof;
 char amino[4],Type[5];
 float f;

 /* assumption - header was already read  and the file rewinded */
 eof=fscanf(Fmod,"%s %d %f",Type,&i,&f);
 if(eof==EOF) return(0);

 /* as reaching this place means that non-generic alphabet is going to
    be used - set the generic addresses to minus one */
 for(i=0;i<MAXS;++i)
   *(Adr+i)=-1;

 for(i=0;i<op->n_type;++i)
   {
    fscanf(Fmod,"%s%d",amino,Imod+i); 
    j=seq2dig(amino);
    if(j<0) return(0); /* if wrong ama in Model */
    else 
      {
	*(Alph+i)=j;
	*(Adr+j)=i;
      }
   }
 /* fix the GAP residue address for quick reference later on */
 /*  CHANGE if(op->gaps) op->gap_adr=*(Adr+seq2dig("GAP")); */
 op->gap_adr=*(Adr+seq2dig("GAP"));
 /*  CHANGE */
 return(1);

/* here is the example of Model file contents:
einm  
300
6.4
ALA 15
ARG 15
ASN 15
ASP 15
CYS 15
GLN 15
GLU 15
GLY 15
HIS 15
ILE 15
LEU 15
LYS 15
MET 15
PHE 15
PRO 15
SER 15
THR 15
TRP 15
TYR 15
VAL 15
*/

}

/* read_head_of_Mod reads type of potential and n_par as given in file Fmod */
/* it also used to read the rest (see read_Model) to redefine n_type */
int read_head_of_Mod(FILE *Fmod,k_opt *op)
{
 int i,eof,n,sum;
 char amino[4],Type[5];

 eof=fscanf(Fmod,"%s %d %f",Type,&op->n_par,&op->r_cut);
 if(eof==EOF) return(0);
 
 /* turn off the default value for einm */
 if(strcmp(Type,"eij")==0) 
   {
     op->eij=1;
     op->ein=0;
     op->einm=0;
     op->evdw=0;
   }
 if(strcmp(Type,"ein")==0) 
   {
     op->ein=1;
     op->einm=0;
     op->evdw=0;
     op->eij=0; 
   }
 if(strcmp(Type,"einm")==0) 
   {
     op->einm=1;
     op->ein=0;
     op->evdw=0;
     op->eij=0; 
   }
 if(strcmp(Type,"evdw")==0) 
   {
     op->evdw=1;
     op->ein=0;
     op->einm=0;
     op->eij=0; /* turn off the default value for eij */
   }

 /* now define how many types of letters (alphabet) are to be used */
 /* op->n_type=NTYPE; -  default if Model contains only header */
 n=0; sum=0;
 while(eof!=EOF)
   {    
     eof=fscanf(Fmod,"%s%d",Type,&i);
     if(eof!=EOF) 
       {
	 ++n;
	 sum+=i;
       }
   }
 /* the section of Model file containing description of letters and number
    of parameters for each letter may be skipped for eij, evdw - however to
    redefine alphabet one should add it (otherwise standard 20 letters are
    assumed */
 if(n!=0) 
   {
     op->n_type=n;
     if(op->ein || op->einm ) /* THOM */
       {
	 if(!(sum==op->n_par)) 
	   {
	     printf(" WARNING: inconsistency in Model \n");
	     printf(" n_par doesn't fit def. of alphabet - taking the latter\n");
	     op->n_par=sum;
	   }
       }
     if(op->eij) /* pairwise */
       {
	 if(op->n_type*(op->n_type+1)/2!=op->n_par)
	   {
	     printf(" WARNING: there is inconsistency in Model file! \n");
	     printf(" n_par doesn't fit def. of alphabet - taking the latter\n");
	     op->n_par=op->n_type*(op->n_type+1)/2;
	   }
       }
     if(op->evdw) /* pairwise - continuous */
       {
	 if(op->n_type*(op->n_type+1)!=op->n_par)
	   {
	     printf(" WARNING: there is inconsistency in Model file! \n");
	     printf(" n_par doesn't fit def. of alphabet - taking the latter\n");
	     op->n_par=op->n_type*(op->n_type+1);
	   }
       }
   }
 else  if(op->ein || op->einm) 
       {
	 printf(" You must specify number of types for each ama in Model ...\n");
	 exit(1);
       }
 
 /* now rewind file for the second attempt to read the remaining part */
 rewind(Fmod);
 
 return(1);

}

/* prepare_files opens all the files except for Model */
void prepare_files(a_file *F, k_opt *op, k_info *Info, float *Potn) 
{

 /* open the files */
 /* some names are fixed here and not in set_file_names ! */
 F->coor=fopen(F->f_coor,"r"); 
 F->seq=fopen(F->f_seq,"r");
 F->type=fopen("cont_by_Type","w");
 F->cont=fopen(F->f_cont,"r"); 
 /* right now op->distance is coupled with op->strucmp !!! */
 if(op->distance) F->rij=fopen(F->f_rij,"r"); 
 if(op->fngrps) F->fngrps=fopen(F->f_fps,"w");
 if(op->freeze && op->wDc) F->rhs=fopen(F->f_rhs,"w"); 
 if(op->evdw) F->vdw_a=fopen(F->f_vdw_a,"r"); 
 if(op->evdw) F->vdw_b=fopen(F->f_vdw_b,"r");

 /* switch from reading to writing CMs if necessary */ 
 if(op->evdw)
   {
     if(F->vdw_a==NULL || F->vdw_b==NULL || !op->r_cont)
       {
	 op->r_cont=0; /* recreate all the other CMs */
	 F->vdw_a=fopen(F->f_vdw_a,"w"); 
	 F->vdw_b=fopen(F->f_vdw_b,"w");
       }     
   }
 if(op->distance)
   {
     if(F->rij==NULL || !op->r_cont)
       {
	 op->r_cont=0; /* recreate all the other CMs */
	 F->rij=fopen(F->f_rij,"w");
       }
   }
 if(F->cont==NULL || !op->r_cont)
   {
     op->r_cont=0; /* redundant when !op->r_cont ... */
     F->cont=fopen(F->f_cont,"w");
   }

 /* open potential file or read potential instead of opening file */
 if(!op->def_pot)
   {
     F->pot=fopen(F->f_pot,"r");
     if(F->pot==NULL) open_err(F->f_pot);
   }
 else get_default_pot(F,op,Potn);;
 
 /* this gives precedence to current.pot */
 /*
 F->pot=fopen(F->f_pot,"r");
 if(F->pot==NULL) 
   {
     if(op->def_pot) get_default_pot(F,op,Potn);
     else open_err(F->f_pot);
   }
   */
 
 F->out1=fopen(F->f_out,"w");  
 if(op->wDc) F->out2=fopen(F->f_dcon,"w");
 if(op->mps) F->mps=fopen(F->f_mps,"w");
 if(op->info) F->info=fopen(F->f_info,"w");
 if(op->info2) F->iprot=fopen(F->f_iprot,"w");
 if(op->list) F->list=fopen(F->f_list,"r");
 if(op->build_db || op->rm_redund) F->newlist=fopen(F->f_newlist,"w");

 /* files for recognition */
 if(op->exam  || op->align || op->strucmp)
   {
     F->exam=fopen(F->f_exam,"r");
     F->best=fopen(F->f_best,"w");
     F->balg=fopen(F->f_balg,"w");
     F->query=fopen(F->f_query,"w"); /* switch to "r" after creating */
     F->qxyz=fopen(F->f_qxyz,"w"); /* switch to "r" after creating */
     F->qcont=fopen(F->f_qcont,"w"); /* switch to "r" after creating */
     F->qrij=fopen(F->f_qrij,"w"); /* switch to "r" after creating */
   }

 /* check the basic files */ 
 if(F->coor==NULL)
   {
     printf("Sorry - cannot find coordinate file: %s \n",F->f_coor);
     exit(1);
   }
 if(F->seq==NULL)
   {
     printf("Sorry - cannot find sequence file: %s \n",F->f_seq);
     exit(1);
   }
 if(F->cont==NULL)
   {
     printf("Cannot open ContMap - I give up \n");
     exit(1);
   }
 if(F->out1==NULL) open_err(F->f_out); 
 if(F->out2==NULL && op->wDc) open_err(F->f_dcon);
 if(F->rhs==NULL && op->wDc && op->freeze) open_err(F->f_rhs);
 if(F->mps==NULL && op->mps) open_err(F->f_mps);
 if(F->info==NULL && op->info) open_err(F->f_info);
 if(F->iprot==NULL && op->info2) open_err(F->f_iprot);
 if(F->list==NULL && op->list) open_err(F->f_list);
 if(F->newlist==NULL && (op->build_db || op->rm_redund)) 
   open_err(F->f_newlist);
 if(op->exam  || op->align || op->strucmp)
   {
     if(F->exam==NULL) open_err(F->f_exam);
     if(F->best==NULL) open_err(F->f_best);
     if(F->balg==NULL) open_err(F->f_balg);
     if(F->query==NULL) open_err(F->f_query); 
     if(F->qxyz==NULL) open_err(F->f_qxyz); 
     if(F->qcont==NULL) open_err(F->f_qcont); 
     if(F->qrij==NULL) open_err(F->f_qrij); 
   }
 
}

/* set the default file names */
void set_file_names(a_file *F, k_opt *op, k_info *Info)
{

 strcpy(F->f_dcon,"current.dcon");
 strcpy(F->f_rhs,"current.rhs");
 strcpy(F->f_mps,"current.mps");
 strcpy(F->f_pot,"current.pot");
 strcpy(F->f_info,"current.info");
 strcpy(F->f_iprot,"proteins.info");
 strcpy(F->f_exam,"seq_to_examine.txt");
 strcpy(F->f_list,"list.txt");
 strcpy(F->f_newlist,"new_list.txt");
 strcpy(F->f_query,"Query");
 strcpy(F->f_qxyz,"Query_xyz");
 strcpy(F->f_qrij,"Query_rij");
 strcpy(F->f_qcont,"Query_CM");
 strcpy(F->f_out,"xscan.log");
 strcpy(F->f_best,"best.log");
 strcpy(F->f_balg,"alignments.log");
 strcpy(F->f_coor,"XYZ");
 strcpy(F->f_seq,"SEQ");
 strcpy(F->f_mod,"Model");
 strcpy(F->f_cont,"ContMap");
 strcpy(F->f_vdw_a,"ContMap_vdw_A");	 
 strcpy(F->f_vdw_b,"ContMap_vdw_B");	 
 strcpy(F->f_rij,"ContMap_rij");
 strcpy(F->f_fps,"FINGERPS");
  
  
}

/* set the default options */
void set_default_opt(a_file *F, k_opt *op, k_info *Info)
{

 /* set up THOM2 as the default model of the potential  */
 op->eij=0; 
 op->ein=0; 
 op->einm=1; 
 op->evdw=0; 

 /* set up the default parameters for the model of the potential */

 /* the default number of letters */
 op->n_type=NTYPE; /* give explicit def. of an alph. in Model to change it */ 
 ++op->n_type; /* now GAP is default letter */

 /* the default number of paramts - consistent with the default eij,n_type */
 /* op->n_par=NTYPE*(NTYPE+1)/2; */
 op->n_par=op->n_type*NTHOM2+1; /* now standard THOM2 is the default */
 /* adjusted in Model for the actual n_type */

 /* default value of the cutoff */
 op->r_cut=R_high; 

 /* set up the default mode of operation as gapless threading for LP learning*/
 op->crd2xyz=0; 
 op->pdb2xyz=0; 
 op->fngrps=0;
 op->learn=0; 
 op->decoy=0;
 op->exam=1;
 op->strucmp=0;
 op->gaps=1; /*  use gaps in seq-stru alignment */
 op->gapless=1; /* this is not opposite to the one above */
 op->align=0; /* do not assume seq-seq alignment */
 op->debug=0;

 /* set up the defaults for all the utilities */

 /* do not print extensive output */
 op->info=0; 
 op->info2=0; 
 op->info3=0; 
 op->info4=0;
 op->ver=0;
 
 /* try to read existing ContMap first */
 op->r_cont=1;
 op->distance=0; /* do not employ the distance map */
 
 /* do not assume matrix format of pot_file */
 op->mat=0; 

 /* pot is (not scaled) "energy" and not scoring function */
 op->resc_pot=1.0; 

 /* Lennard-Jones rep_pow-atr_pow model */
 op->rep_pow=6; 
 op->atr_pow=2; 

 /* default range of threading loop */
 op->beg=1; 
 op->end=0; /* this will be redefined in alloc_mem */ 

 /* default format of XYZ */
 op->col=2; 
 op->n_col=2; 

 /* do not print LP constraints  */
 op->wDc=0; 
 op->mps=0; 

 /* histogram length */
 op->max_histo=2*NTYPE; 

 /* do not use ene_threshold */
 op->th_ene=0.0; 
 op->noth=1; 
 op->up=0;  

 /* do not attempt building or re-building database */
 op->rm_redund=0;
 op->build_db=0;
 /* do not use one-letter ama codes */
 op->onelet=0; 

 /* default gap penalties */
 op->gap_pen=GAP_PNLT; 
 op->pref_pen=GAP_PNLT; 
 op->cdel=0;
 op->strgap=1; /* env dependent gap penalty for seq-seq align */
 
 op->loc_a=0; /* do not perform local align but rather global one */
 op->t_length=0.50; /* skip alignment if its length is less than 50% */
 op->t_zscore=1.5; /* skip alignment if its Zscore is less than 1.5 */
 op->t_rms=12.0; /* skip alignment if its RMS is larger than 12.0 */
 op->t_rms_db=3.0; /* exclude from database if RMS is smaller than 3.0 */
 op->t_seqid=0.50; /* skip structure if seqid larger than 50% */
 op->search_depth=1; /* default search depth */
 op->justq=0; /* do not stop after building query files */
 
 /* assume that there is no GAP res in current alph */
 op->gap_adr=-1; 

 /* do not freeze any variables */
 op->freeze=0;
 op->nfrozen=0;

 /* by default take 20 best alignments and make 100 shufflings for distrib. */
 op->kbest=20;
 op->k_prn=10;
 op->n_shuffle=100;
 strcpy(op->qname,"query");
 op->qfmt=1; /* default format is that of loopp */
 op->list=0; /* query file is not a list of file names */
 op->prn_shuffle=0; /* do not print alignments for shuffled seqs */
 op->mod_search_depth=0; /* check later whether defaults were altered */
 
 op->pick=0; /* do not pick seq or str by name */
 op->pick_all_seq=0; /* do not pick all if pick */
 strcpy(op->pseq_name,"");
 strcpy(op->pstr_name,"");
 op->def_pot=1; /* do not take potential from a file */

 
} /* end of dafaults */  

/* interpret the input line - at most two parameters for each option */
void command_line(int argc, char *argv[], a_file *F, k_opt *op, k_info *Info) 
{
 int plus_1,plus_2,i,itmp;
 char stmp[MAXS];
  
 for(i=0;i<argc;++i)
   {
     plus_1=0;
     plus_2=0;
     if(i+1<argc)
       {
	 if(strncmp(argv[i+1],"-",1)!=0) plus_1=1;
       }
     if(i+2<argc)
       {
	 if( plus_1 && (strncmp(argv[i+2],"-",1)!=0) ) plus_2=1;
       }
     if(strncmp("-h",argv[i],2)==0 ) print_help(argv,Info);
     if(strcmp(argv[i],"-b")==0) 
       {
	 if(plus_1) op->beg=atoi(argv[i+1]);
	 if(plus_2) op->end=atoi(argv[i+2]);
	 if(op->beg<1) op->beg=1;
	 if(op->end<1) op->end=0;
	 else if(op->end<op->beg)
	   {
	     itmp=op->beg;
	     op->beg=op->end;
	     op->end=op->beg;
	   }
       }
     if(strcmp(argv[i],"-eij")==0) 
       {
	 op->eij=1;
	 op->ein=0;
	 op->einm=0;
	 op->evdw=0;
       }
	 if(strcmp(argv[i],"-evdw")==0)
	 {
		 op->evdw=1;
		 op->eij=0;
		 op->einm=0;
		 op->ein=0;
	 }
     if(strcmp(argv[i],"-ein")==0) 
       {
	 op->ein=1;
	 op->eij=0;
	 op->einm=0;
	 op->evdw=0;
       }
     if(strcmp(argv[i],"-einm")==0) 
       {
	 op->einm=1;
	 op->eij=0;
	 op->ein=0;
	 op->evdw=0;
       }
     if(strncmp(argv[i],"-mod",3)==0) 
       {
	 if(plus_1) strcpy(F->f_mod,argv[i+1]);
       }
     
     if(strcmp(argv[i],"-s")==0) /* seq alignm */
       {
	 if(plus_1) strcpy(F->f_exam,argv[i+1]);
	 op->align=1;
	 op->learn=0;
	 op->exam=0;
	 op->decoy=0;
	 op->strucmp=0;
       }
     if(strcmp(argv[i],"-q")==0) /* inequalities */
       {
	 op->learn=1;
	 op->exam=0;
	 op->align=0;
	 op->decoy=0;
	 op->strucmp=0;
       }
     if(strcmp(argv[i],"-x")==0) /* stru to stru alignm */
       {
	 if(plus_1) strcpy(F->f_exam,argv[i+1]);
	 op->strucmp=1;
	 op->distance=1;
	 op->exam=0;
	 op->learn=0;
	 op->align=0;
	 op->decoy=0;
       }
     if(strcmp(argv[i],"-p")==0) 
       {
	 if(plus_1) strcpy(F->f_pot,argv[i+1]);
	 op->def_pot=0; /* to read potential from a file */
       }
     if(strcmp(argv[i],"-xyz")==0) 
       {
	 if(plus_1) strcpy(F->f_coor,argv[i+1]);

       }
     if(strcmp(argv[i],"-seq")==0) 
       {
	 if(plus_1) strcpy(F->f_seq,argv[i+1]);

       }
     if(strcmp(argv[i],"-cm")==0) 
       {
	 if(plus_1) 
	   {
	     strcpy(F->f_cont,argv[i+1]);
	     strcpy(stmp,argv[i+1]);
	     strcpy(F->f_rij,strcat(stmp,"_rij"));
	     strcpy(stmp,argv[i+1]);
	     strcpy(F->f_vdw_a,strcat(stmp,"_vdw_A"));	
	     strcpy(stmp,argv[i+1]);
	     strcpy(F->f_vdw_b,strcat(stmp,"_vdw_B"));	     
	   }
       }
     if(strcmp(argv[i],"-o")==0) 
       {
	 if(plus_1) strcpy(F->f_out,argv[i+1]);
	 if(plus_2) strcpy(F->f_info,argv[i+2]);

       }
     if(strcmp(argv[i],"-log")==0) 
       {
	 if(plus_1) strcpy(F->f_best,argv[i+1]);
	 if(plus_2) strcpy(F->f_balg,argv[i+2]);

       }
     if(strcmp(argv[i],"-w")==0) 
       {
	 if(plus_1) strcpy(F->f_dcon,argv[i+1]);
	 op->wDc=1;
       }
     if(strcmp(argv[i],"-frz")==0) 
       {
	 if(plus_1) op->nfrozen=atoi(argv[i+1]);
	 op->freeze=1;
       } 
     if(strcmp(argv[i],"-best")==0) 
       {
	 if(plus_1) op->kbest=atoi(argv[i+1]);
	 if(plus_2) op->n_shuffle=atoi(argv[i+2]);
	 if(op->kbest<1) op->kbest=1;
	 if(op->n_shuffle<2) op->n_shuffle=2;
	 op->mod_search_depth=1;
       } 
     if(strcmp(argv[i],"-nprn")==0) 
       {
	 if(plus_1) op->k_prn=atoi(argv[i+1]);
	 if(op->k_prn<1) op->k_prn=1;
       } 
     if(strcmp(argv[i],"-pick")==0) 
       {
	 op->pick=1; /* pick seq and/or structure by name */
	 if(plus_1) strcpy(op->pseq_name,argv[i+1]); 
	 if(plus_2) strcpy(op->pstr_name,argv[i+2]);
       } 
     if(strcmp(argv[i],"-prnshf")==0) 
       {
	 op->prn_shuffle=1;
       }
     if(strcmp(argv[i],"-name")==0) 
       {
	 if(plus_1) strcpy(op->qname,argv[i+1]);

       }
     if(strcmp(argv[i],"-fmt")==0) 
       {
	 if(plus_1) 
	   {
	     if(strncmp("loopp",argv[i+1],1)==0) op->qfmt=1; /* loopp format */
	     if(strncmp("fasta",argv[i+1],1)==0) op->qfmt=2; /* >name  1let */
	     if(strncmp("swiss",argv[i+1],1)==0) op->qfmt=3; /* SwissProt */
	     if(strncmp("pdb",argv[i+1],1)==0) op->qfmt=4; /* PDB SEQRES */
	     if(strncmp("1let",argv[i+1],1)==0) op->qfmt=5; /* 1let */
	     if(strncmp("plain",argv[i+1],2)==0) op->qfmt=5; /* 1let */
	     if(strncmp("crd",argv[i+1],1)==0) op->qfmt=6; /* CHARMM */
	   }
       }
     if(strcmp(argv[i],"-mps")==0) 
       {
	 if(plus_1) strcpy(F->f_mps,argv[i+1]);
	 op->mps=1;
       }  
     if(strcmp(argv[i],"-te")==0) 
       {
	 op->th_ene=0.0;
	 if(plus_1) op->th_ene=atof(argv[i+1]);
	 op->noth=0; 
	 if(plus_2) op->th_ene=-op->th_ene;
       }  
     if(strcmp(argv[i],"-tl")==0) 
       {
	 if(plus_1) op->t_length=atof(argv[i+1])/100.0;
	 if(op->t_length<0.01) op->t_length=0.01;
	 if(op->t_length>1.0) op->t_length=1.0;
       }  
     if(strcmp(argv[i],"-tsid")==0) 
       {
	 if(plus_1) op->t_seqid=atof(argv[i+1])/100.0;
	 if(op->t_seqid<0.01) op->t_seqid=0.01;
	 if(op->t_seqid>1.0) op->t_seqid=1.0;
       }  
     if(strcmp(argv[i],"-tz")==0) 
       {
	 if(plus_1) op->t_zscore=atof(argv[i+1]);
	 if(plus_2) op->t_zscore=-op->t_zscore;
       } 
     if(strcmp(argv[i],"-tr")==0) 
       {
	 if(plus_1) op->t_rms=atof(argv[i+1]);
       } 
     if(strcmp(argv[i],"-trmsDB")==0 || strcmp(argv[i],"-trDB")==0) 
       {
	 if(plus_1) op->t_rms_db=atof(argv[i+1]);
       }
     if(strcmp(argv[i],"-justQ")==0 || strcmp(argv[i],"-justq")==0) 
       {
	 op->justq=1;
       }  
     if(strcmp(argv[i],"-u")==0) 
       {
	 op->up=1;
       }
     if(strcmp(argv[i],"-i")==0) 
       {
	 /* initially op->info=0 */
	 if(plus_1) op->info=atoi(argv[i+1]);
	 if(plus_2) op->max_histo=atoi(argv[i+2]);
	 if(op->info>1) op->info2=1; /* additional info when -i 2 */
	 if(op->info>2) op->info3=1; /* additional info when -i 3 */
	 if(op->info>3) op->info4=1; /* additional info when -i 4 */
	 op->info=1; /* basic info always when option -i */
       }
     if(strncmp("-v",argv[i],2)==0) 
       {
	 op->ver=1;
       }
     if(strcmp(argv[i],"-m")==0) 
       {
	 op->mat=1;
       }
     if(strcmp(argv[i],"-newCM")==0) 
       {
	 op->r_cont=0; /* regenerate ContMap anyway */
       }
     if(strcmp(argv[i],"-g")==0) 
       {
	 if(plus_1) op->pref_pen=atof(argv[i+1]);
	 if(plus_2) op->gap_pen=atof(argv[i+2]);
	 op->gaps=1;
       }
     if(strcmp(argv[i],"-nog")==0) 
       {
	 op->gaps=0;
       }
     if(strcmp(argv[i],"-1let")==0) 
       {
	 op->onelet=1;
       }
     if(strcmp(argv[i],"-l")==0) 
       {
	 op->loc_a=1;
       }
     if(strcmp(argv[i],"-sd")==0) 
       {
	 /* initially op->search_depth=1 */
	 if(plus_1) op->search_depth=atoi(argv[i+1]);
       }
     if(strcmp(argv[i],"-cdel")==0) 
       {
	 op->cdel=1;
       }
     if(strcmp(argv[i],"-cgap")==0) 
       {
	 op->strgap=0;
       }
     if(strcmp(argv[i],"-rsc")==0) 
       {
	 if(plus_1) op->resc_pot=atof(argv[i+1]);
	 if(plus_2) op->resc_pot=-op->resc_pot;
       }
     if(strcmp(argv[i],"-fps")==0) 
       {
	 if(plus_1) strcpy(F->f_fps,argv[i+1]);
	 op->fngrps=1;
       }
     if(strcmp(argv[i],"-k")==0) 
       {
	 if(plus_1) op->col=atoi(argv[i+1]);
	 if(plus_2) op->n_col=atoi(argv[i+2]);
	 if(op->col>3) op->col=2; /* go back to the default value */
	 if(op->n_col>3) op->n_col=2; /* go back to the default value */
	 if(op->col>op->n_col) op->col=op->n_col;
       }
     if(strcmp(argv[i],"-t")==0) /* threading */
       {
	 if(plus_1) strcpy(F->f_exam,argv[i+1]);
	 op->exam=1;
	 op->align=0;
	 op->learn=0;
	 op->decoy=0;
	 op->strucmp=0;
       }
     if(strcmp(argv[i],"-d")==0) /* decoys */
       {
	 op->decoy=1;
	 op->exam=0;
	 op->align=0;
	 op->learn=0;
	 op->strucmp=0;
       }
     if(strcmp(argv[i],"-pvdw")==0) 
       {
	 if(plus_1) op->rep_pow=atoi(argv[i+1]);
	 if(plus_2) op->atr_pow=atoi(argv[i+2]);
       }
     if(strcmp(argv[i],"-dbg")==0) 
       {
	 op->debug=1;
       }
     if(strcmp(argv[i],"-list")==0) 
       {
	 op->list=1;
	 if(plus_1) strcpy(F->f_list,argv[i+1]);
       }
     if(strcmp(argv[i],"-nored")==0) 
       {
	 op->rm_redund=1;
	 op->build_db=0;
	 if(plus_1) strcpy(F->f_newlist,argv[i+1]);
       }
     if(strcmp(argv[i],"-newDB")==0) 
       {
	 op->build_db=1;
	 op->rm_redund=0;
	 if(plus_1) strcpy(F->f_newlist,argv[i+1]);
       }

   } /* end of loop over input arguments */
  
} /* end of command_line */

/* interpret Model file to adjust options */
void interpret_Model(a_file *F, k_opt *op, k_info *Info)
{

 /* set the type of the potential and n_par, n_type - read it from file Model */
 /* assuming standard pair potential and n_par=210 if Model not found */
 /* not anymore - einm (THOM2) is default now */

 /* read the header only to check type of the potential here - 
    the rest will be interpreted in set_Model */
 F->mod=fopen(F->f_mod,"r"); 
 if(!(F->mod==NULL)) 
   {
     if(read_head_of_Mod(F->mod,op))
       {
	 if(op->info2)
	   {
	     sprintf(Info->Sout,
		     " Reading type of potential and n_par from file: Model \n");
	     to_stnd_out(Info->Sout);
	   }
       }
     else
       {
	 sprintf(Info->Sout," Problem while reading Model \n");
	 to_stnd_out(Info->Sout);
	 exit(1);
       }
   }
 else
   {
     if(op->align) 
       {
	 op->eij=1;
	 op->einm=0;
	 op->ein=0;
	 op->evdw=0;
       }
     
     if(op->strucmp) 
       {
	 op->eij=0;
	 op->einm=1;
	 op->ein=0;
	 op->evdw=0;
       }
     
     if(op->eij)
       {
	 op->n_par=NTYPE*(NTYPE+1)/2;
	 op->n_type=NTYPE;
       }
     if(op->ein)
       {
	 op->n_par=op->n_type*NTHOM1;
       }
     if(op->einm)
       {
	 op->n_par=op->n_type*NTHOM2+1; /* plus one for GAP */
	 if(op->strucmp) 
	   {
	     op->n_type=1; /* one-letter (ALL) alph. */
	     op->n_par=op->n_type*NTHOM2;
	   }
       }
     
   }
  
}

/* check consistency and adjust options if necessary */
void check_opt(a_file *F, k_opt *op, k_info *Info) 
{

 /* adjust the depth of search for local (threading) alignments */
 if(op->loc_a && op->exam && !op->mod_search_depth)
   {
     op->kbest=250;
     op->n_shuffle=100;
   }

 /* adjust the depth of search for local sequence alignments */
 if(op->loc_a && op->align && !op->mod_search_depth)
   {
     op->kbest=50;
     op->n_shuffle=100;
   }

 /* adjust the depth of search for stru-stru alignments */
 if(op->strucmp && !op->mod_search_depth)
   {
     op->kbest=100;
   }

 /* adjust gap penalts for pairwise models (frozen environment approximation */
 if(op->eij && op->exam && op->gaps)
   {
     op->gap_pen=GAP_PNLT_EIJ; 
     op->pref_pen=GAP_PNLT_EIJ;  
   }

 /* pick all the query sequences using "all" */
 if(op->pick)
   {
     if(strncmp(op->pseq_name,"all",3)==0) op->pick_all_seq=1;
   }

 /* noredundant only for LOOPP format - it assumes SEQ, XYZ files */
 if(op->rm_redund)
   {  
     if(op->qfmt!=1) op->rm_redund=0;
   }

 /* f_exam is irrelevant in case of list, so just open the same file */
 if(op->list) 
   {
     strcpy(F->f_exam,F->f_seq); 
     /* strcpy(F->f_exam,F->f_list); */
   }

 /* if loopp format and stru2stru align then Query=SEQ */
 if(op->strucmp && op->qfmt==1)
   {
     strcpy(F->f_exam,F->f_seq); 
   }
 
}

/* set_Model reads Model file again to setup alphabet that is used */
void set_Model(a_file *F, k_opt *op, k_info *Info, int *Imod,
	       int *Alph, int *Adr, float *Potn)
{
  int i,j,k,sum;

  /* initiate addresses of generic symbols before you go to check the
     actual alphabet in Model */
  for(i=0;i<MAXS;++i) *(Adr+i)=i;

  /* initiate generic alphabet to work with it when Model not found */
  for(i=0;i<op->n_type;++i)
    {
      *(Alph+i)=i;
      if(op->strucmp) /* set Model for stru-stru align */
	{
	  /* assumption: alphabet shorter than MAXS */
	  for(k=0;k<MAXS;++k) *(Adr+k)=-1; 
	  j=seq2dig("ALL");
	  *(Alph+i)=j;
	  *(Adr+j)=i;
	}
      if(op->ein) *(Imod+i)=NTHOM1;
      if(op->einm) *(Imod+i)=NTHOM2;  
    }
  /* increase Imod for GAP to accomodate penalty for site of no neighbors */
  if(op->einm && !op->strucmp) ++(*(Imod+op->n_type-1)); 
  op->gap_adr=*(Adr+seq2dig("GAP"));
 
  /* define the model of the potential - read it from file Model */
  /* assuming THOM2 potential if Model not found */ 
  if(!(F->mod==NULL)) 
    {
      if(read_Mod(F->mod,op,Imod,Alph,Adr))
	{
	  if(op->info2) 
	    {
	      sprintf(Info->Sout," Reading Model finished \n");
	      to_stnd_out(Info->Sout);
	    }
	}
      else
	{
	  sprintf(Info->Sout," Problem while reading Model \n");
	  to_stnd_out(Info->Sout);
	  exit(1);
	}
    }

  /* get potential now */
  if(!op->def_pot) read_pot(op,Potn,F,Imod); 

} /* end of set_Model */


/* set_options interprets the input line to pick up right options */
void set_options(int argc, char *argv[], a_file *F, k_opt *op, k_info *Info) 
{
  
 /* default values of input file names - use appropriate option to change it */
 set_file_names(F,op,Info);
  
 /* default values of input parameters - use appropriate option to change it */
 set_default_opt(F,op,Info);

 /* loop over input line - at most two parameters for each option */
 command_line(argc,argv,F,op,Info);
 
 /* model of the potential may be overwritten using input file Model */
 interpret_Model(F,op,Info);

 /* adjust and check consistency of options */
 check_opt(F,op,Info);

 
} /* end of set_options */


