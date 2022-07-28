/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000 v. 2.000                                      */
/*                                                                */
/*        PRINT STANDARD OUTPUT AND OTHER USEFUL INFORMATIONS     */
 
#include "loopp.h"
   

void to_stnd_out(char *Str_to_out)
{
  
  printf("%s",Str_to_out);
  
}

void to_file(FILE *Fout,char *Str_to_out)
{
  
  fprintf(Fout,"%s",Str_to_out);
  
}

void exit_end(char *Str_to_out,a_file *F,k_prot *prot,k_opt *op,k_info *Info,
	      k_Model *model,solv_sh *ssh,k_dynpt *dpt)
{
  printf("%s",Str_to_out);
  clean_up(F,prot,op,Info,model,ssh,dpt);
  exit(0);
}

void exit_err(char *Str_to_out,a_file *F,k_prot *prot,k_opt *op,k_info *Info,
	      k_Model *model,solv_sh *ssh,k_dynpt *dpt)
{
  printf("%s",Str_to_out);
  clean_up(F,prot,op,Info,model,ssh,dpt);
  exit(1);
}
  
/* print_res writes down result of a gapless threading */
void prn_res(k_opt *op, a_file *F, float *Potn, int *contR, int *contT,
	     unsigned long int *index, float *T_ene, char *R_name,
	     char *T_name, k_info *Info, int *n, int *Ls,
	     float *contRvdw, float *contTvdw, k_mps *var_mps[],
	     k_mps *tail_mps[])
{
  
  if(op->noth) /* if there is no threshold for energy */
    {
      fprintf(F->out1,Info->fmt,*index,*T_ene,R_name,T_name,
	      *n+1,*n+*Ls);
      ++Info->n_prn;
      if(op->wDc || op->mps) 
	write_ineq(F,op,Potn,contR,contT,&Info->n_prn,
		   contRvdw,contTvdw,&var_mps[0],&tail_mps[0]); 
    }
  else
    {
      if(op->up) 
	{
	  if(*T_ene>op->th_ene)
	    {	      
	      fprintf(F->out1,Info->fmt,*index,*T_ene,R_name,T_name,
		      *n+1,*n+*Ls);
	      ++Info->n_prn;
	      if(op->wDc || op->mps) 
		write_ineq(F,op,Potn,contR,contT,&Info->n_prn,
			   contRvdw,contTvdw,&var_mps[0],&tail_mps[0]); 
	    }
	}
      else if(*T_ene<=op->th_ene)
	{
	  fprintf(F->out1,Info->fmt,*index,*T_ene,R_name,T_name,
		  *n+1,*n+*Ls);
	  ++Info->n_prn;
	  if(op->wDc || op->mps) 
	    write_ineq(F,op,Potn,contR,contT,&Info->n_prn,
		       contRvdw,contTvdw,&var_mps[0],&tail_mps[0]); 
	}
      
    }  
}

/* to initiate statistics and histograms coputation */
void init_info(k_opt *op,a_file *F,k_info *Info,
	       k_mps *var_mps[],k_mps *tail_mps[])
{
  int i;
  
  strcpy(Info->fmt,"%lu\t%f\t%s %s:%d %d\n"); /* format for f_out */
  Info->av_ene=0.0;
  Info->z_score=0.0;
  for(i=0;i<op->max_histo;++i) *(Info->histo+i)=0;  
  for(i=0;i<op->n_par;++i) 
    {
      *(Info->con_all+i)=0;
      *(Info->con_All_decoys+i)=0.0;
      *(Info->con_all_dev+i)=0;
    }
  Info->n_prn=0;
  Info->ncont=0;
  Info->nres=0;
  for(i=0;i<op->n_type;++i) *(Info->res_Stat+i)=0;
  Info->ncont_decoys=0.0;
  Info->ntop=0;

  /* add small info to current.info file */
  if(op->info2) 
    {
      sprintf(Info->Sout,
	      "Contrib. to native energies - training set %s and potential %s\n",
	      F->f_coor,F->f_pot); 
      to_file(F->iprot,Info->Sout);
    }

  if(op->mps)
    {
      /* the pointer arrays *_mps are not dynamically allocated */
      if(op->n_par>MAX) 
	{
	  sprintf(Info->Sout,"Increase MAX to %d+1 and recompile \n",
		  op->n_par);
	  to_stnd_out(Info->Sout);
	  exit(0);
	}
 
      for(i=0;i<op->n_par;++i)
	{	
	  *(var_mps+i)=NULL;
	  *(tail_mps+i)=NULL;
	}
    }
  
  
}

/* info computes moments and histograms of the energy diff. distribution */
void prep_info(k_info *Info,float *T_ene,int *max_histo)
{
  Info->av_ene+=*T_ene; 
  Info->z_score+=(*T_ene)*(*T_ene);
  if(*T_ene>*max_histo-1) *T_ene=*max_histo-1;
  if(*T_ene<0.0) *T_ene=0.0;
  ++*(Info->histo+(int)(*T_ene));  
}

/* DB amino acid residues statistics */
void add_res(k_info *Info, int *Seq, int *l_seq, int *Adr, int *Alph)
{

 int i,j;

 for(i=0;i<*l_seq;++i)
   {
     j=*(Alph + *(Seq+i));
     ++(*(Info->res_Stat+*(Adr+j)));
     ++Info->nres;
   }
 
}

/* DB contact statistics */
void add_cont(int *n_para, int *contR, unsigned long int *n_cont,
	      unsigned long int *cont_all, unsigned long int *cont_all2)
{
 int i;

 for(i=0;i<*n_para;++i)
   {
     *(cont_all+i)+=*(contR+i);
     *(cont_all2+i)+=*(contR+i)*(*(contR+i));
     *n_cont+=*(contR+i);
   }
  
}

/* decoys contact statistics */
void add_cont_decoy(int *n_para,int *contT,double *x_cont,double *cont_All)
{
 int i;
 double eps_loc;
 

 /* put 1e-6 for each contact to avoid huge numbers */
 eps_loc=0.000001;
 
 for(i=0;i<*n_para;++i)
   {
     *(cont_All+i)+=(*(contT+i)*eps_loc);
     /* *x_cont+=(*(contT+i)*eps_loc); */ /* moved to info */
   }
  
}

void prn_title_bar(FILE *Finfo, char *Title, k_opt *op, int *Alph)
{ 
  int i;
  
  to_file(Finfo,Title);
  /* now use Title as a buffer */
  for(i=0;i<op->n_type;++i)
    {
      sprintf(Title,"   %s",dig2amino(*(Alph+i)));
      to_file(Finfo,Title);
    }
  sprintf(Title,"\n");
  to_file(Finfo,Title);
}


void info(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, int *Imod,
	  int *Alph, int *Adr)
{
 int i,j,k,l,m,n;
 char stmp[MAXS];
 float diff,av_diff,tot_diff,av_cont,ftmp,f1,f2;
 double norm,dtmp;
 

 /* small piece of info */
 sprintf(Info->Sout," Total number of residues in the training set: %ld \n",
	 Info->nres);
 to_file(Finfo,Info->Sout);
 sprintf(Info->Sout," Residue frequencies (percent of total #res): \n");
 prn_title_bar(Finfo,Info->Sout,op,Alph);
 for(i=0;i<op->n_type;++i) 
   {
     j=*(Adr + *(Alph+i));
     sprintf(Info->Sout," %5.2f",*(Info->res_Stat+j)*100.0/(float)Info->nres);
     to_file(Finfo,Info->Sout);
   }
 sprintf(Info->Sout,"\n\n");
 to_file(Finfo,Info->Sout);

 sprintf(Info->Sout,
	 " Total number of contacts (or sites) in the training set: %ld \n",
	 Info->ncont); 
 to_file(Finfo,Info->Sout);
 sprintf(Info->Sout,
	 " Contact (or site) type statistics (native structures):\n");
 prn_title_bar(Finfo,Info->Sout,op,Alph);

 /* first print the native contacts statistics */
 info_cont_nat(Finfo,Info,Potn,op,Imod,Alph,Adr);

 /* now print statistics for decoys */
 if(op->learn || op->decoy) info_cont_dec(Finfo,Info,Potn,op,Imod,Alph,Adr);

 /* print potential used for threading or alignment */ 
 info_pot(Finfo,Info,Potn,op,Imod,Alph,Adr);

 /* better TOM2 description in terms of effective pair contributions */
 if(op->einm && op->info3) info_effpair_einm(Finfo,Info,Potn,op,Imod,Alph,Adr);

 /* total stabilization (over training set) decomposed into parameter contr */
 info_ene_tot(Finfo,Info,Potn,op,Imod,Alph,Adr);

 /* start threading info */
 /* compute moments and histograms of the energy diff. distribution */
 sprintf(Info->Sout," total number of alignments = %lu \n",Info->index);
 to_file(Finfo,Info->Sout);

 Info->av_ene=Info->av_ene/(double)(Info->index);
 strcpy(stmp,"");
 if(op->learn || op->decoy) strcpy(stmp,"(with respect to the native one)");
 sprintf(Info->Sout," average energy %s = %5.2f \n",stmp,Info->av_ene);
 to_file(Finfo,Info->Sout);

 Info->z_score=Info->z_score/(double)(Info->index);
 dtmp=Info->z_score-(Info->av_ene)*(Info->av_ene);
 if(dtmp<epsilon) dtmp=epsilon;
 Info->z_score=(Info->av_ene)/sqrt(dtmp); 
 sprintf(Info->Sout," aver over sigma ratio = %5.2f \n",Info->z_score);
 to_file(Finfo,Info->Sout);
 
 sprintf(Info->Sout," Histogram (energy gaps or energies) ... \n");
 to_file(Finfo,Info->Sout);
 for(i=0;i<op->max_histo;++i) 
   {
     sprintf(Info->Sout,"%d\t%ld\n",i,*(Info->histo+i));
     to_file(Finfo,Info->Sout);
   }
 
} 

/* get_stat_poten prints statistical potential resulting from the current DB */
void get_stat_poten(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, 
		    int *Imod, int *Alph, int *Adr)
{
 int i,j,k,l,m,n;
 float f1,f2;

 sprintf(Info->Sout,
	 " Statisical potential resulting from the set of templates [kT]\n");
 prn_title_bar(Finfo,Info->Sout,op,Alph);

 for(i=0;i<op->n_type;++i)
   {
     m=*(Adr + *(Alph+i));
     f1=(float)Info->nres/(float)*(Info->res_Stat+m);
     for(j=0;j<op->n_type;++j) 
       {
	 n=*(Adr + *(Alph+j));
	 f2=(float)Info->nres/(float)*(Info->res_Stat+n);
	 if(i>j) k=mat2vec(&j,&i,&op->n_type);
	 else k=mat2vec(&i,&j,&op->n_type);
	 fprintf(Finfo," %5.2f",
		 -log(*(Info->con_all+k)*f1*f2/(float)Info->ncont));
       }
     fprintf(Finfo,"\n");
   } 
 
}

/* info_cont_nat prints statistics regarding native contacts */
void info_cont_nat(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, 
		   int *Imod, int *Alph, int *Adr)
{
 int i,j,k,l,m,n;
 unsigned long int *cont_all,*cont_all_dev;
 float diff,av_diff,tot_diff,av_cont,ftmp;
 
 cont_all=Info->con_all;
 cont_all_dev=Info->con_all_dev;

 if(op->eij) /* first pairwise models */
   {
     for(i=0;i<op->n_type;++i)
       {
	 for(j=0;j<op->n_type;++j) 
	   {
	     if(i>j) k=mat2vec(&j,&i,&op->n_type);
	     else k=mat2vec(&i,&j,&op->n_type);
	     fprintf(Finfo," %5ld",*(cont_all+k));
	   }
	 fprintf(Finfo,"\n");
       }    

     /* now print the same in terms of percentage of all native contacts */
     sprintf(Info->Sout,
	     " Native contacts (or sites) as percent of total #cont=%ld\n",
	     Info->ncont);
     prn_title_bar(Finfo,Info->Sout,op,Alph);

     for(i=0;i<op->n_type;++i)
       {
	 for(j=0;j<op->n_type;++j) 
	   {
	     if(i>j) k=mat2vec(&j,&i,&op->n_type);
	     else k=mat2vec(&i,&j,&op->n_type);
	     fprintf(Finfo," %5.2f",
		     *(Info->con_all+k)*100.0/(float)Info->ncont);
	   }
	 fprintf(Finfo,"\n");
       } 
     
     /* print also statistical potential resulting from the current DB */
     get_stat_poten(Finfo,Info,Potn,op,Imod,Alph,Adr);

   } /* end of pairwise for native structures */

 if(op->ein || op->einm) /* now onion models */
   {
     k=*(Imod); /* take the number of parameters for first ama */
     for(i=1;i<op->n_type;++i)
       {
	 if(*(Imod+i)>k) k=*(Imod+i); /* take the maximum num of para */
       }
     for(i=0;i<k;++i)
       {
	 l=0;
	 for(j=0;j<op->n_type;++j) 
	   {
	     l+=*(Imod+j);
	     m=l-*(Imod+j);
	     if(*(Imod+j)<i+1) fprintf(Finfo," %s","     ");
	     else fprintf(Finfo," %5ld",*(cont_all+m+i));
	   }
	 fprintf(Finfo,"\n");
       }

     /* now print the same in terms of percentage of all native contacts */
     fprintf(Finfo,
	     " Native contacts (or sites) as percent of total #cont=%ld\n",
	     Info->ncont);
     for(i=0;i<op->n_type;++i) fprintf(Finfo,"   %s",dig2amino(*(Alph+i)));
     fprintf(Finfo,"\n"); 
     for(i=0;i<k;++i)
       {
	 l=0;
	 for(j=0;j<op->n_type;++j) 
	   {
	     l+=*(Imod+j);
	     m=l-*(Imod+j);
	     if(*(Imod+j)<i+1) fprintf(Finfo," %s","     ");
	     else fprintf(Finfo," %5.2f",
	       *(Info->con_all+m+i)*100.0/(float)Info->ncont);
	   }
	 fprintf(Finfo,"\n");
       }
     
     /* print also row deviations (variability of contacts) */
     if(op->info2)
       {
	 fprintf(Finfo,
		" How contact (site) statistics varies from prot to prot:\n");
	 fprintf(Finfo,
		" average number of contacts (sites) of a specific type: \n");
	 for(i=0;i<op->n_type;++i) fprintf(Finfo,"   %s",
					   dig2amino(*(Alph+i)));
	 fprintf(Finfo,"\n");

	 for(i=0;i<k;++i)
	   {
	     l=0;
	     for(j=0;j<op->n_type;++j) 
	       {
		 l+=*(Imod+j);
		 m=l-*(Imod+j);
		 av_cont=(float)(*(cont_all+m+i))/(float)op->max_prot;
		 if(*(Imod+j)<i+1) fprintf(Finfo," %s","     ");
		 else fprintf(Finfo," %5.2f",av_cont);
	       }
	     fprintf(Finfo,"\n");
	   }    
	 fprintf(Finfo,
		" st deviations for contacts (sites) of a specific type:\n");
	 for(i=0;i<k;++i)
	   {
	     l=0;
	     for(j=0;j<op->n_type;++j) 
	       {
		 l+=*(Imod+j);
		 m=l-*(Imod+j);
		 av_cont=(float)(*(cont_all+m+i))/(float)op->max_prot;
		 diff=av_cont*av_cont;
		 diff=(float)(*(cont_all_dev+m+i))/(float)op->max_prot - diff;
		 if(*(Imod+j)<i+1) fprintf(Finfo," %s","     ");
		 else fprintf(Finfo," %5.2f",sqrt(diff));
	       }
	     fprintf(Finfo,"\n");
	   }     
       }
     
   }
 
}

/* info_cont_dec prints statistics regarding decoys contacts */
void info_cont_dec(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, 
		   int *Imod, int *Alph, int *Adr)
{
 int i,j,k,l,m,n;
 double norm,dtmp;
 

 fprintf(Finfo," Contact type statistics (decoys) : \n");
 norm=0.0;
 for(i=0;i<op->n_par;++i)
   { 
     norm+=*(Info->con_All_decoys+i);
   }
 Info->ncont_decoys=norm;
 fprintf(Finfo,
	 " Decoy contacts (sites) as percent of tot #ncont_dec=%15.2f [mln]\n",
	 norm);
 for(i=0;i<op->n_type;++i) fprintf(Finfo,"   %s",dig2amino(*(Alph+i)));
 fprintf(Finfo,"\n");

 if(op->eij)
   {
     for(i=0;i<op->n_type;++i)
       {
	 for(j=0;j<op->n_type;++j) 
	   {
	     if(i>j) k=mat2vec(&j,&i,&op->n_type);
	     else k=mat2vec(&i,&j,&op->n_type);
	     fprintf(Finfo," %5.2f",
		     *(Info->con_All_decoys+k)*100.0/norm);
	   }
	 fprintf(Finfo,"\n");
       }    
   }

 if(op->ein || op->einm) 
   {
     k=*(Imod); /* take the number of parameters for first ama */
     for(i=1;i<op->n_type;++i)
       {
	 if(*(Imod+i)>k) k=*(Imod+i); /* take the maximum num of para */
       }
     for(i=0;i<k;++i)
       {
	 l=0;
	 for(j=0;j<op->n_type;++j) 
	   {
	     l+=*(Imod+j);
	     m=l-*(Imod+j);
	     if(*(Imod+j)<i+1) fprintf(Finfo," %s","     ");
	     else fprintf(Finfo," %5.2f",
			  *(Info->con_All_decoys+m+i)*100.0/norm);
	   }
	 fprintf(Finfo,"\n");
       }
     
   } 
}


/* info_pot writes down the parameters of the employed potential */
void info_pot(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, int *Imod,
		  int *Alph, int *Adr)
{
 int i,j,k,l,m,n;

 fprintf(Finfo," Potential that is used: \n");
 for(i=0;i<op->n_type;++i) fprintf(Finfo,"   %s",dig2amino(*(Alph+i)));
 fprintf(Finfo,"\n");
 
 if(op->eij)
   {
     for(i=0;i<op->n_type;++i)
       {
	 for(j=0;j<op->n_type;++j) 
	   {
	     if(i>j) k=mat2vec(&j,&i,&op->n_type);
	     else k=mat2vec(&i,&j,&op->n_type);
	     fprintf(Finfo," %5.2f",*(Potn+k));
	   }
	 fprintf(Finfo,"\n");
       }     
   }

 if(op->ein || op->einm) 
   {
     k=*(Imod); /* take the number of parameters for ALA */
     for(i=1;i<op->n_type;++i)
       {
	 if(*(Imod+i)>k) k=*(Imod+i); /* take the maximum num of para */
       }
     for(i=0;i<k;++i)
       {
	 l=0;
	 for(j=0;j<op->n_type;++j) 
	   {
	     l+=*(Imod+j);
	     m=l-*(Imod+j);
	     if(*(Imod+j)<i+1) fprintf(Finfo," %s","     ");
	     else fprintf(Finfo," %5.2f",*(Potn+m+i));
	   }
	 fprintf(Finfo,"\n");
       }     
   }
 
}

/* info_ene_tot prints the total energy contributions to all the native aligns */
void info_ene_tot(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, int *Imod,
		  int *Alph, int *Adr)
{
 int i,j,k,l,m,n;
 unsigned long int *cont_all;
 
 cont_all=Info->con_all;

 fprintf(Finfo," Total (whole set) stabilization contributions: \n");
 for(i=0;i<op->n_type;++i) fprintf(Finfo,"   %s",dig2amino(*(Alph+i)));
 fprintf(Finfo,"\n");

 if(op->eij)
   {
     for(i=0;i<op->n_type;++i)
       {
	 for(j=0;j<op->n_type;++j) 
	   {
	     if(i>j) k=mat2vec(&j,&i,&op->n_type);
	     else k=mat2vec(&i,&j,&op->n_type);
	     fprintf(Finfo," %5.0f",*(cont_all+k)*(*(Potn+k)));
	   }
	 fprintf(Finfo,"\n");
       }    
   }

 if(op->ein || op->einm) 
   {
     k=*(Imod); /* take the number of parameters for ALA */
     for(i=1;i<op->n_type;++i)
       {
	 if(*(Imod+i)>k) k=*(Imod+i); /* take the maximum num of para */
       }
     for(i=0;i<k;++i)
       {
	 l=0;
	 for(j=0;j<op->n_type;++j) 
	   {
	     l+=*(Imod+j);
	     m=l-*(Imod+j);
	     if(*(Imod+j)<i+1) fprintf(Finfo," %s","     ");
	     else fprintf(Finfo," %5.0f",*(cont_all+m+i)*(*(Potn+m+i)));
	   }
	 fprintf(Finfo,"\n");
       } 
     fprintf(Finfo,"\n");
 
   }
 
}

/* info_effpair_einm prints 5*5 effective pairwise matrices for THOM2 */
void info_effpair_einm(FILE *Finfo, k_info *Info, float *Potn, k_opt *op, 
		       int *Imod, int *Alph, int *Adr)
{
 int i,j,k,l,m,n;
 float diff,av_diff,tot_diff,av_cont;
 
 /* Assumption: first contact layer of 5 classes and second of 3 classes */
 tot_diff=0;
 fprintf(Finfo," Effective pairwise interactions of the employed potential:\n"); 
 for(i=0;i<op->n_type;++i)
   {
     fprintf(Finfo,"   %s block  \n",dig2amino(*(Alph+i)));
     for(j=0;j<op->n_type;++j) 
       {
	 fprintf(Finfo," %s - %s \n",
		 dig2amino(*(Alph+i)),dig2amino(*(Alph+j)));
	 av_diff=0;
	 for(k=0;k<5;++k) 
	   {
	     n=0;
	     if(k>0 && k<3) n=1;
	     if(k>2) n=2;
	     for(l=0;l<5;++l) 
	       { 
		 m=0;
		 if(l>0 && l<3) m=1;
		 if(l>2) m=2;
		 fprintf(Finfo," %5.2f",
			 *(Potn+15*i+3*k+m)+*(Potn+15*j+3*l+n));
		 diff=*(Potn+15*i+3*k+m)+*(Potn+15*j+3*l+n);
		 diff-=*(Potn+15*i+3*l+n)+*(Potn+15*j+3*k+m);
		 if(k>l) 
		   {
		     av_diff+=diff*diff;
		     tot_diff+=diff*diff;
		   }
	       }
	     fprintf(Finfo,"\n");
	   }
	 fprintf(Finfo," d= %5.2f \n",av_diff);
       }
     fprintf(Finfo,"\n");
   }  
 fprintf(Finfo," tot_d= %5.2f \n",tot_diff);
 
}

/* info computes moments and histograms of the energy diff. distribution */
void info_prot(FILE *Finfo,char *R_name,float *ene,float *Potn,
	       k_opt *op,int *Imod,int *contR,int *Alph,int *Adr)
{
 int i,j,k,l,m;

 fprintf(Finfo," Energy of %s = %8.2f and it is due to ...\n",R_name,*ene);
 for(i=0;i<op->n_type;++i) fprintf(Finfo,"   %s",dig2amino(*(Alph+i)));
 fprintf(Finfo,"\n");

 if(op->eij)
   {
     for(i=0;i<op->n_type;++i)
       {
	 for(j=0;j<op->n_type;++j) 
	   {
	     if(i>j) k=mat2vec(&j,&i,&op->n_type);
	     else k=mat2vec(&i,&j,&op->n_type);
	     fprintf(Finfo," %5.1f",*(contR+k)*(*(Potn+k)));
	   }
	 fprintf(Finfo,"\n");
       }    
   }

 if(op->ein || op->einm) 

   {
     k=*(Imod); /* take the number of parameters for ALA */

     /* take into account varying number of para per ama */
     for(i=1;i<op->n_type;++i)
       {
	 if(*(Imod+i)>k) k=*(Imod+i); /* take the maximum num of para */
       }
     for(i=0;i<k;++i)
       {
	 l=0;
	 for(j=0;j<op->n_type;++j) 
	   {
	     l+=*(Imod+j);
	     m=l-*(Imod+j);
	     if(*(Imod+j)<i+1) fprintf(Finfo,"      ");
	     else fprintf(Finfo," %5.1f",*(contR+m+i)*(*(Potn+m+i)));
	   }
	 fprintf(Finfo,"\n");
       } 
     fprintf(Finfo,"\n");
 
   }

 
}


void prn_setup(a_file *F,k_opt *op,int *Alph)
{
  int i;
  char l_seq[MAXS],stmp[MAXS];

  /* conversion of input/output handling not finished */
      
  /* print what was found in Model */
  printf(" Using alphabet:  \n");
  for(i=0;i<op->n_type;++i) printf(" %s",dig2amino(*(Alph+i)));
  printf("  \n");
  if(op->eij) printf(" Using e_ij model \n");
  if(op->ein) printf(" Using e_i(n) model \n");
  if(op->einm) printf(" Using e_i(n,m) model \n");
  if(op->evdw) 
    {
      printf(" Using e_vdw_ij model \n");
      printf(" using Lennard-Jones type %d-%d \n",op->rep_pow,op->atr_pow);
    }
  printf(" n_para = %d \n",op->n_par);
  printf(" contact cutoff R_cut = %5.2f \n",op->r_cut);

  /* print some information about the options used while threading */
  if(op->exam || op->learn || op->decoy)
    {
      if(op->exam) strcpy(l_seq,F->f_exam);
      else strcpy(l_seq,F->f_seq);
      if(op->decoy) printf(" Generating ineq for independent decoys ... \n");
      if(op->learn) printf(" Generating ineq using gapless threading ... \n");
      if(op->exam) printf(" Searching for the best alignment ... \n");
      if(op->end==op->max_prot) 
       printf(" Threading all the seq. from file %s (through all the stru.)\n",
	       l_seq);
      else
	printf(" Threading seq. from #%d to #%d (through all the stru.)\n",
	       op->beg,op->end);
      printf(" using potential file: %s \n",F->f_pot); 
      strcpy(stmp,"");
      if(op->learn || op->decoy) strcpy(stmp," (relative to native)");
      printf(" writing energies%s to file: %s \n",stmp,F->f_out);
      if(!op->noth) 
	{
	  if(op->up) printf(" only those having energy larger than %3.1f \n",
			    op->th_ene);
	  else printf(" only those having energy smaller than %3.1f \n",
		      op->th_ene);
	}
      if(op->wDc) 
	{
	  printf(" writing contact diff. onto file: %s \n", F->f_dcon);
	  /* if(op->freez) printf(" freezing %d par and storing rhs in: %s \n",
			       op->nfrz,F->f_rhs);*/
	}
      if(op->mps) printf(" writing contact diff. onto MPS file: %s \n",
			 F->f_mps);
      if(op->mps) 
	printf(" WARNING: primal is generated and it eats up a lot of RAM \n");
    }
  
} /* end of prn_setup */
 
void print_help(char *argv[], k_info *Info)
{
  int i;
  char *logo[]={"\n\t", 
  "LOOPP (Learning, Observing and Outputting Protein Patterns) ver. 2.000  \n\t",
  "Protein structure recognition by threading and sequence alignment,      \n\t",
  "developing new folding potentials, structure-structure alignment,       \n\t",
  "building non-redundant databases of folds ",
		                         "and other useful functions ... \n\n\t",
  "               by Jarek Meller and Ron Elber \n" };
  char *ex1[]={
    "To get LP inequalities by gapless threading of seqs from the \n",
    "database through all the structures in the database and write\n",
    "energy diff between decoys and natives structures into file scan.log\n" };

  char *ex2[]={
    "To get LP ineqs using pairwise poten my.pot (instead of default THOM2)\n",
    "and write contact diff between decoys and natives to file con_diff\n" };

  char *ex4[]={
    "To thread all the seqs from file SEQ through structures from XYZ and \n",
    "write differences in contacts between decoys and natives in MPS format\n"};

  char *ex5[]={
    "To get LP inequalities for explicitly generated (e.g. using MD) decoys \n", 
    "included in XYZ file (native structure included as first one in XYZ) \n" };

  char *ex6[]={
    "To thread sequences (in SWISS-PROT format) from file my.seq through all \n",
    "database structures, using default threading potential (THOM2) \n" };

  char *ex7[]={
    "To thread sequences from file SEQ (in the default LOOPP format) through \n",
    "all the structures in the database (by default included in XYZ file) \n" };

  char *ex9[]={
    "To perform sequence to sequence alignments for all the sequences \n",
    "from the file my.seq (in FASTA format) to all the database sequences \n"};

  char *ex10[]={
    "To perform structure to structure alignments for a structure from \n",
    "the PDB file 1MBO.pdb to all the database structures \n" };

  char *options[]={
    "\t -t <query_seq_file>       threading\n",
    "\t -s <query_seq_file>       sequence to sequence alignment\n",
    "\t -x <query_struct_file>    structure to structure alignment\n",
    "\t -q                        generating ineq for LP training\n",
    "\t -d                        LP ineq for explicit decoys\n",
    "\t -i <#level>               info and statistics, #level=1,2,3,4\n",
    "\t -l                        local alignments\n",
    "\t -fmt <format>             query seq (struct) file format \n",
    "\t                           that can be one of the following: \n",
    "\t               fasta (f)   FASTA \n",
    "\t               swiss (s)   SWISSPROT \n",
    "\t               pdb   (p)   PROTEIN DATA BANK \n",
    "\t               crd   (c)   CHARMM \n",
    "\t               1let  (1)   plain (simple one-letter code) \n",
    "\t               loopp (l)   LOOPP (default) \n",
    "\t -name <query_name>        query name for query_fmt=1let\n",
    "\t -pick <query> <match>     narrow the search to a given seq/stru\n",
    "\t -sd <#level>              search depth, #level=0,1,2,3\n",
    "\t -list <list_file>         multiple query files specified in list \n",
    "\t -einm                     THOM2 potentials (default) \n",
    "\t -eij                      pairwise potentials \n",
    "\t -ein                      THOM1 potentials \n",
    "\t -evdw                     LJ potentials \n",
    "\t -p <pot_file>             non-default potential\n",
    "\t -mod <model_file>         non-default energy functional\n",
    "\t -best <#best> <#shuff>    Z-score test range and convergence \n",
    "\t -nprn <#best>             print only #best matches \n",
    "\t -g <#pre_g_pen> <#g_pen>  gap penalties \n",
    "\t -cgap                     constant gap penalty (seq-seq) \n",
    "\t -cdel                     constant deletion penalty (threading)\n",
    "\t -nog                      use gapless threading for recognition\n",
    "\t -tl <#align_length_th>    threshold for alignment length \n",
    "\t -tz <#z_sc_th>            threshold for Z-score \n",
    "\t -tr <#rms_th>             threshold for RMS distance \n",
    "\t -tsid <#seqid_th>         threshold for sequence identity \n",
    "\t -te <#ene_th>             threshold for energy \n",
    "\t -u                        take energies larger then ene_th \n",
    "\t -seq <seq_file>           specify the sequence database file\n",
    "\t -xyz <coor_file>          specify the coordinates database \n",
    "\t -newCM                    enforce recreating Contact Maps \n",
    "\t -k <#col> <#n_col>        specify XYZ file format \n",
    "\t -o <out_file>             specify the name of the stdout file \n",
    "\t -log <best> <aligns>      specify names of the log files \n",
    "\t -mps <mps_file>           write LP constraints in MPS format \n",
    "\t -w <dcon_file>            write differences in contacts \n",
    "\t -frz <#n_frz>             freeze (make constant) #n_frz first param\n",
    "\t -m                        read potential in matrix format \n",
    "\t -rsc <#scale> <neg>       scale potential \n",
    "\t -pvdw <#rep> <#atr>       LJ powers specification \n",
    "\t -b <#start> <#end>        range of gapless threading in trainig\n",
    "\t -newDB <new_list>         build (new) database \n",
    "\t -nored <new_list>         remove redundancies from a database\n",
    "\t -trDB <#rms_th>           RMS distance for exclusion from a database\n",
    "\t -justQ                    stop after building query files \n",
    "\t -dbg                      debugging \n",
    "\t -h                        help \n",
    "\t -v                        print version signature \n" };

  for(i=0;i<7;++i) 
    {
      sprintf(Info->Sout,"%s",logo[i]);
      to_stnd_out(Info->Sout); 
    }
  
  sprintf(Info->Sout,"\n");
  to_stnd_out(Info->Sout);
  sprintf(Info->Sout,"Usage: %s [Option1] [Option2] ... \n\n",argv[0]);
  to_stnd_out(Info->Sout);

  /* print the most important options */
  sprintf(Info->Sout,"Options:  \n");
  to_stnd_out(Info->Sout);
  for(i=0;i<57;++i) 
    {
      sprintf(Info->Sout,"%s",options[i]);
      to_stnd_out(Info->Sout);
    }

  /* now print some examples */
  sprintf(Info->Sout,"\n\nExamples: \n\n");
  to_stnd_out(Info->Sout);
  sprintf(Info->Sout,"\t Designing folding potentials ...\n\n");
  to_stnd_out(Info->Sout);

  for(i=0;i<3;++i) 
    {
      sprintf(Info->Sout,"%s",ex1[i]);
      to_stnd_out(Info->Sout);
    }
  sprintf(Info->Sout,"type> %s -q -o scan.log  \n\n",argv[0]);
  to_stnd_out(Info->Sout);

  for(i=0;i<2;++i) 
    {
      sprintf(Info->Sout,"%s",ex2[i]);
      to_stnd_out(Info->Sout);
    }
  sprintf(Info->Sout,
	  "type> %s -q -eij -p my.pot -o scan.log -w con_diff \n\n",argv[0]);
  to_stnd_out(Info->Sout);

  sprintf(Info->Sout,
	  "To filter out inequalities with energy gaps larger than 5  \n");
  to_stnd_out(Info->Sout);
  sprintf(Info->Sout,
	  "type> %s -q -eij -p my.pot -o out -te 5 -w con_diff \n\n",argv[0]);
  to_stnd_out(Info->Sout);

  for(i=0;i<2;++i) 
    {
      sprintf(Info->Sout,"%s",ex4[i]);
      to_stnd_out(Info->Sout);
    }
  sprintf(Info->Sout,
	  "type> %s -q -eij -p my.pot -seq SEQ -xyz XYZ -mps /tmp/out.mps\n\n",
	  argv[0]);
  to_stnd_out(Info->Sout);

  for(i=0;i<2;++i) 
    {
      sprintf(Info->Sout,"%s",ex5[i]);
      to_stnd_out(Info->Sout);
    }
  sprintf(Info->Sout,"type> %s -eij -p my.pot -seq SEQ -xyz XYZ -d \n",argv[0]);
  to_stnd_out(Info->Sout);

  sprintf(Info->Sout,"\n\t Fold recognition ...\n\n");
  to_stnd_out(Info->Sout);

  for(i=0;i<2;++i) 
    {
      sprintf(Info->Sout,"%s",ex6[i]);
      to_stnd_out(Info->Sout);
    }
  sprintf(Info->Sout,"type> %s -t my.seq -fmt swiss \n\n",argv[0]);
  to_stnd_out(Info->Sout);

  for(i=0;i<2;++i) 
    {
      sprintf(Info->Sout,"%s",ex7[i]);
      to_stnd_out(Info->Sout);
    }
  sprintf(Info->Sout,"type> %s -t SEQ  \n\n",argv[0]);
  to_stnd_out(Info->Sout);

  sprintf(Info->Sout,
	  "To refine the Z-scores for a pair (1mba vs 1lh2) of proteins only\n");
  to_stnd_out(Info->Sout);
  sprintf(Info->Sout,"type> %s -t SEQ -pick 1mba 1lh2 -best 1 1000\n\n",argv[0]);
  to_stnd_out(Info->Sout);

  for(i=0;i<2;++i) 
    {
      sprintf(Info->Sout,"%s",ex9[i]);
      to_stnd_out(Info->Sout);
    }
  sprintf(Info->Sout,"type> %s -s my.seq -fmt fasta \n \n",argv[0]);
  to_stnd_out(Info->Sout);

  for(i=0;i<2;++i) 
    {
      sprintf(Info->Sout,"%s",ex10[i]);
      to_stnd_out(Info->Sout);
    }
  sprintf(Info->Sout,"type> %s  -x 1MBO.pdb -fmt pdb  \n\n\n",argv[0]);
  to_stnd_out(Info->Sout);

  sprintf(Info->Sout,"See loopp_doc.html for more details \n");
  to_stnd_out(Info->Sout);
  exit(0);

}

