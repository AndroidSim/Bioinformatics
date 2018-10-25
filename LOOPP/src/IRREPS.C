/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000       v. 2.000                                */
/*                                                                */
/* BUILD INTERNAL REPRESENTATION, CONTACT MAPS AND OTHER BLOCKS   */
 
#include "loopp.h"
   

/* build_irrep computes contact maps and reference energies */ 
void build_irrep(k_opt *op, a_file *F, solv_sh *ssh, k_prot *prot,
		 k_info *Info, k_Model *model)
{
  int i_prot,n,Lstru,Lseq,n_cont,i,j,k,match,to_read;
  int *contR,*IJcont,*Imod,*Alph,*Adr,*Seq;
  float *Potn,*R_ene,*X,*Y,*Z,*contRvdw,*r12ij,*r6ij,*rij;
  char S_name[MAXS],X_name[MAXS];
  
  /* Use S_name for seq names and X_name for those of stru names  */ 

  /* get local shortcuts */
  contR=prot->conR;
  IJcont=prot->IJcont;
  Imod=model->Imod;
  Potn=model->Potn;
  X=prot->X_res;
  Y=prot->Y_res;
  Z=prot->Z_res;
  R_ene=prot->R_ene;
  Seq=prot->seq;
  contRvdw=prot->conRvdw;
  r12ij=prot->IJr12;
  r6ij=prot->IJr6;
  rij=prot->native_IJrij;
  Alph=model->Alph;
  Adr=model->Adr;

 /* check first consistency of Contact Maps with DataBase files  */
 check_CM(op,F,Info,IJcont,Seq,r12ij,r6ij,rij,Alph,Adr);
  
 i_prot=0;
 n=0;
 while(read_seq(F->seq,&Lseq,S_name,Seq,&i_prot,&op->n_type,Alph,Adr))
   {
     if(op->r_cont) read_contMap(F,op,IJcont,&Lstru,X_name,r12ij,r6ij,rij);
     else /* prepare new CMs */
       {
	 read_coor(F->coor,op,&Lstru,X_name,X,Y,Z);
	 if(Lseq!=Lstru || strcmp(S_name,X_name)!=0)
	   {
	     printf("Mismatch between coor & seq at %s\n",S_name);
	     exit(1);
	   }
	 write_contMap(F,IJcont,&Lstru,X_name,X,Y,Z,op,r12ij,r6ij,rij);
       }

     /* prepare contact vectors and energies for each DB template */
     write_contType(F,op,Imod,contR,IJcont,&n,Seq,&Lseq,X_name,
		    ssh->Ncont,ssh,contRvdw,r12ij,r6ij);
     *(R_ene+i_prot-1)=energy(op,contR,Potn,contRvdw);

     /* copy names into Info.Names table */
     strcpy((Info->Names+i_prot-1)->mstr,X_name);
     if(op->debug) printf(" %s \n",(Info->Names+i_prot-1)->mstr);

     if(op->info)
       {
	 add_cont(&op->n_par,contR,&Info->ncont,Info->con_all,
		  Info->con_all_dev);
	 add_res(Info,Seq,&Lseq,Adr,Alph);
	 /* print additional info to a file */
	 if(op->info2) info_prot(F->iprot,X_name,R_ene+i_prot-1,
				 Potn,op,Imod,contR,Alph,Adr);
       }
   }

 /* close and open files to be read only from now on */
 to_read=1;
 switch_rw_CM(op,F,Info,&to_read);
 fclose(F->type);
 F->type=fopen("cont_by_Type","r");
 rewind(F->seq);
 rewind(F->coor);

 /* prepare query files in the loopp format */
 if(op->exam  || op->align || op->strucmp) build_query(F,op,Info,prot);

}

/* check_ContactMap decides whether the existing CM should be rebuilt */
int check_CM(k_opt *op, a_file *F, k_info *Info, int *IJcont, int *Seq,
	     float *r12ij, float *r6ij, float *rij, int *Alph, int *Adr)
{
 int i_prot,k,Lstru,Lseq,match,to_read;
 char S_name[MAXS],X_name[MAXS]; 

 /* compare SEQ and ContMap to decide */ 
 if(op->r_cont)
   {
     i_prot=0;
     match=1;
     if(op->debug) 
       {
	 sprintf(Info->Sout," Checking Contact Map files ...\n");
	 to_stnd_out(Info->Sout);
       }
     while(match && 
	   read_seq(F->seq,&Lseq,S_name,Seq,&i_prot,&op->n_type,Alph,Adr))
       {
	 if(!read_contMap(F,op,IJcont,&Lstru,X_name,r12ij,r6ij,rij)) match=0;
	 if(Lseq!=Lstru || strcmp(S_name,X_name)!=0) match=0;
	 if(op->debug) 
	   {
	     sprintf(Info->Sout,"\n %d stru %s seq %s",i_prot,X_name,S_name);
	     to_stnd_out(Info->Sout);
	   }
       }
     if(!match)
       {
	 op->r_cont=0;
	 to_read=0; /* switch to write */
	 switch_rw_CM(op,F,Info,&to_read);
       }
     else
       {
	 to_read=1; /* rewind to read again */
	 switch_rw_CM(op,F,Info,&to_read);
       }
     rewind(F->seq);
   }
 
 /* add small info if any of the CMs is to be rebuilt */
 if(!op->r_cont)
   {
     sprintf(Info->Sout,"ContMap will be rebuilt \n");
     to_stnd_out(Info->Sout);
     if(op->distance)
       {
	 sprintf(Info->Sout,"ContMap_rij will be rebuilt \n");
	 to_stnd_out(Info->Sout);     
       }
     if(op->evdw)
       {
	 sprintf(Info->Sout,"ContMap_vdw will be rebuilt \n");
	 to_stnd_out(Info->Sout);
       }
     
     k=(op->col-1)*3;
     sprintf(Info->Sout," Reading columns %d to %d of %s to get coordinates\n",
	    k+1,k+3,F->f_coor);
     to_stnd_out(Info->Sout);
   }
 
 return(1);
 
}

/* switch_rw_CM changes the mode of access to ContMap file */
int switch_rw_CM(k_opt *op, a_file *F, k_info *Info, int *to_read)
{
    
 /* close and open (or rewind) files to be written (or read) only */

 if(*to_read)
   {
     if(op->r_cont) /* continue reading after rewinding */
       {    
	 rewind(F->cont);
	 if(op->distance) rewind(F->rij);
	 if(op->evdw)
	   {
	     rewind(F->vdw_a);
	     rewind(F->vdw_b);
	   }	 
       }
     else /* switch from writing to reading */
       {
	 fclose(F->cont);
	 F->cont=fopen(F->f_cont,"r");
	 if(op->distance)
	   {
	     fclose(F->rij);
	     F->rij=fopen(F->f_rij,"r");     
	   }
	 if(op->evdw)
	   {
	     fclose(F->vdw_a);
	     F->vdw_a=fopen(F->f_vdw_a,"r");
	     fclose(F->vdw_b);
	     F->vdw_b=fopen(F->f_vdw_b,"r");
	   }
       }
   }
 else /* if to_write switch from reading to writing */
   {
     
     fclose(F->cont);
     F->cont=fopen(F->f_cont,"w");
     if(op->distance)
       {
	 fclose(F->rij);
	 F->rij=fopen(F->f_rij,"w");     
       }
     if(op->evdw)
       {
	 fclose(F->vdw_a);
	 F->vdw_a=fopen(F->f_vdw_a,"w");
	 fclose(F->vdw_b);
	 F->vdw_b=fopen(F->f_vdw_b,"w");
       }
   }    
 /* add if NULL then return 0 */
 return(1);
 
}

/* get_contMap finds for each residue i of struc R_name its 
   contacts with other residues j_1,j_2 ...  (using coor X,Y,Z) 
   and puts them to vector IJ */
void get_contMap(int *IJ,int *Lstru,char *R_name,float *X,float *Y,float *Z,
		 k_opt *op)
{
  /* not upgraded to use it for vdw - check write_contMap !! */
 int i,j; 
 float Rij2,r0,r1,r2, *jx, *jy, *jz;

 r0=R_low*R_low;
 r1=op->r_cut*op->r_cut;
 if(op->debug) printf("\n struc %s of %d residues",R_name,*Lstru); 
 for(i=0;i<*Lstru-DELTA;++i)
   {
   *IJ=i;  
   jx=X+DELTA; jy=Y+DELTA; jz=Z+DELTA;  
   for(j=i+DELTA;j<*Lstru;++j)
     {
      Rij2=(*X-*jx)*(*X-*jx)+(*Y-*jy)*(*Y-*jy)+(*Z-*jz)*(*Z-*jz);
      if(Rij2>r0 && Rij2<=r1) 
	{
	  ++IJ;
	  *IJ=j;
	}
      ++jx;  ++jy;  ++jz;
     }
   ++X; ++Y; ++Z;
   ++IJ;
   }
 *IJ=i; /* to end somehow the vector - this last value is used later */
}

/* read_contMap simply reads existing contact map from Fcont */
int read_contMap(a_file *F,k_opt *op,int *IJ,int *Lstru,char *R_name,
		 float *r12ij, float *r6ij, float *rij)
{

 int i,j,eof,i1,*j1,k,kflag; 
 float f1,f2;

 j1=malloc(op->max_cont*sizeof(int));
 if(j1==NULL) return(0);

 eof=fscanf(F->cont,"%s%d",R_name,Lstru);
 if(eof==EOF) return(0);

 if(op->distance) 
   {  
     eof=fscanf(F->rij,"%s%d",R_name,Lstru);
     if(eof==EOF) return(0);
   }

 if(op->evdw) 
   {  
     eof=fscanf(F->vdw_a,"%s%d",R_name,Lstru);
     if(eof==EOF) return(0);
     eof=fscanf(F->vdw_b,"%s%d",R_name,Lstru);
     if(eof==EOF) return(0);
   }
 
 if(op->debug) printf("\n starting with %s\t%d\n",R_name,*Lstru); 
 kflag=1;
 for(i=0;i<*Lstru-DELTA;++i)
   { 
    if(kflag)
      {
	fscanf(F->cont,"%d",&i1); 
	*IJ=i1; ++IJ;
	 if(op->distance) 
	   { 
	     fscanf(F->rij,"%f",&f1);
	     *rij=f1; ++rij;
	   }
	if(op->evdw)
	  {
	    fscanf(F->vdw_a,"%f",&f1); 
	    fscanf(F->vdw_b,"%f",&f2); 
	    *r12ij=f1; ++r12ij;
	    *r6ij=f2; ++r6ij;
	  }
      }
    k=0; kflag=1;
    while(kflag && k<op->max_cont)
      {
	fscanf(F->cont,"%d",&j1[k]);
	if(i1+1==j1[k]) 
	  {
	    kflag=0;
	    *IJ=j1[k]; ++IJ; 
	    i1=j1[k]; 
	  }
	else
	  {
	    *IJ=j1[k]; ++IJ; 
	    ++k;
	  }
	if(op->distance) 
	   { 
	     fscanf(F->rij,"%f",&f1); 
	     *rij=f1; ++rij;
	   }
	if(op->evdw)
	  {
	    fscanf(F->vdw_a,"%f",&f1); 
	    fscanf(F->vdw_b,"%f",&f2); 
	    *r12ij=f1; ++r12ij;
	    *r6ij=f2; ++r6ij;
	  }
      }

   }
 free(j1);
 return(1);

}

/* read_contMap simply reads existing contact map from Fcont */
int read_Query_CM(a_file *F, k_opt *op, int *IJ, int *Lstru, char *R_name,
		 float *rij)
{

 int i,j,eof,i1,*j1,k,kflag; 
 float f1,f2;

 j1=malloc(op->max_cont*sizeof(int));
 if(j1==NULL) return(0);

 eof=fscanf(F->qcont,"%s%d",R_name,Lstru);
 if(eof==EOF) return(0);

 if(op->distance) 
   {  
     eof=fscanf(F->qrij,"%s%d",R_name,Lstru);
     if(eof==EOF) return(0);
   }

 if(op->debug) printf("\n starting with %s\t%d\n",R_name,*Lstru); 
 kflag=1;
 for(i=0;i<*Lstru-DELTA;++i)
   { 
    if(kflag)
      {
	fscanf(F->qcont,"%d",&i1); 
	*IJ=i1; ++IJ;
	 if(op->distance) 
	   { 
	     fscanf(F->qrij,"%f",&f1);
	     *rij=f1; ++rij;
	   }
      }
    k=0; kflag=1;
    while(kflag && k<op->max_cont)
      {
	fscanf(F->qcont,"%d",&j1[k]);
	if(i1+1==j1[k]) 
	  {
	    kflag=0;
	    *IJ=j1[k]; ++IJ; 
	    i1=j1[k]; 
	  }
	else
	  {
	    *IJ=j1[k]; ++IJ; 
	    ++k;
	  }
	if(op->distance) 
	   { 
	     fscanf(F->qrij,"%f",&f1); 
	     *rij=f1; ++rij;
	   }
      }

   }
 free(j1);
 return(1);

}

/* write_contMap finds for each residue i of struc R_name its 
   contacts with other residues j_1,j_2 ...  (using coor X,Y,Z) 
   and puts them to both: file Fcont and vector IJ */
void write_contMap(a_file *F, int *IJ, int *Lstru, char *R_name,
		   float *X, float *Y, float *Z, k_opt *op,
		   float *r12ij, float *r6ij, float *rij)
{

 int i,j,k,m; 
 float Rij2,r0,r1,r2,r6,r12,scale,*jx,*jy,*jz,u;

 r0=R_low*R_low;
 r1=op->r_cut*op->r_cut;
 fprintf(F->cont,"%s\t%d\n",R_name,*Lstru); 
 if(op->distance) fprintf(F->rij,"%s\t%d\n",R_name,*Lstru);
 if(op->evdw) 
   {
     scale=(float)RESC_VDW*RESC_VDW;
     fprintf(F->vdw_a,"%s\t%d\n",R_name,*Lstru); 
     fprintf(F->vdw_b,"%s\t%d\n",R_name,*Lstru);      
   }

 for(i=0;i<*Lstru-DELTA;++i)
   {
   *IJ=i;  fprintf(F->cont,"%d\t\t",i); 
   if(op->distance) 
     {
       *rij=(float)i;
       fprintf(F->rij,"%d\t\t",i);
     }
   if(op->evdw)
     { 
       *r12ij=(float)i;
       *r6ij=(float)i;
       fprintf(F->vdw_a,"%d\t\t",i);
       fprintf(F->vdw_b,"%d\t\t",i);
     }
   
   k=0;
   jx=X+DELTA; jy=Y+DELTA; jz=Z+DELTA;  
   for(j=i+DELTA;j<*Lstru;++j)
     {
      Rij2=(*X-*jx)*(*X-*jx)+(*Y-*jy)*(*Y-*jy)+(*Z-*jz)*(*Z-*jz);
      if(Rij2>r0 && Rij2<=r1) 
	{
	  ++IJ;
	  *IJ=j; fprintf(F->cont,"%d\t",j); 
	  if(op->debug) printf("\n Contact between %d and %d",i,j); 
	  ++k; 
	  if(k>=op->max_cont) 
	    {
	      printf(" WARNING: max_cont too small and some contacts");
	      printf(" will not be accounted for \n");
	    }
	  if(op->distance) 
	    {
	      ++rij;
	      fprintf(F->rij,"%e\t",sqrt(Rij2));
	    }
	  if(op->evdw)
	    {
	      ++r12ij;
	      ++r6ij;
	      r6=scale/Rij2;
	      r12=r6;
	      u=r6; 
	      /* get the requested powers for L-J */
	      for(m=0;m<op->atr_pow/2-1;++m) r6*=u;
	      for(m=0;m<op->rep_pow/2-1;++m) r12*=u;
	      *r12ij=r12;
	      *r6ij=r6; 
	      fprintf(F->vdw_a,"%e\t",r12);
	      fprintf(F->vdw_b,"%e\t",r6); 
	    }
	}
      ++jx;  ++jy;  ++jz;
     }
   ++X; ++Y; ++Z;
   ++IJ; fprintf(F->cont,"\n");
   if(op->distance) 
     {
       ++rij;
       fprintf(F->rij,"\n");
     }
   if(op->evdw)
     { 
       ++r12ij;
       ++r6ij;
       fprintf(F->vdw_a,"\n");
       fprintf(F->vdw_b,"\n");
     }
   }
 fprintf(F->cont,"%d\n",i); 
 *IJ=i; /* to end somehow the vector - this last value is used later */ 
 if(op->distance)
   {
     *rij=(float)i;
     fprintf(F->rij,"%d\n",i);
   }
 if(op->evdw)
   { 
     /* to have the same num of entries - some bytes are wasted here */
     *r12ij=(float)i;
     *r6ij=(float)i;
     fprintf(F->vdw_a,"%d\n",i);
     fprintf(F->vdw_b,"%d\n",i);
   }

}

/* write_Query_ContMap finds for each residue i of struc R_name its 
   contacts with other residues j_1,j_2 ...  (using coor X,Y,Z) 
   and puts them to files f_qcont and f_qrij only */
void write_Query_CM(a_file *F, int *Lstru, char *R_name,
		    float *X, float *Y, float *Z, k_opt *op)
{

 int i,j,k,m; 
 float Rij2,r0,r1,r2,r6,r12,scale,*jx,*jy,*jz,u;

 r0=R_low*R_low;
 r1=op->r_cut*op->r_cut;
 fprintf(F->qcont,"%s\t%d\n",R_name,*Lstru); 
 if(op->distance) fprintf(F->qrij,"%s\t%d\n",R_name,*Lstru);

 for(i=0;i<*Lstru-DELTA;++i)
   {
   fprintf(F->qcont,"%d\t\t",i); 
   if(op->distance) fprintf(F->qrij,"%d\t\t",i);

   k=0;
   jx=X+DELTA; jy=Y+DELTA; jz=Z+DELTA;  
   for(j=i+DELTA;j<*Lstru;++j)
     {
      Rij2=(*X-*jx)*(*X-*jx)+(*Y-*jy)*(*Y-*jy)+(*Z-*jz)*(*Z-*jz);
      if(Rij2>r0 && Rij2<=r1) 
	{
	  fprintf(F->qcont,"%d\t",j); 
	  if(op->debug) printf("\n Contact between %d and %d",i,j); 
	  ++k; 
	  if(k>=op->max_cont) 
	    {
	      printf(" WARNING: max_cont too small and some contacts");
	      printf(" will not be accounted for \n");
	    }
	  if(op->distance) fprintf(F->qrij,"%e\t",sqrt(Rij2));
	}
      ++jx;  ++jy;  ++jz;
     }
   ++X; ++Y; ++Z;
   fprintf(F->qcont,"\n");
   if(op->distance) fprintf(F->qrij,"\n");

   }

 fprintf(F->qcont,"%d\n",i); 
 /* to end somehow the vector - this last value is used later */ 
 if(op->distance) fprintf(F->qrij,"%d\n",i);

}


/* get_contType calculates number of contacts according to their types */
void get_contType(k_opt *op,int *Imod,int *contT,int *IJ,int *start,int *seq,
		  int *Lseq,int *start2,int *Icont,solv_sh *ssh,
		  float *contTvdw,float *r12ij,float *r6ij)
{
 int i,j,k,l,m,alpha,beta,nbij,*IJ_loc;


 /* Icont=ssh->Ncont; Icont is local copy for first solv sh */
 if(op->ein || op->einm)
   for(i=0;i<*Lseq;++i)
     *(Icont+i)=0;
 
 if(op->evdw) 
   {
     nbij=op->n_type*(op->n_type+1)/2;
     for(i=0;i<op->n_par;++i)
       *(contTvdw+i)=0;
   }
 else
   {
     for(i=0;i<op->n_par;++i)
       *(contT+i)=0;    
   }
 
 /* first skip everything up to the point where the threading starts */
 /* fixed for a while !!! when JUMP changes it won't work */
 if(op->gapless)
   {
     l=*start;
     IJ+=*start2;
     if(op->evdw) 
       {
	 r12ij+=*start2;
	 r6ij+=*start2;
       }
     if(*IJ==l) ++*start2;
     i=0;
     while(*(IJ+i)!=l+1) { ++*start2; ++i; }
   }
 else
   {
     if(op->evdw)
       {
	 for(i=0;i<*start;++i)
	   {
	     if(*IJ==i) 
	       {
		 ++IJ; 
		 ++r12ij;
		 ++r6ij;
	       }
	     while(*IJ!=i+1)  
	       {
		 ++IJ; 
		 ++r12ij;
		 ++r6ij;
	       }
	   }
       }
     else
       {
	 for(i=0;i<*start;++i)
	   {
	     if(*IJ==i) ++IJ;     
	     while(*IJ!=i+1) ++IJ; 
	   } 
       }

   }

 if(op->einm) IJ_loc=IJ; /* create a copy for einm */

 k=*start; 
 /* if threading is regular (i.e. with fixed JUMP) one may calculate where to 
    start in IJ during next call to avoid the previous loop  */
 for(i=0;i<*Lseq-DELTA;++i)
   {
     alpha=*(seq+i);
     if(*IJ==k)
       {
	 ++IJ; 
	 if(op->evdw) 
	   { 
	     ++r12ij;
	     ++r6ij;
	   }
       } 
     while(*IJ!=k+1) /* k+1 cannot be a neighbour by covalent exclusion */
       {
	 m=*IJ-*start;
	 if(m<*Lseq)
	   {
	     beta=*(seq+m);
	     if(op->ein || op->einm)
	       {
		 ++*(Icont+i);
		 ++*(Icont+m);
	       }
	     else
	       {
		 if(op->evdw)
		   {
		     if(alpha>beta)
		       {
			 *(contTvdw+mat2vec(&beta,&alpha,&op->n_type))
			   +=*r12ij;
			 /* the Bij parameters come after Aij */
			 *(contTvdw+mat2vec(&beta,&alpha,&op->n_type)+nbij)
			   +=*r6ij;
		       }
		     else
		       {
			 *(contTvdw+mat2vec(&alpha,&beta,&op->n_type))
			   +=*r12ij;
			 *(contTvdw+mat2vec(&alpha,&beta,&op->n_type)+nbij)
			   +=*r6ij;
		       }
		   }
		 else
		   {
		     if(alpha>beta)
		       ++*(contT+mat2vec(&beta,&alpha,&op->n_type));
		     else
		       ++*(contT+mat2vec(&alpha,&beta,&op->n_type));
		   }
	       }
	   }
	 ++IJ; 
	 if(op->evdw) 
	   { 
	     ++r12ij;
	     ++r6ij;
	   }
       }
     ++k;
   }

 if(op->einm) get_einm(op,Imod,contT,IJ_loc,start,seq,Lseq,ssh->Ncont,
		       ssh->Mcont0,ssh->Mcont1,ssh->Mcont2);
 
 if(op->ein)
   {
     /* CHANGE */
     if(op->gap_adr>=0)
       {	 
	 for(i=0;i<*Lseq;++i)
	   if((*(Icont+i)>0) || (op->gap_adr==*(seq+i))) 
	     ++*(contT+imod2vec(seq+i,Icont+i,Imod));
       }
     else
     /* CHANGE */  
       {	 
	 for(i=0;i<*Lseq;++i)
	   if(*(Icont+i)>0) ++*(contT+imod2vec(seq+i,Icont+i,Imod));
	 /* fix to avoid e_i(0) - remove if in line above to have it
	    ++*(contT+imod2vec(&seq[i],Icont+i,Imod)); */
       }
   }

}

/* get_score_eij calculates energy od a site given eij and frozen env */
float get_score_eij(k_opt *op,int *Imod,int *contT,int *IJ,int *start,
		    int *nat_Seq,int *Lseq,int *qseq_i_alpha,int *j_beta,
		    float *contTvdw,float *r12ij,float *r6ij,float *Potn,
		    float *gpen)
{
 int i,j,k,l,m,alpha,beta,nbij,l_sum,include,n_cont;


 /* general comment:
    this function is terribly slow and can be easily extended to calculate
    the identities of the neighbors for each site only once (and stored)
    instead of doing this for each DP table entry, which is as stupid as it
    can be (sorry for this quick patch) */

 if(op->evdw) 
   {
     nbij=op->n_type*(op->n_type+1)/2;
     for(i=0;i<op->n_par;++i)
       *(contTvdw+i)=0;
   }
 else
   {
     for(i=0;i<op->n_par;++i)
       *(contT+i)=0;    
   }
 n_cont=0;
 
 /* do not skip anything as some neighbours may be defined before this site */

 alpha=*qseq_i_alpha; /* query ama to be put on a site j_beta */
 k=*start; 
 /* start assumed to be zero to start from the beginning  */
 /* Lseq refers now to the length of the native (tmpl) seq (and thus stru) */
 if((*j_beta+1)>(*Lseq-DELTA)) l_sum=*Lseq-DELTA;
 else l_sum=*j_beta+1;
 
 for(i=0;i<l_sum;++i)
   {
     /* so now k is actually supposed to be equal to i at each step */
     if(*IJ==k)
       {
	 ++IJ; 
	 if(op->evdw) 
	   { 
	     ++r12ij;
	     ++r6ij;
	   }
       } 
     while(*IJ!=k+1) /* k+1 cannot be a neighbour by exclusion */
       {
	 include=0;
	 m=*IJ-*start;
	 if(k<*j_beta) /* among neighs m of site k might be j_beta */ 
	   {
	     if(m==*j_beta)
	       {
		 beta=*(nat_Seq+k);
		 include=1;
		 ++n_cont;
	       }
	   }
	 else /* if it is the site number j_beta take all neighs */
	   {
	     beta=*(nat_Seq+m);
	     include=1;
	     ++n_cont;
	   }
	 
	 if(include)
	   {
	     
	     if(op->evdw)
	       {
		 if(alpha>beta)
		   {
		     *(contTvdw+mat2vec(&beta,&alpha,&op->n_type))
		       +=*r12ij;
		     /* the Bij parameters come after Aij */
		     *(contTvdw+mat2vec(&beta,&alpha,&op->n_type)+nbij)
		       +=*r6ij;
		   }
		 else
		   {
		     *(contTvdw+mat2vec(&alpha,&beta,&op->n_type))
		       +=*r12ij;
		     *(contTvdw+mat2vec(&alpha,&beta,&op->n_type)+nbij)
		       +=*r6ij;
		   }
	       }
	     else
	       {
		 if(alpha>beta)
		   ++*(contT+mat2vec(&beta,&alpha,&op->n_type));
		 else
		   ++*(contT+mat2vec(&alpha,&beta,&op->n_type));
	       }
	   }
	 
	 ++IJ; 
	 if(op->evdw) 
	   { 
	     ++r12ij;
	     ++r6ij;
	   }
       }
     ++k;
   }
 
 /* for pairwise potentials gap penalty is set to be proportional to 
    n_cont+1 without referring to explicit GAP residue in the Potn */
 *gpen=(n_cont+1)*op->gap_pen;
 /* printf("\n site %d   n_cont %d \n",*j_beta,n_cont); */
 
 return(0.5*energy(op,contT,Potn,contTvdw));

}


/* get_env computes environmnet in terms of solvation shells */
void get_env(k_opt *op,int *Imod,int *contT,int *IJ,int *start,int *seq,
		  int *Lstru,int *start2,int *Icont,solv_sh *ssh,
		  float *r12ij,float *r6ij,int *Stru_row)
{
 int i,j,k,l,m,alpha,beta,nbij,*IJ_loc,i_rewind;


 /* Icont=ssh->Ncont; Icont is local copy for first solv sh */
 if(op->ein || op->einm)
   for(i=0;i<*Lstru;++i)
     {
       *(Icont+i)=0;
       *(Stru_row+i)=i;
     }
 
 if(op->einm) IJ_loc=IJ; /* create a copy for einm */

 k=*start; 
 /* get the Icont vector for the whole structure */
 for(i=0;i<*Lstru-DELTA;++i)
   {
     /* alpha=*(seq+i); */
     if(*IJ==k)
       {
	 ++IJ; 
	 if(op->evdw) 
	   { 
	     ++r12ij;
	     ++r6ij;
	   }
       } 
     while(*IJ!=k+1)
       {
	 m=*IJ-*start;
	 if(m<*Lstru) /* so now the whole structure is taken */
	   {
	     /* beta=*(seq+m); */
	     if(op->ein || op->einm)
	       {
		 ++*(Icont+i);
		 ++*(Icont+m);
	       }

	   }
	 ++IJ; 
	 if(op->evdw) 
	   { 
	     ++r12ij;
	     ++r6ij;
	   }
       }
     ++k;
   }

 /* IJ=IJ_loc; */
 
 if(op->einm) get_env_einm(op,Imod,contT,IJ_loc,start,seq,Lstru,ssh->Ncont,
		       ssh->Mcont0,ssh->Mcont1,ssh->Mcont2);

 /* now you should have each site of a structure characterized by N,Mcont */

}


/* get_env computes environmnet in terms of the first contact shell */
void get_env_ein(k_opt *op,int *Imod,int *IJ,int *Lstru,int *Icont)
{
 int i,j,k,l,m,start;

 /* assumption: all the arrays allocated irrespective of choice of op->ein */
 /* also: always start from the beginning i.e. start=0 */

 start=0;
 for(i=0;i<*Lstru;++i)
   {
     *(Icont+i)=0;
   }
 
 k=start; 
 /* get the Icont vector for the whole structure */
 for(i=0;i<*Lstru-DELTA;++i)
   {
     /* alpha=*(seq+i); */
     if(*IJ==k)
       {
	 ++IJ; 
       } 
     while(*IJ!=k+1)
       {
	 m=*IJ-start;
	 if(m<*Lstru) /* so now the whole structure is taken */
	   {
	     /* beta=*(seq+m); */
	     ++*(Icont+i);
	     ++*(Icont+m);
	   }
	 ++IJ; 
       }
     ++k;
   }


 /* now you should have each site of a structure characterized by Ncont */

}



/* get_score computes local energy contribution of a site */
float get_score(k_opt *op,int *Imod,int *n,int *seq,
		int *Lseq,int *m,int *Icont,float *Potn,
		int *Mcont0,int *Mcont1,int *Mcont2,float *gpen)
{
 int i,j,jg,k,l,alpha,beta,nbij,n_i;
 float score;
 
 score=0.0;
 *gpen=0.0;
 
 if(op->einm) 
   {
     /* assumption: n_par=3(sec. shell)*n_type(first shell) !!! fixed */
     j=0; n_i=*(Icont+*m);
     /* for n_type=4 use the next line */
     /* if(n_i>2 && n_i<6) j+=1; if(n_i>5 && n_i<9) j+=2; if(n_i>8) j+=3; */
     /* currently n_type=5 and thus n_par(per ama)=15 is assumed - fixed !!! */
     if(n_i>2) 
       {
	 ++j;
	 if(n_i>4) 
	   {
	     ++j;
	     if(n_i>6) 
	       {
		 ++j;
		 if(n_i>8) ++j;
	       }	     
	   }
       }
     j*=3;
     jg=j; /* to restart pointer for gap penalties */
     /* this loop could be avoided - Imod */
     /* if stru-stru just take the value for the first letter */
     if(!op->strucmp) for(k=0;k<*(seq+*n);++k) j+=*(Imod+k); 
     if(op->debug)
       {
	 printf("#%d n_i %d I0i %d \n",*m,n_i,*(Mcont0+*m));
	 printf(" addresses I0i %d I1i %d I2i %d \n",j,j+1,j+2); 
       }
     score+=*(Mcont0+*m)*(*(Potn+j));
     score+=*(Mcont1+*m)*(*(Potn+j+1));
     score+=*(Mcont2+*m)*(*(Potn+j+2));

     /* take energy terms for GAP res at site m */
     if(op->gap_adr>=0)
       {
	 /* rewind the pointer first - jg points to a proper gap type */
	 for(k=0;k<op->gap_adr;++k) jg+=*(Imod+k); 
	 /* if there is a lonely GAP take the additional (last) term */
	 if(n_i==0) *gpen=*(Potn+jg+*(Imod+op->gap_adr)-1);
	 else
	   {
	     *gpen+=*(Mcont0+*m)*(*(Potn+jg));
	     *gpen+=*(Mcont1+*m)*(*(Potn+jg+1));
	     *gpen+=*(Mcont2+*m)*(*(Potn+jg+2));
	   }
	 
	 /* rescue patch for the TOM2 model without penalty for lonely GAPS */
	 /* if(n_i==0) *gpen+=op->gap_pen; obsolate now */
	 if(op->debug) printf(" gap penalty = %f \n",*gpen);
       }
     
   }
 
 if(op->ein)
   {
     /* fix to avoid e_i(0) - remove if in line above to have it
	++*(contT+imod2vec(&seq[i],Icont+i,Imod)); */
     if(*(Icont+*m)>0)
       {
       
	 if(op->strucmp)
	   score=*(Potn+*(Icont+*m)); /* take the first block */
	 else
	   {
	     score=*(Potn+imod2vec(seq+*n,Icont+*m,Imod));
	     if(op->debug) printf(" site=%d res=%d score=%f \n",*m,*n,score);
	   }
       }
     if(op->gap_adr>=0) *gpen=*(Potn+imod2vec(&op->gap_adr,Icont+*m,Imod));
   }


 return(score);
 
}


/* write_contType calculates number of contacts according to their 
   types and puts them into both: file Ftype and vector Mat */
void write_contType(a_file *F,k_opt *op,int *Imod,int *contR,int *IJ,
		    int *start,int *seq,int *Lseq,char *R_name,int *Icont,
		    solv_sh *ssh,float *contRvdw,float *r12ij,float *r6ij)
{
 int i,j,k,m,alpha,beta,nbij,*IJ_loc;


 /* Icont=ssh->Ncont; Icont is local copy for first solv sh */
 if(op->ein || op->einm)
   for(i=0;i<*Lseq;++i)
     *(Icont+i)=0;
 
 if(op->evdw) 
   {
     nbij=op->n_type*(op->n_type+1)/2;
     for(i=0;i<op->n_par;++i)
       *(contRvdw+i)=0;
   }
 else
   {
     for(i=0;i<op->n_par;++i)
       *(contR+i)=0;
   }
 
 /* first skip everything up to the point where the threading starts */
 if(op->evdw)
   {
     for(i=0;i<*start;++i)
       {
	 if(*IJ==i) 
	   {
	     ++IJ; 
	     ++r12ij;
	     ++r6ij;
	   }
	 while(*IJ!=i+1)  
	   {
	     ++IJ; 
	     ++r12ij;
	     ++r6ij;
	   }
       }
   }
 else
   {
     for(i=0;i<*start;++i)
       {
	 if(*IJ==i) ++IJ;     
	 while(*IJ!=i+1) ++IJ; 
       } 
   }

 if(op->einm) IJ_loc=IJ; /* create a copy for einm */

 k=*start; 
 /* if threading is regular - fixed JUMP - one may calculate where to 
    start in IJ during next call to avoid the previous loop  */
 for(i=0;i<*Lseq-DELTA;++i)
   {
     alpha=*(seq+i);
     if(*IJ==k) 
       {
	 ++IJ; 
	 if(op->evdw) 
	   { 
	     ++r12ij;
	     ++r6ij;
	   }
       } 
     while(*IJ!=k+1)
       {
	 m=*IJ-*start;
	 if(m<*Lseq)
	   {
	     beta=*(seq+m); 
	     if(op->ein || op->einm)
	       {
		 ++*(Icont+i);
		 ++*(Icont+m);
	       }
	     else
	       {
		 if(op->evdw)
		   {
		     if(alpha>beta)
		       {
			 *(contRvdw+mat2vec(&beta,&alpha,&op->n_type))
			   +=*r12ij;
			 /* the Bij parameters come after Aij */
			 *(contRvdw+mat2vec(&beta,&alpha,&op->n_type)+nbij)
			   +=*r6ij;
		       }
		     else
		       {
			 *(contRvdw+mat2vec(&alpha,&beta,&op->n_type))
			   +=*r12ij;
			 *(contRvdw+mat2vec(&alpha,&beta,&op->n_type)+nbij)
			   +=*r6ij;
		       }
		   }
		 else
		   {
		     if(alpha>beta)
		       ++*(contR+mat2vec(&beta,&alpha,&op->n_type));
		     else
		       ++*(contR+mat2vec(&alpha,&beta,&op->n_type));
		   }
	       }
	   }
	 ++IJ; 
	 if(op->evdw) 
	   { 
	     ++r12ij;
	     ++r6ij;
	   }
       }
     ++k;
   }

 if(op->einm) get_einm(op,Imod,contR,IJ_loc,start,seq,Lseq,ssh->Ncont,
		       ssh->Mcont0,ssh->Mcont1,ssh->Mcont2);
	
 /* CHANGE */
 /*
 if(op->ein)
   {
     for(i=0;i<*Lseq;++i) 
     if(*(Icont+i)>0) ++*(contR+imod2vec(seq+i,Icont+i,Imod));
   }
 */
 if(op->ein)
   {
     if(op->gap_adr>=0) /* GAPs without neighbours will be accounted for */
       {	 
	 for(i=0;i<*Lseq;++i)
	   if((*(Icont+i)>0) || (op->gap_adr==*(seq+i))) 
	     ++*(contR+imod2vec(seq+i,Icont+i,Imod));
       }
     else

       {	 
	 for(i=0;i<*Lseq;++i)
	   if(*(Icont+i)>0) ++*(contR+imod2vec(seq+i,Icont+i,Imod));
	 /* fix to avoid e_i(0) - remove if in line above to have it
	    ++*(contT+imod2vec(&seq[i],Icont+i,Imod)); */
       }
   }
 /* CHANGE */  

 fprintf(F->type,"%s\n",R_name);

 if(op->evdw)
   {
     for(i=0;i<op->n_par;++i)
       {
	 fprintf(F->type,"%e ",*contRvdw);
	 ++contRvdw;
       }
   }
 else
   {
     for(i=0;i<op->n_par;++i)
       {
	 fprintf(F->type,"%d ",*contR);
	 ++contR;
       }
   }
 fprintf(F->type,"\n");

}

/* read_contType reads number of contacts stored according to their 
   types onto file Ftype and puts them into vector contR */
int read_contType(a_file *F,k_opt *op,int *contR,char *T_name,float *contRvdw)
{
 int i,j,k,eof;
 char Loc_name[MAXS];

 eof=fscanf(F->type,"%s",Loc_name); 
 if(op->debug) printf("\n Reading %s",Loc_name);
 if(eof==EOF) return(0);
 if(strcmp(Loc_name,T_name)!=0) 
   { printf("\n Mismatch %s %s ",Loc_name,T_name); return(0); } 

 if(op->evdw)
   {
     for(i=0;i<op->n_par;++i)
       {
	 fscanf(F->type,"%f",contRvdw);
	 ++contRvdw;
       }
   }
 else
   {
     for(i=0;i<op->n_par;++i)
       {
	 fscanf(F->type,"%d",contR);
	 ++contR;
       }
   }
 return(1);
}

/* this function gives contacts by type if second shell is used */
void get_einm(k_opt *op,int *Imod,int *contR,int *IJ_loc,int *start,
	      int *seq,int *Lseq,int *Ncont,int *Mcont0,int *Mcont1,
	      int *Mcont2)
{
 int i,j,k,m,alpha,beta,n_i,n_m;

 for(i=0;i<*Lseq;++i)
   {
     *(Mcont0+i)=0;
     *(Mcont1+i)=0;
     *(Mcont2+i)=0;
   }
 /* division into categories is manually fixed!!! */

 /* having all the neighbours in Ncont one can repeat loop and count Icount1,
    Icount2, Icount3 containing only neighbours of a given number of their
    own neighbours (second shell) */

 k=*start; 
 if(op->debug) printf("seq of length %d\n",*Lseq); 
 for(i=0;i<*Lseq-DELTA;++i)
   {
     alpha=seq[i];  
     if(op->debug) printf("%d\t\t",*(Ncont+i)); 
     if(*IJ_loc==k) ++IJ_loc;
     while(*IJ_loc!=k+1)
       {
	 m=*IJ_loc-*start;
	 if(m<*Lseq)
	   {
	     beta=*(seq+m); 
	     if(op->debug) printf("%d\t",*(Ncont+m)); 
	     n_i=*(Ncont+i);
	     n_m=*(Ncont+m);
	     if(n_m<3) ++*(Mcont0+i);
	     else 
	       {
		 if(n_m<7) ++*(Mcont1+i);
		 else ++*(Mcont2+i);
	       }
	     if(n_i<3) ++*(Mcont0+m);
	     else 
	       {
		 if(n_i<7) ++*(Mcont1+m);
		 else ++*(Mcont2+m);
	       }
	     
	   }
	 ++IJ_loc; 
       }
     ++k; 
     if(op->debug) printf("\n"); 
   } /* now second shell is characterized in Mcont0,1,2 */

 /* assumption: n_par=3(sec. shell)*n_type(first shell) !!! fixed */
 for(i=0;i<*Lseq;++i)
   {
     j=0; n_i=*(Ncont+i);
     /* for n_type=4 use the next line */
     /* if(n_i>2 && n_i<6) j+=1; if(n_i>5 && n_i<9) j+=2; if(n_i>8) j+=3; */
     /* currently n_type=5 and thus n_par(per ama)=15 is assumed - fixed !!! */
     if(n_i>2) 
       {
	 ++j;
	 if(n_i>4) 
	   {
	     ++j;
	     if(n_i>6) 
	       {
		 ++j;
		 if(n_i>8) ++j;
	       }	     
	   }
       }
     j*=3;
     /* this loop could be avoided - Imod */
     for(k=0;k<*(seq+i);++k) j+=*(Imod+k); 
     if(op->debug)
       {
	 printf("#%d n_i %d I0i %d \n",i,n_i,*(Mcont0+i));
	 printf(" addresses I0i %d I1i %d I2i %d \n",j,j+1,j+2); 
       }
     /* if GAP residue is used and a site of no-neighbours is occup by GAP */
     if(op->gap_adr>=0) 
       {
	 if((n_i==0) && (op->gap_adr==*(seq+i))) 
	   ++*(contR+j+*(Imod+op->gap_adr)-1);
	 else
	   {
	     *(contR+j)+=*(Mcont0+i);
	     *(contR+j+1)+=*(Mcont1+i);
	     *(contR+j+2)+=*(Mcont2+i);
	   }
       }
     else
       {
	 *(contR+j)+=*(Mcont0+i);
	 *(contR+j+1)+=*(Mcont1+i);
	 *(contR+j+2)+=*(Mcont2+i);
       }
     
   }
 
}

/* this function gives environment if second shell is used */
void get_env_einm(k_opt *op,int *Imod,int *contR,int *IJ_loc,int *start,
	      int *seq,int *Lstru,int *Ncont,int *Mcont0,int *Mcont1,
	      int *Mcont2)
{
 int i,j,k,m,alpha,beta,n_i,n_m;

 for(i=0;i<*Lstru;++i)
   {
     *(Mcont0+i)=0;
     *(Mcont1+i)=0;
     *(Mcont2+i)=0;
   }
 /* division into categories is manually fixed!!! */

 /* having all the neighbours in Ncont one can repeat loop and count Icount1,
    Icount2, Icount3 containing only neighbours of a given number of their
    own neighbours (second shell) */

 k=*start; 
 if(op->debug) printf("stru of length %d\n",*Lstru); 
 /* get N,Mcont vectors for the whole structure! */
 for(i=0;i<*Lstru-DELTA;++i)
   {
     /* alpha=seq[i]; */ 
     if(op->debug) printf("%d(%d)\t\t",i,*(Ncont+i)); 
     if(*IJ_loc==k) ++IJ_loc;
     while(*IJ_loc!=k+1)
       {
	 m=*IJ_loc-*start;
	 if(m<*Lstru)
	   {
	     /* beta=*(seq+m); */
	     if(op->debug) printf("%d(%d)\t",m,*(Ncont+m)); 
	     n_i=*(Ncont+i);
	     n_m=*(Ncont+m);
	     if(n_m<3) ++*(Mcont0+i);
	     else 
	       {
		 if(n_m<7) ++*(Mcont1+i);
		 else ++*(Mcont2+i);
	       }
	     if(n_i<3) ++*(Mcont0+m);
	     else 
	       {
		 if(n_i<7) ++*(Mcont1+m);
		 else ++*(Mcont2+m);
	       }
	     
	   }
	 ++IJ_loc; 
       }
     ++k; 
     if(op->debug) printf("\n"); 
   } /* now second shell is characterized in Mcont0,1,2 */


}

float energy(k_opt *op,int *contR,float *Potn,float *contRvdw)
{
 int i;
 float energ;

 energ=0.0;

 if(op->evdw)
   { 
   for(i=0;i<op->n_par;++i)
       {
	 if(*(contRvdw+i)!=0.0)
	   energ+=*(contRvdw+i)*(*Potn);
	 ++Potn;
       }
   }
 else
   {
     for(i=0;i<op->n_par;++i)
       {
	 if(*(contR+i)!=0)
	   energ+=(float)*(contR+i)*(*Potn);
	 ++Potn;
       }
   }
 return(energ);
}
 

/* mat2vec returns the vector address of an elem of a symmetric matrix */
int mat2vec(int *i,int *j,int *n_type)
{
  return(*j+(int)(*i*0.5*(*n_type*2-*i-1))); 
  /* upper triangle is rolled to a vector */
}

/* imod2vec returns the vector address of a parameter i */
int imod2vec(int *i,int *n,int *Imod)
{
  int j,k;
  /* here *n is the number of neighbours of yet another residue of type *i, 
   *(Imod+*i) is the number of parameters for residues of type *i */

  j=0;
  for(k=0;k<*i;++k) j+=*(Imod+k);
  if(*n>*(Imod+*i)-1) j+=*(Imod+*i)-1;
  else j+=*n;
  return(j);
}

/* skip reads seq and contType files */
void skip(a_file *F,k_opt *op,int *seq,int *i_prot,int *contR,float *contRvdw,
	  int *Alph,int *Adr)
{
  char T_name[MAXS];
  int Lseq;
  
  if(op->exam) read_seq(F->exam,&Lseq,T_name,seq,i_prot,&op->n_type,Alph,Adr);
  else
    {
      read_seq(F->seq,&Lseq,T_name,seq,i_prot,&op->n_type,Alph,Adr);
      read_contType(F,op,contR,T_name,contRvdw);
    }
}

