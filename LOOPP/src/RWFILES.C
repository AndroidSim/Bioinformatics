/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000 v. 2.000                                      */
/*                                                                */
/*        READ SEQ, XYZ, POT FILES, WRITE LP CONSTRAINTS          */
 
#include "loopp.h"
   

/* read the potential (scoring) function */
void read_pot(k_opt *op, float *Potn, a_file *F, int *Imod)
{
 int i,j,k,l,m;
 float any;

 if(op->mat && op->evdw)
   {
     printf("Sorry - matrix format not supported  for evdw - trying vector \n");
     op->mat=0;
   } 
 if(op->mat) /* beware of the fixed order of amino acids  */
   {
     if(op->eij)
       {
	 
	 for(i=0;i<op->n_type;++i)
	   for(j=0;j<op->n_type;++j)
	     {
	       fscanf(F->pot,"%f",&any);
	       if(j>=i)  /* upper triangle */
		 {
		   *Potn=op->resc_pot*any;
		   ++Potn; 
		 } 
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
		 if(*(Imod+j)>i) 
		   {
		     fscanf(F->pot,"%f",&any);
		     *(Potn+m+i)=op->resc_pot*any;
		   }
	       }
	   }     
       }
   }
 else /* if vector */
   {    
     for(i=0;i<op->n_par;++i)
       {
	 fscanf(F->pot,"%f",&any);
	 *Potn=op->resc_pot*any;
	 ++Potn; 
       }
   }
 
} 


void write_ineq(a_file *F,k_opt *op,float *Potn,int *contR,int *contT,
		unsigned long int *index,float *contRvdw,float *contTvdw,
		k_mps *var_mps[],k_mps *tail_mps[])
{
 int i,j,n;
 k_mps *ptr;
 float rhs;
 

 /* produce either f_dcon or f_mps - it makes little sense to get both  */
 if(op->wDc)
   {
     if(op->evdw) /* get pointers to vdw */
       {
	 for(i=0;i<op->n_par;++i)
	   {
	     fprintf(F->out2,"%e ",*contTvdw-*contRvdw);
	     ++contRvdw; ++contTvdw;
	   }
       }
     else
       {
	 if(op->freeze)
	   {
	     rhs=0.0;
	     for(i=0;i<op->nfrozen;++i)
	       {
		 rhs-=*Potn*(*contT-*contR);
		 ++contR; ++contT; ++Potn;
	       }
	     fprintf(F->rhs,"%f \n",rhs);
	     for(i=op->nfrozen;i<op->n_par;++i)
	       {
		 fprintf(F->out2,"%d ",*contT-*contR);
		 ++contR; ++contT;
	       }
	   }
	 else
	   {
	     for(i=0;i<op->n_par;++i)
	       {
		 fprintf(F->out2,"%d ",*contT-*contR);
		 ++contR; ++contT;
	       }
	   }
       }
     fprintf(F->out2,"\n");
   }
 else if(op->mps) /* storing inequlities as linked lists */
   {
     if(op->evdw) 
       {
	 printf("Not built in - use an external converter ... \n");
	 exit(0);
       }
    
     for(i=0;i<op->n_par;++i)
       {
	 n=*contT-*contR;
	 if(n!=0)
	   {
	     if(var_mps[i]==NULL)
	       {
		 var_mps[i]=(k_mps *)malloc(sizeof(k_mps));
		 var_mps[i]->value=n;
		 var_mps[i]->row=*index-1;
		 var_mps[i]->down=NULL;
		 tail_mps[i]=var_mps[i];
	       }
	     else
	       {
		 ptr=tail_mps[i];
		 ptr->down=(k_mps *)malloc(sizeof(k_mps));
		 if(ptr->down==NULL)
		   {
		     printf("Not enough memory - I stop ...\n");
		     exit(1);
		   }
		 ptr->down->value=n;
		 ptr->down->row=*index-1;
		 ptr->down->down=NULL;
		 tail_mps[i]=ptr->down;
	       }
	   }
	 ++contR; ++contT;
       }
   }
 
 

}

/* read_seq reads sequence of struc T_name and returns 0 when no more
   sequences left to read */
int read_seq(k_opt *op,FILE *Fseq, int *Lseq, char *T_name, int *seq, int *i_prot,
	     int *n_type, int *Alph, int *Adr)
{
 int i,j,eof;
 char amino[4];

 eof=fscanf(Fseq,"%s\n%d",T_name,Lseq);
 if(eof==EOF) return(0);
 if (op->fngrps==1 && op->strucmp==1)
 {
	 for(i=0;i<*Lseq;++i)
	 {
		fscanf(Fseq,"%s",amino);
		/* transform string to generic integer repr. */
		j=seq2dig(amino);
		if(j<0) read_err(T_name);
		/* assumption: sequences longer than allowed to be ignored anyway */
		if(i<op->max_length) *(seq+i)=j;
	 }
 }
 else
 {
	 for(i=0;i<*Lseq;++i)
	 {
		fscanf(Fseq,"%s",amino);
		/* transform string to generic integer repr. and then to current alph. */
		j=seq2dig(amino);
		if(j<0) read_err(T_name);
		*(seq+i)=generic2adr(&j,n_type,Alph,Adr);
	 }
 }
 ++(*i_prot);
 return(1);
}

/*int read_seq(FILE *Fseq, int *Lseq, char *T_name, int *seq, int *i_prot,
	     int *n_type, int *Alph, int *Adr)
{
 int i,j,eof;
 char amino[4];

 eof=fscanf(Fseq,"%s\n%d",T_name,Lseq);
 if(eof==EOF) return(0);
 for(i=0;i<*Lseq;++i)
   {
    fscanf(Fseq,"%s",amino);
    /* transform string to generic integer repr. and then to current alph. */
/*    j=seq2dig(amino);
    if(j<0) read_err(T_name);
    *(seq+i)=generic2adr(&j,n_type,Alph,Adr);
   }
 ++(*i_prot);
 return(1);
} */

/* read_seq reads sequences assuming that they are defined in term of
   generic alphabet */
int read_seq_generic(FILE *Fseq, int *Lseq, char *T_name, int *seq, int *i_prot,
		     k_opt *op)
{
 int i,j,eof;
 char amino[4];

 eof=fscanf(Fseq,"%s\n%d",T_name,Lseq);
 if(eof==EOF) return(0);
 for(i=0;i<*Lseq;++i)
   {
    fscanf(Fseq,"%s",amino);
    /* transform string to generic integer repr. */
    j=seq2dig(amino);
    if(j<0) read_err(T_name);
    /* assumption: sequences longer than allowed to be ignored anyway */
    if(i<op->max_length) *(seq+i)=j;
   }
 ++(*i_prot);
 return(1);
}


int read_coor(FILE *Fin, k_opt *op, int *Lstru, char *R_name,
	      float *X, float *Y, float *Z)
{
 int i,eof;
 float f;

 /* for consistency with Dror's database the second 3-col block is taken
    for the representation of side chains - use op->col to change it */
 /* first three coordinates can be used for C_alpha representation */
 /* another 3-column block could be used for C_beta representation */
 
 eof=fscanf(Fin,"%s\n%d",R_name,Lstru);
 if(eof==EOF) return(0);

 if(op->n_col==2)
   for(i=0;i<*Lstru;++i)
     {
       if(op->col==2) fscanf(Fin,"%f%f%f%f%f%f",&f,&f,&f,X+i,Y+i,Z+i);
       if(op->col==1) fscanf(Fin,"%f%f%f%f%f%f",X+i,Y+i,Z+i,&f,&f,&f);
     }

 if(op->n_col==1)
   for(i=0;i<*Lstru;++i)
     fscanf(Fin,"%f%f%f",X+i,Y+i,Z+i);

 if(op->n_col==3)
   for(i=0;i<*Lstru;++i)
     {
       if(op->col==2) 
	 fscanf(Fin,"%f%f%f%f%f%f%f%f%f",&f,&f,&f,X+i,Y+i,Z+i,&f,&f,&f);
       if(op->col==1) 
	 fscanf(Fin,"%f%f%f%f%f%f%f%f%f",X+i,Y+i,Z+i,&f,&f,&f,&f,&f,&f);
       if(op->col==3) 
	 fscanf(Fin,"%f%f%f%f%f%f%f%f%f",&f,&f,&f,&f,&f,&f,X+i,Y+i,Z+i);
     }

 return(1);
}

int read_coor_all3(FILE *Fin, k_opt *op, int *Lstru, char *R_name,
		   float *X, float *Y, float *Z, float *Xa, float *Ya, 
		   float *Za, float *Xb, float *Yb, float *Zb)
{
 int i,eof;
 float f;

 /* for consistency with Dror's database the second 3-col block is taken
    for the representation of side chains - use op->col to change it */
 /* first three coordinates can be used for C_alpha representation */
 /* another 3-column block could be used for C_beta representation */
 
 eof=fscanf(Fin,"%s\n%d",R_name,Lstru);
 if(eof==EOF) return(0);

 if(op->n_col==2)
   for(i=0;i<*Lstru;++i) fscanf(Fin,"%f%f%f%f%f%f",Xa+i,Ya+i,Za+i,X+i,Y+i,Z+i);

 if(op->n_col==1)
   for(i=0;i<*Lstru;++i) fscanf(Fin,"%f%f%f",Xa+i,Ya+i,Za+i);

 if(op->n_col==3)
   for(i=0;i<*Lstru;++i)
     fscanf(Fin,"%f%f%f%f%f%f%f%f%f",Xa+i,Ya+i,Za+i,X+i,Y+i,Z+i,Xb+i,Yb+i,Zb+i);

 return(1);
}


/* write_mps transforms diff in cont vector into MPS format */
void write_mps(a_file *F,k_opt *op,unsigned long int *index,k_mps *var_mps[])
{
 k_mps *ptr;
 int i,j;

 /* only primal problem is written - dual does not need so much memory 
    but may be tricky with new constraints and obj functions */
 /* this function is based on a function written by Dror Tobi */

 printf(" Writting LP constraints in MPS format to file %s \n",F->f_mps);
 printf(" It may take some time and some disk space ... \n");

 fprintf(F->mps,"NAME          unnamed\n");
 
 /* writing rows */
 fprintf(F->mps,"ROWS\n");
 fprintf(F->mps," N  r_0\n");
 for(i=0;i<*index;++i) fprintf(F->mps," G  r_%d\n",i+1);
 /* writing columns */
 fprintf(F->mps,"COLUMNS\n");
 for(i=0;i<op->n_par;++i)
   {
    ptr=*(var_mps+i);
    if(ptr!=NULL)
      {
       while(ptr->down != NULL)
       {
	  fprintf(F->mps,"    var_%d",i+1); 
	  space(F->mps,i+1,4);
	  fprintf(F->mps,"  r_%d",ptr->row+1);
	  space(F->mps,ptr->row+1,6);
	  fprintf(F->mps,"  %d\n",ptr->value);
	  ptr=ptr->down;
       }
       fprintf(F->mps,"    var_%d",i+1); 
       space(F->mps,i+1,4);
       fprintf(F->mps,"  r_%d",ptr->row+1);
       space(F->mps,ptr->row+1,6);
       fprintf(F->mps,"  %d\n",ptr->value);
       ptr=ptr->down;
      }
   }
 /*  writing vector b */
 fprintf(F->mps,"RHS\n");
 for(i=1;i<=*index;++i)
   {
    fprintf(F->mps,"    RHS       ");
    fprintf(F->mps,"r_%d",i);
    space(F->mps,i,6);
    fprintf(F->mps,"  %g\n",epsilon);
   }

 /* writing bounds */
 fprintf(F->mps,"BOUNDS\n");
 for(i=1;i<=op->n_par;++i)
   {
    fprintf(F->mps," UP BND       var_%d",i);
    space(F->mps,i,4);
    fprintf(F->mps,"  %d\n",10);
    fprintf(F->mps," LO BND       var_%d",i);
    space(F->mps,i,4);
    fprintf(F->mps,"  %d\n",-10);

   }
 fprintf(F->mps,"ENDATA\n");
}

/* just small addition to write_mps */
void space(FILE *out,int v,int n)
{
 int i,a;

 a=1;
 while(v>9)
   {
    v=v/10;
    ++a;
   }
 for(i=a;i<n;++i) fprintf(out," ");

}


void read_err(char *name)
{
  printf("Error while reading %s\n",name);
  exit(1);  
}

void open_err(char *name)
{
  printf("Sorry - cannot open file %s\n",name);
  exit(1);  
}

