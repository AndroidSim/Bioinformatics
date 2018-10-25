/* modification to the energy calculation of the fea (frozen
   environment approximation) in the threading section of the
   program LOOPP written by Jarek Meller and Ron Elber

  this section written by Andrew Smith */

#include "loopp.h"


/*struct fea_mod
{
	int *env_eij;
   float *r12_env;
   float *r6_env;
}; 
typedef struct fea_mod f_mod;
f_mod fea_data;
f_mod *fea; */

/* function prototypes */

void get_env_eij(k_opt *op,f_mod *fea,int *Imod,int *contT,int *IJ,int *nat_Seq,
				int *Lnatseq,float *contTvdw,float *r12ij,float *r6ij);
float get_energy_eij(k_opt *op,f_mod *fea,int *Imod,int *contT,int *start,int *j_beta,
					 int *qseq_i_alpha,int *nat_Seq,int *Lnatseq,
                     float *Potn,float *contTvdw,float *r12ij,float *r6ij,
                     float *gpen);


void get_env_eij(k_opt *op,f_mod *fea,int *Imod,int *contT,int *IJ,int *nat_Seq,
				int *Lnatseq,float *contTvdw,float *r12ij,float *r6ij)
{
	int i,j,k,alpha,beta,*env_eij /* **test */;
   float *r12_env,*r6_env;

/*	fea=&fea_data; */ 
   env_eij=fea->env_eij;
   r12_env=fea->r12_env;
   r6_env=fea->r6_env;

/*	test=malloc((*Lnatseq)*sizeof(int));
	if (test==NULL)
	{
		printf("unable to allocate memory for test");
		exit(1);
	}

	for (i=0; i<op->n_type; ++i)
	{
		*test=malloc(op->n_type*sizeof(int));
		if (test==NULL)
		{
			printf("unable to allocate memory for test");
			exit(1);
		}
	} */

   /* allocate memory for env_eij,r12_env,r6_env data structures */

	fea->env_eij=malloc((*Lnatseq*(op->n_type))*sizeof(int));
	if (!fea->env_eij)
	{
		printf("unable to allocate memory for env_eij\n");
		exit(1);
	}
	if (op->evdw)
	{
		fea->r12_env=malloc((*Lnatseq*(op->n_type))*sizeof(float));
		if (!fea->r12_env)
		{
			printf("unable to allocate memory for r12_env\n");
			exit(1);
		}
		fea->r6_env=malloc((*Lnatseq*(op->n_type))*sizeof(float));
		if (!fea->r6_env)
		{
   			printf("unable to allocate memory for r6_env\n");
			exit(1);
		}
	}


   /* initialize env_eij to zero contacts for each site of template structure */

	for (i=0; i<*Lnatseq; ++i)
	{
		for (j=0; j<op->n_type; ++j)
		{
			fea->env_eij[i*op->n_type+j]=0;
		/*	test[i][j]=0; */
		}
	}

   /* begin loop over IJ to get all the contacts and store them in env_eij */

	for (k=0; k<*Lnatseq; ++k)
	{
		if (*IJ==k)
		{
			++IJ; 
			if (op->evdw) 
			{ 
				++r12ij;
				++r6ij;
			}
		}
		
		while (*IJ!=(k+1))
		{
			alpha=*(nat_Seq+k);
			beta=*(nat_Seq+*IJ);
			fea->env_eij[k*(op->n_type)+beta]++;
			fea->env_eij[*IJ*(op->n_type)+alpha]++;
			++IJ;
			if (op->evdw)
			{
				fea->r12_env[k*(op->n_type)+beta]=*r12ij;
				fea->r12_env[*IJ*(op->n_type)+alpha]=*r12ij;
				fea->r6_env[k*(op->n_type)+beta]=*r6ij;
				fea->r6_env[*IJ*(op->n_type)+alpha]=*r6ij;
				++r12ij;
				++r6ij;
			}
		}
	}

/*   if (op->evdw)
   {
   	return (fea->r12_env);
   }
   else
   {
   	return (fea->env_eij);
   }
/*	return(0.5*energy(op,contT,Potn,contTvdw)); */
}

float get_energy_eij(k_opt *op,f_mod *fea,int *Imod,int *contT,int *start,int *j_beta,
					 int *qseq_i_alpha,int *nat_Seq,int *Lnatseq,
                     float *Potn,float *contTvdw,float *r12ij,float *r6ij,
                     float *gpen)
{
	int i,j,k,n_cont,nbij,par_ind,*env_eij,alpha;
	float energ,*r12_env,*r6_env;

/*   fea=&fea_data; */
   env_eij=fea->env_eij;
   r12_env=fea->r12_env;
   r6_env=fea->r6_env;

/*	initialize energy (energ), number of contacts for gaps (n_cont), and
    the number of contacts of a specific type (contT and contTvdw) to 0 */

	energ=0.0;
	n_cont=0;

   if (op->evdw)
   {
      nbij=op->n_type*(op->n_type+1)/2;
      for (i=0;i<op->n_par;++i)
      {
      	*(contTvdw+i)=0.0;
      }
   }
	else
   {
   	for (i=0;i<op->n_par;++i)
      {
      	*(contT+i)=0;
      }
   }

	/* calculate the number of contacts at a position of the template */

	for (j=0; j<op->n_type; ++j)
	{
		n_cont+=env_eij[*j_beta*(op->n_type)+j];
	}

	/* for pairwise potentials gap penalty is set to be proportional to 
       n_cont+1 without referring to explicit GAP residue in the Potn */

	*gpen=(n_cont+1)*op->gap_pen;

	/* use passed sequence identity (now alpha) from dpt to calculate 
	   contT and contTvdw for energy calculation of putting alpha into 
	   structural site j_beta */

	alpha=*qseq_i_alpha;

/*   if(op->evdw)
   {
   	if(alpha>beta)
      {
      	*(contTvdw+mat2vec(&beta,&alpha,&op->n_type))
         +=*r12ij;
         /* the Bij parameters come after Aij */
/*         *(contTvdw+mat2vec(&beta,&alpha,&op->n_type)+nbij)
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
   } */

	/* calculate the data structures contT and contTvdw */
	/* par_ind is used as a parameter index for contT and contTvdw to 
	   place the contacting or interacting aa types */

   if (op->evdw)
   {
   	/* calculate contTvdw using r12ij and r6ij */
      for (i=1; i<alpha; ++i)
      {
      	par_ind+=i;
      }
      for (i=0; i<alpha; ++i)
      {
      	contTvdw[par_ind+i]=fea->r12_env[*j_beta*(op->n_type)+i];
        contTvdw[par_ind+i]-=fea->r6_env[*j_beta*(op->n_type)+i];
      }
      k=alpha;
      while (k<=(op->n_type))
      {
      	 par_ind=par_ind+k;
         contTvdw[par_ind+(alpha-1)]=fea->r12_env[*j_beta*(op->n_type)+(k+1)];
         contTvdw[par_ind+(alpha-1)]-=fea->r6_env[*j_beta*(op->n_type)+(k+1)];
         k++;
      }
   }
   else
   {
   	/* calculate contT for non-continuous or contact potential */
	   par_ind=0;
   	  for (i=1; i<alpha; ++i)
      {
      	par_ind+=i;
      }
      for (i=0; i<alpha; ++i)
      {
      	contT[par_ind+i]=fea->env_eij[*j_beta*(op->n_type)+i];
      }
      k=alpha;
      while (k<=(op->n_type))
      {
      	 par_ind=par_ind+k;
         contT[par_ind+(alpha-1)]=fea->env_eij[*j_beta*(op->n_type)+(k+1)];
         k++;
      }
   }

   if(op->evdw)
	{
		for(i=0;i<op->n_par;++i)
		{
			if(*(contTvdw+i)!=0.0)
			{
				energ+=*(contTvdw+i)*(*Potn);
				++Potn;
			}
		}
	}
	else
	{
		for(i=0;i<op->n_par;++i)
		{
			if(*(contT+i)!=0)
			{
				energ+=(float)*(contT+i)*(*Potn);
				++Potn;
			}
		}
	}

	return(energ);
}