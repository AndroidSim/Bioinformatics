/* modification to the energy calculation of the fea (frozen
   environment approximation) in the threading section of the
   program LOOPP written by Jarek Meller and Ron Elber

  this section written by Andrew Smith */

#include "loopp.h"


void get_env_eij(k_opt *op,f_mod *fea,int *Imod,int *IJ,int *nat_Seq,
				int *Lnatseq,int *nop,float *r12ij,float *r6ij,float *s_p_env)
{
	int i,j,k,x,y,alpha,beta,*env_eij,*tn_cont,*flag,n_flags,n_type;
    float *r12_env,*r6_env;
	
/*	fps=&fps_data; */

   /* allocate memory for env_eij,r12_env,r6_env data structures */

	if (op->fngrps && op->strucmp)
	{
		n_type=fps->fpsn_type;
	}
	else
	{
		n_type=op->n_type;
	}

	if (op->eij || op->strucmp)
	{
		fea->env_eij=malloc((*Lnatseq*(n_type))*sizeof(int));
		if (!fea->env_eij)
		{
			printf("unable to allocate memory for env_eij\n");
			exit(1);
		}

		fea->flag=malloc((*Lnatseq*(n_type))*sizeof(int));
		if (!fea->flag)
		{
			printf("unable to allocate memory for flag\n");
			exit(1);
		}
	}
	if (op->evdw)
	{
		fea->r12_env=malloc((*Lnatseq*(n_type))*sizeof(float));
		if (!fea->r12_env)
		{
			printf("unable to allocate memory for r12_env\n");
			exit(1);
		}

		fea->r6_env=malloc((*Lnatseq*(n_type))*sizeof(float));
		if (!fea->r6_env)
		{
   			printf("unable to allocate memory for r6_env\n");
			exit(1);
		}
	}

	fea->tn_cont=malloc(*Lnatseq*sizeof(int));
	if (!fea->tn_cont)
	{
		printf("unable to allocate memory for tn_cont\n");
		exit(1);
	}
	
	if (op->eij || op->strucmp)
	{
		env_eij=fea->env_eij;
		flag=fea->flag;
	}
	if (op->evdw)
	{
		r12_env=fea->r12_env;
		r6_env=fea->r6_env;
	}
	tn_cont=fea->tn_cont;

   /* initialize env_eij to zero contacts for each site of template structure */

	for (i=0; i<*Lnatseq; ++i)
	{
		for (j=0; j<n_type; ++j)
		{
			if (op->eij || op->strucmp)
			{
				*(env_eij+i*n_type+j)=0;
				*(flag+i*n_type+j)=0;
				(*nop)++;
			}
			if (op->evdw)
			{
				*(r12_env+i*n_type+j)=0;
				*(r6_env+i*n_type+j)=0;
			}
		}
		*(tn_cont+i)=0;
	}
	n_flags=0;

   /* begin loop over IJ to get all the contacts and store them in env_eij */
	
	for (k=0; k<*Lnatseq-DELTA; ++k)
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
			x=k*(n_type)+beta;
			y=*IJ*(n_type)+alpha;
			if (op->eij || op->strucmp)
			{
				env_eij[x]++;
				env_eij[y]++;
				if (flag[x]==0)
				{
					flag[x]=x;
					n_flags++;
				}
				if (flag[y]==0)
				{
					flag[y]=y;
					n_flags++;
				}
				tn_cont[k]++;
				tn_cont[*IJ]++;
			}
			if (op->evdw)
			{
				r12_env[x]+=*r12ij;
				r12_env[y]+=*r12ij;
				r6_env[x]+=*r6ij;
				r6_env[y]+=*r6ij;
				tn_cont[k]++;
				tn_cont[*IJ]++;
				++r12ij;
				++r6ij;
			}
			++IJ;
			(*nop)++;
		}
	}

	*s_p_env+=(float)n_flags/(*Lnatseq*(n_type));

}

float get_energy_eij(k_opt *op,f_mod *fea,int *Imod,int *IJ,int *contT,
					int *start,int *j_beta,int *qseq_i_alpha,int *nat_Seq,
                     int *Lnatseq,float *Potn,float *contTvdw,float *r12ij,
                     float *r6ij,float *gpen,int *nop,int *nope/*,double *e_time*/)
{
	int i,j,n_cont,nbij,par_ind,alpha,*env_eij,*tn_cont;
	float energ,*r12_env,*r6_env;

/*	time_t time1,time2;

	time1=time(NULL);*/

   env_eij=fea->env_eij;
   tn_cont=fea->tn_cont;
   r12_env=fea->r12_env;
   r6_env=fea->r6_env; 

/*	initialize energy (energ), number of contacts for gaps (n_cont), and
    the number of contacts of a specific type (contT and contTvdw) to 0 */

	energ=0.0;
	n_cont=0;
	par_ind=0;

   if (op->evdw)
   {
      nbij=op->n_type*(op->n_type+1)/2;
   }

	/* calculate the number of contacts at a position of the template */

	n_cont=tn_cont[*j_beta];

	/* for pairwise potentials gap penalty is set to be proportional to 
       n_cont+1 without referring to explicit GAP residue in the Potn */

	*gpen=(n_cont+1)*op->gap_pen;

	/* use passed sequence identity (now alpha) from dpt to calculate 
	   contT and contTvdw for energy calculation of putting alpha into 
	   structural site j_beta */

	alpha=*qseq_i_alpha;

 	/* calculate the data structures contT and contTvdw */
	/* par_ind is used as a parameter index for contT and contTvdw to 
	   place the contacting or interacting aa types */

   if (op->evdw)
   {
   	/* calculate contTvdw using r12ij and r6ij */

/*	  for (i=0; i<alpha; ++i)
	  {
		  *(contTvdw+mat2vec(&i,&alpha,&op->n_type))+=*(r12_env+*j_beta*(op->n_type)+i);
		  /* the Bij parameters come after Aij */
/*		  *(contTvdw+mat2vec(&i,&alpha,&op->n_type)+nbij)+=*(r6_env+*j_beta*(op->n_type)+i);
	  }
	  for (j=alpha; j<op->n_type; ++j)
	  {
		  *(contTvdw+mat2vec(&alpha,&j,&op->n_type))+=*(r12_env+*j_beta*(op->n_type)+j);
		  *(contTvdw+mat2vec(&alpha,&j,&op->n_type)+nbij)+=*(r6_env+*j_beta*(op->n_type)+j);
	  }*/

	   for (i=0; i<alpha; ++i)
	   {
		   energ+=*(r12_env+*j_beta*(op->n_type)+i)*(*(Potn+mat2vec(&i,&alpha,&op->n_type)));
		   energ+=*(r6_env+*j_beta*(op->n_type)+i)*(*(Potn+mat2vec(&i,&alpha,&op->n_type)+nbij));
	   }
	   for (j=alpha; j<op->n_type; ++j)
	   {
		   energ+=*(r12_env+*j_beta*(op->n_type)+j)*(*(Potn+mat2vec(&alpha,&j,&op->n_type)));
		   energ+=*(r6_env+*j_beta*(op->n_type)+j)*(*(Potn+mat2vec(&alpha,&j,&op->n_type)+nbij));
	   }
   }
   else
   {
   	/* calculate contT for non-continuous or contact potential */	

	   for (i=0; i<alpha; ++i)
	   {
		   energ+=(float)*(env_eij+*j_beta*(op->n_type)+i)*(*(Potn+mat2vec(&i,&alpha,&op->n_type)));
	   }
	   for (j=alpha; j<op->n_type; ++j)
	   {
		   energ+=(float)*(env_eij+*j_beta*(op->n_type)+j)*(*(Potn+mat2vec(&alpha,&j,&op->n_type)));
	   }
   } 

	return((float)0.5*energ); 

/*   time2=time(NULL);
   *e_time=difftime(time2,time1);*/

}




