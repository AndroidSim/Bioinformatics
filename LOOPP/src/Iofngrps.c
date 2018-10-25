
/* iofngrps.c is a file for reading and writing the fingerprint files
	written by Andrew Smith */

#include "loopp.h"

void write_fngrps(a_file *F,k_opt *op,f_mod *fea,fgrprnt *fps,int *IJ,
					int *seq,int *Lseq,int *Lstru,char *R_name,float *Potn,int quordata)
{
	int i,j,k,*env_eij,*tn_cont,*flag;
	int n_cont,nbij,n_type;
	float *r12_env,*r6_env,*x2x_Potn,*temp_pot;
	double *dev_pij,energ,sum_prob,*s_energy,*profile,*norm_pj,m,n;
	double sum_exp,*entropy;
	
	fps=&fps_data;
	if (op->strucmp)
	{
		n_type=fps->fpsn_type;
	}
	else
	{
		n_type=op->n_type;
	}

	/* allocate memory for s(i,j) matrix = *s_energy */

	if (*Lseq!=*Lstru)
	{
		printf("Lseq != Lstru, Lseq = %d, Lstru = %d \n",*Lseq,*Lstru);
		exit(1);
	}

	fps->s_energy=malloc((*Lstru*(n_type))*sizeof(double));
	if (!fps->s_energy)
	{
		printf("unable to allocate memory for s_energy\n");
		exit(1);
	}
	
	fps->profile=malloc((*Lstru*(n_type))*sizeof(double));
	if (!fps->profile)
	{
		printf("unable to allocate memory for profile\n");
		exit(1);
	}

	dev_pij=malloc((*Lstru*(n_type))*sizeof(double));
	if (!dev_pij)
	{
		printf("unable to allocate memory for dev_pij\n");
		exit(1);
	}

	norm_pj=malloc(*Lstru*sizeof(double));
	if (!norm_pj)
	{
		printf("unable to allocate memory for norm_pj\n");
		exit(1);
	}

	entropy=malloc(*Lstru*sizeof(double));
	if (!entropy)
	{
		printf("unable to allocate memory for entropy\n");
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
	s_energy=fps->s_energy;
	profile=fps->profile;

	if (op->evdw)
	{
		nbij=op->n_type*(op->n_type+1)/2;
	}

	if (op->strucmp)
	{
		x2x_Potn=fps->x2x_Potn; 
		/* switch pointers between stru2stru or threading potential and fingerprnt potential */
		temp_pot=Potn;
		Potn=x2x_Potn;
		x2x_Potn=temp_pot; 
	}

/*	for (i=0; i<*Lstru; ++i)
	{
		for (j=0; j<n_type; ++j)
		{
			if (j<n_type-1)
			{
				printf("%d",env_eij[i*n_type+j]);
			}
			else
			{
				printf("%d\n",env_eij[i*n_type+j]);
			}
		}
	} */

	/* create environment contact types for both eij and evdw using
	   frozen environment approx. environment data in structure f_mod 
	   pointed to by fea */

	/* initialize variables */

	n_cont=0;

	for (i=0; i<*Lstru; ++i)
	{
		for (j=0; j<n_type; ++j)
		{
			*(s_energy+i*(n_type)+j)=(double)0.0;
		}
	}

	for (i=0; i<*Lstru-1; ++i)
	{
		for (j=0; j<n_type; ++j)
		{
			*(profile+i*(n_type)+j)=(double)0.0;
		}
	}

	/* get_env_eij() was here */

/* calculate the number of contacts at a position of the template */

/*	n_cont=tn_cont[*j_beta]; */

	if (op->evdw)
	{
   	/* calculate energy for elements of s(i,j) using r12ij and r6ij */

		for (i=0; i<*Lstru; ++i)
		{
			for (j=0; j<op->n_type; ++j)
			{
				energ=(double)0.0;
				for (k=0; k<op->n_type; ++k)
				{
					if (k>j)
					{
						energ+=(*(r12_env+i*(n_type)+k)*(*(Potn+mat2vec(&j,&k,&n_type))));
						energ+=(*(r6_env+i*(n_type)+k)*(*(Potn+mat2vec(&j,&k,&n_type)+nbij)));
					}
					else
					{
						energ+=(*(r12_env+i*(n_type)+k)*(*(Potn+mat2vec(&k,&j,&n_type))));
						energ+=(*(r6_env+i*(n_type)+k)*(*(Potn+mat2vec(&k,&j,&n_type)+nbij)));
					}
				}
				*(s_energy+i*(n_type)+j)=(double)0.5*energ;
			}
		}
	}
	if (op->eij || op->strucmp)
	{
   	/* calculate energy for elements of s(i,j) for non-continuous or contact potential */
		
		for (i=0; i<*Lstru; ++i)
		{
			for (j=0; j<n_type; ++j)
			{
				energ=(double)0.0;
				for (k=0; k<n_type; ++k)
				{
					if (k>j)
					{
						energ+=(double)(*(env_eij+i*(n_type)+k)*(*(Potn+mat2vec(&j,&k,&n_type))));
					}
					else
					{
						energ+=(double)(*(env_eij+i*(n_type)+k)*(*(Potn+mat2vec(&k,&j,&n_type))));
					}  
				}
				*(s_energy+i*(n_type)+j)=(double)0.5*energ;
				/* debugging */
/*				if ((i==1 && j==4) || (i==1 && j==19))
				{
					printf("for stru pos %d and aa %c, energy = %f\n",i,*dig2one(j),s_energy[i*n_type+j]);
				}
				if (i==1)
				{
					printf("energy for stru pos %d and aa %c = %f\n",i,*dig2one(j),s_energy[i*n_type+j]);
				} */
			}
		}
	} 
	
/*	energ=((float)0.5*energ); */
/* create Boltzmann like probability profile using s_energy */
	for (i=0; i<*Lstru; ++i)
	{
		sum_exp=(double)0.0;
		for (k=0; k<n_type; ++k)
		{
			m=*(s_energy+i*(n_type)+k);
			sum_exp+=(double)exp(((double)-1.0)*m);
/*			printf("normalization constant = %f\n",sum_exp); */ 
		}
		sum_prob=(double)0.0;
		for (j=0; j<n_type; ++j)
		{
			n=*(s_energy+i*(n_type)+j);
			*(profile+i*(n_type)+j)=(double)((exp(((double)-1.0)*n))/sum_exp);
			sum_prob+=*(profile+i*(n_type)+j);
/*			printf("energy = %f\n",*(s_energy+i*(n_type)+j));
			printf("e to the energy = %f\n",exp(-1*(*(s_energy+i*(n_type)+j))));
			printf("fingerprint = %f\n",*(profile+i*(n_type)+j));
			printf("\n"); */
		}
/*		if (sum_prob!=(float)1.0)
		{
			printf("the sum of the probabilities for stru pos %d does not equal 1, but = %f\n",i,sum_prob);
		} */
		/* debugging */
/*		if (i==1)
		{
			printf("sum_exp for stru pos %d = %f\n",i,sum_exp);
		} */
	}

	if (op->strucmp)
	{
		/* switch back potential pointers */
		temp_pot=x2x_Potn;
		x2x_Potn=Potn;
		Potn=temp_pot;
	}
	
	/* calculate norm of each site aa weight profile vector */
	for (j=0; j<*Lstru; ++j)
	{
		*(norm_pj+j)=(double)0.0;
		*(entropy+j)=(double)0.0;
		for (k=0; k<n_type; ++k)
		{
			*(norm_pj+j)+=(double)pow(*(profile+j*(n_type)+k),2);
			*(entropy+j)+=(double)(*(profile+j*(n_type)+k)*log(*(profile+j*(n_type)+k)));
		}
		*(norm_pj+j)=(double)sqrt(*(norm_pj+j));
		*(entropy+j)=(double)(-1)*(*(entropy+j)/(log((double)n_type)));
	}

	/* calculate the deviation from native for each amino acid of each 
		site aa weight profile vector */
	for (i=0; i<*Lstru; ++i)
	{
		for (j=0; j<n_type; ++j)
		{
			*(dev_pij+i*(n_type)+j)=*(profile+i*(n_type)+j)-*(profile+i*(n_type)+*(seq+i));
		}
	}

	/* if quordata==0 then fngrprnts are written to query fngrprnt files, if quordata==1 then
		fngrprnts are written to database fngrprnt files */
	if (quordata==0)
	{
		fprintf(F->q_fngrps,"%s\t%d\t%f\n",R_name,*Lstru,op->resc_pot);
		fprintf(F->q_fngrps2,"%s\t%d\n",R_name,*Lstru);
		fprintf(F->q_fngrps3,"%s\t%d\n",R_name,*Lstru);
		for (i=0; i<*Lstru; ++i)
		{
			fprintf(F->q_fngrps,"%4d  %2d(%c)  %4.3f  %4.3f    ",i,*(seq+i),*dig2one(*(seq+i)),profile[i*n_type+*(seq+i)],entropy[i]);
			fprintf(F->q_fngrps2,"%4d  %2d(%c)  % 05.4f    ",i,*(seq+i),*dig2one(*(seq+i)),dev_pij[i*n_type+*(seq+i)]);
			for (j=0; j<n_type; ++j)
			{
				fprintf(F->q_fngrps,"%4.3f  ",*(profile+i*n_type+j));
				fprintf(F->q_fngrps2,"% 05.4f ",*(dev_pij+i*n_type+j));
			}
			fprintf(F->q_fngrps3,"%4.3f\n",*(norm_pj+i));
			fprintf(F->q_fngrps,"\n");
			fprintf(F->q_fngrps2,"\n");
		} 
/*		fprintf(F->q_fngrps,"\n");
		fprintf(F->q_fngrps2,"\n");
		fprintf(F->q_fngrps3,"\n"); */

		free(fps->s_energy);
		free(fps->profile);
		free(dev_pij);
		free(norm_pj);
		free(entropy);
	}

	if (quordata==1)
	{
		fprintf(F->fngrps,"%s\t%d\t%f\n",R_name,*Lstru,op->resc_pot);
		fprintf(F->fngrps2,"%s\t%d\n",R_name,*Lstru);
		fprintf(F->fngrps3,"%s\t%d\n",R_name,*Lstru);
		fprintf(F->E_fngrps,"%s\t%d\n",R_name,*Lstru);
		fprintf(F->Sfile,"%s\t%d\n",R_name,*Lstru);

/*	fwrite(profile,sizeof profile,1,F->fngrps); */
/*  write profile or fingerprint information to fingerprint files */
		for (i=0; i<*Lstru; ++i)
		{
			fprintf(F->fngrps,"%4d  %2d(%c)  %4.3f  %4.3f    ",i,*(seq+i),*dig2one(*(seq+i)),profile[i*n_type+*(seq+i)],entropy[i]);
			fprintf(F->fngrps2,"%4d  %2d(%c)  % 05.4f    ",i,*(seq+i),*dig2one(*(seq+i)),dev_pij[i*n_type+*(seq+i)]);
			fprintf(F->E_fngrps,"%4d  %2d(%c)  %-+8.2f  ",i,*(seq+i),*dig2one(*(seq+i)),s_energy[i*n_type+*(seq+i)]);
			for (j=0; j<n_type; ++j)
			{
				fprintf(F->fngrps,"%4.3f  ",*(profile+i*n_type+j));
				fprintf(F->fngrps2,"% 05.4f ",*(dev_pij+i*n_type+j)); 
				fprintf(F->E_fngrps,"%-+8.2f",*(s_energy+i*n_type+j));
			}
			fprintf(F->fngrps3,"%4.3f\n",*(norm_pj+i));
			fprintf(F->fngrps,"\n");
			fprintf(F->fngrps2,"\n"); 
			fprintf(F->E_fngrps,"\n");
			fprintf(F->Sfile,"%4.3f\n",*(entropy+i));
		} 
/*		fprintf(F->fngrps,"\n");
		fprintf(F->fngrps2,"\n");
		fprintf(F->fngrps3,"\n");
		fprintf(F->E_fngrps,"\n"); */

		free(fps->s_energy);
		free(fps->profile);
		free(dev_pij);
		free(norm_pj);
		free(entropy);
	}				
}

void init_fngrps(a_file *F,k_opt *op,fgrprnt *fps)
{
	int i;

	for (i=0; i<1; ++i)
	{
		printf("i am crazy\n");
	}
}

void set_fpp(a_file *F,k_opt *op,fgrprnt *fps)
{
	int i;

	fps=&fps_data;
 
	 if (op->strucmp)
	 {	 
	/*	fps_Potn=fps->fps_Potn; */
	/* set fingerprint potential to default eij potential or evdw potential*/
		 if (fps->fps_potn==0)
		 {
			 fps->fpsn_type=NTYPE;
			 fps->x2x_Potn=malloc(((fps->fpsn_type*(fps->fpsn_type+1))/2)*sizeof(float));
			 if (!fps->x2x_Potn)
			 {
				 printf("unable to allocate memory for x2x_Potn\n");
				 exit(1);
			 }
			 tobi_elber_eij(fps->x2x_Potn);
		 }
		 if (fps->fps_potn==1)
		 {
			 fps->fpsn_type=NTYPE/2;
			 fps->x2x_Potn=malloc((fps->fpsn_type*(fps->fpsn_type+1))*sizeof(float));
			 if (!fps->x2x_Potn)
			 {
				 printf("unable to allocate memory for x2x_Potn\n");
				 exit(1);
			 }
			 meller_elber_vdw(fps->x2x_Potn);
		 }
	/* switch pointers between stru2stru or threading potential and fingerprnt potential */
	/*	temp_pot=Potn;
		Potn=fps_Potn;
		fps_Potn=temp_pot; */
		if (op->resc_pot!=(float)1.0)
		{ 
			for (i=0; i<210; ++i)
			{
				fps->x2x_Potn[i]=((fps->x2x_Potn[i])*(op->resc_pot));
/*				++fps->x2x_Potn; */
			}
		} 
/*		*naa_types=fps->fpsn_type; */
	 }
}

void prepare_fps_files(a_file *F,k_opt *op)
{
	strcpy(F->f_fps2,"FINGERPS2");
	strcpy(F->f_fps3,"FINGERPS3");
	strcpy(F->E_fps,"E_FINGERPS");
	strcpy(F->f_qfps,"Q_FINGERPS");
	strcpy(F->f_qfps2,"Q_FINGERPS2");
	strcpy(F->f_qfps3,"Q_FINGERPS3");
	strcpy(F->f_norm,"FINGERPS_NORM");
	strcpy(F->f_dv,"FINGERPS_DV");
	strcpy(F->f_entropy,"FPS_ENTROPY");

	/* if op->fngrps, then open fingerprint files */
/*	if (F->fngrps==NULL)
	{
		F->fngrps=fopen(F->f_fps,"w");
		if (F->fngrps==NULL)
		{
			printf("cannot open file fngrps for writing\n");
			exit(1);
		}
	} */
	if (!(op->r_cont))
	{
		fclose(F->fngrps);
		F->fngrps=fopen(F->f_fps,"w");
		if (F->fngrps==NULL)
		{
			printf("cannot open file FINGERPS for writing in init_fngrps\n");
			exit(1);
		}
	}
	F->fngrps2=fopen(F->f_fps2,"w");
	if (F->fngrps2==NULL)
	{
		printf("cannot open file fngrps2 for writing\n");
		exit(1);
	}
	F->fngrps3=fopen(F->f_fps3,"w");
	if (F->fngrps3==NULL)
	{
		printf("cannot open file fngrps3 for writing\n");
		exit(1);
	}
	F->E_fngrps=fopen(F->E_fps,"w");
	if (F->E_fngrps==NULL)
	{
		printf("cannot open file E_fngrps for writing\n");
		exit(1);
	}
	F->q_fngrps=fopen(F->f_qfps,"w");
	if (F->fngrps2==NULL)
	{
		printf("cannot open file fngrps2 for writing\n");
		exit(1);
	}
	F->q_fngrps2=fopen(F->f_qfps2,"w");
	if (F->fngrps3==NULL)
	{
		printf("cannot open file fngrps3 for writing\n");
		exit(1);
	}
	F->q_fngrps3=fopen(F->f_qfps3,"w");
	if (F->E_fngrps==NULL)
	{
		printf("cannot open file E_fngrps for writing\n");
		exit(1);
	}
	F->fps_norm=fopen(F->f_norm,"w");
	if (F->fps_norm==NULL)
	{
		printf("cannot open file fps_norm for writing\n");
		exit(1);
	}
	F->fps_dv=fopen(F->f_dv,"w");
	if (F->fps_dv==NULL)
	{
		printf("cannot open file fps_dv for writing\n");
		exit(1);
	}
	F->Sfile=fopen(F->f_entropy,"w");
	if (F->Sfile==NULL)
	{
		printf("cannot open file Sfile for writing\n");
		exit(1);
	}
}

void free_up_fea(k_opt *op,f_mod *fea)
{
	if (op->evdw)
	{
		free(fea->r12_env);
		free(fea->r6_env);
	}
	if (op->eij)
	{
		free(fea->env_eij);
		free(fea->flag);
	}
	free(fea->tn_cont);
}

void close_fps_files(a_file *F,k_opt *op,fgrprnt *fps,int quordata)
{
	if (quordata==0)
	{
		fclose(F->q_fngrps); 
		fclose(F->q_fngrps2);
		fclose(F->q_fngrps3);
	}
	if (quordata==1)
	{
		fclose(F->fngrps); 
		fclose(F->fngrps2);
		fclose(F->fngrps3);
		fclose(F->E_fngrps);
		fclose(F->Sfile);
	}
}

int read_fngrps(FILE *fpsf,char *fname,k_opt *op,fgrprnt *fps,float *fngrprnt,float *rsc,int *Lstru,char *X_name,int quordata)
{
	int eof,i,j,naa,pos,naa_type;
	char aa;
	float natfps,S,*array;

	naa=0;
	pos=0;
	aa='0';

	naa_type=fps->fpsn_type;

/*	if (fpsf==NULL)
	{
		fpsf=fopen(fname,"r");
		if (fpsf==NULL)
		{
			printf("cannot open file %s for reading for stru2stru align\n",fname);
			return(0);
		}
	} */

	eof=fscanf(fpsf,"%s\t%d\t%f\n",X_name,Lstru,rsc);
	if(eof==EOF)
	{
		rewind(fpsf);
		fscanf(fpsf,"%s\t%d\t%f\n",X_name,Lstru);
	}

	if (quordata==0)
	{
		fps->C_profile=malloc((*Lstru*naa_type)*sizeof(float));
		if (!fps->C_profile)
		{
			printf("unable to allocate memory for individual protein fingerprints\n");
			return(0);
		}
		fngrprnt=fps->C_profile;

	}
	if (quordata==1)
	{
		fps->R_profile=malloc((*Lstru*naa_type)*sizeof(float));
		if (!fps->R_profile)
		{
			printf("unable to allocate memory for individual protein fingerprints\n");
			return(0);
		}
		fngrprnt=fps->R_profile;
	}
	/* if quordata==2, then read into a dummy storage array (not going to be used) because
		the fxn is being called just to read the protein names (such as in check_CM) */
	if (quordata==2)
	{
		array=malloc((*Lstru*naa_type)*sizeof(float));
		if (!array)
		{
			printf("unable to allocate dummy array fngrprnt in read_fngrps\n");
			return(0);
		}
		fngrprnt=array;
	}

	for (i=0; i<*Lstru; ++i)
	{
		for (j=0; j<naa_type; ++j)
		{
			if (j==0)
			{
				fscanf(fpsf,"%d%d(%c)%f%f%f",&pos,&naa,&aa,&natfps,&S,fngrprnt);
				++fngrprnt;
			}
			if (j<op->n_type-1 && j!=0)
			{
				fscanf(fpsf,"%f",fngrprnt);
				++fngrprnt;
			}
			if (j==op->n_type-1)
			{
				fscanf(fpsf,"%f\n",fngrprnt);
				++fngrprnt;
			}
		}
	} 

/*	for (i=0; i<(*Lstru*op->n_type); ++i)
	{
		fscanf(fpsf,"%f",fngrprnt);
		++fngrprnt;
	} */

	/* debugging */
/*	for (i=0; i<*Lstru; ++i)
	{
		for (j=0; j<op->n_type; ++j)
		{
			if (j<op->n_type-1)
			{
/*				printf("%4.3f ",fngrprnt);
				++fngrprnt; */
/*				printf("%4.3f",fngrprnt[i*op->n_type+j]); 
			}
			else
			{
/*				printf("%4.3f\n",fngrprnt);
				++fngrprnt; */
/*				printf("%4.3f\n",fngrprnt[i*op->n_type+j]); 
			}
		}
	} */
	
	if (quordata==2)
	{
		free(array);
	}
	return(1);
}