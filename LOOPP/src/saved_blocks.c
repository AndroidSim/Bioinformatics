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

/*void get_env_eij(k_opt *op,f_mod *fea,int *Imod,int *contT,int *IJ,int *nat_Seq,
				int *Lnatseq,float *contTvdw,float *r12ij,float *r6ij);
float get_energy_eij(k_opt *op,f_mod *fea,int *Imod,int *contT,int *start,
							int *j_beta,int *qseq_i_alpha,int *nat_Seq,int *Lnatseq,
                     float *Potn,float *contTvdw,float *r12ij,float *r6ij,
                     float *gpen); */

/*	for (i=0; i<*Lnatseq; ++i)
	{
		for (j=0; j<op->n_type; ++j)
		{
			*(tn_cont+i)+=*(env_eij+i*(op->n_type)+j);
			(*nop)++;
		}
	}*/
/*	for (i=0; i<*Lnatseq; ++i)
	{
		for (j=0; j<op->n_type; ++j)
		{
			if (j<op->n_type-1)
			{
				printf("%d",env_eij[i*op->n_type+j]);
			}
			else
			{
				printf("%d\n",env_eij[i*op->n_type+j]);
			}
		}
	} */
 /*     for (i=0;i<op->n_par;++i)
      {
      	*(contTvdw+i)=0.0;
      }*/
/*	else
   {
   	for (i=0;i<op->n_par;++i)
      {
      	*(contT+i)=0;
		(*nop)++;
      }
   }*/
/*	for (j=0; j<op->n_type; ++j)
	{
		n_cont+=env_eij[*j_beta*(op->n_type)+j];
		(*nop)++;
	} */
/*      for (i=0; i<alpha; ++i)
      {
      	par_ind+=i;
      }
      for (i=0; i<alpha; ++i)
      {
         contTvdw[par_ind+i]=fea->r12_env[*j_beta*(op->n_type)+i];
         contTvdw[par_ind+i]+=fea->r6_env[*j_beta*(op->n_type)+i];
      }
      for (k=alpha; k<=op->n_type; ++k)
      {
      	par_ind=par_ind+k;
        contTvdw[par_ind+(alpha-1)]=fea->r12_env[*j_beta*(op->n_type)+(k+1)];
        contTvdw[par_ind+(alpha-1)]+=fea->r6_env[*j_beta*(op->n_type)+(k+1)];
      } */


/*	  for (i=0; i<=alpha; ++i)
      {
      	par_ind+=i;
      }
      for (j=0; j<=alpha; ++j)
      {
      	contT[par_ind+j]=fea->env_eij[*j_beta*(op->n_type)+j];
      }
      for (k=alpha+1; k<op->n_type; ++k)
      {
      	par_ind+=k;
        contT[par_ind+alpha]=fea->env_eij[*j_beta*(op->n_type)+k];
      }  */

/*	  for (i=0; i<alpha; ++i)
	  {
		  *(contT+mat2vec(&i,&alpha,&op->n_type))=0;
		  *(contT+mat2vec(&i,&alpha,&op->n_type))+=*(env_eij+*j_beta*(op->n_type)+i);
		  energ+=(float)*(contT+mat2vec(&i,&alpha,&op->n_type))*(*(Potn+mat2vec(&i,&alpha,&op->n_type)));
		  (*nop)++;
	  }
	  for (j=alpha; j<op->n_type; ++j)
	  {
		  *(contT+mat2vec(&alpha,&j,&op->n_type))=0;
		  *(contT+mat2vec(&alpha,&j,&op->n_type))+=*(env_eij+*j_beta*(op->n_type)+j);
		  energ+=(float)*(contT+mat2vec(&alpha,&j,&op->n_type))*(*(Potn+mat2vec(&alpha,&j,&op->n_type)));
		  (*nop)++;
	  }*/
/*    if(op->evdw)
	{
		for(i=0;i<op->n_par;++i)
		{
			if(*(contTvdw+i)!=0.0)
			{
				energ+=*(contTvdw+i)*(*Potn);
			}
			++Potn;
		}
	}
	else
	{
		for(i=0;i<op->n_par;++i)
		{
			if(*(contT+i)!=0)
			{
				energ+=(float)*(contT+i)*(*Potn);
				(*nope)++;
			}
			++Potn;
		}
	}*/
/*	return((float)0.5*energy(op,contT,Potn,contTvdw));*/
/*	 if (op->fngrps)
	 {
		 if (F->fngrps==NULL)
		 {
			 F->fngrps=fopen(F->f_fps,"r");
			 if (F->fngrps==NULL)
			 {
				 printf("cannot open file fngrps for reading\n");
				 exit(1);
			 }
		 }

		 eof=fscanf(F->fngrps,"%s\t%d\n",name,&len_stru);
		 if(eof==EOF)
		 {
			 rewind(F->fngrps);
			 fscanf(F->fngrps,"%s\t%d\n",name,&len_stru);
		 }		 
		 while ((strcmp(name,C_name)!=0) && len_stru!=l_qstru)
		 {
			 eof=fscanf(F->fngrps,"%s\t%d\n",name,&len_stru);
			 if(eof==EOF)
			 {
				 rewind(F->fngrps);
				 fscanf(F->fngrps,"%s\t%d\n",name,&len_stru);
			 }
		 }

		 fps->C_profile=malloc((l_qstru*op->n_type)*sizeof(float));
		 if (!fps->C_profile)
		 {
			 printf("unable to allocate memory for C_profile\n");
			 exit(1);
		 }
		 C_profile=fps->C_profile;

		 for (i=0; i<(l_qstru*op->n_type); ++i)
		 {
			 fscanf(F->fngrps,"%f",C_profile);
			 ++C_profile;
		 }
	 } */
/*		 if (op->fngrps)
		 {
			 if (F->fngrps==NULL)
			 {
				 F->fngrps=fopen(F->f_fps,"r");
				 if (F->fngrps==NULL)
				 {
					 printf("cannot open file fngrps for reading\n");
					 exit(1);
				 }
			 }

			 eof=fscanf(F->fngrps,"%s\t%d\n",name,&len_stru);
			 if(eof==EOF)
			 {
				 rewind(F->fngrps);
				 fscanf(F->fngrps,"%s\t%d\n",name,&len_stru);
			 }		 
			 while ((strcmp(name,R_name)!=0) && len_stru!=n_col)
			 {
				 eof=fscanf(F->fngrps,"%s\t%d\n",name,&len_stru);
				 if(eof==EOF)
				 {
					 rewind(F->fngrps);
					 fscanf(F->fngrps,"%s\t%d\n",name,&len_stru);
				 }
			 }

			 fps->R_profile=malloc((l_qstru*op->n_type)*sizeof(float));
			 if (!fps->R_profile)
			 {
				 printf("unable to allocate memory for R_profile\n");
				 exit(1);
			 }
			 R_profile=fps->R_profile;

			 for (i=0; i<(l_qstru*op->n_type); ++i)
			 {
				 fscanf(F->fngrps,"%f",R_profile);
				 ++R_profile;
			 }
		 } */

/*				if (j<n_type-1)
				{			
					fprintf(F->q_fngrps,"%4.3f  ",profile[i*n_type+j]);
					fprintf(F->q_fngrps2,"% 05.4f ",dev_pij[i*n_type+j]); 
				}
				else
				{
					fprintf(F->q_fngrps,"%4.3f\n",profile[i*n_type+j]);
					fprintf(F->q_fngrps2,"% 05.4f\n",dev_pij[i*n_type+j]); 
				} */
/*			if (dev_pij[i*n_type+*(seq+i)]>=(float)0.0)
			{
				fprintf(F->fngrps2,"%4d  %2d(%c)   %5.4f    ",i,*(seq+i),*dig2one(*(seq+i)),dev_pij[i*n_type+*(seq+i)]);
			}
			else
			{
				fprintf(F->fngrps2,"%4d  %2d(%c)  %5.4f    ",i,*(seq+i),*dig2one(*(seq+i)),dev_pij[i*n_type+*(seq+i)]);
			}
			if (s_energy[i*n_type+*(seq+i)]>=(float)0.0)
			{
				fprintf(F->E_fngrps,"%4d  %2d(%c)   %6.3f    ",i,*(seq+i),*dig2one(*(seq+i)),s_energy[i*n_type+*(seq+i)]);
			}
			else
			{
				fprintf(F->E_fngrps,"%4d  %2d(%c)  %6.3f    ",i,*(seq+i),*dig2one(*(seq+i)),s_energy[i*n_type+*(seq+i)]);
			} */
/*					if (dev_pij[i*n_type+j]>=(float)0.0)
					{
						fprintf(F->fngrps2," %5.4f ",dev_pij[i*n_type+j]);
					}
					else
					{
						fprintf(F->fngrps2,"%5.4f ",dev_pij[i*n_type+j]);
					}
					if (s_energy[i*n_type+j]>=(float)0.0)
					{
						fprintf(F->E_fngrps," %6.3f ",s_energy[i*n_type+j]);
					}
					else
					{
						fprintf(F->E_fngrps,"%6.3f ",s_energy[i*n_type+j]);
					} */
/*					if (dev_pij[i*n_type+j]>=(float)0.0)
					{
						fprintf(F->fngrps2," %5.4f\n",dev_pij[i*n_type+j]);
					}
					else
					{
						fprintf(F->fngrps2,"%5.4f\n",dev_pij[i*n_type+j]);
					}
					if (s_energy[i*n_type+j]>=(float)0.0)
					{
						fprintf(F->E_fngrps," %6.3f\n",s_energy[i*n_type+j]);
					}
					else
					{
						fprintf(F->E_fngrps,"%6.3f\n",s_energy[i*n_type+j]);
					} */
/*				if (j<n_type-1)
				{			
					fprintf(F->fngrps,"%4.3f  ",profile[i*n_type+j]);
					fprintf(F->fngrps2,"% 05.4f ",dev_pij[i*n_type+j]); 
					fprintf(F->E_fngrps,"%-+8.2f",s_energy[i*n_type+j]);
				}
				else
				{
					fprintf(F->fngrps,"%4.3f\n",profile[i*n_type+j]);
					fprintf(F->fngrps2,"% 05.4f\n",dev_pij[i*n_type+j]); 
					fprintf(F->E_fngrps," %-+8.2f\n",s_energy[i*n_type+j]);
				} */

/* example of trying to use pointer increment in writing to files
	--does not work-- */
/*		for (i=0; i<*Lstru; ++i)
		{
			for (j=0; j<n_type; ++j)
			{
				if (j<n_type-1)
				{
					fprintf(F->fngrps,"%4.3f ",*profile);
	/*				fprintf(F->fngrps2,"%5.4f ",dev_pij[i*n_type+j]); */
	/*				if (*dev_pij>=0)
					{
						fprintf(F->fngrps2," %5.4f ",*dev_pij);
					}
					else
					{
						fprintf(F->fngrps2,"%5.4f ",*dev_pij);
					} 
				}
				else
				{
					fprintf(F->fngrps,"%4.3f\n",*profile);
	/*				fprintf(F->fngrps2,"%5.4f\n",dev_pij[i*n_type+j]); */
	/*				if (*dev_pij>=0)
					{
						fprintf(F->fngrps2," %5.4f\n",*dev_pij);
					}
					else
					{
						fprintf(F->fngrps2,"%5.4f\n",*dev_pij);
					} 
				}
				++profile;
				++dev_pij; 
			}
			fprintf(F->fngrps3,"%4.3f\n",norm_pj);
			++norm_pj;
		} 
		fprintf(F->fngrps,"\n");
		fprintf(F->fngrps2,"\n");
		fprintf(F->fngrps3,"\n"); */

/* get energy of native structure alignment using norm(dis_fingerprints) scoring scheme */
/*	 *(R_ene+i_str)=0.0;
	 for (i=0; i<n_col; ++i)
	 {
		 norm=(double)0.0;
		 for (k=0; k<op->n_type; ++k)
		 {
			 norm+=pow(*(R_profile+i*(op->n_type)+k),2);
		 }
		 norm=sqrt(norm);
		 *(R_ene+i_str)+=norm;
	 }
	 /* or (this is the right one) */
/*	 *(R_ene+i_str)=0.0;
	 for (i=0; i<n_col; ++i)
	 {
		 *(R_ene+i_str)+=(-1*0.5);
	 } */

/* normalize score Tij from -1/2 -> 1/2 using length of alignment 
	as normalization constant, right now using length of longest structure as constant */
/*			 if (l_qstru>lstru)
			 {
				 l_align=l_qstru;
				 Tij=(Tij/l_align);
			 }
			 else if (lstru>l_qstru)
			 {
				 l_align=lstru;
				 Tij=(Tij/l_align);
			 }
			 else
			 {
				 Tij=(Tij/lstru);
			 } */
/* normalize score Tij from -1/2 -> 1/2 using length of alignment 
	as normalization constant, right now using length of longest structure as constant */
/*		 if (*l_qstru>lstru)
		 {
			 l_align=*l_qstru;
			 Tij=(Tij/l_align);
		 }
		 else if (lstru>*l_qstru)
		 {
			 l_align=lstru;
			 Tij=(Tij/l_align);
		 }
		 else
		 {
			 Tij=(Tij/lstru);
		 } */

