/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000    v. 2.000                                   */
/*                                                                */
/*        PREPARE QUERY FILES                                     */
 
#include "loopp.h"
   
/* prepare Query file to be used as a source of query sequences */
int build_query(a_file *F, k_opt *op, k_info *Info, k_prot *prot,k_Model *model)
{
 int i_prot,lseq,n_domain,*Seq,lstru,membrane,not_empty;
 char file_name[MAXS],stmp[MAXS];
 k_res *Res;
 
 char Q_name[MAXS],A_name[MAXS];
 int *IJcont,*Imod,noper,lastru;
 float *rij,*r6ij,*r12ij,*Potn,s_pct_env;
 
 Seq=prot->seq;
 Res=prot->Res;

 rij=prot->query_IJrij;
 IJcont=prot->IJcont;
 Potn=model->Potn;
 Imod=model->Imod;
 r12ij=prot->IJr12;
 r6ij=prot->IJr6;
 rij=prot->native_IJrij;
 noper=0;
 s_pct_env=0.0;
 
 /* assumption: f_exam, f_coord read for the first time here */

 if(op->qfmt==1)  /* loopp format */
   {
     i_prot=0;
     /* begin loop over sequences for structure T_name */
     while(read_seq_generic(F->exam,&lseq,op->qname,Seq,&i_prot,op))
       {  
	 /* to use -b start end */
	 if(!op->pick_all_seq && i_prot>op->end) break;
	 if(op->debug)
	   {
	     sprintf(Info->Sout,"\n seq %s of %d res",op->qname,lseq);
	     to_stnd_out(Info->Sout);
	   }  

	 /* i_prot incremented in read_seq so it goes from 1 to max_query */
	 if(i_prot>op->beg-1 && lseq<=op->max_length) 
	   {
	     /* use standard f_coor to get cordinates */
	     if(op->strucmp) 
	       {
		 xyz2query_CM(F,op,Info,prot,&i_prot,&not_empty);
		 if(not_empty) seq2query(F,op,Info,Seq,&lseq); 
	       }
	     else seq2query(F,op,Info,Seq,&lseq);
	   }
       }   
   }
 /* end of loopp format */

 if(op->qfmt==2)  /* FASTA format: >name seq(1let) */
   { 
     prepare_query_fasta(F,op,Info);     
   }
 /* end of FASTA format */


 if(op->qfmt==3)  /* SWISS-PROT format */
   { 
     if(!prepare_query_swiss(F,op,Info,Seq))
       {
	 sprintf(Info->Sout,
		 " WARNING - there is a problem with the file %s \n",
		 F->f_exam);
	 to_stnd_out(Info->Sout);
       }
   }
 /* end of SWISS-PROT format */

  
 if(op->qfmt==4)  /* PDB format */
   { 
     /* if list then f_list is supposed to be a list of (PDB) file names */
     if(op->list)
       {
	 fclose(F->exam);
	 strcpy(stmp,F->f_exam);
	 rewind(F->list);
	 /* loop over all entries in list */
	 while(read_entry_list(F,op,Info,file_name))
	   {
	     F->exam=fopen(file_name,"r");
	     if(F->exam==NULL) open_err(file_name);
	     strcpy(F->f_exam,file_name);
	     
	     n_domain=prepare_query_pdb(F,op,Info,Seq);
	     if(n_domain==0) query_file_warning(F,op,Info);
	     if(op->strucmp) prepare_xyz_pdb(F,op,Info,Res,Seq,&n_domain);

	     fclose(F->exam);
	   }
	 if(op->strucmp) prepare_query_CM(F,op,Info,prot);
	 F->exam=fopen(stmp,"r"); /* restore original */
       }
     else /* just a single PDB file */
       {

	 n_domain=prepare_query_pdb(F,op,Info,Seq);
	 if(n_domain==0) query_file_warning(F,op,Info);
	 /* call it only in case of stru-stru alignments */
	 if(op->strucmp)
	   {
	     prepare_xyz_pdb(F,op,Info,Res,Seq,&n_domain);
	     prepare_query_CM(F,op,Info,prot);
	   }

       }
  
   }
 /* end of PDB format */
 

 if(op->qfmt==5)  /* 1let (plain) format */
   { 
     prepare_query_1let(F,op,Info);    
   }
 /* end of 1let (plain) format */


 if(op->qfmt==6)  /* CHARMM (MOIL) format */
   { 
     if(op->list)
       {
	 fclose(F->exam);
	 strcpy(stmp,F->f_exam);
	 rewind(F->list);
	 /* loop over all entries in list */
	 while(read_entry_list(F,op,Info,file_name))
	   {
	     F->exam=fopen(file_name,"r");
	     if(F->exam==NULL) open_err(file_name);
	     strcpy(F->f_exam,file_name);
	     
	     if(!prepare_query_crd(F,op,Info,prot,Seq))
	       query_file_warning(F,op,Info);

	     fclose(F->exam);
	   }
	 if(op->strucmp) prepare_query_CM(F,op,Info,prot);
	 F->exam=fopen(stmp,"r"); /* restore original */
       }
     else /* single CRD file */
       {
	 if(prepare_query_crd(F,op,Info,prot,Seq)) 
	   {
	     if(op->strucmp) prepare_query_CM(F,op,Info,prot);
	   }
	 
       }
     
   }
 /* end of CRD format */


 /* now close and reopen f_query */
 fclose(F->query); 
 F->query=fopen(F->f_query,"r"); 
 if(F->query==NULL) open_err(F->f_query);
 if(op->strucmp)
   {
     fclose(F->qrij);
     F->qrij=fopen(F->f_qrij,"r"); 
     if(F->qrij==NULL) open_err(F->f_qrij);
     fclose(F->qcont);
     F->qcont=fopen(F->f_qcont,"r"); 
     if(F->qcont==NULL) open_err(F->f_qcont);
     fclose(F->qxyz);
     F->qxyz=fopen(F->f_qxyz,"r"); 
     if(F->qxyz==NULL) open_err(F->f_qxyz);
   }
 
  /* create query fingerprints from F->qcont and store in query fingerprint
	files */
 if (op->fngrps && op->strucmp)
 {
/*	 fclose(F->qcont);
	 F->qcont=fopen(F->f_qcont,"r");
	 if (F->qcont==NULL)
	 {
	 	 printf("cannot open file fngrps for reading\n");
		 exit(1);
	 }
	 rewind(F->qcont); */
	 while(read_seq_generic(F->query,&lseq,Q_name,Seq,&i_prot,op))
	 {
		 fps=&fps_data;
		 fea=&fea_data;
		 i_prot=0;
		 read_Query_CM(F,op,IJcont,&lastru,A_name,rij);	 
		 get_env_eij(op,fea,Imod,IJcont,Seq,&lastru,&noper,r12ij,r6ij,&s_pct_env);
		 write_fngrps(F,op,fea,fps,IJcont,Seq,&lseq,&lastru,Q_name,Potn,0);
		 free_up_fea(op,fea);
	 }
	 fclose(F->q_fngrps); 
	 fclose(F->q_fngrps2);
	 fclose(F->q_fngrps3);
	 rewind(F->qcont);
	 rewind(F->query);
	 rewind(F->qrij);
 }

 
 /* now screen query sequences */
 i_prot=0;
 membrane=0;
 while(read_seq_generic(F->query,&lseq,op->qname,Seq,&i_prot,op))
   {
     /* check if membrane protein - now only warning */
     if(check_for_membrane(op,F,Info,Seq,&lseq))  membrane=1;
     /* print warning if more than 5 percent in low complexity regions */
     if(low_complexity(op,F,Info,Seq,&lseq)>5.0)
       {
	 sprintf(Info->Sout,
		 "\n WARNING - sequence %s contains LOW COMPLEXITY regions \n",
		 op->qname);
	 to_stnd_out(Info->Sout);
	 to_file(F->best,Info->Sout);
	 sprintf(Info->Sout,
		 "           that may result in spuriously high Z-scores! \n");
	 to_stnd_out(Info->Sout);
	 to_file(F->best,Info->Sout);
       }
     
     /* include here other filters ... */
   }
 if(membrane && op->exam) membrane_warning(F,op,Info);
 

 /* rewind the files */
 rewind(F->exam);
 rewind(F->coor);
 rewind(F->seq);  
 rewind(F->query);

 return(1);

}

void query_file_warning(a_file *F, k_opt *op, k_info *Info)
{

  sprintf(Info->Sout,
	  "\n WARNING - there is a problem with the file %s \n",
	  F->f_exam);
  to_stnd_out(Info->Sout);
  
}

void membrane_warning(a_file *F, k_opt *op, k_info *Info)
{
  
  sprintf(Info->Sout,
    "\n Our threading predictions are not reliable for membrane proteins!\n");
  to_stnd_out(Info->Sout);
  to_file(F->best,Info->Sout);
  
}

/* check how many sites do not have coordinates assigned */
float check_if_empty(int *Lstru, float *X, float *Y, float *Z)
{
  int i,n_empty;
  float perc_empty;
  
  perc_empty=0.0;
  n_empty=0;
  
  for(i=0;i<*Lstru;++i)
    {
      if(*(X+i)>999.0) ++n_empty;
    }
  perc_empty=(100.0*n_empty)/(float)(*Lstru);
  return(perc_empty);
  
}

/* to prepare query f_qcont file if using loopp format */
int xyz2query_CM(a_file *F, k_opt *op, k_info *Info, k_prot *prot, int *i_prot,
		 int *not_empty)
{
  int Lstru;
  float *X, *Y, *Z, *Xa, *Ya, *Za, *Xb, *Yb, *Zb;
  char X_name[MAXS];
  
  /* prepare local shortcuts */
  X=prot->X_res;
  Y=prot->Y_res;
  Z=prot->Z_res;

  /* prepare additional vectors for coordinates */
  Xa=malloc(op->max_length*sizeof(int));
  Ya=malloc(op->max_length*sizeof(int));
  Za=malloc(op->max_length*sizeof(int));
  Xb=malloc(op->max_length*sizeof(int));
  Yb=malloc(op->max_length*sizeof(int));
  Zb=malloc(op->max_length*sizeof(int));
  if(Xa==NULL || Ya==NULL || Za==NULL || Xb==NULL || Yb==NULL || Zb==NULL) 
    alloc_err("local_copy_of_theStructure"); 
  
  if(!read_coor_all3(F->coor,op,&Lstru,X_name,X,Y,Z,Xa,Ya,Za,Xb,Yb,Zb)) 
    return(0);

  /* check the percentage of sites without coor assigned */
  if(check_if_empty(&Lstru,X,Y,Z)>TH_EMPTY) 
    {
      *not_empty=0; 
      sprintf(Info->Sout,
	      " Skipping %s - no coor assigned for more than %5.2f %% of res\n",
	      op->qname,TH_EMPTY);
      to_stnd_out(Info->Sout);     
    }
  else *not_empty=1;

  
  if((*i_prot>op->beg-1) && (strcmp(op->qname,X_name)==0) && *not_empty )
    {
      /* if list then f_list is supposed to be a list of names (LOOPP fmt) */
      if(op->list && op->qfmt==1)
	{
	  rewind(F->list);
	  /* loop over all entries in list */
	  while(read_entry_list(F,op,Info,op->pseq_name))
	    {
	      if(strncmp(op->qname,op->pseq_name,strlen(op->pseq_name))==0)
		{
		  write_Query_CM(F,&Lstru,X_name,X,Y,Z,op);
		  write_query_XYZ(F,op,Info,&Lstru,X,Y,Z,Xa,Ya,Za,Xb,Yb,Zb);
		}
	    }
      
	}
      else
	{
	  if(op->pick) /* if pick by name */
	    {
	      if(op->pick_all_seq) strcpy(op->pseq_name,op->qname);
	      if(strncmp(op->qname,op->pseq_name,strlen(op->pseq_name))==0)
		{
		  write_Query_CM(F,&Lstru,X_name,X,Y,Z,op);
		  write_query_XYZ(F,op,Info,&Lstru,X,Y,Z,Xa,Ya,Za,Xb,Yb,Zb);
		}
	    }
	  else 
	    {
	      write_Query_CM(F,&Lstru,X_name,X,Y,Z,op); 
	      write_query_XYZ(F,op,Info,&Lstru,X,Y,Z,Xa,Ya,Za,Xb,Yb,Zb); 
	    }
	}
      
    }

  free(Xa);
  free(Ya);
  free(Za);
  free(Xb);
  free(Yb);
  free(Zb);
  
  return(1);
}

/* to prepare query f_qcont file if using PDB format */
int prepare_query_CM(a_file *F, k_opt *op, k_info *Info, k_prot *prot)
{
  int Lstru;
  float *X, *Y, *Z;
  char X_name[MAXS];
  
  /* first close and reopen f_qxyz */
  fclose(F->qxyz);
  F->qxyz=fopen(F->f_qxyz,"r");
  /* prepare local shortcuts */
  X=prot->X_res;
  Y=prot->Y_res;
  Z=prot->Z_res;
  /* op->distance is used inside write_Query_CM */
  while(read_coor(F->qxyz,op,&Lstru,X_name,X,Y,Z))
    write_Query_CM(F,&Lstru,X_name,X,Y,Z,op);

  return(1);
}

/* read another name (of file or protein) from f_exam if op->list */
int read_entry_list(a_file *F, k_opt *op, k_info *Info, char *name)
{
  int eof;
  
  eof=fscanf(F->list,"%s",name);
  if(eof==EOF) return(0);   

  return(1);
  
}

/* write sequence into a file f_query (if picked) */
void seq2query(a_file *F, k_opt *op, k_info *Info,int *Seq,int *lseq)
{

  /* if list then f_list is supposed to be a list of names (LOOPP fmt) */
  /* switch this off when SWISS ... */
  /* pick conditions here must be consistent with xyz2query - add function!!! */
  if(op->list && op->qfmt==1)
    {
      rewind(F->list);
      /* loop over all entries in list */
      while(read_entry_list(F,op,Info,op->pseq_name))
	{
	  if(strncmp(op->qname,op->pseq_name,strlen(op->pseq_name))==0)
	    write_down_seq(F,op,Info,Seq,lseq); 	  
	}
      
    }
  else
    {
      /* here use standard pick to include/exclude sequences */
      if(op->pick) /* if pick by name */
	{
	  if(op->pick_all_seq) strcpy(op->pseq_name,op->qname);
	  if(strncmp(op->qname,op->pseq_name,strlen(op->pseq_name))==0)
	    write_down_seq(F,op,Info,Seq,lseq); 
	}
      else write_down_seq(F,op,Info,Seq,lseq);   
    }
  
}

/* write sequence into a file f_query using standard loopp format */
void write_down_seq(a_file *F,k_opt *op,k_info *Info,int *Seq,int *lseq)
{
 int i;

 /* assuming that sequence is defined in terms of generic alphabet */
 
 /* write down the query sequence now */
 sprintf(Info->Sout,"%s\t%d \n",op->qname,*lseq);
 to_file(F->query,Info->Sout);
 
 for(i=0;i<*lseq;++i)
   { 
     sprintf(Info->Sout,"%s ",dig2amino(*(Seq+i)));
     if(((i+1)%20)== 0) strcat(Info->Sout,"\n");
     to_file(F->query,Info->Sout);
   }
 sprintf(Info->Sout,"\n");
 to_file(F->query,Info->Sout);
}

/* write coordinates into a file f_qxyz using standard loopp format */
void write_query_XYZ(a_file *F, k_opt *op, k_info *Info, int *lstru,
		     float *X, float *Y, float *Z, float *Xa, float *Ya, 
		     float *Za, float *Xb, float *Yb, float *Zb)
{
 int i;
 char fmt[MAXS];
 
 /* assuming that only the right 3-column block will be used */
 sprintf(Info->Sout,"%s\t%d \n",op->qname,*lstru);
 to_file(F->qxyz,Info->Sout);

 if(op->n_col==1)
   strcpy(fmt," %10.3f %10.3f %10.3f\n");
 if(op->n_col==2)
   strcpy(fmt," %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n");
 if(op->n_col==3) strcpy(fmt,
      " %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n");

 /* so now copy X, Y, Z into each 3-column block */
 for(i=0;i<*lstru;++i)
   { 
     if(op->n_col==1)
       sprintf(Info->Sout,fmt,*(Xa+i),*(Ya+i),*(Za+i));
     if(op->n_col==2)
       sprintf(Info->Sout,fmt,*(Xa+i),*(Ya+i),*(Za+i),*(X+i),*(Y+i),*(Z+i));
     if(op->n_col==3)
       sprintf(Info->Sout,fmt,*(Xa+i),*(Ya+i),*(Za+i),*(X+i),*(Y+i),*(Z+i),
	       *(Xb+i),*(Yb+i),*(Zb+i)); 
     to_file(F->qxyz,Info->Sout);
   }
 sprintf(Info->Sout,"\n");
 to_file(F->qxyz,Info->Sout);
}

/* to prepare standard f_query file if using swiss format */
int prepare_query_swiss(a_file *F, k_opt *op, k_info *Info, int *Seq)
{
 int n_seq,i,i_length,i_found_let,j,*Lngs,eof,error,w;
 char ident[MAXS],stmp[MAXS],buffer[MAXS],*ctmp; 
 

 /* assuming standard SWISS-PRO format */
 ctmp=buffer;
 n_seq=0;
 eof=fscanf(F->exam,"%s",ident);
 if(eof==EOF) return(0); 
 
 while(eof!=EOF)
   {
     /* read the whole line if not name or length */
     if(strcmp(ident,"ID")!=0 && strcmp(ident,"SQ")!=0 && 
	strcmp(ident,"ORIGIN")!=0)
       {
	 w='\0';
	 while(w!='\n') w=fgetc(F->exam);
       }
     if(strcmp(ident,"ID")==0)
	{
	  eof=fscanf(F->exam,"%s",op->qname);
	  if(eof==EOF) return(0); 
	  w='\0';
	  while(w!='\n') w=fgetc(F->exam);
	}
     if(strcmp(ident,"SQ")==0)
	{
	  eof=fscanf(F->exam,"%s",stmp);
	  if(eof==EOF) return(0); 
	  eof=fscanf(F->exam,"%d",&i_length);
	  if(eof==EOF) return(0); 
	  w='\0';
	  while(w!='\n') w=fgetc(F->exam);
	}
     if(strcmp(ident,"ORIGIN")==0)
       {
	  eof=fscanf(F->exam,"%s",stmp);
	  if(eof==EOF) return(0); 
	  i_found_let=0;
	  error=0;
	  
	  while(strcmp(stmp,"//")!=0) 
	    {
	      w=fgetc(F->exam); 
	      *ctmp=w;
	      while(w!='\n' && w!=EOF) 
		{
		  if(w!=' ' && w!='\n' && i_found_let<op->max_length)
		    {
		      /* j=one2dig(ctmp);
		      *(Seq+i_found_let)=j;
		      if(j<0) error=1; 
		      ++i_found_let; */
		      j=one2dig(ctmp);
		      /* skip wrong symbols */
		      if(j>=0)  
			{
			  *(Seq+i_found_let)=j;
			  ++i_found_let; 
			}
		    }
		  w=fgetc(F->exam);
		  *ctmp=w;
		}	
	      eof=fscanf(F->exam,"%s",stmp); 
	      if(eof==EOF) return(0);
	    }

	 /* skip the sequence if anything wrong */ 
	 if(i_found_let>=op->max_length || error)
	   {
	     if(op->info) skip_query_warn(op,Info); 
	   }
	 else 
	   {
	     if(op->info && i_length<i_found_let)
	       {
		 sprintf(Info->Sout,
			 " WARNING - sequence %s shorter than declared \n",
			 op->qname);
		 to_stnd_out(Info->Sout);		 
	       }
	     seq2query(F,op,Info,Seq,&i_found_let); 
	   }

       }

     eof=fscanf(F->exam,"%s",ident); /* start reading again */
   }
 
 return(1);
 
}

/* to prepare standard f_query file if using PDB format */
int prepare_query_pdb(a_file *F,k_opt *op,k_info *Info,int *Seq)
{
 int n_seq,i,i_length,i_found_let,j,eof,error,length_prv,sum,sum_prv,
     header,read_next,i_domain;
 char ident[MAXS],stmp[MAXS],sym[MAXS],chain[MAXS],chain_prv[MAXS],name[MAXS]; 
 

 /* assuming `standard' PDB format - section SEQRES will be read here
    to collect all the non-identical chains (but only the consecutive
    ones are checked for redundancy), white spaces are assumed between
    residue names */

 /* to count number of non-redundant (nseq) and all (i_domain) chains */
 n_seq=0;
 i_domain=0;
 op->noredund=0;
 /* take the name of the query file as a name of your query */
 get_clean_name(name,F->f_exam);
 /* read first string */ 
 eof=fscanf(F->exam,"%s",ident);
 if(eof==EOF) return(0);
 /* initiate the counters for the subsequent (current and prv) chains */ 
 i_found_let=0;
 error=0; 
 strcpy(chain_prv,"");
 length_prv=0;
 sum=0;
 sum_prv=0;
 header=1;
 read_next=1;
 /* start reading the rest of PDB file */
 while(eof!=EOF && header)
   {

     /* read the whole line if not SEQRES section */
     if(strcmp(ident,"SEQRES")!=0)
       {

	 /* write down whatever was found before */
	 if(i_found_let>0)
	   {
	     if(i_found_let>=op->max_length || i_found_let<=SHORT)
	       {
		 extend_name(op,name,chain);
		 if(op->info) skip_query_warn(op,Info);
	       }
	     else if(i_found_let!=length_prv && sum!=sum_prv)
	       { 
		 extend_name(op,name,chain);
		 seq2query(F,op,Info,Seq,&i_found_let); 
		 ++n_seq;
		 length_prv=i_found_let;
		 sum_prv=sum;
	       }
	     ++i_domain;
	     sum=0;
	     i_found_let=0;
	     error=0;
	   }
	 /* read line to the end */	 
	 fgets(stmp,MAXS,F->exam);

       }
     else
       {
	 /* read row number and chain label first */
	 fscanf(F->exam,"%s %s",sym,chain);
	 /* assuming that MAXS>80 */
	 fgets(stmp,MAXS,F->exam); 
	 /* if another chain has started write down the previous one */
	 /* i_found_let is now length of the previous chain to be written
	    in this cycle, whereas length_prv is length of the chain before
	    the previous one - if new chain then consider two previous ones */
	 if((i_found_let>0) && (strcmp(chain,chain_prv)!=0))
	   {
	     if(i_found_let>=op->max_length || i_found_let<=SHORT)
	       {
		 /* printf("Curr_ch=%s Short prev_ch=%s l_prv=%d l_bef=%d\n",
			chain,chain_prv,i_found_let,length_prv); */
		 extend_name(op,name,chain_prv);
		 if(op->info) skip_query_warn(op,Info);		 
	       }
	     else if(i_found_let!=length_prv && sum!=sum_prv)
	       {
		 /* printf("Curr_ch=%s  Prev_ch=%s l_prv=%d l_before=%d\n",
			chain,chain_prv,i_found_let,length_prv); */
		 extend_name(op,name,chain_prv);
		 seq2query(F,op,Info,Seq,&i_found_let); 
		 ++n_seq;
		 length_prv=i_found_let;
		 sum_prv=sum;
	       }
	     ++i_domain;
	     /* sum_prv=sum; length_prv=i_found_let; */
	     sum=0;
	     i_found_let=0;
	     error=0;
	   }
	 /* get ama residues from SEQRES line */
	 error=line2seq(Seq,&i_found_let,stmp,&sum,&read_next);

       }

     if(i_found_let>0) strcpy(chain_prv,chain);
     eof=fscanf(F->exam,"%s",ident); /* start reading next line */
     if(strcmp(ident,"ATOM")==0) header=0; /* stop when ATOM */
     
   }
 if(n_seq==i_domain) op->noredund=1;
 
 return(n_seq);
 
}

/* to prepare standard f_query file if using PDB format */
int prepare_xyz_pdb(a_file *F, k_opt *op, k_info *Info, k_res *Res,
		    int *Seq, int *n_domain)
{
 int n_seq,i,i_length,i_found_let,j,eof,error,i_domain,n_res,length_prv,
     sum,sum_prv,first_mod,n_sc_last,l_seq,read_next,n_written;
 char ident[MAXS],stmp[MAXS],sym[MAXS],chain[MAXS],chain_prv[MAXS],name[MAXS],
      chain_curr[MAXS],chain_incl[MAXS]; 
 float x,y,z;

 /* assuming fixed (`standard') PDB format */
 /* one passage for each domain to read its sequence (once reached) and 
    atomic parameters of the current domain: n_seq goes from 0 to i_domain */
 /* ignoring all alternatives (for an atom, side chain, model etc)
    except for the first one */

 /* take the name of the query file as a name of your query */
 get_clean_name(name,F->f_exam);
 i_domain=1; 
 n_seq=0;
 n_written=0;
 strcpy(chain_prv,"");

 /* i_domain is allowed to increment only once in each passage */
 /* ignore domain if it appears to be (sum, length) the same as previous one */
 while(n_written<*n_domain)
   {
     /* printf("Begin next domain i_dom=%d n_seq_prv=%d n_dom=%d\n",
	    i_domain,n_seq,*n_domain); */
     n_seq=0;
     read_next=1;
     strcpy(chain_prv,"");
     restart_counters(&sum,&i_found_let,&error);
     sum_prv=0; 
     length_prv=0;
     n_res=0;
     init_Res(Res,&n_res); 
     first_mod=1;
     read_next=1;
     rewind(F->exam);
     /* read first token */
     eof=fscanf(F->exam,"%s",ident);
     if(eof==EOF) return(0); 
     /* start passage through the whole PDB file for each domain */
     /* ignore domains read previously and to be read in the next passage */
     while(eof!=EOF && first_mod)
       {

	 /* read the whole line if not SEQRES section */
	 if(strcmp(ident,"SEQRES")!=0)
	   {

	     /* define the current chain if found before */
	     if(i_found_let>0 && n_seq<i_domain)
	       {
		 if((i_found_let<op->max_length) && (i_found_let>SHORT) &&
		    (i_found_let!=length_prv) && (sum!=sum_prv))
		   {
		     extend_name(op,name,chain);
		     ++n_seq;
		     l_seq=i_found_let;
		     sum_prv=sum;
		     strcpy(chain_incl,chain_prv);
		     length_prv=i_found_let;
		     if(n_seq==i_domain) 
		       {
			 read_next=0; 
			 ++n_written;
		       }
		   }
		 else
		   {
		     ++n_seq;
		     read_next=1;
		   }
		 restart_counters(&sum,&i_found_let,&error);
	       }
	     /* read line to the end */	 
	     fgets(stmp,MAXS,F->exam);
	 
	   }

	 /* read and decompose the whole line to get residues if SEQRES */
	 if(strcmp(ident,"SEQRES")==0)
	   {
	     /* read row number and chain label */
	     fscanf(F->exam,"%s %s",sym,chain);
	     fgets(stmp,MAXS,F->exam); /* assuming that MAXS>80 */
	     /* go on till first domain still not read is reached */
	     /* printf("Trying to reach next chain prv=%s curr=%s n_seq= %d\n",
		    chain_prv,chain,n_seq); */
	     if(n_seq<i_domain)
	       {
		 /* if another chain has started define the current chain
		    as the previous one (the one which just ended) */ 
		 if((i_found_let>0) && (strcmp(chain,chain_prv)!=0))
		   {
		     if((i_found_let<op->max_length) && (i_found_let>SHORT) &&
			(i_found_let!=length_prv) && (sum!=sum_prv))
		       {
			 /* printf("CURR_ch=%s PRV_ch=%s l_prv=%d l_bef=%d\n",
				chain,chain_prv,i_found_let,length_prv); */
			 extend_name(op,name,chain_prv);
			 ++n_seq; 
			 l_seq=i_found_let;
			 /* printf("Chain_to_be_writt=%s n_seq=%d nwr_prv=%d\n",
				chain_prv,n_seq,n_written); */
			 strcpy(chain_incl,chain_prv);
			 strcpy(chain_prv,chain);
			 sum_prv=sum;
			 length_prv=i_found_let;
			 /* next non-ignored domain reached - 
			    write it down in the current cycle */
			 if(n_seq==i_domain) 
			   {
			     read_next=0; 
			     ++n_written;
			   }
			 
		       }
		     else
		       {
			 /* printf("Ignoring chain= %s n_seq= %d i_dom= %d\n",
				chain_prv,n_seq,i_domain); */
			 ++n_seq;
			 read_next=1;
		       }
		     
		     restart_counters(&sum,&i_found_let,&error);
		   }
	 
		 error=line2seq(Seq,&i_found_let,stmp,&sum,&read_next);
		 
	       }	     
	   }
	 
	 /* decompose line to get atom parameters if ATOM */
	 /* now chain_incl is the name of the chain to be included */
	 /* if there is only one chain n_domain=1 */
	 /* assuming fixed format here, rest of the line read already */
	 if((n_seq==i_domain) && (strcmp(ident,"ATOM")==0))
	   {
	     /* if there was only one chain (no chain_name char) */
	     if(*n_domain==1 && op->noredund) strcpy(chain_curr,chain_incl);
	     else
	       {
		 strncpy(chain_curr,&stmp[17],1); 
		 chain_curr[1]='\0'; 
	       }
	     if(strcmp(chain_curr,chain_incl)==0 && !read_next)
	       {
		 /* printf("seq=%d domain=%d \n",n_seq,i_domain); */
		 n_sc_last=get_line_pdb(F,op,Info,Res,Seq,&n_res,stmp);
		 
	       }
	     
	   } 

	 if(i_found_let>0) strcpy(chain_prv,chain); 
	 /* start reading next line */
	 eof=fscanf(F->exam,"%s",ident); 
	 /* stop when second model starts (NMR) */
	 if(strcmp(ident,"ENDMDL")==0) first_mod=0; 

       } /* end of one passage over pdb file */
     
     /* now write down the collected coordinates of a chain (domain) */
     if(n_seq==i_domain && !read_next)
       {
	 check_Res(Res,&n_res,&n_sc_last);
	 write_down_xyz(F,op,Info,Res,Seq,&l_seq,&n_res);
       }
     
     ++i_domain;
     
   } /* end of loop over domains with complete passage for each domain */
 
 
 return(1);
 
}

/* write_down_xyz tries to match assumed and actual seqs and writes coor */
int write_down_xyz(a_file *F,k_opt *op,k_info *Info,k_res *Res,
		   int *Seq,int *l_seq,int *n_res)
{
  int i,i_seq,i_res,s1,s2,s3,r1,r2,r3,i_rmax,next;
  
  /* try to match triplets i,i+1,i+2 where i+1,i+3 may be zeros 
     if the next residue is not included in PDB (i.e. in Res) */

  /* for(i=0;i<*n_res;++i) 
    {
       
      printf("%s %s %7.2f %7.2f %7.2f \n",
	     (Res+i+1)->name,(Res+i+1)->cnum,
	     (Res+i+1)->x_sc,(Res+i+1)->y_sc,(Res+i+1)->z_sc);
    } */ /* one may produce a clean pdb file above */
  sprintf(Info->Sout,"%s\t%d\n",op->qname,*l_seq);
  to_file(F->qxyz,Info->Sout);
  i_rmax=0;
  i_res=1; /* first residue is not in zeroth  */
  for(i_seq=0;i_seq<*l_seq;++i_seq) 
    {
      
      get_triplet(&s1,&s2,&s3,Seq,l_seq,&i_seq);
      /* get also current Res triplet */
      if(i_res>*n_res) r1=-1;
      else r1=seq2dig((Res+i_res)->name);
      if(i_res+1>*n_res) r2=-1;
      else
	{
	  next=(Res+i_res+1)->num;
	  if((Res+i_res)->num == (next-1))
	    r2=seq2dig((Res+i_res+1)->name);
	  else r2=-1;
	}
      if(i_res+2>*n_res || r2==-1) r3=-1;
      else
	{      
	  next=(Res+i_res+2)->num;
	  if((Res+i_res+1)->num == (next-1))
	    r3=seq2dig((Res+i_res+2)->name);
	  else r3=-1;
	}
      
      /*
      printf("s1=%d s2=%d s3=%d r1=%d r2=%d r3=%d #s=%d #r=%d \n",
	     s1,s2,s3,r1,r2,r3,i_seq+1,(Res+i_res)->num); */
	     

      /* check one whether the two triplets match each other */
      if(match_res(&s1,&s2,&s3,&r1,&r2,&r3))
	{
	  write_xyz_line(F,op,Info,Res,&i_res);
	  ++i_res;
	}
      else
	write_xyz_line(F,op,Info,Res,&i_rmax);
	     
    }
	 
	 
 return(1);  
}

int get_triplet(int *s1,int *s2,int *s3,int *Seq,int *l_seq,int *i_seq)
{

  *s1=*(Seq+*i_seq);
  if((*i_seq+1)<*l_seq) *s2=*(Seq+*i_seq+1);
  else *s2=-1;
  if((*i_seq+2)<*l_seq) *s3=*(Seq+*i_seq+2);
  else *s3=-1;      

  return(1);
  
}

int match_res(int *s1,int *s2,int *s3,int *r1,int *r2,int *r3)
{
  int check;
  
  /*
  printf("s1=%d s2=%d s3=%d r1=%d r2=%d r3=%d \n",
	 *s1,*s2,*s3,*r1,*r2,*r3); 
  if(*s1 != *r1) printf("disagree\n"); */
  
  
  if(*s1 != *r1) return(0);

  check=1;
  if((*s2<0) || (*r2<0)) check=0; 
  if((*s2 != *r2) && check)  return(0);

  check=1;
  if((*s3<0) || (*r3<0)) check=0;
  if((*s3 != *r3) && check)  return(0);
  
  return(1);
  
}

void check_Res(k_res *Res,int *n_res,int *n_sc_last)
{
  int i;
  float fn_sc;

  /* put RMAX into the first (zeroth) to be used when convenient */
  (Res)->x_a=RMAX; 
  (Res)->y_a=RMAX;
  (Res)->z_a=RMAX;
  (Res)->x_sc=RMAX;
  (Res)->y_sc=RMAX;
  (Res)->z_sc=RMAX;
  (Res)->x_b=RMAX;
  (Res)->y_b=RMAX;
  (Res)->z_b=RMAX;

  /* now from the second to the last before last */
  for(i=1;i<*n_res+1;++i) 
    {  
      if((Res+i)->x_a==0.0) (Res+i)->x_a=RMAX; 
      if((Res+i)->y_a==0.0) (Res+i)->y_a=RMAX;
      if((Res+i)->z_a==0.0) (Res+i)->z_a=RMAX;
      if((Res+i)->x_sc==0.0)(Res+i)->x_sc=RMAX;
      if((Res+i)->y_sc==0.0) (Res+i)->y_sc=RMAX;
      if((Res+i)->z_sc==0.0) (Res+i)->z_sc=RMAX;
      if((Res+i)->x_b==0.0) (Res+i)->x_b=RMAX;
      if((Res+i)->y_b==0.0) (Res+i)->y_b=RMAX;
      if((Res+i)->z_b==0.0) (Res+i)->z_b=RMAX;      
    }

  /* rescale now the side chain of the last */
  if(*n_sc_last == 0) *n_sc_last=1;
  fn_sc=(float)*n_sc_last;
  (Res+*n_res)->x_sc=((Res+*n_res)->x_sc/fn_sc); 
  (Res+*n_res)->y_sc=((Res+*n_res)->y_sc/fn_sc); 
  (Res+*n_res)->z_sc=((Res+*n_res)->z_sc/fn_sc);

}

void init_Res(k_res *Res,int *n_res)
{
  
  (Res+*n_res)->num=0;
  (Res+*n_res)->n_sc=0;
  (Res+*n_res)->x_a=0.0; 
  (Res+*n_res)->y_a=0.0;
  (Res+*n_res)->z_a=0.0;
  (Res+*n_res)->x_sc=0.0;
  (Res+*n_res)->y_sc=0.0;
  (Res+*n_res)->z_sc=0.0;
  (Res+*n_res)->x_b=0.0;
  (Res+*n_res)->y_b=0.0;
  (Res+*n_res)->z_b=0.0;
  strcpy((Res+*n_res)->name,"");
  strcpy((Res+*n_res)->cnum,"");
  (Res+*n_res)->secstr=0; /* undefined secondary structure */ 
  /* (Res+*n_res)->alt=0; */
}

/* line2seq tokenizes string to get ama residues */
int line2seq(int *Seq,int *i_let,char *line,int *sum,int *read_next)
{
  int error,j;
  char *token_ptr;
  
  error=0;
  /* now the first ama or number of res comes */
  token_ptr=strtok(line," ");
  while(token_ptr!=0)
    {
      j=seq2dig(token_ptr);
      if(j<0) error=1;  
      /* wrong symbol - just skip it this time */
      else
	{
	  if(*read_next) *(Seq+*i_let)=j;
	  ++(*i_let); 
	  *sum+=j;
	}
      token_ptr=strtok(NULL," ");
    }
  return(error);
  
}

/* this is just small function to avoid typing */
void restart_counters(int *sum,int *i_let,int *error)
{
  *sum=0;
  *i_let=0;
  *error=0;  
}

/* get PDB coordinates and tranform them into LOOPP XYZ file */
int get_line_pdb(a_file *F,k_opt *op,k_info *Info,k_res *Res,
		 int *Seq,int *n_res,char *line)
{
 int i,n,beg,n_dig,i_res,pick_line,alternative;
 char xtmp[MAXS],Res_num[10],Res_num_prv[10],Res_name[4],Atom[5],copy;
 float x_sc,y_sc,z_sc,x_a,y_a,z_a,x_b,y_b,z_b,fn_sc;
 
 /* line is supposed to contain complete PDB line, except for ATOM */
 /* n_res points to the residue to be filled and it is incremented when
    it is found that a new residue started except for the first one */

 beg=0; /* could be a parameter */
 pick_line=1;
 alternative=0;
 
 /* get the residue number and name, atom name  */
 strncpy(Res_num,&line[beg+18],4); Res_num[4]='\0';
 strncpy(Res_name,&line[beg+13],3); Res_name[3]='\0';

 /* skip alternative coordinates, wrong residues aso */
 if(seq2dig(Res_name)<0) pick_line=0;
 copy=line[beg+12];
 if(copy!=' ' && copy!='A' && copy!='1') pick_line=0;
 strncpy(Atom,&line[beg+9],4); Atom[4]='\0';
 copy=line[beg+8];
 if(copy!=' ' && copy!='A' && copy!='1') pick_line=0;
 /* check for alternative residues - take the last */
 copy=line[beg+22];
 if(copy!=' ')
   {
     alternative=1;
     if(copy=='A' || copy=='1') strcpy(Res_num,"9998");
     else strcpy(Res_num,"9999");
   }
 

 /* the first coordinate column should start at 26th column */
 beg=26;
 n_dig=8;
 
 if(pick_line) 
   {    
     /* printf("%s #curr=%s #prv=%s\n",Res_name,Res_num,(Res+*n_res)->cnum);*/
     /* if next residue reached finish the previous and re-initiate */
     /* if(strncmp(Res_num_prv,Res_num,4)!=0) */
     if(strncmp(Res_num,(Res+*n_res)->cnum,4)!=0)
       {
	 if((Res+*n_res)->n_sc == 0) (Res+*n_res)->n_sc=1;
	 if(*n_res>0) fn_sc=(float)(Res+*n_res)->n_sc;
	 else fn_sc=1.0;
	 (Res+*n_res)->x_sc=((Res+*n_res)->x_sc/fn_sc); 
	 (Res+*n_res)->y_sc=((Res+*n_res)->y_sc/fn_sc); 
	 (Res+*n_res)->z_sc=((Res+*n_res)->z_sc/fn_sc);

	 /* go to next residue with the current data */
	 ++(*n_res); 
	 (Res+*n_res)->x_sc=0.0;
	 (Res+*n_res)->y_sc=0.0;
	 (Res+*n_res)->z_sc=0.0;
	 (Res+*n_res)->n_sc=0;	 

	 /* now you may substitute what you already know about new res */
	 (Res+*n_res)->num=atoi(Res_num);
	 
	 /* going back dropped of the plan */
	 /*
	 if(*n_res>0) (Res+*n_res-1)->plmin_1=seq2dig(Res_name);
	 if(*n_res>1) (Res+*n_res-2)->plmin_2=seq2dig(Res_name);
	 */
	 strcpy((Res+*n_res)->name,Res_name);
	 strcpy((Res+*n_res)->cnum,Res_num);

       }
     /*
     if(alternative) (Res+*n_res)->alt=1; 
     else
       { */
	 
     /* now you may continue reading new residue */  
     if(pick_a(Res_name,Atom)) 
       {
	 get_xyz_col(line,&beg,&x_a,&y_a,&z_a,&n_dig);
	 (Res+*n_res)->x_a=x_a; 
	 (Res+*n_res)->y_a=y_a;
	 (Res+*n_res)->z_a=z_a;	 
       }
     
     if(pick_b(Res_name,Atom)) 
       {
	 get_xyz_col(line,&beg,&x_b,&y_b,&z_b,&n_dig);   
	 (Res+*n_res)->x_b=x_b;
	 (Res+*n_res)->y_b=y_b;
	 (Res+*n_res)->z_b=z_b;	 
       }

     if(pick(Res_name,Atom))
       {
	 get_xyz_col(line,&beg,&x_sc,&y_sc,&z_sc,&n_dig);
	 (Res+*n_res)->x_sc+=x_sc;
	 (Res+*n_res)->y_sc+=y_sc;
	 (Res+*n_res)->z_sc+=z_sc;
	 ++((Res+*n_res)->n_sc);
       }
     
     if(op->debug)
       {
	 strncpy(xtmp,&line[beg],3*n_dig); xtmp[3*n_dig]='\0';
	 sprintf(Info->Sout,"%s %s %s %s ",Atom,Res_name,Res_num,xtmp);
	 to_stnd_out(Info->Sout);
	 if(pick_a(Res_name,Atom)) 
	   sprintf(Info->Sout," %7.2f %7.2f %7.2f\n",x_a,y_a,z_a);
	 else sprintf(Info->Sout,"\n");
	 to_stnd_out(Info->Sout);
       }

     /*  } */
     
   }

 /* in case of last residue one needs n_sc later on */
 return((Res+*n_res)->n_sc);
	 
}

/* having all residue informations write it down to f_qxyz */
void write_xyz_line(a_file *F,k_opt *op,k_info *Info,k_res *Res,int *i_res)
{
  char fmt[MAXS];
  int i;
  
 /* fix the format of XYZ file */
 if(op->n_col==1)
   strcpy(fmt," %10.3f %10.3f %10.3f\n");
 if(op->n_col==2)
   strcpy(fmt," %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n");
 if(op->n_col==3)
   strcpy(fmt,
	  " %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n");

 /* remember that residues in Res start from 1 and not from 0 */
 /* i=*i_res+1;  */ /* now it is taken care outside */
 i=*i_res;
 
 if(op->n_col==1)
   sprintf(Info->Sout,fmt,(Res+i)->x_a,(Res+i)->y_a,(Res+i)->z_a);
 if(op->n_col==2)
   sprintf(Info->Sout,fmt,(Res+i)->x_a,(Res+i)->y_a,(Res+i)->z_a,
	   (Res+i)->x_sc,(Res+i)->y_sc,(Res+i)->z_sc);
 if(op->n_col==3)
   sprintf(Info->Sout,fmt,(Res+i)->x_a,(Res+i)->y_a,(Res+i)->z_a,
	   (Res+i)->x_sc,(Res+i)->y_sc,(Res+i)->z_sc,
	   (Res+i)->x_b,(Res+i)->y_b,(Res+i)->z_b);
 to_file(F->qxyz,Info->Sout);
  
}

/* get_xyz_col reads x,y,z coordinates from PDB or CRD line */
void get_xyz_col(char *line,int *beg,float *x,float *y,float *z,int *n_dig)
{
  char xtmp[MAXS];
  int start;
  
  start=*beg;
  
  /* assumption: n_dig < MAXS  */
  strncpy(xtmp,&line[start],*n_dig); xtmp[*n_dig]='\0';
  *x=(float)atof(xtmp); 
  start+=(*n_dig);
  strncpy(xtmp,&line[start],*n_dig); xtmp[*n_dig]='\0';
  *y=(float)atof(xtmp); 
  start+=(*n_dig);
  strncpy(xtmp,&line[start],*n_dig); xtmp[*n_dig]='\0';
  *z=(float)atof(xtmp); 
}

void get_clean_name(char *name,char *file_name)
{
 int j;
 char stmp[MAXS];
 
 /* ignore path  */
 if(strrchr(file_name,'/') == NULL) strcpy(stmp,file_name);
 else strcpy(stmp,strrchr(file_name,'/')+1);
 /* if NT path */
 if(strrchr(file_name,'\\') != NULL) 
   strcpy(stmp,strrchr(file_name,'\\')+1);
 /* ignore extension */ 
 j=strcspn(stmp,"."); 
 strncpy(name,stmp,j);   
 name[j]='\0';
 
}

/* add the chain id to the name and pass it to op->qname */
void extend_name(k_opt *op,char *name,char *chain_name)
{
 int j;
 char stmp[MAXS];
 
 strcpy(op->qname,name);
 if(strlen(chain_name)==1)
   {
     strcat(op->qname,"_");
     strcat(op->qname,chain_name);
   }
}

void skip_query_warn(k_opt *op,k_info *Info)
{
  
  sprintf(Info->Sout,
	  " WARNING - skipping (too long or too short) sequence %s \n",
	  op->qname);
  to_stnd_out(Info->Sout);

}

/* to prepare standard f_query file if using 1-let inp */
void prepare_query_1let(a_file *F, k_opt *op, k_info *Info)
{
 int n_seq,i,i_length,i_found_let,j,*Lngs,let;
 char buffer[MAX],*w; 
 

 /* assuming that seqs are using only generic (one-lett.) symbols and are 
    separated by a hash (#) symbol, white spaces ignored */

 w=buffer; 

 /* first - read to get number of sequences */
 n_seq=0;
 j='\0'; 
 while(j!=EOF)
   {
     j=fgetc(F->exam); 
     i_length=0;
     while(j!='#' && j!=EOF) 
       {
	 if(j!=' ' && j!='\n')
	   {
	     ++i_length; 
	   }
	 j=fgetc(F->exam); 
       }
     if(i_length>0) ++n_seq;
     
   }
 /* now one may allocate memory */
 Lngs=malloc(n_seq*sizeof(int));
 if(Lngs==NULL) alloc_err("Lengths"); 


 /* second - read to get lengths */
 rewind(F->exam);
 i=0;
 j='\0'; 
 while(j!=EOF)
   {
     j=fgetc(F->exam); *w=j;
     i_length=0;
     while(j!='#' && j!=EOF) 
       {
	 if(j!=' ' && j!='\n')
	   {
	     let=one2dig(w);
	     if(let>=0) ++i_length; 
	   }
	 j=fgetc(F->exam); *w=j;
       }
     if(i_length>0)
       {
	 /* if(ilength>op->max_length) reallocate(); */
	 /* right now quick patch */
	 if(i_length>op->max_length)
	   {
	     sprintf(Info->Sout,"Sequence %s is too long \t max_len=%d\n",
		     op->qname,op->max_length);
	     to_file(F->best,Info->Sout);
	     to_stnd_out(Info->Sout);
	     exit(1);
	   }
	 *(Lngs+i)=i_length; 
	 ++i;
       }
     
   }

 /* next rewind and read again to write down */
 rewind(F->exam);
 i=0;
 while(i<n_seq)
   {
     i_length=*(Lngs+i);
     if(i_length<op->max_length) 
       {
	 if(n_seq==1)
	   sprintf(Info->Sout,"%s\t%d \n",op->qname,i_length);
	 else
	   sprintf(Info->Sout,"%s_%d\t%d \n",op->qname,i+1,i_length);
	 to_file(F->query,Info->Sout);
	 /* read until next non-empty */
	 i_found_let=0;
	 j='\0';
	 while(i_found_let==0 && j!=EOF) 
	   {
	     j=fgetc(F->exam); *w=j; 
	     while(j!='#' && j!=EOF) 
	       {
		 if(j!=' ' && j!='\n')
		   {
		     let=one2dig(w);
		     /* skip wrong symbols */
		     if(let>=0)
		       {
			 sprintf(Info->Sout,"%s ",dig2amino(let));
			 if(((i_found_let+1)%20)== 0) strcat(Info->Sout,"\n");
			 to_file(F->query,Info->Sout);
			 ++i_found_let; 
		       }
		   }
		 j=fgetc(F->exam); *w=j; 
	       }
	   }
	  ++i;
	  sprintf(Info->Sout,"\n");
	  to_file(F->query,Info->Sout);
       }
     else /* skip sequence if too long */
       {
	 i_found_let=0;
	 j='\0';
	 while(i_found_let==0 && j!=EOF) 
	   {
	     j=fgetc(F->exam); 
	     while(j!='#' && j!=EOF) 
	       {
		 if(j!=' ' && j!='\n')
		   {
		     ++i_found_let; 
		   }
		 j=fgetc(F->exam); 
	       }
	   }
	  ++i;
       }

   }


}


/* to prepare standard f_query file if using name and 1-let inp */
void prepare_query_fasta(a_file *F, k_opt *op, k_info *Info)
{
 int n_seq,i,i_length,i_found_let,j,*Lngs,w,let;
 char buffer[MAXS],*ctmp; 
 

 /* assuming that seqs are using only generic (one-lett.) symbols and are 
    separated by a larger than (>) symbol proceeded with a name, 
    white spaces ignored, name of a sequence ends when end of line  */

 ctmp=buffer; 

 /* first - read to get number of sequences */
 n_seq=0;
 w=fgetc(F->exam); 
 while(w!=EOF)
   {
     
     /* if character > found then read till end of line */
     if(w=='>')
       {
	 i=0;
	 w=fgetc(F->exam);
	 *ctmp=w; /* conversion to unsigned */
	 while(w!='\n')
	   {
	     *(op->qname+i)=*ctmp;
	     ++i;
	     w=fgetc(F->exam);
	     *ctmp=w;
	   }
	 *(op->qname+i)='\0';
       }
     /* printf("%s \n",op->qname); */
     
     i_length=0;
     while(w!='>' && w!=EOF) 
       {
	 if(w!=' ' && w!='\n')
	   {
	     ++i_length; 
	   }
	 w=fgetc(F->exam);  
       }
     if(i_length>0) ++n_seq;
     
   } 
 /* now one may allocate memory */
 Lngs=malloc(n_seq*sizeof(int));
 if(Lngs==NULL) alloc_err("Lengths"); 


 /* second - read to get lengths */
 rewind(F->exam);
 i=0;
 w=fgetc(F->exam); 
 while(w!=EOF)
   {
     /* if character > found then read till end of line */
     if(w=='>')
       {
	 j=0;
	 w=fgetc(F->exam);
	 *ctmp=w;
	 while(w!='\n')
	   {
	     *(op->qname+j)=*ctmp;
	     ++j;
	     w=fgetc(F->exam);
	     *ctmp=w;
	   }
	 *(op->qname+j)='\0';
       }
     
     i_length=0;
     while(w!='>' && w!=EOF) 
       {
	 if(w!=' ' && w!='\n')
	   {
	     let=one2dig(ctmp);
	     if(let>=0) ++i_length; 
	   }
	 w=fgetc(F->exam); *ctmp=w; 
       }
     if(i_length>0)
       {
	 *(Lngs+i)=i_length; 
	 ++i;
       }
     
   }

 /* next rewind and read again to write down */
 rewind(F->exam);
 i=0;
 i_found_let=0;
 w=fgetc(F->exam); 
 *ctmp=w;
 while(w!=EOF) /* read empty seqs too */
   {
	     
     /* if character > found then read till end of line */
     if(w=='>')
       {
	 j=0;
	 w=fgetc(F->exam);
	 while(w!='\n')
	   {
	     if(w==' ') w='_';
	     *ctmp=w;
	     *(op->qname+j)=*ctmp;
	     ++j;
	     w=fgetc(F->exam);
	     *ctmp=w;
	   }
	 *(op->qname+j)='\0';
	 i_length=*(Lngs+i);
	 if(i_length<op->max_length) 
	   {
	     sprintf(Info->Sout,"%s\t%d \n",op->qname,i_length);
	     if(op->pick)
	       {
		 if(op->pick_all_seq) strcpy(op->pseq_name,op->qname);
		 if(strncmp(op->qname,op->pseq_name,strlen(op->pseq_name))==0)
		   to_file(F->query,Info->Sout);
	       }
	     else to_file(F->query,Info->Sout);
	   }
     
       }
     /* now read the sequence and write it down */
     while(w!='>' && w!=EOF) 
       {
	 if(w!=' ' && w!='\n')
	   {
	     if(i_length<op->max_length) 
	       {
		 let=one2dig(ctmp);
		 /* skip wrong symbols */
		 if(let>=0)
		   {
		     sprintf(Info->Sout,"%s ",dig2amino(let));
		     if(((i_found_let+1)%20)== 0) strcat(Info->Sout,"\n");
		     if(op->pick)
		       {
			 if(strncmp(op->qname,op->pseq_name,
				    strlen(op->pseq_name))==0)
			   to_file(F->query,Info->Sout);
		       }
		     else to_file(F->query,Info->Sout);
		     ++i_found_let;
		   }
		 
	       }
	     
	   }
	 w=fgetc(F->exam); /* with this go to the next let or seq */
	 *ctmp=w;
       }
     ++i;
     i_found_let=0;
     if(i_length<op->max_length) 
       {
	 sprintf(Info->Sout,"\n");
	 if(op->pick)
	   {
	     if(strncmp(op->qname,op->pseq_name,strlen(op->pseq_name))==0)
	       to_file(F->query,Info->Sout);
	   }
	 else to_file(F->query,Info->Sout);
       }
	 
   }


}




/* picks the side chain atoms to calculate residue center */
int pick(char *amino,char *atom)
{
 if(strncmp("GLY",amino,3)==0)
   {
    if(atom[0]=='C' && atom[1]=='A') return(1);
    else
      return(0);
   }
 if(strncmp("ALA",amino,3)==0)
   {
    if(atom[0]=='C' && atom[1]=='B') return(1);
    else
      return(0);
   }
 if(atom[0]=='H') return(0);
 if(atom[0]=='C' && atom[1]=='A') return(0);
 if(atom[0]=='N' && atom[1]==' ') return(0);
 if(atom[0]=='C' && atom[1]==' ') return(0);
 if(atom[0]=='O' && atom[1]==' ') return(0);
 return(1);
}

/* this function picks CA atoms */
int pick_a(char *amino,char *atom)
{
  if(atom[0]=='C' && atom[1]=='A') return(1);
  else return(0);

}

/* this function picks CB atoms */
int pick_b(char *amino,char *atom)
{
  if(atom[0]=='C' && atom[1]=='B') return(1);
  else return(0);

}

/*                                                                       */
/* a couple of obsolete functions that may be still useful in the future */
/*                                                                       */

/* prepare standard f_seq file with 3-lett. codes if using 1-lett. inp */
void one2three(a_file *F,k_opt *op)
{
 int m,Fseq_1let,n_seq,i,j,*Lngs;
 char T_name[MAXS],name[MAXS],buffer[MAX],*w; 
 FILE *onelet;
 

 /* assuming that seqs are using only generic (one-lett.) symbols and are 
    stored one seq per line in the f_seq file */

 Fseq_1let=open(F->f_seq,O_RDONLY);
 strcpy(name,F->f_seq);
 strcpy(T_name,strcat(name,"_3let"));
 F->seq=fopen(T_name,"w");
 if(F->seq==NULL) open_err(T_name);
     
 /* assuming that each line contains whole sequence */
 /* assuming CRD-like type of header with number of lines coming first */
 n_seq=get_header_crd(Fseq_1let);
 w=buffer; /* patch - MAX becomes a limitation for a length of seq */
 Lngs=malloc(n_seq*sizeof(int));
 if(Lngs==NULL) alloc_err("Lngs"); 
 i=0;
 while(i<n_seq)
   {
     *w='\0';
     j=0;
     while(*w!='\n') 
       {
	 read(Fseq_1let,w,1); ++j;
	 if(op->debug) printf("%s",w);
       }
     --j;
     *(Lngs+i)=j;
     ++i;
   }

 /* close the file and read it again to add header */
 close(Fseq_1let);
 Fseq_1let=open(F->f_seq,O_RDONLY); 
 /* onelet=fopen(F->f_seq,"r"); */
 n_seq=get_header_crd(Fseq_1let); 
 /* fscanf(onelet,"%d\n",&n_seq); */
 i=0;
 while(i<n_seq)
   {
     fprintf(F->seq,"seq_%d\t%d \n",i+1,*(Lngs+i));
     *w='\0';
     read(Fseq_1let,w,1); 
     /* j=fgetc(onelet); *w=j; */
     while(*w!='\n') 
       {
	 fprintf(F->seq,"%s \n",dig2amino(one2dig(w)));
	 read(Fseq_1let,w,1);  
	 /* j=fgetc(onelet); *w=j; */
       }
     ++i;
   }

 close(Fseq_1let); 
 fclose(F->seq);
 strcpy(F->f_seq,T_name); /* now replace the old file */

}


/* get_red_repr allows using the program as an interface to get reduced repr */
void get_red_repr(a_file *F,k_opt *op,float *X,float *Y,float *Z,int *seq)
{
 int Fcrd,Fpdb,n,i,n_lines; 
 char name_xyz[MAXS],name_seq[MAXS],name[MAXS],fmt[MAXS];


 strcpy(fmt," %10.3f %10.3f %10.3f \n");

 if(op->crd2xyz)
   {
     Fcrd=open(F->f_crd,O_RDONLY);

     /* read the header first to allocate */
     n_lines=get_header_crd(Fcrd);
     seq=malloc(n_lines*sizeof(int)); /* n_mono < n_lines */
     if(seq==NULL) alloc_err("seq"); 
     X=malloc(n_lines*sizeof(float));
     if(X==NULL) alloc_err("X"); 
     Y=malloc(n_lines*sizeof(float));
     if(Y==NULL) alloc_err("Y"); 
     Z=malloc(n_lines*sizeof(float));
     if(Z==NULL) alloc_err("Z"); 
     close(Fcrd);
     Fcrd=open(F->f_crd,O_RDONLY);

     /* fix the names of xyz and seq files */
     strcpy(name,F->f_crd);
     strcpy(name_xyz,strcat(name,"_xyz"));
     strcpy(name,F->f_crd);
     strcpy(name_seq,strcat(name,"_seq"));
     F->coor=fopen(name_xyz,"w");
     F->seq=fopen(name_seq,"w");

     /* read the crd file */
     n=read_crd(1,Fcrd,F->coor,op,X,Y,Z,seq);

     /* start writing to xyz and seq */
     fprintf(F->seq,"%s\n\t%d\n",F->f_crd,n);
     fprintf(F->coor,"%s\n\t%d\n",F->f_crd,n);
     for(i=0;i<n;++i)
       fprintf(F->seq,"%s\n",dig2amino(*(seq+i)));
     /* use proper format depending on the number of 3-col blocks */
     if(op->n_col==1)
       {
	 for(i=0;i<n;++i)
	     fprintf(F->coor,fmt,*(X+i),*(Y+i),*(Z+i));
       }
     else 
       {
	 close(Fcrd);
	 Fcrd=open(F->f_crd,O_RDONLY);
	 n=read_crd(2,Fcrd,F->coor,op,X,Y,Z,seq);
       }     
   }
 
 if(op->pdb2xyz) 
   {
     printf(" Sorry - not implemented - use external interface \n");
     exit(0);
     Fpdb=open(F->f_pdb,O_RDONLY);
     /* you may now call read_pdb */
   }
 
 printf(" Converting into reduced represenation - see files %s, %s\n",
	name_seq,name_xyz);

 /* clean up and exit */
 close(Fcrd);
 free(seq);
 free(X); free(Y); free(Z);
 exit(0);  /* as it is just interface finish now */
 
}

/* this function replaces at present the one above */
/* prepare query files in case of CHARMM format - mixture of old and new code */
int prepare_query_crd(a_file *F, k_opt *op, k_info *Info,
		      k_prot *prot, int *Seq)
{
 int Fcrd,n,i,n_lines; 
 float *X, *Y, *Z;
 char name[MAXS],fmt[MAXS];

 /* prepare local shortcuts */
 X=prot->X_res;
 Y=prot->Y_res;
 Z=prot->Z_res;
 strcpy(fmt," %10.3f %10.3f %10.3f \n");

 fclose(F->exam); /* close to use this strange hybrid ... */
 Fcrd=open(F->f_exam,O_RDONLY);
 
 /* fix the name of the query sequence */
 get_clean_name(name,F->f_exam);
 strcpy(op->qname,name);

 /* check for max_length */
 n_lines=get_header_crd(Fcrd);
 close(Fcrd);
 n_lines+=10; /* add arbitrarily 10 lines for header */
 if(n_lines>op->max_length) return(0);

 /* read the crd file */
 Fcrd=open(F->f_exam,O_RDONLY);
 n=read_crd(1,Fcrd,F->qxyz,op,X,Y,Z,Seq);

 /* start writing to xyz and seq */
 /* fprintf(F->query,"%s\t%d\n",name,n); */

 seq2query(F,op,Info,Seq,&n); 
 
 
 /* for(i=0;i<n;++i) fprintf(F->query,"%s\n",dig2amino(*(seq+i))); */

 fprintf(F->qxyz,"%s\t%d\n",name,n);
 /* use proper format depending on the number of 3-col blocks */
 if(op->n_col==1)
   {
     for(i=0;i<n;++i) fprintf(F->qxyz,fmt,*(X+i),*(Y+i),*(Z+i));
   }
 else 
   {
     close(Fcrd);
     Fcrd=open(F->f_exam,O_RDONLY);
     n=read_crd(2,Fcrd,F->qxyz,op,X,Y,Z,Seq);
   }     

 /* restore previous state of files  */
 close(Fcrd); 
 F->exam=fopen(F->f_exam,"r");
 
 return(1);
 
}

/* read head of the CRD file to determine n_lines */
int get_header_crd(int Fcrd)
{
 int n_lines;
 char lnum[6],buf[83];
  
 /* use read_line for lines with varying length */
 read_line(Fcrd,buf);
 strncpy(lnum,&buf[0],1); lnum[1]='\0';
 if(strcmp(lnum,"*")!=0)
   {
     strncpy(lnum,&buf[0],5); lnum[5]='\0'; 
     n_lines=(int)atoi(lnum);     
   }
 else
   {
     while(strcmp(lnum,"*")==0)
       {
	 read_line(Fcrd,buf);
	 strncpy(lnum,&buf[0],1); lnum[1]='\0';	 
       }
     strncpy(lnum,&buf[0],5); lnum[5]='\0'; 
     n_lines=(int)atoi(lnum);      
   }
 /* printf("number of coor lines in CRD file %d \n",n_lines); */
 return(n_lines);
 
}


/* read_crd reads CRD coordinate file to get SEQ and XYZ */
int read_crd(int turn,int Fcrd,FILE *Fxyz,k_opt *op,
	     float *X,float *Y,float *Z,int *seq)
{
 int i,n,num,n_sc,n_lines;
 char cnum[10],Canum[10],Ranum[10],amino[4],atom[5],buf[83],fmt[MAXS];
 float x_sc,y_sc,z_sc,x_a,y_a,z_a,x_b,y_b,z_b;
 

 /* get the header first */
 n_lines=get_header_crd(Fcrd);

 /* fix the format of XYZ file */
 if(op->n_col==2)
   strcpy(fmt," %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n");
 if(op->n_col==3)
   strcpy(fmt," %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n");

 n=read(Fcrd,buf,71); 
 num=0; i=1; /* loosing last line - should be O */
 while(i<n_lines)
 {
   strncpy(Canum,&buf[7],4); Canum[4]='\0';
   strcpy(Ranum,Canum); 
   strncpy(amino,&buf[11],3); amino[3]='\0';
   if(seq2dig(amino)>=0)
     {
       x_sc=0.0; y_sc=0.0; z_sc=0.0;
       n_sc=0;
       while(strncmp(Ranum,Canum,4)==0)
	 {
	   strncpy(atom,&buf[16],4);atom[4]='\0';
	   if(pick_a(amino,atom))
	     {
	       strncpy(cnum,&buf[22],9); cnum[9]='\0';
	       x_a=(float)atof(cnum);  
	       strncpy(cnum,&buf[32],9); cnum[9]='\0';
	       y_a=(float)atof(cnum); 
	       strncpy(cnum,&buf[42],9); cnum[9]='\0';
	       z_a=(float)atof(cnum); 
	     }
	   if(pick_b(amino,atom))
	     {
	       strncpy(cnum,&buf[22],9); cnum[9]='\0';
	       x_b=(float)atof(cnum);  
	       strncpy(cnum,&buf[32],9); cnum[9]='\0';
	       y_b=(float)atof(cnum); 
	       strncpy(cnum,&buf[42],9); cnum[9]='\0';
	       z_b=(float)atof(cnum); 
	     }
	   if(pick(amino,atom))
	     {
	       strncpy(cnum,&buf[22],9); cnum[9]='\0';
	       x_sc+=(float)atof(cnum);
	       strncpy(cnum,&buf[32],9); cnum[9]='\0';
	       y_sc+=(float)atof(cnum);
	       strncpy(cnum,&buf[42],9); cnum[9]='\0';
	       z_sc+=(float)atof(cnum);
	       ++n_sc;
	     }
	   n=read(Fcrd,buf,71); 
	   ++i; 
	   strncpy(Ranum,&buf[7],4); Ranum[4]='\0';
	   /* patch !!! */
	   if(i==n_lines) strcpy(Ranum,"CTER");
	   
	 }
       if(x_sc==0.0 && y_sc==0.0 && z_sc==0.0)
	 {
	   printf("There was a problem while reading ama num %s\n",Canum);
	   x_sc=y_sc=z_sc=999.9;
	   n_sc=1;
	 }
       *(X+num)=(x_sc/(float)n_sc); 
       *(Y+num)=(y_sc/(float)n_sc); 
       *(Z+num)=(z_sc/(float)n_sc); 
       *(seq+num)=seq2dig(amino); 
       if(turn!=1)
	 {
	   if(op->n_col==2)
	     fprintf(Fxyz,fmt,x_a,y_a,z_a,*(X+num),*(Y+num),*(Z+num));

	   if(op->n_col==3)
	     fprintf(Fxyz,fmt,x_a,y_a,z_a,*(X+num),*(Y+num),*(Z+num),
		     x_b,y_b,z_b);
	 }
       
       ++num;
     }
   else /* if not ama */ 
     {
       n=read(Fcrd,buf,71);
       ++i;
     }
   
 }

 return(num);

}

/* this function simply reads line of a text file */
void read_line(int inp,char *w)
{
  w--;
  while(*w!='\n')
    {
      w++;
      read(inp,w,1);
      /* if(eof(inp)) exit(0); */
    }
  
}
