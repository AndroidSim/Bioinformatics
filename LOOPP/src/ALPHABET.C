/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000      v. 2.000                                 */
/*                                                                */
/*       DEFINITION OF AN ALPHABET AND TRANSLATION FUNCTIONS      */
 
#include "loopp.h"
   

/* One may add new symbols by a straightforward extension of the
   functions included below - there is no external file with the
   definition of a generic alphabet (yet) */

/* this function returns address of a generic ama in the current alphabet */
int generic2adr(int *i, int *n_type, int *Alph, int *Adr)
{
 int j,repl_hyd,repl_pol,repl_chg,repl_chg_neg,repl_all;

 /* check first replacements - use generic order of ama from seq2dig */
 repl_hyd=0;
 repl_pol=0;
 repl_chg=0;
 repl_chg_neg=0;
 repl_all=0;
 for(j=0;j<*n_type;++j)
   {
     if(*(Alph+j)==21) repl_hyd=1;
     if(*(Alph+j)==22) repl_pol=1;
     if(*(Alph+j)==23) repl_chg=1;
     if(*(Alph+j)==24) repl_chg_neg=1;
     if(*(Alph+j)==25) repl_all=1;
   }
  
 /* if symbol of ama was given explicitly in the definition of alphabet
    in Model, then its address is already in Adr, otherwise it has to be
    found */
 j=-1;
 if(*(Adr+*i)!=-1) j=*(Adr+*i);
 else
   {
     if(repl_all) j=*(Adr+25); /* if ALL is used */
     
     if(repl_hyd) /* ama that are hydrophobic in first place */
       {
	 if(*i==0) j=*(Adr+21); /* ALA */
	 if(*i==4) j=*(Adr+21); /* CYS */
	 if(*i==8) j=*(Adr+21); /* HIS */
	 if(*i==9) j=*(Adr+21); /* ILE */
	 if(*i==10) j=*(Adr+21); /* LEU */
	 if(*i==12) j=*(Adr+21); /* MET */
	 if(*i==13) j=*(Adr+21); /* PHE */
	 if(*i==14) j=*(Adr+21); /* PRO */
	 if(*i==17) j=*(Adr+21); /* TRP */
	 if(*i==18) j=*(Adr+21); /* TYR */
	 if(*i==19) j=*(Adr+21); /* VAL */
       }
     if(repl_pol) /* non-hydrophobic */
       {
	 if(*i==1) j=*(Adr+22); /* ARG */
	 if(*i==2) j=*(Adr+22); /* ASN */
	 if(*i==3) j=*(Adr+22); /* ASP */
	 if(*i==5) j=*(Adr+22); /* GLN */
	 if(*i==6) j=*(Adr+22); /* GLU */
	 if(*i==7) j=*(Adr+22); /* GLY */
	 if(*i==11) j=*(Adr+22); /* LYS */
	 if(*i==15) j=*(Adr+22); /* SER */
	 if(*i==16) j=*(Adr+22); /* THR */
       }
     if(repl_chg) /* charged */
       {
	 if(*i==1) j=*(Adr+23); /* ARG */
	 if(*i==3) j=*(Adr+23); /* ASP */
	 if(*i==6) j=*(Adr+23); /* GLU */
	 if(*i==11) j=*(Adr+23); /* LYS */	 
       } 
     if(repl_chg_neg) /* negatively charged */
       {
	 if(*i==3) j=*(Adr+24); /* ASP */
	 if(*i==6) j=*(Adr+24); /* GLU */
       }  
     /* if there is no replacement e.g. both LEU and HYD missing */
     if(j<0) 
       {
	 printf(" Sorry, symbol %s does not have replacement in Model\n",
		dig2amino(*i));
	 exit(1);
       }
   }

 /* special symbols are expected to be used without an attempt to replace */ 
 
 return(j);
}


/* this function converts amino acid code into (generic) int */
int seq2dig(char *s)
{
 char amino[4];

 strncpy(amino,s,3);
 amino[3]='\0';
 
 /* the symbols as given below define in fact the generic alphabet
    used implictly throughout the code */

 /* to change !!! - this bunch of if's is not very elegant and one
    should define an array of symbols to go through it when desired */

 if(strncmp("ALA",amino,3)==0)return(0);
 if(strncmp("ARG",amino,3)==0)return(1);
 if(strncmp("ASN",amino,3)==0)return(2);
 if(strncmp("ASP",amino,3)==0)return(3);
 if(strncmp("CYS",amino,3)==0)return(4);
 if(strncmp("GLN",amino,3)==0)return(5);
 if(strncmp("GLU",amino,3)==0)return(6);
 if(strncmp("GLY",amino,3)==0)return(7);
 if(strncmp("HIS",amino,3)==0)return(8);
 if(strncmp("ILE",amino,3)==0)return(9);
 if(strncmp("LEU",amino,3)==0)return(10);
 if(strncmp("LYS",amino,3)==0)return(11);
 if(strncmp("MET",amino,3)==0)return(12);
 if(strncmp("PHE",amino,3)==0)return(13);
 if(strncmp("PRO",amino,3)==0)return(14);
 if(strncmp("SER",amino,3)==0)return(15);
 if(strncmp("THR",amino,3)==0)return(16);
 if(strncmp("TRP",amino,3)==0)return(17);
 if(strncmp("TYR",amino,3)==0)return(18);
 if(strncmp("VAL",amino,3)==0)return(19);
 /* now non-standard symbols */
 if(strncmp("GAP",amino,3)==0)return(20);
 if(strncmp("HYD",amino,3)==0)return(21);
 if(strncmp("POL",amino,3)==0)return(22);
 if(strncmp("CHG",amino,3)==0)return(23);
 if(strncmp("CRG",amino,3)==0)return(23); /* two names */
 if(strncmp("CH-",amino,3)==0)return(24);
 if(strncmp("CHN",amino,3)==0)return(24); /* two names */
 if(strncmp("ALL",amino,3)==0)return(25);
 if(strncmp("INS",amino,3)==0)return(26);
 if(strncmp("DEL",amino,3)==0)return(27);
 /* even more fancy symbols start after first 30 symbols */
 if(strncmp("CST",amino,3)==0)return(30); /* CYS in bridge */
 if(strncmp("HST",amino,3)==0)return(31); /* HIS in bridge */
 if(strncmp("USR",amino,3)==0)return(40); /* user defined */
 /* you may surely add more here, but add it in dig2amino as well */
 /* MAXS is currently used to allocate arrays for generic alph. */
 if(MAXS<50) printf(" WARNING - MAXS small and alphabet may be wrong \n");
 
 return(-1);
}

/* this function converts one letter amino code to (generic) int */
int one2dig(char *s)
{
 char amino[2];

 strncpy(amino,s,1);
 amino[1]='\0';

 if(strncmp("A",amino,1)==0)return(0);
 if(strncmp("R",amino,1)==0)return(1);
 if(strncmp("N",amino,1)==0)return(2);
 if(strncmp("D",amino,1)==0)return(3);
 if(strncmp("C",amino,1)==0)return(4);
 if(strncmp("Q",amino,1)==0)return(5);
 if(strncmp("E",amino,1)==0)return(6);
 if(strncmp("G",amino,1)==0)return(7);
 if(strncmp("H",amino,1)==0)return(8);
 if(strncmp("I",amino,1)==0)return(9);
 if(strncmp("L",amino,1)==0)return(10);
 if(strncmp("K",amino,1)==0)return(11);
 if(strncmp("M",amino,1)==0)return(12);
 if(strncmp("F",amino,1)==0)return(13);
 if(strncmp("P",amino,1)==0)return(14);
 if(strncmp("S",amino,1)==0)return(15);
 if(strncmp("T",amino,1)==0)return(16);
 if(strncmp("W",amino,1)==0)return(17);
 if(strncmp("Y",amino,1)==0)return(18);
 if(strncmp("V",amino,1)==0)return(19);

 if(strncmp("a",amino,1)==0)return(0);
 if(strncmp("r",amino,1)==0)return(1);
 if(strncmp("n",amino,1)==0)return(2);
 if(strncmp("d",amino,1)==0)return(3);
 if(strncmp("c",amino,1)==0)return(4);
 if(strncmp("q",amino,1)==0)return(5);
 if(strncmp("e",amino,1)==0)return(6);
 if(strncmp("g",amino,1)==0)return(7);
 if(strncmp("h",amino,1)==0)return(8);
 if(strncmp("i",amino,1)==0)return(9);
 if(strncmp("l",amino,1)==0)return(10);
 if(strncmp("k",amino,1)==0)return(11);
 if(strncmp("m",amino,1)==0)return(12);
 if(strncmp("f",amino,1)==0)return(13);
 if(strncmp("p",amino,1)==0)return(14);
 if(strncmp("s",amino,1)==0)return(15);
 if(strncmp("t",amino,1)==0)return(16);
 if(strncmp("w",amino,1)==0)return(17);
 if(strncmp("y",amino,1)==0)return(18);
 if(strncmp("v",amino,1)==0)return(19);

 return(-1);
}

/* this function converts back int into one-lett. amino acid code  */
char *dig2one(int i)
{ 
 /* assuming implictly generic alphabet */
 if(i==0)return("A");
 if(i==1)return("R");
 if(i==2)return("N");
 if(i==3)return("D");
 if(i==4)return("C");
 if(i==5)return("Q");
 if(i==6)return("E");
 if(i==7)return("G");
 if(i==8)return("H");
 if(i==9)return("I");
 if(i==10)return("L");
 if(i==11)return("K");
 if(i==12)return("M");
 if(i==13)return("F");
 if(i==14)return("P");
 if(i==15)return("S");
 if(i==16)return("T");
 if(i==17)return("W");
 if(i==18)return("Y");
 if(i==19)return("V");
 /* now non-standard symbols */
 /* I had to invent something ... */
 if(i==20)return("-"); /* GAP */
 if(i==21)return("B"); /* hydrophoBic */
 if(i==22)return("O"); /* pOlar */
 if(i==23)return("U"); /* plUs-minUs */
 if(i==24)return("X"); /* negative character */
 if(i==25)return("*"); /* ALL */
 if(i==26)return("+"); /* insertion */
 if(i==27)return("-"); /* deletion */
 if(i==30)return("J"); /* sulphur bridge (joint CYSs) */
 if(i==31)return("Z"); /* attach ions with some Z e.g. Z=2 for Fe */
 if(i==31)return("="); /* user defined */
 /* MAXS is currently used to allocate arrays for generic alph. */
 return("?");
}

/* this function converts back int into amino acid code  */
char *dig2amino(int i)
{ 
 /* assuming implictly generic alphabet */
 if(i==0)return("ALA");
 if(i==1)return("ARG");
 if(i==2)return("ASN");
 if(i==3)return("ASP");
 if(i==4)return("CYS");
 if(i==5)return("GLN");
 if(i==6)return("GLU");
 if(i==7)return("GLY");
 if(i==8)return("HIS");
 if(i==9)return("ILE");
 if(i==10)return("LEU");
 if(i==11)return("LYS");
 if(i==12)return("MET");
 if(i==13)return("PHE");
 if(i==14)return("PRO");
 if(i==15)return("SER");
 if(i==16)return("THR");
 if(i==17)return("TRP");
 if(i==18)return("TYR");
 if(i==19)return("VAL");
 /* now non-standard symbols */
 if(i==20)return("GAP");
 if(i==21)return("HYD");
 if(i==22)return("POL");
 if(i==23)return("CHG");
 if(i==24)return("CH-");
 if(i==25)return("ALL");
 if(i==26)return("INS");
 if(i==27)return("DEL");
 if(i==30)return("CST");
 if(i==31)return("HST");
 if(i==40)return("USR");
 /* MAXS is currently used to allocate arrays for generic alph. */
 return("ERR");
}

/* check if the query sequence contains low complexity regions */
float low_complexity(k_opt *op, a_file *F, k_info *Info, int *Seq, int *l_seq)
{
  int i,j,k,add,curr,n_four,n_eight,n_twelve,n_all,*in_lowCompl;
  float low_compl;
  
  
  in_lowCompl=malloc((*l_seq)*sizeof(int));
  if(in_lowCompl==NULL) alloc_err("local_pointer");   
  for(i=0;i<*l_seq;++i) *(in_lowCompl+i)=0;
  low_compl=0.0;

  /* now slide three windows (of length 4,8 and 12) through the sequence */
  for(i=12;i<*l_seq;++i)
    {
      curr=*(Seq+i);
      n_four=1;
      n_eight=0;
      n_twelve=0;
      for(j=i-1;j>i-4;--j) 
	{
	  add=1;
	  for(k=j+1;k<=i;++k) 
	    {
	      if(*(Seq+k)==*(Seq+j)) 
		{
		  add=0;
		  break;
		}
	    }
	  if(add) ++n_four;
	}
      n_eight+=n_four;
      for(j=i-4;j>i-8;--j) 
	{
	  add=1;
	  for(k=j+1;k<=i;++k) 
	    {
	      if(*(Seq+k)==*(Seq+j)) 
		{
		  add=0;
		  break;
		}
	    }
	  if(add) ++n_eight;
	}
      n_twelve+=n_eight;
      for(j=i-8;j>i-12;--j) 
	{
	  add=1;
	  for(k=j+1;k<=i;++k) 
	    {
	      if(*(Seq+k)==*(Seq+j)) 
		{
		  add=0;
		  break;
		}
	    }
	  if(add) ++n_twelve;
	}
      /*
      printf(" ires=%d  n_4=%d  n_8=%d  n_12=%d \n",i,n_four,n_eight,n_twelve);
      */
      if(n_twelve<4) 
	{
	  for(j=i;j>i-12;--j) *(in_lowCompl+j)=1; 
	}
      else
	{
	  if(n_eight<3) 
	    {
	      for(j=i;j>i-8;--j) *(in_lowCompl+j)=1;
	    }
	  else
	    {
	      if(n_four<2) 
		{
		  for(j=i;j>i-4;--j) *(in_lowCompl+j)=1;
		}
	    }
	}
      
    }
  
  n_all=0;
  for(i=0;i<*l_seq;++i)
    {
      if(*(in_lowCompl+i)==1) ++n_all;
    }  
  /* return percent of the residues in low complexity regions */
  low_compl=(100.0*n_all)/(float)(*l_seq);
  /*
  if(low_compl>2.0)
    printf(" %f percent in low complexity regions \n",low_compl); */
  
  free(in_lowCompl);
  return(low_compl);
  
}

/* check_for_membrane computes average hydrophobicity index */
int check_for_membrane(k_opt *op, a_file *F, k_info *Info, int *Seq, int *l_seq)
{
  int likely_membrane;
  float h_ind,h_th_likely,h_th_very_likely,h_th_sure;
  char stmp[MAXS];
    
  likely_membrane=0;
  /* set the thresholds */
  h_th_likely=0.4;
  h_th_very_likely=0.6; /* this will do the job for the time being */
  h_th_sure=0.8;
  
  h_ind=get_hydrophob_ind(op,Info,Seq,l_seq);
  if( h_ind > h_th_likely )
    {
      strcpy(stmp,"may be ");
      if( h_ind > h_th_very_likely ) strcpy(stmp,"is likely ");
      if( h_ind > h_th_sure ) strcpy(stmp,"is very likely ");
      
      sprintf(Info->Sout,"\n WARNING: average hydrophobicity= %5.2f \n",h_ind);
      to_stnd_out(Info->Sout);
      to_file(F->best,Info->Sout);
      sprintf(Info->Sout,
	      " Sequence %s %sof a MEMBRANE protein \n",op->qname,stmp);
      to_stnd_out(Info->Sout);
      to_file(F->best,Info->Sout);
      likely_membrane=1;
    }
  
  return(likely_membrane);
  
}

float get_hydrophob_ind(k_opt *op, k_info *Info, int *Seq, int *l_seq)
{
  int i;
  float h_ind;
  /* hydropathy index of J. Kyte and R.F. Doolittle */
  float hydrophob_ind[20]={ 1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,
			    3.8,-3.9,1.9,2.8,-1.6,-0.8,-0.7,-0.9,-1.3,4.2 };
  
  /* assumption: sequence is given in the generic alphabet */
  h_ind=0.0;
  for(i=0;i<*l_seq;++i)
    {  
      h_ind+=hydrophob_ind[*(Seq+i)];
    }
  
  return(h_ind/(float)(*l_seq));
  
}


