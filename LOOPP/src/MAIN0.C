/* Learning, Observing and Outputting Protein Patterns (LOOPP)    */
/*        by Jarek Meller and Ron Elber                           */
/* Jerusalem (Hebrew University) and Ithaca (Cornell University)  */
/*        1999/2000     v. 2.000                                  */
/*                                                                */
/*                   MAIN                                         */
 
#include "loopp.h"
   
/* 
   Unlike the previous version of LOOPP (called then Linear Optimization
   Of Protein Potentials), which was primarily supposed to be used to 
   develop protein folding potentials using linear programming,
   the present version of the program is meant to be used for FOLD
   RECOGNITION in the first place. Thus, one may use LOOPP to perform:
   a) standard sequence to sequence alignments
   b) structurally biased sequence to sequence alignments
   c) sequence to structure (THREADING) alignments
   d) structure to structure alignemnts.
   One may also use this program to generate non-redundant libraries of
   folds, both for the training and recognition. Our novel structure to
   structure alignment method is used for this purpose.

   Several models of (off-lattice, reduced representation) threading 
   potentials are implemented in the current version. These include:
   a) standard pairwise potentials 
   b) continuous pairwise models (of Lennard-Jones type) 
   c) profile Threading Onion Models, 
   as introduced and optimized (using LP training) by JM and RE 
   (see loopp_doc.html for details). If you don't like the potentials
   (or scoring functions) trained by us, you may still use LOOPP to design
   improved folding potentials of your own.

   In the training mode (options -q or -d) the LP constraints, resulting
   from the requirement of perfect recognition of native shapes within
   a database, are generated. In order to obtain LP constraints the 
   gapless threading is performed and the differencies between native 
   structures and the decoys are stored in files to be read by LP solvers.
   Energy differences between native and decoys are printed to a file
   xscan.log. Differences in contacts are put into a file current.dcon/mps 
   (to write differences in contacts while scaning use option -w/-mps).
    
   Sequences and structures of proteins are read from files 
   SEQ and XYZ and potential from current.pot , respectively 
   (if not specified otherwise). The default format of XYZ file is:
   name num_of_res 
   x_ca_i  y_ca_i  z_ca_i  x_sc_i  y_sc_i  z_sc_i 
   ...
   where sc stands for side chain coordinates. The default format of 
   SEQ file is:
   name num_of_res 
   sym_i
   ...
   where sym_i denotes identity of the subsequent residues in terms of
   an alphabet defined by the user. The definition of an alphabet is very 
   flexible and ranges from H/P two-letter alphabet to extended alphabets 
   with user defined symbols.

   Most of the options are taken from the standard input, although 
   the descriptions of the scoring function (threading potential) is 
   defined in the file Model. It is used to specify parameters, as 
   type of potential, alphabet, cutoff distance for non-continuous models.
   Check function set_options for Model format. Options specified in the
   Model file overwrite command line specifications. For most of the standard
   recognition functions, though, the default setup should be fine and
   the Model file does not need to be used. Use option -h to get small help.


   Please send comments and report bugs to: 
   meller@cs.cornell.edu
   meller@phys.uni.torun.pl
   http://www.cs.cornell.edu/home/meller/  

*/


/* main main main main main main */

main(int argc, char *argv[])
{
 int i;
 /* define structures and pointers to structures */
 k_opt opt,*op; 
 a_file Fs,*F; 
 solv_sh s_sh,*ssh; 
 k_dynpt dynpt_stru,*dpt;
 k_prot prot_stru,*prot; 
 k_Model Mod_stru,*model; 
 k_info Info_stru,*Info;
 k_mps *var_mps[MAX],*tail_mps[MAX];
 time_t begin,end;

 /* init pointers to structures */
 op=&opt; ssh=&s_sh; F=&Fs; Info=&Info_stru; 
 model=&Mod_stru; prot=&prot_stru; dpt=&dynpt_stru;  
 /* so each pointer points to an address where the corresponding
    structure is kept and then each array (pointer) that is a member 
    of this structure points to an address of the first element of
    this array as allocated by malloc in alloc_mem */
 
 begin=time(NULL);
 /* print small hint */
 if(argc<2) 
   {
     sprintf(Info->Sout," Perhaps try option -help first ... \n");
     to_stnd_out(Info->Sout);
   }
  
 /* set the default values of input parameters, interpret command line */
 set_options(argc,argv,F,op,Info);

 /* allocate memory for the arrays and other stuff */
 if(!alloc_mem(F,prot,op,Info,model,ssh,dpt))
   {
     sprintf(Info->Sout,"Allocation of memory failed ... \n"); 
     to_stnd_out(Info->Sout);
     exit(1); 
   }

 /* prepare all the necessary files - default or modified names from options */
 prepare_files(F,op,Info,model->Potn); 
 
 /* read Model description and get potential (scoring function) */
 set_Model(F,op,Info,model->Imod,model->Alph,model->Adr,model->Potn);

 /* print information about options and files used later on */
 if(op->info2) prn_setup(F,op,model->Alph);  

 /* init info and store_mps */
 init_info(op,F,Info,&var_mps[0],&tail_mps[0]); 

 /* compute contact maps and reference energies */
 build_irrep(op,F,ssh,prot,Info,model);

 /* skip the rest if only the query files are to be created */
 if(op->justq)
   {
     clean_up(F,prot,op,Info,model,ssh,dpt);
     exit(0);    
   }
 
 /* if align perform sequence to sequence alignment  */
 if(op->align) align_seq(op,F,ssh,prot,Info,model,dpt);
 
 /* compare (align) two structures using averaged THOM2 potential  */
 if(op->strucmp) stru2stru(op,F,ssh,prot,Info,model,dpt);
 
 /* if exam perform sequence to structure alignment (threading) */ 
 if(op->exam) thread_seq(op,F,ssh,prot,Info,model,dpt,&var_mps[0],&tail_mps[0]);

 /* in learning phase, in turn, use gapless threading of database seqeunces into
    database structures to compare native energies with those of decoys and to 
    get LP constraints in terms of contacts */
 if(op->learn) 
   get_ineq(op,F,ssh,prot->conR,prot->conT,prot->IJcont,model->Imod,
	    model->Potn,prot->R_ene,prot->seq,Info,prot->conRvdw,
	    prot->conTvdw,prot->IJr12,prot->IJr6,prot->native_IJrij,
	    model->Alph,model->Adr,&var_mps[0],&tail_mps[0]);

 /* if you have independently generated decoys use this function 
    to get inequalities */
 if(op->decoy) 
   get_ineq_dec(op,F,ssh,prot->conR,prot->conT,prot->IJcont,model->Imod,
		model->Potn,prot->R_ene,prot->seq,Info,prot->conRvdw,
		prot->conTvdw,model->Alph,model->Adr,&var_mps[0],&tail_mps[0]);

 /* get additional informations and contact statistics */
 if(op->info) info(F->info,Info,model->Potn,op,model->Imod,model->Alph,
		   model->Adr);

 /* write down the inequalities in MPS format */
 if(op->mps) write_mps(F,op,&Info->n_prn,&var_mps[0]);

 /* now one should nicely close all the files and free memory */ 
 if(op->ver) prn_signature(argc,argv,F,op,Info);
 clean_up(F,prot,op,Info,model,ssh,dpt);
 end=time(NULL);
 printf("loopp ran for %f seconds\n",difftime(end,begin));
 return(0); 
 
}  
/* end of main   end of main   end of main   end of main   end of main */

