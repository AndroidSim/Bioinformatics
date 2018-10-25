/*	i=1;
	while((eof=fscanf(seq_file,"%c",&sequence[i-1]))!=EOF)
	{
/*		printf("int of sequence = %d\n",(int)sequence[i-1]); */
/*		if(sequence[i-1]=='\n') /* iscntrl(sequence[i-1])) */
/*		{
			sequence=realloc(sequence,--i);
/*			printf("seq char = %c",sequence[i-1]);
			printf("i = %d\n",i); */
/*		}
		if (isdigit(sequence[i-1]))
		{
			sequence=realloc(sequence,--i);
			l_seq=i;
			++n_seq;
			fscanf(seq_file,".\n");
			break; 
		}
		++i;
		sequence=realloc(sequence,i);
	}
	if (eof==EOF)
	{
		return(0);
	} */
/*		if(ch=='\n') 
		{
			sequence=realloc(sequence,--i);
		} */
	/* open scoring matrix file and read into memory the number of aas */
/*	strcpy(file->f_sm_name,"Blosum50.txt");
	file->f_sm=fopen(file->f_sm_name,"r");
	if (file->f_sm==NULL)
	{
		printf("can not open file %s\n",file->f_sm_name);
		return(0);
	} */
	/* open sequence file if number of sequences equals zero 
	if (n_seq==0)
	{
		strcpy(file->f_seq_name,"sequences.txt");
		file->f_seq=fopen(file->f_seq_name,"r");
		if (file->f_seq==NULL)
		{
			printf("can not open file %s\n",file->f_seq_name);
			return(0);
		}
	} */
