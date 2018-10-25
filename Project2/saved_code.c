/*	fgets(string,81,pdb_file);
	line_num=0;
	while (string)
	{
		fgets(string,81,pdb_file);
		if (strncmp(string,"SEQRES",6)==0)
		{
			printf("sequence residue\n");
		}
		if (strncmp(string,"ATOM",4)==0)
		{
			printf("atom\n");
		}
		++line_num;
	} */

/*	for (i=0; i<2; ++i)
	{
		if (i==0)
		{
			fprintf(a_file,"query sequence:\t");
		}
		if (i==1)
		{
			fprintf(a_file,"datab sequence:\t");
		}
		for (j=0; j<l_alignment; ++j)
		{
			fprintf(a_file,"%c",integer2aachar(alignment[i][j],acids));
		}
		fprintf(a_file,"\n");
	}
	fprintf(a_file,"\n\n"); */