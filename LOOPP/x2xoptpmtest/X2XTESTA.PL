# the following program is to test the x2x (structure-to-structure) alignment module in the program loopp
# by aligning certain proteins structures to a database of protein structures and analyzing the output files
# for different (beta,d0) parameter values
# the program was started on 1/23/2001
# written by andrew smith, graduate student at cornell university
#!/Perl/bin/perl -w

# START OF PROGRAM

# set i/o file names
$f_seqhl = "SEQ_hl";
$f_seqte = "SEQ";
$f_xyzhl = "XYZ_hl";
$f_xyzte = "XYZ";
$f_best = "best.log";
$f_data = ">data.log";

# set options
$query = "1mba"; # the query protein structure to align into database of protein structures
$local = "-l";	 # local or global x2x alignments
$cgporngp = 0;	 # if cgporngp = 1, then (beta,d0) points will be searched with varying constant gap penalties 
				 # if cgporngp = 0, then (beta,d0) points will be searched with environment dependent gappen
$save = 10;		 # the number of saved parameters from optimize.pl
$Tzs = 5.0;		 # z-score threshold for determining families using loopp s2s alignments
$n2prnt = 20;    # the number of matches to print
$n_best = 20;	 # number of best matches to print	

# expand single query into a set a query structures to loopp over
@query = ("1mba","1ccr","4cms");
# store $save best (beta,d0) and gap penalties from optimize.pl
@beta = (6,6,6,5,2,3,8,6,10,5);
@d0 = (0.2,0.2,0.2,0.1,0.1,0.1,0.3,0.2,0.4,0.2);
@gp = (0.1,0.3,0.2,0.1,0.1);
# if environment dependent gap penalties (ngp where n=number of contacts), then two params gpara and gparb
if ($cgporngp == 0)
{
#	$gpa = $gpb = 0.0;
#	@gp = (@gpar0,@gpar1,@gpar2,@gpar3,@gpar4) = ([$gpa,$gpb],[$gpa,$gpb],[$gpa,$gpb],[$gpa,$gpb],[$gpa,$gpb]);
	@gp = (@gpar0,@gpar1,@gpar2,@gpar3,@gpar4) = ([0.03,0.3],[0.02,0.3],[0.01,0.3],[0.01,0.1],[0.01,0.2],[0.01,0.1],[0.05,0.1],[0.04,0.3],[0.01,0.1],[0.02,0.1]);
}
# debugging
#$m = 0;
#$k = 0;
#$f_seqhl = "SEQ_globins";
#$f_xyzhl = "XYZ_globins";

# open any output files to write to
open(DATA,$f_data) or die "can not open $f_data for writing: $!\n";
# loop over (beta,d0) pairs to generate output files
for ($m=0; $m<3; ++$m)
{
for ($k=0; $k<$save; ++$k)
{
	# clear previous arrays before next parameter loop
	@xnames = ();
	@x_ener = ();
	@x_rmsd = ();
	@xfmnames = ();
	@xfm_ener = ();
	@xfm_rmsd = ();
	@nxfmnames = ();
	@nxfm_ener = ();
	@nxfm_rmsd = ();
	@snames = ();
	@s_ener = ();
	@zscore = ();
	@family = ();

	if ($cgporngp == 1)
	{
		# call loopp with beta, d0, and constant gap penalty parameters
		@command_line = "loopp.exe -x -seq $f_seqhl -xyz $f_xyzhl -i 2 -fps -rsc $beta[$k] -fd0 $d0[$k] -g $gp[$k] $gp[$k] -pick $query[$m] -nprn $n2prnt -best $n_best -tr 50 $local\n";
		print @command_line;
		system(@command_line);
	}
	if ($cgporngp == 0)
	{
		# call loopp with beta, do, and environment dependent gap penalty parameters
		@command_line = "loopp.exe -x -seq $f_seqhl -xyz $f_xyzhl -i 2 -fps -rsc $beta[$k] -fd0 $d0[$k] -g 0.3 0.3 -gpar $gp[$k][0] $gp[$k][1] -pick $query[$m] -nprn $n2prnt -best $n_best -tr 50 $local\n";
		print @command_line;
		system(@command_line);
	}

	# read best.log file to get data from loopp x2x alignment
	# store matching protein names in @xnames, their energies in @x_ener, and their rmsds in @x_rmsd
	open(BESTLOG,$f_best) or die "can not open file $f_best for reading: $!\n";
	$n = 0;
	while ($line = <BESTLOG>)
	{
		if ($line =~ /ene=/)
		{
			($value1,$value2,$value3,$value4,$value5,$rest) = split(" ",$line);
			$xnames[$n] = $value1;
			$x_ener[$n] = $value3;
			$x_rmsd[$n] = $value5;
			$n++;
		}
	}
#	prnt_array(\@xnames,"xnames");
#	prnt_array(\@x_ener,"x_ener");
#	prnt_array(\@x_rmsd,"x_rmsd");

	# calculate sum of the rmsds
	$sumrms = 0.0;
	for ($i=0; $i<@x_rmsd; ++$i)
	{
		$sumrms += @x_rmsd[$i];
	}

	print "\n";
#	print "sumrms = $sumrms\n";
	print "\n";
	
	# get % seq identity of best matched proteins to query of x2x alignment
	($qprot,$match,$seqid)= &get_seqid;

#	prnt_array(\@$qprot,"qprot");
#	prnt_array(\@$match,"match");
#	prnt_array(\@$seqid,"seqid");
	
	# get family members determined by loopp z-score > 5.0 for local seq-seq alignments
	@command_line2 = "loopp.exe -s $f_seqhl -seq $f_seqhl -xyz $f_xyzhl -i 2 -pick $query[$m] -nprn $n2prnt -best $n_best $local\n";
	print @command_line2;
	system(@command_line2);
	print "\n";

	# read best.log file to get data from loopp s2s alignment
	# store matching protein names in @snames, their energies in @s_ener, and their z scores in @zscore
	open(BESTLOG,$f_best) or die "can not open file $f_best for reading: $!\n";
	$n2 = 0;
	while ($line = <BESTLOG>)
	{
		if ($line =~ /ene=/)
		{
			($value1,$value2,$value3,$value4,$value5,$rest) = split(" ",$line);
			$snames[$n2] = $value1;
			$s_ener[$n2] = $value3;
			$zscore[$n2] = $value5;
			$n2++;
		}
	}

#	prnt_array(\@snames,"snames");
#	prnt_array(\@s_ener,"s_ener");
#	prnt_array(\@zscore,"zscore");

#	print "n = $n\n";
#	print "n2 = $n2\n";

	# get % seq identity of best matched proteins to query of s2s alignment
	($sqprot,$smatch,$sseqid)= &get_seqid;
	
	# find protein matches with z scores > Tzs, this is definition of family member
	# assuming z-scores are already order from highest->lowest as index goes from 0->n
	for ($i=0; $i<@zscore; ++$i)
	{
		if ($zscore[$i] >= $Tzs)
		{
			$family[$i] = $snames[$i];
		}
	}
	
#	prnt_array(\@family,"family");
	
	# determine how many of the xnames for x2x are defined family members from the database and save those
	# names and their data (x2x energies and rmsds) to other arrays
	# also save names and data of non-family members
	$n_xfm = 0;
	$n_nxfm = 0;
	for ($i=0; $i<$n; ++$i)  # original: $n = @xnames, but n diff for diff queries and parameters
	{
		$pick = 0;
		for ($j=0; $j<@family; ++$j)
		{
			if ($xnames[$i] eq $family[$j])
			{	
				$xfmnames[$n_xfm] = $xnames[$i];
				$xfm_ener[$n_xfm] = $x_ener[$i];
				$xfm_rmsd[$n_xfm] = $x_rmsd[$i];
				++$n_xfm;
				$pick = 1;
				last;
			}
		}
		if ($pick == 0)	
		{
			$nxfmnames[$n_nxfm] = $xnames[$i];
			$nxfm_ener[$n_nxfm] = $x_ener[$i];
			$nxfm_rmsd[$n_nxfm] = $x_rmsd[$i];
			++$n_nxfm;
		}
	}
	$fracrec = $n_xfm/@family;
	
#	print "n_xfm = $n_xfm\n";
#	print "n_nxfm = $n_nxfm\n";
#	print "fraction of family members recognized = $fracrec\n";
#	prnt_array(\@xfmnames,"xfmnames");
#	prnt_array(\@xfm_ener,"xfm_ener");
#	prnt_array(\@xfm_rmsd,"xfm_rmsd");
#	prnt_array(\@nxfmnames,"nxfmnames");
#	prnt_array(\@nxfm_ener,"nxfm_ener");
#	prnt_array(\@nxfm_rmsd,"nxfm_rmsd");
	
	# determine highest rmsd for the family members of the x2x alignment (xfmnames and xfm_rmsd)
	$pos = 0;
	$max = $xfm_rmsd[$pos];
	for ($i=0; $i<$n_xfm; ++$i)
	{
		if ($xfm_rmsd[$i] > $max)
		{
			$max = $xfm_rmsd[$i];
			$pos = $i;
		}
	}
	$highrmsd = $max;	
#	print "maximum x2x family member rmsd = $max\n";
	
	# determine the number of false positives
	# find maximum energy of all x2x family members and then determine how many of the non-family members have
	# energies lower than that maximum value
	$pos = 0;
	$max = $xfm_ener[$pos];
	for ($i=0; $i<@xfm_ener; ++$i)
	{
		if ($xfm_ener[$i] > $max)
		{
			$max = $xfm_ener[$i];
			$pos = $i;
		}
	}
#	print "max energy = $max\n";
	$n_fp = 0;
	for ($i=0; $i<@nxfm_ener; ++$i)
	{
		if ($nxfm_ener[$i] < $max)
		{
			++$n_fp;
		}
	}
#	print "n_fp = $n_fp\n";
	
	if ($cgporngp == 0)
	{
		print DATA "run: beta = $beta[$k], d0 = $d0[$k], gpara = $gp[$k][0], gparb = $gp[$k][1]\n";
	}
	if ($cgporngp == 1)
	{
		print DATA "run: beta = $beta, d0 = $d0, cgp = $gp[$k]\n";
	}
	print DATA "query = $query[$m] and database = $f_seqhl and $f_xyzhl\n";
	print DATA "command line for x2x = @command_line\n";
	print DATA "command line for s2s = @command_line2\n";
	print DATA "xnames\tx_ener\tx_rmsd\tfamily\tnon-family\tmatch\t\t%seqid\tsmatch\t%sseqid\n";
	for ($i=0; $i<$n; ++$i)
	{
		print DATA "$xnames[$i]\t\t$x_ener[$i]\t$x_rmsd[$i]\t\t$family[$i]\t\t$nxfmnames[$i]\t\t$$match[$i]\t\t$$seqid[$i]\t\t$$smatch[$i]\t\t$$sseqid[$i]\n";
	}
	print DATA "\n";
	print DATA "number of family members pulled out of database = $n_xfm\n";
	print DATA "fraction of family members recognized = $fracrec\n";
	print DATA "number of false positives = $n_fp\n";
	print DATA "highest family rmsd = $highrmsd\n";
	print DATA "\n\n";	 	  
}
}

close BESTLOG;
close ALIGN;
close DATA;

die "exiting x2xtest program\n";

#END OF PROGRAM

#=================================================================================================================
# &prnt_array prints an array to stdout
sub prnt_array
{
	my ($array,$name) = @_;
	my ($i,$n_elem);

	print "$name = ";
	for ($i=0; $i<@$array; ++$i)
	{
		print "$$array[$i]";
		print " ";
	}
	print "\n";
	$n_elem = @$array;
	print "number of elements in array = $n_elem\n";
	print "\n";
#	print "daliz = ";
#	for ($i=0; $i<$n; ++$i)
#	{
#		print "$$daliz[$i]";
#		print " ";
#	}
#	print "\n";
}

#=================================================================================================================
# &read_famset reads database file and scop family file(s) finds all matching proteins and puts them into an array
sub read_famset
{
	my ($dbseqfile,$dbxyzfile,$famfile) = @_;
	my ($line);
	
	open (DBNAME,$dbseqfile) or die "can not open $dbseqfile for reading: $!\n";
	open (FAMILY,$famfile) or die "can not open $famfile for reading: $!\n";
	while ($line = <$dbseqfile>)
	{
		if ($line !~ /##/ || $line !~ /chain/)
		{
		}
	}
}

#=================================================================================================================
# &get_seqid reads the protein match and % sequence identity from the file "alignments.log" and puts them into the 
# arrays xname and seqid
sub get_seqid
{
	my ($line,$f_align,@query,@match,@seqid,$term1,$term2,$term3,$term4,$term5,$term6,@term,$n,$n2);

	$f_align = "alignments.log";
	open(ALIGN,$f_align) or die "can not open $f_align for reading seqid: $!\n";
	
	$n = 0;
	while ($line = <ALIGN>)
	{
		if ($line =~ /query=/)
		{
			($term1,$term2,$term3,$term4,$term5,$term6) = split(" ",$line);
#			print "term2 = $term2\n";
#			print "term4 = $term4\n";
			$query[$n] = $term2;
			$match[$n] = $term4;
			next;
		}
		if ($line =~ /ident=/)
		{
			@term = split("=",$line);
#			prnt_array(\@term,"term");
			for ($n2=0; $n2<@term; ++$n2)
			{
				if ($term[$n2] =~ /.ident/)
				{
#					print "term = $term[$n2]\n";
					chop $term[++$n2];
					chop $term[$n2];
					$seqid[$n] = $term[$n2];
					last;
				}
			}
			++$n;					
		}	
	}
#	print "n = $n\n";
#	prnt_array(\@query,"query");
#	prnt_array(\@match,"match");
#	prnt_array(\@seqid,"seqid");

	return(\@query,\@match,\@seqid);
}


