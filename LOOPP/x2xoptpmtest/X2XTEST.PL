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
$save = 3;		 # the number of saved parameters from optimize.pl
$Tzs = 5.0;		 # z-score threshold for determining families using loopp s2s alignments
$n2prnt = 20;    # the number of matches to print
$n_best = 20;	 # number of best matches to print	

# expand single query into a set a query structures to loopp over
@query = ("1mba","1ccr","4cms");
# store $save best (beta,d0) and gap penalties from optimize.pl
@beta = (3,3,3);
@d0 = (0.1,0.1,0.1);
@gp = (0.1,0.3,0.2,0.1,0.1);
@dn0 = (1,2,3,4,5,6,7,8,9,10);
@lambda = (0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0);
# if environment dependent gap penalties (ngp where n=number of contacts), then two params gpara and gparb
if ($cgporngp == 0)
{
#	$gpa = $gpb = 0.0;
#	@gp = (@gpar0,@gpar1,@gpar2,@gpar3,@gpar4) = ([$gpa,$gpb],[$gpa,$gpb],[$gpa,$gpb],[$gpa,$gpb],[$gpa,$gpb]);
#	@gp = (@gpar0,@gpar1,@gpar2,@gpar3,@gpar4) = ([0.04,0.1],[0.05,0.15],[0.05,0.1],[0.035,0.1],[0.045,0.15],[0.04,0.05],[0.045,0.1],[0.04,0.15],[0.015,0.15],[0.045,0.05]);
	@gp = (@gpar0,@gpar1,@gpar2,@gpar3,@gpar4) = ([0.05,0.1],[0.045,0.05],[0.01,0.1]);
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
for ($k=0; $k<=$save; ++$k)
{
	for ($r=0; $r<=9; ++$r)
	{
	for ($s=0; $s<=10; ++$s)
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
#		@command_line = "loopp.exe -x -seq $f_seqhl -xyz $f_xyzhl -i 2 -fps -rsc $beta[$k] -fd0 $d0[$k] -g 0.3 0.3 -gpar $gp[$k][0] $gp[$k][1] -pick $query[$m] -nprn $n2prnt -best $n_best -tr 50 $local\n";
#		print @command_line;
#		system(@command_line);
		# call loopp with (beta,d0,ngpa,ngpb,dn0,lambda) for combined metric
		@command_line = "loopp.exe -x -seq $f_seqhl -xyz $f_xyzhl -i 2 -fps -ncx2x $dn0[$r] $lambda[$s] -rsc $beta[$k] -fd0 $d0[$k] -g 0.3 0.3 -gpar $gp[$k][0] $gp[$k][1] -pick $query[$m] -nprn $n2prnt -best $n_best -tr 50 $local\n";
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

	# calculate the spearman rank-order correlation coefficient between the structure scores (x_ener) and their 
	# rmsds
	($d,$zd,$probd,$rs,$probrs) = (0.0,0.0,0.0,0.0,0.0);
	($d,$zd,$probd,$rs,$probrs) = spear(\@x_ener,\@x_rmsd,$n,$d,$zd,$probd,$rs,$probrs);

	print "\n";
	print "for spearman rank-order correlation coefficient:\n";
	print "rs = $rs\n";
	print "d = $d\n";
	print "zd = $zd\n";
	print "probd = $probd\n";
	print "probrs = $probrs\n";
	print "\n"; 

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
		print DATA "run: beta = $beta[$k], d0 = $d0[$k], gpara = $gp[$k][0], gparb = $gp[$k][1], dn0 = $dn0[$r], lambda = $lambda[$s]\n";
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
	print DATA "rs of x_ener and x_rmsd = $rs\n";
	print DATA "\n\n";
	}
	}	 	  
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

# ============================================================================================
# the following 6 or so subroutines are taken from "Numerical Recipes in C: The Art of Scientific Computing"
# and rewritten in Perl script

# &spear calculates the Spearman Rank-Order Correlation Coefficient between @loopps and @daliz
sub spear 
{
	my ($array1,$array2,$n,$d,$zd,$probd,$rs,$probrs) = @_;
	my ($i,$j,$vard,$t,$sg,$sf,$fac,$en3n,$en,$df,$aved,@wksp1,@wksp2);

	# debugging 
#	print "n = $n\n";
#	prnt_array(\@$array1,"array1");
#	print "\n";
#	prnt_array(\@$array2,"array2");
#	print "\n"; 

	for ($j=0; $j<$n; $j++)
	{
		$wksp1[$j] = $$array1[$j];
		$wksp2[$j] = $$array2[$j];
	}

	# i wrote the following block
	if ($n < 50)
	{
		piksr2($n,\@wksp1,\@wksp2);
	}
	else
	{
		sort2($n,\@wksp1,\@wksp2);
	}

	# debugging
#	prnt_array(\@wksp1,"wksp1");
#	print "\n";
#	prnt_array(\@wksp2,"wksp2");
#	print "\n";

	$sf = crank($n,\@wksp1,$sf);
#	print "after ranking wksp1: \n";
#	prnt_array(\@wksp1,"wksp1");
#	print "sf = $sf\n";
#	print "\n";
	
	# i wrote the following block
	if ($n < 50)
	{
		piksr2($n,\@wksp2,\@wksp1);
	}
	else
	{
		sort2($n,\@wksp2,\@wksp1);
	}

	# debugging
#	prnt_array(\@wksp1,"wksp1");
#	print "\n";
#	prnt_array(\@wksp2,"wksp2");
#	print "\n";

	$sg = crank($n,\@wksp2,$sg);
#	print "after ranking wksp2: \n";
#	prnt_array(\@wksp2,"wksp2");
#	print "sg = $sg\n";
#	print "\n";

#	pf_array(SCORES,\@wksp2,"ranked daliz");
#	printf SCORES "\n";
#	pf_array(SCORES,\@wksp1,"ranked loopps w/r ranked daliz");
	
	$d = 0.0;
	for ($j=0; $j<$n; $j++)
	{
		$d += ($wksp1[$j]-$wksp2[$j])**2;
	}

	$rs = spearman(\@wksp1,\@wksp2);
#	print "rs from spearman = $rs\n";

#	print "d = $d\n";
	$en = $n;
	$en3n = $en*$en*$en-$en;
#	print "en3n = $en3n\n";
	$aved = $en3n/6.0-($sf+$sg)/12.0;
#	print "aved = $aved\n";
	$fac = (1.0-$sf/$en3n)*(1.0-$sg/$en3n);
#	print "fac = $fac\n";
	$vard = (($en-1.0)*$en*$en*($en+1.0)**2/36.0)*$fac;
#	print "vard = $vard\n";
	$zd = ($d-$aved)/sqrt $vard;
	$probd = erfcc(abs $zd/1.4142136);
#	$rs = (1.0-(6.0/$en3n)*($d+0.5*($sf+$sg)))/$fac;
#	print "rs = $rs\n";
#	$t = ($rs)*sqrt(($en-2.0)/(($rs+1.0)*(1.0-$rs)));
	$df = $en-2.0;
	$probrs = betai(0.5*$df,0.5,$df/($df+$t*$t));

#	print "within spear fxn the values are: \n";
#	print "rs = $rs\n";
#	print "d = $d\n";
#	print "zd = $zd\n";
#	print "probd = $probd\n";
#	print "probrs = $probrs\n";
#	print "\n";

	return($d,$zd,$probd,$rs,$probrs);
}

# ============================================================================================
# &crank replaces the elements of a sorted array by their ranks, including midranking of ties.
# the script is slightly modified from the num. rec. to take into account the indices of the 
# array starts with zero and the ranked array should not contain zeros
sub crank
{
	my ($n,$w,$s) = @_;
	my ($j,$ji,$jt,$t,$rank);

	$s = 0.0;
	$j = 0;
	while ($j<$n-1) # original: while (j<n)
	{
		if ($$w[$j+1] != $$w[$j])
		{
			$$w[$j]=$j+1; # original: w[j]=j, but rank=j+1
			++$j;
		}
		else
		{
	LINE:	for ($jt=$j+1; $jt<$n; $jt++) #original: for (jt=j+1; jt<=n; jt++)
			{
				if ($$w[$jt] != $$w[$j])
				{
					last LINE;
				}
			}
			$rank = 0.5*($j+$jt+1); # original: rank=0.5(j+jt-1),now j=j+1,jt=jt+1, so rank=0.5(j+jt+1)
			for ($ji=$j; $ji<$jt; $ji++)  # original: for ($ji=$j; $ji<=($jt-1); $ji++)
			{
				$$w[$ji] = $rank;
			}
			$t = $jt-$j;
			$s += $t*$t*$t-$t;
			$j = $jt;
		}
	}
	if ($j==$n-1) # original: if (j==n)
	{
		$$w[$j] = $n; # original: w[n]=n
	}

	return $s;
}

# ==============================================================================================
# &sort2 sorts an array into ascending numerical order using the Heapsort algorithm, while making
# the corresponding rearrangement of another array
sub sort2
{
	my ($n,$ra,$rb) = @_;
	my ($l,$j,$ir,$i,$rra,$rrb);

	print "in sort2 subroutine, n = $n\n";
	print "array ra in sort2 = ";
	prnt_array(\@$ra,"ra");
	print "array rb in sort2 = ";
	prnt_array(\@$rb,"rb");

#	$l = ($n>>1)+1;
	$l = ($n>>1);
#	$l = ($n-1)>>1;
#	$ir = $n;
	$ir = $n-1;
	for (;;)
	{
		if ($l>0) # original: if ($l>1)
		{
			$rra = $$ra[--$l];
			$rrb = $$rb[$l];
		}
		else
		{
			$rra = $$ra[$ir];
			$rrb = $$rb[$ir];
			$$ra[$ir] = $$ra[0]; # original: ra[ir]=ra[1]
			$$rb[$ir] = $$rb[0]; # original: rb[ir]=rb[1]
			if (--$ir == 0) # original: if (--ir == 1)
			{
				$$ra[0] = $rra; # original: ra[1]=rra
				$$rb[0] = $rrb; # original: rb[1]=rrb
				return;
			}
		}
		$i = $l;
		$j = ($l<<1)+1; # original: j=(l<<1)
		while ($j <= $ir)
		{
			if ($j < $ir && $$ra[$j] < $$ra[$j+1])
			{
				++$j;
			}
			if ($rra < $$ra[$j])
			{
				$$ra[$i] = $$ra[$j];
				$$rb[$i] = $$rb[$j];
				$j += ($i = $j);
			}
			else
			{
				$j = $ir+1;
			}
		}
		$$ra[$i] = $rra;
		$$rb[$i] = $rrb;
	}
}

# ================================================================================================
# &piksr2 is an alternative sort routine similar to sort2 but is O(N**2) whereas sort2 = O(Nlog base 2 N)
# therefore this is only practical for n<50
sub piksr2
{
	my ($n,$arr,$brr) = @_;
	my ($i,$j,$a,$b);

	for ($j=1; $j<$n; $j++)
	{
		$a = $$arr[$j];
		$b = $$brr[$j];
		$i = $j-1;
		while ($$arr[$i] > $a)
		{
			$$arr[$i+1] = $$arr[$i];
			$$brr[$i+1] = $$brr[$i];
			$i--;
		}
		$$arr[$i+1] = $a;
		$$brr[$i+1] = $b;
	}
}
# ================================================================================================
# &betai returns the incomplete beta function
sub betai
{
	my ($a,$b,$x) = @_;
	my ($bt);

	if ($x < 0.0 || $x > 1.0)
	{
		die "bad x in subroutine betai: $!\n"
	}
	if ($x == 0.0 || $x == 1.0)
	{
		$bt = 0.0;
	}
	else
	{
		$bt = exp(gammln($a+$b)-gammln($a)-gammln($b)+$a*log($x)+$b*log(1.0-$x));
	}
	if ($x < ($a+1.0)/($a+$b+2.0))
	{
		return $bt*betacf($a,$b,$x)/$a;
	}
	else
	{
		return 1.0-$bt*betacf($b,$a,1.0-$x)/$b;
	}
}

# =================================================================================================
# &betacf is for continued fraction for incomplete beta function, used by betai
sub betacf
{
	my ($a,$b,$x) = @_;
	my ($qap,$qam,$qab,$em,$tem,$d,$bz,$bm,$bp,$bpp,$az,$am,$ap,$app,$aold,$m);

	$bm = 1.0;
	$az = 1.0;
	$am = 1.0;
	$ITMAX = 100;
	$EPS = 0.0000003;

	$qab = $a+$b;
	$qap = $a+1.0;
	$qam = $a-1.0;
	$bz = 1.0-$qab*$x/$qap;
	for ($m=1; $m<=$ITMAX; $m++)
	{
		$em = $m;
		$tem = $em+$em;
		$d = $em*($b-$em)*$x/(($qam+$tem)*($a+$tem));
		$ap = $az+$d*$am;
		$bp = $bz+$d*$bm;
		$d = -($a+$em)*($qab+$em)*$x/(($qap+$tem)*($a+$tem));
		$app = $ap+$d*$az;
		$bpp = $bp+$d*$bz;
		$aold = $az;
		$am = $ap/$bpp;
		$bm = $bp/$bpp;
		$az = $app/$bpp;
		$bz = 1.0;
		if (abs($az-$aold) < ($EPS*(abs $az)))
		{
			return $az;
		}
	}
	print "a or b too big, or ITMAX too small in betacf\n";
}

# ================================================================================================
# &gammln returns the value ln(gamma(xx)) for xx>0.  full accuracy is obtained for xx>1.
# for 0<xx<1, the reflection formula can be used first
sub gammln
{
	my ($xx) = @_;
	my ($x,$tmp,$ser,$j);
	my @cof = (76.18009173,-86.50532033,24.01409822,-1.231739516,0.00120858003,-0.00000536382);

	$x = $xx-1.0;
	$tmp = $x+5.5;
	$tmp -= ($x+0.5)*log $tmp;
	$ser = 1.0;
	for ($j=0; $j<=5; $j++)
	{
		$x += 1.0;
		$ser += $cof[$j]/$x;
	}
	return -tmp+log(2.50662827465*$ser);
}

# ===============================================================================================
# &erfcc returns the complementary error function erfc(x) with fractional error everywhere < 1.2x10-7
sub erfcc
{
	my $x = @_;
	my ($t,$z,$ans);

	$z = abs $x;
	$t = 1.0/(1.0+0.5*$z);
	$ans = $t*exp(-$z*$z-1.26551223+$t*(1.00002368+$t*(0.37409196+$t*(0.09678418+$t*(-0.18628806+$t*(0.27886807+$t*(-1.13520398+$t*(1.48851587+$t*(-0.82215223+$t*0.17087277)))))))));
	return ($x >= 0.0 ? $ans : 2.0-$ans);
}

#=====================================================================================================
#=====================================================================================================
# &spearman is another way to calculate spearman rank-order correlation coefficient
sub spearman
{
	my ($data1,$data2) = @_;
	my ($i,$n,$sum1,$avg1,$sum2,$avg2,$sqrootsum1,$var1,$var2,$sd1,$sd2,$covar,$rocc);
	
	if (@$data1 != @$data2)
	{
		die "length_data1 != length_data2: length_data1 = @$data1, length_data2 = @$data2: $!\n"
	}
	else
	{
		$n = @$data1;  # or =@$data2, does not matter
	}
	
#	print "within spearman: \n";
#	print "loopps = ";
#	prnt_array(\@$data1);
#	print "daliz = ";
#	prnt_array(\@$data2);  

	# calculate average of rank data 1
	$sum1 = 0.0;
	for ($i=0; $i<$n; ++$i)
	{
		$sum1 += $$data1[$i];
	}
	$avg1 = $sum1/$n;
#	print "avg of loopps rank = $avg1\n";

	# calculate average of rank data 2
	$sum2 = 0.0;
	for ($i=0; $i<$n; ++$i)
	{
		$sum2 += $$data2[$i];
	}
	$avg2 = $sum2/$n;
#	print "avg of daliz rank = $avg2\n";

	# calculate linear correlation coefficient between array1 and array2
	for ($i=0; $i<$n; ++$i)
	{
		$var1 += ($$data1[$i]-$avg1)**2;
		$var2 += ($$data2[$i]-$avg2)**2;
	}
	$sd1 = sqrt($var1);
	$sd2 = sqrt($var2);
#	print "sd loopps rank = $sd1\n";
#	print "sd daliz rank = $sd2\n";
	$covar = 0.0;
	for ($i=0; $i<$n; ++$i)
	{
		$covar += ($$data1[$i]-$avg1)*($$data2[$i]-$avg2);
	}
#	print "covar = $covar\n";
	$rocc = $covar/($sd1*$sd2);

	return $rocc;
}

