# this program is for the optimization of the parameters, d0 and beta, used in stru2stru module of loopp
# program started on 12/8/00
# written by Andrew Smith
#!/Perl/bin/perl -w

#START OF PROGRAM
# initialize major boolean parameters
$true = 1;
$false = 0;

# initialize other options
$local = 'l'; # local or global x2x alignments
$save = 5; # number of the top (beta,d0) coordinates to save with best correlation coefficients
$gpornot = 1; # if gapsornot = 1, then the best (beta,d0) points will be searched with varying
			  # gap penalties and if gapsornot =0, then above will not happen

# set step size and maximum/minimum values for optimization parameters d0 and beta (and gap search if selected)
# also set query protein 
$stepbeta = 1.0;		# original: 0.5
$stepd0 = 0.1;			# original: 0.1		
$minbeta = 1.0;			# original: 0.1
$maxbeta = 20.0;		# original: 50
$mind0 = 0.1;			# original: 0.1
$maxd0 = 2.0;			# original: 10
$stepgp = 0.1;
$mingp = 0.1;
$maxgp = 0.5;
$query = "1mba";	

# input/output sequence and structure files 
$f_seq = "SEQ_globins";
$f_xyz = "XYZ_globins";
$f_best = "best.log";

# initialize output arrays
#@sprs = (0.0,0.0,0.0,0.0,0.0,);
#@oneoverkt = (0.0,0.0,0.0,0.0,0.0);
#@deenot = (0.0,0.0,0.0,0.0,0.0);
#@sumlrms0 = (0.0,0.0,0.0,0.0,0.0);

# debugging
#$maxbeta = $minbeta;
#$maxd0 = $mind0;

$start = (times)[0];
($ustart,$sstart,$custart,$csstart) = times;

$f_scores = ">scores.log";
open(SCORES,$f_scores) or die "can not open file $f_scores for writing: $!\n";

print "\n";
print "starting search\n\n";

$n_points = 0;
$beta = 10;
$d0 = 0.3;
#for ($beta=$minbeta; $beta<=$maxbeta; $beta+=$stepbeta)
#{
#LOOP2:	for ($d0=$mind0; $d0<=$maxd0; $d0+=$stepd0)
#	{
		if ($n_points != 0)
		{
			print "searching next point: \n";
		}
		print "beta = $beta\n";
		print "d0 = $d0\n";
		# call loopp program with do and kt as command line arguments
		@command_line = "loopp.exe -x -seq $f_seq -xyz $f_xyz -fps -rsc $beta -fd0 $d0 -g 0.3 0.3 -pick $query -nprn 16 -best 16 -tr 50 -l";
#		print @command_line;
		print "loopp.exe -x -seq $f_seq -xyz $f_xyz -fps -rsc $beta -fd0 $d0 -g 0.3 0.3 -pick $query -nprn 16 -best 16 -tr 50 -l\n";
		system("loopp.exe -x -seq $f_seq -xyz $f_xyz -fps -rsc $beta -fd0 $d0 -g 0.3 0.3 -pick $query -nprn 16 -best 16 -tr 50 -$local");
		
		# read best.log file to get scores from loopp x2x alignment
		# store loopp scores in @loopps
		$index = 0;
		open(BESTLOG,$f_best) or die "can not open file $f_best for reading: $!\n";
		$line = <BESTLOG>;
#		print "line1 = $line\n";
		$line = <BESTLOG>;
#		print "line2 = $line\n";
		$line = <BESTLOG>;
#		print "line3 = $line\n";
		$line = <BESTLOG>;
#		print "line4 = $line\n";
		$line = <BESTLOG>;
#		print "line5 = $line\n";
		while ($line = <BESTLOG>)
		{
			if ($line =~ /ene=/)
			{
				($value1,$value2,$value3,$value4,$value5,$rest) = split(" ",$line);
				$loopps[$index] = $value3;
				$loopp_rms[$index] = $value5;
				$index++;
			}
		}

		# if loopp determines that with the given parameters of beta and d0 some alignments do not meet the
		# criteria for a good match or alignment (such as number of aligned residues,low energy, low rms
		# such that not all of the best matched proteins are printed to best.log then skip this (beta,d0)
		# coordinate
		if ($index < 16)
		{
			next LOOP2;
		}

#		print @loopps;
#		print "number of data values for loopps = $index\n";
#		print "loopp_rms = ";
#		prnt_array(\@loopp_rms);

		# store DALI Z_scores and rmsds in an array
		@daliz = (30.9,18.6,18.5,18.1,18.0,17.8,17.6,16.7,16.6,15.7,15.5,15.5,14.8,14.3,13.3,12.2);
		@dali_rms = (0.0,1.9,1.8,1.9,2.0,2.0,2.2,2.0,1.9,2.3,2.3,2.5,2.4,2.4,2.8,2.7);		

		printf SCORES "run: beta = $beta, d0 = $d0\n";
		printf SCORES "command line = @command_line\n";
		printf SCORES "dali_scores\t\t\tloopp_scores\n";
		for ($k=0; $k<$index; ++$k)
		{
			printf SCORES "$daliz[$k]\t\t\t\t$loopps[$k]\n";
		}
		printf SCORES "\n";

		$sumlrms = 0.0;
		for ($k=0; $k<@loopp_rms; ++$k)
		{
			$sumlrms += @loopp_rms[$k];
		}
#		$variable = @loopp_rms;
#		print "$variable\n";
#		print "@loopp_rms\n";
#		print "sum of loopp rms = $sumlrms\n";

		# optimization of the parameters involves saving the ($beta,$d0) coordinate pair the has the largest 
		# negative spearman rank-order correlation coefficient and/or lowest sum of rmsds
		($d,$zd,$probd,$rs,$probrs) = (0.0,0.0,0.0,0.0,0.0);
		($d,$zd,$probd,$rs,$probrs) = spear(\@loopps,\@daliz,$index,$d,$zd,$probd,$rs,$probrs);
		
#		printf SCORES "loopps = ";
#		pf_array(SCORES,\@loopps);
#		printf SCORES "\n";
#		printf SCORES "daliz = ";
#		pf_array(SCORES,\@daliz);
		printf SCORES "\n";
		printf SCORES "spearman rank-order correlation coefficient = $rs\n";
		printf SCORES "sum of the rmsds for loopp scores = $sumlrms\n\n";

#		print "\n";
#		print "rs = $rs\n";
#		print "d = $d\n";
#		print "zd = $zd\n";
#		print "probd = $probd\n";
#		print "probrs = $probrs\n";
		print "\n";

#		if ($n_points == 0)
#		{
#			$rs0 = $rs;
#			$sprs[$n_points] = $rs;
#			$beta0 = $beta;
#			$1overkt[$n_points] = $beta;
#			$d00 = $d0;
#			$deenot[$n_points] = $d0;
#			$sumlrms0[$n_points] = $sumlrms;
#		}
#		else 
#		{
#			if ($rs < $rs0 && $sumlrms < $sumlrms0)
#			{
#				# same as for n_points == 0			
#			}
#		}
		if ($n_points < $save)
		{
			$sprs[$n_points] = $rs;
			$oneoverkt[$n_points] = $beta;
			$deenot[$n_points] = $d0;
			$sumlrms0[$n_points] = $sumlrms;
		}
		else
		{
			$pos = 0;
			$max = $sprs[$pos];
			for ($k=1; $k<$save; ++$k)
			{
				if ($sprs[$k] > $max)
				{
					$max = $sprs[$k];
					$pos = $k;
				}
			}
			if ($rs < $sprs[$pos] && $sumlrms < $sumlrms0[$pos])
			{
				$sprs[$pos] = $rs;
				$oneoverkt[$pos] = $beta;
				$deenot[$pos] = $d0;
				$sumlrms0[$pos] = $sumlrms;
			}
		}
		++$n_points;

#		print "rs0 = $rs0\n";
#		print "beta0 = $beta0\n";
#		print "d00 = $d00\n";
#		print "\n";
#	}
#}

close BESTLOG;
close SCORES;

# debugging or test run
@oneoverkt = (3,7,5,4,7);
@deenot = (0.1,0.2,0.2,0.1,0.3);

if ($gpornot == 1)
{
	$f_gpscores = ">gpscores.log";
	open(GPSCORES,$f_gpscores) or die "can not open file $f_gpscores for writing: $!\n";

	$n_gppoints = 0;
	# searching different gap penalties to find optimal one for best (beta,d0) coordinate pair
	print "searching gap penalty range\n\n";
	for ($k=0; $k<$save; ++$k)
	{
		$n_gp = 0;
LOOP3:	for ($gappen=$mingp; $gappen<=$maxgp; $gappen+=$stepgp)
		{
			# call loopp with specified gappens
			print "loopp.exe -x -seq $f_seq -xyz $f_xyz -fps -rsc $oneoverkt[$k] -fd0 $deenot[$k] -g $gappen $gappen -pick $query -nprn 16 -best 16 -tr 50 -l\n";
			system("loopp.exe -x -seq $f_seq -xyz $f_xyz -fps -rsc $oneoverkt[$k] -fd0 $deenot[$k] -g $gappen $gappen -pick $query -nprn 16 -best 16 -tr 50 -$local");

			# read best.log file to get scores from loopp x2x alignment
			# store loopp scores in @loopps
			$index = 0;
			open(BESTLOG,$f_best) or die "can not open file $f_best for reading: $!\n";
			$line = <BESTLOG>;
#			print "line1 = $line\n";
			$line = <BESTLOG>;
#			print "line2 = $line\n";
			$line = <BESTLOG>;
#			print "line3 = $line\n";
			$line = <BESTLOG>;
#			print "line4 = $line\n";
			$line = <BESTLOG>;
#			print "line5 = $line\n";
			while ($line = <BESTLOG>)
			{
				if ($line =~ /ene=/)
				{
					($value1,$value2,$value3,$value4,$value5,$rest) = split(" ",$line);
					$loopps[$index] = $value3;
					$loopp_rms[$index] = $value5;
					$index++;
				}
			}

			# if loopp determines that with the given parameters of beta and d0 some alignments do not meet the
			# criteria for a good match or alignment (such as number of aligned residues,low energy, low rms
			# such that not all of the best matched proteins are printed to best.log then skip this (beta,d0)
			# coordinate
			if ($index < 16)
			{
				next LOOP3;
			}

			$sumrms = 0.0;
			for ($m=0; $m<@loopp_rms; ++$m)
			{
				$sumrms += @loopp_rms[$m];
			}

			# store DALI Z_scores and rmsds in an array
			@daliz = (30.9,18.6,18.5,18.1,18.0,17.8,17.6,16.7,16.6,15.7,15.5,15.5,14.8,14.3,13.3,12.2);
			@dali_rms = (0.0,1.9,1.8,1.9,2.0,2.0,2.2,2.0,1.9,2.3,2.3,2.5,2.4,2.4,2.8,2.7);
			
			# optimization of the parameters involves saving the ($beta,$d0) coordinate pair the has the largest 
			# spearman rank-order correlation coefficient
			($dgp,$zdgp,$probdgp,$rsgp,$probrsgp) = (0.0,0.0,0.0,0.0,0.0);
			($dgp,$zdgp,$probdgp,$rsgp,$probrsgp) = spear(\@loopps,\@daliz,$index,$dgp,$zdgp,$probdgp,$rsgp,$probrsgp);

#			$sprmn[$n_gp] = $rsgp;
#			$sumrms0[$n_gp] = $sumrms;
			
			if ($n_gp == 0)
			{
				printf GPSCORES "gap penalty run of ($oneoverkt[$k],$deenot[$k])\n";
				printf GPSCORES "gappen\t\t\tspearman_rocc\t\t\tsum_rms\n";
			} 
			printf GPSCORES "$gappen\t\t\t\t$sprmn[$n_gp]\t\t$sumgprms0[$n_gp]\n";
			printf GPSCORES "\n";

			# optimize gappen by finding the gappen with the maximum negative spearman r-o cc and/or lowest sum of
			# rmsds
			if ($n_gp == 0)
			{
				$gp[$k] = $gappen;
				$gpsprmn = $sprmn;
				$sumgprms = $sumrms;
			}
			else 
			{
				if ($sprmn < $gpsprmn && $sumrms < $sumgprms)
				{
					$gp[$k] = $gappen;
					$gpsprmn = $sprmn;
					$sumgprms = $sumrms;		
				}
			}
#			$gpbeta[$k] = $oneoverkt[$k];
#			$gpd0[$k] = $deenot[$k];

			++$n_gp;
			++$n_gppoints;
		}	
	}
}


$end = (times)[0];
($uend,$send,$cuend,$csend) = times;

print "the search is finished\n";
print "number points searched = $n_points\n";
# print "the optimal parameters are: beta = $beta0 and d0 = $d00\n";	
# print "with the maximum spearman rank-order correlation coefficient = $rs0\n";
print "the top $save best parameters are: \n\n";
for ($k=0; $k<$save; ++$k)
{
#	$varb = $oneoverkt[$k];
#	$vard = $deenot[k];
#	$vars = $sprs[$k];
#	print "beta = $varb and d0 = $vard with rs = $vars\n";
	print "beta = $oneoverkt[$k] and d0 = $deenot[$k] with rs = $sprs[$k] and sumrms = $sumlrms0[$k]\n";
}
print "\n";

if (gpornot == 1)
{
	print "number of points searched for best gap penalty per (beta,d0) = $n_gp\n";
	print "the optimal gap penalty for the best $save parameters are:\n\n";
	for $k (0..($save-1))
	{
		print "for beta = $oneoverkt[$k] and d0 = $deenot[$k], best gappen = $gp[$k]\n";
	}
}

($difftime,$udifftime,$sdifftime) = (($end-$start),($uend-$ustart),($send-$sstart));
print "program took $difftime CPU seconds\n";
print "program took $udifftime user seconds and $sdifftime system seconds\n";

die "exiting optimize program $!\n";
# END OF PROGRAM
				
# ============================================================================================
# the following 6 or so subroutines are taken from "Numerical Recipes in C: The Art of Scientific Computing"
# and rewritten in Perl script

# &spear calculates the Spearman Rank-Order Correlation Coefficient between @loopps and @daliz
sub spear 
{
	my ($loopps,$daliz,$n,$d,$zd,$probd,$rs,$probrs) = @_;
	my ($i,$j,$vard,$t,$sg,$sf,$fac,$en3n,$en,$df,$aved,@wksp1,@wksp2);

	# debugging 
#	print "n = $n\n";
#	print "loopps = ";
#	prnt_array(\@$loopps);
#	print "\n";
#	print "daliz = ";
#	prnt_array(\@$daliz);
#	print "\n"; 

	for ($j=0; $j<$n; $j++)
	{
		$wksp1[$j] = $$loopps[$j];
		$wksp2[$j] = $$daliz[$j];
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
#	print "wksp1 = ";
#	prnt_array(\@wksp1);
#	print "\n";
#	print "wksp2 = ";
#	prnt_array(\@wksp2);
#	print "\n";

	$sf = crank($n,\@wksp1,$sf);
#	print "after ranking wksp1: \n";
#	print "wksp1 = ";
#	prnt_array(\@wksp1);
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
#	print "wksp1 = ";
#	prnt_array(\@wksp1);
#	print "\n";
#	print "wksp2 = ";
#	prnt_array(\@wksp2);
#	print "\n";

	$sg = crank($n,\@wksp2,$sg);
#	print "after ranking wksp2: \n";
#	print "wksp2 = ";
#	prnt_array(\@wksp2);
#	print "sg = $sg\n";
#	print "\n";

	printf SCORES "ranked daliz = ";
	pf_array(SCORES,\@wksp2);
	printf SCORES "\n";
	printf SCORES "ranked loopps w/r ranked daliz = ";
	pf_array(SCORES,\@wksp1);
	
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
	prnt_array(\@$ra);
	print "array rb in sort2 = ";
	prnt_array(\@$rb);

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

# =================================================================================================
#end of numerical recipes subroutines
# &prnt_array prints an array to stdout
sub prnt_array
{
	my ($array) = @_;
	my ($i,$n_elem);

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

#=====================================================================================================
# &prntf_array prints an array to a file in a column
sub prntf_array
{
	my ($filehandle,$array,$beglineoffset) = @_;
	my ($i,$filepos,$filepos2,$filepos3);

	$filepos = tell $filehandle;
	for ($i=0; $i<@$array; ++$i)
	{
		printf $filehandle "$$array[$i]\n";
		$filepos2 = tell $filehandle;
		$filepos3 = $filepos2+($filepos-$beglineoffset);
		seek($filehandle,$filepos3,0);
	}
}

# =====================================================================================================
# &pf_array prints an array to a file horizontally
sub pf_array
{
	my ($filehandle,$array) = @_;
	my ($i);

	for ($i=0; $i<@$array; ++$i)
	{
		printf $filehandle "$$array[$i]";
		printf $filehandle " ";
	}
	printf $filehandle "\n";
} 

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

# =================================================================================================
# &kendl1 calculates the nonparametric correlation statistic called kendall's tau.
# this subroutine was taken from Numerical Recipes in C.
sub kendl1
{
	my ($data1,$data2,$n,$tau,$z,$prob) = @_;
	my ($n1,$n2,$k,$j,$is,$svar,$aa,$a1,$a1);

	$n1 = 0;
	$n2 = 0;
	$is = 0;

	for ($j=0; $j<$n-1; $j++)
	{
		for ($k=($j+1); $k<$n; $k++)
		{
			$a1 = $$data1[$j]-$$data1[$k];
			$a2 = $$data2[$j]-$$data2[$k];
			$aa = $a1*$a2;
			if ($aa)
			{
				++$n1;
				++$n2;
				$aa > 0.0 ? ++$is : --$is;
			}
			else
			{
				if ($a1) 
				{
					++$n1;
				}
				if ($a2)
				{
					++$n2;
				}
			}
		}
	}
	$tau = $is/(sqrt($n1)*sqrt($n2));
	$svar = (4.0*$n+10.0)/(9.0*$n*($n-1.0));
	$z = ($tau)/sqrt($svar);
	$prob = erfcc(abs($z)/1.4142136);
}
		
		
		 




