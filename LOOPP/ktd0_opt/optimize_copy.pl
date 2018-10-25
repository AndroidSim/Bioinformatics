# this program is for the optimization of the parameters, d0 and beta (formerly kt), used in stru2stru module of loopp
# program started on 12/8/00
# written by Andrew Smith
#!/Perl/bin/perl -w

# initialize major boolean parameters
$true = 1;
$false = 0;

# initialize other options
$local = 'l';

# set step size and maximum/minimum values for optimization parameters d0 and beta
# also set query protein 
$stepbeta = 0.5;
$stepd0 = 0.1;
$minbeta = 0.1;
$maxbeta = 50;
$mind0 = 0.1;
$maxd0 = 10;
$query = "1mba";

# input/output sequence and structure files 
$f_seq = "SEQ_globins";
$f_xyz = "XYZ_globins";
$f_best = "best.log";

# debugging
#$maxbeta = $minbeta;
#$maxd0 = $mind0;

$n_runs = 0;
$beta = 10;
$d0 = 0.5;
#for ($beta=$minbeta; $beta<=$maxbeta; $kt+=$stepbeta)
#{
#	for ($d0=$mind0; $d0<=$maxd0; $d0+=$stepd0)
#	{
		# call loopp program with do and kt as command line arguments
		print "loopp.exe -x -seq $f_seq -xyz $f_xyz -fps -rsc $beta -fd0 $d0 -g 0.5 0.5 -pick $query -nprn 16 -best 16 -tr 50";
		system("loopp.exe -x -seq $f_seq -xyz $f_xyz -fps -rsc $beta -fd0 $d0 -g 0.5 0.5 -pick $query -nprn 16 -best 16 -tr 50");
		
		# read best.log file to get scores from loopp x2x alignment
		# store loopp scores in @loopps
		$index = 0;
		open(BESTLOG,$f_best) or die "can't open file $f_best: $!\n";
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
				($value1,$value2,$value3,$rest) = split(" ",$line);
				$loopps[$index] = $value3;
				$index++;
			}
		}
#		print @loopps;
		print "number of data values for loopps = $index\n";

# optimization of the parameters involves saving the ($beta,$d0) coordinate pair the has the largest 
# spearman rank-order correlation coefficient
#store DALI Z_scores in an array
		@daliz = (30.9,18.6,18.5,18.1,18.0,17.8,17.6,16.7,16.6,15.7,15.5,15.5,14.8,14.3,13.3,12.2);		
		($d,$zd,$probd,$rs,$probrs) = (0.0,0.0,0.0,0.0,0.0);
		($d,$zd,$probd,$rs,$probrs) = spear(\@loopps,\@daliz,$index,$d,$zd,$probd,$rs,$probrs);

		print "rs = $rs\n";
		print "d = $d\n";
		print "zd = $zd\n";
		print "probd = $probd\n";
		print "probrs = $probrs\n";
		print "\n";

		if (n_runs == 0)
		{
			$rs0 = $rs;
			$beta0 = $beta;
			$d00 = $d0;
		}
		else 
		{
			if ($rs > $rs0)
			{
				$rs0 =$rs;
				$beta0 = $beta;
				$d00 = $d0;		
			}
		}
		++$n_runs;
#	}
#}

print "n_runs = $n_runs\n";
print "The search is finished\n";
print "The optimal parameter are: beta = $beta0 and d0 = $d0\n";
print "\n";
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
	print "n = $n\n";
	print "loopps = ";
	prnt_array(\@$loopps);
	print "\n";
	print "daliz = ";
	prnt_array(\@$daliz);
	print "\n"; 

	for ($j=0; $j<$n; $j++)
	{
		$wksp1[$j] = $$loopps[$j];
		$wksp2[$j] = $$daliz[$j];
	}

	if ($n < 50)
	{
		piksr2($n,\@wksp1,\@wksp2);
	}
	else
	{
		sort2($n,\@wksp1,\@wksp2);
	}

	# debugging
	print "wksp1 = ";
	prnt_array(\@wksp1);
	print "\n";
	print "wksp2 = ";
	prnt_array(\@wksp2);
	print "\n";

	$sf = crank($n,\@wksp1,$sf);
	print "after ranking wksp1: \n";
	print "wksp1 = ";
	prnt_array(\@wksp1);
	print "sf = $sf\n";
	print "\n";
	
	if ($n < 50)
	{
		piksr2($n,\@wksp2,\@wksp1);
	}
	else
	{
		sort2($n,\@wksp2,\@wksp1);
	}

	# debugging
	print "wksp1 = ";
	prnt_array(\@wksp1);
	print "\n";
	print "wksp2 = ";
	prnt_array(\@wksp2);
	print "\n";

	$sg = crank($n,\@wksp2,$sg);
	print "after ranking wksp2: \n";
	print "wksp2 = ";
	prnt_array(\@wksp2);
	print "sg = $sg\n";
	print "\n";

	$d = 0.0;
	for ($j=0; $j<$n; $j++)
	{
		$d += ($wksp1[$j]-$wksp2[$j])**2;
	}

	$en = $n;
	$en3n = $en*$en*$en-$en;
	$aved = $en3n/6.0-($sf+$sg)/12.0;
	$fac = (1.0-$sf/$en3n)*(1.0-$sg/$en3n);
	$vard = (($en-1.0)*$en*$en*($en+1.0)**2/36.0)*$fac;
	$zd = ($d-$aved)/sqrt $vard;
	$probd = erfcc(abs $zd/1.4142136);
	$rs = (1.0-(6.0/$en3n)*($d+0.5*($sf+$sg)))/$fac;
	$t = ($rs)*sqrt(($en-2.0)/(($rs+1.0)*(1.0-$rs)));
	$df = $en-2.0;
	$probrs = betai(0.5*$df,0.5,$df/($df+$t*$t));

	print "within spear fxn the values are: \n";
	print "rs = $rs\n";
	print "d = $d\n";
	print "zd = $zd\n";
	print "probd = $probd\n";
	print "probrs = $probrs\n";
	print "\n";

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
	LINE:	for ($jt=$j+1; $jt<$n; $jt++)
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
#	$ir = $n;
	$ir = $n-1;
	for (;;)
	{
		if ($l>1)
		{
			$rra = $$ra[--$l];
			$rrb = $$rb[$l];
		}
		else
		{
			$rra = $$ra[$ir];
			$rrb = $$rb[$ir];
			$$ra[$ir] = $$ra[0];
			$$rb[$ir] = $$rb[0];
			if (--$ir == 0)
			{
				$$ra[0] = $rra;
				$$rb[0] = $rrb;
				return;
			}
		}
		$i = $l;
		$j = $l << 1;
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





