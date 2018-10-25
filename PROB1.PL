#!/Perl/bin/perl -w

$answer = 0;
$answer = N3(5,5,5);
print "answer = $answer\n";

#================================================
sub N3
{
	my ($n,$m,$l) = @_;
	my ($ans);

#	$ans = 0;
	if ($n==0 && $m==0 && $l==0)
	{
		$ans = 1;
		return $ans;
	}
	if ($n==0 && $m==0 && $l!=0)
	{
		$ans = 1;
		return $ans;
	}
	if ($n==0 && $l==0 && $m!=0)
	{
		$ans = 1;
		return $ans;
	}
	if ($m==0 && $l==0 && $n!=0)
	{
		$ans = 1;
		return $ans;
	}
	if ($n==0 && $m!=0 && $l!=0)
	{
		$ans = N2($m,$l);
		return $ans;
	}
	if ($m==0 && $l!=0 && $n!=0)
	{
		$ans = N2($n,$l);
		return $ans;
	}
	if ($l==0 && $m!=0 && $n!=0)
	{
		$ans = N2($m,$n);
		return $ans;
	}
	if ($n!=0 && $m!=0 && $l!=0)
	{
		$ans = N3($n-1,$m-1,$l-1)+N3($n,$m,$l-1)+N3($n,$m-1,$l)+N3($n-1,$m,$l)+N3($n-1,$m-1,$l)+N3($n-1,$m,$l-1)+N3($n,$m-1,$l-1);
	}
	return $ans;
}

#===================================
sub N2
{
	my ($i,$j) = @_;
	my ($sol);

	if ($i==0 && $j==0)
	{
		$sol = 1;
		return $sol;
	}
	if ($i==0 && $j!=0)
	{
		$sol = 1;
		return $sol;
	}
	if ($i!=0 && $j==0)
	{
		$sol = 1;
		return $sol;
	}
	if ($i!=0 && $j!=0)
	{
		$sol = N2($i-1,$j-1)+N2($i,$j-1)+N2($i-1,$j);
	}
	return $sol;
}
	
