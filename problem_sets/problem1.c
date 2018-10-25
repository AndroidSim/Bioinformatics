
#include "stdio.h"
#include "math.h"

int N3(int n, int m, int l);
int N2(int i, int j);

void main (void)
{
	int answer;

	answer = 0;
	answer = N3(5,5,5);
	printf("answer = %d\n",answer);
}

/*================================================================================*/
int N3(int n, int m, int l)
{
	int ans;

	if (n==0 && m==0 && l==0)
	{
		ans = 1;
		return(ans);
	}
	if (n==0 && m==0 && l!=0)
	{
		ans = 1;
		return(ans);
	}
	if (n==0 && l==0 && m!=0)
	{
		ans = 1;
		return(ans);
	}
	if (m==0 && l==0 && n!=0)
	{
		ans = 1;
		return(ans);
	}
	if (n==0 && m!=0 && l!=0)
	{
		ans = N2(m,l);
		return(ans);
	}
	if (m==0 && l!=0 && n!=0)
	{
		ans = N2(n,l);
		return(ans);
	}
	if (l==0 && m!=0 && n!=0)
	{
		ans = N2(m,n);
		return(ans);
	}
	if (n!=0 && m!=0 && l!=0)
	{
		ans = N3(n-1,m-1,l-1)+N3(n,m,l-1)+N3(n,m-1,l)+N3(n-1,m,l)+N3(n-1,m-1,l)+N3(n-1,m,l-1)+N3(n,m-1,l-1);
	}
	return(ans);
}

/*==============================================================================*/
int N2(int i, int j)
{
	int sol;

	if (i==0 && j==0)
	{
		sol = 1;
		return(sol);
	}
	if (i==0 && j!=0)
	{
		sol = 1;
		return(sol);
	}
	if (i!=0 && j==0)
	{
		sol = 1;
		return(sol);
	}
	if (i!=0 && j!=0)
	{
		sol = N2(i-1,j-1)+N2(i,j-1)+N2(i-1,j);
	}
	return(sol);
}
