// I'm just playing around on this to see if I can make the Halton Sequence in C++.
// Next phase is to build this into the Master Sample.
//

#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

//' Internal function for log base b, since it isn't obvious in C++.
// [[Rcpp::export(rng = false)]]
double log_b(double x, int & base)
{
	double tmp;
	tmp = log(x) / log(base);
	return tmp;
}

//' Draw Halton Sequence values for a single dimension.
//' Note that this was borrowed from the Internet and is not my implmementation.
//'
//'
//' @param x An integer for the starting index k >= 0.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param n Number of samples to draw.
//'
//' @examples
//' HaltonSeq(k = 0, base = 2, n = 10)
//'
// [[Rcpp::export(rng = false)]]
NumericVector HaltonSeq(const int & k, double & base, int & n) 
{
	NumericVector xk(n);
	int index = k;
	double f;
	for (int i = 0; i < n; i++) {
		index = k + i; 
		f = 1;
		while(index > 0){
			f = f/base;
			xk[i] = xk[i] + f*(fmod(index,base));
			index = index/base;
		}
	}
    return xk;
}

//' Fast Implementation of finding x and y order numbers to feed into the linear congruence equation.
//' This solves equation solves for a_i in the HIP paper.
//' This is essentially an internal function and shouldn't be worried about.
//'
//'
//' @param lxy A Matrix of lower x y coordinates of a Halton Box inside the unit box.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param J Integer of 2 values that represent the numbers 2^J1, 3^J2.
//'
//'
// [[Rcpp::export(rng = false)]]
NumericMatrix GetBoxIndices(NumericMatrix& lxy, IntegerVector& base, IntegerVector J)
{
int n = lxy.nrow();
NumericVector ai(2);
NumericMatrix results(n,2);

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			ai[j] = 0;
			for(int k = 1; k <= J[j]; k++)
			{
				ai[j] += (int(floor(lxy(i,j) * pow( base[j], k))) % base[j]) * pow( base[j], k - 1);
			}
		results(i,j) = ai[j];
		}
	}
	return results;
}

//' Solve system of linear congruence from HIP paper to order HIP boxes.
//' See page 5 of Robertson et al. 2018 Halton Iterative Partitioning.
//' This is essentially an internal function and shouldn't be worried about.
//'
//'
//' @param A Matrix that is in numeric for computational reasons but is the a_i solutions for all HIP boxes.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param J Integer of 2 values that represent the numbers 2^J1, 3^J2.
//'
//'
// [[Rcpp::export(rng = false)]]
NumericVector SolveCongruence(NumericMatrix& A, NumericVector& base, NumericVector J)
{
	NumericVector b = NumericVector::create( pow(base[0], J[0]), pow(base[1], J[1]) );
	double B = b[0]*b[1];

	double j, ll, jj, x;
	bool test;
	NumericVector Index(A.nrow());
	double tmp;
	for(int i = 0; i < A.nrow(); i++)
	{	
		j = 0; test = false;

		ll = (A(i,0) > A(i,1)) ? 0 : 1;		// for speed use the larger value to loop through.
		jj = (ll == 1) ? 0 : 1;				// Track jj as the one to test against.
		while(!test && (j < B))
		{
			x = A(i,ll) + (b[ll])*j;		// Find multiple
			tmp = floor(x / b[jj]);	// Does that particular multiple match with the other?
			test = (x - b[jj]*tmp) == A(i,jj);
			j ++;
		}
		if(test) Index(i) = x;			// Make sure there is an actual solution
		if(!test) Index(i) = -j;		// If no solution make it obvious.
	}
	return Index;
}