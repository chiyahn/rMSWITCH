#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rowSumC(NumericMatrix mat) 
{
	int nrow = mat.nrow();
	int ncol = mat.ncol();
	NumericVector vec = NumericVector(nrow);
	for (int i = 0; i < nrow; i++)
	  for (int j = 0; j < ncol; j++)
		vec[i] += mat(i,j);
	return vec;
}