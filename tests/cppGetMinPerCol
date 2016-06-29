#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
SEXP GetMinPerCol (Rcpp::NumericMatrix matrix_rcpp)
{
	arma::mat    matrix(matrix_rcpp.begin(),
								matrix_rcpp.nrow(), matrix_rcpp.ncol(), false);

	int n = matrix.n_rows;
	int m = matrix.n_cols;
	arma::colvec min_col(m);

	for (int j = 0; j < m; j++)
	{
		double min_value = std::numeric_limits<double>::infinity();
		for (int i = 0; i < n; i++)
		{
			if (min_value > matrix.at(i,j))
				{
					min_value = matrix.at(i,j);
					min_col(j) = min_value;
				}
		}
	}

	return wrap(min_col);
}
