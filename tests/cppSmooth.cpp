#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
SEXP Smooth (Rcpp::NumericMatrix xi_k_rcpp,
			Rcpp::NumericMatrix transition_probs_rcpp)
{
	int n = xi_k_rcpp.nrow();
	int M = xi_k_rcpp.ncol();

	arma::mat 	xi_n_t(M, n);
	arma::mat 	xi_k(xi_k_rcpp.begin(),
					 n, M, false);
	arma::mat 	transition_probs(transition_probs_rcpp.begin(),
								M, M, false);
	arma::mat 	xi_k_t = xi_k.t(); // produces a copy of a transpose

	xi_n_t.col(n-1) = xi_k_t.col(n-1);
	for (int k = (n-2); k >= 0; k--)
 		xi_n_t.col(k) = xi_k_t.col(k) %
 			(transition_probs * (xi_n_t.col(k+1) / (transition_probs * xi_k_t.col(k))));

	return (wrap(xi_n_t.t()));
}
