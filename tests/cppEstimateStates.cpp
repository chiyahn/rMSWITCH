// has been already included in cppFilterIndep.cpp.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SQRT2PI = 2.50662827463; // sqrt(2*pi)

SEXP EtaIndep  (Rcpp::NumericVector y_rcpp,
				Rcpp::NumericMatrix y_lagged_rcpp,
				Rcpp::NumericMatrix z_dependent_rcpp,
				Rcpp::NumericMatrix z_independent_rcpp,
				Rcpp::NumericVector beta_rcpp,
				Rcpp::NumericVector mu_rcpp,
				Rcpp::NumericVector sigma_rcpp,
				Rcpp::NumericMatrix gamma_dependent_rcpp,
				Rcpp::NumericVector gamma_independent_rcpp)
{
	int M = mu_rcpp.size();
	int n = y_rcpp.size();
	arma::mat eta(n, M);

	arma::colvec y(y_rcpp.begin(), y_rcpp.size(), false);
	arma::mat    y_lagged(y_lagged_rcpp.begin(),
								y_lagged_rcpp.nrow(), y_lagged_rcpp.ncol(), false);
	arma::mat    z_dependent(z_dependent_rcpp.begin(),
								z_dependent_rcpp.nrow(), z_dependent_rcpp.ncol(), false);
	arma::mat    z_independent(z_independent_rcpp.begin(),
								z_independent_rcpp.nrow(), z_independent_rcpp.ncol(), false);
	arma::colvec beta(beta_rcpp.begin(), beta_rcpp.size(), false);
	arma::colvec mu(mu_rcpp.begin(), mu_rcpp.size(), false);
	arma::colvec sigma(sigma_rcpp.begin(), sigma_rcpp.size(), false);
	arma::mat    gamma_dependent(gamma_dependent_rcpp.begin(),
								gamma_dependent_rcpp.nrow(), gamma_dependent_rcpp.ncol(), false);
	arma::colvec gamma_independent(gamma_independent_rcpp.begin(),
								gamma_independent_rcpp.size(), false);


	for (int j = 0; j < M; j++)
	{
	  eta.col(j) = y - y_lagged * beta - z_dependent * gamma_dependent.col(j)
									- z_independent * gamma_independent - mu(j);
		eta.col(j) = eta.col(j) % eta.col(j); // element-wise multiplication
		eta.col(j) = exp(-eta.col(j) / (2 * (sigma(j) * sigma(j))));
		eta.col(j) = eta.col(j) / (SQRT2PI * sigma(j));
	}

	return (wrap(eta));
}

// [[Rcpp::export]]
SEXP EstimateStates  (Rcpp::NumericVector y_rcpp,
				Rcpp::NumericMatrix y_lagged_rcpp,
				Rcpp::NumericMatrix z_dependent_rcpp,
				Rcpp::NumericMatrix z_independent_rcpp,
				Rcpp::NumericVector beta_rcpp,
				Rcpp::NumericVector mu_rcpp,
				Rcpp::NumericVector sigma_rcpp,
				Rcpp::NumericMatrix gamma_dependent_rcpp,
				Rcpp::NumericVector gamma_independent_rcpp)
{
	SEXP eta_rcpp = EtaIndep(y_rcpp, y_lagged_rcpp,
						z_dependent_rcpp, z_independent_rcpp,
						beta_rcpp, mu_rcpp, sigma_rcpp,
						gamma_dependent_rcpp, gamma_independent_rcpp);
	arma::mat eta = Rcpp::as<arma::mat>(eta_rcpp);
	int n = eta.n_rows;
	int M = eta.n_cols;

	arma::colvec states(n);

	for (int i = 0; i < n; i++)
	{
		double best = -std::numeric_limits<double>::infinity();
		int best_index = -1;
		for (int j = 0; j < M; j++)
			if (eta(i,j) > best)
			{
				best = eta(i,j);
				best_index = j;
			}
		states(i) = best_index; // CAUTION: the first state is numerized as zero.
	}

	return (wrap(states));
}
