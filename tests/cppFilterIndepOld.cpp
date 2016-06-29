#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SQRT2PI = 2.50662827463; // sqrt(2*pi)

// [[Rcpp::export]]
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
SEXP FilterIndepOld  (Rcpp::NumericVector y_rcpp,
					Rcpp::NumericMatrix y_lagged_rcpp,
					Rcpp::NumericMatrix z_dependent_rcpp,
					Rcpp::NumericMatrix z_independent_rcpp,
					Rcpp::NumericVector beta_rcpp,
					Rcpp::NumericVector mu_rcpp,
					Rcpp::NumericVector sigma_rcpp,
					Rcpp::NumericMatrix gamma_dependent_rcpp,
					Rcpp::NumericVector gamma_independent_rcpp,
					Rcpp::NumericMatrix transition_probs_rcpp,
					Rcpp::NumericVector initial_dist_rcpp
					)
{
	int n = y_rcpp.size();
	int M = mu_rcpp.size();
	arma::mat xi_k_t(M, n); // make a transpose first for easier column operations.

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
	arma::mat    transition_probs(transition_probs_rcpp.begin(),
								transition_probs_rcpp.nrow(), transition_probs_rcpp.ncol(), false);
	arma::colvec initial_dist(initial_dist_rcpp.begin(), initial_dist_rcpp.size(), false);

	double likelihood = 0;

	SEXP eta_rcpp = EtaIndep(y_rcpp, y_lagged_rcpp,
						z_dependent_rcpp, z_independent_rcpp,
						beta_rcpp, mu_rcpp, sigma_rcpp,
						gamma_dependent_rcpp, gamma_independent_rcpp);
	arma::mat eta_t = (Rcpp::as<arma::mat>(eta_rcpp)).t();

	xi_k_t.col(0) = eta_t.col(0) % initial_dist;
	double total = sum(xi_k_t.col(0));
	xi_k_t.col(0) = xi_k_t.col(0) / total;
	likelihood += log(total);

	for (int k = 1; k < n; k++)
	{
	  	xi_k_t.col(k) = eta_t.col(k) % (transition_probs * xi_k_t.col(k-1));
	 	total = sum(xi_k_t.col(k));
		xi_k_t.col(k) = xi_k_t.col(k) / total;
		likelihood += log(total);
	}

	return Rcpp::List::create(Named("xi.k") = wrap(xi_k_t.t()),
														Named("likelihood") = wrap(likelihood));
}

// // [[Rcpp::export]]
// SEXP EstimateStates  (Rcpp::NumericVector y_rcpp,
// 				Rcpp::NumericMatrix y_lagged_rcpp,
// 				Rcpp::NumericMatrix z_dependent_rcpp,
// 				Rcpp::NumericMatrix z_independent_rcpp,
// 				Rcpp::NumericVector beta_rcpp,
// 				Rcpp::NumericVector mu_rcpp,
// 				Rcpp::NumericVector sigma_rcpp,
// 				Rcpp::NumericMatrix gamma_dependent_rcpp,
// 				Rcpp::NumericVector gamma_independent_rcpp)
// {
// 	SEXP eta_rcpp = EtaIndep(y_rcpp, y_lagged_rcpp,
// 						z_dependent_rcpp, z_independent_rcpp,
// 						beta_rcpp, mu_rcpp, sigma_rcpp,
// 						gamma_dependent_rcpp, gamma_independent_rcpp);
// 	arma::mat eta = Rcpp::as<arma::mat>(eta_rcpp);
//  int n = eta.n_cols;
//  int M = eta.n_rows;
//
// 	arma::colvec states(n);
//
// 	for (int i = 0; i < n; i++)
// 	{
// 		double best = -std::numeric_limits<double>::infinity();
// 		int best_index = -1;
// 		for (int j = 0; j < M; j++)
// 			if (eta(i,j) > best)
// 				best_index = j;
// 		states(i) = best_index; // CAUTION: the first state is numerized as zero.
// 	}
//
// 	return (wrap(states));
// }
