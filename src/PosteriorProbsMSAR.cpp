#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SQRT2PI = 2.50662827463; // sqrt(2*pi)
const double LOG2PI_OVERTWO = 0.91893853320467274178; // (log(2*pi) / 2)

// Returns n by M matrix that represents a posterior probs of observations for
//  univariate MS-AR model where beta is switching.
// Note that even if beta is non-switching, setting beta as a s by M matrix with
// repeated column of the original beta will give you the posterior probs
// matrix for MS-AR model with non-switching beta.
// [[Rcpp::export]]
SEXP PosteriorProbsMSAR (Rcpp::NumericVector y_rcpp,
					Rcpp::NumericMatrix y_lagged_rcpp,
					Rcpp::NumericMatrix z_dependent_rcpp,
					Rcpp::NumericMatrix z_independent_rcpp,
					Rcpp::NumericMatrix transition_probs_rcpp,
					Rcpp::NumericVector initial_dist_rcpp,
					Rcpp::NumericMatrix beta_rcpp,
					Rcpp::NumericVector mu_rcpp,
					Rcpp::NumericVector sigma_rcpp,
					Rcpp::NumericMatrix gamma_dependent_rcpp,
					Rcpp::NumericVector gamma_independent_rcpp
					)
{
	int n = y_rcpp.size();
	int M = mu_rcpp.size();
	arma::mat xi_k_t(M, n); // make a transpose first for easier column operations.

	arma::colvec y(y_rcpp.begin(), y_rcpp.size(), false);
	arma::mat    y_lagged(y_lagged_rcpp.begin(),
								y_lagged_rcpp.nrow(),
								y_lagged_rcpp.ncol(), false);
	arma::mat    z_dependent(z_dependent_rcpp.begin(),
								z_dependent_rcpp.nrow(),
								z_dependent_rcpp.ncol(), false);
	arma::mat    z_independent(z_independent_rcpp.begin(),
								z_independent_rcpp.nrow(),
								z_independent_rcpp.ncol(), false);
	arma::mat    transition_probs(transition_probs_rcpp.begin(),
								transition_probs_rcpp.nrow(),
								transition_probs_rcpp.ncol(), false);
	arma::colvec initial_dist(initial_dist_rcpp.begin(),
								initial_dist_rcpp.size(), false);
	arma::colvec beta(beta_rcpp.begin(), 
                   beta_rcpp.nrow(), beta_rcpp.ncol(), false);
	arma::colvec mu(mu_rcpp.begin(), mu_rcpp.size(), false);
	arma::colvec sigma(sigma_rcpp.begin(), sigma_rcpp.size(), false);
	arma::mat    gamma_dependent(gamma_dependent_rcpp.begin(),
								gamma_dependent_rcpp.nrow(),
								gamma_dependent_rcpp.ncol(), false);
	arma::colvec gamma_independent(gamma_independent_rcpp.begin(),
								gamma_independent_rcpp.size(), false);

	for (int k = 0; k < n; k++)
	{
		// initial setting; keep track of minimum value and its index
		// to divide everything by the min. value in order to prevent
		// possible numerical errors when computing posterior probs.
		int min_index = 0;
		double min_value = std::numeric_limits<double>::infinity();
		double* ratios = new double[M];

		arma::colvec xi_past;
		if (k > 0)
			xi_past = transition_probs * xi_k_t.col(k-1);
		else
			xi_past = transition_probs * initial_dist;

		for (int j = 0; j < M; j++)
		{
			arma::colvec xi_k_t_jk = y.row(k) - y_lagged.row(k) * beta.col(j) -
	      z_dependent.row(k) * gamma_dependent.col(j) -
	      z_independent.row(k) * gamma_independent - mu(j);
			xi_k_t(j,k) = xi_k_t_jk(0); // explicit gluing
	    xi_k_t(j,k) *= xi_k_t(j,k);
			xi_k_t(j,k) = xi_k_t(j,k) / (2 * (sigma(j) * sigma(j)));
			if (min_value > xi_k_t(j,k))
			{
				min_value = xi_k_t(j,k);
				min_index = j;
			}
			// SQRT2PI only matters in calculation of eta;
			// you can add it in the final log-likelihood.
			ratios[j] = xi_past(j) / sigma(j);
		}

		for (int j = 0; j < M; j++)
		{
			if (j == min_index)
				xi_k_t(j,k) = 1.0;
			else
				xi_k_t(j,k) = (ratios[j] / ratios[min_index]) *
				exp(min_value - xi_k_t(j,k));
		}
		delete[] ratios; // clear memory
	}

	return (wrap(xi_k_t.t()));
}
