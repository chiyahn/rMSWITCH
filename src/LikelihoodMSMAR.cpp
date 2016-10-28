#include <RcppArmadillo.h>
#include "Utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SQRT2PI = 2.50662827463; // sqrt(2*pi)
const double LOG2PI_OVERTWO = 0.91893853320467274178; // (log(2*pi) / 2)

// Returns likelihood of univariate MSM-AR model where beta is switching.
// Note that even if beta is non-switching, setting beta as a s by M matrix with
// repeated column of the original beta will give you the likelihood for
// MS-AR model with non-switching beta.
// TODO: Implement the version where z_dependent/z_independent exist(lagged variables should be implemented for MSM models)
// [[Rcpp::export]]
SEXP LikelihoodMSMAR (Rcpp::NumericVector y_rcpp,
					Rcpp::NumericMatrix y_lagged_rcpp,
					Rcpp::NumericMatrix z_dependent_rcpp,
					Rcpp::NumericMatrix z_independent_rcpp,
					Rcpp::NumericMatrix transition_probs_extended_rcpp,
					Rcpp::NumericVector initial_dist_extended_rcpp,
					Rcpp::NumericMatrix beta_rcpp,
					Rcpp::NumericVector mu_rcpp,
					Rcpp::NumericVector sigma_rcpp,
					Rcpp::NumericMatrix gamma_dependent_rcpp,
					Rcpp::NumericVector gamma_independent_rcpp,
					Rcpp::IntegerMatrix state_conversion_mat_rcpp
					)
{
	int n = y_rcpp.size();

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
	arma::mat    transition_probs_extended(transition_probs_extended_rcpp.begin(),
								transition_probs_extended_rcpp.nrow(),
								transition_probs_extended_rcpp.ncol(), false);
	arma::colvec initial_dist_extended(initial_dist_extended_rcpp.begin(),
								initial_dist_extended_rcpp.size(), false);
	arma::mat    beta(beta_rcpp.begin(),
                   beta_rcpp.nrow(), beta_rcpp.ncol(), false);
	arma::colvec mu(mu_rcpp.begin(), mu_rcpp.size(), false);
	arma::colvec sigma(sigma_rcpp.begin(), sigma_rcpp.size(), false);
	arma::mat    gamma_dependent(gamma_dependent_rcpp.begin(),
								gamma_dependent_rcpp.nrow(),
								gamma_dependent_rcpp.ncol(), false);
	arma::colvec gamma_independent(gamma_independent_rcpp.begin(),
								gamma_independent_rcpp.size(), false);
	arma::imat    state_conversion_mat(state_conversion_mat_rcpp.begin(),
								state_conversion_mat_rcpp.nrow(),
								state_conversion_mat_rcpp.ncol(), false);

  arma::mat transition_probs_extended_t = transition_probs_extended.t();
  double likelihood = 0;
	int M_extended = transition_probs_extended_t.n_rows;
	int M = gamma_dependent.n_cols;
	int s = beta.n_rows;
	int s_plus_one = s + 1;
	int M_extended_block = IntPower(M, s);
	arma::mat xi_k_t(M_extended, n, arma::fill::zeros); // make a transpose first for col operations.

	for (int k = 0; k < n; k++)
	{
		// initial setting; keep track of minimum value and its index
		// to divide everything by the min. value in order to prevent
		// possible numerical errors when computing posterior probs.
		int min_index = -1;
		double min_value = std::numeric_limits<double>::infinity();
		double* ratios = new double[M_extended];
		double row_sum = 0;

		arma::colvec xi_past;
		if (k > 0)
			xi_past = transition_probs_extended_t * xi_k_t.col(k-1);
		else
			xi_past = initial_dist_extended;
    xi_past /= arma::sum(xi_past);

		for (int j_M = 0; j_M < M; j_M++)
		{
			for (int j_extra = 0; j_extra < M_extended_block; j_extra++)
			{
				int j = j_M * M_extended_block + j_extra;
				arma::colvec xi_k_t_jk = y.row(k) -
		      z_dependent.row(k) * gamma_dependent.col(j_M) -
		      z_independent.row(k) * gamma_independent - mu(j_M);
				for (int lag = 0; lag < s; lag++)
				{
					int lagged_index = state_conversion_mat.at((lag + 1), j);
					xi_k_t_jk -= beta.at(lag, lagged_index) *
												(y_lagged.at(k, lag) - mu(lagged_index));
				}
				xi_k_t(j,k) = xi_k_t_jk(0); // explicit gluing
		    xi_k_t(j,k) *= xi_k_t(j,k);
				xi_k_t(j,k) = xi_k_t(j,k) / (2 * (sigma(j_M) * sigma(j_M)));

				if (min_value > xi_k_t(j,k))
				{
					min_value = xi_k_t(j,k);
					min_index = j;
				}
				// SQRT2PI only matters in calculation of eta;
				// you can remove it in the final log-likelihood.
				ratios[j] = xi_past(j) / sigma(j_M);
			}
		}

		for (int j = 0; j < M_extended; j++)
		{
			if (j == min_index)
				xi_k_t(j,k) = 1.0;
			else
				xi_k_t(j,k) = (ratios[j] / ratios[min_index]) *
											exp(min_value - xi_k_t(j,k));
			row_sum += xi_k_t(j,k);
		}
		xi_k_t.col(k) /= row_sum;

		likelihood += log(row_sum) - min_value + log(ratios[min_index]);

		delete[] ratios; // clear memory
	}
	likelihood -= n * LOG2PI_OVERTWO;

	return (wrap(likelihood));
}
