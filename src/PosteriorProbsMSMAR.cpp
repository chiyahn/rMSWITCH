/*
// Computes posterior probabilities given a AR-MSI model.
// Written by Chiyoung Ahn
*/
//#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include "Utilities.h"
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
SEXP PosteriorProbsMSMAR (Rcpp::NumericVector y_rcpp,
					Rcpp::NumericMatrix y_lagged_rcpp,
					Rcpp::NumericMatrix z_dependent_rcpp,
					Rcpp::NumericMatrix z_independent_rcpp,
					Rcpp::NumericMatrix z_dependent_lagged_rcpp,
					Rcpp::NumericMatrix z_independent_lagged_rcpp,
					Rcpp::NumericMatrix transition_probs_rcpp,
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
	arma::mat    z_dependent_lagged(z_dependent_lagged_rcpp.begin(),
								z_dependent_lagged_rcpp.nrow(),
								z_dependent_lagged_rcpp.ncol(), false);
	arma::mat    z_independent(z_independent_rcpp.begin(),
								z_independent_rcpp.nrow(),
								z_independent_rcpp.ncol(), false);
	arma::mat    z_independent_lagged(z_independent_lagged_rcpp.begin(),
								z_independent_lagged_rcpp.nrow(),
								z_independent_lagged_rcpp.ncol(), false);
	arma::mat    transition_probs(transition_probs_rcpp.begin(),
								transition_probs_rcpp.nrow(),
								transition_probs_rcpp.ncol(), false);
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

  arma::mat transition_probs_extended = GetExtendedTransitionProbs(
                              transition_probs, state_conversion_mat);
  arma::mat transition_probs_extended_t = transition_probs_extended.t();
	int M = transition_probs_extended_t.n_rows;
	int M_reduced = gamma_dependent.n_cols;
	int s = beta.n_rows;
	int M_extended_block = IntPower(M_reduced, s);
	arma::mat xi_k_t(M, n, arma::fill::zeros); // make a transpose first for col operations.
	arma::mat xi_past_t(M, n);

	// partition blocks
	int p = gamma_dependent.n_rows;
	int q = gamma_independent.size();
	arma::mat* z_dependent_lagged_blocks = new arma::mat[s];
	arma::mat* z_independent_lagged_blocks = new arma::mat[s];
	for (int i = 0; i < s; i++)
	{
		int z_dependent_block_first = i * p;
		int z_independent_block_first = i * q;
		z_dependent_lagged_blocks[i] = z_dependent_lagged.cols(z_dependent_block_first,
			z_dependent_block_first + p - 1);
		z_independent_lagged_blocks[i] = z_independent_lagged.cols(z_independent_block_first,
			z_independent_block_first + q - 1);
	}

	for (int k = 0; k < n; k++)
	{
		// initial setting; keep track of minimum value and its index
		// to divide everything by the min. value in order to prevent
		// possible numerical errors when computing posterior probs.
		int min_index = -1;
		double min_value = std::numeric_limits<double>::infinity();
		double* ratios = new double[M];
		double row_sum = 0;

		if (k > 0)
			xi_past_t.col(k) = transition_probs_extended_t * xi_k_t.col(k-1);
		else
			xi_past_t.col(k) = initial_dist_extended;
		xi_past_t.col(k) /= arma::sum(xi_past_t.col(k));

		for (int j_M = 0; j_M < M_reduced; j_M++)
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
					xi_k_t_jk -= beta.at(lag, j_M) *
						(y_lagged.at(k, lag) -
						z_dependent_lagged_blocks[lag].row(k) * gamma_dependent.col(lagged_index) -
						z_independent_lagged_blocks[lag].row(k) * gamma_independent -
						mu(lagged_index));
				}

				xi_k_t(j,k) = xi_k_t_jk(0); // explicit gluing
				xi_k_t(j,k) *= xi_k_t(j,k);
				xi_k_t(j,k) = xi_k_t(j,k) /
					(2 * sigma(j_M) * sigma(j_M)); // WATCH: sigma switching
				if (min_value > xi_k_t(j,k) && xi_past_t(j, k) > 0)
				{
					min_value = xi_k_t(j,k);
					min_index = j;
				}

				// SQRT2PI only matters in calculation of eta;
				// you can remove it in the final log-likelihood.
				ratios[j] = xi_past_t(j, k) / sigma(j_M); // WATCH: sigma switching

			}
		}

		for (int j = 0; j < M; j++)
		{
			if (j == min_index)
				xi_k_t(j,k) = 1.0;
			else
				xi_k_t(j,k) = (ratios[j] / ratios[min_index]) *
											exp(min_value - xi_k_t(j,k));
			row_sum += xi_k_t(j,k);
		}

		xi_k_t.col(k) /= row_sum;

		delete[] ratios; // clear memory
	}

	// clear memory for blocks
	delete[] z_dependent_lagged_blocks;
	delete[] z_independent_lagged_blocks;

	// smoothed probabilities
	arma::mat xi_n_t(M, n);
  xi_n_t.col(n-1) = xi_k_t.col(n-1);
  for (int k = (n-2); k >= 0; k--)
  {
    xi_n_t.col(k) = xi_k_t.col(k) %
  (transition_probs_extended_t * (xi_n_t.col(k+1) / (xi_past_t.col(k+1))));
    xi_n_t.col(k) /= sum(xi_n_t.col(k)); // normalize
  }

	// make 'simplified' xi_n based on current states
	arma::mat xi_k = xi_k_t.t();
	arma::mat xi_n = xi_n_t.t();
	arma::mat xi_k_simplified(n, M_reduced, arma::fill::zeros);
	arma::mat xi_n_simplified(n, M_reduced, arma::fill::zeros);

	for (int j = 0; j < M_reduced; j++)
		for (int j_block = 0; j_block < M_extended_block; j_block++)
		{
			xi_k_simplified.col(j) += xi_k.col(j * M_extended_block + j_block);
			xi_n_simplified.col(j) += xi_n.col(j * M_extended_block + j_block);
		}

	return Rcpp::List::create(Named("xi.k") = wrap(xi_k_simplified),
														Named("xi.n") = wrap(xi_n_simplified),
														Named("xi.k.extended") = wrap(xi_k),
														Named("xi.n.extended") = wrap(xi_n));
}
