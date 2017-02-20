/*
// Computes estimates for AR-MSMAH models with EM-algorithm.
// Written by Chiyoung Ahn
*/
#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include "Smooth.h"
#include "Theta.h"
#include "Utilities.h"
#include "Xi.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SQRT2PI = 2.50662827463100050241; // sqrt(2*pi)
const double LOG2PI_OVERTWO = 0.91893853320467274178; // (log(2*pi) / 2)

// Computes xi_k and likelihood at this stage, and returns
// an instance of Xi that contains xi_k, likelihood, and empty xi_n that needs
// to be computed in the smooth step.
Xi FilterARMSMAH (arma::colvec* py,
                arma::mat* py_lagged,
                arma::mat* pz_dependent,
                arma::mat* pz_independent,
                Theta* ptheta,
                arma::imat* pstate_conversion_mat,
                arma::mat* ptransition_probs_extended)
{
  double likelihood = 0;

  int n = py->n_rows;
  int M_reduced = ptheta->gamma_dependent.n_cols;
	int M = ptheta->initial_dist.size();
	int s = ptheta->beta.n_rows;
	int M_extended_block = IntPower(M_reduced, s);

  arma::mat xi_k_t(M, n, arma::fill::zeros); // make transpose first for easier column operations.
  arma::mat xi_past_t(M, n, arma::fill::zeros);
  arma::mat transition_probs_extended_t = ptransition_probs_extended->t();

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
			xi_past_t.col(k) = ptheta->initial_dist;
    xi_past_t.col(k) /= arma::sum(xi_past_t.col(k));

    for (int j_M = 0; j_M < M_reduced; j_M++)
		{
			for (int j_extra = 0; j_extra < M_extended_block; j_extra++)
			{
				int j = j_M * M_extended_block + j_extra;
				arma::colvec xi_k_t_jk = py->row(k) -
		      pz_dependent->row(k) * ptheta->gamma_dependent.col(j_M) -
		      pz_independent->row(k) * ptheta->gamma_independent - ptheta->mu(j_M);
				for (int lag = 0; lag < s; lag++)
				{
					int lagged_index = pstate_conversion_mat->at((lag + 1), j);
					xi_k_t_jk -= ptheta->beta.at(lag, j_M) *  // WATCH: beta switching
												(py_lagged->at(k, lag) - ptheta->mu(lagged_index));
				}

				xi_k_t(j,k) = xi_k_t_jk(0); // explicit gluing
		    xi_k_t(j,k) *= xi_k_t(j,k);
				xi_k_t(j,k) = xi_k_t(j,k) /
          (2 * (ptheta->sigma(j_M) * ptheta->sigma(j_M))); // WATCH: sigma switching

				if (min_value > xi_k_t(j,k) && xi_past_t(j, k) > 0)
				{
					min_value = xi_k_t(j,k);
					min_index = j;
				}

				// SQRT2PI only matters in calculation of eta;
				// you can remove it in the final log-likelihood.
				ratios[j] = xi_past_t(j, k) / ptheta->sigma(j_M); // WATCH: sigma switching
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

		likelihood += log(row_sum) - min_value + log(ratios[min_index]);

		delete[] ratios; // clear memory
	}


  likelihood -= n * LOG2PI_OVERTWO;
  ptheta->likelihood = likelihood;
  Xi xi(xi_k_t.t(), xi_past_t, arma::mat());
  return xi;
}

// Returns Xi, which contains xi_k, xi_n, and likelihood at this stage.
Xi ExpectationStepARMSMAH	(arma::colvec* py,
                          arma::mat* py_lagged,
                          arma::mat* pz_dependent,
                          arma::mat* pz_independent,
                          Theta* ptheta,
                          arma::imat* pstate_conversion_mat)
{
  arma::mat transition_probs_extended = GetExtendedTransitionProbs(
                              ptheta->transition_probs,
                             *pstate_conversion_mat);

  Xi filter = FilterARMSMAH(py, py_lagged, pz_dependent, pz_independent,
                          ptheta, pstate_conversion_mat,
                          &transition_probs_extended);
  filter.xi_n = Smooth(&filter.xi_k, &filter.xi_past_t,
                        &transition_probs_extended);
  return filter;
}

// Returns an maximized theta based on computed xi_k and xi_n from an E-step.
Theta MaximizationStepARMSMAH (arma::colvec* py,
                            arma::mat* py_lagged,
                            arma::mat* pz_dependent,
                            arma::mat* z_dependent_lagged_blocks,
                            arma::mat* pz_independent,
                            arma::mat* z_independent_lagged_blocks,
                            arma::mat* py_lagged_t,
                            arma::mat* pz_dependent_t,
                            arma::mat* pz_independent_t,
                            Theta* ptheta0,
                            arma::mat* pxi_k,
                            arma::mat* pxi_past_t,
                            arma::mat* pxi_n,
                            double transition_probs_min,
                            double transition_probs_max,
                            double sigma_min,
                            arma::imat* pstate_conversion_mat)
{
  int M_reduced = ptheta0->transition_probs.n_cols;
  int M = ptheta0->initial_dist.size();
  int s = ptheta0->beta.n_rows;
  int s1 = s + 1;
  int p_dep = ptheta0->gamma_dependent.n_rows;
  int p_indep = ptheta0->gamma_independent.n_rows;
  int n = py->n_rows;
  int n_minus_one = n - 1;

  arma::mat 		transition_probs(M_reduced, M_reduced, arma::fill::zeros);
  arma::mat 		transition_probs_extended(M, M, arma::fill::zeros);
  arma::mat 		beta(s, M_reduced, arma::fill::zeros); // s by M_reduced matrix WATCH: A case?
  arma::colvec 	mu(M_reduced, arma::fill::zeros);  // M_reduced-length vec
  arma::colvec 	sigma(M_reduced, arma::fill::zeros); // M_reduced-length vec WATCH: H case?
  arma::mat 		gamma_dependent(p_dep, M_reduced, arma::fill::zeros); // p_dep by M_reduced mat
  arma::colvec	gamma_independent(p_indep, arma::fill::zeros); // p_indep-length vec
  arma::colvec	initial_dist(M, arma::fill::zeros);

  // 0. do some computation first.
  arma::mat xi_n_sigmasq = *pxi_n; // pi_{kx}/sigma_(x_1)^2
  arma::mat remainders(n, M, arma::fill::zeros);
  for (int j = 0; j < M; j++)
  {
    // fill pi_jk / sigma_j^2
    int current_index = pstate_conversion_mat->at(0,j);
    xi_n_sigmasq.row(j) /= ptheta0->sigma(current_index) *
                            ptheta0->sigma(current_index); // WATCH: H case?
  }

  // 1. Estimation for transition_probs
  double* totals_i = new double[M_reduced]();

  for (int j = 0; j < M; j++)
  {
    int current_index = pstate_conversion_mat->at(0, j);
    int next_index = pstate_conversion_mat->at(1, j);

    for (int k = 0; k < n_minus_one; k++)
    {
      transition_probs(current_index, next_index) += pxi_n->at(k+1,j);
      totals_i[current_index] += pxi_n->at(k, j);
    }
  }

  for (int i = 0; i < M_reduced; i++)
  {
    transition_probs.row(i) /= totals_i[i];

    // hard bound
    for (int j = 0; j < M_reduced; j++)
    {
      transition_probs(i,j) = std::max(transition_probs(i,j),
                       transition_probs_min);
      transition_probs(i,j) = std::min(transition_probs(i,j),
                       transition_probs_max);
    }
    transition_probs.row(i) /= arma::sum(transition_probs.row(i));
  }
  delete[] totals_i;

  // 2. Estimation for beta, mu, sigma, gamma
  // 2-1. beta (switching)
  arma::cube 	beta_divisors(s, s, M_reduced, arma::fill::zeros); // WATCH: A case?
  arma::mat 	beta_numerators(s, M_reduced, arma::fill::zeros); // WATCH: A case?
  arma::mat   beta_j_jM_remainders(M, M_reduced, arma::fill::zeros); // WATCH: A case?

  for (int j = 0; j < M; j++)
  {
    int current_index = pstate_conversion_mat->at(0, j);
    for (int k = 0; k < n; k++)
    {
      // arma::colvec beta_jk_coeff = beta_skj_coeffs.slice(j).col(k);
      arma::colvec beta_jk_coeff(s, arma::fill::zeros);
      for (int lag = 0; lag < s; lag++)
      {
        int lagged_index = pstate_conversion_mat->at((lag + 1), j);
        beta_jk_coeff(lag) = py_lagged->at(k, lag) - ptheta0->mu(lagged_index);
      }

      beta_divisors.slice(current_index) += xi_n_sigmasq(k,j) *
        beta_jk_coeff * beta_jk_coeff.t();
      beta_numerators.col(current_index) += xi_n_sigmasq(k,j) * beta_jk_coeff *
        (py->at(k) - ptheta0->mu(current_index)); // FIXIT: assume no regressors
    }
  }

  for (int j_M = 0; j_M < M_reduced; j_M++)
    beta.col(j_M) = arma::solve(beta_divisors.slice(j_M), beta_numerators.col(j_M));


  // update remainders
  for (int j = 0; j < M; j++)
  {
    int current_index = pstate_conversion_mat->at(0, j);
    // fill remainders
    for (int k = 0; k < n; k++)
    {
      arma::colvec remainder_kj = py->row(k) -
        pz_dependent->row(k) * ptheta0->gamma_dependent.col(current_index) -
        pz_independent->row(k) * ptheta0->gamma_independent -
        ptheta0->mu(current_index); // explicit glowing
      for (int lag = 0; lag < s; lag++)
      {
        int lagged_index = pstate_conversion_mat->at((lag + 1), j);
        double beta_skj_coeff = py_lagged->at(k, lag) - ptheta0->mu(lagged_index);
        remainder_kj -= beta(lag, current_index) * beta_skj_coeff;
      }
      remainders(k,j) = remainder_kj(0);
    }
  }


  // 2-2. mu
  arma::colvec 	mu_divisors(M_reduced, arma::fill::zeros);
  arma::colvec 	mu_numerators(M_reduced, arma::fill::zeros);

  for (int j = 0; j < M; j++)
  {
    int current_index = pstate_conversion_mat->at(0, j);
    double xi_nj_sigmasq_sum = sum(xi_n_sigmasq.col(j));

    for (int j_M = 0; j_M < M_reduced; j_M++)
    {
      double mu_j_coeff = 0;

      if (current_index == j_M)
        mu_j_coeff--;

      for (int lag = 0; lag < s; lag++)
      {
        int lagged_index = pstate_conversion_mat->at((lag + 1), j);
        if (lagged_index == j_M)
          mu_j_coeff -= ptheta0->beta.at(lag, current_index); // WATCH: A case?
      }
      mu_divisors(j_M) += mu_j_coeff * mu_j_coeff * xi_nj_sigmasq_sum;
      mu_numerators(j_M) += mu_j_coeff * sum(xi_n_sigmasq.col(j) %
                            (mu_j_coeff * ptheta0->mu(j_M) - remainders.col(j)));
    }
  }
  mu = mu_numerators / mu_divisors;

  // update remainders
  for (int j = 0; j < M; j++)
  {
    int current_index = pstate_conversion_mat->at(0, j);
    remainders.col(j) += ptheta0->mu(current_index);
    remainders.col(j) -= mu(current_index);
    for (int lag = 0; lag < s; lag++)
    {
      int lagged_index = pstate_conversion_mat->at((lag + 1), j);
      remainders.col(j) += beta(lag, current_index) *
                            ptheta0->mu.at(lagged_index); // WATCH: A case?
      remainders.col(j) -= beta(lag, current_index) *
                            mu.at(lagged_index); // WATCH: A case?
    }
  }

  // // 2-2. gamma_dependent
  // if (!SetToZeroIfAlmostZero(&ptheta0->gamma_dependent)) // validity check
  //   for (int j = 0; j < M; j++)
  //   {
  //     arma::mat     gamma_dependent_part_one(p_dep, p_dep, arma::fill::zeros);
  //     arma::colvec  gamma_dependent_part_two(p_dep, arma::fill::zeros);
  //     for (int k = 0; k < n; k++)
  //     {
  //       gamma_dependent_part_one += pxi_n->at(k,j) *
  //         (pz_dependent_t->col(k) * pz_dependent->row(k));
  //       gamma_dependent_part_two += pxi_n->at(k,j) * pz_dependent_t->col(k) *
  //         (py->at(k) - py_lagged->row(k) * ptheta0->beta -
  //         pz_independent->row(k) * ptheta0->gamma_independent - mu(j));
  //     }
  //
  //     gamma_dependent.col(j) = inv(gamma_dependent_part_one) *
  //                               gamma_dependent_part_two;
  //   }
  //

  // // 2-4. gamma_independent
  // if (!SetToZeroIfAlmostZero(&ptheta0->gamma_independent)) // validity check
  // {
  //   arma::mat 		gamma_independent_part_one(p_indep, p_indep, arma::fill::zeros);
  //   arma::colvec 	gamma_independent_part_two(p_indep, arma::fill::zeros);
  //   for (int k = 0; k < n; k++)
  //   {
  //     double prop_sum = 0;
  //     for (int j = 0; j < M; j++)
  //     {
  //       double prop = pxi_n->at(k,j) / (ptheta0->sigma(0) * ptheta0->sigma(0));
  //       prop_sum += prop;
  //       gamma_independent_part_two += prop * pz_independent_t->col(k) *
  //         (py->at(k) -
  //           py_lagged->row(k) * beta -
  //           pz_dependent->row(k) * gamma_dependent.col(j) - mu(j));
  //     }
  //     gamma_independent_part_one += prop_sum *
  //       (pz_independent_t->col(k) * pz_independent->row(k));
  //   }
  //   gamma_independent = inv(gamma_independent_part_one) *
  //                         gamma_independent_part_two;
  // }
  //

  // 2-5. sigma (switching) : WATCH: H case?
  arma::colvec 	sigma_divisors(M_reduced, arma::fill::zeros); // WATCH: H case?
  arma::colvec 	sigma_numerators(M_reduced, arma::fill::zeros); // WATCH: H case?
  for (int j = 0; j < M; j++)
  {
    int current_index = pstate_conversion_mat->at(0,j);
    sigma_divisors(current_index) += sum(pxi_n->col(j));
    sigma_numerators(current_index) += sum(pxi_n->col(j) %
                          (remainders.col(j) % remainders.col(j)));
  }
  sigma = arma::sqrt(sigma_numerators / sigma_divisors);

  // impose the hard constraint
  for (int j = 0; j < M_reduced; j++)
    sigma(j) = std::max(sigma(j), sigma_min);

  // 2-6. initial_dist
  initial_dist = (pxi_n->row(0)).t();

  gamma_dependent = ptheta0->gamma_dependent; // FIXIT: assume no regressors.
  gamma_independent = ptheta0->gamma_independent; // FIXIT: assume no regressors.

  Theta theta = Theta(beta, mu, sigma, gamma_dependent, gamma_independent,
                     transition_probs, initial_dist);

  return theta;
}

// [[Rcpp::export]]
SEXP EMcppARMSMAH (Rcpp::NumericVector y_rcpp,
                  Rcpp::NumericMatrix y_lagged_rcpp,
                  Rcpp::NumericMatrix z_dependent_rcpp,
        					Rcpp::NumericMatrix z_independent_rcpp,
        					Rcpp::NumericMatrix z_dependent_lagged_rcpp,
        					Rcpp::NumericMatrix z_independent_lagged_rcpp,
                  Rcpp::NumericMatrix beta0_rcpp,
                  Rcpp::NumericVector mu0_rcpp,
                  Rcpp::NumericVector sigma0_rcpp,
                  Rcpp::NumericMatrix gamma_dependent0_rcpp,
                  Rcpp::NumericMatrix gamma_independent0_rcpp,
                  Rcpp::NumericMatrix transition_probs0_rcpp,
                  Rcpp::NumericVector initial_dist0_rcpp,
                  int maxit,
                  double epsilon,
                  double transition_probs_min,
                  double transition_probs_max,
                  double sigma_min,
                  Rcpp::IntegerMatrix state_conversion_mat_rcpp
        					)
{
  // conversion first
  arma::colvec y(y_rcpp.begin(), y_rcpp.size(), false);
  arma::mat    y_lagged(y_lagged_rcpp.begin(),
                y_lagged_rcpp.nrow(), y_lagged_rcpp.ncol(), false);
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
  arma::mat		 beta0(beta0_rcpp.begin(),
                beta0_rcpp.nrow(), beta0_rcpp.ncol(), false);
  arma::colvec mu0(mu0_rcpp.begin(), mu0_rcpp.size(), false);
  arma::colvec sigma0(sigma0_rcpp.begin(), sigma0_rcpp.size(), false);
  arma::mat    gamma_dependent0(gamma_dependent0_rcpp.begin(),
                 gamma_dependent0_rcpp.nrow(),
                 gamma_dependent0_rcpp.ncol(), false);
  arma::mat    gamma_independent0(gamma_independent0_rcpp.begin(),
                 gamma_independent0_rcpp.nrow(),
                 gamma_independent0_rcpp.ncol(), false);
  arma::mat    transition_probs0(transition_probs0_rcpp.begin(),
                transition_probs0_rcpp.nrow(),
                transition_probs0_rcpp.ncol(), false);
  arma::colvec initial_dist0(initial_dist0_rcpp.begin(),
                initial_dist0_rcpp.size(), false);
  arma::imat    state_conversion_mat(state_conversion_mat_rcpp.begin(),
                state_conversion_mat_rcpp.nrow(),
                state_conversion_mat_rcpp.ncol(), false);
  arma::mat y_lagged_t = y_lagged.t();
  arma::mat z_dependent_t = z_dependent.t();
  arma::mat z_independent_t = z_independent.t();

  arma::colvec likelihoods(maxit, arma::fill::zeros);
  int index_exit = -1;

  Theta theta_original(beta0, mu0, sigma0,
             gamma_dependent0, gamma_independent0,
             transition_probs0, initial_dist0);
  Theta theta_updated = theta_original;
  Theta theta = theta_updated;
  theta.likelihood = -std::numeric_limits<double>::infinity();
  likelihoods(0) = -std::numeric_limits<double>::infinity();

  // partition blocks
  int s = beta0.n_rows;
  int p = gamma_dependent0.n_rows;
  int q = gamma_independent0.size();
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

  // 1. EM iterations
  for (int i = 1; i < maxit; i++)
  {
    index_exit++;
    Xi 		e_step = ExpectationStepARMSMAH(&y, &y_lagged, &z_dependent, &z_independent,
                          &theta_updated, &state_conversion_mat);

    likelihoods(i) = theta_updated.likelihood;

    // Stop if 1. the difference in likelihoods is small enough
    // or 2. likelihood decreases (a decrease is due to locating local max.
    // out of hard constraints in maximization step, which suggests that this
    // is not a good candidate anyway.)
    // Use negative condition to stop if 3. the updated likelihood is NaN.
    if (!((theta_updated.likelihood - theta.likelihood) >= epsilon))
      break;

    theta = theta_updated;
    theta_updated = MaximizationStepARMSMAH(&y, &y_lagged,
                                  &z_dependent, z_dependent_lagged_blocks,
                                  &z_independent, z_independent_lagged_blocks,
                                  &y_lagged_t, &z_dependent_t, &z_independent_t,
                                  &theta,
                                  &(e_step.xi_k), &(e_step.xi_past_t),
                                  &(e_step.xi_n),
                                  transition_probs_min,
                                  transition_probs_max,
                                  sigma_min,
                                  &state_conversion_mat);
  }

  // remove unused elements
  likelihoods = likelihoods.subvec(0, std::min((index_exit + 1), (maxit - 1)));

  // clear memory for blocks
  delete[] z_dependent_lagged_blocks;
  delete[] z_independent_lagged_blocks;


  Rcpp::List theta_R = Rcpp::List::create(Named("beta") = wrap(theta.beta),
                            Named("mu") = wrap(theta.mu),
                            Named("sigma") = wrap(theta.sigma),
                            Named("gamma.dependent") =
                              wrap(theta.gamma_dependent),
                            Named("gamma.independent") =
                              wrap(theta.gamma_independent),
                            Named("transition.probs") =
                              wrap(theta.transition_probs),
                            Named("initial.dist") =
                              wrap(theta.initial_dist)
                            );

  return Rcpp::List::create(Named("theta") = wrap(theta_R),
                            Named("likelihood") = wrap(theta.likelihood),
                            Named("likelihoods") = wrap(likelihoods));
}
