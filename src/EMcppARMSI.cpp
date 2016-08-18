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
Xi FilterARMSI (arma::colvec* py,
                arma::mat* py_lagged,
                arma::mat* pz_dependent,
                arma::mat* pz_independent,
                Theta* ptheta)
{
  int n = py->n_rows;
  int M = ptheta->transition_probs.n_cols;
  arma::mat xi_k_t(M, n); // make transpose first for easier column operations.
  arma::mat xi_past_t(M, n);
  double likelihood = 0;

  for (int k = 0; k < n; k++)
  {
    // initial setting; keep track of minimum value and its index
    // to divide everything by the min. value in order to prevent
    // possible numerical errors when computing posterior probs.
    int min_index = -1;
    double min_value = std::numeric_limits<double>::infinity();
    double* ratios = new double[M];
    double row_sum = 0;
    double sigma = ptheta->sigma(0); // non-switching sigma

    if (k > 0)
      xi_past_t.col(k) = ptheta->transition_probs * xi_k_t.col(k-1);
    else
      xi_past_t.col(k) = ptheta->initial_dist;

    for (int j = 0; j < M; j++)
    {
      arma::colvec xi_k_t_jk = py->row(k) -
        py_lagged->row(k) * ptheta->beta -
        pz_dependent->row(k) * ptheta->gamma_dependent.col(j) -
        pz_independent->row(k) * ptheta->gamma_independent - ptheta->mu(j);
      xi_k_t(j,k) = xi_k_t_jk(0) * xi_k_t_jk(0); // explicit gluing
      xi_k_t(j,k) = xi_k_t(j,k) / (2 * (sigma * sigma));
      if (min_value > xi_k_t(j,k))
      {
        min_value = xi_k_t(j,k);
        min_index = j;
      }
      // SQRT2PI only matters in calculation of eta;
      // you can add it in the final log-likelihood.
      ratios[j] = xi_past_t(j, k) / sigma;
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

  Xi xi(xi_k_t.t(), xi_past_t, arma::mat(), likelihood);
  return xi;
}

// Returns Xi, which contains xi_k, xi_n, and likelihood at this stage.
Xi ExpectationStepARMSI	(arma::colvec* py,
                          arma::mat* py_lagged,
                          arma::mat* pz_dependent,
                          arma::mat* pz_independent,
                          Theta* ptheta)
{
  Xi filter = FilterARMSI(py, py_lagged, pz_dependent, pz_independent,
                            ptheta);
  filter.xi_n = Smooth(&filter.xi_k, &filter.xi_past_t,
                        &(ptheta->transition_probs));
  return filter;
}

// Returns an maximized theta based on computed xi_k and xi_n from an E-step.
Theta MaximizationStepARMSI (arma::colvec* py,
                            arma::mat* py_lagged,
                            arma::mat* pz_dependent,
                            arma::mat* pz_independent,
                            arma::mat* py_lagged_t,
                            arma::mat* pz_dependent_t,
                            arma::mat* pz_independent_t,
                            Theta* ptheta0,
                            arma::mat* pxi_k,
                            arma::mat* pxi_past_t,
                            arma::mat* pxi_n,
                            arma::colvec* sigma_first_step,
                            double transition_probs_min,
                            double transition_probs_max)
{
  int M = ptheta0->initial_dist.size();
  int s = ptheta0->beta.n_rows;
  int p_dep = ptheta0->gamma_dependent.n_rows;
  int p_indep = ptheta0->gamma_independent.n_rows;
  int n = py->n_rows;
  int n_minus_one = n - 1;
  double sigma0 = ptheta0->sigma(0); // non-switching sigma

  arma::mat 		transition_probs(M, M);
  arma::mat 		beta(s, M);
  arma::colvec 	mu(M);  // M-length vec
  arma::colvec 	sigma(1); // 1-length vec
  arma::mat 		gamma_dependent(p_dep, M); // p_dep by M mat
  arma::colvec	gamma_independent(p_indep); // p_indep-length vec
  arma::colvec	initial_dist(M);

  // 1. Estimation for transition_probs
  for (int i = 0; i < M; i++)
   {
    double total = 0;
    for (int k = 0; k < n_minus_one; k++)
      total += pxi_n->at(k,i);
    for (int j = 0; j < M; j++)
    {
      double prob_ij = 0;
      for (int k = 1; k < n; k++)
        prob_ij += pxi_n->at(k,j) *
                    ptheta0->transition_probs(i,j) * pxi_k->at((k-1),i) /
                    pxi_past_t->at(j,k);
      transition_probs(i,j) = prob_ij / total;
      // enforce ub/lb (hard constraints)
      transition_probs(i,j) = std::max(transition_probs(i,j),
                                        transition_probs_min);
      transition_probs(i,j) = std::min(transition_probs(i,j),
                                        transition_probs_max);
    }
    // normalize
    transition_probs.row(i) = transition_probs.row(i) /
                              sum(transition_probs.row(i));
  }

  // 2. Estimation for beta, mu, sigma, gamma
  // 2-1. mu
  for (int j = 0; j < M; j++)
    mu(j) = sum(pxi_n->col(j) %
            (*py - *py_lagged * ptheta0->beta -
              *pz_dependent * (ptheta0->gamma_dependent).col(j) -
              *pz_independent * ptheta0->gamma_independent)) /
              sum(pxi_n->col(j));


  // 2-2. gamma_dependent
  if (!SetToZeroIfAlmostZero(&ptheta0->gamma_dependent)) // validity check
    for (int j = 0; j < M; j++)
    {
      arma::mat     gamma_dependent_part_one(p_dep, p_dep);
      arma::colvec  gamma_dependent_part_two(p_dep);
      for (int k = 0; k < n; k++)
      {
        gamma_dependent_part_one += pxi_n->at(k,j) *
          (pz_dependent_t->col(k) * pz_dependent->row(k));
        gamma_dependent_part_two += pxi_n->at(k,j) * pz_dependent_t->col(k) *
          (py->at(k) - py_lagged->row(k) * ptheta0->beta -
          pz_independent->row(k) * ptheta0->gamma_independent - mu(j));
      }

      gamma_dependent.col(j) = inv(gamma_dependent_part_one) *
                                gamma_dependent_part_two;
    }

  // 2-3. beta (non-switching)
  arma::mat 		beta_part_one(s, s);
  arma::colvec 	beta_part_two(s);
  for (int k = 0; k < n; k++)
  {
    double prop_sum = 0;
    for (int j = 0; j < M; j++)
    {
      double prop = pxi_n->at(k,j) / (sigma0 * sigma0);
      prop_sum += prop;
      beta_part_two += prop * py_lagged_t->col(k) *
        (py->at(k) -
        pz_independent->row(k) * ptheta0->gamma_independent -
        pz_dependent->row(k) * gamma_dependent.col(j) - mu(j));
    }
    beta_part_one += prop_sum *
      (py_lagged_t->col(k) * (py_lagged->row(k)));
  }
  beta = inv(beta_part_one) * beta_part_two;

  // 2-4. gamma_independent
  if (!SetToZeroIfAlmostZero(&ptheta0->gamma_independent)) // validity check
  {
    arma::mat 		gamma_independent_part_one(p_indep, p_indep);
    arma::colvec 	gamma_independent_part_two(p_indep);
    for (int k = 0; k < n; k++)
    {
      double prop_sum = 0;
      for (int j = 0; j < M; j++)
      {
        double prop = pxi_n->at(k,j) / (sigma0 * sigma0);
        prop_sum += prop;
        gamma_independent_part_two += prop * pz_independent_t->col(k) *
          (py->at(k) -
            py_lagged->row(k) * beta -
            pz_dependent->row(k) * gamma_dependent.col(j) - mu(j));
      }
      gamma_independent_part_one += prop_sum *
        (pz_independent_t->col(k) * pz_independent->row(k));
    }
    gamma_independent = inv(gamma_independent_part_one) *
                          gamma_independent_part_two;
  }

  // 2-5. sigma (non-switching)
  for (int j = 0; j < M; j++)
    for (int k = 0; k < n; k++)
    {
      // res is a scalar, with one element
      arma::colvec res = py->at(k) -
        py_lagged->row(k) * beta -
        pz_independent->row(k) * gamma_independent -
        pz_dependent->row(k) * gamma_dependent.col(j) -
        mu(j);
      sigma(0) += pxi_n->at(k,j) * res(0) * res(0);
    }
  // impose the hard constraint; sigma >= 0.5 * sigma.hat
  sigma(0) = std::max(sqrt(sigma(0)), 0.5 * sigma_first_step->at(0));

  // 2-6. initial_dist
  initial_dist = (pxi_n->row(0)).t();

  Theta theta = Theta(beta, mu, sigma, gamma_dependent, gamma_independent,
                     transition_probs, initial_dist);
  return theta;
}

// [[Rcpp::export]]
SEXP EMcppARMSI (Rcpp::NumericVector y_rcpp,
                  Rcpp::NumericMatrix y_lagged_rcpp,
                  Rcpp::NumericMatrix z_dependent_rcpp,
                  Rcpp::NumericMatrix z_independent_rcpp,
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
                  double transition_probs_max)
{
  // conversion first
  arma::colvec y(y_rcpp.begin(), y_rcpp.size(), false);
  arma::mat    y_lagged(y_lagged_rcpp.begin(),
                y_lagged_rcpp.nrow(), y_lagged_rcpp.ncol(), false);
  arma::mat    z_dependent(z_dependent_rcpp.begin(),
                z_dependent_rcpp.nrow(), z_dependent_rcpp.ncol(), false);
  arma::mat    z_independent(z_independent_rcpp.begin(),
                z_independent_rcpp.nrow(), z_independent_rcpp.ncol(), false);
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
  arma::mat y_lagged_t = y_lagged.t();
  arma::mat z_dependent_t = z_dependent.t();
  arma::mat z_independent_t = z_independent.t();

  arma::colvec likelihoods(maxit);
  double likelihood;
  int index_exit = 0;

  Theta theta_original(beta0, mu0, sigma0,
             gamma_dependent0, gamma_independent0,
             transition_probs0, initial_dist0);
  Theta theta_temp = theta_original;
  Theta theta = theta_temp;

  likelihoods(0) = -std::numeric_limits<double>::infinity();
  likelihood = likelihoods(0);

  // 1. EM iterations
  for (int i = 1; i < maxit; i++)
  {
    index_exit++;
    Xi 		e_step = ExpectationStepARMSI(&y, &y_lagged, &z_dependent, &z_independent,
                          &theta_temp);

    likelihoods(i) = e_step.likelihood;

    // stop if 1. the difference in likelihoods is small enough
    // or 2. likelihood decreases (a decrease is due to locating local max.
    // out of hard constraints in maximization step, which suggests that this
    // is not a good candidate anyway.)
    if ((likelihoods(i) - likelihoods(i-1)) < epsilon)
      break;

    theta = theta_temp;
    theta_temp = MaximizationStepARMSI(&y, &y_lagged,
                                  &z_dependent, &z_independent,
                                  &y_lagged_t, &z_dependent_t, &z_independent_t,
                                  &theta,
                                  &(e_step.xi_k), &(e_step.xi_past_t),
                                  &(e_step.xi_n),
                                  &theta_original.sigma,
                                  transition_probs_min,
                                  transition_probs_max);

  }
  likelihood = likelihoods(index_exit); // copy the most recent likelihood
  likelihoods = likelihoods.subvec(0, index_exit); // remove unused elements

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
                            Named("likelihood") = wrap(likelihood),
                            Named("likelihoods") = wrap(likelihoods));
}
