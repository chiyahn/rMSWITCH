// WARNING: THIS HAS NOT BEEN TESTED YET

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double EPS_ALMOST_ZERO = 0.0000001; // used to check almost zero
const double SQRT2PI = 2.50662827463; // sqrt(2*pi)

// Check if a given vector is a zero matrix.
bool IsZero (arma::mat matrix)
{
	int m = matrix.n_rows;
	int n = matrix.n_cols;

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (std::abs(matrix(i,j)) > EPS_ALMOST_ZERO)
				return FALSE;

	return TRUE;
}

bool SetToZeroIfAlmostZero (arma::mat* pmatrix)
{
  int m = pmatrix->n_rows;
  int n = pmatrix->n_cols;

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (std::abs(pmatrix->at(i,j)) > EPS_ALMOST_ZERO)
        return FALSE;

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      pmatrix->at(i,j) = 0;

  return TRUE;
}

// Check if a given vector is a zero vector.
bool IsZero (arma::vec vector)
{
	int n = vector.n_elem;

	for (int i = 0; i < n; i++)
	  if (std::abs(vector(i)) > EPS_ALMOST_ZERO)
			return FALSE;

	return TRUE;
}

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


// [[Rcpp::export]]
SEXP MaximizationStepIndep  (Rcpp::NumericVector y_rcpp,
															Rcpp::NumericMatrix y_lagged_rcpp,
															Rcpp::NumericMatrix z_dependent_rcpp,
															Rcpp::NumericMatrix z_independent_rcpp,
															Rcpp::NumericMatrix beta0_rcpp,
															Rcpp::NumericVector mu0_rcpp,
															Rcpp::NumericVector sigma0_rcpp,
															Rcpp::NumericMatrix gamma_dependent0_rcpp,
															Rcpp::NumericVector gamma_independent0_rcpp,
															Rcpp::NumericMatrix transition_probs0_rcpp,
															Rcpp::NumericVector initial_dist0_rcpp,
															Rcpp::NumericMatrix xi_k_rcpp,
															Rcpp::NumericMatrix xi_n_rcpp,
															Rcpp::NumericVector sigma0_origin_rcpp
															)
{
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
								gamma_dependent0_rcpp.nrow(), gamma_dependent0_rcpp.ncol(), false);
	arma::colvec gamma_independent0(gamma_independent0_rcpp.begin(),
								gamma_independent0_rcpp.size(), false);
	arma::mat    transition_probs0(transition_probs0_rcpp.begin(),
								transition_probs0_rcpp.nrow(),
								transition_probs0_rcpp.ncol(), false);
	arma::colvec initial_dist0(initial_dist0_rcpp.begin(),
	initial_dist0_rcpp.size(), false);
	arma::mat    xi_k(xi_k_rcpp.begin(),
								xi_k_rcpp.nrow(), xi_k_rcpp.ncol(), false);
	arma::mat    xi_n(xi_n_rcpp.begin(),
								xi_n_rcpp.nrow(), xi_n_rcpp.ncol(), false);
	arma::colvec sigma0_origin(sigma0_origin_rcpp.begin(),
															sigma0_origin_rcpp.size(), false);

	arma::mat y_lagged_t = y_lagged.t(); // transpose of y_lagged
	arma::mat z_dependent_t = z_dependent.t(); // transpose of z_dependent
	arma::mat z_independent_t = z_independent.t(); // transpose of z_independent


	int M = initial_dist0.size();
	int s = beta0.n_rows;
	int p_dep = gamma_dependent0.n_rows;
	int p_indep = gamma_independent0.n_rows;
	int n = y.n_rows;

	arma::mat 		transition_probs(M, M);
	arma::mat 		beta(s, M); // s by M matrix (ind. case => beta is a single col.)
	arma::colvec 	mu(M);  // M-length vec
	arma::colvec 	sigma(M); // M-length vec
	arma::mat 		gamma_dependent(p_dep, M); // p_dep by M mat
	arma::colvec	gamma_independent(p_indep); // p_indep-length vec
	arma::colvec	initial_dist(M);

	// 1. Estimation for transition_probs
	for (int i = 0; i < M; i++)
   {
    double total = 0;
		for (int k = 0; k < n; k++)
    	total += xi_n(k,i);
  	for (int j = 0; j < M; j++)
  	{
    	double prob_ij = 0;
	    for (int k = 1; k < n; k++)
	    {
	      arma::colvec prob_ij_k = (transition_probs0 * (xi_k.row(k-1).t()));
	      prob_ij   += xi_n(k,j) * transition_probs0(i,j) * xi_k((k-1),i) / prob_ij_k(j);
	    }
	    transition_probs(i,j) = prob_ij / total;
	    transition_probs(i,j) = std::max(transition_probs(i,j), 0.02); // hard constraint
	    transition_probs(i,j) = std::min(transition_probs(i,j), 0.98); // hard constraint
    }
		// normalize
  	transition_probs.row(i) = transition_probs.row(i) / sum(transition_probs.row(i));
  }

	// 2. Estimation for beta, mu, sigma, gamma
	// 2-1. mu WATCH: this is specific for independent beta and scalar AR
	for (int j = 0; j < M; j++)
	  mu(j) = sum(xi_n.col(j) % (y - y_lagged * beta0)) / sum(xi_n.col(j));

	// 2-2. gamma_dependent WATCH: this is specific for independent beta and scalar AR
	if (!SetToZeroIfAlmostZero(&gamma_dependent)) // validity check
  	for (int j = 0; j < M; j++)
  	{
  	  arma::mat     gamma_dependent_part_one(p_dep, p_dep);
  	  arma::colvec  gamma_dependent_part_two(p_dep);
  	  for (int k = 0; k < n; k++)
  	  {
  	    gamma_dependent_part_one += xi_n(k,j) *
  	      (z_dependent_t.col(k) * z_dependent.row(k));
  	    gamma_dependent_part_two += xi_n(k,j) * z_dependent_t.col(k) *
  	      (y(k) - y_lagged.row(k) * beta0 - z_independent.row(k) * gamma_independent0 - mu(j));
  	  }

  	  gamma_dependent.col(j) = inv(gamma_dependent_part_one) * gamma_dependent_part_two;
  	}

	// 2-3. beta WATCH: this is specific for independent beta and scalar AR
	arma::mat 		beta_part_one(s, s);
	arma::colvec 	beta_part_two(s);
	for (int k = 0; k < n; k++)
	{
	  double prop_sum = 0;
	  for (int j = 0; j < M; j++)
	  {
			double prop = xi_n(k,j) / (sigma0(j) * sigma0(j));
	    prop_sum += prop;
	    beta_part_two += prop * y_lagged_t.col(k) *
				(y(k) -
				z_independent.row(k) * gamma_independent0 -
				z_dependent.row(k) * gamma_dependent.col(j) - mu[j]);
	  }
	  beta_part_one += prop_sum * (y_lagged_t.col(k) * (y_lagged.row(k)));
	}
	beta = inv(beta_part_one) * beta_part_two;

	// 2-4. gamma_independent WATCH: this is specific for independent beta and scalar AR
	if (!SetToZeroIfAlmostZero(&gamma_independent)) // validity check
	{
		arma::mat 		gamma_independent_part_one(p_indep, p_indep);
		arma::colvec 	gamma_independent_part_two(p_indep);
	  for (int k = 0; k < n; k++)
	  {
	    double prop_sum = 0;
	    for (int j = 0; j < M; j++)
	    {
	      double prop = xi_n(k,j) / (sigma0(j) * sigma0(j));
	      prop_sum += prop;
	      gamma_independent_part_two += prop * z_independent_t.col(k) *
					(y[k] -
						y_lagged.row(k) * beta -
						z_dependent.row(k) * gamma_dependent.col(j) - mu(j));
	    }
	    gamma_independent_part_one += prop_sum *
				(z_independent_t.col(k) * z_independent.row(k));
	  }
	  gamma_independent = inv(gamma_independent_part_one) * gamma_independent_part_two;
	}

	// 2-5. sigma (switching) WATCH: this is specific for independent beta and scalar AR
	for (int j = 0; j < M; j++)
	{
	  sigma(j) = 0;
	  for (int k = 0; k < n; k++)
	  {
	    // res is a scalar, with one element
	    arma::colvec res = y(k) -
	      y_lagged.row(k) * beta -
	      z_independent.row(k) * gamma_independent -
	      z_dependent.row(k) * gamma_dependent.col(j) -
	      mu(j);
	    sigma(j) += (xi_n(k,j) / sum(xi_n.col(j))) * res(0) * res(0);
	  }
		// impose the hard constraint; sigma >= 0.01 * sigma.hat
	  sigma(j) = std::max(sqrt(sigma(j)), 0.01 * sigma0_origin(j));
	}

	// 2-6. initial_dist
	SEXP states_rcpp = EstimateStates(y_rcpp, y_lagged_rcpp,
						z_dependent_rcpp, z_independent_rcpp,
						wrap(beta), wrap(mu), wrap(sigma),
						wrap(gamma_dependent), wrap(gamma_independent));

	arma::colvec states = Rcpp::as<arma::colvec>(states_rcpp);

	for (int i = 0; i < n; i++)
	 	initial_dist(((int)states(i)))++; // WATCH: might truncate.
	initial_dist /= n;

	return Rcpp::List::create(Named("beta") = wrap(beta),
														Named("mu") = wrap(mu),
														Named("sigma") = wrap(sigma),
														Named("gamma.dependent") = wrap(gamma_dependent),
														Named("gamma.independent") = wrap(gamma_independent),
														Named("transition.probs") = wrap(transition_probs),
														Named("initial.dist") = wrap(initial_dist)
														);
}
