#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double EPS_ALMOST_ZERO = 0.0000001; // used to check almost zero
const double SQRT2PI = 2.50662827463100050241; // sqrt(2*pi)
const double LOG2PI_OVERTWO = 0.91893853320467274178; // (log(2*pi) / 2)

// A parameter that describes a model
struct Theta
{
	arma::mat 		beta;
	arma::colvec 	mu;								 // WATCH: can be in mat in VAR models
	arma::colvec	sigma;						 // WATCH: can be in mat in VAR models
	arma::mat 		gamma_dependent;
	arma::mat 		gamma_independent;
	arma::mat 		transition_probs;
	arma::colvec 	initial_dist;
  Theta() {}
	Theta(arma::mat _beta, arma::mat _mu, arma::mat _sigma,
				arma::mat _gamma_dependent, arma::mat _gamma_independent,
				arma::mat _transition_probs, arma::colvec _initial_dist)
	{
			beta = _beta;
			mu = _mu;
			sigma = _sigma;
			gamma_dependent = _gamma_dependent;
			gamma_independent = _gamma_independent;
			transition_probs = _transition_probs;
			initial_dist = _initial_dist;
	}
};

// Collection of xi_k, xi_n where each (i,j)th member of xi_k and xi_n
// represents p_{i|k}(j) and p_{i|n}(j).
struct Xi
{
	arma::mat xi_k;				// n by M matrix
	arma::mat xi_past_t; // n by M matrix
	arma::mat xi_n;				// n by M matrix
	double		likelihood;

	Xi (arma::mat xi_k_, arma::mat xi_past_t_, arma::mat xi_n_,
		double likelihood_)
	{
		xi_k = xi_k_;
		xi_past_t = xi_past_t_;
		xi_n = xi_n_;
		likelihood = likelihood_;
	}
};

// Check if every element of a given matrix is almost zero,
// determined by criterion EPS_ALMOST_ZERO. if it is, set it to a zero matrix.
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

// Returns n by 1 column that contains a minimum of each column in m by n mat
arma::colvec GetMinPerCol (arma::mat* pmatrix)
{
	int m = pmatrix->n_rows;
	int n = pmatrix->n_cols;
	arma::colvec min_col(n);

	for (int j = 0; j < n; j++)
	{
		double min_value = std::numeric_limits<double>::infinity();
		for (int i = 0; i < m; i++)
		{
			if (min_value > pmatrix->at(i,j))
				{
					min_value = pmatrix->at(i,j);
					min_col(j) = min_value;
				}
		}
	}

	return min_col;
}

// Returns an eta matrix
arma::mat EtaIndep (arma::colvec* py,
                    arma::mat* py_lagged,
                    arma::mat* pz_dependent,
                    arma::mat* pz_independent,
                    Theta* ptheta)
{
  int M = ptheta->transition_probs.n_rows;
  int n = py->n_rows;
  arma::mat eta(n, M);

  for (int j = 0; j < M; j++)
  {
    eta.col(j) = *py - *py_lagged * ptheta->beta -
      *pz_dependent * ptheta->gamma_dependent.col(j) -
      *pz_independent * ptheta->gamma_independent - ptheta->mu(j);
    eta.col(j) = eta.col(j) % eta.col(j); // element-wise multiplication
    eta.col(j) = exp(-eta.col(j) / (2 * (ptheta->sigma(j) * ptheta->sigma(j))));
    eta.col(j) = eta.col(j) / (SQRT2PI * ptheta->sigma(j));
  }

  return eta;
}

// Computes xi_k and likelihood at this stage, and returns
// an instance of Xi that contains xi_k, likelihood, and empty xi_n that needs
// to be computed in the smooth step.
Xi FilterIndep (arma::colvec* py,
                arma::mat* py_lagged,
                arma::mat* pz_dependent,
                arma::mat* pz_independent,
                Theta* ptheta)
{
  int n = py->n_rows;
  int M = ptheta->transition_probs.n_cols;
  arma::mat xi_k_t(M, n); // make a transpose first for easier column operations.
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

		if (k > 0)
			xi_past_t.col(k) = ptheta->transition_probs * xi_k_t.col(k-1);
		else
			xi_past_t.col(k) = ptheta->initial_dist;

		for (int j = 0; j < M; j++)
		{
			arma::colvec xi_k_t_jk = py->row(k) - py_lagged->row(k) * ptheta->beta -
	      pz_dependent->row(k) * ptheta->gamma_dependent.col(j) -
	      pz_independent->row(k) * ptheta->gamma_independent - ptheta->mu(j);
			xi_k_t(j,k) = xi_k_t_jk(0); // explicit gluing
	    xi_k_t(j,k) *= xi_k_t(j,k);
			xi_k_t(j,k) = xi_k_t(j,k) / (2 * (ptheta->sigma(j) * ptheta->sigma(j)));
			if (min_value > xi_k_t(j,k))
			{
				min_value = xi_k_t(j,k);
				min_index = j;
			}
			// SQRT2PI only matters in calculation of eta;
			// you can add it in the final log-likelihood.
			ratios[j] = xi_past_t(j, k) / ptheta->sigma(j);
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

// Returns xi_n.
arma::mat Smooth (arma::mat* xi_k,
									arma::mat* xi_past_t,
                  arma::mat* transition_probs)
{
  int n = xi_k->n_rows;
  int M = xi_k->n_cols;

  arma::mat xi_n_t(M, n);
  arma::mat xi_k_t = xi_k->t(); // produces a copy of a transpose
	arma::mat	transition_probs_t = transition_probs->t();
  xi_n_t.col(n-1) = xi_k_t.col(n-1);
  for (int k = (n-2); k >= 0; k--)
    xi_n_t.col(k) = xi_k_t.col(k) %
  (transition_probs_t * (xi_n_t.col(k+1) / (xi_past_t->col(k+1))));

  return xi_n_t.t();
}

// Returns Xi, which contains xi_k, xi_n, and likelihood at this stage.
Xi ExpectationStep(arma::colvec* py,
									arma::mat* py_lagged,
									arma::mat* pz_dependent,
									arma::mat* pz_independent,
									Theta* ptheta)
{
	Xi filter = FilterIndep(py, py_lagged, pz_dependent, pz_independent,
														ptheta);
	filter.xi_n = Smooth(&filter.xi_k, &filter.xi_past_t,
												&(ptheta->transition_probs));
	return filter;
}

// Returns estimates for states of every observation based on posterior prob.
arma::colvec EstimateStates (arma::colvec* py,
								arma::mat* py_lagged,
								arma::mat* pz_dependent,
								arma::mat* pz_independent,
								Theta* ptheta)
{
	arma::mat eta = EtaIndep(py, py_lagged, pz_dependent, pz_independent, ptheta);
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

	return states;
}

// Returns an estimate for initial_dist based on posterior probabilities.
arma::colvec ComputeInitialDist (arma::colvec* py,
								arma::mat* py_lagged,
								arma::mat* pz_dependent,
								arma::mat* pz_independent,
								Theta* ptheta)
{
	arma::mat eta = EtaIndep(py, py_lagged, pz_dependent, pz_independent, ptheta);
	int n = eta.n_rows;
	int M = eta.n_cols;

	arma::colvec initial_dist(M);

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
		initial_dist(best_index)++;
	}

	return (initial_dist / n);
}

// Compute a stationary distribution given a transition matrix.
// (Since we force every member of the matrix to be strictly positive,
// the corresponding MC is irreducible thus positive-recurrent as it is finite.)
arma::colvec ComputeStationaryDist (arma::mat* transition_probs, int M)
{
	arma::cx_vec eigval;
	arma::cx_mat eigvec;
	arma::eig_gen(eigval, eigvec, transition_probs->t()); // left eigenv, so take t.
	arma::vec stationary_dist(M);

	// Need to find which eigenvector has a eigenval. of one and extract real parts
	// find index
	int stationary_dist_index = -1;
	for (int i = 0; i < M; i++)
	  if (std::abs(eigval(i).real() - 1) < EPS_ALMOST_ZERO)
	    stationary_dist_index = i;
	// extract real
	for (int i = 0; i < M; i++)
	  stationary_dist(i) = eigvec(i,stationary_dist_index).real();

	stationary_dist = abs(stationary_dist) / sum(abs(stationary_dist));
	return stationary_dist;
}

// Returns an maximized theta based on computed xi_k and xi_n from an E-step.
Theta MaximizationStepIndep (arma::colvec* py,
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
														arma::colvec* sigma_first_step)
{
	int M = ptheta0->initial_dist.size();
	int s = ptheta0->beta.n_rows;
	int p_dep = ptheta0->gamma_dependent.n_rows;
	int p_indep = ptheta0->gamma_independent.n_rows;
	int n = py->n_rows;
	int n_minus_one = n - 1;

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
			// enforce ub/lb.
			transition_probs(i,j) = std::max(transition_probs(i,j), 0.02); // hard constraint
			transition_probs(i,j) = std::min(transition_probs(i,j), 0.98); // hard constraint
		}
		// normalize
		transition_probs.row(i) = transition_probs.row(i) /
															sum(transition_probs.row(i));
	}

	// 2. Estimation for beta, mu, sigma, gamma
	// 2-1. mu WATCH: this is specific for independent beta and scalar AR
	for (int j = 0; j < M; j++)
	  mu(j) = sum(pxi_n->col(j) %
						(*py - *py_lagged * ptheta0->beta -
							*pz_dependent * (ptheta0->gamma_dependent).col(j) -
							*pz_independent * ptheta0->gamma_independent)) / sum(pxi_n->col(j));


	// 2-2. gamma_dependent WATCH: specific for independent beta and scalar AR
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

	// 2-3. beta WATCH: this is specific for independent beta and scalar AR
	arma::mat 		beta_part_one(s, s);
	arma::colvec 	beta_part_two(s);
	for (int k = 0; k < n; k++)
	{
	  double prop_sum = 0;
	  for (int j = 0; j < M; j++)
	  {
			double prop = pxi_n->at(k,j) / (ptheta0->sigma(j) * ptheta0->sigma(j));
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

	// 2-4. gamma_independent WATCH: specific for independent beta and scalar AR
	if (!SetToZeroIfAlmostZero(&ptheta0->gamma_independent)) // validity check
	{
		arma::mat 		gamma_independent_part_one(p_indep, p_indep);
		arma::colvec 	gamma_independent_part_two(p_indep);
	  for (int k = 0; k < n; k++)
	  {
	    double prop_sum = 0;
	    for (int j = 0; j < M; j++)
	    {
	      double prop = pxi_n->at(k,j) / (ptheta0->sigma(j) * ptheta0->sigma(j));
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

	// 2-5. sigma (switching) WATCH: specific for independent beta and scalar AR
	for (int j = 0; j < M; j++)
	{
	  sigma(j) = 0;
	  for (int k = 0; k < n; k++)
	  {
	    // res is a scalar, with one element
	    arma::colvec res = py->at(k) -
	      py_lagged->row(k) * beta -
	      pz_independent->row(k) * gamma_independent -
	      pz_dependent->row(k) * gamma_dependent.col(j) -
	      mu(j);
	    sigma(j) += (pxi_n->at(k,j) / sum(pxi_n->col(j))) * res(0) * res(0);
	  }
		// impose the hard constraint; sigma >= 0.01 * sigma.hat
	  sigma(j) = std::max(sqrt(sigma(j)), 0.01 * sigma_first_step->at(j));
	}

	// 2-6. initial_dist
	initial_dist = (pxi_n->row(0)).t();

	Theta theta = Theta(beta, mu, sigma, gamma_dependent, gamma_independent,
                     transition_probs, initial_dist);
	return theta;
}

// [[Rcpp::export]]
SEXP EMIndepCPP  (Rcpp::NumericVector y_rcpp,
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
									double epsilon)
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

	Theta* thetas = new Theta[maxit]; // keep track, only the last will stay
	Theta theta;
	arma::colvec likelihoods(maxit);
	double likelihood;
	int index_exit = 0;

	Theta theta0(beta0, mu0, sigma0,
								gamma_dependent0, gamma_independent0,
								transition_probs0, initial_dist0);
	thetas[0] = theta0;
	likelihoods(0) = -std::numeric_limits<double>::infinity();
	theta = theta0;
	likelihood = likelihoods(0);

	// 1. EM iterations
	for (int i = 1; i < maxit; i++)
	{
		index_exit++;
		Xi 		e_step = ExpectationStep(&y, &y_lagged, &z_dependent, &z_independent,
													&thetas[i-1]);
		thetas[i] = MaximizationStepIndep(&y, &y_lagged, &z_dependent, &z_independent,
																	&y_lagged_t, &z_dependent_t, &z_independent_t,
																	&thetas[i-1],
																	&(e_step.xi_k), &(e_step.xi_past_t),
																	&(e_step.xi_n),
																	&theta0.sigma);

		likelihoods(i) = e_step.likelihood;

		if (std::abs(likelihoods(i) - likelihoods(i-1)) < epsilon)
			break;
	}
	theta = thetas[index_exit]; // copy the best theta
	likelihood = likelihoods(index_exit); // copy the most recent likelihood
	likelihoods = likelihoods.subvec(0,index_exit); // remove unused elements
	delete[] thetas; // clear memory

	// 2. state estimation
	arma::colvec states = EstimateStates(&y, &y_lagged,
                                      &z_dependent, &z_independent,
																			&theta); // XXX: the first state is zero.
	int n = states.n_elem;
	for (int i = 0; i < n; i++)
		states(i)++; // add one to make the first state one.

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
														Named("likelihoods") = wrap(likelihoods),
														Named("states") = wrap(states)
														);
}
