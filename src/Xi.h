#ifndef __XI__
#define __XI__

#include <RcppArmadillo.h>

// Collection of xi_k, xi_n where each (i,j)th member of xi_k and xi_n
// represents p_{i|k}(j) and p_{i|n}(j).
struct Xi
{
  arma::mat xi_k;				// n by M matrix
  arma::mat xi_past_t; // n by M matrix
  arma::mat xi_n;				// n by M matrix
  double		likelihood;

  Xi (arma::mat xi_k_, arma::mat xi_past_t_,
      arma::mat xi_n_, double likelihood_)
  {
    xi_k = xi_k_;
    xi_past_t = xi_past_t_;
    xi_n = xi_n_;
    likelihood = likelihood_;
  }
};

#endif
