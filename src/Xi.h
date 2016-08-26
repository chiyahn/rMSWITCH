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

  Xi (arma::mat _xi_k, arma::mat _xi_past_t,
      arma::mat _xi_n);

  Xi (arma::mat _xi_k, arma::mat _xi_past_t,
      arma::mat _xi_n, double _likelihood);
};

#endif
