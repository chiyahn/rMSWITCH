#ifndef __SMOOTH__
#define __SMOOTH__

#include <RcppArmadillo.h>

// Returns xi_n.
inline arma::mat Smooth (arma::mat* xi_k,
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
  {
    xi_n_t.col(k) = xi_k_t.col(k) %
      (transition_probs_t * (xi_n_t.col(k+1) / (xi_past_t->col(k+1))));
    xi_n_t.col(k) /= sum(xi_n_t.col(k)); // normalize
  }

  return xi_n_t.t();
}

#endif
