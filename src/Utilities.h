#ifndef __UTILITIES__
#define __UTILITIES__

#include <RcppArmadillo.h>
#include "Theta.h"

const double EPS_ALMOST_ZERO = 0.0000001; // used to check almost zero

// Check if every element of a given matrix is almost zero,
// determined by criterion EPS_ALMOST_ZERO. if it is, set it to a zero matrix.
inline bool SetToZeroIfAlmostZero (arma::mat* pmatrix)
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

// Assume exp > 0.
inline int IntPower(int base, int exp)
{
  int result = 1;
  while (exp)
  {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }

  return result;
}

inline arma::mat GetExtendedTransitionProbs(arma::mat transition_probs,
                                            arma::imat state_conversion_mat)
{
  int s = state_conversion_mat.n_rows - 1;
  int M = transition_probs.n_cols;
  int M_extended = state_conversion_mat.n_cols;
  int M_to_s = IntPower(M, s);
  arma::mat transition_probs_extended(M_extended, M_extended, arma::fill::zeros);
  for (int j = 0; j < M_extended; j++)
  {
    int sub_index = j % M_to_s;
    int last_state = state_conversion_mat.at(0, j);

    for (int i = 0; i < M; i++)
      transition_probs_extended(j, (M_to_s * i + sub_index)) = transition_probs.at(last_state, i);
  }
  return (transition_probs_extended);
}

inline arma::mat ExtendedTransMatToReduced(arma::mat transition_probs_extended,
                                        arma::imat state_conversion_mat,
                                        int M_reduced)
{
  int s = state_conversion_mat.n_rows - 1;
  int M_extended = state_conversion_mat.n_cols;
  int M_to_s = IntPower(M_reduced, s);
  arma::mat transition_probs(M_reduced, M_reduced, arma::fill::zeros);
  for (int j = 0; j < M_extended; j++)
  {
    int sub_index = j % M_to_s;
    int last_state = state_conversion_mat.at(0, j);

    for (int i = 0; i < M_reduced; i++)
      transition_probs(last_state, i) += transition_probs_extended.at(j, (M_to_s * i + sub_index));
  }
  transition_probs /= M_to_s; // smooth it by averaging out
  return (transition_probs);
}

#endif
