#ifndef __UTILITIES__
#define __UTILITIES__

#include <RcppArmadillo.h>

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

#endif
