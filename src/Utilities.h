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

// TODO: DO I IMPLEMENT THIS? DECIDE IT.
// void MatrixToThetas (int theta_count arma::mat* ptheta_matrix, Theta* pthetas);
inline Theta ColumnToTheta (arma::colvec* ptheta_col)
{
  // TODO: FINISH THIS
  return Theta(); // stub
}
inline arma::colvec ThetaToColumn (Theta* ptheta, int col_length)
{
  arma::colvec col(col_length);
  arma::conv_to< arma::colvec >::from(ptheta->transition_probs);

  // TODO: FINISH THIS
  return arma::colvec(1); // stub
}

#endif
