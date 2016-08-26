#include "Theta.h"

Theta::Theta() {}
Theta::Theta(arma::mat _beta, arma::mat _mu, arma::mat _sigma,
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
    likelihood = 0; // set to zero as it is impossible.
}
