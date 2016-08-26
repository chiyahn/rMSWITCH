#ifndef __THETA__
#define __THETA__

#include <RcppArmadillo.h>

// A parameter that describes a model
struct Theta
{
  arma::mat 		beta;
  arma::colvec 	mu;								// WATCH: can be in mat in VAR models
  arma::colvec	sigma;						// WATCH: can be in mat in VAR models
  arma::mat 		gamma_dependent;
  arma::mat 		gamma_independent;
  arma::mat 		transition_probs;
  arma::colvec 	initial_dist;
  double likelihood;              // initial value is zero

  Theta();
  Theta(arma::mat _beta, arma::mat _mu, arma::mat _sigma,
        arma::mat _gamma_dependent, arma::mat _gamma_independent,
        arma::mat _transition_probs, arma::colvec _initial_dist);
};

#endif
