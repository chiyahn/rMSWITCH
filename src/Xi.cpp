#include "Xi.h"

Xi::Xi (arma::mat _xi_k, arma::mat _xi_past_t,
    arma::mat _xi_n, double _likelihood)
{
  xi_k = _xi_k;
  xi_past_t = _xi_past_t;
  xi_n = _xi_n;
  likelihood = _likelihood;
}
