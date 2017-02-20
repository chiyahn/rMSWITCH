library(rMSWITCH)
library(testthat)

test_that("OrderTransitionMatrix", {
  mu1.order <- order(c(1,-1))   # 2, 1
  mu2.order <- order(c(-1,1))   # 1, 2
  mu3.order <- order(c(3,-2,1)) # 2, 3, 1
  mu4.order <- order(c(-1,3,2)) # 1, 3, 2
  mu5.order <- order(c(1,2,3))  # 1, 2, 3
  transition.matrix1 <- matrix(c(0.4,0.7,0.6,0.3), ncol = 2)
  transition.matrix2 <- matrix(c(0.4,0.5,0.2,0.6,0.3,0.3,0,0.2,0.5), ncol = 3)
  expect_equal(OrderTransitionMatrix(transition.matrix1, mu1.order),
               matrix(c(0.3,0.6,0.7,0.4), ncol = 2))
  expect_equal(OrderTransitionMatrix(transition.matrix1, mu2.order),
               transition.matrix1)
  expect_equal(OrderTransitionMatrix(transition.matrix2, mu3.order),
               matrix(c(0.3,0.3,0.6,0.2,0.5,0,0.5,0.2,0.4), ncol = 3))
  expect_equal(OrderTransitionMatrix(transition.matrix2, mu4.order),
               matrix(c(0.4,0.2,0.5,0,0.5,0.2,0.6,0.3,0.3), ncol = 3))
  expect_equal(OrderTransitionMatrix(transition.matrix2, mu5.order),
               transition.matrix2)
})


test_that("GetExtendedTransitionProbs", {
  library(inline)
  GetExtendedTransitionProbsC <- cxxfunction  (signature(
    transition_probs_rcpp="matrix", state_conversion_mat_rcpp="matrix"),
    plugin="RcppArmadillo", body='
    arma::imat state_conversion_mat = Rcpp::as<arma::imat>(state_conversion_mat_rcpp);
    arma::mat transition_probs = Rcpp::as<arma::mat>(transition_probs_rcpp);
    int s = state_conversion_mat.n_rows - 1;
    int M = transition_probs.n_cols;
    int M_extended = state_conversion_mat.n_cols;
    
    int base = M; // calling IntPower(M, s)
    int exp = s;
    int M_to_s = 1;
    while (exp)
    {
      if (exp & 1)
        M_to_s *= base;
      exp >>= 1;
      base *= base;
    }

    arma::mat transition_probs_extended(M_extended, M_extended, arma::fill::zeros);
    for (int j = 0; j < M_extended; j++)
    {
    int sub_index = (j - j % M) / M;
    int last_state = state_conversion_mat.at(0, j);
    
    for (int i = 0; i < M; i++)
      transition_probs_extended(j, (M_to_s * i + sub_index)) = 
        transition_probs.at(last_state, i);
    }
    return wrap(transition_probs_extended);')
  
  transition.matrix1 <- matrix(c(0.4,0.7,0.6,0.3), ncol = 2)
  transition.matrix2 <- matrix(c(0.4,0.5,0.2,0.6,0.3,0.3,0,0.2,0.5), ncol = 3)
  s1 <- 2
  s2 <- 3
  M1 <- nrow(transition.matrix1)
  M2 <- nrow(transition.matrix2)
  state.conversion.mat1 <- GetStateConversionMat(M1, s1)
  state.conversion.mat2 <- GetStateConversionMat(M2, s2)
  state.conversion.mat1R <- GetStateConversionMatForR(M1, s1)
  state.conversion.mat2R <- GetStateConversionMatForR(M2, s2)
  
  extended1 <- GetExtendedTransitionProbsC(transition.matrix1, 
                                           state.conversion.mat1)
  extended1.expected <- GetExtendedTransitionProbs(transition.matrix1, 
                                                   state.conversion.mat1R)
  extended2 <- GetExtendedTransitionProbsC(transition.matrix2, 
                                           state.conversion.mat2)
  extended2.expected <- GetExtendedTransitionProbs(transition.matrix2, 
                                                   state.conversion.mat2R)
  expect_equal(extended1, extended1.expected)
  expect_equal(extended2, extended2.expected)
})
