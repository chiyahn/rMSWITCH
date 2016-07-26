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

test_that("ReducedColumn&TakeItBack", {

  CheckValidityReducedColumn <- function(theta.one, theta.two)
  {
    theta.one.vectorized <- ThetaToReducedColumn(theta.one)
    theta.two.vectorized <- ThetaToReducedColumn(theta.two)
    theta.onetwo.list <- ReducedColumnsToThetas(cbind(theta.one.vectorized,
                                                theta.two.vectorized),
                                                theta0 = theta.one)
    expect_equal(theta.onetwo.list[[1]], theta.one)
    expect_equal(theta.onetwo.list[[2]], theta.two)
  }

  theta.one <- RandomTheta()
  theta.two <- RandomTheta()
  CheckValidityReducedColumn(theta.one, theta.two)

  theta.one <- RandomTheta(M = 3, s = 2, p.dep = 2, p.indep = 1)
  theta.two <- RandomTheta(M = 3, s = 2, p.dep = 2, p.indep = 1)
  CheckValidityReducedColumn(theta.one, theta.two)

  theta.one <- RandomTheta(M = 5, s = 4, p.dep = 2)
  theta.two <- RandomTheta(M = 5, s = 4, p.dep = 2)
  CheckValidityReducedColumn(theta.one, theta.two)

  theta.one <- RandomTheta(M = 5, s = 4, p.indep = 3)
  theta.two <- RandomTheta(M = 5, s = 4, p.indep = 3)
  CheckValidityReducedColumn(theta.one, theta.two)

})


test_that("ReducedStataColumn&TakeItBack", {
  #' ThetaToReducedColumn, as is presented in Stata
  ThetaToReducedStataColumn <- function(theta)
  {
    # transition.probs could have been listed in the order of
    # p11, p21, ..., pM1, p12, ..., pM2, ..., p1(M-1), ..., pM(M-1)
    # taking a transpose will make it listed as
    # p11, p12, ..., p1(M-1), p21, ..., p2(M-1), ..., pM(M-1)
    M <- ncol(theta$transition.probs)
    reduced.transition.probs <- theta$transition.probs[,1:(M-1)]
    return (c(c(theta$beta),
              c(theta$gamma.dependent),
              c(theta$gamma.independent),
              c(theta$mu), c(theta$sigma), c(t(reduced.transition.probs))))

  }

  CheckValidityReducedStataColumn <- function(theta.one, theta.two)
  {
    theta.one.vectorized.stata <- ThetaToReducedStataColumn(theta.one)
    theta.two.vectorized.stata <- ThetaToReducedStataColumn(theta.two)
    theta.onetwo.vectorized <- ReducedStataColumnsToReducedColumns(
      cbind(theta.one.vectorized.stata, theta.two.vectorized.stata),
      theta0 = theta.one)
    theta.onetwo.list <- ReducedColumnsToThetas(theta.onetwo.vectorized,
                                               theta.one)

    # information about initial.dist is not retreived in experiments.
    theta.onetwo.list <- lapply(theta.onetwo.list, function(theta)
    {
      theta$initial.dist <- NULL
      return (theta)
    })
    theta.onetwo.list.original <- lapply(list(theta.one, theta.two),
    function(theta)
    {
      theta$initial.dist <- NULL
      return (theta)
    })
    expect_equal(theta.onetwo.list[[1]], theta.onetwo.list.original[[1]])
    expect_equal(theta.onetwo.list[[2]], theta.onetwo.list.original[[2]])
  }

  theta.one <- RandomTheta()
  theta.two <- RandomTheta()
  CheckValidityReducedStataColumn(theta.one, theta.two)

  theta.one <- RandomTheta(M = 3, s = 2, p.dep = 2, p.indep = 1)
  theta.two <- RandomTheta(M = 3, s = 2, p.dep = 2, p.indep = 1)
  CheckValidityReducedStataColumn(theta.one, theta.two)

  theta.one <- RandomTheta(M = 5, s = 4, p.dep = 2)
  theta.two <- RandomTheta(M = 5, s = 4, p.dep = 2)
  CheckValidityReducedStataColumn(theta.one, theta.two)

  theta.one <- RandomTheta(M = 5, s = 4, p.indep = 3)
  theta.two <- RandomTheta(M = 5, s = 4, p.indep = 3)
  CheckValidityReducedStataColumn(theta.one, theta.two)

})
