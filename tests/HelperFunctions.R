library(testthat)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sourceCpp("cppGetMinPerCol.cpp")
expect_that(as.matrix(c(1,2)), equals(GetMinPerCol(matrix(c(1,4,5,2), ncol = 2))))

test_that("Check ThetaToColumn & ColumnToTheta", {
  ## ThetaToColumn & ColumnToTheta (Case 1)
  p12 <- 0.6
  p21 <- 0.7
  transition.probs <- matrix(c((1-p12),p12,p21,(1-p21)), ncol = 2, byrow = T)
  initial.dist <- c(0.3,0.7)
  beta <- c(0.3, 0.4)
  mu <- c(-2, 2)
  sigma <- c(1, 2)
  gamma.dependent <- NULL
  gamma.independent <- NULL
  
  
  M <- ncol(transition.probs)
  s <- length(beta)
  is.beta.switching <- FALSE
  is.sigma.switching <- TRUE
  p.dep <- 0
  p.indep <- 0
  
  initial.dist.index <- M * M + 1
  beta.index <- M + initial.dist.index
  mu.index <- s * ifelse(is.beta.switching, M, 1) + beta.index
  sigma.index <- M + mu.index # TODO: where does sigma begin in a vectorized theta?
  gamma.dep.index <- ifelse(is.sigma.switching, M, 1) + sigma.index
  gamma.indep.index <- p.dep * M + gamma.dep.index
  
  theta <- list(transition.probs = transition.probs, initial.dist = initial.dist,
                beta = beta, mu = mu, sigma = sigma)
  
  # Transform a vectorized theta back to a list form
  ColumnToTheta <- function(theta.vectorized)
  {
    gamma.dependent <- NULL
    gamma.independent <- NULL
    if (!gamma.dep.index == gamma.indep.index)
      gamma.dependent <- matrix(theta.vectorized[gamma.dep.index:
                                                   (gamma.indep.index - 1)], ncol = M)
    if (!gamma.indep.index == length(theta.vectorized) + 1)
      gamma.independent <- theta.vectorized[gamma.indep.index:
                                              length(theta.vectorized)]
    return (list
            (transition.probs = matrix(theta.vectorized[1:(M*M)], ncol = M, byrow = T),
            initial.dist = theta.vectorized[initial.dist.index:(beta.index - 1)],
            beta = theta.vectorized[beta.index:(mu.index - 1)],
            mu = theta.vectorized[mu.index:(sigma.index - 1)],
            sigma = theta.vectorized[sigma.index:(gamma.dep.index - 1)],
            gamma.dependent = gamma.dependent,
            gamma.independent = gamma.independent
            ))
  }
  
  
  theta.back <- ColumnToTheta(ThetaToColumn(theta))
  expect_that(theta$transition.probs, equals(theta.back$transition.probs))
  expect_that(theta$initial.dist, equals(theta.back$initial.dist))
  expect_that(theta$beta, equals(theta.back$beta))
  expect_that(theta$mu, equals(theta.back$mu))
  expect_that(theta$sigma, equals(theta.back$sigma))
  
  ## ThetaToColumn & ColumnToTheta (Case 2)
  p12 <- 0.6
  p21 <- 0.7
  transition.probs <- matrix(c((1-p12),p12,p21,(1-p21)), ncol = 2, byrow = T)
  initial.dist <- c(0.3,0.7)
  mu <- c(-2, 2)
  sigma <- c(1, 2)
  beta <- c(0.3, 0.4)
  
  gamma.dependent <- c(0.5, 0.3)
  gamma.independent <- NULL
  M <- ncol(transition.probs)
  s <- length(beta)
  is.beta.switching <- FALSE
  is.sigma.switching <- TRUE
  p.dep <- 2
  p.indep <- 0
  
  initial.dist.index <- M * M + 1
  beta.index <- M + initial.dist.index
  mu.index <- s * ifelse(is.beta.switching, M, 1) + beta.index
  sigma.index <- M + mu.index # TODO: where does sigma begin in a vectorized theta?
  gamma.dep.index <- ifelse(is.sigma.switching, M, 1) + sigma.index
  gamma.indep.index <- p.dep * M + gamma.dep.index
  
  theta <- list(transition.probs = transition.probs, initial.dist = initial.dist,
                beta = beta, mu = mu, sigma = sigma,
                gamma.dependent = gamma.dependent, gamma.independent = gamma.independent)
  
  
  # Transform a vectorized theta back to a list form
  ColumnToTheta <- function(theta.vectorized)
  {
    gamma.dependent <- NULL
    gamma.independent <- NULL
    if (!gamma.dep.index == gamma.indep.index)
      gamma.dependent <- matrix(theta.vectorized[gamma.dep.index:
                                                   (gamma.indep.index - 1)], ncol = M)
    if (!gamma.indep.index == length(theta.vectorized) + 1)
      gamma.independent <- theta.vectorized[gamma.indep.index:
                                              length(theta.vectorized)]
    return (list
            (transition.probs = matrix(theta.vectorized[1:(M*M)], ncol = M, byrow = T),
            initial.dist = theta.vectorized[initial.dist.index:(beta.index - 1)],
            beta = theta.vectorized[beta.index:(mu.index - 1)],
            mu = theta.vectorized[mu.index:(sigma.index - 1)],
            sigma = theta.vectorized[sigma.index:(gamma.dep.index - 1)],
            gamma.dependent = gamma.dependent,
            gamma.independent = gamma.independent
            ))
  }
  
  
  theta.back <- ColumnToTheta(ThetaToColumn(theta))
  expect_that(theta$transition.probs, equals(theta.back$transition.probs))
  expect_that(theta$initial.dist, equals(theta.back$initial.dist))
  expect_that(theta$beta, equals(theta.back$beta))
  expect_that(theta$mu, equals(theta.back$mu))
  expect_that(theta$sigma, equals(theta.back$sigma))
  
  
  ## ThetaToColumn & ColumnToTheta (Case 3)
  p12 <- 0.4
  p13 <- 0.3
  p21 <- 0.3
  p23 <- 0.2
  p31 <- 0.3
  p32 <- 0.4
  transition.probs <- matrix(c((1-p12-p13),p12,p13,
                               p21,(1-p21-p23),p23,
                               p31,p32,(1-p31-p32)), 
                             ncol = 3, byrow =T)
  initial.dist <- c(0.3,0.7,0.1)
  mu <- c(-2, 2, 1)
  sigma <- c(1, 2, 0.5)
  beta <- c(0.3, 0.4)
  gamma.dependent <- NULL
  gamma.independent <- c(0.5, 0.3, 0.2)
  M <- ncol(transition.probs)
  s <- length(beta)
  is.beta.switching <- FALSE
  is.sigma.switching <- TRUE
  p.dep <- 0
  p.indep <- 3
  
  initial.dist.index <- M * M + 1
  beta.index <- M + initial.dist.index
  mu.index <- s * ifelse(is.beta.switching, M, 1) + beta.index
  sigma.index <- M + mu.index # TODO: where does sigma begin in a vectorized theta?
  gamma.dep.index <- ifelse(is.sigma.switching, M, 1) + sigma.index
  gamma.indep.index <- p.dep * M + gamma.dep.index
  
  theta <- list(transition.probs = transition.probs, initial.dist = initial.dist,
                beta = beta, mu = mu, sigma = sigma,
                gamma.dependent = gamma.dependent, gamma.independent = gamma.independent)
  
  
  # Transform a vectorized theta back to a list form
  ColumnToTheta <- function(theta.vectorized)
  {
    gamma.dependent <- NULL
    gamma.independent <- NULL
    if (!gamma.dep.index == gamma.indep.index)
      gamma.dependent <- matrix(theta.vectorized[gamma.dep.index:
                                                   (gamma.indep.index - 1)], ncol = M)
    if (!gamma.indep.index == length(theta.vectorized) + 1)
      gamma.independent <- theta.vectorized[gamma.indep.index:
                                              length(theta.vectorized)]
    return (list
            (transition.probs = matrix(theta.vectorized[1:(M*M)], ncol = M, byrow = T),
            initial.dist = theta.vectorized[initial.dist.index:(beta.index - 1)],
            beta = theta.vectorized[beta.index:(mu.index - 1)],
            mu = theta.vectorized[mu.index:(sigma.index - 1)],
            sigma = theta.vectorized[sigma.index:(gamma.dep.index - 1)],
            gamma.dependent = gamma.dependent,
            gamma.independent = gamma.independent
            ))
  }
  
  
  theta.back <- ColumnToTheta(ThetaToColumn(theta))
  expect_that(theta$transition.probs, equals(theta.back$transition.probs))
  expect_that(theta$initial.dist, equals(theta.back$initial.dist))
  expect_that(theta$beta, equals(theta.back$beta))
  expect_that(theta$mu, equals(theta.back$mu))
  expect_that(theta$sigma, equals(theta.back$sigma))
  
  ## ThetaToColumn & ColumnToTheta (Case 4)
  beta <- c(0.2,0.3,0.5)
  gamma.dependent <- c(0.2, 0.3)
  gamma.independent <- c(0.5, 0.3, 0.2)
  M <- ncol(transition.probs)
  s <- length(beta)
  is.beta.switching <- FALSE
  is.sigma.switching <- TRUE
  p.dep <- 2
  p.indep <- 3
  
  initial.dist.index <- M * M + 1
  beta.index <- M + initial.dist.index
  mu.index <- s * ifelse(is.beta.switching, M, 1) + beta.index
  sigma.index <- M + mu.index # TODO: where does sigma begin in a vectorized theta?
  gamma.dep.index <- ifelse(is.sigma.switching, M, 1) + sigma.index
  gamma.indep.index <- p.dep * M + gamma.dep.index
  
  theta <- list(transition.probs = transition.probs, initial.dist = initial.dist,
                beta = beta, mu = mu, sigma = sigma,
                gamma.dependent = gamma.dependent, gamma.independent = gamma.independent)
  
  
  # Transform a vectorized theta back to a list form
  ColumnToTheta <- function(theta.vectorized)
  {
    gamma.dependent <- NULL
    gamma.independent <- NULL
    if (!gamma.dep.index == gamma.indep.index)
      gamma.dependent <- matrix(theta.vectorized[gamma.dep.index:
                                                   (gamma.indep.index - 1)], ncol = M)
    if (!gamma.indep.index == length(theta.vectorized) + 1)
      gamma.independent <- theta.vectorized[gamma.indep.index:
                                              length(theta.vectorized)]
    return (list
            (transition.probs = matrix(theta.vectorized[1:(M*M)], ncol = M, byrow = T),
            initial.dist = theta.vectorized[initial.dist.index:(beta.index - 1)],
            beta = theta.vectorized[beta.index:(mu.index - 1)],
            mu = theta.vectorized[mu.index:(sigma.index - 1)],
            sigma = theta.vectorized[sigma.index:(gamma.dep.index - 1)],
            gamma.dependent = gamma.dependent,
            gamma.independent = gamma.independent
            ))
  }
  
  
  theta.back <- ColumnToTheta(ThetaToColumn(theta))
  expect_that(theta$transition.probs, equals(theta.back$transition.probs))
  expect_that(theta$initial.dist, equals(theta.back$initial.dist))
  expect_that(theta$beta, equals(theta.back$beta))
  expect_that(theta$mu, equals(theta.back$mu))
  expect_that(theta$sigma, equals(theta.back$sigma))
  
  
  ## ThetaToColumn & ColumnToTheta (Case 5: switching beta)
  beta <- matrix(c(0.2,0.3,0.5,0.2,0.1,0.3), ncol = M)
  s <- nrow(beta)
  is.beta.switching <- TRUE
  
  initial.dist.index <- M * M + 1
  beta.index <- M + initial.dist.index
  mu.index <- s * ifelse(is.beta.switching, M, 1) + beta.index
  sigma.index <- M + mu.index # TODO: where does sigma begin in a vectorized theta?
  gamma.dep.index <- ifelse(is.sigma.switching, M, 1) + sigma.index
  gamma.indep.index <- p.dep * M + gamma.dep.index
  
  theta <- list(transition.probs = transition.probs, initial.dist = initial.dist,
                beta = beta, mu = mu, sigma = sigma,
                gamma.dependent = gamma.dependent, gamma.independent = gamma.independent)
  
  
  # Transform a vectorized theta back to a list form
  ColumnToTheta <- function(theta.vectorized)
  {
    gamma.dependent <- NULL
    gamma.independent <- NULL
    if (!gamma.dep.index == gamma.indep.index)
      gamma.dependent <- matrix(theta.vectorized[gamma.dep.index:
                                                   (gamma.indep.index - 1)], ncol = M)
    if (!gamma.indep.index == length(theta.vectorized) + 1)
      gamma.independent <- theta.vectorized[gamma.indep.index:
                                              length(theta.vectorized)]
    return (list
            (transition.probs = matrix(theta.vectorized[1:(M*M)], ncol = M, byrow = T),
            initial.dist = theta.vectorized[initial.dist.index:(beta.index - 1)],
            beta = matrix(theta.vectorized[beta.index:(mu.index - 1)], 
                          ncol = ifelse(is.beta.switching, M, 1)),
            mu = theta.vectorized[mu.index:(sigma.index - 1)],
            sigma = theta.vectorized[sigma.index:(gamma.dep.index - 1)],
            gamma.dependent = gamma.dependent,
            gamma.independent = gamma.independent
            ))
  }
  
  
  theta.back <- ColumnToTheta(ThetaToColumn(theta))
  expect_that(theta$transition.probs, equals(theta.back$transition.probs))
  expect_that(theta$initial.dist, equals(theta.back$initial.dist))
  expect_that(theta$beta, equals(theta.back$beta))
  expect_that(theta$mu, equals(theta.back$mu))
  expect_that(theta$sigma, equals(theta.back$sigma))
})