# This applies only to univariate time series.
ThetaToColumn <- function(theta)
{
  # transition.probs will be listed in the order of
  # p11, p21, ..., pM1, p12, ..., pM2, ..., pMM if we don't use transpose.
  # taking a transpose will make it listed as
  # p11, p12, ..., p1M, p21, ..., p2M, ..., pMM
  return (c(c(t(theta$transition.probs)),
            c(theta$initial.dist),
            c(theta$beta), c(theta$mu), c(theta$sigma),
            c(theta$gamma.dependent),
            c(theta$gamma.independent)))

}

# returns a list of (theta likelihood) where theta is a list of parameters
# and likelihood is a likelihood of the data using the parameters in theta;
# this applies for univariate time series only.
MaximizeLongStep <- function(candidates, y, y.lagged, z.dependent, z.independent)
{
  # use the first candidate to save the information about dimensions
  theta <- candidates[[1]]
  n <- length(y)
  M <- ncol(theta$transition.probs)
  s <- nrow(as.matrix(theta$beta))
  is.beta.switching <- (ncol(as.matrix(theta$beta)) > 1)
  is.sigma.switching <- (length(theta$sigma) > 1)
  p.dep <- 0
  p.indep <- 0
  if (!is.null(theta$gamma.dependent))
    p.dep <- nrow(theta$gamma.dependent)
  else
    z.dependent <- as.matrix(rep(0,n))
  if (!is.null(theta$gamma.independent))
    p.dep <- nrow(theta$gamma.independent)
  else
    z.independent <- as.matrix(rep(0,n))
  
  # if beta is not switching, make it as a switching par;
  # replicate the column to every M to form length(beta) by M matrix
  if (!is.beta.switching)
    candidates <- lapply(candidates,
                          function (theta)
                          { theta$beta <- matrix(rep(theta$beta, M), ncol = M)
                            return (theta)})


  # this holds only for univariate time series.
  initial.dist.index <- M * M + 1
  beta.index <- M + initial.dist.index
  mu.index <- M * s + beta.index
  # mu.index <- s * ifelse(is.beta.switching, M, 1) + beta.index
  sigma.index <- M + mu.index 

  # if gamma.dependent does not exist, should have the same value as gamma.indep.index
  gamma.dep.index <- ifelse(is.sigma.switching, M, 1) + sigma.index
  # if gamma.independent does not exist, should have the same value as length(theta.vectorized) + 1
  gamma.indep.index <- p.dep * M + gamma.dep.index

  # Transform a vectorized theta back to a list form
  ColumnToTheta <- function(theta.vectorized)
  {
    beta <- theta.vectorized[beta.index:(mu.index - 1)]
    if (is.beta.switching)
      beta <- matrix(beta, ncol = M)
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
            beta = beta,
            mu = theta.vectorized[mu.index:(sigma.index - 1)],
            sigma = theta.vectorized[sigma.index:(gamma.dep.index - 1)],
            gamma.dependent = gamma.dependent,
            gamma.independent = gamma.independent
            ))
  }

  # define it dynamically (for indices)
  SLSQPNonSwitchingAR <- function(theta.vectorized,
                                  y, y.lagged, z.dependent, z.independent)
  {
    # dynamically define a constraint function
    ConstraintMCTransition <- function(theta.vectorized)
    {
      constraint.vectorized <- vector()

      # transition.probs
      for (i in 1:M)
      {
        index <- 1 + (i - 1) * M
        sum.ith <- sum(theta.vectorized[index:(index + M - 1)])
        constraint.vectorized <- c(constraint.vectorized,
                                   1-sum.ith)
      }

      # initial.dist
      constraint.vectorized <- c(constraint.vectorized,
                                 1-sum(theta.vectorized[(M*M+1):(M*M+M)]))

      return (constraint.vectorized)
    }

    ObjectiveLogLikelihood <- function(theta.vectorized)
    {
      gamma.dependent <- t(rep(0,M))
      gamma.independent <- 0
      if (gamma.dep.index != gamma.indep.index) # i.e. gamma.dependent exists
        gamma.dependent   <- matrix(theta.vectorized[gamma.dep.index:
                                                       (gamma.indep.index - 1)], ncol = M)
      if (gamma.indep.index <= length(theta.vectorized)) # i.e. gamma.independent exists
        gamma.independent <- theta.vectorized[gamma.indep.index:length(theta.vectorized)]

      # slsqp solves a minimization problem;
      # take a negative value to turn the problem into max. problem
      -LikelihoodMSAR(y, y.lagged, z.dependent, z.independent,
                      matrix(theta.vectorized[1:(M*M)], ncol = M, byrow = T),   # trans.probs
                      c(theta.vectorized[initial.dist.index:(beta.index - 1)]),  # initial.dist
                      beta = matrix(theta.vectorized[beta.index:(mu.index - 1)],
                                    ncol = M),  # beta
                      theta.vectorized[mu.index:(sigma.index - 1)],           # mu
                      theta.vectorized[sigma.index:(gamma.dep.index - 1)],    # sigma
                      gamma.dependent,
                      gamma.independent) # gamma.indep
    }

    # hard constraints to prevent values from bounding off
    transition.probs.lb <- rep(0.02, M*M)
    transition.probs.ub <- rep(0.98, M*M)
    beta <- theta.vectorized[beta.index:(mu.index - 1)]
    beta.lb <- beta * (1 - 0.2 * sign(beta))
    beta.ub <- beta * (1 + 0.2 * sign(beta))
    mu <- theta.vectorized[mu.index:(sigma.index - 1)]
    mu.lb <- mu * (1 - 0.5 * sign(mu))
    mu.ub <- mu * (1 + 0.5 * sign(mu))
    sigma <- theta.vectorized[sigma.index:(gamma.dep.index - 1)]
    sigma.lb <- sigma * 0.5
    sigma.ub <- sigma * 1.5

    # hard constraints on
    # transtion.probs (ub & lb), initial.dist (ub & lb), and sigma (lb)
    lb <- c(transition.probs.lb, rep(0.02, M),
            beta.lb,
            mu.lb,
            sigma.lb)
    if (!length(theta.vectorized) == length(lb)) # if gamma.dep and gamma.indep remain
      lb <- c(lb, rep(-Inf, (length(theta.vectorized) - length(lb))))
    ub <- c(transition.probs.ub, rep(1, M),
            beta.ub,
            mu.ub,
            sigma.ub)
    ub <- c(ub, rep(Inf, (length(theta.vectorized) - length(ub))))

    result <- slsqp(theta.vectorized,
                    fn = ObjectiveLogLikelihood,
                    lower = lb, upper = ub, heq = ConstraintMCTransition)
    result$value <- -result$value # take negative back to make it actual
    return (result)
  }


  candidates.matrix <- sapply(candidates, function (theta) ThetaToColumn(theta))
  
  # perform non-linear optimization in each candidate
  long.results <- apply(candidates.matrix, 2,
                        SLSQPNonSwitchingAR,
                        y = y, y.lagged = y.lagged,
                        z.dependent = z.dependent, z.independent = z.independent)

  long.likelihoods <- unlist(lapply(long.results, "[[", "value"))
  
  # extract the one that returns the best likelihood
  long.result <- long.results[[(which(long.likelihoods==max(long.likelihoods))[1])]]

  return (list(theta = ColumnToTheta(long.result$par),
              likelihood = long.result$value))
}
