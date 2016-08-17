# returns a list of (theta log.likelihood) where theta is a list of parameters
# and log.likelihood is a log.likelihood of the data using
# the parameters in theta; this applies for univariate time series only.
MaximizeLongStep <- function(candidates, y, y.lagged,
                            z.dependent, z.independent,
                            epsilon, maxit,
                            transition.probs.min = 0.01,
                            lb.prob.density = 10e-6,
                            ub.prob.density = (1-10e-6))
{
  # use the first candidate to save the information about dimensions
  theta <- candidates[[1]]
  n <- length(y)
  M <- ncol(theta$transition.probs)
  s <- nrow(as.matrix(theta$beta))
  is.beta.switching <- (ncol(as.matrix(theta$beta)) > 1)
  is.sigma.switching <- (length(theta$sigma) > 1)
  p.dep <- 1 # even if gamma.dependent is NULL, use a zero vector
  p.indep <- 1 # same reason.
  if (!is.null(z.dependent))
    p.dep <- nrow(as.matrix(theta$gamma.independent))
  else
    z.dependent <- as.matrix(rep(0,n))
  if (!is.null(z.independent))
    p.indep <- nrow(as.matrix(theta$gamma.independent))
  else
    z.independent <- as.matrix(rep(0,n))

  # this holds only for univariate time series.
  initial.dist.index <- M * (M-1) + 1 # reduced case
  beta.index <- (M-1) + initial.dist.index # reduced case
  mu.index <- s * ifelse(is.beta.switching, M, 1) + beta.index
  sigma.index <- M + mu.index

  # if gamma.dependent does not exist,
  # should have the same value as (gamma.indep.index - M)
  gamma.dep.index <- ifelse(is.sigma.switching, M, 1) + sigma.index
  # if gamma.independent does not exist,
  # should have the same value as length(theta.vectorized) + 1
  gamma.indep.index <- p.dep * M + gamma.dep.index

  # Transform a vectorized theta back to a list form
  ColumnToTheta <- function(theta.vectorized)
  {
    beta <- theta.vectorized[beta.index:(mu.index - 1)]
    if (is.beta.switching)
      beta <- matrix(beta, ncol = M)
    gamma.dependent <- NULL
    gamma.independent <- NULL
    if (gamma.dep.index != (gamma.indep.index - M))
      gamma.dependent <- matrix(theta.vectorized[gamma.dep.index:
                                                   (gamma.indep.index - 1)],
                                ncol = M)
    if (gamma.indep.index != length(theta.vectorized) + 1)
      gamma.independent <- theta.vectorized[gamma.indep.index:
                                              length(theta.vectorized)]
    return (list
            (transition.probs = matrix(theta.vectorized[1:(M*M)],
                                        ncol = M, byrow = T),
            initial.dist = theta.vectorized[initial.dist.index:
                                            (beta.index - 1)],
            beta = beta,
            mu = theta.vectorized[mu.index:(sigma.index - 1)],
            sigma = theta.vectorized[sigma.index:(gamma.dep.index - 1)],
            gamma.dependent = gamma.dependent,
            gamma.independent = gamma.independent
            ))
  }

  ReducedColumnToTheta <- function(theta.vectorized)
  {
    transition.probs <- matrix(theta.vectorized[1:(M*(M-1))],
                              ncol = (M-1), byrow = T)
    # revive the original from the reduced form.
    transition.probs <- t(apply(transition.probs, 1,
                                function (row) c(row, (1-sum(row)))))
    initial.dist <- theta.vectorized[initial.dist.index:(beta.index - 1)]
    initial.dist <- c(initial.dist, (1 - sum(initial.dist)))
    beta <- theta.vectorized[beta.index:(mu.index - 1)]
    if (is.beta.switching)
      beta <- matrix(beta, ncol = M)
    gamma.dependent <- NULL
    gamma.independent <- NULL
    if (gamma.dep.index != (gamma.indep.index - M))
      gamma.dependent <- matrix(theta.vectorized[gamma.dep.index:
                                                   (gamma.indep.index - 1)],
                                ncol = M)
    if (gamma.indep.index != length(theta.vectorized) + 1)
      gamma.independent <- theta.vectorized[gamma.indep.index:
                                              length(theta.vectorized)]

    return (list
            (transition.probs = transition.probs,
            initial.dist = initial.dist,
            beta = beta,
            mu = theta.vectorized[mu.index:(sigma.index - 1)],
            sigma = theta.vectorized[sigma.index:(gamma.dep.index - 1)],
            gamma.dependent = gamma.dependent,
            gamma.independent = gamma.independent
            ))
  }

  # dynamically define a constraint function
  ConstraintMCTransition <- function(theta.vectorized)
  {
    constraint.vectorized <- vector()

    # transition.probs
    for (i in 1:M)
    {
      index <- 1 + (i - 1) * (M-1)
      sum.ith <- sum(theta.vectorized[index:(index + M - 2)])
      constraint.vectorized <- c(constraint.vectorized,
                                 (1-transition.probs.min)-sum.ith)
    }

    # initial.dist (indices are M*(M-1)+1:M*(M-1)+(M-1))
    constraint.vectorized <- c(constraint.vectorized,
                               ub.prob.density-
                                 sum(theta.vectorized[(M*(M-1)+1):(M*M-1)]))

    return (constraint.vectorized)
  }

  # define it dynamically (for indices)
  SLSQPNonSwitchingAR <- function(theta.vectorized,
                                  y, y.lagged, z.dependent, z.independent)
  {
    # sanity check; if a candidate contains a singularity, you must not use it.
    if (anyNA(theta.vectorized) || is.null(theta.vectorized))
      return (list(convergence = -3, value = -Inf))

    ObjectiveLogLikelihood <- function(theta.vectorized)
    {
      transition.probs <- matrix(theta.vectorized[1:(M*(M-1))],
                                ncol = (M-1), byrow = T)
      initial.dist <- c(theta.vectorized[initial.dist.index:(beta.index - 1)])
      # retrieve the original from the reduced form.
      transition.probs <- t(apply(transition.probs, 1,
                                  function (row) c(row, (1-sum(row)))))
      initial.dist <- c(initial.dist, (1-initial.dist))

      beta <- theta.vectorized[beta.index:(mu.index - 1)]
      if (!is.beta.switching) # make it as a switching parameter if not.
        beta <- rep(beta, M)
      beta <- matrix(beta, ncol = M)

      sigma <- theta.vectorized[sigma.index:(gamma.dep.index - 1)]
      if (!is.sigma.switching) # make it as a switching parameter if not.
        sigma <- rep(sigma, M)

      gamma.dependent <- t(rep(0,M))
      gamma.independent <- 0
      # i.e. gamma.dependent exists
      if (gamma.dep.index != (gamma.indep.index - M))
        gamma.dependent   <- matrix(theta.vectorized[gamma.dep.index:
                                                       (gamma.indep.index - 1)],
                                    ncol = M)
      # i.e. gamma.independent exists
      if (gamma.indep.index <= length(theta.vectorized))
        gamma.independent <- theta.vectorized[gamma.indep.index:
                                              length(theta.vectorized)]

      # slsqp solves a minimization problem;
      # take a negative value to turn the problem into max. problem
      -LikelihoodMSIAR(y, y.lagged, z.dependent, z.independent,
                      transition.probs,
                      initial.dist,  # initial.dist
                      beta = beta,  # beta
                      theta.vectorized[mu.index:(sigma.index - 1)],  # mu
                      sigma,    # sigma
                      gamma.dependent,
                      gamma.independent) # gamma.indep
    }

    # hard constraints to prevent values from bounding off
    transition.probs.lb <- rep(transition.probs.min, M*(M-1))
    transition.probs.ub <- rep((1-(M-1)*transition.probs.min), M*(M-1))
    theta.vectorized[1:(M*(M-1))] <- pmax(transition.probs.lb,
                                          theta.vectorized[1:(M*(M-1))])
    theta.vectorized[1:(M*(M-1))] <- pmin(transition.probs.ub,
                                          theta.vectorized[1:(M*(M-1))])
    initial.dist.lb <- rep(lb.prob.density, (M-1))
    initial.dist.ub <- rep(ub.prob.density, (M-1))
    theta.vectorized[(M*(M-1)+1):(M*M-1)] <- pmax(initial.dist.lb,
                                                  theta.vectorized[(M*(M-1)+1):
                                                  (M*M-1)])
    theta.vectorized[(M*(M-1)+1):(M*M-1)] <- pmin(initial.dist.ub,
                                                  theta.vectorized[(M*(M-1)+1):
                                                  (M*M-1)])
    beta <- theta.vectorized[beta.index:(mu.index - 1)]
    beta.lb <- pmin(-1, beta * (1 - 0.5 * sign(beta)))
    beta.ub <- pmax(1, beta * (1 + 0.5 * sign(beta)))
    mu <- theta.vectorized[mu.index:(sigma.index - 1)]
    mu.lb <- pmin(-1, mu * (1 - 0.8 * sign(mu)))
    mu.ub <- pmax(1, mu * (1 + 0.8 * sign(mu)))
    sigma <- theta.vectorized[sigma.index:(gamma.dep.index - 1)]
    sigma.lb <- sigma * 0.9
    sigma.ub <- sigma * 8

    # hard constraints on
    # transtion.probs (ub & lb), initial.dist (ub & lb), and sigma (lb)
    lb <- c(transition.probs.lb, initial.dist.lb,
            beta.lb,
            mu.lb,
            sigma.lb)
    # if gamma.dep & gamma.indep remain
    if (!length(theta.vectorized) == length(lb))
      lb <- c(lb, rep(-Inf, (length(theta.vectorized) - length(lb))))
    ub <- c(transition.probs.ub, initial.dist.ub,
            beta.ub,
            mu.ub,
            sigma.ub)
    ub <- c(ub, rep(Inf, (length(theta.vectorized) - length(ub))))

    # sanity check for derivatives.
    if (!LongStepSanityCheck(x0 = theta.vectorized, fn = ObjectiveLogLikelihood))
      return (list (convergence = -Inf, value = -Inf))

    result <- slsqp(theta.vectorized,
                    fn = ObjectiveLogLikelihood,
                    lower = lb, upper = ub, hin = ConstraintMCTransition,
                    control = list(maxeval = maxit, ftol_abs = epsilon))
    result$value <- -result$value # take negative back to make it actual
    return (result)
  }

  candidates.matrix <- sapply(candidates,
                              function (theta) ThetaToReducedColumn(theta))

  # perform non-linear optimization in each candidate
  long.results <- apply(candidates.matrix, 2,
                        SLSQPNonSwitchingAR,
                        y, y.lagged, z.dependent, z.independent)
  long.convergence <- unlist(lapply(long.results, "[[", "convergence"))
  long.likelihoods <- unlist(lapply(long.results, "[[", "value"))
  long.likelihoods[!is.finite(long.likelihoods)] <- -Inf # abnormal values
  long.likelihoods[long.convergence < 0] <- -Inf # non-convergence

  # extract the one that returns the best log.likelihood
  long.result <- long.results[[(which(long.likelihoods==
                                      max(long.likelihoods))[1])]]
  if (!is.finite(long.result$value))
    stop("Estimation failed. Try different settings for EM-algorithm.")
  return (list(theta = ReducedColumnToTheta(long.result$par),
               log.likelihood = long.result$value,
               long.results = long.results))
}

LongStepSanityCheck <- function (x0, fn)
{
  # sanity check for derivatives
  heps <- .Machine$double.eps^(1/3) # epsilon used in nloptr package
  n <- length(x0)
  hh <- diag(heps, n)
  gr <- numeric(n)
  for (i in 1:n)
    if (is.na(fn(x0 + hh[,i])) || is.na(fn(x0 - hh[,i])))
      return (FALSE)
  return (TRUE)

}
