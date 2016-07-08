#' Generates a sample of length n from theta based on observed initial values
#' and the state of the very last observation given. This applies only for
#' univariate time series.
#' @export
#' @title GenerateSample
#' @name GenerateSample
#' @param theta A list that represents the parameters of a model with items:
#' \item{beta}{s by 1 column for state-independent coefficients on AR(s)}
#' \item{mu}{M by 1 column that contains state-dependent mu}
#' \item{sigma}{M by 1 column that contains state-dependent sigma}
#' \item{gamma.dependent}{p_dep by M matrix that contains switching
#' coefficients for state-dependent exogenous variables}
#' \item{gamma.independent}{p_indep by 1 column that contains non-switching
#' coefficients for state-independent exogenous variables}
#' \item{transition.probs}{M by M matrix that contains transition probabilities}
#' \item{initial.dist}{M by 1 column that represents an initial distribution}
#' @param n The number of observations to be created.
#' @param initial.y.set n_initial by 1 column that represents previous samples;
#' n_initial must be larger than/equal to s, the number of autoregressive terms.
#' @param initial.state An integer in {1, 2, ..., M} that represents the state of
#' the very last observation of initial.y.set.
#' @param z_dependent n by p_dep matrix of data for exogenous variables that
#' have switching coefficeints.
#' @param z_independent n by p_indep matrix of data for exogenous variables that
#' have non-switching coefficeints.
GenerateSample <- function(theta, n = 100, initial.y.set, initial.state,
  z.dependent = NULL, z.independent = NULL)
{
  M <- nrow(theta$transition.probs)
  mu <- as.matrix(theta$mu)
  sigma <- as.matrix(theta$sigma)
  beta <- as.matrix(0)
  if (!is.null(theta$beta))
    beta <- as.matrix(theta$beta)
  s <- nrow(beta)
  gamma.dependent <- matrix(rep(0,M), ncol = M)
  gamma.independent <- as.matrix(0)

  # reformatting parameters & safety check
  if (ncol(beta) < 2) # even if beta is not switching make it a matrix
  {
    s <- length(beta)
    beta <- sapply(rep(0,M), function(col) beta)
    if (s == 1)
      beta <- t(as.matrix(beta)) # if s == 1, sapply results in a column.
  }
  if (nrow(beta) > length(initial.y.set))
    stop("EXCEPTION: the initial y set must have a length greater than equal to
        the number of regressive terms in the model.")
  if (length(mu) < 2) # even if mu is not switching make it a vector
    mu <- as.matrix(rep(mu[1,1], M))
  if (length(sigma) < 2) # even if sigma is not switching make it a vector
    mu <- as.matrix(rep(sigma[1,1], M))
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,length(initial.y.set) + n),ncol=1)
  else
  {
    z.dependent <- as.matrix(z.dependent)
    if (n != nrow(z.dependent))
      stop("EXCEPTION: the length of samples should match the number of
           observations for z.dependent.")
    z.dependent <- rbind(matrix(rep(Inf, s*ncol(z.dependent)), 
                                ncol = ncol(z.dependent)),
                         as.matrix(z.dependent))
    gamma.dependent <- as.matrix(theta$gamma.dependent)
  }
  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,length(initial.y.set) + n),ncol=1)
  else
  {
    z.independent <- as.matrix(z.independent)
    if (n != nrow(z.independent))
      stop("EXCEPTION: the length of samples should match the number of
        observations for z.independent.")
    z.independent <- rbind(matrix(rep(Inf, s*ncol(z.independent)), 
                                  ncol = ncol(z.independent)),
                           as.matrix(z.independent))
    gamma.independent <- as.matrix(theta$gamma.independent)
  }

  # initialization
  y <- c(initial.y.set, rep(-1, n))
  states <- c(rep(-1, (length(initial.y.set) - 1)),
              initial.state,
              rep(0, n))
  initial.index <- length(initial.y.set) + 1
  last.index <- length(initial.y.set) + n
  for (k in initial.index:last.index)
  {
    previous.state <- states[(k-1)]
    trans.cumsum <- cumsum(theta$transition.probs[previous.state,])
    prob <- runif(1) # decision to switch
    state = 1
    for (j in 2:M)
      if (prob > trans.cumsum[j-1] && prob <= trans.cumsum[j]) {
        state <- j
        break
      }
    states[k] <- state
    y[k] <- mu[state,1] +
      t(y[(k-s):(k-1)]) %*% as.numeric(beta[,state]) +
      z.dependent[k,] %*% as.matrix(gamma.dependent[,state]) +
      z.independent[k,] %*% as.matrix(gamma.independent) +
      rnorm(1,sd=sigma[state,1])
  }

  return (list(sample = y[initial.index:length(y)],
  states = states[initial.index:length(states)]))
}
