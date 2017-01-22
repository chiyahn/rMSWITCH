#' Estimates MLE of a Markov regime-switching autoregressive (MS-AR) model on
#' AR(s) with M states from sample data
#' @export
#' @title EstimateMSAR
#' @name EstimateMSAR
#' @param y n by 1 vector of data for y
#' @param z.dependent n by p.dep matrix of data for switching exogenous
#' variables
#' @param z.independent n by p.indep matrix of data for non-switching exogenous
#' variables
#' @param M The number of states in the model
#' @param s The number of terms used for AR(s)
#' @param is.beta.switching Specifies whether autoregressive terms are switching
#' @param is.sigma.swithcing Specifies whether the model is heteroscedastic.
#' @param is.MSM Specifies whether the model is switching in mean (MSM) or
#' intercept (MSI).
#' @param initial.theta An initial guess for MS-AR model
#' @param epsilon Epsilon used as convergence criterion.
#' @param maxit The maximum number of iterations.
#' @param short.n Number of initial draws for EM estimation
#' @param short.iterations Maximum iteration used to perform short EMs
#' @param transition.probs.min Minimum set for transition prob. matrix
#' @param sigma.min Minimum set for variance.
#' @param nloptr Determines whether nonlinear optimization package is used.
#' @param estimate.fisher Determines whether the variance of each estimate is going
#' to be computed.
#' @return  A list with items:
#' \item{beta}{s by 1 column for state-independent coefficients on AR(s)}
#' \item{mu}{M by 1 column that contains state-dependent mu}
#' \item{sigma}{M by 1 column that contains state-dependent sigma}
#' \item{gamma.dependent}{p_dep by M matrix that contains switching
#' coefficients for state-dependent exogenous variables}
#' \item{gamma.independent}{p_indep by 1 column that contains non-switching
#' coefficients for state-independent exogenous variables}
#' \item{transition.probs}{M by M matrix that contains transition probabilities}
#' \item{initial.dist}{M by 1 column that represents an initial distribution}
#' \item{nlotpr}{Determines whether nonlinear optimization package is used
#' in estimation.}
#' \item{estimate.fisher}{Determines whether the variance of each estimate is going
#' to be computed.}
#' @examples
#' theta <- RandomTheta()
#' y <- GenerateSample(theta = theta)$y
#' EstimateMSAR(y = y, M = 2, s = 1,
#'              is.beta.switching = FALSE,
#'              is.sigma.switching = TRUE)
EstimateMSAR <- function(y = y, z.dependent = NULL, z.independent = NULL,
                        M = 2, s = 2,
                        is.beta.switching = FALSE,
                        is.sigma.switching = TRUE,
                        is.MSM = FALSE,
                        initial.theta = NULL,
                        epsilon = 1e-08, maxit = 2000,
                        short.n = 10, short.epsilon = 1e-03,
                        short.iterations = 200,
                        transition.probs.min = 0.01,
                        sigma.min = 0.02,
                        nloptr = FALSE,
                        estimate.fisher = TRUE) {
  if (test.on) # initial values controlled by test.on
    set.seed(test.seed)

  p.dependent <- 0
  if (s + 1 > length(y))
  {
    print ("EXCEPTION: The length of observations must be greater than s.")
    return (NULL)
  }
  if (is.MSM)
    stop ("MSM models are currently not supported.")

  # formatting dataset
  transition.probs.max <- 1 - (M-1)*transition.probs.min
  
  lagged.and.sample <- GetLaggedAndSample(y, s)
  y.lagged <- lagged.and.sample$y.lagged
  y.sample <- lagged.and.sample$y.sample
  n <- length(y.sample)
  initial.params <- NULL

  # remove first s rows of z.dependent and z.independent
  if (!is.null(z.dependent))
    z.dependent <- as.matrix(as.matrix(z.dependent[(s+1):length(y),]))
  if (!is.null(z.independent))
    z.independent <- as.matrix(as.matrix(z.independent[(s+1):length(y),]))

  # 1. Get initial parameter using regmix if initial.theta is not given
  if (is.null(initial.theta) || M == 1)
  {
    initial.params <- GetInitialParams(y.sample, y.lagged,
                                       z.dependent, z.independent,
                                       M = M, s = s, p.dependent = p.dependent,
                                       is.beta.switching = is.beta.switching,
                                       is.sigma.switching = is.sigma.switching)
    initial.theta <- initial.params$theta
  }
  else
  {
    initial.theta$beta <- as.matrix(initial.theta$beta)
    if (!is.null(initial.theta$gamma.dependent))
      initial.theta$gamma.dependent <- as.matrix(initial.theta$gamma.dependent)
  }

  if (M > 1)
  {
    # 2. Run short EM
    # how many candidates would you like to find?
    short.n.candidates <- max(floor(sqrt(log(n)*(1+s)*M)*2*short.n), 400)
    short.thetas <- lapply(1:short.n.candidates,
                          function(j) 
                            EstimateMSARInitShort(theta = initial.theta,
                                                  transition.probs.min =
                                                    transition.probs.min,
                                                  transition.probs.max =
                                                    transition.probs.max))
    # For compatibility with cpp codes, change gammas to
    # appropriate zero vectors. After computation, they will be returned NULL.
    if (is.null(z.dependent))
      initial.theta$gamma.dependent <- matrix(rep(0,M), ncol = M)
    if (is.null(z.independent))
      initial.theta$gamma.independent <- as.matrix(0)
    # include the original theta
    short.thetas[[length(short.thetas) + 1]] <- initial.theta
    short.results <- MaximizeShortStep(short.thetas = short.thetas,
                          y = y.sample, y.lagged = y.lagged,
                          z.dependent = z.dependent,
                          z.independent = z.independent,
                          is.beta.switching = is.beta.switching,
                          is.sigma.switching = is.sigma.switching,
                          maxit = short.iterations, epsilon = short.epsilon,
                          transition.probs.min = transition.probs.min,
                          transition.probs.max = transition.probs.max,
                          sigma.min = sigma.min)

    short.likelihoods <- sapply(short.results, "[[", "likelihood")
    
    # 3. Run long step
    long.thetas <- lapply(short.results, "[[", "theta")
    long.thetas <- long.thetas[order(short.likelihoods,decreasing=T)[1:
                                           min(length(long.thetas), short.n)]] 

    if (nloptr)
      long.result <- MaximizeLongStepNLOPTR(long.thetas,
                                      y = y.sample, y.lagged = y.lagged,
                                      z.dependent = z.dependent,
                                      z.independent = z.independent,
                                      is.beta.switching = is.beta.switching,
                                      is.sigma.switching = is.sigma.switching,
                                      epsilon = epsilon, maxit = maxit,
                                      transition.probs.min = transition.probs.min,
                                      transition.probs.max = transition.probs.max,
                                      sigma.min = sigma.min)    
    else
      long.result <- MaximizeLongStep(long.thetas,
                                      y = y.sample, y.lagged = y.lagged,
                                      z.dependent = z.dependent,
                                      z.independent = z.independent,
                                      is.beta.switching = is.beta.switching,
                                      is.sigma.switching = is.sigma.switching,
                                      epsilon = epsilon, maxit = maxit,
                                      transition.probs.min = transition.probs.min,
                                      transition.probs.max = transition.probs.max,
                                      sigma.min = sigma.min)
    if (!long.result$succeeded)
    {
      print("Estimation failed. Try different settings for EM-algorithm.")
      return (NULL)
    }
  }
  else
  {
    long.result <- list(log.likelihood = initial.params$log.likelihood,
                        theta = initial.params$theta)
  }

  # 4. Final formatting
  theta <- long.result$theta
  theta$initial.dist <- theta$initial.dist / sum(theta$initial.dist)
  if (is.null(z.dependent))
    theta$gamma.dependent <- NULL
  if (is.null(z.independent))
    theta$gamma.independent <- NULL

  # 4.1. Sort them based on mu
  if (M > 1)
  {
    mu.order  <- order(theta$mu)
    theta$transition.probs <- OrderTransitionMatrix(theta$transition.probs,
                                                    mu.order)
    theta$initial.dist <- theta$initial.dist[mu.order]
    if (is.beta.switching)
      theta$beta <- matrix(theta$beta[,mu.order], ncol = M)
    theta$mu        <- theta$mu[mu.order]
    if (is.sigma.switching)
      theta$sigma     <- theta$sigma[mu.order]
    if (!is.null(theta$gamma.dependent))
      theta$gamma.dependent <- theta$gamma.dependent[,mu.order]
  }

  # 4.2. Basic information
  log.likelihood <- long.result$log.likelihood

  # 4-3. Based on formatted dataset, get posterior probabilities
  posterior.probs <- EstimatePosteriorProbs(theta = theta,
                      y = y.sample, y.lagged = y.lagged,
                      z.dependent = z.dependent, z.independent = z.independent)
  states <- EstimateStates(posterior.probs$xi.n) # use smoothed probabilities
  fisher.estimated <- NULL
  
  if (estimate.fisher)
    fisher.estimated <- EstimateFisherInformation(theta = theta,
                        y = y.sample, y.lagged = y.lagged,
                        z.dependent = z.dependent, z.independent = z.independent)

  msar.model <- list(theta = theta,
                     log.likelihood = log.likelihood,
                     posterior.probs.filtered = posterior.probs$xi.k,
                     posterior.probs.smoothed = posterior.probs$xi.n,
                     fisher.estimated = fisher.estimated,
                     states = states,
                     call = match.call(),
                     is.MSM = is.MSM,
                     label = "msar.model")
  class(msar.model) <- "msar.model"

  return (msar.model)
}

# Returns a list of lagged y (y.lagged) and corresponding y sample (y.sample)
GetLaggedAndSample <- function(y, s)
{
  y <- as.numeric(y)
  y.lagged <- sapply(seq(0,s), GetLaggedColumn, y, s) # (n-s) by s matrix
  y.sample <- as.matrix(y.lagged[,1])
  y.lagged <- as.matrix(y.lagged[,-1])
  return (list (y.lagged = y.lagged, y.sample = y.sample))
}

# Get initial theta to run EM algorithm, using normalregMix package.
GetInitialParams <- function (y.sample, y.lagged, z.dependent, z.independent,
                              M, s, p.dependent,
                              is.beta.switching, is.sigma.switching)
{
  regmix.result <- regmixPMLE(y = y.sample, x = cbind(y.lagged, z.dependent),
                              z = z.independent, m = M, vcov.method="none")
  regmix.theta <- regmix.result$parlist
  regmix.transition.probs <- StatesToTransitionProbs(states =
                                                regmix.result$components, M = M)
  regmix.beta <- mean(regmix.theta$mubeta[2:(s+1),])
  regmix.sigma <- 1
  regmix.gamma.dependent <- NULL
  regmix.gamma.independent <- regmix.theta$gamma

  regmix.initial.dist <- StatesToInitialDist(states = regmix.result$components,
                                            M = M)

  if (is.beta.switching) # estimate for state-dependent. beta
    regmix.beta <- matrix(regmix.theta$mubeta[2:(s+1),], ncol = M)
  else if (s > 1)# estimate for state-indep. beta
    regmix.beta <- apply(as.matrix(regmix.theta$mubeta[2:(s+1),]), 1, mean)
  if (is.sigma.switching)
    regmix.sigma <- regmix.theta$sigma
  else
    regmix.sigma <- sqrt(sum(regmix.theta$alpha * regmix.theta$alpha *
                            regmix.theta$sigma * regmix.theta$sigma))
  if (p.dependent > 0)
  {
    # estimate for state-dep. gamma
    regmix.gamma.dependent <- regmix.theta$mubeta[(s+2):(s+1+p.dependent),]
    regmix.gamma.dependent <- as.matrix(regmix.gamma.dependent)
    # if is one-dim, the extracted is a seq (becomes a col); take a transpose.
    if (p.dependent == 1)
      regmix.gamma.dependent <- t(regmix.gamma.dependent)
  }
  if (!is.null(regmix.gamma.independent))
    regmix.gamma.independent <- as.matrix(regmix.gamma.independent)

  regmix.theta <- list(beta = as.matrix(regmix.beta),
                       mu = regmix.theta$mubeta[1,],
                       sigma = regmix.sigma,
                       gamma.dependent = regmix.gamma.dependent,
                       gamma.independent = regmix.gamma.independent,
                       transition.probs = regmix.transition.probs,
                       initial.dist = regmix.initial.dist)

  return(list(log.likelihood = regmix.result$loglik, theta = regmix.theta))
}

# Get a column of 0 \leq j \leq s lagged variable
GetLaggedColumn <- function (j, col, s) {
  if (j != s)
    col <- col[-(0:(s-j))] # destroy first s-j elements
  return (col[1:(length(col)-j)])
}

# Gives variation in theta given
EstimateMSARInitShort <- function(theta, 
                                  transition.probs.min,
                                  transition.probs.max) {
  beta0 <- theta$beta
  mu0 <- theta$mu
  sigma0 <- theta$sigma
  gamma.dependent0 <- theta$gamma.dependent
  gamma.independent0 <- theta$gamma.independent
  transition.probs0 = theta$transition.probs
  initial.dist0 = theta$initial.dist
  M <- ncol(transition.probs0)

  sigma.epsilon <- 0.6 # minimum sigma

  transition.probs <- matrix(0, ncol = M, nrow = M)
  for (i in 1:M)
  {
    transition.probs[i,i] <- runif(1, (transition.probs.max * 0.7), transition.probs.max)
    transition.probs[i,][-i] <- runif((M-1), transition.probs.min, transition.probs[i,i])  
    transition.probs[i,] <- transition.probs[i,] / sum(transition.probs[i,])
  }

  initial.dist = VariationInRow(initial.dist0, 
                                transition.probs.min,
                                transition.probs.max)
  beta <- beta0 # beta estimate should be accurate enough
  mu <- mu0 + rnorm(length(mu0))
  sigma <- sapply(sigma0, function(sig) max(sig + rnorm(1), sigma.epsilon))

  # gamma estimates should be accurate enough
  gamma.dependent <- gamma.dependent0
  gamma.independent <- gamma.independent0

  # For compatibility with cpp codes, change gammas to
  # appropriate zero vectors. After computation, they will be returned NULL.
  if (!is.null(gamma.dependent))
    gamma.dependent = as.matrix(gamma.dependent)
  else
    gamma.dependent = matrix(rep(0, M), ncol = M)
  if (!is.null(gamma.independent))
    gamma.independent = as.matrix(gamma.independent)
  else
    gamma.independent = as.matrix(0)

  theta <- list(beta = as.matrix(beta),
                mu = mu,
                sigma = sigma,
                gamma.dependent = gamma.dependent,
                gamma.independent = gamma.independent,
                transition.probs = transition.probs,
                initial.dist = initial.dist)
  return (theta)
}

# Produces a rough estimate for transition matrix from estimated states
StatesToTransitionProbs <- function(states, M)
{
  transition.probs <- matrix(rep(0,M*M), ncol = M)
  n <- length(states) - 1
  for (i in 1:M)
    for (j in 1:M)
      for (k in 1:n) # n has been already subtracted by one
        if (states[k] == i && states[(k+1)] == j)
          transition.probs[i,j] <- 1 + transition.probs[i,j]
  transition.probs <- t(apply(transition.probs, 1,
                              function(row)
                              {
                                if (sum(row) != 0)
                                  return (row / sum(row))
                                else # state does not appear/appears at last
                                  return (rep(1/M, M))
                              }))
  return (transition.probs)
}

# Produces a rough estimate for initial distribution from estimated states
StatesToInitialDist <- function(states, M)
{
  initial.dist <- seq(1,M)
  n <- length(states)
  initial.dist <- sapply(initial.dist, function(i) sum(states == i) / n)
  return (initial.dist)
}

# Gives a random variation in a row of a transition probability matrix
VariationInRow <- function(row, transition.probs.min,
                           transition.probs.max)
{
  row <- runif(length(row),transition.probs.min,
               transition.probs.max)
  return (row / sum(row))
}

EstimateFisherInformation <- function(theta, y, y.lagged,
                                      z.dependent, z.independent,
                                      eps = 1e-6)
{
  # use the first candidate to save the information about dimensions
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

  LogLikelihoods <- function(theta.vectorized)
  {
    transition.probs <- as.matrix(1)
    if (M > 1)
    {
      transition.probs <- matrix(theta.vectorized[1:(M*(M-1))],
                                ncol = (M-1), byrow = T)
      transition.probs <- t(apply(transition.probs, 1,
                                  function (row) c(row, (1-sum(row)))))
    }
    initial.dist <- c(theta.vectorized[initial.dist.index:(beta.index - 1)])
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

    LikelihoodsMSIAR(y, y.lagged, z.dependent, z.independent,
                    transition.probs,
                    initial.dist,  # initial.dist
                    beta = beta,  # beta
                    theta.vectorized[mu.index:(sigma.index - 1)],  # mu
                    sigma,    # sigma
                    gamma.dependent,
                    gamma.independent) # gamma.indep
  }

  x <- ThetaToReducedColumn(theta)
  # Defines a step (make sure it does not bind with the ub/lb)
  h <- pmax(eps, abs(x)) * eps ^ {2/3}
  xh <- x + h
  h <- xh - x
  h.diag <- diag(h)

  G <- matrix(0, nrow = length(x), ncol = n)
  H <- matrix(0, nrow = length(x), ncol = length(x))

  for (i in 1:length(x))
    G[i,] <- (LogLikelihoods(x + h.diag[i,]) - LogLikelihoods(x - h.diag[i,])) /
              (2 * h[i])

  for (k in 1:n)
    H <- H + G[,k] %*% t(G[,k])

  return (H / n)
}
