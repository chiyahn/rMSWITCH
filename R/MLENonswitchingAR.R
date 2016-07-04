#' Estimates MLE of Markov Regime Switching (MRS) Model on AR(s) with M states from sample data
#' where regressors on the autoregression are independent of states
#' @export
#' @title MLENonswitchingAR
#' @name MLENonswitchingAR
#' @param y n by 1 vector of data for y
#' @param z n by p matrix of data for z (exogenous variables)
#' @param z.is.switching p by 1 vector of booleans that indicate
#' whether jth column (1 \leq j \leq p) of z has a corresponding
#' coefficient that is state-dependent
#' @param M The number of states in the model
#' @param s The number of terms used for AR(s)
#' @param phi.initial The mrs model
#' @param epsilon Epsilon used as convergence criterion.
#' @param maxit The maximum number of iterations.
#' @param short.n Number of short EMs
#' @param short.iterations Maximum iteration used to perform short EMs
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
MLENonswitchingAR <- function(y = y, z = NULL, z.is.switching = NULL, M = 3, s = 2, theta.initial = NULL,
                        epsilon = 1e-08, maxit = 2000, short.n = 100, short.epsilon = 1e-02,
                        short.iterations = 30) {
  # TODO: change name to MLENonswitchingAR
  p.dependent <- 0
  if (s + 1 > length(y))
  {
    print ("EXCEPTION: The length of observations must be greater than s.")
    return (NULL)
  }
  if (!is.null(z) && ncol(z) != length(z.is.switching))
  {
    print ("EXCEPTION: You must specify which terms of coefficients for z are switching.")
    return (NULL)
  }

  short.n.candidates <- max(short.n*5*(s+length(z.is.switching))*M, 100)

  # formatting dataset
  y <- as.numeric(y)
  y.lagged <- sapply(seq(s,0), GetLaggedColumn, y, s) # (n-s) by s matrix
  y.sample <- as.matrix(y.lagged[,1])
  y.lagged <- as.matrix(y.lagged[,-1])
  n <- length(y.sample)
  z.dependent <- NULL
  z.independent <- NULL

  # divide z into two parts: z that is state-dependent and independent.
  if (!is.null(z))
  {
    if (is.null(z.is.switching[1]))
    {
      print ("WARNING: You must specify which terms of coefficients for z are switching.")
      print ("By default, all coefficients are going to be assumed to be dependent on states.")
      z.is.switching <- rep(TRUE, ncol(z))
    }
    z.lagged <- apply(z, 2, GetLaggedColumn, j = s, s = s) # remove the first s terms
    z.dependent <- as.matrix(z.lagged[,z.is.switching])
    z.independent <- as.matrix(z.lagged[,!z.is.switching])
    p.dependent <- ncol(z.dependent)
    # if one is a column of length 0, transform it into just NULL.
    if (length(z.dependent) == 0)
      z.dependent <- NULL
    if (length(z.independent) == 0)
      z.independent <- NULL
  }



  # 1. Get initial parameter using regmix if theta.initial is not given
  if (is.null(theta.initial))
    theta.initial <- GetInitialTheta(y.sample, y.lagged, z.dependent, z.independent, s, p.dependent, M)
  if (is.null(z.dependent))
  {
    z.dependent <- as.matrix(rep(0,n))
    theta.initial$gamma.dependent <- t(as.matrix(rep(0,M)))
  }
  if (is.null(z.independent))
  {
    z.independent <- as.matrix(rep(0,n))
    theta.initial$gamma.independent <- as.matrix(0)
  }

  ## TODO: On short & run EM, if dim(z) > 0, use ExpectationMaximizationIndepExo instead.
  # 2. Run short EM
  short.thetas <- lapply(1:short.n.candidates, function(j) MLENonswitchingARInitShort(theta.initial))
  short.thetas[[length(short.thetas) + 1]] <- theta.initial # include the original theta


  short.results <- lapply(short.thetas, ExpectationMaximizationIndep,
                          y = y.sample, y.lagged = y.lagged, z_dependent = z.dependent, z_independent = z.independent,
                          maxit = short.iterations, epsilon = short.epsilon)
  short.likelihoods <- sapply(short.results, "[[", "likelihood")

  # 3. Run long EM
  long.thetas <- short.thetas[order(short.likelihoods,decreasing=T)[1:short.n]] # pick best short.n thetas
  long.results <- lapply(long.thetas, ExpectationMaximizationIndep,
                         y = y.sample, y.lagged = y.lagged, z_dependent = z.dependent, z_independent = z.independent,
                         maxit = maxit, epsilon = epsilon)
  long.likelihoods <- sapply(long.results, "[[", "likelihood")
  em.result <- long.results[[(which(long.likelihoods==max(long.likelihoods))[1])]] # gives you the best result

  if (is.null(z.is.switching))
  {
    em.result$theta$gamma.dependent <- NULL
    em.result$theta$gamma.independent <- NULL
  }
  else if (!is.element(TRUE, z.is.switching)) # i.e. none of z is switching
    em.result$theta$gamma.dependent <- NULL
  else if (!is.element(FALSE, z.is.switching)) # i.e. all z values are switching
    em.result$theta$gamma.independent <- NULL


  return (em.result)
}

ExpectationMaximizationIndep <- function (y, y.lagged,
                                          z_dependent, z_independent, theta,
                                          maxit, epsilon)
{
  EMIndepCPP(y, y.lagged, z_dependent, z_independent,
        theta$beta, theta$mu, theta$sigma,
        theta$gamma.dependent, theta$gamma.independent,
        theta$transition.probs, theta$initial.dist,
        maxit, epsilon)
}

# Get initial theta to run EM algorithm, using normalregMix package.
GetInitialTheta <- function (y.sample, y.lagged, z.dependent, z.independent, s, p.dependent, M)
{
  regmix.result <- regmixPMLE(y = y.sample, x = cbind(y.lagged, z.dependent), z = z.independent, m = M, vcov.method="OPG")
  regmix.theta <- regmix.result$parlist
  regmix.transition.probs <- StatesToTransitionProbs(states = regmix.result$indices, M = M)
  regmix.beta <- mean(regmix.theta$mubeta[2,]) # WATCH: state-independent beta case.
  regmix.gamma.dependent <- NULL
  regmix.gamma.independent <- regmix.theta$gamma

  regmix.initial.dist <- StatesToInitialDist(states = regmix.result$indices, M = M)

  if (s > 1) # WATCH: state-independent beta case.
    regmix.beta <- apply(regmix.theta$mubeta[2:(s+1),], 1, mean) # estimate for state-ind. beta
  if (p.dependent > 0)
  {
    regmix.gamma.dependent <- regmix.theta$mubeta[(s+2):(s+1+p.dependent),] # estimate for state-dep. gamma
    regmix.gamma.dependent <- as.matrix(regmix.gamma.dependent)
    if (p.dependent == 1) # if is one-dim, the extracted is a seq (becomes a col), so take a transpose.
      regmix.gamma.dependent <- t(regmix.gamma.dependent)
  }
  if (!is.null(regmix.gamma.independent))
    regmix.gamma.independent <- as.matrix(regmix.gamma.independent)

  regmix.theta <- list(beta = as.matrix(regmix.beta), # WATCH: you might want to take a transpose for s=1 if state-dependent
                       mu = regmix.theta$mubeta[1,],
                       sigma = regmix.theta$sigma,
                       gamma.dependent = regmix.gamma.dependent,
                       gamma.independent = regmix.gamma.independent,
                       transition.probs = regmix.transition.probs,
                       initial.dist = regmix.initial.dist)

  return(regmix.theta)
}

# Get a column of 0 \leq j \leq s lagged variable
GetLaggedColumn <- function (j, y, s) {
  if (j == 0)
    return (y[(s+1):(length(y) - s)])
  col <- y[-(0:(s+j))] # destroy first s+j elements
  col <- col[1:(length(col) - s + j)]
  return (col)
}

# Gives variation in theta given
MLENonswitchingARInitShort <- function(theta) {
  beta0 <- theta$beta
  mu0 <- theta$mu
  sigma0 <- theta$sigma
  gamma.dependent0 <- theta$gamma.dependent
  gamma.independent0 <- theta$gamma.independent
  transition.probs0 = theta$transition.probs
  initial.dist0 = theta$initial.dist

  sigma.epsilon <- 0.2 # minimum sigma

  beta <- beta0 # beta estimate should be accurate enough
  gamma.dependent <- gamma.dependent0 # gamma estimate should be accurate enough
  gamma.independent <- gamma.independent0 # gamma estimate should be accurate enough
  mu <- mu0 + rnorm(length(mu0))
  sigma <- sapply(sigma0, function(sig) max(sig + rnorm(1), sigma.epsilon))
  transition.probs = t(apply(transition.probs0, 1, VariationInRow))
  initial.dist = VariationInRow(initial.dist0)

  theta <- list(beta = as.matrix(beta),
                mu = mu,
                sigma = sigma,
                gamma.dependent = as.matrix(gamma.dependent),
                gamma.independent = as.matrix(gamma.independent),
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
  transition.probs <- transition.probs / n
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
VariationInRow <- function(row)
{
  variations <- runif(length(row), min = -min(row), max = 1 - max(row))
  row <- sapply((row + variations), function(i) min(0.98, max(0.02, i)))
  return (row / sum(row))
}

# Determines which state each observation belongs to based on theta
EstimateStatesIndep <- function(y, y.lagged, theta) {
  nu <- NuIndep(y, y.lagged, theta$beta, theta$mu, theta$sigma)
  return (apply(nu, 1, function(i) (which(i==max(i)))))
}

# MLENonswitchingAR using EM algorithm written in R
MLENonswitchingARR <- function(y = y, z = NULL, z.is.switching = FALSE, M = 3, s = 2, theta.initial = NULL,
                        epsilon = 1e-08, maxit = 60, short.n = 200, short.iterations = 10) {

  p.dependent <- 0
  if (s + 1 > length(y))
  {
    print ("EXCEPTION: The length of observations must be greater than s.")
    return (NULL)
  }
  if (!is.null(z) && ncol(z) != length(z.is.switching))
  {
    print ("EXCEPTION: You must specify which terms of coefficients for z are switching.")
    return (NULL)
  }


  short.n <- max(short.n, 20)

  # formatting dataset
  y <- as.numeric(y)
  y.lagged <- sapply(seq(s,0), GetLaggedColumn, y, s) # (n-s) by s matrix
  y.sample <- y.lagged[,1]
  y.lagged <- as.matrix(y.lagged[,-1])
  n <- length(y.sample)
  z.dependent <- NULL
  z.independent <- NULL

  # divide z into two parts: z that is state-dependent and independent.
  if (!is.null(z))
  {
    if (is.null(z.is.switching[1]))
    {
      print ("WARNING: You must specify which terms of coefficients for z are switching.")
      print ("By default, all coefficients are going to be assumed to be dependent on states.")
      z.is.switching <- rep(TRUE, ncol(z))
    }
    z.lagged <- apply(z, 2, GetLaggedColumn, j = s, s = s) # remove the first s terms
    z.dependent <- as.matrix(z.lagged[,z.is.switching])
    z.independent <- as.matrix(z.lagged[,!z.is.switching])
    p.dependent <- ncol(z.dependent)
    # if one is a column of length 0, transform it into just NULL.
    if (length(z.dependent) == 0)
      z.dependent <- NULL
    if (length(z.independent) == 0)
      z.independent <- NULL
  }



  # 1. Get initial parameter using regmix if theta.initial is not given
  if (is.null(theta.initial))
    theta.initial <- GetInitialTheta(y.sample, y.lagged,
                                     z.dependent, z.independent, s, p.dependent, M)
  if (is.null(z.dependent))
  {
    z.dependent <- as.matrix(rep(0,n))
    theta.initial$gamma.dependent <- t(as.matrix(rep(0,M)))
  }
  if (is.null(z.independent))
  {
    z.independent <- as.matrix(rep(0,n))
    theta.initial$gamma.independent <- as.matrix(0)
  }

  ## TODO: On short & run EM, if dim(z) > 0, use ExpectationMaximizationIndepExo instead.
  # 2. Run short EM
  short.thetas <- lapply(1:short.n, function(j) MLENonswitchingARInitShort(theta.initial))
  short.thetas[[length(short.thetas) + 1]] <- theta.initial # include the original theta


  short.results <- lapply(short.thetas, ExpectationMaximizationIndepR,
                          y = y.sample, y.lagged = y.lagged,
                          z_dependent = z.dependent, z_independent = z.independent,
                          maxit = short.iterations, epsilon = epsilon)
  short.likelihoods <- sapply(short.results, "[[", "likelihood")

  # 3. Run long EM
  long.thetas <- short.thetas[order(short.likelihoods,decreasing=T)[1:20]] # pick best 20 thetas
  long.results <- lapply(long.thetas, ExpectationMaximizationIndepR,
                         y = y.sample, y.lagged = y.lagged,
                         z_dependent = z.dependent, z_independent = z.independent,
                         maxit = maxit, epsilon = epsilon)
  long.likelihoods <- sapply(long.results, "[[", "likelihood")
  em.result <- long.results[[(which(long.likelihoods==max(long.likelihoods))[1])]] # gives you the best result

  return (em.result)
}
