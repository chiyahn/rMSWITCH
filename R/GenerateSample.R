#' Generates a sample of length n from theta based on observed initial values
#' and the state of the very last observation given. This applies only for
#' univariate time series.
#' @export
#' @title GenerateSample
#' @name GenerateSample
#' @param theta A list that represents the parameters of a model with items:
#' \item{transition.probs}{M by M matrix that contains transition probabilities}
#' \item{initial.dist}{M by 1 column that represents an initial distribution}
#' \item{beta}{s by 1 column for state-independent coefficients on AR(s)}
#' \item{mu}{M by 1 column that contains state-dependent mu}
#' \item{sigma}{M by 1 column that contains state-dependent sigma}
#' \item{gamma.dependent}{p_dep by M matrix that contains switching
#' coefficients for state-dependent exogenous variables}
#' \item{gamma.independent}{p_indep by 1 column that contains non-switching
#' coefficients for state-independent exogenous variables}
#' @param n The number of sample observations to be created.
#' @param initial.y.set n_initial by 1 column that represents previous samples;
#' n_initial must be larger than/equal to s, the number of autoregressive terms.
#' By default, initial.y.set is going to be determined by rnorm(s).
#' @param initial.state An integer in {1, 2, ..., M} that represents the state of
#' the very last observation of initial.y.set. The default value is 1.
#' @param z_dependent n by p_dep matrix of data for exogenous variables that
#' have switching coefficeints.
#' @param z_independent n by p_indep matrix of data for exogenous variables that
#' have non-switching coefficeints.
#' @param is.MSM Determines whether the model follows MSM-AR. If it is set to be
#' TRUE, the model is assumed to be MSM-AR. 
#' @return  A list with items:
#' \item{y}{(n + length(initial.y.set)) by 1 column that represents a sample
#' appended with previous values used to estimate autoregressive terms}
#' \item{y.sample}{n by 1 column that represents a sample of the model}
#' \item{y.lagged}{n by s matrix that represents a lagged sample of the model,
#' where kth column represents a kth lagged column}
#' \item{states}{M by 1 column that contains state-dependent sigma}
#' \item{msar.model}{An instance in msar.model that represents the actual model
#' used to create the sample.}
#' @examples
#' GenerateSample()
#' theta <- RandomTheta(M = 2, s = 3)
#' GenerateSample(theta)
#' GenerateSample(theta, n = 200)
GenerateSample <- function(theta = NULL, n = 100,
                           initial.y.set = NULL, initial.state = 1,
                           z.dependent = NULL, z.independent = NULL,
                           is.MSM = FALSE)
{
  if (is.null(theta))
    theta <- RandomTheta()
  M <- nrow(theta$transition.probs)
  mu <- as.matrix(theta$mu)
  sigma <- as.matrix(theta$sigma)
  beta <- as.matrix(0)
  if (!is.null(theta$beta))
    beta <- as.matrix(theta$beta)
  s <- nrow(beta)
  gamma.dependent <- matrix(rep(0,M), ncol = M)
  gamma.independent <- as.matrix(0)
  is.beta.switching = (ncol(as.matrix(beta)) > 1)
  is.sigma.switching = (length(sigma) > 1)
  
  if (is.null(initial.y.set))
    initial.y.set <- rnorm(s)
  
  if (length(initial.y.set) < s)
    stop ("EXCEPTION: The length of initial.y.set cannot be smaller than s.")
  if (length(initial.y.set) > s)
    warning ("The length of initial.y.set is greater than s;
             only the last s observations are going to be used for
             sample generation.")
  
  # reformatting parameters & safety check
  if (ncol(beta) < M) # even if beta is not switching make it a matrix
  {
    s <- length(beta)
    beta <- sapply(rep(0,M), function(col) beta)
    if (s == 1)
      beta <- t(as.matrix(beta)) # if s == 1, sapply results in a column.
  }
  if (length(sigma) < M) # even if sigma is not switching make it a matrix
    sigma <- matrix(rep(theta$sigma, M), nrow = M)
  
  if (nrow(beta) > length(initial.y.set))
    stop("EXCEPTION: the initial y set must have a length greater than equal to
         the number of regressive terms in the model.")
  if (length(mu) < 2) # even if mu is not switching make it a vector
    mu <- matrix(rep(mu[1,1], M), nrow = M)
  
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,(s + n)),ncol=1)
  else
  {
    z.dependent <- as.matrix(z.dependent)
    if (n > nrow(z.dependent))
      stop("EXCEPTION: the number of observations for z.dependent cannot be
           smaller than the length of samples, n.")
    z.dependent <- rbind(matrix(rep(NaN, (s*ncol(z.dependent))),
                                ncol = ncol(z.dependent)),
                         as.matrix(z.dependent[(nrow(z.dependent) - n + 1):
                                                 (nrow(z.dependent)),]))
    gamma.dependent <- as.matrix(theta$gamma.dependent)
  }
  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,(s + n)),ncol=1)
  else
  {
    z.independent <- as.matrix(z.independent)
    if (n > nrow(z.independent))
      stop("EXCEPTION: the number of observations for z.independent cannot be
           smaller than the length of samples, n.")
    z.independent <- rbind(matrix(rep(NaN, (s*ncol(z.independent))),
                                  ncol = ncol(z.independent)),
                           as.matrix(z.independent[(nrow(z.independent) - n + 1):
                                                     (nrow(z.independent)),]))
    gamma.independent <- as.matrix(theta$gamma.independent)
  }
  
  # initialization
  initial.y.set <- as.numeric(initial.y.set)
  initial.y.set <- initial.y.set[(length(initial.y.set) - s + 1):
                                   length(initial.y.set)] # only last s
  y <- c(initial.y.set, rep(-Inf, n))
  states <- vector()
  if (is.MSM)
  {
    if (length(initial.state) == (s + 1))
      states <- c(initial.state,
                  rep(0, n))
    else
    {
      # determine initial.state by initial.dist
      initial.dist.cumsum <- cumsum(theta$initial.dist)
      if (length(theta$initial.dist) == M) # if in short form, uniformly draw.
      {
        states <- sample(1:M, s, replace = T)
        states <- c(states, 1)
        prob <- runif(1)
        for (j in 2:M)
          if (prob > initial.dist.cumsum[j-1] &&
              prob <= initial.dist.cumsum[j])
            states[length(states)] <- j
      }
      else
      {
        state.conversion.mat <- GetStateConversionMatForR(M = M, s = s)
        M.extended <- ncol(state.conversion.mat)
        prob <- runif(1)
        states <- rev(state.conversion.mat[,1])
        for (j in 2:M.extended)
          if (prob > initial.dist.cumsum[j-1] &&
              prob <= initial.dist.cumsum[j])
          {
            states <- rev(state.conversion.mat[,j])
            break
          }
      }
      states <- c(states, rep(0, n))
    }
  }
  else
    states <- c(rep(-1, (length(initial.y.set) - 1)),
                initial.state,
                rep(0, n))
  initial.index <- length(initial.y.set) + 1
  last.index <- length(initial.y.set) + n
  if (is.MSM)
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
        z.dependent[k,] %*% as.matrix(gamma.dependent[,state]) +
        z.independent[k,] %*% as.matrix(gamma.independent) +
        rnorm(1,sd=sigma[state,1])
      for (lagged.index in 1:s)
      {
        # TODO: add lagged z variables in ar terms.
        lagged.state <- states[(k-lagged.index)]
        y[k] <- as.numeric(beta[lagged.index,lagged.state]) *
          (y[(k - lagged.index)] - mu[lagged.state]) +
          y[k]
      }
    }
  else
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
        t(rev(y[(k-s):(k-1)])) %*% as.numeric(beta[,state]) +
        z.dependent[k,] %*% as.matrix(gamma.dependent[,state]) +
        z.independent[k,] %*% as.matrix(gamma.independent) +
        rnorm(1,sd=sigma[state,1])
    }
  states <- states[initial.index:last.index]
  
  posterior.probs <- matrix(rep(0,M*n), ncol = M)
  for (i in initial.index:n)
    posterior.probs[i,states[i]] = 1
  msar.model <- list(theta = theta,
                     log.likelihood = Inf,
                     aic = Inf, bic = Inf,
                     posterior.probs.filtered = posterior.probs,
                     posterior.probs.smoothed = posterior.probs,
                     states = states,
                     call = match.call(),
                     M = M,
                     s = s,
                     is.beta.switching = is.beta.switching,
                     is.sigma.switching = is.sigma.switching,
                     is.MSM = is.MSM,
                     label = "msar.model")
  
  lagged.and.sample <- GetLaggedAndSample(y, s)
  y.sample <- lagged.and.sample$y.sample
  y.lagged <- lagged.and.sample$y.lagged
  
  
  if (is.MSM)
    return (list(y = y,
                 y.sample = y.sample,
                 y.lagged = y.lagged,
                 states = states,
                 msar.model = msar.model))
  else
    return (list(y = y,
                 y.sample = y.sample,
                 y.lagged = y.lagged,
                 states = states[initial.index:length(states)],
                 msar.model = msar.model))
}

#' Returns an (n + s) by replications matrix that represents samples of
#' length (n + s) observations, where s is a length of autoregressive terms.
#' @export
#' @title GenerateSamples
#' @name GenerateSamples
#' @param theta A list that represents the parameters of a model with items:
#' \item{transition.probs}{M by M matrix that contains transition probabilities}
#' \item{initial.dist}{M by 1 column that represents an initial distribution}
#' \item{beta}{s by 1 column for state-independent coefficients on AR(s)}
#' \item{mu}{M by 1 column that contains state-dependent mu}
#' \item{sigma}{M by 1 column that contains state-dependent sigma}
#' \item{gamma.dependent}{p_dep by M matrix that contains switching
#' coefficients for state-dependent exogenous variables}
#' \item{gamma.independent}{p_indep by 1 column that contains non-switching
#' coefficients for state-independent exogenous variables}
#' @param n The number of sample observations to be created.
#' @param replications The number of replications for samples.
#' @param is.MSM Determines whether the model follows MSM-AR. If it is set to be
#' TRUE, the model is assumed to be MSM-AR. 
#' @param burn.in counts for burn-in samples.
#' @return  A list with items:
#' \item{samples}{(n + length(initial.y.set)) by 1 column that represents
#' a sample appended with previous values used to estimate autoregressive terms}
#' \item{states}{n by 1 column that represents a sample of the model}
#' @examples
#' theta <- RandomTheta(M = 2, s = 3)
#' GenerateSamples(theta)
#' theta <- RandomTheta(M = 3, s = 2)
#' GenerateSamples(theta, n = 800)
GenerateSamples <- function(theta, n = 200, replications = 200,
                            initial.y.set = NULL,
                            is.MSM = FALSE, burn.in = 800)
{
  M <- ncol(theta$transition.probs)
  s <- nrow(as.matrix(theta$beta))
  probs <- runif(replications)
  states <- rep(1, replications)
  n.plus.burn.in <- n + burn.in
  
  
  # Format it as a switching model if not.
  if (ncol(as.matrix(theta$beta)) < M)
    theta$beta <- matrix(rep(theta$beta, M), ncol = M)
  if (length(theta$sigma) < M)
    theta$sigma <- rep(theta$sigma, M)
  theta$beta <- as.matrix(theta$beta)
  theta$mu <- as.matrix(theta$mu)
  theta$sigma <- as.matrix(theta$sigma)
  
  initial.dist.cumsum <- cumsum(theta$initial.dist)
  for (j in 2:(M^(s+1)))
    states[which(probs > initial.dist.cumsum[j-1] &&
                   probs <= initial.dist.cumsum[j])] <- j
  if (is.MSM)
  {
    state.conversion.mat <- GetStateConversionMatForR(M = M, s = s)
    return (sapply(states, function (state)
      tail(GenerateMSMSampleQuick(initial.states = rev(state.conversion.mat[,state]),
                                  theta = theta, n = n.plus.burn.in,
                                  initial.y.set = initial.y.set, 
                                  M = M, s = s), n)))
  }
  else
    return (sapply(states, function (state) 
      tail(GenerateMSISampleQuick(initial.state = state,
             theta = theta, n = n.plus.burn.in,
             initial.y.set = initial.y.set, 
             M = M, s = s), n)))
  
}

GenerateMSISampleQuick <- function(initial.state, theta, n,
                                   initial.y.set, M, s)
{
  if (is.null(initial.y.set))
    initial.y.set <- rnorm(s)
  y <- c(initial.y.set, rep(-Inf, n))
  states <- c(rep(-1, (s - 1)),
              initial.state,
              rep(0, n))
  
  initial.index <- s + 1
  last.index <- length(initial.y.set) + n
  
  for (k in initial.index:last.index)
  {
    previous.state <- states[(k-1)]
    trans.cumsum <- cumsum(theta$transition.probs[previous.state,])
    prob <- runif(1) # decision to switch
    state <- 1
    
    if (M > 1)
      for (j in 2:M)
        if (prob > trans.cumsum[j-1] && prob <= trans.cumsum[j]) {
          state <- j
          break
        }
    
    states[k] <- state
    y[k] <- theta$mu[state,1] +
      t(rev(y[(k-s):(k-1)])) %*% as.numeric(theta$beta[,state]) +
      rnorm(1,sd=theta$sigma[state,1])
  }
  
  return (y)
}

GenerateMSMSampleQuick <- function(initial.states, theta, n,
                                   initial.y.set, M, s)
{
  if (is.null(initial.y.set))
    initial.y.set <- rnorm(s)
  y <- c(initial.y.set, rep(-Inf, n))
  states <- c(initial.states,
              rep(0, n))
  
  initial.index <- s + 1
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
    y[k] <- theta$mu[state,1] +
      rnorm(1,sd=theta$sigma[state,1])
    for (lagged.index in 1:s)
    {
      lagged.state <- states[(k-lagged.index)]
      y[k] <- as.numeric(theta$beta[lagged.index,lagged.state]) *
        (y[(k - lagged.index)] - theta$mu[lagged.state]) +
        y[k]
    }
  }
  
  
  return (y)
}

