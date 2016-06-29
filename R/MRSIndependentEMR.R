ExpectationMaximizationIndepR <- function(y, y.lagged, z_dependent, z_independent, theta0, maxit, epsilon)
{
  # load c++ code for EtaIndep
  # setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  # sourceCpp("cppEtaIndep.cpp")
  # sourceCpp("cppFilterIndep.cpp")
  # sourceCpp("cppSmooth.cpp")
  # sourceCpp("cppEstimateStates.cpp")
  # sourceCpp("cppMaximizationStepIndep.cpp")

  thetas <- list()
  thetas[[1]] <- theta0
  theta <- theta0
  likelihoods <- list(-99999)
  likelihood <- likelihoods[[1]]
  for (i in 2:maxit)
  {
    e.step <- ExpectationStepR(y, y.lagged, z_dependent, z_independent, thetas[[i-1]])
    thetas[[i]] <- MaximizationStepIndepR(y, y.lagged, z_dependent, z_independent,
                                         theta$beta, theta$mu, theta$sigma,
                                         theta$gamma.dependent, theta$gamma.independent,
                                         theta$transition.probs, theta$initial.dist,
                                         e.step$xi.k, e.step$xi.n, theta0$sigma)
    theta <- thetas[[i]]

    likelihoods[[i]] <- e.step$likelihood
    likelihood <- likelihoods[[i]]
    if (abs(likelihoods[[i]] - likelihoods[[i-1]]) < epsilon)
      break
  }

  # !! CHECK: if EstimateStates is called from C++, add one.
  states <- EstimateStatesR(y, y.lagged, z_dependent, z_independent,
                          theta$beta, theta$mu, theta$sigma,
                          theta$gamma.dependent, theta$gamma.independent)
  #states <- states + 1
  # !! CHECK: if EstimateStates is called from C++, add one.
  return (list(theta = theta, likelihood = likelihood,
    likelihoods = unlist(likelihoods), states = states))
}

ExpectationStepR <- function(y, y.lagged, z_dependent, z_independent, theta) {
  filter <- FilterIndepR(y, y.lagged, z_dependent, z_independent,
                        theta$beta, theta$mu, theta$sigma,
                        theta$gamma.dependent, theta$gamma.independent,
                        theta$transition.probs, theta$initial.dist)
  xi.n <- SmoothR(filter$xi.k, theta$transition.probs)
  return (list(xi.k = filter$xi.k, xi.n = xi.n, likelihood = filter$likelihood))
}


MaximizationStepIndepR <- function(y, y.lagged, z_dependent, z_independent,
                                  beta, mu, sigma,
                                  gamma.dependent, gamma.independent,
                                  transition.probs, initial.dist,
                                  xi.k, xi.n, sigma0.origin) {
  beta0 <- as.matrix(beta)
  mu0 <- mu
  sigma0 <- sigma
  gamma_dependent0 <- gamma.dependent
  gamma_independent0 <- gamma.independent

  transition.probs0 = transition.probs
  initial.dist0 = initial.dist

  M <- nrow(transition.probs0)
  s <- nrow(as.matrix(beta))
  n <- length(y)

  transition.probs <- matrix(rep(0, M*M), ncol = M)

  # 1. Estimates for transition probs
  for (i in 1:M)
  {
    total <- 0
    for (k in 1:n)
      total <- xi.n[k,i] + total
    for (j in 1:M)
    {
      prob.ij <- 0
      for (k in 2:n)
      {
        prob.ij.k <- xi.n[k,j] * transition.probs0[i,j] * xi.k[(k-1),][i] /
                    (transition.probs0 %*% xi.k[(k-1),])[j]
        prob.ij <- prob.ij + prob.ij.k
      }
      transition.probs[i,j] <- prob.ij / total
      transition.probs[i,j] <- max(transition.probs[i,j], 0.02) # hard constraint
      transition.probs[i,j] <- min(transition.probs[i,j], 0.98) # hard constraint
    }
    transition.probs[i,] <- transition.probs[i,] /
                            sum(transition.probs[i,]) # normalize again
  }

  # 2. Estimates for beta, mu, sigma
  beta <- beta0 # s-length vec
  mu <- mu0  # M-length vec
  sigma <- sigma0 # M-length vec
  beta <- beta0 # s-length vec
  gamma_dependent <- gamma_dependent0 # p_dep by M mat
  gamma_independent <- gamma_independent0 # p_indep-length vec

  # mu
  for (j in 1:M)
    mu[j] <- sum(xi.n[,j] * (y - y.lagged %*% beta0)) / sum(xi.n[,j])

  # gamma_dependent
  if (!sum(gamma_dependent == 0) == length(gamma_dependent)) # validity check
    for (j in 1:M)
    {
      gamma_dependent_part_one <- 0
      gamma_dependent_part_two <- 0
      for (k in 1:n)
      {
        gamma_dependent_part_one <- xi.n[k,j] *
          (z_dependent[k,] %*% t(z_dependent[k,])) + gamma_dependent_part_one
        gamma_dependent_part_two <- xi.n[k,j] * z_dependent[k,] *
          (y[k] -
            y.lagged[k,] %*% beta0 -
            z_independent[k,] %*% gamma_independent0 -
            mu[j]) + gamma_dependent_part_two
      }
      gamma_dependent[,j] <- solve(gamma_dependent_part_one) %*%
                            gamma_dependent_part_two
    }

  beta.part.one <- 0
  beta.part.two <- 0
  for (k in 1:n)
  {
    prop.sum <- 0
    for (j in 1:M)
    {
      prop <- xi.n[k,j] / sigma0[j]^2
      prop.sum <- prop + prop.sum
      beta.part.two <- prop * y.lagged[k,] *
        (y[k] - z_independent[k,] %*% as.matrix(gamma_independent) -
        z_dependent[k,] %*% as.matrix(gamma_dependent[,j]) - mu[j]) +
        beta.part.two
    }
    beta.part.one <- prop.sum * (y.lagged[k,] %*% t(y.lagged[k,])) + beta.part.one
  }
  beta <- solve(beta.part.one) %*% beta.part.two

  # gamma_independent
  if (!sum(gamma_independent == 0) == length(gamma_independent)) # validity check
  {
    gamma_independent_part_one <- 0
    gamma_independent_part_two <- 0
    for (k in 1:n)
    {
      prop.sum <- 0
      for (j in 1:M)
      {
        prop <- xi.n[k,j] / sigma0[j]^2
        prop.sum <- prop + prop.sum
        gamma_independent_part_two <- prop * z_independent[k,] *
          (y[k] - y.lagged[k,] %*% as.matrix(beta) -
          z_dependent[k,] %*% as.matrix(gamma_dependent)[,j] - mu[j]) +
          gamma_independent_part_two
      }
      gamma_independent_part_one <- prop.sum *
        (z_independent[k,] %*% t(z_independent[k,])) + gamma_independent_part_one
    }
    gamma_independent <- solve(gamma_independent_part_one) %*%
                          gamma_independent_part_two
  }

  # sigma (dependent)
  for (j in 1:M)
  {
    sigma[j] <- 0
    for (k in 1:n)
      sigma[j] <- (xi.n[k,j] / sum(xi.n[,j])) *
        (y[k] - y.lagged[k,] %*% beta -
          z_independent[k,] %*% as.matrix(gamma_independent) -
          z_dependent[k,] %*% as.matrix(gamma_dependent)[,j] - mu[j])^2 + sigma[j]
    sigma[j] <- max(sqrt(sigma[j]), 0.01 * sigma0.origin) # hard constraint; sigma >= 0.01 * sigma.hat
  }

  # !! CHECK: if EstimateStates is called from C++, add one.
  states <- EstimateStatesR(y, y.lagged, z_dependent, z_independent,
                           beta, mu, sigma,
                           gamma_dependent, gamma_independent)
  # states <- states + 1
  # !! CHECK: if EstimateStates is called from C++, add one.

  initial.dist <- sapply(seq(1,M), function(j, states, n) sum(states==j)/n, states, n)

  return (list(beta = beta, mu = mu, sigma = sigma,
               gamma.dependent = gamma_dependent,
               gamma.independent = gamma_independent,
               transition.probs = transition.probs, initial.dist = initial.dist))
}


# Returns n by M matrix that whose element eta_kj determines
# the probability of kth observation belonging to jth regime
EtaIndepR <- function(y, y.lagged, z_dependent, z_independent,
                      beta, mu, sigma,
                      gamma_dependent, gamma_independent) {
  M <- length(mu)
  n <- length(y)
  eta <- matrix(rep(0,n*M), ncol = M)


  for (j in 1:M)
  {
    eta[,j] <- -(y - (y.lagged %*% as.matrix(beta)) -
                   (z_dependent %*% as.matrix(gamma_dependent)[,j]) -
                   (z_independent %*% as.matrix(gamma_independent)) - mu[j])^2
    eta[,j] <- exp(eta[,j] / (2 * (sigma[j] ^ 2))) / (sqrt(2*pi) * sigma[j])
  }

  return (eta)
}

# Estimate xi_{k|k}; xi.k is an n by M matrix.
FilterIndepR <- function(y, y.lagged, z_dependent, z_independent,
                        beta, mu, sigma,
                        gamma_dependent, gamma_independent,
                        transition.probs, initial.dist)
{
  n <- length(y)
  M <- ncol(transition.probs)
  xi.k <- matrix(rep(0,n*M), ncol = M)

  likelihood <- 0

  eta <- EtaIndepR(y, y.lagged, z_dependent, z_independent,
                  beta, mu, sigma,
                  gamma_dependent, gamma_independent) # n by M matrix

  # k=1 uses initial distribution
  xi.k[1,] <- eta[1,] * initial.dist
  total <- sum(xi.k[1,])
  xi.k[1,] <- xi.k[1,] / total
  likelihood <- log(total) + likelihood
  for (k in 2:n)
  {
    xi.k[k,] <- eta[k,] * (transition.probs %*% xi.k[(k-1),])
    total <- sum(xi.k[k,])
    xi.k[k,] <- xi.k[k,] / total
    likelihood <- log(total) + likelihood
  }

  return(list(xi.k = xi.k, likelihood = likelihood))
}

# Returns a smooth inference by backwards recursion; xi.n is an n by M matrix.
SmoothR <- function(xi.k, transition.probs)
{
  n <- nrow(xi.k)
  M <- ncol(transition.probs)
  xi.n <- matrix(rep(0,n*M), ncol = M)

  xi.n[n,] <- xi.k[n,]

  for (k in (n-1):1)
    xi.n[k,] <- xi.k[k,] *
    (transition.probs %*% (xi.n[(k+1),] / (transition.probs %*% xi.k[k,])))

  return (xi.n)
}

# Determines which state each observation belongs to based on theta
EstimateStatesR <- function(y, y.lagged, z_dependent, z_independent,
                            beta, mu, sigma,
                            gamma.dependent, gamma.independent) {
  eta <- EtaIndepR(y, y.lagged, z_dependent, z_independent,
                  beta, mu, sigma,
                  gamma.dependent, gamma.independent)
  return (apply(eta, 1, function(i) (which(i==max(i)))[1]))
}
