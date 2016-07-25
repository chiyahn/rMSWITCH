ExpectationMaximizationIndepR <- function(y, y.lagged, z.dependent, z.independent, theta0, maxit, epsilon)
{
  # load c++ code for EtaIndep
  # setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  # sourceCpp("cppEtaIndep.cpp")
  # sourceCpp("cppFilterIndep.cpp")
  # sourceCpp("cppSmooth.cpp")
  # sourceCpp("cppEstimateStates.cpp")
  # sourceCpp("cppMaximizationStepIndep.cpp")

  # if z.dependent/z.independent = NULL, set it as a zero vector.
  n <- length(y)
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,n))

  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,n))


  thetas <- list()
  thetas[[1]] <- theta0
  theta <- theta0
  likelihoods <- list(-Inf)
  likelihood <- likelihoods[[1]]
  for (i in 2:(maxit+1))
  {
    e.step <- ExpectationStepR(y, y.lagged, z.dependent, z.independent, thetas[[i-1]])

    likelihoods[[i]] <- e.step$likelihood


    # stop if 1. the difference in likelihoods is small enough
    # or 2. likelihood decreases (a decrease is due to locating local max.
    # out of hard constraints in maximization step, which suggests that this
    # is not a good candidate anyway.)
    if ((likelihoods[[i]] - likelihoods[[i-1]]) < epsilon)
      break

    likelihood <- likelihoods[[i]]
    theta <- thetas[[(i-1)]]
    thetas[[i]] <- MaximizationStepIndepR(y, y.lagged, z.dependent, z.independent,
                                         theta$beta, theta$mu, theta$sigma,
                                         theta$gamma.dependent, theta$gamma.independent,
                                         theta$transition.probs, theta$initial.dist,
                                         e.step$xi.k, e.step$xi.past,
                                         e.step$xi.n, theta0$sigma)



  }

  # !! CHECK: if EstimateStates is called from C++, add one.
  states <- EstimateStatesR(y, y.lagged, z.dependent, z.independent,
                          theta$beta, theta$mu, theta$sigma,
                          theta$gamma.dependent, theta$gamma.independent)
  #states <- states + 1
  # !! CHECK: if EstimateStates is called from C++, add one.
  return (list(theta = theta, likelihood = likelihood,
    likelihoods = unlist(likelihoods), states = states, thetas = thetas))
}

ExpectationStepR <- function(y, y.lagged, z.dependent, z.independent, theta) {
  filter <- FilterIndepR(y, y.lagged, z.dependent, z.independent,
                        theta$beta, theta$mu, theta$sigma,
                        theta$gamma.dependent, theta$gamma.independent,
                        theta$transition.probs, theta$initial.dist)
  xi.n <- SmoothR(filter$xi.k, theta$transition.probs)
  return (list(xi.k = filter$xi.k, xi.past = filter$xi.past,
               xi.n = xi.n, likelihood = filter$likelihood))
}


MaximizationStepIndepR <- function(y, y.lagged, z.dependent, z.independent,
                                  beta, mu, sigma,
                                  gamma.dependent, gamma.independent,
                                  transition.probs, initial.dist,
                                  xi.k, xi.past, xi.n, sigma0.origin) {
  transition.probs0 <- transition.probs
  initial.dist0 <- initial.dist
  beta0 <- as.matrix(beta)
  mu0 <- mu
  sigma0 <- sigma
  gamma_dependent0 <- gamma.dependent
  gamma_independent0 <- gamma.independent

  transition.probs <- matrix(rep(0, M*M), ncol = M)
  beta <- beta0 # s-length vec
  mu <- mu0  # M-length vec
  sigma <- sigma0 # M-length vec
  beta <- beta0 # s-length vec
  gamma_dependent <- gamma_dependent0 # p_dep by M mat
  gamma_independent <- gamma_independent0 # p_indep-length vec

  M <- nrow(transition.probs0)
  s <- nrow(as.matrix(beta))
  n <- length(y)
  n.minus.one <- n - 1

  # 1. Estimates for transition probs
  for (i in 1:M)
  {
    total <- 0
    for (k in 1:n.minus.one)
      total <- xi.n[k,i] + total
    for (j in 1:M)
    {
      prob.ij <- 0
      for (k in 2:n)
      {
        prob.ij.k <- xi.n[k,j] * transition.probs0[i,j] * xi.k[(k-1),i] /
                    (xi.past[k,j])
        prob.ij <- prob.ij + prob.ij.k
      }
      transition.probs[i,j] <- prob.ij / total
      # Enforce ub/lb.
      transition.probs[i,j] <- max(transition.probs[i,j], 0.05) # hard constraint
      transition.probs[i,j] <- min(transition.probs[i,j], 0.95) # hard constraint
    }
    transition.probs[i,] <- transition.probs[i,] /
                            sum(transition.probs[i,]) # normalize again
  }

  # 2. Estimates for beta, mu, sigma

  # mu
  for (j in 1:M)
    mu[j] <- sum(xi.n[,j] * (y - y.lagged %*% beta0 -
      z.dependent %*% gamma_dependent[,j] -
      z.independent %*% gamma.independent)) /
      sum(xi.n[,j])

  # gamma_dependent
  if (!sum(abs(gamma_dependent) < 0.00001) == length(gamma_dependent)) # validity check
    for (j in 1:M)
    {
      gamma_dependent_part_one <- 0
      gamma_dependent_part_two <- 0
      for (k in 1:n)
      {
        gamma_dependent_part_one <- xi.n[k,j] *
          (z.dependent[k,] %*% t(z.dependent[k,])) + gamma_dependent_part_one
        gamma_dependent_part_two <- xi.n[k,j] * z.dependent[k,] *
          (y[k] -
            y.lagged[k,] %*% beta0 -
            z.independent[k,] %*% gamma_independent0 -
            mu[j]) + gamma_dependent_part_two
      }
      gamma_dependent[,j] <- solve(gamma_dependent_part_one) %*%
                            gamma_dependent_part_two
    }

  # beta
  beta.part.one <- 0
  beta.part.two <- 0
  for (k in 1:n)
  {
    prop.sum <- 0
    for (j in 1:M)
    {
      prop <- xi.n[k,j] / (sigma0[j] * sigma0[j])
      prop.sum <- prop + prop.sum
      beta.part.two <- prop * y.lagged[k,] *
        (y[k] - z.independent[k,] %*% as.matrix(gamma_independent) -
        z.dependent[k,] %*% as.matrix(gamma_dependent[,j]) - mu[j]) +
        beta.part.two
    }
    beta.part.one <- prop.sum * (y.lagged[k,] %*% t(y.lagged[k,])) + beta.part.one
  }
  beta <- solve(beta.part.one) %*% beta.part.two

  # gamma_independent
  if (!sum(abs(gamma_independent) < 0.00001) == length(gamma_independent)) # validity check
  {
    gamma_independent_part_one <- 0
    gamma_independent_part_two <- 0
    for (k in 1:n)
    {
      prop.sum <- 0
      for (j in 1:M)
      {
        prop <- xi.n[k,j] / (sigma0[j] * sigma0[j])
        prop.sum <- prop + prop.sum
        gamma_independent_part_two <- prop * z.independent[k,] *
          (y[k] - y.lagged[k,] %*% as.matrix(beta) -
          z.dependent[k,] %*% as.matrix(gamma_dependent)[,j] - mu[j]) +
          gamma_independent_part_two
      }
      gamma_independent_part_one <- prop.sum *
        (z.independent[k,] %*% t(z.independent[k,])) + gamma_independent_part_one
    }
    gamma_independent <- solve(gamma_independent_part_one) %*%
                          gamma_independent_part_two
  }

  # sigma (dependent)
  for (j in 1:M)
  {
    sigma[j] <- 0
    for (k in 1:n)
    {
      res <- (y[k] - y.lagged[k,] %*% beta -
        z.independent[k,] %*% as.matrix(gamma_independent) -
        z.dependent[k,] %*% as.matrix(gamma_dependent)[,j] - mu[j])
      sigma[j] <- (xi.n[k,j] / sum(xi.n[,j])) * res * res + sigma[j]
    }
    sigma[j] <- max(sqrt(sigma[j]), 0.5 * sigma0.origin) # hard constraint; sigma >= 0.2 * sigma.hat
  }

  # initial.dist
  initial.dist <- xi.n[1,]

  return (list(beta = beta, mu = mu, sigma = sigma,
               gamma.dependent = gamma_dependent,
               gamma.independent = gamma_independent,
               transition.probs = transition.probs, initial.dist = initial.dist))
}


# Returns n by M matrix that whose element eta_kj determines
# the probability of kth observation belonging to jth regime
EtaIndepR <- function(y, y.lagged, z.dependent, z.independent,
                      beta, mu, sigma,
                      gamma_dependent, gamma_independent) {
  M <- length(mu)
  n <- length(y)
  eta <- matrix(rep(0,n*M), ncol = M)


  for (j in 1:M)
  {
    eta[,j] <- -(y - (y.lagged %*% as.matrix(beta)) -
                   (z.dependent %*% as.matrix(gamma_dependent)[,j]) -
                   (z.independent %*% as.matrix(gamma_independent)) - mu[j])^2
    eta[,j] <- exp(eta[,j] / (2 * (sigma[j] ^ 2))) / (sqrt(2*pi) * sigma[j])
  }

  return (eta)
}

# Estimate xi_{k|k}; xi.k is an n by M matrix.
FilterIndepR <- function(y, y.lagged, z.dependent, z.independent,
                        beta, mu, sigma,
                        gamma_dependent, gamma_independent,
                        transition.probs, initial.dist)
{
  n <- length(y)
  M <- ncol(transition.probs)
  xi.k <- matrix(rep(0,n*M), ncol = M)
  xi.past <- matrix(rep(0,n*M), ncol = M)
  likelihood <- 0

  for (k in 1:n)
  {
    min.index <- -1
    min.value <- Inf
    ratios <- rep(-Inf, M)
    row.sum <- 0

    if (k > 1)
      xi.past[k,] <- transition.probs %*% xi.k[(k-1),]
    else
      xi.past[k,] <- initial.dist

    for (j in 1:M)
    {
      xi.k[k,j] <- (y[k] - (y.lagged[k,] %*% as.matrix(beta)) -
                (z.dependent[k,] %*% as.matrix(gamma_dependent)[,j]) -
                (z.independent[k,] %*% as.matrix(gamma_independent)) - mu[j])^2
      xi.k[k,j] <- xi.k[k,j] / (2 * (sigma[j] * sigma[j]))
      if (min.value > xi.k[k,j])
      {
        min.value <- xi.k[k,j]
        min.index <- j
      }
      ratios[j] <- xi.past[k,j] / sigma[j]
    }

    for (j in 1:M)
    {
      if (j == min.index)
        xi.k[k,j] <- 1.0
      else
        xi.k[k,j] <- (ratios[j] / ratios[min.index]) * exp(min.value - xi.k[k,j])
      row.sum <- xi.k[k,j] + row.sum
    }

    xi.k[k,] = xi.k[k,] / sum(xi.k[k,])

    likelihood = likelihood + log(row.sum) - min.value + log(ratios[min.index])
  }
  likelihood = likelihood - n * log(2*pi) / 2

  return(list(xi.k = xi.k, xi.past = xi.past, likelihood = likelihood))
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
    (t(transition.probs) %*% (xi.n[(k+1),] / (transition.probs %*% xi.k[k,])))

  return (xi.n)
}

# Determines which state each observation belongs to based on theta
EstimateStatesR <- function(y, y.lagged, z.dependent, z.independent,
                            beta, mu, sigma,
                            gamma.dependent, gamma.independent) {
  eta <- EtaIndepR(y, y.lagged, z.dependent, z.independent,
                  beta, mu, sigma,
                  gamma.dependent, gamma.independent)
  return (apply(eta, 1, function(i) (which(i==max(i)))[1]))
}

ComputeStationaryDist <- function(transition.probs, M)
{
  eigen.list <- eigen(t(transition.probs))
  eigen.values <- eigen.list$values
  eigen.vectors <- eigen.list$vectors
  stationary.dist.index <- -1
  for (i in 1:M)
    if (eigen.values[i] == 1)
      stationary.dist.index = i;
  stationary.dist <- eigen.vectors[,stationary.dist.index]
  stationary.dist <- abs(stationary.dist) / sum(abs(stationary.dist));
  return (stationary.dist)
}
