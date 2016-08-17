GenerateSamplesStationary <- function (theta, n = 200, replications = 200,
                                      is.MSM = FALSE, burn.in = 400)
{
  if (is.MSM)
    stop("MSM models are currently not supported.")
  
  M <- ncol(theta$transition.probs)
  if (sum(abs(theta$beta)) >= 1)
    stop("The corresponding model is not stationary as the autoregressive
         coefficients are not located inside a unit circle.")
  s <- nrow(as.matrix(theta$beta))
  probs <- runif(replications)
  states <- rep(1, replications)
  initial.y.values <- matrix(0, ncol = replications, nrow = s)
  
  # Format it as a switching model if not.
  if (ncol(as.matrix(theta$beta)) < M)
    theta$beta <- matrix(rep(theta$beta, M), ncol = M)
  if (length(theta$sigma) < M)
    theta$sigma <- rep(theta$sigma, M)
  theta$beta <- as.matrix(theta$beta)
  theta$mu <- as.matrix(theta$mu)
  theta$sigma <- as.matrix(theta$sigma)
  
  if (M == 1)
  {
    initial.y.values <- matrix(rnorm(replications*s, 
                                mean = 0, 
                                sd = as.numeric(theta$sigma) / 
                                      ((1 - sum(theta$beta * 
                                                theta$beta)))),
                               ncol = replications)
    states <- rep(1,replications)
    samples <- matrix(0, ncol = replications, nrow = (n + s))
    for (i in 1:replications)
      samples[,i] <- GenerateSampleQuick(initial.state = states[i], 
                                         initial.y.set = initial.y.values[,i],
                                         theta = theta,
                                         n = n, M = M, s = s,
                                         is.MSM = is.MSM)
    
    return (samples)
  }
  
  # case M > 1
  initial.y.values <- matrix(0,
                             ncol = replications,
                             nrow = s)
  stationary.dist <- ComputeStationaryDist(theta$transition.probs)
  stationary.dist.cumsum <- cumsum(theta$stationary.dist)
  for (j in 2:M)
    states[which(probs > stationary.dist.cumsum[j-1] &&
                   probs <= stationary.dist.cumsum[j])] <- j
  
  samples <- matrix(0, ncol = replications, nrow = (n + burn.in + s))
  for (i in 1:replications)
    samples[,i] <- GenerateSampleQuick(initial.state = states[i], 
                                       initial.y.set = initial.y.values[,i],
                                       theta = theta,
                                       n = (n + burn.in), M = M, s = s,
                                       is.MSM = is.MSM)
  
  return (samples[(burn.in + 1):(nrow(samples)),])
}