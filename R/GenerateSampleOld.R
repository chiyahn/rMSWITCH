# Generates a sample for M = 2
GenerateSampleM2 <- function(n, beta, mu1, mu2, sigma1, sigma2, p12, p21,
                             z.dependent = NULL,
                             z.independent = NULL,
                             gamma.dependent = matrix(rep(0,2), ncol = 2),
                             gamma.independent = as.matrix(0))
{
  y <- as.list(rnorm(s))
  s <- length(beta)
  states <- as.list(rep(1,s))
  M <- 2
  
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,(n+s)),ncol=1)
  else 
    z.dependent <- rbind(matrix(rep(Inf, s*ncol(z.dependent)), 
                                ncol = ncol(z.dependent)),
                         as.matrix(z.dependent))
  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,(n+s)),ncol=1)
  else
    z.independent <- rbind(matrix(rep(Inf, s*ncol(z.independent)), 
                                  ncol = ncol(z.independent)),
                           as.matrix(z.independent))
  
  initial.index = (s+1)
  last.index = (s+n)
  for (i in initial.index:last.index)
  {
    prob <- runif(1) # decision to switch
    if (states[[i-1]] == 1)
      states[[i]] = (prob > (1-p12)) + 1
    else
      states[[i]] = 2 - (prob < p21)
    
    if (states[[i]] == 1)
      y[[i]] <- mu1 + t(unlist(y[(i-s):(i-1)])) %*% as.numeric(beta) +
        z.dependent[i,] %*% as.matrix(gamma.dependent[,1]) +
        z.independent[i,] %*% as.matrix(gamma.independent) +
        rnorm(1,sd=sigma1)
    else
      y[[i]] <- mu2 + t(unlist(y[(i-s):(i-1)])) %*% as.numeric(beta) +
        z.dependent[i,] %*% as.matrix(gamma.dependent[,2]) +
        z.independent[i,] %*% as.matrix(gamma.independent) +
        rnorm(1,sd=sigma2)
  }
  return (as.numeric(y[initial.index:last.index]))
}

# Generates a sample for M = 3
GenerateSampleM3 <- function(n, beta, mu1, mu2, mu3, sigma1, sigma2, sigma3,
                             p12, p13, p21, p23, p31, p32,
                             z.dependent = NULL,
                             z.independent = NULL,
                             gamma.dependent = matrix(rep(0,3), ncol = 3),
                             gamma.independent = as.matrix(0))
{
  y <- as.list(rnorm(s))
  s <- length(beta)
  states <- as.list(rep(1,s))
  transition.probs <- matrix(c((1-p12-p13),p12,p13,
                               p21,(1-p21-p23),p23,
                               p31,p32,(1-p31-p32)), 
                             ncol = 3, byrow =T)
  M <- ncol(transition.probs)
  
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,(n+s)),ncol=1)
  else 
    z.dependent <- rbind(matrix(rep(Inf, s*ncol(z.dependent)), 
                                ncol = ncol(z.dependent)),
                         as.matrix(z.dependent))
  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,(n+s)),ncol=1)
  else
    z.independent <- rbind(matrix(rep(Inf, s*ncol(z.independent)), 
                                  ncol = ncol(z.independent)),
                           as.matrix(z.independent))
  initial.index = (s+1)
  last.index = (s+n)
  for (i in initial.index:last.index)
  {
    prob <- runif(1) # decision to switch
    trans.cumsum <- cumsum(theta$transition.probs[states[[i-1]],])
    states[[i]] = 1
    for (j in 2:M)
      if (prob > trans.cumsum[j-1] && prob <= trans.cumsum[j]) {
        states[[i]] <- j
        break
      }
    
    if (states[[i]] == 1)
      y[[i]] <- mu1 + t(unlist(y[(i-s):(i-1)])) %*% as.numeric(beta) +
        z.dependent[i,] %*% as.matrix(gamma.dependent[,1]) +
        z.independent[i,] %*% as.matrix(gamma.independent) +
        rnorm(1,sd=sigma1)
    else if (states[[i]] == 2)
      y[[i]] <- mu2 + t(unlist(y[(i-s):(i-1)])) %*% as.numeric(beta) +
        z.dependent[i,] %*% as.matrix(gamma.dependent[,2]) +
        z.independent[i,] %*% as.matrix(gamma.independent) +
        rnorm(1,sd=sigma2)
    else
      y[[i]] <- mu3 + t(unlist(y[(i-s):(i-1)])) %*% as.numeric(beta) +
        z.dependent[i,] %*% as.matrix(gamma.dependent[,3]) +
        z.independent[i,] %*% as.matrix(gamma.independent) +
        rnorm(1,sd=sigma3)
    
  }
  return (as.numeric(y[initial.index:last.index]))
}