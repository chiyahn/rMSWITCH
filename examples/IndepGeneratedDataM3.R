library(rmrs)
library(normalregMix)

# Generates a sample for M = 3
GenerateSampleM3 <- function(n, beta, mu1, mu2, mu3, sigma1, sigma2, sigma3, p12, p13, p21, p23, p31, p32)
{
  y <- as.list(rnorm(s))
  s <- length(beta)
  states <- as.list(rep(1,s))
  for (i in (s+1):n)
  {
    prob <- runif(1,0,1) # decision to switch
    if (states[[i-1]] == 1)
      states[[i]] = 1 + (prob < p12) + 2 * (p12 < prob && prob < (p12 + p23))
    else if (states[[i-1]] == 2)
      states[[i]] = 2 + (prob < p23) - (p23 < prob && prob < (p23 + p21))
    else
      states[[i]] = 3 - (prob < p32) - 2 * (p32 < prob && prob < (p32 + p31))
    
    if (states[[i]] == 1)
      y[[i]] <- mu1 + t(as.numeric(y[(i-s):(i-1)])) %*% as.numeric(beta) + rnorm(1,sd=sigma1)
    else if (states[[i]] == 2)
      y[[i]] <- mu2 + t(as.numeric(y[(i-s):(i-1)])) %*% as.numeric(beta) + rnorm(1,sd=sigma2)
    else
      y[[i]] <- mu3 + t(as.numeric(y[(i-s):(i-1)])) %*% as.numeric(beta) + rnorm(1,sd=sigma3)
    
  }
  return (as.numeric(y))
}

# Get elementwise-average of list of matrices
GetMatrixAvg <- function (list.of.matrices)
{
  Reduce("+", list.of.matrices) / length(list.of.matrices)
}

set.seed(123456)
sample.count <- 1
sample.set <- list()

# model specification
n <- 200
p12 <- 0.4
p13 <- 0.3
p21 <- 0.3
p23 <- 0.2
p31 <- 0.3
p32 <- 0.4
mu1 <- -1
mu2 <- 1
mu3 <- 0.5
sigma1 <- 1.3
sigma2 <- 1.1
sigma3 <- 0.5
beta <- c(0.5, 0.3)
s <- length(beta)

# Creates samples
for (sample.index in 1:sample.count)
  sample.set[[sample.index]] <- GenerateSampleM3(n, beta, mu1, mu2, mu3, 
                                                 sigma1, sigma2, sigma3, 
                                                 p12, p13, p21, p23, p31, p32)

# Performs estimation
results <- lapply(sample.set, MRSMLEIndep, M = 3, s = s, maxit = 30)
thetas <- lapply(results, "[[", "theta")
betas <- lapply(thetas, "[[", "beta")
mus <- lapply(thetas, "[[", "mu")
sigmas <- lapply(thetas, "[[", "sigma")
tranprobs <- lapply(thetas, "[[", "transition.probs")
likelihoods <- lapply(results, "[[", "likelihood")

betas.avg <- GetMatrixAvg(lapply(thetas, "[[", "beta"))
mus.avg <- GetMatrixAvg(lapply(thetas, "[[", "mu"))
sigmas.avg <- GetMatrixAvg(lapply(thetas, "[[", "sigma"))
tranprobs.avg <- GetMatrixAvg(lapply(thetas, "[[", "transition.probs"))

betas.avg
mus.avg
sigmas.avg
tranprobs.avg
likelihoods
