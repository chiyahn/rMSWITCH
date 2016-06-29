library(rmrs)
library(normalregMix)

# Generates a sample for M = 2
GenerateSampleM2 <- function(n, beta, mu1, mu2, sigma1, sigma2, p12, p21)
{
  y <- as.list(rnorm(s))
  s <- length(beta)
  states <- as.list(rep(1,s))
  for (i in (s+1):n)
  {
    prob <- runif(1,0,1) # decision to switch
    if (states[[i-1]] == 1)
      states[[i]] = (prob < p12) + 1
    else
      states[[i]] = 2 - (prob < p21)
    
    if (states[[i]] == 1)
      y[[i]] <- mu1 + t(as.numeric(y[(i-s):(i-1)])) %*% as.numeric(beta) + rnorm(1,sd=sigma1)
    else
      y[[i]] <- mu2 + t(as.numeric(y[(i-s):(i-1)])) %*% as.numeric(beta) + rnorm(1,sd=sigma2)
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
p12 <- 0.6
p21 <- 0.7
mu1 <- -1
mu2 <- 1
sigma1 <- 1
sigma2 <- 2
beta <- 0.7
s <- length(beta)

# Creates samples
for (sample.index in 1:sample.count)
  sample.set[[sample.index]] <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21)

# Performs estimation
results <- lapply(sample.set, MRSMLEIndep, M = 2, s = s)
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
