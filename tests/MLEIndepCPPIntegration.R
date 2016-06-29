
#install.packages("MswM")
library(MSwM)
library(normalregMix)
library(rMSWITCH)
library(Rcpp)
library(RcppArmadillo)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Generates a sample for M = 2
GenerateSampleM2 <- function(n, beta, mu1, mu2, sigma1, sigma2, p12, p21, 
                             z.dependent = NULL, 
                             z.independent = NULL,
                             gamma.dependent = as.matrix(rep(0,2), ncol = 2), 
                             gamma.independent = as.matrix(0))
{
  y <- as.list(rnorm(s))
  s <- length(beta)
  states <- as.list(rep(1,s))
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,n),ncol=1)
  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,n),ncol=1) 
  
  for (i in (s+1):n)
  {
    prob <- runif(1,0,1) # decision to switch
    if (states[[i-1]] == 1)
      states[[i]] = (prob < p12) + 1
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
  return (as.numeric(y))
}

# Generates a sample for M = 3
GenerateSampleM3 <- function(n, beta, mu1, mu2, mu3, sigma1, sigma2, sigma3, 
                             p12, p13, p21, p23, p31, p32,
                             z.dependent = NULL, 
                             z.independent = NULL,
                             gamma.dependent = as.matrix(rep(0,3), ncol = 3), 
                             gamma.independent = as.matrix(0))
{
  y <- as.list(rnorm(s))
  s <- length(beta)
  states <- as.list(rep(1,s))
  
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,n),ncol=1)
  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,n),ncol=1) 
  
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
  return (as.numeric(y))
}
# Generates exogeneous variable sample
GenerateExo <- function(n, p)
{
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # changes wd; only works running in RStudio
  return (matrix(as.numeric(read.csv("interestRates.csv")[1:(n*p),1]), ncol = p))
}

# Get elementwise-average of list of matrices
GetMatrixAvg <- function (list.of.matrices)
{
  Reduce("+", list.of.matrices) / length(list.of.matrices)
}


MLEIndepCPPComparison <- function(y = y, z = NULL, iterations = 2, z.is.switching = FALSE, M = 3, s = 2, theta.initial = NULL, 
                        epsilon = 1e-06, maxit = 60, short.n = 200, short.iterations = 10) {
  
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
  theta <- theta.initial
  short.result <- ExpectationMaximizationIndepR(y = y.sample, y.lagged = y.lagged, 
                                                z_dependent = z.dependent, z_independent = z.independent,
                                                theta = theta, 
                                                maxit = iterations, epsilon = epsilon) 
  short.result2 <- ExpectationMaximizationIndep (y.sample, y.lagged, z.dependent, z.independent,
                                                 theta, iterations, epsilon)
  short.result2 <- ExpectationMaximizationIndep (y.sample, y.lagged, z.dependent, z.independent,
                                                 short.result2$theta, iterations, epsilon)
  print(short.result$likelihood)
  print(short.result2$likelihood)
  print(short.result$theta)
  print(short.result2$theta)
  print(length(y.sample))
  
  # 2. Check if EtaIndep is close
  sourceCpp("cppEtaIndep.cpp")
  eta1 <- EtaIndep(y.sample, y.lagged, z.dependent, z.independent,
                  theta$beta, theta$mu, theta$sigma, 
                  theta$gamma.dependent, theta$gamma.independent)
  eta2 <- EtaIndepR(y.sample, y.lagged, z.dependent, z.independent,
                    theta$beta, theta$mu, theta$sigma, 
                    theta$gamma.dependent, theta$gamma.independent)
  print(max(eta1-eta2))
  
  
  # 3. Check if FilterIndep is close
  sourceCpp("cppFilterIndep.cpp")
  sourceCpp("cppFilterIndepOld.cpp")
  filter1 <- FilterIndep(y.sample, y.lagged, z.dependent, z.independent,
                   theta$beta, theta$mu, theta$sigma, 
                   theta$gamma.dependent, theta$gamma.independent,
                   theta$transition.probs, theta$initial.dist)
  filter2 <- FilterIndepR(y.sample, y.lagged, z.dependent, z.independent,
                          theta$beta, theta$mu, theta$sigma, 
                          theta$gamma.dependent, theta$gamma.independent,
                          theta$transition.probs, theta$initial.dist)
  filter3 <- FilterIndepOld(y.sample, y.lagged, z.dependent, z.independent,
                         theta$beta, theta$mu, theta$sigma, 
                         theta$gamma.dependent, theta$gamma.independent,
                         theta$transition.probs, theta$initial.dist)
  print(max(filter1$likelihood-filter2$likelihood))
  print(max(filter1$xi.k-filter2$xi.k))
  print(sum(abs(filter1$xi.k - filter2$xi.k)))
  print(sum(abs(filter1$xi.k - filter3$xi.k)))
  # 4. Check if Smooth is close
  sourceCpp("cppSmooth.cpp")
  smooth1 <- Smooth(filter1$xi.k, theta$transition.probs)
  smooth2 <- SmoothR(filter1$xi.k, theta$transition.probs)
  print(max(smooth1-smooth2))
}

# model specification
n <- 50
p12 <- 0.6
p21 <- 0.7
mu1 <- -2
mu2 <- 2
sigma1 <- 1
sigma2 <- 2
beta <- c(0.4)
M <- 2
s <- 1
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# generates data
set.seed(123456)
y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21)
result.rMRS <- MLEIndepCPPComparison(y, M = M, s = s) #rMRS
result.rMRS
model=lm(y ~ 1)
msmFit(model, k=M, p=s, sw=c(T,F,T)) # MSwM (dependent mu, independent beta1, dependent sigma)

## 2. M = 2, s = 1, p.indep = 1
n <- 100
p12 <- 0.6
p21 <- 0.7
mu1 <- -2
mu2 <- 2
sigma1 <- 1
sigma2 <- 2
beta <- 0.8
M <- 2
s <- 1
p.indep <- 1
gamma.independent <- 0.7

# generates data
z.independent <- GenerateExo(n, p.indep) 

y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21, 
                      z.independent = z.independent, gamma.independent = gamma.independent)

# comparison
result.rMRS <- MLEIndepCPPComparison(y, z = z.independent, z.is.switching = FALSE, M = M, s = s, 
                                     maxit = 2) #rMRS
result.rMRS$theta
result.rMRS$likelihood
model=lm(y ~ z.independent)
msmFit(model, k=M, p=s, sw=c(T,F,F,T)) # MSwM (dependent mu, independent beta1 beta2, dependent sigma)
