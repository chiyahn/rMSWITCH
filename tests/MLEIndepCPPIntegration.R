
#install.packages("MswM")
library(MSwM)
library(normalregMix)
library(rMSWITCH)
library(Rcpp)
library(RcppArmadillo)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
  lagged.and.sample <- GetLaggedAndSample(y, s)
  y.lagged <- lagged.and.sample$y.lagged
  y.sample <- lagged.and.sample$y.sample
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
  theta <- theta.initial
  

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
  filter.cpp <- FilterIndep(y.sample, y.lagged, z.dependent, z.independent,
                   theta$beta, theta$mu, theta$sigma,
                   theta$gamma.dependent, theta$gamma.independent,
                   theta$transition.probs, theta$initial.dist)
  filter.R <- FilterIndepR(y.sample, y.lagged, z.dependent, z.independent,
                          theta$beta, theta$mu, theta$sigma,
                          theta$gamma.dependent, theta$gamma.independent,
                          theta$transition.probs, theta$initial.dist)
  filter.cpp.old <- FilterIndepOld(y.sample, y.lagged, z.dependent, z.independent,
                            theta$beta, theta$mu, theta$sigma,
                            theta$gamma.dependent, theta$gamma.independent,
                            theta$transition.probs, theta$initial.dist)
  print(theta)
  print(max(filter.cpp$likelihood-filter.R$likelihood))
  print(max(filter.cpp$xi.k-filter.R$xi.k))
  print(sum(abs(filter.cpp$xi.k - filter.R$xi.k)))
  print(sum(abs(filter.cpp$xi.k - filter.cpp.old$xi.k)))
  # 4. Check if Smooth is close
  sourceCpp("cppSmooth.cpp")
  smooth.cpp <- Smooth(filter.cpp$xi.k, filter.cpp$xi.past.t, theta$transition.probs)
  smooth.R <- SmoothR(filter.cpp$xi.k, theta$transition.probs)
  print(max(smooth.cpp-smooth.R))
  return (list(filter.cpp = filter.cpp, filter.R = filter.R,
               smooth.cpp = smooth.cpp, smooth.R = smooth.R))
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
theta <- list(beta = beta, mu = c(mu1, mu2), sigma = c(sigma1, sigma2),
              transition.probs = matrix(c((1-p12), p21, p12, (1-p21)), ncol = 2),
              initial.dist = c(1,0))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# generates data
set.seed(123456)
y <- GenerateSample(n = n, theta = theta)
msar.model <- MLEIndepCPPComparison(y, M = M, s = s) #rMRS
print(msar.model$filter.R$xi.k)
apply(msar.model$filter.R$xi.k, 1, function(row) sum(row))
print(msar.model$filter.cpp$xi.k)
apply(msar.model$filter.cpp$xi.k, 1, function(row) sum(row))
print(msar.model$smooth.R)
apply(msar.model$smooth.R, 1, function(row) sum(row))
print(msar.model$smooth.cpp)
apply(msar.model$smooth.cpp, 1, function(row) sum(row))
beepr::beep(2)
model=lm(y ~ 1)
msmFit(model, k=M, p=s, sw=c(T,F,T)) # MSwM (dependent mu, independent beta1, dependent sigma)
# 
# ## 2. M = 2, s = 1, p.indep = 1
# n <- 100
# p12 <- 0.6
# p21 <- 0.7
# mu1 <- -2
# mu2 <- 2
# sigma1 <- 1
# sigma2 <- 2
# beta <- 0.8
# M <- 2
# s <- 1
# p.indep <- 1
# gamma.independent <- 0.7
# 
# # generates data
# z.independent <- GenerateExo(n, p.indep)
# 
# y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21,
#                       z.independent = z.independent, gamma.independent = gamma.independent)
# 
# # comparison
# msar.model <- MLEIndepCPPComparison(y, z = z.independent, z.is.switching = FALSE, M = M, s = s,
#                                      maxit = 2) #rMRS
# msar.model$theta
# msar.model$likelihood
# model=lm(y ~ z.independent)
# msmFit(model, k=M, p=s, sw=c(T,F,F,T)) # MSwM (dependent mu, independent beta1 beta2, dependent sigma)
