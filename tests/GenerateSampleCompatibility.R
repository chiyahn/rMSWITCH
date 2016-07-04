library(MSwM)
library(normalregMix)
library(rMSWITCH)
library(Rcpp)
library(RcppArmadillo)
library(testthat)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Generates exogeneous variable sample
GenerateExo <- function(n, p)
{
  # changes wd; only works running in RStudio
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  return (matrix(as.numeric(read.csv("interestRates.csv")[1:(n*p),1]), ncol = p))
}

# Get elementwise-average of list of matrices
GetMatrixAvg <- function (list.of.matrices)
{
  Reduce("+", list.of.matrices) / length(list.of.matrices)
}

## 1. M = 2, s = 1
# model specification
n <- 50
p12 <- 0.6
p21 <- 0.7
mu1 <- -2
mu2 <- 2
sigma1 <- 1
sigma2 <- 2
beta <- 0.8
M <- 2
s <- 1

# generates data
set.seed(123456)
y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21)
set.seed(123456)
transition.probs <- matrix(c((1-p12),p12,p21,(1-p21)), ncol = 2, byrow = T)
mu = c(mu1,mu2)
sigma = c(sigma1,sigma2)
theta <- list(beta = beta, mu = mu, sigma = sigma, 
              transition.probs = transition.probs, initial.dist = c(0,1))
y2 <- GenerateSample(theta = theta, n = n, initial.y.set = rnorm(length(theta$beta)), 
                     initial.state = 1)$sample
expect_that(y, equals(y2))


## 2. M = 2, s = 2
# model specification
beta <- c(0.4, 0.5)
s <- 2

# generates data
set.seed(123456)
y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21)
set.seed(123456)
theta$beta <- beta
y2 <- GenerateSample(theta = theta, n = n, initial.y.set = rnorm(length(theta$beta)), 
                     initial.state = 1)$sample
expect_that(y, equals(y2))


## 3. M = 3, s = 2
# model specification
n <- 200
beta <- c(0.5, 0.3)
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
M <- 3
s <- length(beta)

# generates data
set.seed(123456)
y <- GenerateSampleM3(n, beta, mu1, mu2, mu3,
                      sigma1, sigma2, sigma3,
                      p12, p13, p21, p23, p31, p32)
set.seed(123456)
transition.probs <- matrix(c((1-p12-p13),p12,p13,
                             p21,(1-p21-p23),p23,
                             p31,p32,(1-p31-p32)), 
                           ncol = 3, byrow =T)
mu = c(mu1,mu2,mu3)
sigma = c(sigma1,sigma2,sigma3)
theta <- list(beta = beta, mu = mu, sigma = sigma, 
              transition.probs = transition.probs, initial.dist = c(0,1))
y2 <- GenerateSample(theta = theta, n = n, initial.y.set = rnorm(length(theta$beta)), 
                     initial.state = 1)$sample
expect_that(y, equals(y2))


## 4. M = 2, s = 1, p.indep = 1
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
set.seed(123456)
z.independent <- GenerateExo(n, p.indep)
y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21,
                      z.independent = z.independent, 
                      gamma.independent = gamma.independent)
set.seed(123456)
mu = c(mu1,mu2)
sigma = c(sigma1,sigma2)
transition.probs <- matrix(c((1-p12),p12,p21,(1-p21)), 
                           ncol = 2, byrow = T)
theta <- list(beta = beta, mu = mu, sigma = sigma, 
              transition.probs = transition.probs, initial.dist = c(0,1),
              gamma.independent = gamma.independent)
y2 <- GenerateSample(theta = theta, n = n, initial.y.set = rnorm(length(theta$beta)), 
                     z.independent = z.independent,
                     initial.state = 1)$sample
expect_that(y, equals(y2))


## 5. M = 2, s = 1, p.dep = 1
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
p.dep <- 1
gamma.dependent <- matrix(c(0.3,0.7), ncol = 2)

# generates data
set.seed(123456)
z.dependent <- GenerateExo(n, p.dep)
y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21,
                      z.dependent = z.dependent, 
                      gamma.dependent = gamma.dependent)
set.seed(123456)
mu = c(mu1,mu2)
sigma = c(sigma1,sigma2)
transition.probs <- matrix(c((1-p12),p12,p21,(1-p21)), 
                           ncol = 2, byrow = T)
theta <- list(beta = beta, mu = mu, sigma = sigma, 
              transition.probs = transition.probs, initial.dist = c(0,1),
              gamma.dependent = gamma.dependent)
y2 <- GenerateSample(theta = theta, n = n, initial.y.set = rnorm(length(theta$beta)), 
                     z.dependent = z.dependent,
                     initial.state = 1)$sample
expect_that(y, equals(y2))


## 6. M = 3, s = 2, p.dep = 3, p.indep = 2
# model specification
n <- 140
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
M <- 3
s <- length(beta)
p.dep <- 3
p.indep <- 2
gamma.dependent <- matrix(c(0.3,0.1,0.1,0.1,0.3,0.1,0.1,0.1,0.3), ncol = 3)
gamma.independent <- matrix(c(0.2,0.5), ncol = 1)

# generates data
set.seed(123456)
z <- GenerateExo(n, p=(p.dep + p.indep))
z.dependent <- z[,1:p.dep]
z.independent <- z[,(p.dep+1):(p.dep+p.indep)]
y <- GenerateSampleM3(n, beta, mu1, mu2, mu3,
                      sigma1, sigma2, sigma3,
                      p12, p13, p21, p23, p31, p32,
                      z.dependent = z.dependent,
                      z.independent = z.independent,
                      gamma.dependent = gamma.dependent,
                      gamma.independent = gamma.independent)
set.seed(123456)
transition.probs <- matrix(c((1-p12-p13),p12,p13,
                             p21,(1-p21-p23),p23,
                             p31,p32,(1-p31-p32)), 
                           ncol = 3, byrow =T)
mu = c(mu1,mu2,mu3)
sigma = c(sigma1,sigma2,sigma3)
theta <- list(beta = beta, mu = mu, sigma = sigma, 
              transition.probs = transition.probs, initial.dist = c(0,1),
              gamma.independent = gamma.independent,
              gamma.dependent = gamma.dependent)
y2 <- GenerateSample(theta = theta, n = n, initial.y.set = rnorm(length(theta$beta)), 
                     z.dependent = z.dependent, z.independent = z.independent,
                     initial.state = 1)$sample
expect_that(y, equals(y2))

