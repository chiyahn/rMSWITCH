###########################################################
# Compares the result with MSwM
###########################################################

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
                             gamma.dependent = matrix(rep(0,2), ncol = 2),
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
                             gamma.dependent = matrix(rep(0,3), ncol = 3),
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

# comparison
result.rMRS <- MRSMLEIndep(y, M = M, s = s) #rMRS
result.rMRS$theta
result.rMRS$likelihood
model=lm(y ~ 1)
msmFit(model, k=M, p=s, sw=c(T,F,T)) # MSwM (dependent mu, independent beta, dependent sigma)

## 2. M = 2, s = 2
beta <- c(0.4, 0.5)
s <- 2

# generates data
y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21)

# comparison
result.rMRS <- MRSMLEIndep(y, M = M, s = s) #rMRS
result.rMRS$theta
result.rMRS$likelihood
model=lm(y ~ 1)
msmFit(model, k=M, p=s, sw=c(T,F,F,T)) # MSwM (dependent mu, independent beta1 beta2, dependent sigma)

## 3. M = 3, s = 2
# model specification
set.seed(654321)
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
M <- 3
s <- length(beta)

# generates data
y <- GenerateSampleM3(n, beta, mu1, mu2, mu3,
                      sigma1, sigma2, sigma3,
                      p12, p13, p21, p23, p31, p32)


# comparison
result.rMRS <- MRSMLEIndep(y, M = M, s = s) #rMRS
result.rMRS$theta
result.rMRS$likelihood
model=lm(y ~ 1)
msmFit(model, k=M, p=s, sw=c(T,F,F,T)) # MSwM (dependent mu, independent beta1 beta2, dependent sigma)

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
z.independent <- GenerateExo(n, p.indep)

y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21,
                      z.independent = z.independent, gamma.independent = gamma.independent)

# comparison
result.rMRS <- MRSMLEIndep(y, z = z.independent, z.is.switching = FALSE, M = M, s = s) #rMRS
result.rMRS$theta
result.rMRS$likelihood
model=lm(y ~ z.independent)
msmFit(model, k=M, p=s, sw=c(T,F,F,T)) # MSwM (dependent mu, independent beta1 beta2, dependent sigma)

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
z.dependent <- GenerateExo(n, p.indep)
y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21,
                      z.dependent = z.dependent, gamma.dependent = gamma.dependent)

# comparison
result.rMRS <- MRSMLEIndep(y, z = z.dependent, z.is.switching = TRUE, M = M, s = s) #rMRS
result.rMRS$theta
result.rMRS$likelihood
model=lm(y ~ z.dependent)
msmFit(model, k=M, p=s, sw=c(T,T,F,T)) # MSwM (dependent mu, z, independent beta1, dependent sigma)

## 6. M = 2, s = 2, p.dep = 2, p.indep = 2
# model specification
n <- 100
p12 <- 0.6
p21 <- 0.7
mu1 <- -2
mu2 <- 2
sigma1 <- 1
sigma2 <- 2
beta <- c(0.3, 0.5)
M <- 2
s <- 2
p.dep <- 2
p.indep <- 2
gamma.dependent <- matrix(c(0.3,0.6,0.6,0.3), ncol = 2)
gamma.independent <- matrix(c(0.4,0.5), ncol = 1)

# generates data
z.dependent <- z[,1:p.dep]
z.independent <- z[,(p.dep+1):(p.dep+p.indep)]
y <- GenerateSampleM2(n, beta, mu1, mu2, sigma1, sigma2, p12, p21,
                      z.dependent = z.dependent, gamma.dependent = gamma.dependent,
                      z.independent = z.independent, gamma.independent = gamma.independent)
z <- GenerateExo(n, p=(p.dep + p.indep))

# comparison
result.rMRS <- MRSMLEIndep(y, z = z,  z.is.switching = c(T,T,F,F), M = M, s = s) #rMRS
result.rMRS$theta
result.rMRS$likelihood
model=lm(y ~ z.independent + z.dependent)
msmFit(model, k=M, p=s, sw=c(T,T,T,F,F,F,F,T)) # MSwM (dependent mu, z1, z2, independent z3, z4, beta1, beta2, dependent sigma)


## 7. M = 3, s = 2, p.dep = 3, p.indep = 2
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
y <- (y / 10)[3:nrow(z.independent)]
z  <- z[3:nrow(z.independent),]
z.dependent <- z.dependent[3:nrow(z.independent),]
z.independent <- z.independent[3:nrow(z.independent),]



# comparison
result.rMRS <- MRSMLEIndep(y, z = z, z.is.switching = c(T,T,T,F,F), M = M, s = s) #rMRS
result.rMRS$theta
result.rMRS$likelihoods
model=lm(y ~ z.independent + z.dependent)
msmFit(model, k=M, p=s, sw=c(T,T,T,T,F,F,F,F,T))  # MSwM (dependent mu, z1, z2, z3, independent z4, z5, beta1, beta2, dependent sigma)
