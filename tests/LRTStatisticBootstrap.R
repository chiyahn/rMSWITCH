#################################################################
# Test the number of regimes based on bootstrapped LRT statistics
#################################################################
library(normalregMix)
library(doParallel)
library(rMSWITCH)
library(nloptr)
set.seed(123456)
n <- 100
transition.probs <- matrix(c(0.9,0.2,0.1,0.8), ncol = 2)
beta <- matrix(c(0.6,0.3), ncol = 2)
mu = c(-0.6,0.6)
sigma = c(0.6,1.4)
theta <- list(beta = beta, mu = mu, sigma = sigma, 
              transition.probs = transition.probs, 
              initial.dist = c(1,0))
M <- ncol(transition.probs)
s <- nrow(as.matrix(beta))

# generates data
sample <- GenerateSample(theta, n = n)
y <- sample$y

# comparison
set.seed(123456)
test <- TestMSAR(y = y, M = (M-1), s = 1, 
                 is.beta.switching = TRUE,
                 is.sigma.switching = TRUE,
                 crit.method = 'bootstrap',
                 parallel = TRUE)
