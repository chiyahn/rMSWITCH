
## 1. M = 2, s = 1
# model specification
set.seed(123456)
n <- 400
transition.probs <- matrix(c(0.8,0.15,0.2,0.85), ncol = 2)
beta <- matrix(c(0.9,0.4), ncol = 2)
mu = c(-0.7,0.7)
sigma = c(0.6)
theta <- list(beta = beta, mu = mu, sigma = sigma,
              transition.probs = transition.probs, initial.dist = c(1,0))
M <- ncol(transition.probs)
s <- nrow(as.matrix(beta))
is.beta.switching <- TRUE
is.sigma.switching <- FALSE
p.dependent <- 0
p.independent <- 0

sample <- GenerateSample(theta, n = n, is.MSM = TRUE)
epsilon <- 1e-08
maxit <- 2000
short.iterations = 200
short.epsilon <- 1e-03
transition.probs.min <- 0.02
transition.probs.max <- 0.98

# data
y <- sample$y
y.lagged <- sample$y.lagged
y.sample <- sample$y.sample
z.dependent <- as.matrix(rep(0,n))
z.independent <- as.matrix(rep(0,n))
n <- length(y.sample)
initial.theta <- RandomTheta(is.sigma.switching = FALSE, is.beta.switching = TRUE, is.MSM = TRUE)
initial.theta$beta <- as.matrix(initial.theta$beta)
initial.theta$gamma.dependent <- t(as.matrix(rep(0, 2)))
initial.theta$gamma.independent <- as.matrix(rep(0, 1))
thetas <- list(initial.theta)

library(nloptr)
library(rMSWITCH)
result.nloptr <- MaximizeLongStepNLOPTR(thetas,
                                  y = y.sample, y.lagged = y.lagged,
                                  z.dependent = z.dependent,
                                  z.independent = z.independent,
                                  is.beta.switching = is.beta.switching,
                                  is.sigma.switching = is.sigma.switching, maxit = maxit,
                                  epsilon = epsilon,
                                  transition.probs.min = transition.probs.min,
                                  transition.probs.max = transition.probs.max,
                                  sigma.min = 0.05,
                                  is.MSM = TRUE)
result.nloptr