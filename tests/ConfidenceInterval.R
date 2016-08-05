####################################################
# Check if estimated parameters are consistent with
# confidence interval created from the estimated
# Fisher information matrix.
####################################################
## Configuration
n <- 1600
replications <- 2500

## Case 1
# model specification
library(rMSWITCH)
library(normalregMix)
library(nloptr)
library(doParallel)

set.seed(123456)
initial.y.set <- 0

transition.probs <- matrix(c(0.8,0.3,0.2,0.7), ncol = 2)
beta <- matrix(c(0.6,0.3), ncol = 2)
mu = c(-0.8,0.8)
sigma = c(0.5)
theta <- list(beta = beta, mu = mu, sigma = sigma,
              transition.probs = transition.probs,
              initial.dist = c(0.95,0.05))
M <- ncol(transition.probs)
s <- nrow(as.matrix(beta))

# generates data
sample <- GenerateSample(theta, n = n,
                         initial.y.set = initial.y.set)
y <- sample$y
model <- EstimateMSAR(y = y, z.dependent = NULL, z.independent = NULL,
                     M, s,
                     is.beta.switching = TRUE,
                     is.sigma.switching = FALSE)
fisher.estimated <- model$fisher.estimated
theta.estimated <- model$theta

samples <- GenerateSamples(theta = theta.estimated, n = n,
                           replications = replications,
                           initial.y.set = initial.y.set)
models <- SamplesToModels(samples, M, s,
                          is.beta.switching = TRUE,
                          is.sigma.switching = FALSE,
                          parallel = TRUE)
thetas <- t(sapply(models$thetas,
                 function (theta) ThetaToReducedColumn(theta)))
thetas <- apply(thetas, 2, sort)
write.csv(as.data.frame(sample$y),file="y1.csv")
write.csv(as.data.frame(thetas),file="thetas1.csv")

fisher.estimated <- EstimateFisherInformation(theta = theta.estimated,
                                              y = sample$y.sample, y.lagged = sample$y.lagged,
                                              z.dependent = NULL, z.independent = NULL)
sd <- sqrt(diag(solve(fisher.estimated)))

lb <- ThetaToReducedColumn(theta.estimated) - 1.96*sd/sqrt(n)
ub <- ThetaToReducedColumn(theta.estimated) + 1.96*sd/sqrt(n)

print(lb)
print(ThetaToReducedColumn(theta.estimated))
print(ub)

freqs <- rep(0, length(lb))
for (i in 1:length(lb))
  freqs[i] <- mean((thetas[,i] < ub[i]) * (thetas[,i] > lb[i]))
print(freqs)

## Case 2
set.seed(123456)
initial.y.set <- 0
transition.probs <- matrix(c(0.5,0.5,0.5,0.5), ncol = 2)
beta <- matrix(c(0.1,0.2), ncol = 2)
mu = c(-0.7,0.7)
sigma = c(0.5)
theta <- list(beta = beta, mu = mu, sigma = sigma,
              transition.probs = transition.probs,
              initial.dist = c(1,0))
M <- ncol(transition.probs)
s <- nrow(as.matrix(beta))

# generates data
sample <- GenerateSample(theta, n = n,
                         initial.y.set = initial.y.set)
y <- sample$y
model <- EstimateMSAR(y = y, z.dependent = NULL, z.independent = NULL,
                      M, s,
                      is.beta.switching = TRUE,
                      is.sigma.switching = FALSE)
fisher.estimated <- model$fisher.estimated
theta.estimated <- model$theta

samples <- GenerateSamples(theta = theta.estimated, n = n,
                           replications = replications,
                           initial.y.set = initial.y.set)
models <- SamplesToModels(samples, M, s,
                          is.beta.switching = TRUE,
                          is.sigma.switching = FALSE,
                          parallel = TRUE)
thetas <- t(sapply(models$thetas,
                   function (theta) ThetaToReducedColumn(theta)))
thetas <- apply(thetas, 2, sort)
write.csv(as.data.frame(sample$y),file="y2.csv")
write.csv(as.data.frame(thetas),file="thetas2.csv")

fisher.estimated <- EstimateFisherInformation(theta = theta.estimated,
                                              y = sample$y.sample, y.lagged = sample$y.lagged,
                                              z.dependent = NULL, z.independent = NULL)
sd <- sqrt(diag(solve(fisher.estimated)))

lb <- ThetaToReducedColumn(theta.estimated) - 1.96*sd/sqrt(n)
ub <- ThetaToReducedColumn(theta.estimated) + 1.96*sd/sqrt(n)

print(lb)
print(ThetaToReducedColumn(theta.estimated))
print(ub)

freqs <- rep(0, length(lb))
for (i in 1:length(lb))
  freqs[i] <- mean((thetas[,i] < ub[i]) * (thetas[,i] > lb[i]))
print(freqs)

## Case 3
set.seed(123456)
initial.y.set <- 0
transition.probs <- matrix(c(0.8,0.3,0.2,0.7), ncol = 2)
beta <- matrix(c(0.1,0.2), ncol = 2)
mu = c(-0.7,0.7)
sigma = c(0.5)
theta <- list(beta = beta, mu = mu, sigma = sigma,
              transition.probs = transition.probs,
              initial.dist = c(0.95,0.05))
M <- ncol(transition.probs)
s <- nrow(as.matrix(beta))

# generates data
sample <- GenerateSample(theta, n = n,
                         initial.y.set = initial.y.set)
y <- sample$y
model <- EstimateMSAR(y = y, z.dependent = NULL, z.independent = NULL,
                      M, s,
                      is.beta.switching = TRUE,
                      is.sigma.switching = FALSE)
fisher.estimated <- model$fisher.estimated
theta.estimated <- model$theta

samples <- GenerateSamples(theta = theta.estimated, n = n,
                           replications = replications,
                           initial.y.set = initial.y.set)
models <- SamplesToModels(samples, M, s,
                          is.beta.switching = TRUE,
                          is.sigma.switching = FALSE,
                          parallel = TRUE)
thetas <- t(sapply(models$thetas,
                   function (theta) ThetaToReducedColumn(theta)))
thetas <- apply(thetas, 2, sort)
write.csv(as.data.frame(sample$y),file="y3.csv")
write.csv(as.data.frame(thetas),file="thetas3.csv")

fisher.estimated <- EstimateFisherInformation(theta = theta.estimated,
                                              y = sample$y.sample, y.lagged = sample$y.lagged,
                                              z.dependent = NULL, z.independent = NULL)
sd <- sqrt(diag(solve(fisher.estimated)))

lb <- ThetaToReducedColumn(theta.estimated) - 1.96*sd/sqrt(n)
ub <- ThetaToReducedColumn(theta.estimated) + 1.96*sd/sqrt(n)

print(lb)
print(ThetaToReducedColumn(theta.estimated))
print(ub)

freqs <- rep(0, length(lb))
for (i in 1:length(lb))
  freqs[i] <- mean((thetas[,i] < ub[i]) * (thetas[,i] > lb[i]))
print(freqs)


## Case 4
set.seed(123456)
initial.y.set <- 0
transition.probs <- matrix(c(0.5,0.5,0.5,0.5), ncol = 2)
beta <- matrix(c(0.6,0.3), ncol = 2)
mu = c(-0.7,0.7)
sigma = c(0.5)
theta <- list(beta = beta, mu = mu, sigma = sigma,
              transition.probs = transition.probs,
              initial.dist = c(0.95,0.05))
M <- ncol(transition.probs)
s <- nrow(as.matrix(beta))

# generates data
sample <- GenerateSample(theta, n = n,
                         initial.y.set = initial.y.set)
y <- sample$y
model <- EstimateMSAR(y = y, z.dependent = NULL, z.independent = NULL,
                      M, s,
                      is.beta.switching = TRUE,
                      is.sigma.switching = FALSE)
fisher.estimated <- model$fisher.estimated
theta.estimated <- model$theta

samples <- GenerateSamples(theta = theta.estimated, n = n,
                           replications = replications,
                           initial.y.set = initial.y.set)
models <- SamplesToModels(samples, M, s,
                          is.beta.switching = TRUE,
                          is.sigma.switching = FALSE,
                          parallel = TRUE)
thetas <- t(sapply(models$thetas,
                   function (theta) ThetaToReducedColumn(theta)))
thetas <- apply(thetas, 2, sort)
write.csv(as.data.frame(sample$y),file="y4.csv")
write.csv(as.data.frame(thetas),file="thetas4.csv")

fisher.estimated <- EstimateFisherInformation(theta = theta.estimated,
                                              y = sample$y.sample, y.lagged = sample$y.lagged,
                                              z.dependent = NULL, z.independent = NULL)
sd <- sqrt(diag(solve(fisher.estimated)))

lb <- ThetaToReducedColumn(theta.estimated) - 1.96*sd/sqrt(n)
ub <- ThetaToReducedColumn(theta.estimated) + 1.96*sd/sqrt(n)

print(lb)
print(ThetaToReducedColumn(theta.estimated))
print(ub)

freqs <- rep(0, length(lb))
for (i in 1:length(lb))
  freqs[i] <- mean((thetas[,i] < ub[i]) * (thetas[,i] > lb[i]))
print(freqs)
