test.on <- FALSE
test.seed <- 8888577

#' Turns on/off the test mode.
#' @export
#' @title testMode
#' @name testMode
#' @description When the modified EM-algorithm is run, initial values are randomly created
#' based on the data given. If the test mode is turned on, these initial values
#' are going to be created with the random seed provided. This method would be suitable
#' for users who would like to replicate experiments. By default, the test mode is turned off.
#' @param on Option to turn on the test mode
#' @param seed The random seed to be used for initialization
#' @param hide.message Determines whether to print the current seed and status
testMode <- function(on = FALSE, seed = 8888577, hide.message = TRUE)
{
  unlockBinding("test.on", getNamespace("normalregMix"))
  unlockBinding("test.seed", getNamespace("normalregMix"))
  assign("test.on", on, getNamespace("normalregMix"))
  assign("test.seed", seed, getNamespace("normalregMix"))
  lockBinding("test.on", getNamespace("normalregMix"))
  lockBinding("test.seed", getNamespace("normalregMix"))

  if (!hide.message)
    print(paste("The test mode is currently",
                switch(as.character(test.on), "TRUE" = "ON", "FALSE" = "OFF"),
                "with seed",
                as.character(test.seed)))
}

#' Draws two diagnosis plots that visually present 1. changes in posterior
#' probabilities for each regime across observations and 2. estimated regimes
#' based on posterior probabilities.
#' @export
#' @title DiagPlot
#' @name DiagPlot
#' @param msar.model An instance in msar.model represents one MSM-AR/MSI-AR model
#' @param y n by 1 column that represents a time series
#' @examples
#' theta <- RandomTheta(M = 2, s = 3)
#' sample.meta <- GenerateSample(theta)
#' y <- sample.meta$y
#' msar.model <- sample.meta$msar.model
#' DiagPlot(msar.model, y)
DiagPlot <- function(msar.model, y)
{
  library(ggplot2)
  library(reshape2)
  s <- nrow(as.matrix(msar.model$theta$beta))
  y <- as.matrix(y)
  y.original <- y
  y <- y[(s+1):(length(y)),]
  n <- length(y)
  M <- ncol((msar.model$theta)$transition.probs)


  states <- msar.model$states
  posterior.probs.smoothed <- msar.model$posterior.probs.smoothed # use smoothed probabilities

  # 2-1. generate a plot of posterior probabilities (smoothed)
  df.posterior.probs.smoothed <- as.data.frame(cbind(posterior.probs.smoothed, k = seq(1:n)))
  colnames(df.posterior.probs.smoothed) <- c(sapply(seq(1:M), function (x) paste ("State", x)),
                                             "k")


  molten.posterior.probs.smoothed <- melt(df.posterior.probs.smoothed, id="k",
                                          value.name="posterior.probability.smoothed",
                                          variable.name="State")
  print(ggplot(data=molten.posterior.probs.smoothed,
               aes(x=k, y=posterior.probability.smoothed, colour=State)) +
          geom_line())

  # 2-2. generate a plot of y data and shade a region for each stte
  plot.states.beginnings <- c(1)
  plot.states.endings <- vector()
  plot.states.values <- states[1]
  for (k in 2:n)
    if (plot.states.values[length(plot.states.values)] != states[k])
    {
      if (k != n)
      {
        plot.states.beginnings <- c(plot.states.beginnings, k)
        plot.states.values <- c(plot.states.values, states[k])
      }
      plot.states.endings <- c(plot.states.endings, k)

    }
  if (states[(n-1)] == states[n])
    plot.states.endings <- c(plot.states.endings, n)


  df.states <- as.data.frame(cbind(y = y.original, k = seq(0,(length(y.original)-1))))
  colnames(df.states) <- c("y", "k")
  ymin <- min(y)
  ymax <- max(y)
  # setting both y1 & y2 -Inf/Inf does not work.
  plot.states.shades.df <- data.frame(beginnings = plot.states.beginnings,
                                      endings = plot.states.endings,
                                      States = plot.states.values,
                                      y1 = -Inf,
                                      y2 = abs(max(y)) * 4)

  print(ggplot() +
    geom_rect(data=plot.states.shades.df,
                      mapping=aes(xmin=plot.states.beginnings,
                                  xmax=plot.states.endings,
                                  ymin=y1,
                                  ymax=y2, fill=factor(States)),
                      color="gray", alpha=0.3) +
    geom_line(aes(x=k, y=y), color='black',data=df.states) +
      coord_cartesian(ylim = c(-abs(min(y)) * 1.2, abs(max(y)) * 1.2)))

}

# Given theta that represents a list of parameters, return the number of
# parameters that can be used to calculate AIC/BIC.
NumberOfParameters <- function(theta)
{
  M <- ncol(theta$transition.probs)
  count <- 0

  # 1. from transition.probs
  count <- M * (M-1) + count

  # 2. from initial.dist
  count <- (M-1) + count

  # 3. from beta, mu, sigma, gamma.dependent, gamma.independent
  count <- length(theta$beta) + length(theta$mu) + length(theta$sigma) +
    length(theta$gamma.dependent) + length(theta$gamma.independent) +
    count

  return (count)
}

# Given a list of parameters (theta) and data, return n by M matrix
# that contains posterior probabilities for each observation.
EstimatePosteriorProbs <- function(theta, y, y.lagged,
                                   z.dependent, z.independent, is.MSM = FALSE)
{
  if (is.MSM)
    stop ("MSM models are currently not supported.")

  M <- ncol(theta$transition.probs)
  n <- length(y)
  is.beta.switching <- (ncol(as.matrix(theta$beta)) > 1)
  is.sigma.switching <- (length(theta$sigma) > 1)

  beta <- theta$beta
  sigma <- theta$sigma
  gamma.dependent <- matrix(rep(0, M), ncol = M)
  gamma.independent <- as.matrix(0)

  # even if beta/sigma is not switching, make it like a switching parameter
  # by creating a matrix with a duplicated column so that we can use a single
  # code to estimate posterior probabilities.
  if (!is.beta.switching)
    beta <- matrix(replicate(M, beta), ncol = M)
  if (!is.sigma.switching)
    sigma <- replicate(M, sigma)

  if (!is.null(z.dependent))
  {
    z.dependent <- as.matrix(z.dependent)
    gamma.dependent <- theta$gamma.dependent
  }
  else
    z.dependent <- as.matrix(rep(0,n))

  if (!is.null(z.independent))
  {
    z.independent <- as.matrix(z.independent)
    gamma.independent <- theta$gamma.independent
  }
  else
    z.independent <- as.matrix(rep(0,n))


  return (PosteriorProbsMSAR(y, y.lagged, z.dependent, z.independent,
                             theta$transition.probs,
                             theta$initial.dist,
                             beta,
                             theta$mu,
                             sigma,
                             gamma.dependent,
                             gamma.independent))
}

#' Returns a valid theta that represents parameters of a random MS-AR model;
#' ideal for creating a sample for testing.
#' @export
#' @title RandomTheta
#' @name RandomTheta
#' @param M The number of states in the model for the null hypothesis
#' @param s The number of terms used for AR(s)
#' @param p.dep The dimension of exogeneous variables that are switching
#' @param p.indep The dimension of exogeneous variables that are non-switching
#' @return A list that represents the parameters of a model with items:
#' \item{transition.probs}{M by M matrix that contains transition probabilities}
#' \item{initial.dist}{M by 1 column that represents an initial distribution}
#' \item{beta}{s by 1 column for state-independent coefficients on AR(s)}
#' \item{mu}{M by 1 column that contains state-dependent mu}
#' \item{sigma}{M by 1 column that contains state-dependent sigma}
#' \item{gamma.dependent}{p_dep by M matrix that contains switching
#' coefficients for state-dependent exogenous variables}
#' \item{gamma.independent}{p_indep by 1 column that contains non-switching
#' coefficients for state-independent exogenous variables}
#' @examples
#' RandomTheta()
#' RandomTheta(M = 3, s = 2)
#' RandomTheta(M = 3, s = 3, p.dep = 1)
RandomTheta <- function(M = 2, s = 1, p.dep = 0, p.indep = 0)
{
  transition.probs <- matrix(runif(M*M, 0.3, 0.5), ncol = M)
  diag(transition.probs) <- 2 # stay in your current state longer.
  transition.probs <- t(apply(transition.probs, 1,
                              function(row) (row / sum(row))))
  initial.dist <- c(0.9, rep(0.1/(M-1), (M-1)))
  beta <- runif(s, -0.4, 0.4)
  beta <- runif(1, 0.4, 0.7) * (beta / sum(abs(beta)))
  mu <- runif(M, -0.5, 0.5)
  sigma <- runif(M, 0.2, 3)

  gamma.dependent <- NULL
  gamma.independent <- NULL
  if (p.dep > 0)
  {
    gamma.dependent <- matrix(runif(p.dep * M, -0.4, 0.4), ncol = M)
    gamma.dependent <- 0.8 * (gamma.dependent / sum(abs(gamma.dependent)))
  }
  if (p.indep > 0)
  {
    gamma.independent <- runif(p.indep, -0.4, 0.4)
    gamma.independent <- 0.8 * (gamma.independent / sum(abs(gamma.independent)))
  }
  return (list(transition.probs = transition.probs,
               initial.dist = initial.dist,
               beta = beta, mu = mu, sigma = sigma,
               gamma.dependent = gamma.dependent,
               gamma.independent = gamma.independent))
}

# Given M by n matrix of posterior.probs, returns an n by 1 integer vector that
# has the index of the regime based on posterior.probs for each observation.
EstimateStates <- function(posterior.probs)
{
  posterior.probs <- as.matrix(posterior.probs)
  # for each row, returns an indice of a coulmn (state) that is max in the row
  apply(posterior.probs, 1, function(i) (which(i==max(i))))
}

#' Order a transition matrix given an order of a column.
#' Example:
#' mu <- c(1,-1)
#' mu.order <- order(mu)   # 2, 1
#' transition.matrix <- matrix(c(0.4,0.7,0.6,0.3), ncol = 2)
#' OrderTransitionMatrix(transition.matrix, mu.order)
OrderTransitionMatrix <- function(transition.matrix, new.order)
{
  M <- ncol(transition.matrix)
  ordered.matrix <- matrix(rep(0,(M*M)), ncol = M)
  for (i in 1:M)
    for (j in 1:M)
      ordered.matrix[i,j] <- transition.matrix[new.order[i],new.order[j]]
  return (ordered.matrix)
}

#' Transform a theta to a reduced column, ordered in
#' transition.probs, initial.dist, beta, mu, sigma, gamma.dep, gamma.indep.
#' where each row of transition.probs is reduced to have a length of M-1,
#' initial.dist is reduced to have a length of M-1.
ThetaToReducedColumn <- function(theta)
{
  # transition.probs could have been listed in the order of
  # p11, p21, ..., pM1, p12, ..., pM2, ..., p1(M-1), ..., pM(M-1)
  # taking a transpose will make it listed as
  # p11, p12, ..., p1(M-1), p21, ..., p2(M-1), ..., pM(M-1)
  M <- ncol(theta$transition.probs)
  reduced.transition.probs <- theta$transition.probs[,1:(M-1)]
  reduced.initial.dist <- theta$initial.dist[1:(M-1)]
  return (c(c(t(reduced.transition.probs)),
            c(reduced.initial.dist),
            c(theta$beta), c(theta$mu), c(theta$sigma),
            c(theta$gamma.dependent),
            c(theta$gamma.independent)))

}

#' Given a |theta.reduced| by n matrix, whose column represents a reduced form of
#' parameters and a theta0 that retains information about parameters,
#' return an equivalent list of thetas (theta is in a list).
ReducedColumnsToThetas <- function(theta.matrix, theta0)
{
  M <- ncol(theta0$transition.probs)
  s <- nrow(as.matrix(theta0$beta))
  is.beta.switching <- (ncol(as.matrix(theta0$beta)) > 1)
  is.sigma.switching <- (length(theta0$sigma) > 1)
  p.dep <- 0
  p.indep <- 0
  if (!is.null(theta0$gamma.dependent))
    p.dep <- nrow(as.matrix(theta0$gamma.dependent))
  if (!is.null(theta0$gamma.independent))
    p.indep <- nrow(as.matrix(theta0$gamma.independent))

  # this holds only for univariate time series.
  initial.dist.index <- M * (M-1) + 1 # reduced case
  beta.index <- (M-1) + initial.dist.index # reduced case
  mu.index <- s * ifelse(is.beta.switching, M, 1) + beta.index
  sigma.index <- M + mu.index
  gamma.dep.index <- ifelse(is.sigma.switching, M, 1) + sigma.index
  gamma.indep.index <- p.dep * M + gamma.dep.index


  ReducedColumnToThetaDynamic <- function (theta.vectorized)
  {
    transition.probs <- matrix(theta.vectorized[1:(M*(M-1))],
                                ncol = (M-1), byrow = T)
    # revive the original from the reduced form.
    transition.probs <- t(apply(transition.probs, 1,
                                function (row) c(row, (1-sum(row)))))
    initial.dist <- theta.vectorized[initial.dist.index:(beta.index - 1)]
    initial.dist <- c(initial.dist, (1 - sum(initial.dist)))
    beta <- theta.vectorized[beta.index:(mu.index - 1)]
    if (is.beta.switching)
      beta <- matrix(beta, ncol = M)
    gamma.dependent <- NULL
    gamma.independent <- NULL
    if (!gamma.dep.index == gamma.indep.index)
      gamma.dependent <- matrix(theta.vectorized[gamma.dep.index:
                                (gamma.indep.index - 1)], ncol = M)
    if (!gamma.indep.index == length(theta.vectorized) + 1)
      gamma.independent <- theta.vectorized[gamma.indep.index:
                                              length(theta.vectorized)]

    return (list
            (transition.probs = transition.probs,
            initial.dist = initial.dist,
            beta = beta,
            mu = theta.vectorized[mu.index:(sigma.index - 1)],
            sigma = theta.vectorized[sigma.index:(gamma.dep.index - 1)],
            gamma.dependent = gamma.dependent,
            gamma.independent = gamma.independent
            ))
  }

  return (apply(theta.matrix, 2, ReducedColumnToThetaDynamic))
}

#' Given a |theta.reduced| by n matrix, whose column represents a reduced form of
#' parameters from Stata and a theta0 that retains information about parameters,
#' return an equivalent list of thetas (theta is in a list).
#' Stata saves parameters in the order of
#' beta, gamma.dependent, gamma.independent, mu, sigma, transition.probs.
ReducedStataColumnsToReducedColumns <- function(theta.matrix.stata, theta0)
{
  M <- ncol(theta0$transition.probs)
  s <- nrow(as.matrix(theta0$beta))
  is.beta.switching <- (ncol(as.matrix(theta0$beta)) > 1)
  is.sigma.switching <- (length(theta0$sigma) > 1)
  p.dep <- 0
  p.indep <- 0
  if (!is.null(theta0$gamma.dependent))
    p.dep <- nrow(as.matrix(theta0$gamma.dependent))
  if (!is.null(theta0$gamma.independent))
    p.indep <- nrow(as.matrix(theta0$gamma.independent))
  no.gamma <- (p.dep + p.indep == 0)

  # this holds only for univariate time series.
  beta.end.index <- ifelse(is.beta.switching, (M * s), s)
  gamma.dep.index <- beta.end.index + 1
  gamma.indep.index <- p.dep * M + (beta.end.index + 1)
  gamma.indep.end.index <- p.indep + p.dep * M + beta.end.index
  mu.index <- gamma.indep.end.index + 1
  transition.probs.index <- nrow(theta.matrix.stata) - (M * (M - 1)) + 1

  StataColumnToReducedThetaColumn <- function (stata.column)
  {
    transition.probs.column <- stata.column[transition.probs.index:
                                            length(stata.column)]
    mu.sigma.column <- stata.column[mu.index:(transition.probs.index - 1)]
    beta.column <- stata.column[1:beta.end.index]
    gamma.dep.indep.column <- NULL
    if (!no.gamma)
      gamma.dep.indep.column <- stata.column[gamma.dep.index:
                                            gamma.indep.end.index]

    return (c(transition.probs.column,
              rep(-Inf, (M-1)),
              beta.column,
              mu.sigma.column,
              gamma.dep.indep.column))
  }

  theta.matrix <- apply(theta.matrix.stata, 2, StataColumnToReducedThetaColumn)
  return (theta.matrix)
}
