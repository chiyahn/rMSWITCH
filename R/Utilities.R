# given msar instance, draw a diagnosis plot that consists of states probabilities
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
  posterior.probs <- msar.model$posterior.probs

  # 2-1. generate a plot of posterior probabilities
  df.posterior.probs <- as.data.frame(cbind(posterior.probs, k = seq(1:n)))
  colnames(df.posterior.probs) <- c(sapply(seq(1:M), function (x) paste ("State", x)),
                                    "k")


  molten.posterior.probs <- melt(df.posterior.probs, id="k",
                                 value.name="posterior.probability",
                                 variable.name="State")
  print(ggplot(data=molten.posterior.probs,
               aes(x=k, y=posterior.probability, colour=State)) +
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

#' Returns a valid theta that represents parameters of a random MS-AR model
#' Ideal for creating a sample for testing.
#' @export
#' @title RandomTheta
#' @name RandomTheta
#' @param M The number of states in the model for the null hypothesis
#' @param s The number of terms used for AR(s)
#' @param nbtsp The number of bootstrap observations; by default, it is set to be 199
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @return A list of class \code{modelMSwitch} with items:
#' \item{theta}{The parameter estimates as a list containing alpha, mu, and sigma (and gamma if z is included in the model).}
#' \item{loglik}{The maximized value of the log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations}
#' \item{states}{n by 1 vector of integers that indicates the indices of components
#' each observation belongs to based on computed posterior probabilities}
#' \item{m}{The number of components in the mixture.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' normalmixMEMtest(y = eruptions, m = 1)
#' normalmixMEMtest(y = eruptions, m = 2)
RandomTheta <- function(M = 2, s = 1, p.dep = 0, p.indep = 0)
{
  transition.probs <- matrix(runif(M*M), ncol = M)
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
               gamma.dependent = gamma.dependent, gamma.independent = gamma.independent))
}

# Given M by n matrix of posterior.probs, returns an n by 1 integer vector that
# has the index of the regime based on posterior.probs for each observation.
EstimateStates <- function(posterior.probs)
{
  posterior.probs <- as.matrix(posterior.probs)
  # for each row, returns an indice of a coulmn (state) that is max in the row
  apply(posterior.probs, 1, function(i) (which(i==max(i))))
}
