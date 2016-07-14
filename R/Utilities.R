# given msar instance, draw a diagnosis plot that consists of states probabilities
DiagPlot <- function(msar.model, y)
{
  library(ggplot2)
  library(reshape2)
  
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


# Given M by n matrix of posterior.probs, returns an n by 1 integer vector that
# has the index of the regime based on posterior.probs for each observation.
EstimateStates <- function(posterior.probs)
{
  posterior.probs <- as.matrix(posterior.probs)
  # for each row, returns an indice of a coulmn (state) that is max in the row
  apply(posterior.probs, 1, function(i) (which(i==max(i))))
}
