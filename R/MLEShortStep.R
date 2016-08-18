# Given short.thetas, a list of thetas, process EM algorithm to each theta
# and return them as a list with the same length.
# If gamma.dependent or gamma.independent = NULL, change them to appropriate 
# zero vectors since NULL values will make cpp codes throw exceptions.
MaximizeShortStep <- function(short.thetas, 
                              y, y.lagged, z.dependent, z.independent,
                              is.beta.switching,
                              is.sigma.switching,
                              maxit, epsilon,
                              transition.probs.min,
                              transition.probs.max)
{
  n <- length(y)
  
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,n))

  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,n))
  
  
  if (is.beta.switching)
    if (is.sigma.switching)
      return (lapply(short.thetas, EM.MSIAH.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max))
    else
      return (lapply(short.thetas, EM.MSIA.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max))
  else
    if (is.sigma.switching)
      return (lapply(short.thetas, EM.MSIH.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max))
    else
      return (lapply(short.thetas, EM.MSI.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max))
}
EM.MSIA.AR <- function (theta, y, y.lagged,
                        z.dependent, z.independent, 
                        maxit, epsilon, 
                        transition.probs.min, transition.probs.max)
{
  
  EMcppARMSIA(y, y.lagged, z.dependent, z.independent,
              theta$beta, theta$mu, theta$sigma,
              theta$gamma.dependent, theta$gamma.independent,
              theta$transition.probs, theta$initial.dist,
              maxit, epsilon,
              transition.probs.min, transition.probs.max)
}

EM.MSIAH.AR <- function (theta, y, y.lagged,
                         z.dependent, z.independent, 
                         maxit, epsilon, 
                         transition.probs.min, transition.probs.max)
{
  
  EMcppARMSIAH(y, y.lagged, z.dependent, z.independent,
              theta$beta, theta$mu, theta$sigma,
              theta$gamma.dependent, theta$gamma.independent,
              theta$transition.probs, theta$initial.dist,
              maxit, epsilon,
              transition.probs.min, transition.probs.max)
}

EM.MSI.AR <- function (theta, y, y.lagged,
                        z.dependent, z.independent, 
                        maxit, epsilon, 
                       transition.probs.min, transition.probs.max)
{
  
  EMcppARMSI(y, y.lagged, z.dependent, z.independent,
              theta$beta, theta$mu, theta$sigma,
              theta$gamma.dependent, theta$gamma.independent,
              theta$transition.probs, theta$initial.dist,
              maxit, epsilon,
              transition.probs.min, transition.probs.max)
}


EM.MSIH.AR <- function (theta, y, y.lagged,
                       z.dependent, z.independent, 
                       maxit, epsilon, 
                       transition.probs.min, transition.probs.max)
{

  EMcppARMSIH(y, y.lagged, z.dependent, z.independent,
             theta$beta, theta$mu, theta$sigma,
             theta$gamma.dependent, theta$gamma.independent,
             theta$transition.probs, theta$initial.dist,
             maxit, epsilon,
             transition.probs.min, transition.probs.max)
}
