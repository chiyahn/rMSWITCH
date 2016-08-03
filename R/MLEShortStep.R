# Given short.thetas, a list of thetas, process EM algorithm to each theta
# and return them as a list with the same length.
# If gamma.dependent or gamma.independent = NULL, change them to appropriate 
# zero vectors since NULL values will make cpp codes throw exceptions.
MaximizeShortStep <- function(short.thetas, 
                              y, y.lagged, z.dependent, z.independent,
                              is.beta.switching,
                              is.sigma.switching,
                              maxit, epsilon)
{
  n <- length(y)
  
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,n))

  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,n))
  
  if (is.sigma.switching)
    if (is.beta.switching)
      return (lapply(short.thetas, EM.MSMAH.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon))
    else
      return (lapply(short.thetas, EM.MSMH.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon))
  else
    return (NULL)
}

EM.MSMAH.AR <- function (theta, y, y.lagged,
                         z.dependent, z.independent, 
                         maxit, epsilon)
{
  
  EMcppARMSMAH(y, y.lagged, z.dependent, z.independent,
              theta$beta, theta$mu, theta$sigma,
              theta$gamma.dependent, theta$gamma.independent,
              theta$transition.probs, theta$initial.dist,
              maxit, epsilon)
}

EM.MSMH.AR <- function (theta, y, y.lagged,
                       z.dependent, z.independent, 
                       maxit, epsilon)
{

  EMcppARMSMH(y, y.lagged, z.dependent, z.independent,
             theta$beta, theta$mu, theta$sigma,
             theta$gamma.dependent, theta$gamma.independent,
             theta$transition.probs, theta$initial.dist,
             maxit, epsilon)
}
