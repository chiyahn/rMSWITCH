# Given short.thetas, a list of thetas, process EM algorithm to each theta
# and return them as a list with the same length.
# If gamma.dependent or gamma.independent = NULL, change them to appropriate 
# zero vectors since NULL values will make cpp codes throw exceptions.
MaximizeShortStep <- function(short.thetas, 
                              y, y.lagged, z.dependent, z.independent,
                              maxit, epsilon)
{
  n <- length(y)
  
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,n))

  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,n))

  return (lapply(short.thetas, ExpectationMaximizationIndep,
                 y = y, y.lagged = y.lagged,
                 z.dependent = z.dependent, z.independent = z.independent,
                 maxit = maxit, epsilon = epsilon))
}

ExpectationMaximizationIndep <- function (theta, y, y.lagged,
                                          z.dependent, z.independent, 
                                          maxit, epsilon)
{
  # TODO: Remove the following five lines if there is no function calling this
  # other than MaximizeShortStep.
  n <- length(y)
  if (is.null(z.dependent))
    z.dependent <- as.matrix(rep(0,n))

  if (is.null(z.independent))
    z.independent <- as.matrix(rep(0,n))
  
  EMIndepCPP(y, y.lagged, z.dependent, z.independent,
             theta$beta, theta$mu, theta$sigma,
             theta$gamma.dependent, theta$gamma.independent,
             theta$transition.probs, theta$initial.dist,
             maxit, epsilon)
}
