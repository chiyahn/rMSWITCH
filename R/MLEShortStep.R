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
                              transition.probs.max,
                              sigma.min, 
                              z.dependent.lagged = NULL,
                              z.independent.lagged = NULL,
                              is.MSM = FALSE)
{
  n <- length(y)
  M <- ncol(as.matrix(short.thetas[[1]]$transition.probs))
  s <- nrow(as.matrix(short.thetas[[1]]$beta)) # WARNING: assume beta exists
  
  if (is.null(z.dependent))
  {
    z.dependent <- as.matrix(rep(0,n))
    z.dependent.lagged <- matrix(0, ncol = s, nrow = n)
  }
  if (is.null(z.independent))
  {
    z.independent <- as.matrix(rep(0,n))
    z.independent.lagged <- matrix(0, ncol = s, nrow = n)
  }
  
  if (is.MSM)
  {
    state.conversion.mat <- GetStateConversionMat(M, s)
    if (is.beta.switching)
      if (is.sigma.switching)
        return (lapply(short.thetas, EM.MSMAH.AR,
                       y = y, y.lagged = y.lagged,
                       z.dependent = z.dependent, z.independent = z.independent,
                       z.dependent.lagged = z.dependent.lagged, 
                       z.independent.lagged = z.independent.lagged,
                       maxit = maxit, epsilon = epsilon,
                       transition.probs.min = transition.probs.min,
                       transition.probs.max = transition.probs.max,
                       sigma.min = sigma.min,
                       state.conversion.mat = state.conversion.mat))
    else
      return (lapply(short.thetas, EM.MSMA.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     z.dependent.lagged = z.dependent.lagged, 
                     z.independent.lagged = z.independent.lagged,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max,
                     sigma.min = sigma.min,
                     state.conversion.mat = state.conversion.mat))
    else
      if (is.sigma.switching)
        return (lapply(short.thetas, EM.MSMH.AR,
                       y = y, y.lagged = y.lagged,
                       z.dependent = z.dependent, z.independent = z.independent,
                       z.dependent.lagged = z.dependent.lagged, 
                       z.independent.lagged = z.independent.lagged,
                       maxit = maxit, epsilon = epsilon,
                       transition.probs.min = transition.probs.min,
                       transition.probs.max = transition.probs.max,
                       sigma.min = sigma.min,
                       state.conversion.mat = state.conversion.mat))
    else
      return (lapply(short.thetas, EM.MSM.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     z.dependent.lagged = z.dependent.lagged, 
                     z.independent.lagged = z.independent.lagged,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max,
                     sigma.min = sigma.min,
                     state.conversion.mat = state.conversion.mat))
  }
  
  
  if (is.beta.switching)
    if (is.sigma.switching)
      return (lapply(short.thetas, EM.MSIAH.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max,
                     sigma.min = sigma.min))
    else
      return (lapply(short.thetas, EM.MSIA.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max,
                     sigma.min = sigma.min))
  else
    if (is.sigma.switching)
      return (lapply(short.thetas, EM.MSIH.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max,
                     sigma.min = sigma.min))
    else
      return (lapply(short.thetas, EM.MSI.AR,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max,
                     sigma.min = sigma.min))
}