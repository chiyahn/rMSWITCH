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
                              sigma0, 
                              z.dependent.lagged = NULL,
                              z.independent.lagged = NULL,
                              is.MSM = FALSE,
                              force.persistence = FALSE,
                              penalty.term = 0)
{
  n <- length(y)
  M <- ncol(as.matrix(short.thetas[[1]]$transition.probs))
  s <- 0
  if (!is.null(short.thetas[[1]]$beta))
  {
    s <- nrow(as.matrix(short.thetas[[1]]$beta)) 
  }
  if (is.null(z.dependent))
  {
    z.dependent <- as.matrix(rep(0,n))
    z.dependent.lagged <- matrix(0, ncol = max(0, s), nrow = n)
  }
  if (is.null(z.independent))
  {
    z.independent <- as.matrix(rep(0,n))
    z.independent.lagged <- matrix(0, ncol = max(0, s), nrow = n)
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
  
  if (force.persistence)
    if (is.sigma.switching)
      return (lapply(short.thetas, EM.MSIH.AR.QuZhuo,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max,
                     sigma.min = sigma.min))
    else
      return (lapply(short.thetas, EM.MSI.AR.QuZhuo,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max,
                     sigma.min = sigma.min))
  
  
  if (penalty.term > 0)
    if (is.sigma.switching)
      return (lapply(short.thetas, EM.MSIH.AR.Penalized,
                     y = y, y.lagged = y.lagged,
                     z.dependent = z.dependent, z.independent = z.independent,
                     maxit = maxit, epsilon = epsilon,
                     transition.probs.min = transition.probs.min,
                     transition.probs.max = transition.probs.max,
                     sigma.min = sigma.min, sigma0 = sigma0,
                     penalty.term = penalty.term))
  
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