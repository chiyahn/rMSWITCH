EM.MSIA.AR <- function (theta, y, y.lagged,
                        z.dependent, z.independent, 
                        maxit, epsilon, 
                        transition.probs.min, transition.probs.max,
                        sigma.min)
{
  
  EMcppARMSIA(y, y.lagged, z.dependent, z.independent,
              theta$beta, theta$mu, theta$sigma,
              theta$gamma.dependent, theta$gamma.independent,
              theta$transition.probs, theta$initial.dist,
              maxit, epsilon,
              transition.probs.min, transition.probs.max,
              sigma.min)
}

EM.MSIAH.AR <- function (theta, y, y.lagged,
                         z.dependent, z.independent, 
                         maxit, epsilon, 
                         transition.probs.min, transition.probs.max,
                         sigma.min)
{
  
  EMcppARMSIAH(y, y.lagged, z.dependent, z.independent,
               theta$beta, theta$mu, theta$sigma,
               theta$gamma.dependent, theta$gamma.independent,
               theta$transition.probs, theta$initial.dist,
               maxit, epsilon,
               transition.probs.min, transition.probs.max,
               sigma.min)
}

EM.MSI.AR <- function (theta, y, y.lagged,
                       z.dependent, z.independent, 
                       maxit, epsilon, 
                       transition.probs.min, transition.probs.max,
                       sigma.min)
{
  
  EMcppARMSI(y, y.lagged, z.dependent, z.independent,
             theta$beta, theta$mu, theta$sigma,
             theta$gamma.dependent, theta$gamma.independent,
             theta$transition.probs, theta$initial.dist,
             maxit, epsilon,
             transition.probs.min, transition.probs.max,
             sigma.min)
}


EM.MSIH.AR <- function (theta, y, y.lagged,
                        z.dependent, z.independent, 
                        maxit, epsilon, 
                        transition.probs.min, transition.probs.max,
                        sigma.min)
{
  
  EMcppARMSIH(y, y.lagged, z.dependent, z.independent,
              theta$beta, theta$mu, theta$sigma,
              theta$gamma.dependent, theta$gamma.independent,
              theta$transition.probs, theta$initial.dist,
              maxit, epsilon,
              transition.probs.min, transition.probs.max,
              sigma.min)
}

EM.MSIH.AR.Penalized <- function (theta, y, y.lagged,
                                  z.dependent, z.independent, 
                                  maxit, epsilon, 
                                  transition.probs.min, transition.probs.max,
                                  sigma.min, sigma0, penalty.term)
{
  EMcppARMSIHPenalized(y, y.lagged, z.dependent, z.independent,
              theta$beta, theta$mu, theta$sigma,
              theta$gamma.dependent, theta$gamma.independent,
              theta$transition.probs, theta$initial.dist,
              maxit, epsilon,
              transition.probs.min, transition.probs.max,
              sigma.min, sigma0, penalty.term)
}

EM.MSI.AR.QuZhuo <- function (theta, y, y.lagged,
                       z.dependent, z.independent, 
                       maxit, epsilon, 
                       transition.probs.min, transition.probs.max,
                       sigma.min)
{
  
  EMcppARMSIQuZhuo(y, y.lagged, z.dependent, z.independent,
             theta$beta, theta$mu, theta$sigma,
             theta$gamma.dependent, theta$gamma.independent,
             theta$transition.probs, theta$initial.dist,
             maxit, epsilon,
             transition.probs.min, transition.probs.max,
             sigma.min)
}

EM.MSIH.AR.QuZhuo <- function (theta, y, y.lagged,
                        z.dependent, z.independent, 
                        maxit, epsilon, 
                        transition.probs.min, transition.probs.max,
                        sigma.min)
{
  
  EMcppARMSIHQuZhuo(y, y.lagged, z.dependent, z.independent,
              theta$beta, theta$mu, theta$sigma,
              theta$gamma.dependent, theta$gamma.independent,
              theta$transition.probs, theta$initial.dist,
              maxit, epsilon,
              transition.probs.min, transition.probs.max,
              sigma.min)
}

EM.MSMAH.AR <- function (theta, y, y.lagged,
                         z.dependent, z.independent, 
                         z.dependent.lagged, z.independent.lagged,
                         maxit, epsilon, 
                         transition.probs.min, transition.probs.max,
                         sigma.min, state.conversion.mat)
{
  
  EMcppARMSMAH(y, y.lagged, z.dependent, z.independent,
               z.dependent.lagged, z.independent.lagged,
               theta$beta, theta$mu, theta$sigma,
               theta$gamma.dependent, theta$gamma.independent,
               theta$transition.probs, theta$initial.dist,
               maxit, epsilon,
               transition.probs.min, transition.probs.max,
               sigma.min, state.conversion.mat)
}

EM.MSMA.AR <- function (theta, y, y.lagged,
                        z.dependent, z.independent, 
                        z.dependent.lagged, z.independent.lagged,
                        maxit, epsilon, 
                        transition.probs.min, transition.probs.max,
                        sigma.min, state.conversion.mat)
{
  EMcppARMSMA(y, y.lagged, z.dependent, z.independent,
               z.dependent.lagged, z.independent.lagged,
               theta$beta, theta$mu, theta$sigma,
               theta$gamma.dependent, theta$gamma.independent,
               theta$transition.probs, theta$initial.dist,
               maxit, epsilon,
               transition.probs.min, transition.probs.max,
               sigma.min, state.conversion.mat)
}

EM.MSM.AR <- function (theta, y, y.lagged,
                       z.dependent, z.independent, 
                       z.dependent.lagged, z.independent.lagged,
                       maxit, epsilon, 
                       transition.probs.min, transition.probs.max,
                       sigma.min, state.conversion.mat)
{
  
  EMcppARMSM(y, y.lagged, z.dependent, z.independent,
             z.dependent.lagged, z.independent.lagged,
             theta$beta, theta$mu, theta$sigma,
             theta$gamma.dependent, theta$gamma.independent,
             theta$transition.probs, theta$initial.dist,
             maxit, epsilon,
             transition.probs.min, transition.probs.max,
             sigma.min, state.conversion.mat)
}


EM.MSMH.AR <- function (theta, y, y.lagged,
                        z.dependent, z.independent, 
                        z.dependent.lagged, z.independent.lagged,
                        maxit, epsilon, 
                        transition.probs.min, transition.probs.max,
                        sigma.min, state.conversion.mat)
{
  EMcppARMSMH(y, y.lagged, z.dependent, z.independent,
              z.dependent.lagged, z.independent.lagged,
              theta$beta, theta$mu, theta$sigma,
              theta$gamma.dependent, theta$gamma.independent,
              theta$transition.probs, theta$initial.dist,
              maxit, epsilon,
              transition.probs.min, transition.probs.max,
              sigma.min, state.conversion.mat)
}

