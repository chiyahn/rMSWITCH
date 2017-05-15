TestMSARSequence <- function(y, z.dependent = NULL, z.independent = NULL,
                             s = 1, max.M = 3,
                             is.beta.switching = FALSE,
                             is.sigma.switching = TRUE,
                             is.MSM = FALSE,
                             cl = NULL, parallel = TRUE,
                             short.n = 5,
                             transition.probs.min = 0.01,
                             sigma.min = 0.02,
                             nloptr = FALSE,
                             bootstrap.count = 199)
{
  initial.theta <- NULL
  aic <- bic <- double(max.M)
  pvals <- LRT.statistics <- double(max.M)
  log.likelihoods <- double(max.M)
  model.chosen <- NULL
  
  for (M in 1:max.M)
  {
    msar.model0 <- EstimateMSAR(y = y, 
                                z.dependent = z.dependent, z.independent = z.independent,
                                M = M, s = s,
                                is.beta.switching = is.beta.switching,
                                is.sigma.switching = is.sigma.switching,
                                is.MSM = is.MSM,
                                initial.theta = initial.theta,
                                short.n = short.n,
                                transition.probs.min = transition.probs.min,
                                sigma.min = sigma.min,
                                nloptr = nloptr,
                                crit.method = NULL)
    
    msar.model1 <- EstimateMSAR(y = y, 
                                z.dependent = z.dependent, z.independent = z.independent,
                                M = (M+1), s = s,
                                is.beta.switching = is.beta.switching,
                                is.sigma.switching = is.sigma.switching,
                                is.MSM = is.MSM,
                                short.n = short.n,
                                transition.probs.min = transition.probs.min,
                                sigma.min = sigma.min,
                                nloptr = nloptr,
                                crit.method = NULL)
    
    aic[M] <- -2*x$log.likelihood + 2*NumberOfParameters(msar.model0$theta)
    bic[M] <- -2*x$log.likelihood + log(n)*NumberOfParameters(msar.model0$theta)
    
    LRT.statistics[M] <- 2*(msar.model1$log.likelihood - msar.model0$log.likelihood)
    initial.theta <- msar.model1$theta
    print(sprintf("%d-regime switching model estimate:\n", M))
    print(msar.model0)
    
    cat("LRT statistic ", sprintf('%.3f', LRT.statistics[M]), "\n")
    
    crit.result <- TestMSARCritBoot(LRT.statistic0 = LRT.statistics[M],
                                    msar.model0 = msar.model0,
                                    s = s,
                                    is.beta.switching,
                                    is.sigma.switching,
                                    y = y,
                                    z.dependent = z.dependent,
                                    z.independent = z.independent,
                                    cl = cl, parallel = parallel,
                                    bootstrap.count = bootstrap.count,
                                    short.n = short.n,
                                    transition.probs.min = transition.probs.min,
                                    sigma.min = sigma.min,
                                    nloptr = nloptr)
    pvals[M] <- crit.result$pval
  }
  
  for (M in 1:max.M)
    if (pvals[M] >= 0.05) {
      cat(sprintf("\nThe number of components selected by Sequential Hypothesis Testing (alpha=0.05) = %.i", M), 
          " \n")
      cat(sprintf("The number of components selected by AIC = %.i", 
                  which.min(aics)), " \n")
      cat(sprintf("The number of components selected by BIC = %.i", 
                  which.min(bics)), " \n")
      
      model.chosen   <-  EstimateMSAR(y = y, 
                                      z.dependent = z.dependent, z.independent = z.independent,
                                      M = M, s = s,
                                      is.beta.switching = is.beta.switching,
                                      is.sigma.switching = is.sigma.switching,
                                      is.MSM = is.MSM,
                                      crit.method = NULL)
      
      cat(sprintf("\nThe summary of estimated Markov %.i", M), "regime switching model is as follows \n")
      print(summary(model.chosen))
      break
    }
  
  return (model.chosen)
}

#' Returns a LRT statistic given the data for y, z, s based on
#' the null hypothesis H_0: M = m_0 with the alternative H_1: M = m_0 + 1
#' @export
#' @title TestMSAR
#' @name TestMSAR
#' @param y n by 1 vector of data
#' @param z.dependent n by p.dep matrix of data for
#' switching exogenous variables
#' @param z.independent n by p.indep matrix of data for
#' non-switching exogenous variables
#' @param M The number of states in the model in a null hypothesis;
#' the alternative will be the model with (M+1) regimes.
#' @param s The number of terms used for AR(s)
#' @param is.beta.switching Specifies whether autoregressive terms are switching
#' @param is.sigma.swithcing Specifies whether the model is heteroscedastic.
#' @param is.MSM Specifies whether the model is switching in mean (MSM) or
#' intercept (MSI).
#' @param cl Cluster used for parallelization; if it is \code{NULL},
#' the system will automatically create a new one for computation accordingly.
#' @param parallel Determines whether package \code{doParallel} is used for
#' calculation in parallelized process.
#' @param crit.method Method used to compute p-values,
#' one of \code{"none"} and \code{"boot"}. The default option is
#' '\code{"none"}. When \code{method = "boot"},
#' the p-values are generated by bootstrapping.
#' @param short.n Number of initial draws for EM estimation
#' @param transition.probs.min Minimum set for transition prob. matrix
#' @param sigma.min Minimum set for variance.
#' @param nloptr Determines whether nonlinear optimization package is used.
#' @param bootstrap.count The number of bootstrap samples;
#' by default, it is set to be 199
#' @param msar.model0 Estimated model for null hypothesis
#' @param msar.model1 Estimated model for alternative hypothesis
#' @return A list of class \code{msar.model} with items:
#' \item{LRT.statistic}{An LRT statistic computed from the null M = M0 and
#' alternative M = (M+1)}
#' \item{crit}{The critical value to reject the null hypothesis}
#' \item{pvals}{p-value at 10%, 5%, and 1%}
#' \item{msar.model0}{A list of class \code{msar.model} estimated
#' from a null hypothesis M = M0}
#' \item{msar.model1}{A list of class \code{msar.model} estimated
#' from an alternative hypothesis M = M0+1}
#' @examples
TestMSAR <- function(y, z.dependent = NULL, z.independent = NULL,
                      M = 1, s = 1,
                      is.beta.switching = FALSE,
                      is.sigma.switching = TRUE,
                      is.MSM = FALSE,
                      cl = NULL, parallel = TRUE,
                      crit.method = c('none', 'bootstrap'),
                      estimate.fisher = TRUE,
                      short.n = 5,
                      transition.probs.min = 0.01,
                      sigma.min = 0.02,
                      nloptr = FALSE,
                      bootstrap.count = 199,
                      msar.model0 = NULL, msar.model1 = NULL)
{
  crit.method <- match.arg(crit.method)

  if (is.null(msar.model0))
    msar.model0 <- EstimateMSAR(y = y,
                                z.dependent = z.dependent,
                                z.independent = z.independent,
                                M = M, s = s,
                                is.beta.switching = is.beta.switching,
                                is.sigma.switching = is.sigma.switching,
                                is.MSM = is.MSM,
                                estimate.fisher = estimate.fisher,
                                short.n = short.n,
                                transition.probs.min = transition.probs.min,
                                sigma.min = sigma.min,
                                nloptr = nloptr)
  if (is.null(msar.model1))
    msar.model1 <- EstimateMSAR(y = y,
                                z.dependent = z.dependent,
                                z.independent = z.independent,
                                M = (M + 1), s = s,
                                is.beta.switching = is.beta.switching,
                                is.sigma.switching = is.sigma.switching,
                                is.MSM = is.MSM,
                                estimate.fisher = estimate.fisher,
                                short.n = short.n,
                                transition.probs.min = transition.probs.min,
                                sigma.min = sigma.min,
                                nloptr = nloptr)

  LRT.statistic <- 2*(msar.model1$log.likelihood - msar.model0$log.likelihood)

  if (crit.method == "bootstrap") {
    crit.result <- TestMSARCritBoot(LRT.statistic0 = LRT.statistic,
                                    msar.model0 = msar.model0, 
                                    s = s,
                                    is.beta.switching,
                                    is.sigma.switching,
                                    y = y,
                                    z.dependent = z.dependent,
                                    z.independent = z.independent,
                                    cl = cl, parallel = parallel,
                                    bootstrap.count = bootstrap.count,
                                    short.n = short.n,
                                    transition.probs.min = transition.probs.min,
                                    sigma.min = sigma.min,
                                    nloptr = nloptr)
  }

  else {
    crit.result <- list()
    crit.result$crit <- rep(NA, 3)
    crit.result$pval <- NA
    crit.result$bootstrap.statistics <- NULL
    crit.result$bootstrap.estimates.null <- NULL
  }
  msar.test <- list(LRT.statistic = LRT.statistic,
                    crit = crit.result$crit,
                    pval = crit.result$pval,
                    msar.model0 = msar.model0, msar.model1 = msar.model1,
                    bootstrap.statistics = crit.result$bootstrap.statistics,
                    bootstrap.estimates.null = crit.result$bootstrap.estimates.null,
                    crit.method = crit.method)
  class(msar.test) <- "msar.test"
  return (msar.test)
}

#' Calculate p-values and critical value given a LRT.statistic based on
#' a bootstrapping method.
TestMSARCritBoot <- function (LRT.statistic0, 
                              msar.model0, s,
                              is.beta.switching,
                              is.sigma.switching,
                              y, z.dependent, z.independent,
                              cl = NULL, parallel = TRUE,
                              bootstrap.count = 199,
                              short.n = 5,
                              transition.probs.min = 0.01,
                              sigma.min = 0.02,
                              nloptr = FALSE)
{
  theta0 <- msar.model0$theta
  M <- nrow(theta0$transition.probs)
  n <- nrow(msar.model0$posterior.probs.smoothed)
  bootstrap.samples <- GenerateSamples(theta = theta0, n = n,
                                      replications = bootstrap.count,
                                      initial.y.set = y[1:s],
                                      is.MSM = msar.model0$is.MSM)

  test.results <- list()
  if (parallel) {
    if (is.null(cl))
      cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    test.results <- foreach (col.index = 1:bootstrap.count) %dopar% {
        library(normalregMix)
        library(nloptr)
        test.result <-  TestMSAR (y = bootstrap.samples[,col.index],
                  z.dependent = z.dependent, z.independent = z.independent,
                  M = M, s = s,
                  is.beta.switching = is.beta.switching,
                  is.sigma.switching = is.sigma.switching,
                  is.MSM = msar.model0$is.MSM,
                  cl = NULL, parallel = FALSE,
                  crit.method = 'none', estimate.fisher = FALSE,
                  short.n = short.n,
                  transition.probs.min = transition.probs.min,
                  sigma.min = sigma.min,
                  nloptr = nloptr)
        return (list(LRT.statistic = test.result$LRT.statistic,
                     bootstrap.estimate.null = 
                       ThetaToReducedColumn((test.result$msar.model0)$theta)))
      }
    on.exit(cl)
  }
  else
    test.results <- apply(bootstrap.samples, 2, function(col)
    {
      test.result <-TestMSAR (y = col, 
                              z.dependent = z.dependent, z.independent = z.independent,
                              M = M, s = s,
                              is.beta.switching = is.beta.switching,
                              is.sigma.switching = is.sigma.switching,
                              is.MSM = msar.model0$is.MSM,
                              cl = NULL, parallel = FALSE,
                              crit.method = 'none', estimate.fisher = FALSE,
                              short.n = short.n,
                              transition.probs.min = transition.probs.min,
                              sigma.min = sigma.min,
                              nloptr = nloptr)
      return (list(LRT.statistic = test.result$LRT.statistic,
                   bootstrap.estimate.null = 
                     ThetaToReducedColumn((test.result$msar.model0)$theta)))  
    })

  LRT.statistics <- sort(sapply(test.results, "[[", "LRT.statistic")) # sort it
  qs <- ceiling(bootstrap.count * c(0.90,0.95,0.99))
  crit <- LRT.statistics[qs]
  pval <- mean(LRT.statistics > LRT.statistic0)
  
  bootstrap.statistics <- sapply(test.results, "[[", "LRT.statistic")
  bootstrap.estimates.null <- t(sapply(test.results, "[[", "bootstrap.estimate.null"))

  return (list(crit = crit, pval = pval, 
               bootstrap.estimates.null = bootstrap.estimates.null,
               bootstrap.statistics = bootstrap.statistics))
}
