print.msar.model <- function(x, ...) {
  n <- nrow(x$posterior.probs.filtered)
  theta <- x$theta
  

  coeffs <- ThetaToReducedColumn(theta)
  if (x$is.MSM)
    coeffs <- ThetaToEssentialColumn(theta)
  
  if (!is.null(x$fisher.estimated))
  {
    ses <- sqrt(diag(solve(x$fisher.estimated)))/sqrt(n)
    tvals <- coeffs / ses
    
    coeff.matrix <- data.frame(cbind(
                          estimate = coeffs, 
                          se = ses, 
                          t.value = tvals, 
                          lb.95 = coeffs - ses * qnorm(0.975) / sqrt(n),
                          ub.95 = coeffs + ses * qnorm(0.975) / sqrt(n),
                          p.value = 2*pnorm(-abs(tvals))))
    rownames(coeff.matrix) <- GetColumnNames(theta, essential = x$is.MSM)
    printCoefmat(coeff.matrix, P.values = TRUE, has.Pvalue = TRUE)
  }
  else
  {
    coeff.matrix <- data.frame(cbind(
      estimate = coeffs))
    rownames(coeff.matrix) <- GetColumnNames(theta, essential = x$is.MSM)
    print(coeff.matrix)
  }
  aic <- -2*x$log.likelihood + 2*NumberOfParameters(theta)
  bic <- -2*x$log.likelihood + log(n)*NumberOfParameters(theta)
  cat(sprintf("\nlog.likelihood at estimate:  %.3f\n", x$log.likelihood))
  cat(sprintf("AIC: %.3f\n", aic))
  cat(sprintf("BIC: %.3f\n", bic))
  
}

print.summary.msar.model <- function(x, ...)
{
  print.msar.model(x)
}

print.msar.test <- function(x, ...) {
  M0 <- ncol(((x$msar.model0)$theta)$transition.probs)
  
  cat(sprintf("1. Estimated model under the null hypothesis of %d regime(s):\n", M0))
  print(x$msar.model0)
  cat(sprintf("\n2. Testing the null hypothesis of %d regime(s):\n", M0))
  cat(c("LRT test statistic: ",sprintf('%.3f ', x$LRT.statistic)),"\n")
  if (x$crit.method == "asy") {
    cat(c("asymptotic p-value: ",sprintf('%.3f ', x$pval)),"\n")
  } else if (x$crit.method == "bootstrap") {
    cat(c("bootstrap p-value: ",sprintf('%.3f ', x$pval)),"\n")
  }
}

print.summary.msar.test <- function(x, ...)
{
  print.msar.test(x)
}