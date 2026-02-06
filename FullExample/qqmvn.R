library(mvtnorm)
## QQ plot function to assess if data are multivariate normal
## idea: compute Mahalanobis distances, md, using *sample* mean and *sample* cov.  
## If data have sample size n and dimension p
## Then md * n / (n-1)^2 are Beta (p/2, (n -p -1)/2) 
## Contains code to generate CI intervals around the qqline 
## But this is a bit intensive to run as involves Monte Carlo steps
qqmvn <- function(y, CI = TRUE, ci_lty = 2, ci_lwd = 2, ci_col = 2, 
                  line_col = 1, line_lwd = 1, line_lty =1, 
                  xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", 
                  ...){
  n <- nrow(y)
  p <- ncol(y)
  ybar <- colMeans(y)
  
  S <- cov(y)
  ## normalising constant of MVN dist.
  c1 <- p/2*log(2*pi) + 0.5*log(det(S)) 
  ## cheap/lazy way to compute Mahalanobis distance
  md1 <- apply(y, 1, mvtnorm::dmvnorm, mean = ybar, sigma = S, log = TRUE)
  ## gets rid of normalising constant term
  md1 <- -2*(md1 + c1) 
  
  qp <- qbeta(ppoints(n), p/2, (n -p -1)/2)
  sm <- md1*n/(n-1)^2
  plot(qp, sort(sm), xlab = xlab, ylab = ylab, ...)
  abline(0, 1, col = line_col, lwd = line_lwd, lty = line_lty)
  
  if(CI){
    rbmat <- replicate(1000, sort(rbeta(n, p/2, (n -p -1)/2)))
    rbci <- apply(rbmat, 1, quantile, probs = c(0.025, 0.975))
    lines(qbeta(ppoints(n), p/2, (n -p -1)/2), rbci[1, ], col = ci_col, lwd = ci_lwd, lty = ci_lty)
    lines(qbeta(ppoints(n), p/2, (n -p -1)/2), rbci[2, ], col = ci_col, lwd = ci_lwd, lty = ci_lty)
  }
}

## simulated examples, ignore for now ----
# y1 <- rmvnorm(200, c(0, 0), sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2))
# y2 <- rmvt(200, sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2), df = 3)
# 
# qqmvn(y1, CI = TRUE, main = "Simulated MVN")
# qqmvn(y2, CI = TRUE, main = "Simulated MVT")
# qqmvn(exp(y1), CI = TRUE, main = "Simulated MV log-Normal")

## function that makes our IS sampler outputs compatible with qqmvn ----
## one way this is inefficient: 
## it generates a new set of samples every time it is called 
## should probably construct a separate function to generate samples
## This can then be called by the plotting function 
check_mvn <- function(alpha, 
                      iter, 
                      coeff,
                      cov,
                      dist,
                      ex_info,
                      CI=FALSE, 
                      main = "", ...){
  ## Calls IS function, for (presumably optimal) alpha and iter settings
  get_IS_params_ws <- IS(alpha = alpha, 
                         iter = iter, 
                         coeff = coeff,
                         cov = cov,
                         dist=dist,
                         ex_info = ex_info)
  
  ## Generate set of parameter samples
  ## This uses sampling with replacement, weighted by importance weights 
  samples_weighted <- sample(1:iter, 
                             replace = TRUE,
                             size = 5000,
                             prob = get_IS_params_ws$weights)
  
  ## generate qq plots -- straightforward if using exponential as 1-dim
  if(length(coeff) == 1){
    qqnorm(get_IS_params_ws$samples[samples_weighted, ], main = main, ...)
    qqline(get_IS_params_ws$samples[samples_weighted, ], col = 2, lty = 2, lwd = 2)
  } else{
    ## otherwise use qqmvn function
    qqmvn(get_IS_params_ws$samples[samples_weighted, ], CI = CI, main = main, ...)
    ## another option -- but can't easily manipulate titles...
    ## install.packages("MVN")
    ## library("MVN")
    ## MVN::multivariate_diagnostic_plot(get_IS_params_ws$samples[samples_weighted, ])
  }
}
## wrapper function to facilitate, e.g., using lapply with is.models
check_mvn_wrapper <- function(x, ...) {
  check_mvn(alpha = x$alpha_star, 
            iter = x$iter_star, 
            coeff=x$orig$coefficients,
            cov=x$orig$cov,
            dist=x$orig$dist,
            ex_info = x$ex_info,
            CI=TRUE,
            main = paste("Q-Q Plot for", x$orig$dist, "distribution", sep = " "),
            ...
  )
}
