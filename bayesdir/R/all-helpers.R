# Generated from _main.Rmd: do not edit by hand

#' Compute the function A(kappa) = I_{p/2}(kappa) / I_{p/2-1}(kappa)
#' for the von Mises–Fisher distribution.
#' @param kappa concentration parameter (non‑negative scalar)
#' @param p space dimension (default = 2)
#' @return scalar value of A(kappa)
#' @export
get_Ap <- function(kappa, p = 2) {
  if (kappa < 1e-6) return(0)
  besselI(kappa, p/2, expon.scaled = TRUE) /
    besselI(kappa, p/2 - 1, expon.scaled = TRUE)
}

#' Compute the score function s(eta, x) = x - E[X] for vMF.
#' @param eta natural parameter vector (length p)
#' @param x unit vector observation (length p)
#' @param p space dimension (default = 2)
#' @return score vector of length p
#' @export
score_fn <- function(eta, x, p = 2) {
  kappa <- sqrt(sum(eta^2))
  if (kappa < 1e-6) return(x)
  expected <- get_Ap(kappa, p) * eta / kappa
  x - expected
}



#' Compute the sum of weights sum_{t = start}^{end} 1/(t + N)^2.
#' @param start starting index (t >= 1)
#' @param end ending index
#' @param N number of observed data points
#' @return scalar value of the sum
#' @export
compute_weight_sum <- function(start, end, N) {
  sum(1 / ((start:end) + N)^2)
}

#' Compute MLE given data X 
#' @param X nrow = p ncol = number of observation
#' @return MLE as a p-dim vector
#' @export
mle_vmf <- function(X){
  n <- ncol(X)
  
  sum_vec <- rowSums(X)  # 2-dim vector sum
  Robs <- sqrt(sum(sum_vec^2))
  Rbar <- Robs / n
  mu_hat <- atan2(sum_vec[2], sum_vec[1])
  
  if (Rbar < 0.001) {
      kappa_hat <- 0
    } else if (Rbar > 0.999) {
      kappa_hat <- 1 / (2 * (1 - Rbar))
    } else {
      kappa_hat <- uniroot(function(k) get_Ap(k, 2) - Rbar, c(0.001, 100))$root
    }
  
  c_hat <- kappa_hat * c(cos(mu_hat), sin(mu_hat))
  
  return(c_hat)
}

