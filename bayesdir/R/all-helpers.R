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

#' Compute the score function s(eta, x) for vMF.
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


#' Fisher information matrix for the vMF model
#'
#' @param C  natural parameter vector (length p)
#' @param p  dimension (defaults to length(C))
#' @return   p x p Fisher information matrix
#' @export
get_FIM_vmf <- function(C, p = length(C)) {
  kappa <- sqrt(sum(C^2))
  I_p <- diag(p)
  if (kappa < 1e-6) return(I_p / p)

  mu <- C / kappa
  Ap <- get_Ap(kappa, p)
  Ap_prime <- 1 - Ap^2 - ((p - 1) / kappa) * Ap

  mu_mu <- mu %*% t(mu)
  Ap_prime * mu_mu + (Ap / kappa) * (I_p - mu_mu)
}


#' Compute the sum of weights sum_{t = start}^{end} 1/(t + N)^2.
#' @param n number of observed data points
#' @param M number of iterations
#' @export
#' Helper to compute weights sums
compute_weights <- function(n, M) {
  t_early <- 1:M
  W_nm <- sum(1 / (n + t_early)^2)
  # Infinite sum R_nm = sum_{t=M+1}^\infty 1/(n+t)^2 approx 1/(n+M)
  R_nm <- 1 / (n + M) 
  list(W_nm = W_nm, R_nm = R_nm)
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
      kappa_hat <- stats::uniroot(function(k) get_Ap(k, 2) - Rbar, c(0.001, 100))$root
    }
  
  c_hat <- kappa_hat * c(cos(mu_hat), sin(mu_hat))
  
  return(c_hat)
}
