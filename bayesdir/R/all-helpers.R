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


#' functions for simulations study
#' @param n sample size
#' @param n_sim sim times
#' @param mu_ture,kappa_true true param
#' @param B number of chains
#' @param M number of iteration for each chain
#' @param seed for reproducability
#' @param n_cores number of cores you can use
#' @param save_all save data or not
#' @return all the data with a summary we gonna used
#' @export
sim_martingale_study <- function(n, n_sim = 500,
                                 mu_true = 2, kappa_true = 4,
                                 B = 5000,M = 50,
                                 seed = 1234,
                                 n_cores = 1,
                                 save_all = FALSE) {

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  # ---- Storage ----
  cover_mu <- cover_kappa <- numeric(n_sim)
  len_mu <- len_kappa <- numeric(n_sim)
  mean_mu_vec <- mean_kappa_vec <- numeric(n_sim)
  times <- numeric(n_sim)
  if (save_all) {results <- vector("list", n_sim)}
  
  eta_true <- kappa_true * c(cos(mu_true), sin(mu_true))
  
  
  # ---- Simulation loop ----
  for (i in 1:n_sim) {
    # Generate observations
    X <- replicate(n,rstiefel::rmf.vector(eta_true))
    
    
    # Run B independent chains
    t0 <- proc.time()[3]
    chains<- run_predresamp_vmf(X, B = 5000, M=50,n_cores = n_cores,
                                return_chains = save_all)
    t1 <- proc.time()[3]
    times[i] <- t1 - t0
    
    # Extract samples
    sample_list <- chains$sample
    
    mu_sam <- sample_list$mu_sam
    kappa_sam <- sample_list$kappa_sam
    
    mean_mu_vec[i] <- mean(mu_sam)
    mean_kappa_vec[i] <- mean(kappa_sam)
    
    # 95% credible intervals for mu and kappa
    ci_mu <- quantile(mu_sam, probs = c(0.025, 0.975))
    cover_mu[i] <- (mu_true >= ci_mu[1] && mu_true <= ci_mu[2])
    len_mu[i] <- ci_mu[2] - ci_mu[1]
    
    ci_kappa <- quantile(kappa_sam, c(0.025, 0.975))
    cover_kappa[i] <- (kappa_true >= ci_kappa[1]) & (kappa_true <= ci_kappa[2])
    len_kappa[i] <- ci_kappa[2] - ci_kappa[1]
    
    if (save_all) {results[[i]]<- chains}
    
  }
  
  # ---- Build output table ----
  summarys <- data.frame(
    n = n,
    Coverage_mu = mean(cover_mu) * 100,
    Length_mu = mean(len_mu),
    Coverage_kappa = mean(cover_kappa) * 100,
    Length_kappa = mean(len_kappa),
    Mean_SE_mu      = sd(mean_mu_vec),
    Mean_SE_kappa   = sd(mean_kappa_vec),
    Mean_Time = mean(times)
  )
  
  if (save_all) {
    list(summary = summarys, data_all = results)
  } else {
    summarys
  }
}


