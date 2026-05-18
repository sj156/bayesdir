# Generated from _main.Rmd: do not edit by hand

#' functions for simulations study
#' @param n sample size
#' @param n_sim sim times
#' @param mu_ture,kappa_true true param
#' @param B number of chains
#' @param M number of iteration for each chain
#' @param seed 
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

