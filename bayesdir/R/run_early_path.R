# Generated from _main.Rmd: do not edit by hand

#' Run the early path for a single chain (Step 1 of Hybrid Method)
#' @param c0 initial natural parameter (MLE)
#' @param n observed sample size
#' @param M predictive steps
#' @return A list containing displacement D and quadratic variation Q
#' @export
run_early_path <- function(c0, n, M) {
  p <- length(c0)
  c_curr <- c0
  Q_mat <- matrix(0, p, p)
  
  for (t in 1:M) {
    gamma_nt <- 1 / (n + t)
    # Sample from predictive distribution
    x_new <- rstiefel::rmf.vector(c_curr)
    # Compute score
    s <- score_fn(c_curr, x_new, p)
    # Accumulate quadratic variation
    Q_mat <- Q_mat + (gamma_nt^2) * (s %*% t(s))
    # Update parameter
    c_curr <- c_curr + gamma_nt * s
  }
  
  list(
    D_nm = c_curr - c0,
    Q_nm = Q_mat
  )
}

#' Run Hybrid Predictive Resampling for vMF
#' @param X data matrix (p x n)
#' @param B number of independent chains
#' @param M number of early predictive steps
#' @param n_cores number of cores for parallel processing
#' @export
run_hybrid_vmf <- function(X, B = 5000, M = 50, n_cores = 1) {
  n <- ncol(X)
  p <- nrow(X)
  
  # 0. Initial Estimate (MLE)
  c0 <- mle_vmf(X) 
  weights <- compute_weights(n, M)
  
  # 1. Run B chains to get Early Path data (Parallelized)
  early_results <- parallel::mclapply(1:B, function(b) {
    run_early_path(c0, n, M)
  }, mc.cores = n_cores)
  
  # Extract Ds and Qs
  D_list <- lapply(early_results, `[[`, "D_nm")
  Q_list <- lapply(early_results, `[[`, "Q_nm")
  
  # 2. Pooling (Section 1 & 5 of methodology)
  bar_Q_nm <- Reduce("+", Q_list) / B
  hat_I_nm <- bar_Q_nm / weights$W_nm
  
  # Inversion helper with stabilization
  safe_inv <- function(Mat) {
    eig <- eigen((Mat + t(Mat)) / 2)
    vals_inv <- 1 / pmax(eig$values, 1e-10)
    eig$vectors %*% diag(vals_inv) %*% t(eig$vectors)
  }
  hat_I_inv <- safe_inv(hat_I_nm)
  
  # 3. Tail Generation and Post-hoc Correction (Section 8)
  # Tail variance: R_nm * hat_I_nm
  tail_cov <- weights$R_nm * hat_I_nm
  # Pre-calculate matrix square root for tail sampling: Z ~ N(0, tail_cov)
  tail_eig <- eigen(tail_cov)
  tail_sqrt <- tail_eig$vectors %*% diag(sqrt(pmax(tail_eig$values, 0))) %*% t(tail_eig$vectors)
  
  c_final_samples <- sapply(1:B, function(b) {
    # Generate Gaussian Tail Z
    Z_nm <- tail_sqrt %*% rnorm(p)
    # Apply Correction: c_hyb = c0 + I_inv * (D + Z)
    c_hyb <- c0 + hat_I_inv %*% (D_list[[b]] + Z_nm)
    return(c_hyb)
  })
  
  # 4. Transform to Polar for output
  kappa_sam <- apply(c_final_samples, 2, function(v) sqrt(sum(v^2)))
  if (p == 2) {
    # Convert to angles if dimension is 2
    mu_sam <- apply(c_final_samples, 2, function(v) atan2(v[2], v[1]))
  }else{ mu_sam <- sweep(c_final_samples, 2, kappa_sam, "/")}
  
  list(
    samples = list(
      c_sam = c_final_samples,
      mu_sam = mu_sam,
      kappa_sam = kappa_sam
    ),
    info_matrix = hat_I_nm
  )
}
