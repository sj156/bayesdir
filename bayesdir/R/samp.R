# Generated from _main.Rmd: do not edit by hand

#' Run a single martingale posterior chain (no preconditioner).
#' @param c0 initial natural parameter vector(usually use MLE)
#' @param N number of observed data points
#' @param M number of iteration (default = 20)
#' @return A list with components:
#'   \item{c_final}{final natural parameter vector (after possible correction)}
#'   \item{mu_final}{mean direction (angle in radians if p=2)}
#'   \item{kappa_final}{concentration parameter}
#'   \item{x_chain}{all new x generated}
#'   \item{c_chain}{all new c generated}
#'@export
run_single_martingale <- function(c0, N, M = 20,save_path = FALSE) {
                                  
  p <- length(c0)
  c_curr <- c0
  V_mat <- matrix(0, p, p)          # accumulated un‑preconditioned score outer products
  

  if (save_path) {
    x_chain <- matrix(0, p, M)
    c_chain <- matrix(0, p, M)
  }
  
  w_cur <- compute_weight_sum(1, M, N)       
  w_rem <- compute_weight_sum(M + 1, M + 10000, N)  
  total_w <- w_cur + w_rem
  
 # ---- Hybrid acceleration: M steps + Gaussian correction ----
  w_cur <- compute_weight_sum(1, M, N)
  for (t in 1:M) {
      x_new <- rstiefel::rmf.vector(c_curr)
      s <- score_fn(c_curr, x_new, p)
      wt <- 1 / (t + N)
      V_mat <- V_mat + wt^2 * tcrossprod(s)
      c_curr <- c_curr + wt * s
      
      if (save_path) {
      x_chain[, t] <- x_new
      c_chain[, t] <- c_curr
    }
      
  }
  

  r <- w_rem / w_cur
  Cov_rem <- r * V_mat
  Sigma_used <- V_mat + Cov_rem   # total covariance estimate
  
  # Gaussian correction
  eig <- eigen(Cov_rem, symmetric = TRUE)
  eig$values[eig$values < 0] <- 0
  sqrt_Cov_rem <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  epsilon <- rnorm(p)
  c_curr <- c_curr + sqrt_Cov_rem %*% epsilon
    
    
    
  
  # ---- Post‑hoc variance correction ----
  delta <- c_curr - c0
  
  # Use eigendecomposition for stable inversion of the covariance matrix
  eigS <- eigen(Sigma_used, symmetric = TRUE)
  eigS$values[eigS$values < 1e-12] <- 0
  inv_Sigma <- eigS$vectors %*% diag(ifelse(eigS$values > 0, 1/eigS$values, 0)) %*% t(eigS$vectors)
  c_curr <- c0 + total_w * (inv_Sigma %*% delta)

  # ---- Convert to polar coordinates (mu, kappa) for p=2----
  kappa_final <- sqrt(sum(c_curr^2))
  if (p == 2) {
    mu_final <- atan2(c_curr[2], c_curr[1])
    if (mu_final < 0) mu_final <- mu_final + 2 * pi
  } else {
    mu_final <-(1/kappa_final) * c_curr
  }
  
  res <- list(
    c_final     = c_curr,
    mu_final    = mu_final,
    kappa_final = kappa_final
  )
  if (save_path) {
    res$x_chain <- x_chain
    res$c_chain <- c_chain
  }
  
  res
}

#' run B chains for empirical density
#' @param X data dim = c(p,n)
#' @param B expected sample size
#' @param M number of iteration
#' @param n_cores number of cores used in mclapply
#' @param save_chain whether to save the entire chain.
#' @export
run_predresamp_vmf <- function(X,B = 5000, M = 50, n_cores = 1,return_chains = FALSE){
  
  n  <- ncol(X)
  p  <- nrow(X)
  c0 <- mle_vmf(X)
  
  chains <- mclapply(1:B, function(b) run_single_martingale(c0 = c0, N = n, M = M,save_path = return_chains),
                     mc.cores = n_cores)
  
  
  c_sam <- sapply(chains, `[[`, "c_final")
  kappa_sam <-sapply(chains,`[[`, "kappa_final")
  mu_sam <-sapply(chains,`[[`, "mu_final")
  
  sample_list <- list(
    c_sam     = c_sam,
    mu_sam    = mu_sam,
    kappa_sam = kappa_sam
  )

  out <- list(sample = sample_list)
  if (return_chains) {
    out$chain <- chains
  }
  
  out
}
