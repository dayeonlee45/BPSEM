rm(list=ls())
setwd("C:/Users/daylee/Desktop/SplineSEM/Simulation")

library(MASS)
library(rstan)
library(lavaan)
library(foreach)
library(doParallel)

# ==============================================================================
# Simulation Parameters
# ==============================================================================
N <- 100
n_reps <- 100

# [Data Gen] 참값 설정 (데이터 생성용일 뿐, Stan에는 주지 않습니다)
lambda_gen_vec <- c(1.0, 0.5, 0.5, 0.5)
beta_P2  <- 0.3; tau_P2   <- 0.0; delta_P2 <- 0.5
beta_13  <- 0.4; tau_13   <- -0.67; delta_13 <- 0.1
beta_23  <- 0.4; tau_23   <- 0.0; delta_23 <- 0.0
psi_fixed_sd  <- sqrt(0.2) # 생성용 참값
k_softplus <- 20

output_dir <- "sim_1"
if(!dir.exists(output_dir)) dir.create(output_dir)
final_csv_path <- file.path(output_dir, "sim_1.csv")
if(file.exists(final_csv_path)) file.remove(final_csv_path)

# ==============================================================================
# Helper Function
# ==============================================================================
smooth_hinge <- function(x, knot, k) {
  log1p(exp(k * (x - knot))) / k
}

# ==============================================================================
# [Stan Model] - FIXED VERSION
# ==============================================================================
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_code <- "
functions {
  vector smooth_hinge_vec(vector x, real knot, real k) {
    return log1p_exp(k * (x - knot)) / k;
  }
}
data {
  int<lower=1> N;            
  int<lower=1> P;            
  int<lower=1> M;            
  matrix[N, P] Y;            
  real<lower=0> sharpness_k;
  
  real<lower=0> slab_scale; 
  real<lower=0> slab_df;
  
  vector[9] lambda_prior_mu;
  vector[9] lambda_prior_sigma;
}
parameters {
  vector[9] lambda_free;
  vector<lower=0>[P] theta_sd;
  
  vector[M] alpha;
  vector[N] p_vec;
  real beta_p2; real tau_p2;  
  real delta_raw_p2; real<lower=0> lambda_p2;
  
  real beta13; real tau13; 
  real beta23; real tau23;  
  real delta_raw_13; real<lower=0> lambda_13; 
  real delta_raw_23; real<lower=0> lambda_23; 
  
  real<lower=0> tau_global; 
  real<lower=0> caux; 
  
  vector<lower=0>[2] psi_e_sd;
  vector<lower=0>[1] psi_endo_sd;
  
  matrix[N, M] eta; 
}
transformed parameters {
  matrix[P, M] Lambda;
  real c2; real tau2;
  real lt_p2; real lt_13; real lt_23;
  real delta_p2; real delta13; real delta23; 
  
  Lambda = rep_matrix(0.0, P, M);
  {
    int k = 0; 
    for (m in 1:M) {
      int start_row = (m - 1) * 4 + 1;
      Lambda[start_row, m] = 1.0; 
      for (i in 1:3) {
        k += 1;
        Lambda[start_row + i, m] = lambda_free[k];
      }
    }
  }
  
  c2 = square(slab_scale) * caux; 
  tau2 = square(tau_global);
  lt_p2 = sqrt( c2 * square(lambda_p2) / (c2 + tau2 * square(lambda_p2)) );
  lt_13 = sqrt( c2 * square(lambda_13) / (c2 + tau2 * square(lambda_13)) );
  lt_23 = sqrt( c2 * square(lambda_23) / (c2 + tau2 * square(lambda_23)) );
  
  delta_p2 = tau_global * lt_p2 * delta_raw_p2;
  delta13  = tau_global * lt_13 * delta_raw_13;
  delta23  = tau_global * lt_23 * delta_raw_23;
}
model {
  // --- 1. Priors ---
  p_vec ~ std_normal(); 
  lambda_free ~ normal(lambda_prior_mu, lambda_prior_sigma);
  theta_sd ~ cauchy(0, 0.5); 
  alpha ~ normal(0, 1);
  
  // Structural slope parameters
  beta_p2 ~ normal(0, 0.5);
  beta13 ~ normal(0, 0.5);
  beta23 ~ normal(0, 0.5);
  
  // [FIX 1] Knot locations: Wider prior to allow better estimation
  // Changed from normal(0, 0.5) to normal(0, 1.5) to allow knots in wider range
  tau_p2  ~ normal(0, 1.5); 
  tau13   ~ normal(0, 1.5); 
  tau23   ~ normal(0, 1.5); 
  
  // Horseshoe
  tau_global ~ cauchy(0, 1.0); 
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df); 
  
  delta_raw_p2 ~ std_normal(); lambda_p2 ~ cauchy(0, 1);
  delta_raw_13 ~ std_normal(); lambda_13 ~ cauchy(0, 1);
  delta_raw_23 ~ std_normal(); lambda_23 ~ cauchy(0, 1); 
  
  // [FIX 2] Residual variances: Slightly informative prior for better identification
  // Using half-normal with scale 0.5, but could be more informative
  psi_e_sd ~ normal(0, 0.5);
  psi_endo_sd ~ normal(0, 0.5);
  
  // --- 2. Likelihood ---
  {
    vector[N] eta1_vec = eta[, 1];
    vector[N] eta2_vec = eta[, 2];
    vector[N] eta3_vec = eta[, 3];
    
    // Anchor: P -> eta1
    target += normal_lpdf(eta1_vec | alpha[1] + p_vec, psi_e_sd[1]);
    
    // Free Path: P -> eta2
    vector[N] h_p2 = smooth_hinge_vec(p_vec, tau_p2, sharpness_k); 
    target += normal_lpdf(eta2_vec | alpha[2] + beta_p2 * p_vec + delta_p2 * h_p2, psi_e_sd[2]); 
    
    // Endogenous Model
    vector[N] h_13 = smooth_hinge_vec(eta1_vec, tau13, sharpness_k);
    vector[N] h_23 = smooth_hinge_vec(eta2_vec, tau23, sharpness_k); 
    
    vector[N] mu_eta3 = alpha[3] + 
                        (beta13 * eta1_vec + delta13 * h_13) + 
                        (beta23 * eta2_vec + delta23 * h_23); 
                          
    target += normal_lpdf(eta3_vec | mu_eta3, psi_endo_sd[1]); 
    
    // Measurement Model
    matrix[N, P] mu_Y = eta * Lambda';
    for (p in 1:P) {
      target += normal_lpdf(Y[, p] | mu_Y[, p], theta_sd[p]);
    }
  }
}
"

cat(">>> Compiling Stan Model... (This takes a minute)\n")
stan_mod <- stan_model(model_code = stan_code)
cat(">>> Model Compiled!\n\n")

# ==============================================================================
# [Single Replication Function]
# ==============================================================================
run_single_replication <- function(rep_id, stan_mod, N, lambda_gen_vec, 
                                   beta_P2, tau_P2, delta_P2,
                                   beta_13, tau_13, delta_13,
                                   beta_23, tau_23, delta_23,
                                   psi_fixed_sd, k_softplus) {
  
  library(MASS)
  library(rstan)
  library(lavaan)
  
  current_seed <- 10000 + rep_id
  set.seed(current_seed)
  
  # --- A. Data Generation (FIXED: Standardize latents to match model assumptions) ---
  p_vec_sim <- rnorm(N, 0, 1)
  
  # Generate raw latents
  eta1_raw <- 1.0 * p_vec_sim + rnorm(N, 0, psi_fixed_sd)
  
  h_P2_sim <- smooth_hinge(p_vec_sim, tau_P2, k_softplus)
  eta2_raw <- beta_P2 * p_vec_sim + delta_P2 * h_P2_sim + rnorm(N, 0, psi_fixed_sd)
  
  # [FIX 3] Standardize exogenous latents BEFORE computing endogenous model
  # This ensures tau_13 is on the correct scale
  eta1_std <- as.numeric(scale(eta1_raw))
  eta2_std <- as.numeric(scale(eta2_raw))
  
  # Endogenous model using STANDARDIZED latents
  h_13_sim <- smooth_hinge(eta1_std, tau_13, k_softplus)
  h_23_sim <- smooth_hinge(eta2_std, tau_23, k_softplus)
  eta3_mean <- (beta_13 * eta1_std + delta_13 * h_13_sim) +
    (beta_23 * eta2_std + delta_23 * h_23_sim)
  eta3_raw  <- eta3_mean + rnorm(N, 0, psi_fixed_sd)
  eta3_std <- as.numeric(scale(eta3_raw))
  
  # Measurement model using STANDARDIZED latents
  data_obs <- matrix(0, nrow=N, ncol=12)
  for(j in 1:12) {
    if(j <= 4) { factor_scores <- eta1_std; factor_idx <- j; } 
    else if(j <= 8) { factor_scores <- eta2_std; factor_idx <- j - 4; } 
    else { factor_scores <- eta3_std; factor_idx <- j - 8; }
    
    loading <- lambda_gen_vec[factor_idx]
    error_var <- 0.2 
    data_obs[,j] <- loading * factor_scores + rnorm(N, 0, sqrt(error_var))
  }
  
  data_obs <- scale(data_obs, center=TRUE, scale=FALSE)
  df_obs <- as.data.frame(data_obs)
  colnames(df_obs) <- paste0("y", 1:12)
  
  # --- B. Empirical Priors (Lavaan) ---
  lavaan_model_syntax <- '
    eta1 =~ y1 + y2 + y3 + y4
    eta2 =~ y5 + y6 + y7 + y8
    eta3 =~ y9 + y10 + y11 + y12
  '
  
  lav_fit_attempt <- try({ cfa(lavaan_model_syntax, data = df_obs, std.lv=TRUE, warn=FALSE) }, silent = TRUE)
  
  emp_mu_vec <- rep(0.5, 9); emp_sigma_vec <- rep(0.4, 9)
  if (!inherits(lav_fit_attempt, "try-error") && inspect(lav_fit_attempt, "converged")) {
    pe <- parameterEstimates(lav_fit_attempt)
    loadings_est <- pe[pe$op == "=~" & pe$lhs %in% c("eta1", "eta2", "eta3"), ]
    markers <- c("y1", "y5", "y9")
    free_loadings <- loadings_est[!loadings_est$rhs %in% markers, ]
    if (nrow(free_loadings) == 9) {
      emp_mu_vec <- free_loadings$est
      emp_sigma_vec <- free_loadings$se
      emp_sigma_vec[is.na(emp_sigma_vec)] <- 0.5
      emp_sigma_vec[emp_sigma_vec < 0.05] <- 0.05
    }
  }
  
  # --- C. Stan Estimation ---
  stan_data_list <- list(
    N = N, P = 12, M = 3,
    Y = data_obs,
    sharpness_k = k_softplus,
    slab_scale = 5.0, slab_df = 2,
    lambda_prior_mu = emp_mu_vec,
    lambda_prior_sigma = emp_sigma_vec
  )
  
  fit_attempt <- try({
    sampling(stan_mod, data = stan_data_list, chains = 4, iter = 8000, warmup = 4000, 
             cores = 1, refresh=0, seed = current_seed, 
             init = 0, 
             control = list(adapt_delta = 0.95, max_treedepth = 12))
  }, silent = TRUE)
  
  if (inherits(fit_attempt, "try-error")) {
    return(NULL)
  }
  
  # --- D. Extract Results ---
  pars_to_extract <- c(
    "beta_p2", "delta_p2", "tau_p2",
    "beta13", "delta13", "tau13",
    "beta23", "delta23", "tau23"
  )
  
  summary_obj <- summary(fit_attempt, pars = pars_to_extract, probs = c(0.025, 0.975))$summary
  res_row <- data.frame(rep_id = rep_id, converged = 1)
  
  for (par_name in pars_to_extract) {
    if (par_name %in% rownames(summary_obj)) {
      res_row[[paste0(par_name, "_mean")]] <- summary_obj[par_name, "mean"]
      res_row[[paste0(par_name, "_sd")]]   <- summary_obj[par_name, "sd"]
      res_row[[paste0(par_name, "_2.5%")]]  <- summary_obj[par_name, "2.5%"]
      res_row[[paste0(par_name, "_97.5%")]] <- summary_obj[par_name, "97.5%"]
      
      ci_lower <- summary_obj[par_name, "2.5%"]
      ci_upper <- summary_obj[par_name, "97.5%"]
      res_row[[paste0(par_name, "_sig")]] <- ifelse(ci_lower * ci_upper > 0, 1, 0)
    }
  }
  
  return(res_row)
}

# ==============================================================================
# [Main Simulation Loop - PARALLEL]
# ==============================================================================

n_cores <- parallel::detectCores()
cat(sprintf(">>> Setting up parallel processing with %d cores...\n", n_cores))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat(">>> Exporting objects to workers...\n")
clusterExport(cl, c("stan_mod", "N", "lambda_gen_vec",
                    "beta_P2", "tau_P2", "delta_P2",
                    "beta_13", "tau_13", "delta_13",
                    "beta_23", "tau_23", "delta_23",
                    "psi_fixed_sd", "k_softplus", "smooth_hinge"))

clusterEvalQ(cl, {
  library(MASS)
  library(rstan)
  library(lavaan)
  options(mc.cores = 1)
  rstan_options(auto_write = TRUE)
})

cat(">>> Starting parallel simulation...\n")
cat(sprintf(">>> Running %d replications across %d cores...\n\n", n_reps, n_cores))

start_time <- Sys.time()

results_list <- foreach(
  r = 1:n_reps,
  .packages = c("MASS", "rstan", "lavaan"),
  .combine = rbind,
  .errorhandling = "remove"
) %dopar% {
  run_single_replication(
    rep_id = r,
    stan_mod = stan_mod,
    N = N,
    lambda_gen_vec = lambda_gen_vec,
    beta_P2 = beta_P2,
    tau_P2 = tau_P2,
    delta_P2 = delta_P2,
    beta_13 = beta_13,
    tau_13 = tau_13,
    delta_13 = delta_13,
    beta_23 = beta_23,
    tau_23 = tau_23,
    delta_23 = delta_23,
    psi_fixed_sd = psi_fixed_sd,
    k_softplus = k_softplus
  )
}

end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")

stopCluster(cl)

# ==============================================================================
# [Save Results]
# ==============================================================================
cat("\n>>> Simulation completed!\n")
cat(sprintf(">>> Elapsed time: %.2f minutes\n", elapsed_time))
cat(sprintf(">>> Successful replications: %d / %d\n\n", 
            ifelse(is.null(results_list), 0, nrow(results_list)), n_reps))

if (!is.null(results_list) && nrow(results_list) > 0) {
  write.table(results_list, file = final_csv_path, sep = ",", 
              col.names = TRUE, row.names = FALSE)
  cat(sprintf(">>> Results saved to: %s\n", final_csv_path))
} else {
  cat(">>> WARNING: No successful replications to save!\n")
}

cat("\n>>> Done!\n")
