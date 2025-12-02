# ==============================================================================
# Monte Carlo Simulation for Piecewise Structural Relations
# Three Latent Variables Model with Piecewise Relationships
# ==============================================================================

# ==============================================================================
# [User Input: Manipulation Parameters]
# ==============================================================================

# 1. Sample Size
N <- 200 

# 2. Factor Loadings (Standardized)
lambda_val <- 0.8  # 모든 loading을 0.8로 가정 (첫번째 indicator는 marker로 1.0 변환됨)

# 3. Phantom Variable Structure (P -> eta1, P -> eta2)
# Common cause creating correlation between eta1 and eta2
beta_P  <- 0.3   # Baseline slope
tau_P   <- 0.0   # Knot location
delta_P <- 0.2   # Slope difference (Non-linearity)

# 4. eta3 ~ eta1 Path (Piecewise)
beta_13  <- 0.4  # Baseline slope
tau_13   <- 0.5  # Knot location for eta1 -> eta3
delta_13 <- 0.3  # Slope difference

# 5. eta3 ~ eta2 Path (Linear, fixed slope = 0.4)
beta_23  <- 0.4  # Linear slope fixed
tau_23   <- 0.0  # (Not used since linear)
delta_23 <- 0.0  # 0.0 implies Linear relationship

# 6. Replication Settings
n_reps <- 500
output_dir <- "sim_data_1" # CSV 저장 폴더명

# ==============================================================================
# [Helper Functions & Setup]
# ==============================================================================

if(!dir.exists(output_dir)) dir.create(output_dir)

sharp_hinge <- function(x, knot) {
  pmax(0, x - knot)
}

# --- A. Large Sample Calculation (Standardization Logic) ---
# 입력된 파라미터(beta, delta)가 표준화 계수가 되도록 잔차 분산(psi)을 역산

cat(">>> Calculating Residual Variances for Standardization...\n")

n_large <- 500000
set.seed(12345) # 파라미터 계산용 고정 시드

# [Step 1] Exogenous Latent Variables (eta1, eta2)
# Phantom Variable P ~ N(0,1)
P_large <- rnorm(n_large, 0, 1)

# Path from Phantom
path_P_large <- beta_P * P_large + delta_P * sharp_hinge(P_large, tau_P)
var_path_P   <- var(path_P_large)

# eta1, eta2의 분산이 1이 되기 위한 잔차 분산 (psi11, psi22)
psi11 <- 1 - var_path_P
psi22 <- 1 - var_path_P

if(psi11 <= 0) stop("Error: Phantom path coefficients are too large (Explains > 100% variance)")

# Generate large eta1, eta2 for next step
e1_large <- rnorm(n_large, 0, sqrt(psi11))
e2_large <- rnorm(n_large, 0, sqrt(psi22))
eta1_large <- path_P_large + e1_large
eta2_large <- path_P_large + e2_large

# [Step 2] Endogenous Latent Variable (eta3)
# eta3 = beta1*eta1 + delta1*hinge(eta1) + beta2*eta2 + delta2*hinge(eta2) + error
# Note: delta_23 is 0, so the second part is linear
eta3_pred_large <- (beta_13 * eta1_large + delta_13 * sharp_hinge(eta1_large, tau_13)) +
  (beta_23 * eta2_large + delta_23 * sharp_hinge(eta2_large, tau_23))

var_eta3_expl <- var(eta3_pred_large)
psi33 <- 1 - var_eta3_expl

if(psi33 <= 0) stop("Error: Structural path coefficients for eta3 are too large")

cat("   Psi11 (eta1 resid):", round(psi11, 4), "\n")
cat("   Psi22 (eta2 resid):", round(psi22, 4), "\n")
cat("   Psi33 (eta3 resid):", round(psi33, 4), "\n")
cat(">>> Standardization setup complete.\n\n")

# ==============================================================================
# [Measurement Model Setup]
# ==============================================================================

# Measurement Model Parameters
# 3 Factors, 4 Indicators each = 12 items
lambda_vec <- rep(lambda_val, 4) 

Lambda_mat <- matrix(0, nrow=12, ncol=3)
# Factor 1 (eta1): Items 1-4
Lambda_mat[1:4, 1] <- lambda_vec
# Factor 2 (eta2): Items 5-8
Lambda_mat[5:8, 2] <- lambda_vec
# Factor 3 (eta3): Items 9-12
Lambda_mat[9:12, 3] <- lambda_vec

# Measurement Error Variances (Theta)
# Var(X) = L^2 * Var(F) + Theta. Assuming Var(X)=1, Var(F)=1.
theta_vec <- 1 - lambda_vec^2
Theta_diag <- rep(theta_vec, 3) # 12 items
Theta_mat  <- diag(Theta_diag)

# ==============================================================================
# [Load Libraries]
# ==============================================================================

library(MASS)
library(rstan)

# 병렬 처리 설정 (PC 사양에 맞게 조정)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ==============================================================================
# [Stan Model Definition & Compilation]
# ==============================================================================

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
  matrix[N, 1] P_mat; 
  
  real<lower=0> slab_scale; 
  real<lower=0> slab_df;
}

parameters {
  // Measurement
  vector[9] lambda_free;
  vector<lower=0>[P] theta_sd;
  real beta30; 
  
  // Structural
  real beta_p; real tau_p;  
  real beta13; real tau13; 
  real beta23;             
  
  // RHS / Horseshoe
  real delta_raw_p; real<lower=0> lambda_p;
  real delta_raw_13; real<lower=0> lambda_13; 

  real<lower=0> tau_global; 
  real<lower=0> caux; 
  
  vector<lower=0>[2] psi_e_sd;    
  vector<lower=0>[1] psi_endo_sd; 
  
  matrix[N, M] eta; 
}

transformed parameters {
  matrix[P, M] Lambda;
  real c2; real tau2;
  real lt_p; real lt_13;
  real delta_p_shared;
  real delta13;

  // Lambda Construction
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

  // RHS Logic
  c2 = square(slab_scale) * caux; 
  tau2 = square(tau_global);

  lt_p  = sqrt( c2 * square(lambda_p)  / (c2 + tau2 * square(lambda_p)) );
  lt_13 = sqrt( c2 * square(lambda_13) / (c2 + tau2 * square(lambda_13)) );
  
  delta_p_shared = tau_global * lt_p  * delta_raw_p;
  delta13        = tau_global * lt_13 * delta_raw_13;
}

model {
  real k_softplus = sharpness_k;
  
  // Priors
  lambda_free ~ normal(0.8, 0.5);
  theta_sd ~ cauchy(0, 0.5);
  beta30 ~ normal(0, 1);
  
  beta_p ~ normal(0, 0.5);
  beta13 ~ normal(0, 0.5);
  beta23 ~ normal(0, 0.5);
  
  tau_p  ~ normal(0, 1); 
  tau13  ~ normal(0, 1); 
  
  tau_global ~ cauchy(0, 0.5); 
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df); 
  delta_raw_p  ~ std_normal(); lambda_p  ~ cauchy(0, 1);
  delta_raw_13 ~ std_normal(); lambda_13 ~ cauchy(0, 1);
  
  psi_e_sd ~ cauchy(0, 0.5);
  psi_endo_sd ~ cauchy(0, 0.5);
  
  // Likelihood
  {
    vector[N] eta1_vec = eta[, 1];
    vector[N] eta2_vec = eta[, 2];
    vector[N] eta3_vec = eta[, 3];
    vector[N] p_vec = P_mat[, 1];
    
    // Exogenous
    vector[N] h_p = smooth_hinge_vec(p_vec, tau_p, k_softplus); 
    vector[N] mu_eta1 = beta_p * p_vec + delta_p_shared * h_p; 
    vector[N] mu_eta2 = beta_p * p_vec + delta_p_shared * h_p; 

    target += normal_lpdf(eta1_vec | mu_eta1, psi_e_sd[1]);
    target += normal_lpdf(eta2_vec | mu_eta2, psi_e_sd[2]); 
    
    // Endogenous
    vector[N] h_13 = smooth_hinge_vec(eta1_vec, tau13, k_softplus);
    vector[N] mu_eta3 = rep_vector(beta30, N) + 
                        (beta13 * eta1_vec + delta13 * h_13) + 
                        (beta23 * eta2_vec); 
                        
    target += normal_lpdf(eta3_vec | mu_eta3, psi_endo_sd[1]); 
    
    // Measurement
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
# [Main Simulation Loop]
# ==============================================================================

cat(">>> Starting Monte Carlo Loop (", n_reps, " replications)...\n")
cat(">>> Results will be saved to:", output_dir, "/\n\n")

# Set seed for reproducibility (but allow variation across replications)
set.seed(999)

# Track failed replications
failed_reps <- 0
successful_reps <- 0

for (r in 1:n_reps) {
  
  # --- A. Data Generation ---
  P_mat_sim <- as.matrix(rnorm(N, 0, 1)) # Input for Stan
  path_P <- beta_P * P_mat_sim + delta_P * sharp_hinge(P_mat_sim, tau_P)
  
  eta1 <- path_P + rnorm(N, 0, sqrt(psi11))
  eta2 <- path_P + rnorm(N, 0, sqrt(psi22))
  
  # eta3 generation
  eta3_str <- (beta_13 * eta1 + delta_13 * sharp_hinge(eta1, tau_13)) +
    (beta_23 * eta2 + delta_23 * sharp_hinge(eta2, tau_23))
  eta3 <- eta3_str + rnorm(N, 0, sqrt(psi33))
  
  eta_true <- cbind(eta1, eta2, eta3)
  errors <- mvrnorm(N, mu = rep(0, 12), Sigma = Theta_mat)
  data_obs <- eta_true %*% t(Lambda_mat) + errors
  
  # --- B. Stan Estimation ---
  
  stan_data_list <- list(
    N = N, P = 12, M = 3,
    Y = data_obs,
    sharpness_k = 20,
    P_mat = P_mat_sim,  # Pass the generated Phantom variable
    slab_scale = 2.5, slab_df = 4
  )
  
  # Try-Catch for stability (skip failed runs)
  fit_attempt <- try({
    sampling(
      stan_mod, 
      data = stan_data_list,
      chains = 2,       # Speed up: 2 chains
      iter = 2000,      # Speed up: 2000 iter (1000 warmup)
      warmup = 1000,
      cores = 2,        # Parallel chains
      refresh = 0,      # Silence output
      control = list(adapt_delta = 0.90, max_treedepth = 12)
    )
  }, silent = TRUE)
  
  if (inherits(fit_attempt, "try-error")) {
    cat("   [Rep", r, "] Stan Error. Skipping...\n")
    failed_reps <- failed_reps + 1
    next
  }
  
  # Check for convergence issues
  if (is.null(fit_attempt)) {
    cat("   [Rep", r, "] Fit is NULL. Skipping...\n")
    failed_reps <- failed_reps + 1
    next
  }
  
  # --- C. Extract Results ---
  # 주요 파라미터 추출
  target_pars <- c("beta_p", "tau_p", "delta_p_shared",
                   "beta13", "tau13", "delta13",
                   "beta23")
  
  summary_fit <- try({
    summary(fit_attempt, pars = target_pars)$summary
  }, silent = TRUE)
  
  if (inherits(summary_fit, "try-error") || is.null(summary_fit)) {
    cat("   [Rep", r, "] Summary extraction failed. Skipping...\n")
    failed_reps <- failed_reps + 1
    next
  }
  
  # Check convergence (Rhat < 1.1)
  rhat_values <- summary_fit[, "Rhat"]
  converged <- ifelse(all(rhat_values < 1.1, na.rm = TRUE), 1, 0)
  
  # 결과 데이터프레임 생성
  res_row <- data.frame(
    rep_id = r,
    converged = converged,
    # Posterior means
    beta_p_mean = summary_fit["beta_p", "mean"],
    tau_p_mean = summary_fit["tau_p", "mean"],
    delta_p_shared_mean = summary_fit["delta_p_shared", "mean"],
    beta13_mean = summary_fit["beta13", "mean"],
    tau13_mean = summary_fit["tau13", "mean"],
    delta13_mean = summary_fit["delta13", "mean"],
    beta23_mean = summary_fit["beta23", "mean"],
    # Posterior SDs
    beta_p_sd = summary_fit["beta_p", "sd"],
    tau_p_sd = summary_fit["tau_p", "sd"],
    delta_p_shared_sd = summary_fit["delta_p_shared", "sd"],
    beta13_sd = summary_fit["beta13", "sd"],
    tau13_sd = summary_fit["tau13", "sd"],
    delta13_sd = summary_fit["delta13", "sd"],
    beta23_sd = summary_fit["beta23", "sd"],
    # 95% CI Lower bounds
    beta_p_lower = summary_fit["beta_p", "2.5%"],
    tau_p_lower = summary_fit["tau_p", "2.5%"],
    delta_p_shared_lower = summary_fit["delta_p_shared", "2.5%"],
    beta13_lower = summary_fit["beta13", "2.5%"],
    tau13_lower = summary_fit["tau13", "2.5%"],
    delta13_lower = summary_fit["delta13", "2.5%"],
    beta23_lower = summary_fit["beta23", "2.5%"],
    # 95% CI Upper bounds
    beta_p_upper = summary_fit["beta_p", "97.5%"],
    tau_p_upper = summary_fit["tau_p", "97.5%"],
    delta_p_shared_upper = summary_fit["delta_p_shared", "97.5%"],
    beta13_upper = summary_fit["beta13", "97.5%"],
    tau13_upper = summary_fit["tau13", "97.5%"],
    delta13_upper = summary_fit["delta13", "97.5%"],
    beta23_upper = summary_fit["beta23", "97.5%"],
    # Rhat values
    beta_p_rhat = summary_fit["beta_p", "Rhat"],
    tau_p_rhat = summary_fit["tau_p", "Rhat"],
    delta_p_shared_rhat = summary_fit["delta_p_shared", "Rhat"],
    beta13_rhat = summary_fit["beta13", "Rhat"],
    tau13_rhat = summary_fit["tau13", "Rhat"],
    delta13_rhat = summary_fit["delta13", "Rhat"],
    beta23_rhat = summary_fit["beta23", "Rhat"],
    stringsAsFactors = FALSE
  )
  
  # Add true parameter values for comparison
  res_row$beta_p_true <- beta_P
  res_row$tau_p_true <- tau_P
  res_row$delta_p_true <- delta_P
  res_row$beta13_true <- beta_13
  res_row$tau13_true <- tau_13
  res_row$delta13_true <- delta_13
  res_row$beta23_true <- beta_23
  
  # --- D. Save Individual Replication Result ---
  # 각 replication마다 결과를 저장 (사용자 요청사항)
  rep_file <- file.path(output_dir, paste0("rep_", sprintf("%04d", r), ".csv"))
  write.csv(res_row, rep_file, row.names = FALSE)
  
  successful_reps <- successful_reps + 1
  
  if (r %% 10 == 0) {
    cat("   ... replicate", r, "completed. (Success:", successful_reps, 
        ", Failed:", failed_reps, ")\n")
  }
}

# ==============================================================================
# [Combine and Save Final Results]
# ==============================================================================

cat("\n>>> Combining all replication results...\n")

# Read all individual files and combine
rep_files <- list.files(output_dir, pattern = "^rep_.*\\.csv$", full.names = TRUE)

if (length(rep_files) > 0) {
  all_results <- do.call(rbind, lapply(rep_files, read.csv))
  
  # Sort by rep_id
  all_results <- all_results[order(all_results$rep_id), ]
  
  # Save combined results
  final_file <- file.path(output_dir, "mc_simulation_results_combined.csv")
  write.csv(all_results, final_file, row.names = FALSE)
  
  cat(">>> Simulation Finished!\n")
  cat("   Total replications attempted:", n_reps, "\n")
  cat("   Successful replications:", successful_reps, "\n")
  cat("   Failed replications:", failed_reps, "\n")
  cat("   Individual results saved to:", output_dir, "/\n")
  cat("   Combined results saved to:", final_file, "\n\n")
  
  # Print summary statistics
  cat(">>> Summary Statistics (Posterior Means):\n")
  print(summary(all_results[, grep("_mean$", colnames(all_results))]))
  
  # Convergence rate
  if ("converged" %in% colnames(all_results)) {
    conv_rate <- mean(all_results$converged, na.rm = TRUE)
    cat("\n>>> Convergence Rate:", round(conv_rate * 100, 2), "%\n")
  }
  
  print(head(all_results))
} else {
  cat(">>> ERROR: No successful replications found!\n")
}
