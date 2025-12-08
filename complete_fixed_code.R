rm(list=ls())
setwd("~/Library/CloudStorage/GoogleDrive-daylee@umd.edu/My Drive/Research/★Spline SEM/Simulation")

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

lambda_gen_vec <- c(1.0, 0.8, 0.8, 0.8)

beta_P2  <- 0.3; tau_P2   <- 0.0; delta_P2 <- 0.5
beta_13  <- 0.4; tau_13   <- -0.67; delta_13 <- 0.5
beta_23  <- 0.4; tau_23   <- 0.0; delta_23 <- 0.0

psi_fixed_sd  <- sqrt(0.2)
smooth_eps    <- 0.1

output_dir <- "sim_6"
if(!dir.exists(output_dir)) dir.create(output_dir)
final_csv_path <- file.path(output_dir, "sim_6.csv")
if(file.exists(final_csv_path)) file.remove(final_csv_path)

# ==============================================================================
# Helper Function: Zero-Centered Harring Spline (R Version)
# ==============================================================================
harring_spline_gen <- function(x, intercept_at_0, beta, delta, knot, eps) {
  alpha2 <- beta + delta / 2.0
  alpha3 <- delta / 2.0
  
  raw_x <- alpha2 * (x - knot) + alpha3 * sqrt((x - knot)^2 + eps)
  raw_0 <- alpha2 * (0 - knot) + alpha3 * sqrt((0 - knot)^2 + eps)
  
  y <- intercept_at_0 + (raw_x - raw_0)
  return(y)
}

# ==============================================================================
# [Stan Model] - IDENTIFICATION FIXES APPLIED
# ==============================================================================
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_code <- "
functions {
  vector harring_spline_centered(vector x, real alpha2, real alpha3, real gamma, real eps) {
    int N = num_elements(x);
    vector[N] raw_vals;
    real val_at_zero;
    
    raw_vals = alpha2 * (x - gamma) + alpha3 * sqrt(square(x - gamma) + eps);
    val_at_zero = alpha2 * (0 - gamma) + alpha3 * sqrt(square(0 - gamma) + eps);
    
    return raw_vals - val_at_zero;
  }
}

data {
  int<lower=1> N;             
  int<lower=1> P;             
  int<lower=1> M;             
  matrix[N, P] Y;             
  real<lower=0> smooth_eps;   
  
  real<lower=0> slab_scale; 
  real<lower=0> slab_df;
  
  vector[9] lambda_prior_mu;
  vector[9] lambda_prior_sigma;
}

parameters {
  vector[9] lambda_free;
  vector<lower=0>[P] theta_sd; 
  
  // [FIX] alpha_zero[1] 제거 - eta1의 intercept는 0으로 고정
  real alpha_zero_2;  
  real alpha_zero_3;  
  
  vector[N] p_vec; 
  
  real alpha2_p2; 
  real alpha2_13;
  real alpha2_23;

  real tau_p2;  
  real tau13;   
  real tau23;   
  
  real alpha3_raw_p2; real<lower=0> lambda_p2;
  real alpha3_raw_13; real<lower=0> lambda_13; 
  real alpha3_raw_23; real<lower=0> lambda_23; 
  
  real<lower=0> tau_global; 
  real<lower=0> caux; 
  
  vector<lower=0>[2] psi_e_sd;    
  vector<lower=0>[1] psi_endo_sd; 
  
  // [FIX] eta_raw: 원시 잠재변수 (평균 제약 전)
  matrix[N, M] eta_raw;  
}

transformed parameters {
  matrix[P, M] Lambda;
  // [FIX] eta: 각 열의 평균이 0인 잠재변수
  matrix[N, M] eta;  
  
  real c2 = square(slab_scale) * caux; 
  real tau_sq = square(tau_global);
  
  real lt_p2 = sqrt( c2 * square(lambda_p2) / (c2 + tau_sq * square(lambda_p2)) );
  real lt_13 = sqrt( c2 * square(lambda_13) / (c2 + tau_sq * square(lambda_13)) );
  real lt_23 = sqrt( c2 * square(lambda_23) / (c2 + tau_sq * square(lambda_23)) );
  
  real alpha3_p2 = tau_global * lt_p2 * alpha3_raw_p2;
  real alpha3_13 = tau_global * lt_13 * alpha3_raw_13;
  real alpha3_23 = tau_global * lt_23 * alpha3_raw_23;
  
  // [FIX] eta의 각 열(잠재변수)의 평균을 0으로 제약
  for (m in 1:M) {
    real mean_eta_m = mean(eta_raw[, m]);
    eta[, m] = eta_raw[, m] - mean_eta_m;
  }
  
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
}

model {
  // [FIX] p_vec의 평균을 0에 가깝도록 soft constraint
  p_vec ~ std_normal();
  sum(p_vec) ~ normal(0, 0.1 * sqrt(N));
  
  lambda_free ~ normal(lambda_prior_mu, lambda_prior_sigma);
  target += normal_lpdf(theta_sd | 0.45, 0.15);

  alpha_zero_2 ~ normal(0, 1);
  alpha_zero_3 ~ normal(0, 1);
  
  alpha2_p2 ~ normal(0, 0.5);
  alpha2_13 ~ normal(0, 0.5);
  alpha2_23 ~ normal(0, 0.5);
  
  tau_p2  ~ normal(0, 0.5); 
  tau13   ~ normal(0, 0.5); 
  tau23   ~ normal(0, 0.5); 
  
  tau_global ~ cauchy(0, 1.0); 
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df); 
  
  alpha3_raw_p2 ~ std_normal(); lambda_p2 ~ cauchy(0, 1);
  alpha3_raw_13 ~ std_normal(); lambda_13 ~ cauchy(0, 1);
  alpha3_raw_23 ~ std_normal(); lambda_23 ~ cauchy(0, 1); 
  
  psi_e_sd ~ normal(0.4, 0.2);
  psi_endo_sd ~ normal(0.4, 0.2);
  
  // eta_raw에 대한 prior
  for (m in 1:M) {
    eta_raw[, m] ~ std_normal();
  }
  
  // --- Likelihood ---
  {
    vector[N] eta1_vec = eta[, 1];
    vector[N] eta2_vec = eta[, 2];
    vector[N] eta3_vec = eta[, 3];
    
    // [FIX] eta1: intercept = 0으로 고정
    target += normal_lpdf(eta1_vec | p_vec, psi_e_sd[1]);
    
    vector[N] spline_p2 = harring_spline_centered(p_vec, alpha2_p2, alpha3_p2, tau_p2, smooth_eps);
    target += normal_lpdf(eta2_vec | alpha_zero_2 + spline_p2, psi_e_sd[2]); 
    
    vector[N] spline_13 = harring_spline_centered(eta1_vec, alpha2_13, alpha3_13, tau13, smooth_eps);
    vector[N] spline_23 = harring_spline_centered(eta2_vec, alpha2_23, alpha3_23, tau23, smooth_eps);
    
    vector[N] mu_eta3 = alpha_zero_3 + spline_13 + spline_23;
    target += normal_lpdf(eta3_vec | mu_eta3, psi_endo_sd[1]); 
    
    matrix[N, P] mu_Y = eta * Lambda';
    for (p in 1:P) {
      target += normal_lpdf(Y[, p] | mu_Y[, p], theta_sd[p]);
    }
  }
}

generated quantities {
  real beta_p2 = alpha2_p2 - alpha3_p2;
  real delta_p2 = 2 * alpha3_p2;
  
  real beta13 = alpha2_13 - alpha3_13;
  real delta13 = 2 * alpha3_13;
  
  real beta23 = alpha2_23 - alpha3_23;
  real delta23 = 2 * alpha3_23;
}
"

cat(">>> Compiling Stan Model (with Identification Fixes)...\n")
stan_mod <- stan_model(model_code = stan_code)
cat(">>> Model Compiled!\n\n")

# ==============================================================================
# [Single Replication Function]
# ==============================================================================
run_single_replication <- function(rep_id, stan_mod, N, lambda_gen_vec, 
                                   beta_P2, tau_P2, delta_P2,
                                   beta_13, tau_13, delta_13,
                                   beta_23, tau_23, delta_23,
                                   psi_fixed_sd, smooth_eps, harring_gen_func) {
  
  library(MASS)
  library(rstan)
  library(lavaan)
  
  current_seed <- 10000 + rep_id
  set.seed(current_seed)
  
  # --- A. Data Generation ---
  p_vec_sim <- rnorm(N, 0, 1)
  
  # Eta1: intercept = 0 (고정)
  eta1_raw <- 0 + 1.0 * p_vec_sim + rnorm(N, 0, psi_fixed_sd)
  
  # Eta2
  eta2_mean <- harring_gen_func(p_vec_sim, intercept_at_0=0, 
                                beta=beta_P2, delta=delta_P2, knot=tau_P2, eps=smooth_eps)
  eta2_raw <- eta2_mean + rnorm(N, 0, psi_fixed_sd)
  
  # Eta3
  part1 <- harring_gen_func(eta1_raw, intercept_at_0=0, 
                            beta=beta_13, delta=delta_13, knot=tau_13, eps=smooth_eps)
  part2 <- harring_gen_func(eta2_raw, intercept_at_0=0, 
                            beta=beta_23, delta=delta_23, knot=tau_23, eps=smooth_eps)
  eta3_raw <- (part1 + part2) + rnorm(N, 0, psi_fixed_sd)
  
  # Measurement
  data_obs <- matrix(0, nrow=N, ncol=12)
  for(j in 1:12) {
    if(j <= 4) { factor_scores <- eta1_raw; factor_idx <- j; } 
    else if(j <= 8) { factor_scores <- eta2_raw; factor_idx <- j - 4; } 
    else { factor_scores <- eta3_raw; factor_idx <- j - 8; }
    
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
    smooth_eps = smooth_eps,
    slab_scale = 5.0, slab_df = 2,
    lambda_prior_mu = emp_mu_vec,
    lambda_prior_sigma = emp_sigma_vec
  )
  
  fit_attempt <- try({
    sampling(stan_mod, data = stan_data_list, chains = 4, iter = 6000, warmup = 3000, 
             cores = 1, refresh=0, seed = current_seed, 
             init = 0, 
             control = list(adapt_delta = 0.95, max_treedepth = 12))
  }, silent = TRUE)
  
  if (inherits(fit_attempt, "try-error")) return(NULL)
  
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
      res_row[[paste0(par_name, "_sd")]]    <- summary_obj[par_name, "sd"]
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
# [Main Simulation Loop]
# ==============================================================================
n_cores <- parallel::detectCores()
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c("stan_mod", "N", "lambda_gen_vec",
                    "beta_P2", "tau_P2", "delta_P2",
                    "beta_13", "tau_13", "delta_13",
                    "beta_23", "tau_23", "delta_23",
                    "psi_fixed_sd", "smooth_eps", "harring_spline_gen"))

clusterEvalQ(cl, {
  library(MASS)
  library(rstan)
  library(lavaan)
  options(mc.cores = 1)
  rstan_options(auto_write = TRUE)
})

cat(">>> Starting simulation with IDENTIFICATION FIXES...\n")
cat(">>> Fixes: alpha_zero[1]=0, eta means=0, p_vec sum constraint\n\n")

results_list <- foreach(
  r = 1:n_reps,
  .packages = c("MASS", "rstan", "lavaan"),
  .combine = rbind,
  .errorhandling = "remove"
) %dopar% {
  run_single_replication(
    rep_id = r, stan_mod = stan_mod, N = N,
    lambda_gen_vec = lambda_gen_vec,
    beta_P2 = beta_P2, tau_P2 = tau_P2, delta_P2 = delta_P2,
    beta_13 = beta_13, tau_13 = tau_13, delta_13 = delta_13,
    beta_23 = beta_23, tau_23 = tau_23, delta_23 = delta_23,
    psi_fixed_sd = psi_fixed_sd, smooth_eps = smooth_eps,
    harring_gen_func = harring_spline_gen
  )
}

stopCluster(cl)

if (!is.null(results_list) && nrow(results_list) > 0) {
  write.table(results_list, file = final_csv_path, sep = ",", col.names = TRUE, row.names = FALSE)
  cat(sprintf(">>> Results saved to: %s\n", final_csv_path))
} else {
  cat(">>> WARNING: No successful replications!\n")
}

cat("\n>>> Done!\n")
