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
# [Main Simulation Loop]
# ==============================================================================

cat(">>> Starting Monte Carlo Simulation (", n_reps, " replications)...\n")

# Measurement Model Parameters
# 3 Factors, 4 Indicators each = 12 items
# Lambda: First item is marker (1.0), others are lambda_val/lambda1 ratio? 
# Usually in simulation for standardized factors:
# If Latent Var=1, and we want loading=lambda_val, then Loading = lambda_val.
# Residual (Theta) = 1 - lambda_val^2.

lambda_vec <- rep(lambda_val, 4) 
# For CFA generation, usually we set the first loading to 1.0 later for estimation, 
# but for data generation, we can generate with true standardized loadings.
# Or if you want to strictly follow "First indicator = 1.0 (Marker)" in generation:
# Since Latent Var=1, if L1=1.0, then Item Var = 1^2*1 + Theta. This makes item not standardized.
# To keep Items standardized (Var=1):
# We use uniform loadings (e.g., 0.8) for generation. When analyzing, we fix one to 1.0.

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
Theta_diag <- rep(c(theta_vec, theta_vec, theta_vec)) # 12 items
Theta_mat  <- diag(Theta_diag)


library(MASS)
library(rstan)

# 병렬 처리 설정 (PC 사양에 맞게 조정)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ==============================================================================
# 1. Define Stan Model (Compile ONCE outside the loop)
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
# 2. Simulation Settings & Setup
# ==============================================================================

# User Inputs
N <- 200
n_reps <- 500  # 총 반복 횟수

# True Parameters (Standardized Generation)
lambda_val <- 0.8 
beta_P <- 0.3; tau_P <- 0.0; delta_P <- 0.2
beta_13 <- 0.4; tau_13 <- 0.5; delta_13 <- 0.3
beta_23 <- 0.4; delta_23 <- 0.0; tau_23 <- 0.0 # Linear

# Helper for Hinge
sharp_hinge <- function(x, knot) pmax(0, x - knot)

# --- Psi Calculation for Standardization ---
set.seed(123)
n_large <- 100000
P_large <- rnorm(n_large, 0, 1)
path_P_large <- beta_P * P_large + delta_P * sharp_hinge(P_large, tau_P)
var_path_P <- var(path_P_large)
psi11 <- 1 - var_path_P
psi22 <- 1 - var_path_P
if(psi11 < 0) stop("Phantom effect too large")

eta1_l <- path_P_large + rnorm(n_large, 0, sqrt(psi11))
eta2_l <- path_P_large + rnorm(n_large, 0, sqrt(psi22))
eta3_pred_l <- (beta_13 * eta1_l + delta_13 * sharp_hinge(eta1_l, tau_13)) +
  (beta_23 * eta2_l + delta_23 * sharp_hinge(eta2_l, tau_23))
psi33 <- 1 - var(eta3_pred_l)
if(psi33 < 0) stop("Structural effect too large")

# Measurement setup
lambda_vec <- rep(lambda_val, 4)
Theta_diag <- rep(1 - lambda_vec^2, 3)
Theta_mat <- diag(Theta_diag)
Lambda_mat <- matrix(0, 12, 3)
Lambda_mat[1:4,1] <- lambda_vec
Lambda_mat[5:8,2] <- lambda_vec
Lambda_mat[9:12,3] <- lambda_vec

# Result Storage
results_list <- list()

# ==============================================================================
# 3. Main Loop: Generate -> Fit -> Save
# ==============================================================================

cat(">>> Starting Monte Carlo Loop (", n_reps, " reps)...\n")
set.seed(999)

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
      control = list(adapt_delta = 0.90)
    )
  }, silent = TRUE)
  
  if (inherits(fit_attempt, "try-error")) {
    cat("   [Rep", r, "] Stan Error. Skipping...\n")
    next
  }
  
  # --- C. Extract Results ---
  # 주요 파라미터만 추출하여 저장
  summary_fit <- summary(fit_attempt, 
                         pars = c("beta_p", "tau_p", "delta_p_shared",
                                  "beta13", "tau13", "delta13",
                                  "beta23"))$summary
  
  # 결과 데이터프레임 생성 (Posterior Mean, 2.5%, 97.5%)
  # Transpose to single row
  res_row <- as.data.frame(t(summary_fit[, "mean"]))
  colnames(res_row) <- paste0(colnames(res_row), "_mean")
  
  # Add simulation info
  res_row$rep_id <- r
  res_row$converged <- ifelse(all(summary_fit[, "Rhat"] < 1.1), 1, 0)
  
  results_list[[r]] <- res_row
  
  if (r %% 10 == 0) cat("   ... replicate", r, "completed.\n")
}

# ==============================================================================
# 4. Save Final Results
# ==============================================================================

final_results <- do.call(rbind, results_list)
write.csv(final_results, "mc_simulation_results.csv", row.names = FALSE)

cat(">>> Simulation Finished! Results saved to 'mc_simulation_results.csv'\n")
print(head(final_results))
