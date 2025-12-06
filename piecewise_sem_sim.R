rm(list = ls())
# setwd("YOUR_PATH_HERE")

library(MASS)
library(rstan)
library(lavaan)

# ==============================================================================
# Simulation Parameters
# ==============================================================================
N <- 1000
n_reps <- 1

lambda_marker <- 1.0
lambda_other  <- 0.8

# --- Exogenous Relations (P -> eta1, P -> eta2) ---
# Path 1 (P -> eta1) is the linear anchor (implicitly beta = 1, delta = 0)
# Path 2 (P -> eta2) keeps the piecewise signal
beta_P2  <- 0.3
tau_P2   <- 0.0
delta_P2 <- 0.8

# --- Endogenous Relations (eta1/eta2 -> eta3) ---
beta_13  <- 0.4
tau_13   <- -0.67
delta_13 <- 0.1
beta_23  <- 0.4
tau_23   <- 0.0
delta_23 <- 0.0

output_dir <- "sim_5"
if (!dir.exists(output_dir)) dir.create(output_dir)
final_csv_path <- file.path(output_dir, "sim_5.csv")

if (file.exists(final_csv_path)) {
  file.remove(final_csv_path)
}

# ==============================================================================
# Helper Functions & Population Variance Calculation
# ==============================================================================
sharp_hinge <- function(x, knot) {
  pmax(0, x - knot)
}

n_large <- 1000000
set.seed(73)

P_large <- rnorm(n_large, 0, 1)

# Path 1: Pure Linear Anchor (Assumed 1-to-1 with P)
path_P_eta1 <- 1.0 * P_large
# Path 2: Piecewise relation
path_P_eta2 <- beta_P2 * P_large + delta_P2 * sharp_hinge(P_large, tau_P2)

var_signal_eta1 <- var(path_P_eta1)
var_signal_eta2 <- var(path_P_eta2)

# Residual variances keep total variance close to 1 while avoiding zeros
psi11 <- max(1 - var_signal_eta1, 1e-6)
psi22 <- max(1 - var_signal_eta2, 1e-6)

# Generate Latents
eta1_large <- path_P_eta1 + rnorm(n_large, 0, sqrt(psi11))
eta2_large <- path_P_eta2 + rnorm(n_large, 0, sqrt(psi22))

# Generate Eta3
eta3_pred_large <- (beta_13 * eta1_large + delta_13 * sharp_hinge(eta1_large, tau_13)) +
  (beta_23 * eta2_large + delta_23 * sharp_hinge(eta2_large, tau_23))
var_signal_eta3 <- var(eta3_pred_large)
psi33 <- max(1 - var_signal_eta3, 1e-6)

cat(sprintf(">>> Fixed Psi: psi11=%.4f, psi22=%.4f, psi33=%.4f\n\n", psi11, psi22, psi33))

lambda_gen_vec <- c(lambda_marker, lambda_other, lambda_other, lambda_other)
theta_vec <- 1 - lambda_gen_vec^2
theta_vec[theta_vec < 0.05] <- 0.05

# ==============================================================================
# [Stan Model]
# ==============================================================================
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_code <- "
functions {
  vector hinge_vec(vector x, real knot) {
    int N = num_elements(x);
    vector[N] out;
    for (n in 1:N) {
      out[n] = fmax(0, x[n] - knot);
    }
    return out;
  }
}
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> M;
  matrix[N, P] Y;

  real<lower=0> slab_scale;
  real<lower=0> slab_df;
  vector[9] lambda_prior_mu;
  vector[9] lambda_prior_sigma;
}
parameters {
  vector[9] lambda_free;
  vector<lower=0>[P] theta_sd;
  real beta30;

  vector[N] p_vec;

  // Phantom path parameters
  real beta_p2;
  real tau_p2;
  real delta_raw_p2;
  real<lower=0> lambda_p2;

  // Endogenous parameters
  real beta13;
  real tau13;
  real beta23;
  real tau23;

  real delta_raw_13;
  real<lower=0> lambda_13;
  real delta_raw_23;
  real<lower=0> lambda_23;

  real<lower=0> tau_global;
  real<lower=0> caux;

  vector<lower=0>[2] psi_e_sd;
  vector<lower=0>[1] psi_endo_sd;

  matrix[N, M] eta;
}
transformed parameters {
  matrix[P, M] Lambda;
  real c2;
  real tau2;
  real lt_p2;
  real lt_13;
  real lt_23;

  real delta_p2;
  real delta13;
  real delta23;

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

  lt_p2 = sqrt(c2 * square(lambda_p2) / (c2 + tau2 * square(lambda_p2)));
  lt_13 = sqrt(c2 * square(lambda_13) / (c2 + tau2 * square(lambda_13)));
  lt_23 = sqrt(c2 * square(lambda_23) / (c2 + tau2 * square(lambda_23)));

  delta_p2 = tau_global * lt_p2 * delta_raw_p2;
  delta13  = tau_global * lt_13 * delta_raw_13;
  delta23  = tau_global * lt_23 * delta_raw_23;
}
model {
  p_vec ~ std_normal();

  lambda_free ~ normal(lambda_prior_mu, lambda_prior_sigma);
  theta_sd ~ cauchy(0, 0.5);
  beta30 ~ normal(0, 1);

  beta_p2 ~ normal(0, 0.5);
  beta13 ~ normal(0, 0.5);
  beta23 ~ normal(0, 0.5);

  tau_p2 ~ normal(0, 0.5);
  tau13  ~ normal(0, 0.5);
  tau23  ~ normal(0, 0.5);

  tau_global ~ cauchy(0, 1.0);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);

  delta_raw_p2 ~ std_normal();
  lambda_p2 ~ cauchy(0, 1);
  delta_raw_13 ~ std_normal();
  lambda_13 ~ cauchy(0, 1);
  delta_raw_23 ~ std_normal();
  lambda_23 ~ cauchy(0, 1);

  psi_e_sd ~ cauchy(0, 0.5);
  psi_endo_sd ~ cauchy(0, 0.5);

  {
    vector[N] eta1_vec = eta[, 1];
    vector[N] eta2_vec = eta[, 2];
    vector[N] eta3_vec = eta[, 3];

    vector[N] mu_eta1 = p_vec;
    vector[N] h_p2 = hinge_vec(p_vec, tau_p2);
    vector[N] mu_eta2 = beta_p2 * p_vec + delta_p2 * h_p2;

    target += normal_lpdf(eta1_vec | mu_eta1, psi_e_sd[1]);
    target += normal_lpdf(eta2_vec | mu_eta2, psi_e_sd[2]);

    vector[N] h_13 = hinge_vec(eta1_vec, tau13);
    vector[N] h_23 = hinge_vec(eta2_vec, tau23);
    vector[N] mu_eta3 = rep_vector(beta30, N) +
                        (beta13 * eta1_vec + delta13 * h_13) +
                        (beta23 * eta2_vec + delta23 * h_23);

    target += normal_lpdf(eta3_vec | mu_eta3, psi_endo_sd[1]);

    matrix[N, P] mu_Y = eta * Lambda';
    for (p in 1:P) {
      target += normal_lpdf(Y[, p] | mu_Y[, p], theta_sd[p]);
    }
  }
}
"

stan_mod <- stan_model(model_code = stan_code)

# ==============================================================================
# [Main Simulation Loop]
# ==============================================================================
lavaan_model_syntax <- '
  eta1 =~ y1 + y2 + y3 + y4
  eta2 =~ y5 + y6 + y7 + y8
  eta3 =~ y9 + y10 + y11 + y12
'

for (r in 1:n_reps) {
  cat(sprintf("--- Rep %d / %d ---\n", r, n_reps))

  current_seed <- 10000 + r
  set.seed(current_seed)

  # --- A. Data Generation ---
  P_mat_sim <- rnorm(N, 0, 1)

  path_P_eta1 <- 1.0 * P_mat_sim
  path_P_eta2 <- beta_P2 * P_mat_sim + delta_P2 * sharp_hinge(P_mat_sim, tau_P2)

  eta1 <- path_P_eta1 + rnorm(N, 0, sqrt(psi11))
  eta2 <- path_P_eta2 + rnorm(N, 0, sqrt(psi22))

  eta3_signal_from_eta1 <- beta_13 * eta1 + delta_13 * sharp_hinge(eta1, tau_13)
  eta3_signal_from_eta2 <- beta_23 * eta2 + delta_23 * sharp_hinge(eta2, tau_23)
  eta3 <- eta3_signal_from_eta1 + eta3_signal_from_eta2 + rnorm(N, 0, sqrt(psi33))

  eta_true <- cbind(eta1, eta2, eta3)

  data_obs <- matrix(0, nrow = N, ncol = 12)
  for (j in 1:12) {
    if (j <= 4) {
      factor_scores <- eta1
      factor_idx <- j
    } else if (j <= 8) {
      factor_scores <- eta2
      factor_idx <- j - 4
    } else {
      factor_scores <- eta3
      factor_idx <- j - 8
    }
    loading <- lambda_gen_vec[factor_idx]
    error_var <- theta_vec[factor_idx]
    data_obs[, j] <- loading * factor_scores + rnorm(N, 0, sqrt(error_var))
  }

  data_obs <- scale(data_obs, center = TRUE, scale = FALSE)
  df_obs <- as.data.frame(data_obs)
  colnames(df_obs) <- paste0("y", 1:12)

  # --- B. Empirical Priors ---
  emp_mu_vec <- rep(0.5, 9)
  emp_sigma_vec <- rep(0.4, 9)

  lav_fit_attempt <- try({
    cfa(lavaan_model_syntax, data = df_obs, warn = FALSE)
  }, silent = TRUE)

  use_lavaan <- FALSE
  if (!inherits(lav_fit_attempt, "try-error")) {
    if (inspect(lav_fit_attempt, "converged")) {
      est_vars <- parameterEstimates(lav_fit_attempt)
      variances <- est_vars[est_vars$op == "~~" & est_vars$lhs == est_vars$rhs, "est"]
      if (all(variances > 0.001)) {
        use_lavaan <- TRUE
      }
    }
  }

  if (use_lavaan) {
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
    slab_scale = 5.0, slab_df = 2,
    lambda_prior_mu = emp_mu_vec,
    lambda_prior_sigma = emp_sigma_vec
  )

  fit_attempt <- try({
    sampling(
      stan_mod,
      data = stan_data_list,
      chains = 4,
      iter = 6000,
      warmup = 3000,
      cores = 4,
      refresh = 0,
      seed = current_seed,
      init = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    )
  }, silent = TRUE)

  if (inherits(fit_attempt, "try-error")) {
    cat("  > [Error] Stan init failed. Skipping.\n")
    next
  }

  # --- D. Extract & Save ---
  pars_to_extract <- c(
    "beta_p2", "delta_p2", "tau_p2",
    "beta13", "delta13", "tau13",
    "beta23", "delta23", "tau23"
  )

  summary_obj <- summary(fit_attempt, pars = pars_to_extract, probs = c(0.025, 0.975))$summary
  res_row <- data.frame(rep_id = r, converged = 1)

  for (par_name in pars_to_extract) {
    if (par_name %in% rownames(summary_obj)) {
      res_row[[paste0(par_name, "_mean")]] <- summary_obj[par_name, "mean"]
      res_row[[paste0(par_name, "_sd")]] <- summary_obj[par_name, "sd"]
      res_row[[paste0(par_name, "_2.5%")]] <- summary_obj[par_name, "2.5%"]
      res_row[[paste0(par_name, "_97.5%")]] <- summary_obj[par_name, "97.5%"]

      ci_lower <- summary_obj[par_name, "2.5%"]
      ci_upper <- summary_obj[par_name, "97.5%"]
      res_row[[paste0(par_name, "_sig")]] <- ifelse(ci_lower * ci_upper > 0, 1, 0)
    }
  }

  write.table(
    res_row,
    file = final_csv_path,
    sep = ",",
    col.names = !file.exists(final_csv_path),
    row.names = FALSE,
    append = TRUE
  )

  cat(sprintf("  > Saved Rep %d to %s\n", r, final_csv_path))
}
