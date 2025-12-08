# ==============================================================================
# IDENTIFICATION FIXES
# ==============================================================================
# 
# 문제 1: alpha_zero[1]과 p_vec이 함께 식별되지 않음
#   - eta1 = alpha_zero[1] + p_vec + error
#   - 해결: alpha_zero[1] = 0으로 고정 (p_vec이 이미 N(0,1)이므로)
#
# 문제 2: eta의 평균이 제약되지 않음
#   - 해결: transformed parameters에서 eta의 평균을 0으로 제약
#
# 문제 3: p_vec의 평균이 명시적으로 제약되지 않음
#   - 해결: p_vec의 합을 0으로 제약 (soft constraint via prior는 이미 있음)

# 수정된 Stan 모델 코드:

stan_code_fixed <- "
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
  
  // [FIX 1] alpha_zero[1] 제거 (고정값 0 사용)
  // alpha_zero[1]은 0으로 고정되어야 함 (p_vec이 N(0,1)이므로)
  real alpha_zero_2;  // eta2의 intercept만 추정
  real alpha_zero_3;  // eta3의 intercept만 추정
  
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
  
  matrix[N, M] eta_raw;  // [FIX 2] 원시 eta (나중에 평균 0으로 변환)
}

transformed parameters {
  matrix[P, M] Lambda;
  matrix[N, M] eta;  // [FIX 2] 평균이 0인 eta
  
  real c2 = square(slab_scale) * caux; 
  real tau_sq = square(tau_global);
  
  real lt_p2 = sqrt( c2 * square(lambda_p2) / (c2 + tau_sq * square(lambda_p2)) );
  real lt_13 = sqrt( c2 * square(lambda_13) / (c2 + tau_sq * square(lambda_13)) );
  real lt_23 = sqrt( c2 * square(lambda_23) / (c2 + tau_sq * square(lambda_23)) );
  
  real alpha3_p2 = tau_global * lt_p2 * alpha3_raw_p2;
  real alpha3_13 = tau_global * lt_13 * alpha3_raw_13;
  real alpha3_23 = tau_global * lt_23 * alpha3_raw_23;
  
  // [FIX 2] eta의 각 열(잠재변수)의 평균을 0으로 제약
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
  // [FIX 3] p_vec의 평균을 0으로 제약 (soft constraint)
  // std_normal() prior가 이미 있지만, 추가로 합이 0에 가깝도록 제약
  p_vec ~ std_normal();
  // Soft constraint: sum(p_vec) should be close to 0
  sum(p_vec) ~ normal(0, 0.1 * sqrt(N));  // N이 클수록 더 엄격
  
  lambda_free ~ normal(lambda_prior_mu, lambda_prior_sigma);
  target += normal_lpdf(theta_sd | 0.45, 0.15);

  // Structural Intercepts
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
  
  // eta_raw에 대한 prior (평균은 나중에 제거됨)
  for (m in 1:M) {
    eta_raw[, m] ~ std_normal();
  }
  
  // --- Likelihood ---
  {
    vector[N] eta1_vec = eta[, 1];
    vector[N] eta2_vec = eta[, 2];
    vector[N] eta3_vec = eta[, 3];
    
    // [FIX 1] eta1: alpha_zero[1] = 0으로 고정
    target += normal_lpdf(eta1_vec | 0 + p_vec, psi_e_sd[1]);
    
    // eta2: alpha_zero_2 사용
    vector[N] spline_p2 = harring_spline_centered(p_vec, alpha2_p2, alpha3_p2, tau_p2, smooth_eps);
    target += normal_lpdf(eta2_vec | alpha_zero_2 + spline_p2, psi_e_sd[2]); 
    
    // eta3: alpha_zero_3 사용
    vector[N] spline_13 = harring_spline_centered(eta1_vec, alpha2_13, alpha3_13, tau13, smooth_eps);
    vector[N] spline_23 = harring_spline_centered(eta2_vec, alpha2_23, alpha3_23, tau23, smooth_eps);
    
    vector[N] mu_eta3 = alpha_zero_3 + spline_13 + spline_23;
    target += normal_lpdf(eta3_vec | mu_eta3, psi_endo_sd[1]); 
    
    // Measurement Model
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

cat("=== Identification Fixes Applied ===\n")
cat("1. alpha_zero[1] = 0으로 고정 (p_vec과의 식별 문제 해결)\n")
cat("2. eta의 각 열 평균을 0으로 제약 (transformed parameters에서)\n")
cat("3. p_vec의 합을 0에 가깝도록 soft constraint 추가\n")
cat("4. eta_raw를 사용하여 평균 제약 구현\n\n")
