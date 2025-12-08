# =========================================================================
# Monte Carlo Simulation: BDI Index for Growth Mixture Models
# BDI = 1 - Cohen's Kappa (Boundary Discrepancy Index)
# =========================================================================

rm(list=ls())
library(MplusAutomation)
library(e1071)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)

# =========================================================================
# Simulation Settings
# =========================================================================
N_REP <- 100        # 반복 횟수
N_SUB <- 200        # 샘플 사이즈
K_FIT <- 2          # Mplus 적합 클래스 수
AMBIGUOUS_PCT <- 0.20  # 모호한 사례 비율 (하위 20%)
MIN_BOUNDARY_N <- 10   # 최소 경계 사례 수

# 작업 경로 설정 (본인 환경에 맞게)
setwd("C:/Users/daylee/Desktop/SVM index in GMM/Replication")

# 결과 저장용 데이터프레임
results_df <- data.frame(
  Rep_ID = integer(),
  Accuracy = numeric(),
  BDI = numeric(),
  SV_Ratio = numeric(),
  Entropy = numeric(),
  Converged = logical(),
  Warning = logical(),
  Boundary_N = integer(),      # 추가: 경계 사례 수
  Class_Balance = numeric(),    # 추가: 경계 사례의 클래스 균형도
  stringsAsFactors = FALSE
)

# Mplus 결과 파일 저장 폴더
sim_dir <- "Sim_Output_K2_n=200_unconstrained"
if(!dir.exists(sim_dir)) dir.create(sim_dir)
setwd(sim_dir) 

# =========================================================================
# Helper Functions
# =========================================================================

# Cohen's Kappa 계산 함수
calculate_kappa <- function(predicted, actual) {
  # 예측값과 실제값을 factor로 변환 (같은 levels 보장)
  all_levels <- unique(c(levels(factor(predicted)), levels(factor(actual))))
  predicted <- factor(predicted, levels = all_levels)
  actual <- factor(actual, levels = all_levels)
  
  # 혼동 행렬
  cm <- table(Predicted = predicted, Actual = actual)
  
  n_obs <- sum(cm)
  if(n_obs == 0) return(NA)
  
  # 관찰된 일치 비율
  p_o <- sum(diag(cm)) / n_obs
  
  # 기대 일치 비율 (chance agreement)
  p_e <- sum(rowSums(cm) * colSums(cm)) / (n_obs^2)
  
  # Kappa 계산
  if((1 - p_e) == 0 || is.na(p_e)) {
    # 완벽한 기대 일치 (모든 사례가 한 클래스) → kappa = 0
    return(0)
  } else {
    kappa <- (p_o - p_e) / (1 - p_e)
    return(kappa)
  }
}

# 클래스 균형도 계산 (0.5 = 완벽한 균형, 1.0 = 완전 불균형)
calculate_class_balance <- function(class_vector) {
  if(length(class_vector) == 0) return(NA)
  class_counts <- table(class_vector)
  max_prop <- max(class_counts) / sum(class_counts)
  return(max_prop)
}

# =========================================================================
# Start Simulation Loop
# =========================================================================
cat("Starting Simulation with", N_REP, "replications...\n")
cat("Sample size:", N_SUB, ", Classes:", K_FIT, ", Ambiguous %:", AMBIGUOUS_PCT, "\n\n")

for (r in 1:N_REP) {
  
  if(r %% 10 == 0) cat("Processing Replication:", r, "...\n")
  
  # -----------------------------------------------------------------------
  # Step 1: Data Generation (수정: N_SUB에 맞게 생성)
  # -----------------------------------------------------------------------
  set.seed(1234 + r) 
  
  Time_Points <- 0:4
  
  # 클래스 비율 설정 (예: 60/40)
  prop1 <- 0.6
  N1 <- round(N_SUB * prop1)
  N2 <- N_SUB - N1
  
  # Class 1
  Mu1 <- c(2.0, 0.5)
  Sigma1 <- matrix(c(0.2, 0.05, 0.05, 0.1), 2, 2)
  RE1 <- mvrnorm(N1, Mu1, Sigma1)
  
  # Class 2
  Mu2 <- c(9.0, 1.5)
  Sigma2 <- matrix(c(0.5, 0.1, 0.1, 0.2), 2, 2)
  RE2 <- mvrnorm(N2, Mu2, Sigma2)
  
  Alpha <- c(RE1[,1], RE2[,1])
  Beta <- c(RE1[,2], RE2[,2])
  
  # 종단 데이터 생성
  df_list <- list()
  for(i in 1:N_SUB) {
    y_values <- Alpha[i] + Beta[i] * Time_Points + rnorm(5, 0, 0.5)
    df_list[[i]] <- data.frame(ID = i, Time = Time_Points, Y = y_values)
  }
  df_long <- do.call(rbind, df_list)
  
  # Wide Format 변환
  df_wide <- df_long %>%
    mutate(Time_Label = paste0("y", Time + 1)) %>%
    pivot_wider(id_cols = "ID", names_from = Time_Label, values_from = Y) %>%
    dplyr::select(ID, y1, y2, y3, y4, y5)
  
  # -----------------------------------------------------------------------
  # Step 2: Run Mplus
  # -----------------------------------------------------------------------
  class_spec <- paste0("%c#", 1:K_FIT, "%\n [i s];\n i s;\n i WITH s;\n", collapse = "\n")
  
  mplus_model <- mplusObject(
    TITLE = paste0("Rep ", r),
    VARIABLE = paste0("NAMES = ID y1-y5;\n USEVARIABLES = y1-y5;\n CLASSES = c(", K_FIT, ");"),
    ANALYSIS = "TYPE = MIXTURE; STARTS = 800 400;", 
    MODEL = paste0("%OVERALL%\n i s | y1@0 y2@1 y3@2 y4@3 y5@4;\n\n", class_spec),
    OUTPUT = "TECH1;", 
    SAVEDATA = "FILE IS saved_prob.dat; SAVE = CPROB;",
    rdata = df_wide,
    usevariables = c("y1", "y2", "y3", "y4", "y5")
  )
  
  # Mplus 실행
  res <- mplusModeler(mplus_model, modelout = paste0("model_rep", r, ".inp"), run = 1L, quiet = TRUE)
  
  # -----------------------------------------------------------------------
  # Step 3: Calculate BDI & Kappa (Robust Exception Handling)
  # -----------------------------------------------------------------------
  
  # [예외처리 1] Mplus 실행 실패
  if(is.null(res) || is.null(res$results)) {
    results_df[r, ] <- list(r, NA, NA, NA, NA, FALSE, TRUE, NA, NA)
    next
  }
  
  # [예외처리 2] Savedata 없음
  bc_data <- res$results$savedata
  if(is.null(bc_data) || !is.data.frame(bc_data) || nrow(bc_data) == 0) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error=function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, FALSE, TRUE, NA, NA)
    next
  }
  
  # [예외처리 3] CPROB 컬럼 없음
  prob_cols <- grep("CPROB", names(bc_data), value = TRUE)
  if(length(prob_cols) == 0) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error=function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, FALSE, TRUE, NA, NA)
    next
  }
  
  # 최대 사후 확률 계산
  bc_data$Max_Prob <- apply(bc_data[, prob_cols, drop=FALSE], 1, max)
  
  # 모호한 사례 선택 (하위 AMBIGUOUS_PCT%)
  boundary_cutoff <- quantile(bc_data$Max_Prob, probs = AMBIGUOUS_PCT, na.rm = TRUE)
  boundary_indices <- which(bc_data$Max_Prob <= boundary_cutoff)
  boundary_data <- bc_data[boundary_indices, ]
  
  # [예외처리 4] 경계 사례 수 부족
  if(nrow(boundary_data) < MIN_BOUNDARY_N) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error=function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, TRUE, TRUE, nrow(boundary_data), NA)
    next
  }
  
  # 클래스 할당을 factor로 변환
  boundary_data$C <- as.factor(boundary_data$C)
  
  # 클래스 균형도 계산
  class_balance <- calculate_class_balance(boundary_data$C)
  
  # [예외처리 5] 단일 클래스만 있는 경우
  if(length(unique(boundary_data$C)) < 2) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error=function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, TRUE, TRUE, nrow(boundary_data), class_balance)
    next
  }
  
  # Y 변수 이름 확인 (Mplus가 대문자로 저장할 수 있음)
  y_cols <- grep("^[Yy][1-5]$", names(boundary_data), value = TRUE)
  if(length(y_cols) < 5) {
    # 대소문자 문제일 수 있음 - 직접 확인
    y_cols <- c("Y1", "Y2", "Y3", "Y4", "Y5")
    if(!all(y_cols %in% names(boundary_data))) {
      y_cols <- c("y1", "y2", "y3", "y4", "y5")
    }
  }
  
  if(length(y_cols) < 5 || !all(y_cols %in% names(boundary_data))) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error=function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, TRUE, TRUE, nrow(boundary_data), class_balance)
    next
  }
  
  # SVM 학습
  svm_formula <- as.formula(paste("C ~", paste(y_cols, collapse = " + ")))
  
  svm_fit <- tryCatch({
    svm(svm_formula, 
        data = boundary_data, 
        kernel = "radial", 
        cost = 0.5, 
        cross = 5)
  }, error = function(e) {
    warning(paste("SVM fitting failed in rep", r, ":", e$message))
    NULL
  })
  
  if(is.null(svm_fit)) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error=function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, TRUE, TRUE, nrow(boundary_data), class_balance)
    next
  }
  
  # 지표 계산
  svm_acc <- svm_fit$tot.accuracy / 100
  pred_class <- predict(svm_fit, boundary_data)
  
  # Cohen's Kappa 계산
  kappa <- calculate_kappa(pred_class, boundary_data$C)
  
  # BDI 계산 (1 - kappa)
  BDI <- 1 - kappa
  
  # Support Vector Ratio 계산
  # 주의: tot.nSV가 올바른 속성인지 확인 필요
  n_sv <- ifelse(is.null(svm_fit$tot.nSV), 
                 length(svm_fit$index),  # 대안: index 길이 사용
                 svm_fit$tot.nSV)
  sv_ratio <- n_sv / nrow(boundary_data)
  
  # 경고 메시지 여부 확인
  has_warning <- length(res$results$errors) > 0 || length(res$results$warnings) > 0
  
  # Entropy 추출
  entropy_val <- tryCatch(res$results$summaries$Entropy, error=function(e) NA)
  
  # 결과 저장
  results_df[r, ] <- list(
    r, 
    svm_acc, 
    BDI, 
    sv_ratio, 
    entropy_val,
    TRUE, 
    has_warning,
    nrow(boundary_data),
    class_balance
  )
}

# =========================================================================
# Step 4: Analyze Results
# =========================================================================
setwd("..") # 상위 폴더로 복귀

cat("\nSimulation Completed!\n")
cat("Total replications:", nrow(results_df), "\n")

total_reps <- nrow(results_df)

# 전체 요약 통계
summary_stats <- results_df %>%
  filter(Converged == TRUE & Warning == FALSE) %>% 
  summarise(
    # [1. 수렴 및 유효 해 비율]
    N_Valid = n(),
    Valid_Rate = n() / total_reps,
    
    # [2. 성과 지표 평균 (유효한 해들만의 평균)]
    Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    SD_Accuracy = sd(Accuracy, na.rm = TRUE),
    
    # [3. BDI 지표]
    Mean_BDI = mean(BDI, na.rm = TRUE),
    SD_BDI = sd(BDI, na.rm = TRUE),
    Median_BDI = median(BDI, na.rm = TRUE),
    Min_BDI = min(BDI, na.rm = TRUE),
    Max_BDI = max(BDI, na.rm = TRUE),
    
    # [4. SV Ratio 지표]
    Mean_SV_Ratio = mean(SV_Ratio, na.rm = TRUE),
    SD_SV_Ratio = sd(SV_Ratio, na.rm = TRUE),
    Median_SV_Ratio = median(SV_Ratio, na.rm = TRUE),
    
    # [5. Entropy]
    Mean_Entropy = mean(Entropy, na.rm = TRUE),
    SD_Entropy = sd(Entropy, na.rm = TRUE),
    
    # [6. 추가 진단 정보]
    Mean_Boundary_N = mean(Boundary_N, na.rm = TRUE),
    Mean_Class_Balance = mean(Class_Balance, na.rm = TRUE)
  )

print(summary_stats)

# 경고가 있는 경우와 없는 경우 비교
if(sum(results_df$Warning, na.rm = TRUE) > 0) {
  cat("\n=== Comparison: With vs Without Warnings ===\n")
  comparison <- results_df %>%
    filter(Converged == TRUE) %>%
    group_by(Warning) %>%
    summarise(
      N = n(),
      Mean_BDI = mean(BDI, na.rm = TRUE),
      Mean_SV_Ratio = mean(SV_Ratio, na.rm = TRUE),
      Mean_Entropy = mean(Entropy, na.rm = TRUE)
    )
  print(comparison)
}

# 결과 저장
write.csv(results_df, "simulation_results.csv", row.names = FALSE)
write.csv(summary_stats, "summary_statistics.csv", row.names = FALSE)

cat("\nResults saved to simulation_results.csv and summary_statistics.csv\n")
