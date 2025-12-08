rm(list = ls())

library(MplusAutomation)
library(e1071)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)

# -----------------------------------------------------------------------------
# User Settings
# -----------------------------------------------------------------------------
setwd("C:/Users/daylee/Desktop/SVM index in GMM/Replication")

N_REP <- 100                  # Number of Monte Carlo replications
N_SUB <- 200                  # Sample size per replication
K_FIT <- 2                    # Number of classes to fit in the GMM
CLASS_WEIGHTS <- c(0.6, 0.4)  # True class prevalence used at simulation

AMBIGUITY_CEILING <- 0.65     # Posterior max-probability threshold for ambiguity
MIN_BOUNDARY_PER_CLASS <- 8   # Require at least this many per class in boundary set
TRAIN_PROP <- 0.7             # Split for SVM train/test inside boundary pool

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
compute_kappa <- function(pred, actual) {
  if (length(pred) == 0 || length(actual) == 0) return(NA_real_)
  all_levels <- union(levels(factor(pred)), levels(factor(actual)))
  pred <- factor(pred, levels = all_levels)
  actual <- factor(actual, levels = all_levels)
  cm <- table(Predicted = pred, Actual = actual)
  n_obs <- sum(cm)
  if (n_obs == 0) return(NA_real_)
  p_o <- sum(diag(cm)) / n_obs
  row_tot <- rowSums(cm)
  col_tot <- colSums(cm)
  p_e <- sum(row_tot * col_tot) / (n_obs^2)
  if ((1 - p_e) == 0) {
    kappa <- 0
  } else {
    kappa <- (p_o - p_e) / (1 - p_e)
  }
  kappa <- max(min(kappa, 1), -1)
  return(kappa)
}

# -----------------------------------------------------------------------------
# Containers
# -----------------------------------------------------------------------------
results_df <- data.frame(
  Rep_ID = integer(),
  Accuracy = numeric(),
  BDI = numeric(),
  SV_Ratio = numeric(),
  Entropy = numeric(),
  Kappa_Mplus_Truth = numeric(),
  Kappa_SVM_Truth = numeric(),
  Converged = logical(),
  Warning = logical(),
  stringsAsFactors = FALSE
)

sim_dir <- "Sim_Output_K2_n=200_unconstrained"
if (!dir.exists(sim_dir)) dir.create(sim_dir)
setwd(sim_dir)

cat("Starting Simulation with", N_REP, "replications...\n")

for (r in seq_len(N_REP)) {
  if (r %% 10 == 0) cat("Processing Replication:", r, "...\n")

  set.seed(1234 + r)
  Time_Points <- 0:4

  true_classes <- sample(seq_len(K_FIT), N_SUB, replace = TRUE, prob = CLASS_WEIGHTS)

  mu_list <- list(
    c(2.0, 0.5),
    c(9.0, 1.5)
  )
  sigma_list <- list(
    matrix(c(0.2, 0.05, 0.05, 0.1), 2, 2),
    matrix(c(0.5, 0.1, 0.1, 0.2), 2, 2)
  )

  Alpha <- numeric(N_SUB)
  Beta <- numeric(N_SUB)
  for (k in seq_len(K_FIT)) {
    idx <- which(true_classes == k)
    if (length(idx) > 0) {
      rand_eff <- mvrnorm(length(idx), mu = mu_list[[k]], Sigma = sigma_list[[k]])
      Alpha[idx] <- rand_eff[, 1]
      Beta[idx] <- rand_eff[, 2]
    }
  }

  df_list <- vector("list", N_SUB)
  for (i in seq_len(N_SUB)) {
    y_values <- Alpha[i] + Beta[i] * Time_Points + rnorm(length(Time_Points), 0, 0.5)
    df_list[[i]] <- data.frame(ID = i, Time = Time_Points, Y = y_values)
  }
  df_long <- bind_rows(df_list)

  df_wide <- df_long %>%
    mutate(Time_Label = paste0("y", Time + 1)) %>%
    pivot_wider(id_cols = "ID", names_from = Time_Label, values_from = Y) %>%
    arrange(ID)

  df_wide$True_Class <- factor(true_classes)
  df_wide <- df_wide %>% select(ID, True_Class, y1, y2, y3, y4, y5)

  class_spec <- paste0("%c#", 1:K_FIT, "%\n  [i s];\n  i s;\n  i WITH s;\n", collapse = "\n")
  variable_section <- paste(
    "NAMES = ID True_Class y1-y5;",
    "USEVARIABLES = y1-y5;",
    paste0("CLASSES = c(", K_FIT, ");"),
    "IDVARIABLE = ID;",
    "AUXILIARY = True_Class;",
    sep = "\n"
  )

  mplus_model <- mplusObject(
    TITLE = paste0("Rep ", r),
    VARIABLE = variable_section,
    ANALYSIS = "TYPE = MIXTURE; STARTS = 800 400;",
    MODEL = paste0("%OVERALL%\n  i s | y1@0 y2@1 y3@2 y4@3 y5@4;\n\n", class_spec),
    OUTPUT = "TECH1;",
    SAVEDATA = "FILE IS saved_prob.dat; SAVE = CPROB;",
    rdata = df_wide,
    usevariables = c("y1", "y2", "y3", "y4", "y5", "True_Class")
  )

  res <- tryCatch(
    mplusModeler(mplus_model, modelout = paste0("model_rep", r, ".inp"), run = 1L, quiet = TRUE),
    error = function(e) NULL
  )

  if (is.null(res) || is.null(res$results)) {
    results_df[r, ] <- list(r, NA, NA, NA, NA, NA, NA, FALSE, TRUE)
    next
  }

  bc_data <- res$results$savedata
  if (is.null(bc_data) || !is.data.frame(bc_data) || nrow(bc_data) == 0) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error = function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, NA, NA, FALSE, TRUE)
    next
  }

  prob_cols <- grep("CPROB", names(bc_data), value = TRUE)
  if (length(prob_cols) == 0) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error = function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, NA, NA, FALSE, TRUE)
    next
  }

  bc_data <- bc_data %>%
    mutate(Max_Prob = apply(select(., all_of(prob_cols)), 1, max)) %>%
    left_join(df_wide %>% select(ID, True_Class), by = "ID")

  boundary_data <- bc_data %>%
    filter(Max_Prob <= AMBIGUITY_CEILING)

  class_counts <- table(boundary_data$C)
  if (
    nrow(boundary_data) < (2 * MIN_BOUNDARY_PER_CLASS) ||
      any(class_counts < MIN_BOUNDARY_PER_CLASS)
  ) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error = function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, NA, NA, TRUE, TRUE)
    next
  }

  boundary_data <- boundary_data %>%
    mutate(
      C = factor(C),
      True_Class = factor(True_Class)
    )

  set.seed(9000 + r)
  train_idx <- sample(seq_len(nrow(boundary_data)), size = max(1, floor(TRAIN_PROP * nrow(boundary_data))))
  if (length(train_idx) == nrow(boundary_data)) train_idx <- head(train_idx, -1)
  if (length(train_idx) == 0) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error = function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, NA, NA, TRUE, TRUE)
    next
  }

  boundary_train <- boundary_data[train_idx, ]
  boundary_test <- boundary_data[-train_idx, ]

  svm_fit <- tryCatch(
    svm(C ~ Y1 + Y2 + Y3 + Y4 + Y5,
        data = boundary_train,
        kernel = "radial",
        cost = 0.5),
    error = function(e) NULL
  )

  if (is.null(svm_fit) || nrow(boundary_test) == 0) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error = function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, NA, NA, TRUE, TRUE)
    next
  }

  pred_class <- predict(svm_fit, boundary_test)
  pred_class <- factor(pred_class, levels = levels(boundary_data$C))

  cm <- table(Predicted = pred_class, Actual = boundary_test$C)
  if (sum(cm) == 0) {
    entropy_val <- tryCatch(res$results$summaries$Entropy, error = function(e) NA)
    results_df[r, ] <- list(r, NA, NA, NA, entropy_val, NA, NA, TRUE, TRUE)
    next
  }

  kappa_boundary <- compute_kappa(pred_class, boundary_test$C)
  BDI <- (1 - kappa_boundary) / 2
  accuracy <- mean(pred_class == boundary_test$C)

  kappa_truth_mplus <- compute_kappa(boundary_test$C, boundary_test$True_Class)
  kappa_truth_svm <- compute_kappa(pred_class, boundary_test$True_Class)

  sv_ratio <- svm_fit$tot.nSV / nrow(boundary_data)

  has_warning <- length(res$results$errors) > 0 || length(res$results$warnings) > 0
  entropy_val <- tryCatch(res$results$summaries$Entropy, error = function(e) NA)

  results_df[r, ] <- list(
    r,
    accuracy,
    BDI,
    sv_ratio,
    entropy_val,
    kappa_truth_mplus,
    kappa_truth_svm,
    TRUE,
    has_warning
  )
}

setwd("..")

cat("\nSimulation Completed!\n")

total_reps <- nrow(results_df)
summary_stats <- results_df %>%
  filter(Converged == TRUE & Warning == FALSE) %>%
  summarise(
    N_Valid = n(),
    Valid_Rate = n() / total_reps,
    Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    Mean_BDI = mean(BDI, na.rm = TRUE),
    SD_BDI = sd(BDI, na.rm = TRUE),
    Mean_SV_Ratio = mean(SV_Ratio, na.rm = TRUE),
    SD_SV_Ratio = sd(SV_Ratio, na.rm = TRUE),
    Mean_Entropy = mean(Entropy, na.rm = TRUE),
    Mean_Kappa_Mplus_Truth = mean(Kappa_Mplus_Truth, na.rm = TRUE),
    Mean_Kappa_SVM_Truth = mean(Kappa_SVM_Truth, na.rm = TRUE)
  )

summary_stats
