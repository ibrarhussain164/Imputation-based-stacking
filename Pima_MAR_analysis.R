# =============================================================================
# Pima Indians Diabetes Analysis Pipeline -- MAR Missingness
# Imputation vs. Imputation-Based Stacking, across 10-50% missingness
# =============================================================================

rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(caret); library(pROC); library(mice); library(dplyr); library(tidyr)
  library(readr); library(mlbench); library(data.table); library(irr); library(rpart)
  library(missForest)
})

options(warn = -1)
select <- dplyr::select   # prevent MASS::select conflict

# ---------------------------------------------------------------------------
# 1) CONFIGURATION
# ---------------------------------------------------------------------------

output_dir <- file.path(getwd(), "Pima_MAR_Results")
dir.create(file.path(output_dir, "post_analysis"), recursive = TRUE, showWarnings = FALSE)

train_fraction  <- 0.80
m               <- 5
main_seed       <- 123
n_random_splits <- 50
missing_props   <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50)

methods           <- c("pmm", "rf", "cart", "norm", "midastouch")
method_order      <- c("Complete", "CCA", "Mean", "MissInd", "MissForest",
                       "PMM", "RF", "CART", "NORM", "MIDASTOUCH",
                       "MIE_SE_GLM", "MIE_SE_KNN", "MIE_SE_DT", "Stacking")
classifier_models <- c("GLM", "KNN", "DT")

cv_folds_standard <- 5
cv_folds_stacking <- 5
cv_folds_meta     <- 5

fixed_knn_grid <- data.frame(k = seq(3, 17, by = 2))
get_dt_grid    <- function() expand.grid(cp = c(0.0005, 0.001, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.1))

safe_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(as.data.frame(df), path)
}

# ---------------------------------------------------------------------------
# 2) DATA PREPARATION
# ---------------------------------------------------------------------------

prepare_pima_data <- function() {
  data("PimaIndiansDiabetes2", package = "mlbench")
  dat <- na.omit(PimaIndiansDiabetes2)
  dat$y <- factor(ifelse(dat$diabetes == "pos", "pos", "neg"), levels = c("neg", "pos"))
  dat$diabetes <- NULL
  dat <- dat[, c("y", "pregnant", "glucose", "pressure", "triceps", "insulin", "mass", "pedigree", "age")]
  names(dat) <- c("y", paste0("x", seq_len(ncol(dat) - 1)))
  dat
}

dat <- prepare_pima_data()
cat(sprintf("Pima: %d rows, %d predictors, prevalence = %.3f\n", nrow(dat), ncol(dat) - 1, mean(dat$y == "pos")))

# ---------------------------------------------------------------------------
# 3) METRIC FUNCTIONS
# ---------------------------------------------------------------------------

compute_auc <- function(true_y, pred_prob) {
  if (length(unique(na.omit(droplevels(as.factor(true_y))))) < 2) return(NA_real_)
  tryCatch(
    as.numeric(pROC::auc(pROC::roc(true_y, pred_prob, quiet = TRUE))),
    error = function(e) NA_real_
  )
}

brier_score <- function(true_y, pred_prob) {
  if (length(unique(na.omit(droplevels(as.factor(true_y))))) < 2) return(NA_real_)
  mean((as.numeric(true_y == levels(as.factor(true_y))[2]) - pred_prob)^2, na.rm = TRUE)
}

confusion_metrics <- function(true_y, pred_class) {
  if (length(unique(na.omit(droplevels(as.factor(true_y))))) < 2) return(list(Sensitivity = NA, Specificity = NA, Precision = NA, F1 = NA))
  cm <- caret::confusionMatrix(pred_class, true_y, positive = "pos")
  prec <- as.numeric(cm$byClass["Precision"]); sens <- as.numeric(cm$byClass["Sensitivity"]); spec <- as.numeric(cm$byClass["Specificity"])
  f1 <- if (is.na(prec) || is.na(sens) || (prec + sens) == 0) NA_real_ else 2 * prec * sens / (prec + sens)
  list(Sensitivity = sens, Specificity = spec, Precision = prec, F1 = f1)
}

compute_calibration <- function(obs_y, pred_prob) {
  if (is.factor(obs_y)) obs_y_num <- as.numeric(obs_y == "pos") else obs_y_num <- as.numeric(obs_y)
  pred_prob <- pmin(pmax(as.numeric(pred_prob), 1e-6), 1 - 1e-6)
  if (length(unique(na.omit(obs_y_num))) < 2 || sum(!is.na(pred_prob)) < 5)
    return(list(cal_slope = NA_real_, cal_intercept = NA_real_))
  logit_pred <- log(pred_prob / (1 - pred_prob))
  cal_fit <- tryCatch(glm(obs_y_num ~ logit_pred, family = binomial()), error = function(e) NULL)
  list(
    cal_slope     = if (!is.null(cal_fit)) as.numeric(coef(cal_fit)[2]) else NA_real_,
    cal_intercept = if (!is.null(cal_fit)) as.numeric(coef(cal_fit)[1]) else NA_real_
  )
}

find_optimal_threshold <- function(obs, pred_prob) {
  obs <- factor(obs, levels = c("neg", "pos"))
  grid <- seq(0.01, 0.99, by = 0.01)
  res <- bind_rows(lapply(grid, function(th) {
    pc <- factor(ifelse(pred_prob >= th, "pos", "neg"), levels = c("neg", "pos"))
    cm <- confusion_metrics(obs, pc)
    data.frame(threshold = th, Youden = cm$Sensitivity + cm$Specificity - 1, F1 = cm$F1)
  }))
  if (nrow(res) == 0) return(0.50)
  res <- res[order(-res$Youden, -res$F1, res$threshold), ]
  as.numeric(res$threshold[1])
}

# ---------------------------------------------------------------------------
# 4) MISSINGNESS AND IMPUTATION HELPERS
# ---------------------------------------------------------------------------

# MAR missingness via mice::ampute; predictor x1 (pregnant) and x8 (age) act
# as fully-observed anchors and are protected from missingness themselves.
induce_mar_ampute <- function(data_x, prop, seed = 123) {
  if (prop <= 0) return(as.data.frame(data_x))
  stopifnot(ncol(data_x) == 8)
  set.seed(seed)
  pattern <- matrix(c(1, 0, 0, 0, 0, 0, 0, 1), nrow = 1)
  weights <- matrix(c(1, 0, 0, 0, 0, 0, 0, 1), nrow = 1)
  amp <- mice::ampute(as.data.frame(data_x), prop = prop, patterns = pattern,
                      weights = weights, bycases = FALSE, mech = "MAR")
  out <- as.data.frame(amp$amp); names(out) <- names(data_x); out
}

winsorize_limits <- function(train_x_complete) {
  lapply(train_x_complete, function(x) {
    st <- boxplot.stats(x)$stats
    list(low = st[1], high = st[5])
  })
}

apply_winsorization <- function(data_x, limits) {
  for (col in names(limits)) {
    if (col %in% names(data_x)) {
      lim <- limits[[col]]
      data_x[[col]] <- pmin(pmax(data_x[[col]], lim$low), lim$high)
    }
  }
  data_x
}

summarize_missingness <- function(df_x, dataset_label, split_id, target_prop) {
  var_tbl <- data.frame(
    Split = split_id, Dataset = dataset_label, Missing_Proportion = target_prop, Variable = names(df_x),
    Missing_Count = sapply(df_x, function(x) sum(is.na(x))), Total_Count = nrow(df_x),
    Missing_Percent = sapply(df_x, function(x) mean(is.na(x)))
  )
  overall_tbl <- data.frame(
    Split = split_id, Dataset = dataset_label, Missing_Proportion = target_prop,
    Total_Cells = nrow(df_x) * ncol(df_x), Missing_Cells = sum(is.na(as.matrix(df_x))),
    Overall_Missing_Percent = mean(is.na(as.matrix(df_x)))
  )
  list(variable_level = var_tbl, overall_level = overall_tbl)
}

impute_train_test <- function(train_x, test_x, method_name, m = 5, seed = 123, maxit = 10) {
  n_train <- nrow(train_x); n_test <- nrow(test_x)
  set.seed(seed)
  mids_train <- mice::mice(train_x, method = method_name, m = m, maxit = maxit, printFlag = FALSE)
  
  train_list <- vector("list", m); test_list <- vector("list", m)
  for (k in seq_len(m)) {
    train_completed <- mice::complete(mids_train, action = k)
    combined <- bind_rows(train_completed, test_x)
    set.seed(seed + k)
    mids_test <- mice::mice.mids(mids_train, newdata = combined, maxit = 1, printFlag = FALSE)
    completed_full <- mice::complete(mids_test, action = k)
    completed_full[1:n_train, ] <- train_completed
    train_list[[k]] <- completed_full[1:n_train, , drop = FALSE]
    test_list[[k]]  <- completed_full[(n_train + 1):(n_train + n_test), , drop = FALSE]
  }
  list(train = train_list, test = test_list)
}

# ---------------------------------------------------------------------------
# 5) MODELING HELPERS
# ---------------------------------------------------------------------------

fit_model_cv <- function(train_df, test_df, model_type, seed, cv_folds = 5) {
  ctrl <- caret::trainControl(method = "cv", number = cv_folds, classProbs = TRUE,
                              summaryFunction = caret::twoClassSummary, savePredictions = "final")
  set.seed(seed)
  fit <- tryCatch({
    if (model_type == "GLM") {
      caret::train(y ~ ., data = train_df, method = "glm", family = "binomial", trControl = ctrl, metric = "ROC")
    } else if (model_type == "KNN") {
      caret::train(y ~ ., data = train_df, method = "knn", trControl = ctrl, metric = "ROC",
                   preProcess = c("center", "scale"), tuneGrid = fixed_knn_grid)
    } else if (model_type == "DT") {
      caret::train(y ~ ., data = train_df, method = "rpart", trControl = ctrl, metric = "ROC", tuneGrid = get_dt_grid())
    } else stop("Unsupported model_type")
  }, error = function(e) {
    message("  [fit_model_cv] caret training failed (", model_type, "): ", e$message)
    NULL
  })
  
  if (is.null(fit)) {
    n_train <- nrow(train_df)
    dummy_cv <- data.frame(rowIndex = seq_len(n_train),
                           obs      = train_df$y,
                           pos      = rep(mean(train_df$y == "pos"), n_train))
    dummy_test <- rep(mean(train_df$y == "pos"), nrow(test_df))
    return(list(cv_pred = dummy_cv, test_prob = dummy_test))
  }
  
  cv_pred <- fit$pred
  for (nm in names(fit$bestTune)) cv_pred <- cv_pred[cv_pred[[nm]] == fit$bestTune[[nm]], ]
  cv_pred <- cv_pred %>% arrange(rowIndex) %>% dplyr::select(rowIndex, obs, pos)
  test_prob <- predict(fit, newdata = test_df, type = "prob")$pos
  list(cv_pred = cv_pred, test_prob = test_prob)
}

fit_meta_learner <- function(meta_train_df, meta_test_df, seed, cv_folds = 5) {
  meta_train_df$y <- factor(meta_train_df$y, levels = c("neg", "pos"))
  base_rate <- mean(meta_train_df$y == "pos")
  
  set.seed(seed)
  folds <- caret::createFolds(meta_train_df$y, k = cv_folds, returnTrain = FALSE)
  oof <- rep(NA_real_, nrow(meta_train_df))
  
  for (f in folds) {
    tr <- meta_train_df[-f, ]; va <- meta_train_df[f, ]
    fit <- tryCatch(glm(y ~ ., data = tr, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) { oof[f] <- mean(tr$y == "pos"); next }
    pred <- tryCatch(as.numeric(predict(fit, newdata = va, type = "response")),
                     error = function(e) rep(mean(tr$y == "pos"), nrow(va)))
    pred[!is.finite(pred)] <- mean(tr$y == "pos")
    oof[f] <- pred
  }
  oof[!is.finite(oof)] <- base_rate
  
  final_fit <- tryCatch(glm(y ~ ., data = meta_train_df, family = binomial()), error = function(e) NULL)
  test_pred <- if (is.null(final_fit)) rep(base_rate, nrow(meta_test_df)) else
    tryCatch(as.numeric(predict(final_fit, newdata = meta_test_df, type = "response")),
             error = function(e) rep(base_rate, nrow(meta_test_df)))
  test_pred[!is.finite(test_pred)] <- base_rate
  list(train_oof_pred = oof, test_pred = test_pred)
}

# ---------------------------------------------------------------------------
# 6) MAIN EXPERIMENT LOOP
# ---------------------------------------------------------------------------

missingness_variable_all <- data.frame()
missingness_overall_all  <- data.frame()
timing_results <- data.frame()

# ---------------------------------------------------------------------------
# CHECKPOINT: load previous progress if it exists
# ---------------------------------------------------------------------------
checkpoint_file <- file.path(output_dir, "checkpoint_raw_results.csv")

if (file.exists(checkpoint_file)) {
  cat("=== CHECKPOINT FOUND — resuming from saved progress ===\n")
  scenario_results <- readr::read_csv(checkpoint_file, show_col_types = FALSE)
  scenario_results$Imputation_Method <- as.character(scenario_results$Imputation_Method)
  cat("Checkpoint: loaded", nrow(scenario_results), "rows covering",
      length(unique(paste(scenario_results$Missing_Proportion, scenario_results$Split))),
      "completed splits\n\n")
  completed_keys <- unique(paste(scenario_results$Missing_Proportion,
                                 scenario_results$Split, sep = "_"))
} else {
  scenario_results <- data.frame()
  completed_keys   <- character(0)
  cat("No checkpoint found — starting fresh.\n\n")
}

for (prop_idx in seq_along(missing_props)) {
  prop <- missing_props[prop_idx]
  cat("Missing proportion:", prop, "\n")
  
  for (split_id in seq_len(n_random_splits)) {
    
    # ---- Skip if already completed ----------------------------------------
    split_key <- paste(prop, split_id, sep = "_")
    if (split_key %in% completed_keys) {
      cat("  Split", split_id, "- skipping (checkpoint)\n")
      next
    }
    
    cat("  Split", split_id, "of", n_random_splits, "\n")
    data_seed <- main_seed + (round(prop * 10) + 1) * 1000 + split_id * 10000
    
    n_before_split   <- nrow(scenario_results)
    split_start_time <- proc.time()
    
    set.seed(data_seed + 1)
    idx <- caret::createDataPartition(dat$y, p = train_fraction, list = FALSE)
    train <- dat[idx, ]; test <- dat[-idx, ]
    train_y <- factor(train$y, levels = c("neg", "pos")); test_y <- factor(test$y, levels = c("neg", "pos"))
    train_x_complete <- train[, -1, drop = FALSE]; test_x_complete <- test[, -1, drop = FALSE]
    
    # Winsorize on COMPLETE data first, then induce missingness
    limits <- winsorize_limits(train_x_complete)
    train_x <- apply_winsorization(train_x_complete, limits)
    test_x  <- apply_winsorization(test_x_complete, limits)
    
    if (prop > 0) {
      train_x <- induce_mar_ampute(train_x, prop, seed = data_seed + 2)
      test_x  <- induce_mar_ampute(test_x,  prop, seed = data_seed + 3)
    }
    
    train_miss <- summarize_missingness(train_x, "Train", split_id, prop)
    test_miss  <- summarize_missingness(test_x, "Test", split_id, prop)
    missingness_variable_all <- bind_rows(missingness_variable_all, train_miss$variable_level, test_miss$variable_level)
    missingness_overall_all  <- bind_rows(missingness_overall_all, train_miss$overall_level, test_miss$overall_level)
    
    # ---- Complete-case baseline (prop == 0) ------------------------------
    if (prop == 0) {
      complete_train <- data.frame(train_x, y = train_y)
      for (clf in classifier_models) {
        f <- fit_model_cv(complete_train, test_x, clf, seed = data_seed + 10 + match(clf, classifier_models))
        thresh <- find_optimal_threshold(f$cv_pred$obs, f$cv_pred$pos)
        train_class <- factor(ifelse(f$cv_pred$pos >= thresh, "pos", "neg"), levels = c("neg", "pos"))
        train_cm <- confusion_metrics(f$cv_pred$obs, train_class)
        test_class <- factor(ifelse(f$test_prob >= thresh, "pos", "neg"), levels = c("neg", "pos"))
        test_cm <- confusion_metrics(test_y, test_class)
        
        cal_complete <- compute_calibration(test_y, f$test_prob)
        scenario_results <- bind_rows(scenario_results, data.frame(
          Split = split_id, Missing_Proportion = prop, Imputation_Method = "Complete", Model = clf,
          Train_AUC = compute_auc(f$cv_pred$obs, f$cv_pred$pos), Test_AUC = compute_auc(test_y, f$test_prob),
          Train_Accuracy = mean(train_class == f$cv_pred$obs), Test_Accuracy = mean(test_class == test_y),
          Train_Brier = brier_score(f$cv_pred$obs, f$cv_pred$pos), Test_Brier = brier_score(test_y, f$test_prob),
          Train_Precision = train_cm$Precision, Test_Precision = test_cm$Precision,
          Train_Sensitivity = train_cm$Sensitivity, Test_Sensitivity = test_cm$Sensitivity,
          Train_Specificity = train_cm$Specificity, Test_Specificity = test_cm$Specificity,
          Train_F1 = train_cm$F1, Test_F1 = test_cm$F1, Threshold = thresh,
          Cal_Slope = cal_complete$cal_slope, Cal_Intercept = cal_complete$cal_intercept
        ))
      }
      # Record timing for complete-case split
      split_elapsed <- (proc.time() - split_start_time)[["elapsed"]]
      timing_results <- bind_rows(timing_results, data.frame(
        Split = split_id, Missing_Proportion = prop, Elapsed_Seconds = split_elapsed
      ))
      
      # Save checkpoint
      new_rows <- if (nrow(scenario_results) > n_before_split)
        scenario_results[(n_before_split + 1):nrow(scenario_results), ] else data.frame()
      if (nrow(new_rows) > 0) {
        readr::write_csv(new_rows, checkpoint_file,
                         append = file.exists(checkpoint_file))
        completed_keys <- c(completed_keys, split_key)
        cat("  Checkpoint saved (prop=", prop, ", split=", split_id, ")\n", sep = "")
      }
      next
    }
    
    # ---- Multiple imputation methods -------------------------------------
    imputed <- lapply(methods, function(meth) {
      impute_train_test(train_x, test_x, meth, m = m, seed = data_seed + 100 + match(meth, methods))
    })
    names(imputed) <- methods
    
    train_imp <- list(); test_imp <- list()
    for (meth in methods) {
      train_imp[[meth]] <- lapply(imputed[[meth]]$train, function(d) data.frame(d, y = train_y))
      test_imp[[meth]]  <- imputed[[meth]]$test
    }
    
    for (meth in methods) {
      for (clf in classifier_models) {
        cv_list <- list(); test_list <- list()
        for (k in seq_len(m)) {
          f <- fit_model_cv(train_imp[[meth]][[k]], test_imp[[meth]][[k]], clf,
                            seed = data_seed + 1000 + match(meth, methods) * 100 + match(clf, classifier_models) * 10 + k,
                            cv_folds = cv_folds_standard)
          cv_list[[k]] <- f$cv_pred %>% rename(!!paste0("p", k) := pos)
          test_list[[k]] <- f$test_prob
        }
        pooled_cv <- Reduce(function(x, y) full_join(x, y, by = c("rowIndex", "obs")), cv_list) %>% arrange(rowIndex)
        pooled_cv$prob <- rowMeans(as.matrix(pooled_cv[, grep("^p", names(pooled_cv))]), na.rm = TRUE)
        thresh <- find_optimal_threshold(pooled_cv$obs, pooled_cv$prob)
        train_class <- factor(ifelse(pooled_cv$prob >= thresh, "pos", "neg"), levels = c("neg", "pos"))
        train_cm <- confusion_metrics(pooled_cv$obs, train_class)
        
        test_prob <- rowMeans(do.call(cbind, test_list), na.rm = TRUE)
        test_class <- factor(ifelse(test_prob >= thresh, "pos", "neg"), levels = c("neg", "pos"))
        test_cm <- confusion_metrics(test_y, test_class)
        cal_mice <- compute_calibration(test_y, test_prob)
        
        scenario_results <- bind_rows(scenario_results, data.frame(
          Split = split_id, Missing_Proportion = prop, Imputation_Method = toupper(meth), Model = clf,
          Train_AUC = compute_auc(pooled_cv$obs, pooled_cv$prob), Test_AUC = compute_auc(test_y, test_prob),
          Train_Accuracy = mean(train_class == pooled_cv$obs), Test_Accuracy = mean(test_class == test_y),
          Train_Brier = brier_score(pooled_cv$obs, pooled_cv$prob), Test_Brier = brier_score(test_y, test_prob),
          Train_Precision = train_cm$Precision, Test_Precision = test_cm$Precision,
          Train_Sensitivity = train_cm$Sensitivity, Test_Sensitivity = test_cm$Sensitivity,
          Train_Specificity = train_cm$Specificity, Test_Specificity = test_cm$Specificity,
          Train_F1 = train_cm$F1, Test_F1 = test_cm$F1, Threshold = thresh,
          Cal_Slope = cal_mice$cal_slope, Cal_Intercept = cal_mice$cal_intercept
        ))
      }
    }
    
    # ---- Simple baselines: CCA, Mean, MissInd, MissForest (R1-C6 / R2-C1) ---
    
    store_single_result <- function(method_name, clf, cv_p, test_prob_vec) {
      thresh  <- find_optimal_threshold(cv_p$obs, cv_p$pos)
      tr_cl   <- factor(ifelse(cv_p$pos      >= thresh, "pos", "neg"), levels = c("neg","pos"))
      te_cl   <- factor(ifelse(test_prob_vec >= thresh, "pos", "neg"), levels = c("neg","pos"))
      tr_cm   <- confusion_metrics(cv_p$obs, tr_cl)
      te_cm   <- confusion_metrics(test_y,   te_cl)
      cal     <- compute_calibration(test_y, test_prob_vec)
      data.frame(
        Split = split_id, Missing_Proportion = prop,
        Imputation_Method = method_name, Model = clf,
        Train_AUC      = compute_auc(cv_p$obs, cv_p$pos),
        Test_AUC       = compute_auc(test_y, test_prob_vec),
        Train_Accuracy = mean(tr_cl == cv_p$obs),
        Test_Accuracy  = mean(te_cl == test_y),
        Train_Brier    = brier_score(cv_p$obs, cv_p$pos),
        Test_Brier     = brier_score(test_y, test_prob_vec),
        Train_Precision = tr_cm$Precision,  Test_Precision  = te_cm$Precision,
        Train_Sensitivity=tr_cm$Sensitivity, Test_Sensitivity=te_cm$Sensitivity,
        Train_Specificity=tr_cm$Specificity, Test_Specificity=te_cm$Specificity,
        Train_F1 = tr_cm$F1, Test_F1 = te_cm$F1, Threshold = thresh,
        Cal_Slope = cal$cal_slope, Cal_Intercept = cal$cal_intercept
      )
    }
    
    # 1) CCA
    cat("    Method: CCA\n")
    cca_complete <- complete.cases(train_x)
    if (sum(cca_complete) >= 10) {
      train_x_cca  <- train_x[cca_complete, , drop = FALSE]
      train_y_cca  <- train_y[cca_complete]
      train_means_cca <- colMeans(train_x, na.rm = TRUE)
      test_x_cca   <- as.data.frame(lapply(names(test_x), function(col) {
        v <- test_x[[col]]; v[is.na(v)] <- train_means_cca[[col]]; v
      })); names(test_x_cca) <- names(test_x)
      train_df_cca <- data.frame(train_x_cca, y = train_y_cca)
      for (clf in classifier_models) {
        f <- fit_model_cv(train_df_cca, test_x_cca, clf,
                          seed = data_seed + 30000 + match(clf, classifier_models) * 10,
                          cv_folds = cv_folds_standard)
        scenario_results <- bind_rows(scenario_results,
                                      store_single_result("CCA", clf, f$cv_pred, f$test_prob))
      }
    }
    
    # 2) Mean imputation
    cat("    Method: Mean Imputation\n")
    train_means_m <- colMeans(train_x, na.rm = TRUE)
    train_x_mean  <- train_x; test_x_mean <- test_x
    for (col in names(train_x_mean)) {
      train_x_mean[[col]][is.na(train_x_mean[[col]])] <- train_means_m[[col]]
      test_x_mean[[col]][is.na(test_x_mean[[col]])]   <- train_means_m[[col]]
    }
    train_df_mean <- data.frame(train_x_mean, y = train_y)
    for (clf in classifier_models) {
      f <- fit_model_cv(train_df_mean, test_x_mean, clf,
                        seed = data_seed + 31000 + match(clf, classifier_models) * 10,
                        cv_folds = cv_folds_standard)
      scenario_results <- bind_rows(scenario_results,
                                    store_single_result("Mean", clf, f$cv_pred, f$test_prob))
    }
    
    # 3) Missing indicator ─ binary flag per missing predictor + mean-fill (Sterne et al. 2009)
    cat("    Method: Missing Indicator\n")
    has_miss_cols <- names(train_x)[sapply(train_x, function(x) any(is.na(x)))]
    train_x_ind <- train_x; test_x_ind <- test_x
    for (col in has_miss_cols) {
      train_x_ind[[paste0(col, "_miss")]] <- as.numeric(is.na(train_x[[col]]))
      test_x_ind[[paste0(col,  "_miss")]] <- as.numeric(is.na(test_x[[col]]))
      train_x_ind[[col]][is.na(train_x_ind[[col]])] <- train_means_m[[col]]
      test_x_ind[[col]][is.na(test_x_ind[[col]])]   <- train_means_m[[col]]
    }
    train_df_ind <- data.frame(train_x_ind, y = train_y)
    for (clf in classifier_models) {
      f <- fit_model_cv(train_df_ind, test_x_ind, clf,
                        seed = data_seed + 32000 + match(clf, classifier_models) * 10,
                        cv_folds = cv_folds_standard)
      scenario_results <- bind_rows(scenario_results,
                                    store_single_result("MissInd", clf, f$cv_pred, f$test_prob))
    }
    
    # 4) MissForest
    cat("    Method: MissForest\n")
    n_train_mf    <- nrow(train_x)
    combined_x_mf <- rbind(as.data.frame(train_x), as.data.frame(test_x))
    set.seed(data_seed + 33000)
    mf_out <- tryCatch(
      missForest::missForest(combined_x_mf, verbose = FALSE),
      error = function(e) { message("MissForest failed: ", e$message); NULL }
    )
    if (!is.null(mf_out)) {
      cimp         <- as.data.frame(mf_out$ximp)
      train_x_mf  <- cimp[seq_len(n_train_mf), , drop = FALSE]
      test_x_mf   <- cimp[(n_train_mf + 1):nrow(cimp), , drop = FALSE]
      train_df_mf <- data.frame(train_x_mf, y = train_y)
      for (clf in classifier_models) {
        f <- fit_model_cv(train_df_mf, test_x_mf, clf,
                          seed = data_seed + 33000 + match(clf, classifier_models) * 10,
                          cv_folds = cv_folds_standard)
        scenario_results <- bind_rows(scenario_results,
                                      store_single_result("MissForest", clf, f$cv_pred, f$test_prob))
      }
    }
    
    # ---- Imputation-based stacking ----------------------------------------
    for (clf in classifier_models) {
      stack_train_list <- list(); stack_test_mat <- matrix(NA_real_, nrow(test_x), m)
      
      for (k in seq_len(m)) {
        oof_list <- list(); test_pred_list <- list()
        for (meth in methods) {
          f <- fit_model_cv(train_imp[[meth]][[k]], test_imp[[meth]][[k]], clf,
                            seed = data_seed + 5000 + match(meth, methods) * 100 + match(clf, classifier_models) * 1000 + k,
                            cv_folds = cv_folds_stacking)
          oof_list[[meth]] <- f$cv_pred %>% transmute(row_id = rowIndex, !!paste0(meth, "_pred") := pos)
          test_pred_list[[meth]] <- f$test_prob
        }
        combined_oof <- Reduce(function(x, y) full_join(x, y, by = "row_id"), oof_list) %>% arrange(row_id)
        meta_cols <- grep("_pred$", names(combined_oof), value = TRUE)
        meta_train <- data.frame(row_id = combined_oof$row_id, y = train_y[combined_oof$row_id], combined_oof[, meta_cols])
        meta_test  <- as.data.frame(do.call(cbind, test_pred_list)); colnames(meta_test) <- paste0(methods, "_pred")
        
        meta_fit <- fit_meta_learner(meta_train[, c("y", meta_cols)], meta_test[, meta_cols],
                                     seed = data_seed + 9000 + match(clf, classifier_models) * 100 + k, cv_folds = cv_folds_meta)
        stack_train_list[[k]] <- data.frame(row_id = meta_train$row_id, obs = meta_train$y, pred = meta_fit$train_oof_pred)
        stack_test_mat[, k] <- meta_fit$test_pred
      }
      
      pooled_stack <- Reduce(function(x, y) full_join(x, y, by = c("row_id", "obs")),
                             lapply(seq_along(stack_train_list), function(i) {
                               out <- stack_train_list[[i]]; names(out)[3] <- paste0("p", i); out
                             })) %>% arrange(row_id)
      pooled_stack$prob <- rowMeans(as.matrix(pooled_stack[, grep("^p", names(pooled_stack))]), na.rm = TRUE)
      thresh <- find_optimal_threshold(pooled_stack$obs, pooled_stack$prob)
      train_class <- factor(ifelse(pooled_stack$prob >= thresh, "pos", "neg"), levels = c("neg", "pos"))
      train_cm <- confusion_metrics(pooled_stack$obs, train_class)
      
      test_prob <- rowMeans(stack_test_mat, na.rm = TRUE)
      test_class <- factor(ifelse(test_prob >= thresh, "pos", "neg"), levels = c("neg", "pos"))
      test_cm <- confusion_metrics(test_y, test_class)
      cal_stack <- compute_calibration(test_y, test_prob)
      
      scenario_results <- bind_rows(scenario_results, data.frame(
        Split = split_id, Missing_Proportion = prop, Imputation_Method = "Stacking", Model = clf,
        Train_AUC = compute_auc(pooled_stack$obs, pooled_stack$prob), Test_AUC = compute_auc(test_y, test_prob),
        Train_Accuracy = mean(train_class == pooled_stack$obs), Test_Accuracy = mean(test_class == test_y),
        Train_Brier = brier_score(pooled_stack$obs, pooled_stack$prob), Test_Brier = brier_score(test_y, test_prob),
        Train_Precision = train_cm$Precision, Test_Precision = test_cm$Precision,
        Train_Sensitivity = train_cm$Sensitivity, Test_Sensitivity = test_cm$Sensitivity,
        Train_Specificity = train_cm$Specificity, Test_Specificity = test_cm$Specificity,
        Train_F1 = train_cm$F1, Test_F1 = test_cm$F1, Threshold = thresh,
        Cal_Slope = cal_stack$cal_slope, Cal_Intercept = cal_stack$cal_intercept
      ))
    }
    
    # ---- MIE_SE: Multiple Imputation Ensemble Stacking (Aleryani et al.) -------
    # MICE_SE variant: PMM (m=5 draws) x GLM+KNN+DT = 15 OOF predictions as meta-features.
    # Meta-learner runs once per classifier (GLM/KNN/DT) → 3 rows: MIE_SE_GLM, MIE_SE_KNN, MIE_SE_DT.
    # Base model seeds: 20000 range. Meta-learner seeds: data_seed + 25000 range.
    cat("    Method: MIE_SE (PMM x GLM+KNN+DT, per-clf meta-learner)\n")
    
    mie_se_oof_list  <- list()
    mie_se_test_list <- list()
    mie_se_idx <- 1
    
    for (k in seq_len(m)) {
      for (clf in classifier_models) {
        f_mie <- fit_model_cv(
          train_imp[["pmm"]][[k]],
          test_imp[["pmm"]][[k]],
          clf,
          seed     = data_seed + 20000 + match(clf, classifier_models) * 10 + k,
          cv_folds = cv_folds_standard
        )
        col_name <- paste0("pmm_k", k, "_", clf, "_pred")
        mie_se_oof_list[[mie_se_idx]]  <- f_mie$cv_pred %>%
          transmute(row_id = rowIndex, !!col_name := pos)
        mie_se_test_list[[mie_se_idx]] <- f_mie$test_prob
        mie_se_idx <- mie_se_idx + 1
      }
    }
    
    mie_se_oof <- Reduce(function(x, y) full_join(x, y, by = "row_id"), mie_se_oof_list) %>%
      arrange(row_id)
    mie_se_meta_cols <- grep("_pred$", names(mie_se_oof), value = TRUE)
    
    meta_train_mie <- data.frame(
      y = train_y[mie_se_oof$row_id],
      mie_se_oof[, mie_se_meta_cols]
    )
    mie_se_test_df <- as.data.frame(do.call(cbind, mie_se_test_list))
    colnames(mie_se_test_df) <- mie_se_meta_cols
    
    # One meta-learner per classifier
    for (clf_meta in classifier_models) {
      mie_label <- paste0("MIE_SE_", clf_meta)
      f_meta <- fit_model_cv(
        meta_train_mie,
        mie_se_test_df,
        clf_meta,
        seed     = data_seed + 25000 + match(clf_meta, classifier_models) * 100,
        cv_folds = cv_folds_meta
      )
      meta_oof <- f_meta$cv_pred %>% arrange(rowIndex)
      mie_se_train_prob <- meta_oof$pos
      mie_se_train_obs  <- meta_oof$obs
      mie_se_train_prob[!is.finite(mie_se_train_prob)] <- mean(train_y == "pos")
      thresh_mie_se      <- find_optimal_threshold(mie_se_train_obs, mie_se_train_prob)
      mie_se_train_class <- factor(ifelse(mie_se_train_prob >= thresh_mie_se, "pos", "neg"), levels = c("neg", "pos"))
      mie_se_train_cm    <- confusion_metrics(mie_se_train_obs, mie_se_train_class)
      
      mie_se_test_prob  <- f_meta$test_prob
      mie_se_test_prob[!is.finite(mie_se_test_prob)] <- mean(train_y == "pos")
      mie_se_test_class <- factor(ifelse(mie_se_test_prob >= thresh_mie_se, "pos", "neg"), levels = c("neg", "pos"))
      mie_se_test_cm    <- confusion_metrics(test_y, mie_se_test_class)
      cal_mie           <- compute_calibration(test_y, mie_se_test_prob)
      
      scenario_results <- bind_rows(scenario_results, data.frame(
        Split = split_id, Missing_Proportion = prop, Imputation_Method = mie_label, Model = mie_label,
        Train_AUC      = compute_auc(mie_se_train_obs, mie_se_train_prob),
        Test_AUC       = compute_auc(test_y, mie_se_test_prob),
        Train_Accuracy = mean(mie_se_train_class == mie_se_train_obs),
        Test_Accuracy  = mean(mie_se_test_class == test_y),
        Train_Brier    = brier_score(mie_se_train_obs, mie_se_train_prob),
        Test_Brier     = brier_score(test_y, mie_se_test_prob),
        Train_Precision    = mie_se_train_cm$Precision,
        Test_Precision     = mie_se_test_cm$Precision,
        Train_Sensitivity  = mie_se_train_cm$Sensitivity,
        Test_Sensitivity   = mie_se_test_cm$Sensitivity,
        Train_Specificity  = mie_se_train_cm$Specificity,
        Test_Specificity   = mie_se_test_cm$Specificity,
        Train_F1   = mie_se_train_cm$F1,
        Test_F1    = mie_se_test_cm$F1,
        Threshold  = thresh_mie_se,
        Cal_Slope  = cal_mie$cal_slope, Cal_Intercept = cal_mie$cal_intercept
      ))
    }
    
    # ---- Record per-split computational cost --------------------------------
    split_elapsed <- (proc.time() - split_start_time)[["elapsed"]]
    timing_results <- bind_rows(timing_results, data.frame(
      Split = split_id, Missing_Proportion = prop, Elapsed_Seconds = split_elapsed
    ))
    
    # ---- Save checkpoint ----------------------------------------------------
    new_rows <- if (nrow(scenario_results) > n_before_split)
      scenario_results[(n_before_split + 1):nrow(scenario_results), ] else data.frame()
    if (nrow(new_rows) > 0) {
      readr::write_csv(new_rows, checkpoint_file,
                       append = file.exists(checkpoint_file))
      completed_keys <- c(completed_keys, split_key)
      cat("  Checkpoint saved (prop=", prop, ", split=", split_id,
          ", elapsed=", round(split_elapsed, 1), "s)\n", sep = "")
    }
  }
}

scenario_results <- scenario_results %>%
  mutate(
    Overfit_Gap_AUC         = Train_AUC - Test_AUC,
    Overfit_Gap_Accuracy    = Train_Accuracy - Test_Accuracy,
    Overfit_Gap_Brier       = Test_Brier - Train_Brier,
    Overfit_Gap_Precision   = Train_Precision - Test_Precision,
    Overfit_Gap_Sensitivity = Train_Sensitivity - Test_Sensitivity,
    Overfit_Gap_Specificity = Train_Specificity - Test_Specificity,
    Overfit_Gap_F1          = Train_F1 - Test_F1
  )

scenario_results$Imputation_Method <- factor(scenario_results$Imputation_Method, levels = method_order)
safe_write_csv(scenario_results, file.path(output_dir, "Pima_MAR_raw_results.csv"))
safe_write_csv(missingness_variable_all, file.path(output_dir, "post_analysis", "missingness_variable_level.csv"))
safe_write_csv(missingness_overall_all,  file.path(output_dir, "post_analysis", "missingness_overall.csv"))

# ---- Save per-split timing results ----------------------------------------
timing_summary <- timing_results %>%
  group_by(Missing_Proportion) %>%
  summarise(
    Mean_Elapsed_Seconds = mean(Elapsed_Seconds, na.rm = TRUE),
    SD_Elapsed_Seconds   = sd(Elapsed_Seconds, na.rm = TRUE),
    Min_Elapsed_Seconds  = min(Elapsed_Seconds, na.rm = TRUE),
    Max_Elapsed_Seconds  = max(Elapsed_Seconds, na.rm = TRUE),
    n_splits             = n(),
    .groups = "drop"
  )
safe_write_csv(timing_results,  file.path(output_dir, "post_analysis", "per_split_timing.csv"))
safe_write_csv(timing_summary,  file.path(output_dir, "post_analysis", "timing_summary.csv"))

missingness_summary <- missingness_overall_all %>%
  filter(Missing_Proportion > 0) %>%
  group_by(Missing_Proportion, Dataset) %>%
  summarise(Mean_Missing = mean(Overall_Missing_Percent), SD_Missing = sd(Overall_Missing_Percent), n_splits = n(), .groups = "drop")
safe_write_csv(missingness_summary, file.path(output_dir, "post_analysis", "missingness_mean_sd_summary.csv"))

overfitting_summary <- scenario_results %>%
  filter(!is.na(Missing_Proportion)) %>%
  group_by(Model, Missing_Proportion, Method = Imputation_Method) %>%
  summarise(across(starts_with("Overfit_Gap_"), list(Mean = ~mean(.x, na.rm = TRUE), SD = ~sd(.x, na.rm = TRUE))),
            n_splits = n(), .groups = "drop")
safe_write_csv(overfitting_summary, file.path(output_dir, "post_analysis", "overfitting_mean_sd_summary.csv"))

# ---------------------------------------------------------------------------
# 7) SUMMARY: MEAN +/- SD, 95% CI, KENDALL'S W RANKING
# ---------------------------------------------------------------------------

final_results <- scenario_results %>%
  filter(!is.na(Missing_Proportion)) %>%
  group_by(Model, Missing_Proportion, Method = Imputation_Method) %>%
  summarise(across(c(Test_AUC, Test_Accuracy, Test_Brier, Test_Precision,
                     Test_Sensitivity, Test_Specificity, Test_F1),
                   list(Mean = ~mean(.x, na.rm = TRUE), SD = ~sd(.x, na.rm = TRUE))),
            n_splits = n(), .groups = "drop")

safe_write_csv(final_results, file.path(output_dir, "post_analysis", "final_results_by_model.csv"))

# ---- 95% confidence intervals ---------------------------------------------
metric_info_ci <- data.frame(
  Metric = c("AUC", "Accuracy", "Brier", "Precision", "Sensitivity", "Specificity", "F1"),
  Column = c("Test_AUC", "Test_Accuracy", "Test_Brier", "Test_Precision",
             "Test_Sensitivity", "Test_Specificity", "Test_F1")
)

make_ci_table <- function(results, metric_name, metric_col) {
  results %>%
    filter(!is.na(Missing_Proportion)) %>%
    group_by(Model, Missing_Proportion, Method = Imputation_Method) %>%
    summarise(Metric = metric_name, n = sum(!is.na(.data[[metric_col]])),
              Mean = mean(.data[[metric_col]], na.rm = TRUE), SD = sd(.data[[metric_col]], na.rm = TRUE),
              SE = SD / sqrt(n),
              CI_Lower = Mean - qt(0.975, df = n - 1) * SE, CI_Upper = Mean + qt(0.975, df = n - 1) * SE,
              .groups = "drop")
}

ci_95_results <- bind_rows(lapply(seq_len(nrow(metric_info_ci)),
                                  function(i) make_ci_table(scenario_results, metric_info_ci$Metric[i], metric_info_ci$Column[i])))
safe_write_csv(ci_95_results, file.path(output_dir, "post_analysis", "test_metrics_95CI.csv"))

# ---- Figures: mean +/- SD vs. missingness, one 7-panel figure per classifier
suppressPackageStartupMessages({ library(ggplot2); library(ggpubr); library(scales) })

method_colors <- c(Stacking = "blue", PMM = "grey", RF = "lightgreen", CART = "lightpink",
                   NORM = "lightblue", MIDASTOUCH = "brown",
                   MIE_SE_GLM = "mediumpurple", MIE_SE_KNN = "purple", MIE_SE_DT = "darkviolet",
                   CCA = "darkred", Mean = "orange", MissInd = "darkorange",
                   MissForest = "darkgreen")
breaks_seq <- sort(unique(final_results$Missing_Proportion[final_results$Missing_Proportion > 0]))

plot_metric_model <- function(data, mean_col, sd_col, y_label, title_letter, model_name) {
  pd <- data %>% filter(Method != "Complete", Missing_Proportion > 0)
  ggplot(pd, aes(x = Missing_Proportion, y = .data[[mean_col]], color = Method, group = Method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = .data[[mean_col]] - .data[[sd_col]], ymax = .data[[mean_col]] + .data[[sd_col]]), width = 0.01) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = breaks_seq, labels = percent_format(accuracy = 1)) +
    theme_minimal() +
    labs(title = paste0(title_letter, ". ", y_label, " (", model_name, ")"),
         x = "Missing Percentage", y = paste0("Average ", y_label, " (\u00b1 SD)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold", size = 10),
          legend.position = "bottom")
}

make_classifier_figure <- function(model_name, data_model, output_dir) {
  specs <- list(
    list(mean = "Test_AUC_Mean", sd = "Test_AUC_SD", label = "AUC", letter = "A"),
    list(mean = "Test_Accuracy_Mean", sd = "Test_Accuracy_SD", label = "Accuracy", letter = "B"),
    list(mean = "Test_Brier_Mean", sd = "Test_Brier_SD", label = "Brier", letter = "C"),
    list(mean = "Test_Precision_Mean", sd = "Test_Precision_SD", label = "Precision", letter = "D"),
    list(mean = "Test_Sensitivity_Mean", sd = "Test_Sensitivity_SD", label = "Sensitivity", letter = "E"),
    list(mean = "Test_Specificity_Mean", sd = "Test_Specificity_SD", label = "Specificity", letter = "F"),
    list(mean = "Test_F1_Mean", sd = "Test_F1_SD", label = "F1", letter = "G")
  )
  plots <- lapply(specs, function(s) plot_metric_model(data_model, s$mean, s$sd, s$label, s$letter, model_name))
  fig <- ggarrange(plotlist = plots, ncol = 2, nrow = 4, common.legend = TRUE, legend = "bottom")
  ggsave(file.path(output_dir, "post_analysis", paste0("Fig_", model_name, "_ALL.pdf")),
         plot = fig, width = 14, height = 16, device = cairo_pdf)
  fig
}

for (model_name in c("GLM", "KNN", "DT")) {
  make_classifier_figure(model_name, final_results %>% filter(Model == model_name), output_dir)
}

create_rank_table <- function(dt, metric, higher_better = TRUE) {
  dt_copy <- copy(dt); dt_copy[[metric]] <- as.numeric(dt_copy[[metric]])
  wide <- dcast(dt_copy, Method ~ Missing_Proportion, value.var = metric)
  prop_cols <- setdiff(names(wide), "Method")
  wide[, Mean_Value := rowMeans(.SD, na.rm = TRUE), .SDcols = prop_cols]
  wide[, Overall_Rank := if (higher_better) rank(-Mean_Value, ties.method = "min") else rank(Mean_Value, ties.method = "min")]
  rank_cols <- paste0("Rank_", prop_cols)
  if (higher_better) wide[, (rank_cols) := lapply(.SD, function(x) rank(-x, ties.method = "average")), .SDcols = prop_cols]
  else wide[, (rank_cols) := lapply(.SD, function(x) rank(x, ties.method = "average")), .SDcols = prop_cols]
  kendall_res <- suppressWarnings(irr::kendall(as.matrix(wide[, ..rank_cols])))
  list(rank_table = wide[, c("Method", rank_cols, "Overall_Rank"), with = FALSE], kendall = kendall_res)
}

rank_results <- final_results %>%
  filter(Method != "Complete") %>%
  rename(AUC = Test_AUC_Mean, Accuracy = Test_Accuracy_Mean, Brier = Test_Brier_Mean,
         Precision = Test_Precision_Mean, Sensitivity = Test_Sensitivity_Mean,
         Specificity = Test_Specificity_Mean, F1 = Test_F1_Mean)
setDT(rank_results)

extract_metric_table <- function(dt, metric_name, higher_better) {
  res <- create_rank_table(dt, metric_name, higher_better)
  out <- as.data.frame(res$rank_table)
  out$Metric <- metric_name
  out$W <- round(as.numeric(res$kendall$value), 3)
  out$p_value <- signif(as.numeric(res$kendall$p.value), 4)
  out
}

metric_specs <- list(AUC = TRUE, Accuracy = TRUE, Brier = FALSE, Precision = TRUE,
                     Sensitivity = TRUE, Specificity = TRUE, F1 = TRUE)

make_paper_rank_table <- function(combined_df) {
  # Rename Rank_0.1 -> 10%, Rank_0.2 -> 20%, etc.
  rank_col_pattern <- grep("^Rank_", names(combined_df), value = TRUE)
  new_names <- sapply(rank_col_pattern, function(col) {
    pct <- as.numeric(sub("Rank_", "", col)) * 100
    paste0(round(pct), "%")
  }, USE.NAMES = FALSE)
  names(combined_df)[match(rank_col_pattern, names(combined_df))] <- new_names
  
  # Sort by Metric then Overall_Rank
  combined_df <- combined_df[order(combined_df$Metric, combined_df$Overall_Rank), ]
  
  # Move Metric to first column
  other_cols <- setdiff(names(combined_df), "Metric")
  combined_df <- combined_df[, c("Metric", other_cols)]
  
  # Show W and p_value only on first row of each metric group
  combined_df$p_value <- ifelse(combined_df$p_value < 0.001, "<0.001",
                                as.character(combined_df$p_value))
  first_row <- !duplicated(combined_df$Metric)
  combined_df$W[!first_row] <- NA
  combined_df$p_value[!first_row] <- NA
  
  combined_df
}

for (model_name in c("GLM", "KNN", "DT")) {
  dt_model <- rank_results[Model == model_name]
  combined <- rbind(fill = TRUE, do.call(rbind, lapply(names(metric_specs), function(mm) {
    extract_metric_table(dt_model, mm, metric_specs[[mm]])
  })))
  safe_write_csv(combined, file.path(output_dir, "post_analysis", paste0(model_name, "_combined_rank_table.csv")))
  paper_tbl <- make_paper_rank_table(as.data.frame(combined))
  safe_write_csv(paper_tbl, file.path(output_dir, "post_analysis", paste0(model_name, "_paper_rank_table.csv")))
}

# ---------------------------------------------------------------------------
# PAIRED DIFFERENCES + WIN PROPORTIONS (R1-C5)
# ---------------------------------------------------------------------------

individual_methods <- c("PMM", "RF", "CART", "NORM", "MIDASTOUCH",
                        "CCA", "Mean", "MissInd", "MissForest")

paired_diff_results <- data.frame()

for (model_name in classifier_models) {
  for (prop_val in missing_props[missing_props > 0]) {
    sub <- scenario_results %>%
      filter(Model == model_name, Missing_Proportion == prop_val,
             Imputation_Method %in% c("Stacking", individual_methods))
    if (nrow(sub) == 0) next
    wide <- sub %>%
      select(Split, Imputation_Method, Test_AUC) %>%
      pivot_wider(names_from = Imputation_Method, values_from = Test_AUC)
    if (!"Stacking" %in% names(wide)) next
    indiv_cols <- intersect(individual_methods, names(wide))
    if (length(indiv_cols) == 0) next
    best_indiv <- apply(wide[, indiv_cols, drop = FALSE], 1, max, na.rm = TRUE)
    diff_vec   <- wide$Stacking - best_indiv
    n          <- sum(!is.na(diff_vec))
    if (n < 2) next
    mean_diff <- mean(diff_vec, na.rm = TRUE)
    sd_diff   <- sd(diff_vec, na.rm = TRUE)
    se_diff   <- sd_diff / sqrt(n)
    paired_diff_results <- bind_rows(paired_diff_results, data.frame(
      Model = model_name, Missing_Proportion = prop_val, Metric = "AUC",
      n_splits = n, Mean_Diff = mean_diff, SD_Diff = sd_diff,
      CI_Lower = mean_diff - qt(0.975, df = n-1) * se_diff,
      CI_Upper = mean_diff + qt(0.975, df = n-1) * se_diff,
      Win_Proportion = mean(diff_vec > 0, na.rm = TRUE)
    ))
  }
}

safe_write_csv(paired_diff_results,
               file.path(output_dir, "post_analysis", "stacking_paired_diff_win_prop.csv"))

cat("\nDone. Outputs saved to:", output_dir, "\n")
