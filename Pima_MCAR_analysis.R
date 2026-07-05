# =============================================================================
# Pima Indians Diabetes Analysis Pipeline -- MCAR Missingness
# Imputation vs. Imputation-Based Stacking, across 10-50% missingness
# =============================================================================

rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(caret); library(pROC); library(mice); library(dplyr); library(tidyr)
  library(readr); library(mlbench); library(data.table); library(irr); library(rpart)
  library(missMethods)
})

options(warn = -1)

# ---------------------------------------------------------------------------
# 1) CONFIGURATION
# ---------------------------------------------------------------------------

output_dir <- file.path(getwd(), "Pima_MCAR_Results")
dir.create(file.path(output_dir, "post_analysis"), recursive = TRUE, showWarnings = FALSE)

train_fraction  <- 0.80
m               <- 5
main_seed       <- 123
n_random_splits <- 50
missing_props   <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50)

methods           <- c("pmm", "rf", "cart", "norm", "midastouch")
method_order      <- c("Complete", "PMM", "RF", "CART", "NORM", "MIDASTOUCH", "Stacking")
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
  if (length(unique(true_y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(true_y, pred_prob, quiet = TRUE)))
}

brier_score <- function(true_y, pred_prob) {
  if (length(unique(true_y)) < 2) return(NA_real_)
  mean((as.numeric(true_y == levels(true_y)[2]) - pred_prob)^2)
}

confusion_metrics <- function(true_y, pred_class) {
  if (length(unique(true_y)) < 2) return(list(Sensitivity = NA, Specificity = NA, Precision = NA, F1 = NA))
  cm <- caret::confusionMatrix(pred_class, true_y, positive = "pos")
  prec <- as.numeric(cm$byClass["Precision"]); sens <- as.numeric(cm$byClass["Sensitivity"]); spec <- as.numeric(cm$byClass["Specificity"])
  f1 <- if (is.na(prec) || is.na(sens) || (prec + sens) == 0) NA_real_ else 2 * prec * sens / (prec + sens)
  list(Sensitivity = sens, Specificity = spec, Precision = prec, F1 = f1)
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

# MCAR missingness via missMethods::delete_MCAR, applied cell-wise across all
# eight predictors at the target overall proportion.
induce_mcar_missingness <- function(data_x, prop, seed = 123) {
  if (prop <= 0) return(as.data.frame(data_x))
  set.seed(seed)
  out <- missMethods::delete_MCAR(as.data.frame(data_x), p = prop, cols_mis = names(data_x),
                                   n_mis_stochastic = FALSE, p_overall = TRUE)
  out <- as.data.frame(out); names(out) <- names(data_x); out
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
  fit <- if (model_type == "GLM") {
    caret::train(y ~ ., data = train_df, method = "glm", family = "binomial", trControl = ctrl, metric = "ROC")
  } else if (model_type == "KNN") {
    caret::train(y ~ ., data = train_df, method = "knn", trControl = ctrl, metric = "ROC",
                 preProcess = c("center", "scale"), tuneGrid = fixed_knn_grid)
  } else if (model_type == "DT") {
    caret::train(y ~ ., data = train_df, method = "rpart", trControl = ctrl, metric = "ROC", tuneGrid = get_dt_grid())
  } else stop("Unsupported model_type")

  cv_pred <- fit$pred
  for (nm in names(fit$bestTune)) cv_pred <- cv_pred[cv_pred[[nm]] == fit$bestTune[[nm]], ]
  cv_pred <- cv_pred %>% arrange(rowIndex) %>% select(rowIndex, obs, pos)
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

scenario_results <- data.frame()
missingness_variable_all <- data.frame()
missingness_overall_all  <- data.frame()

for (prop_idx in seq_along(missing_props)) {
  prop <- missing_props[prop_idx]
  cat("Missing proportion:", prop, "\n")

  for (split_id in seq_len(n_random_splits)) {
    cat("  Split", split_id, "of", n_random_splits, "\n")
    data_seed <- main_seed + prop_idx * 1000 + split_id * 10000

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
      train_x <- induce_mcar_missingness(train_x, prop, seed = data_seed + 2)
      test_x  <- induce_mcar_missingness(test_x,  prop, seed = data_seed + 3)
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

        scenario_results <- bind_rows(scenario_results, data.frame(
          Split = split_id, Missing_Proportion = prop, Imputation_Method = "Complete", Model = clf,
          Train_AUC = compute_auc(f$cv_pred$obs, f$cv_pred$pos), Test_AUC = compute_auc(test_y, f$test_prob),
          Train_Accuracy = mean(train_class == f$cv_pred$obs), Test_Accuracy = mean(test_class == test_y),
          Train_Brier = brier_score(f$cv_pred$obs, f$cv_pred$pos), Test_Brier = brier_score(test_y, f$test_prob),
          Train_Precision = train_cm$Precision, Test_Precision = test_cm$Precision,
          Train_Sensitivity = train_cm$Sensitivity, Test_Sensitivity = test_cm$Sensitivity,
          Train_Specificity = train_cm$Specificity, Test_Specificity = test_cm$Specificity,
          Train_F1 = train_cm$F1, Test_F1 = test_cm$F1, Threshold = thresh
        ))
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

        scenario_results <- bind_rows(scenario_results, data.frame(
          Split = split_id, Missing_Proportion = prop, Imputation_Method = toupper(meth), Model = clf,
          Train_AUC = compute_auc(pooled_cv$obs, pooled_cv$prob), Test_AUC = compute_auc(test_y, test_prob),
          Train_Accuracy = mean(train_class == pooled_cv$obs), Test_Accuracy = mean(test_class == test_y),
          Train_Brier = brier_score(pooled_cv$obs, pooled_cv$prob), Test_Brier = brier_score(test_y, test_prob),
          Train_Precision = train_cm$Precision, Test_Precision = test_cm$Precision,
          Train_Sensitivity = train_cm$Sensitivity, Test_Sensitivity = test_cm$Sensitivity,
          Train_Specificity = train_cm$Specificity, Test_Specificity = test_cm$Specificity,
          Train_F1 = train_cm$F1, Test_F1 = test_cm$F1, Threshold = thresh
        ))
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

      scenario_results <- bind_rows(scenario_results, data.frame(
        Split = split_id, Missing_Proportion = prop, Imputation_Method = "Stacking", Model = clf,
        Train_AUC = compute_auc(pooled_stack$obs, pooled_stack$prob), Test_AUC = compute_auc(test_y, test_prob),
        Train_Accuracy = mean(train_class == pooled_stack$obs), Test_Accuracy = mean(test_class == test_y),
        Train_Brier = brier_score(pooled_stack$obs, pooled_stack$prob), Test_Brier = brier_score(test_y, test_prob),
        Train_Precision = train_cm$Precision, Test_Precision = test_cm$Precision,
        Train_Sensitivity = train_cm$Sensitivity, Test_Sensitivity = test_cm$Sensitivity,
        Train_Specificity = train_cm$Specificity, Test_Specificity = test_cm$Specificity,
        Train_F1 = train_cm$F1, Test_F1 = test_cm$F1, Threshold = thresh
      ))
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
safe_write_csv(scenario_results, file.path(output_dir, "Pima_MCAR_raw_results.csv"))
safe_write_csv(missingness_variable_all, file.path(output_dir, "post_analysis", "missingness_variable_level.csv"))
safe_write_csv(missingness_overall_all,  file.path(output_dir, "post_analysis", "missingness_overall.csv"))

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
                    NORM = "lightblue", MIDASTOUCH = "brown")
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

for (model_name in c("GLM", "KNN", "DT")) {
  dt_model <- rank_results[Model == model_name]
  combined <- rbind(fill = TRUE, do.call(rbind, lapply(names(metric_specs), function(mm) {
    extract_metric_table(dt_model, mm, metric_specs[[mm]])
  })))
  safe_write_csv(combined, file.path(output_dir, "post_analysis", paste0(model_name, "_combined_rank_table.csv")))
}

cat("\nDone. Outputs saved to:", output_dir, "\n")
