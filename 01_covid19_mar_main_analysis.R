# =============================================================================
# COVID-19 MAR imputation stacking analysis
# Outcome: Death (0 = survivor, 1 = death)
# Purpose: compare individual imputation methods and imputation-based stacking
# under increasing MAR missingness across repeated train/test splits.
# =============================================================================

# =============================================================================
# 0) CLEAN SESSION AND LOAD PACKAGES
# =============================================================================

rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(readxl)
  library(caret)
  library(pROC)
  library(mice)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(rlang)
  library(flextable)
  library(officer)
  library(irr)
  library(rpart)
  library(knitr)
  library(kableExtra)
  library(cowplot)
  library(grid)
  library(magrittr)
  library(missForest)
})

options(warn = -1)
select <- dplyr::select   # prevent MASS::select conflict


# =============================================================================
# File helpers
# =============================================================================

safe_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(path, ".tmp")
  
  if (file.exists(tmp_path)) {
    unlink(tmp_path, force = TRUE)
  }
  
  readr::write_csv(as.data.frame(df), tmp_path)
  
  if (file.exists(path)) {
    unlink(path, force = TRUE)
  }
  
  ok <- file.rename(tmp_path, path)
  
  if (!ok) {
    stop(paste0(
      "Could not write file: ", path, "\n",
      "Close the file if it is open and run the script again."
    ))
  }
  
  invisible(path)
}

safe_write_text <- function(text, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, con = path, useBytes = TRUE)
  invisible(path)
}


# =============================================================================
# 1) USER CONFIGURATION
# =============================================================================

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_root <- Sys.getenv("COVID19_OUTPUT_ROOT", unset = file.path(getwd(), "outputs"))
output_dir <- Sys.getenv(
  "COVID19_OUTPUT_DIR",
  unset = file.path(output_root, paste0("COVID19_MAR_Stacking_", run_stamp))
)
output_dir <- path.expand(output_dir)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "scenario_raw_results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "scenario_test_summaries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "scenario_overfitting_summaries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "post_analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "post_analysis", "latex_tables"), recursive = TRUE, showWarnings = FALSE)

safe_write_text(output_dir, file.path(getwd(), "latest_covid19_output_dir.txt"))

cat("Working directory:", getwd(), "\n")
cat("Output directory:", output_dir, "\n")

# The data file can be supplied with the first command-line argument or with
# Sys.setenv(COVID19_DATA_FILE = "path/to/file.xlsx").
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[1])) {
  Sys.setenv(COVID19_DATA_FILE = args[1])
}

data_file_path <- Sys.getenv(
  "COVID19_DATA_FILE",
  unset = "~/Downloads/untitled folder/Mortality_incidence_sociodemographic_and_clinical_data_in_COVID19_patients 7.xlsx"
)

train_fraction   <- 0.80
m                <- 5
main_seed        <- 123
n_random_splits  <- 20

missing_props_default <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50)

missing_props_by_scenario <- list(
  COVID19 = missing_props_default
)

methods <- c("pmm", "rf", "cart", "norm", "midastouch")
method_order <- c("Complete", "CCA", "Mean", "MissInd", "MissForest",
                  "PMM", "RF", "CART", "NORM", "MIDASTOUCH",
                  "MIE_SE_GLM", "MIE_SE_KNN", "MIE_SE_DT", "Stacking")
classifier_models <- c("GLM", "KNN", "DT")

cv_folds_standard <- 5
cv_folds_stacking <- 5
cv_folds_meta     <- 5

fixed_knn_grid <- data.frame(k = seq(3, 17, by = 2))


# =============================================================================
# 2) DATA PREPARATION
# =============================================================================

read_covid_data <- function(file_path) {
  if (exists("COVID19", envir = .GlobalEnv)) {
    message("Using existing object: COVID19")
    return(as.data.frame(get("COVID19", envir = .GlobalEnv)))
  }
  
  file_path <- path.expand(file_path)
  
  if (!file.exists(file_path)) {
    stop(paste0(
      "COVID-19 data file not found at: ", file_path, "\n",
      "Place the Excel file at data/covid19_data.xlsx, pass the file path as an argument, ",
      "or set COVID19_DATA_FILE."
    ))
  }
  
  ext <- tolower(tools::file_ext(file_path))
  
  if (ext %in% c("xlsx", "xls")) {
    out <- readxl::read_excel(file_path)
  } else if (ext == "csv") {
    out <- readr::read_csv(
      file_path,
      na = c("", "NA", "NaN", ".", "?", " "),
      show_col_types = FALSE
    )
  } else {
    stop("Unsupported file type. Use .xlsx, .xls, or .csv.")
  }
  
  as.data.frame(out)
}

prepare_covid_data <- function(file_path) {
  covid_raw <- read_covid_data(file_path)
  
  required_source_cols <- c(
    "Death", "Age...27", "LOS", "Severity", "OsSats", "MAP", "Ddimer", "Plts",
    "BUN", "Creatinine", "Sodium", "Glucose", "AST", "ALT", "WBC", "IL6",
    "Ferritin", "CrctProtein", "Procalcitonin"
  )
  
  missing_source <- setdiff(required_source_cols, names(covid_raw))
  if (length(missing_source) > 0) {
    stop(paste0(
      "These required columns were not found in the input file: ",
      paste(missing_source, collapse = ", ")
    ))
  }
  
  dat <- covid_raw %>%
    dplyr::select(
      Death,
      Age = Age...27,
      LOS,
      Severity,
      OsSats,
      MAP,
      Ddimer,
      Plts,
      BUN,
      Creatinine,
      Sodium,
      Glucose,
      AST,
      ALT,
      WBC,
      IL6,
      Ferritin,
      CrctProtein,
      Procalcitonin
    )
  
  required_cols <- c(
    "Death", "Age", "LOS", "Severity", "OsSats", "MAP", "Ddimer", "Plts",
    "BUN", "Creatinine", "Sodium", "Glucose", "AST", "ALT", "WBC", "IL6",
    "Ferritin", "CrctProtein", "Procalcitonin"
  )
  
  for (col in required_cols) {
    dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  }
  
  dat <- dat %>%
    dplyr::filter(Death %in% c(0, 1)) %>%
    stats::na.omit()
  
  dat$y <- factor(ifelse(dat$Death == 1, "pos", "neg"), levels = c("neg", "pos"))
  dat$Death <- NULL
  
  dat <- dat[, c(
    "y", "Age", "LOS", "Severity", "OsSats", "MAP", "Ddimer", "Plts",
    "BUN", "Creatinine", "Sodium", "Glucose", "AST", "ALT", "WBC", "IL6",
    "Ferritin", "CrctProtein", "Procalcitonin"
  )]
  
  variable_map <- data.frame(
    Original_Name = names(dat)[-1],
    Model_Name = paste0("x", seq_len(ncol(dat) - 1)),
    stringsAsFactors = FALSE
  )
  
  names(dat) <- c("y", variable_map$Model_Name)
  attr(dat, "variable_map") <- variable_map
  
  cat("COVID-19 dataset loaded successfully.\n")
  cat("Rows after removing original missing values:", nrow(dat), "\n")
  cat("Outcome counts:\n")
  print(table(dat$y))
  cat("\n")
  
  dat
}

source_datasets <- list(
  COVID19 = prepare_covid_data(data_file_path)
)

covid_variable_map <- attr(source_datasets$COVID19, "variable_map")
safe_write_csv(covid_variable_map, file.path(output_dir, "covid19_variable_map.csv"))


cat("Prepared source datasets:\n")
for (nm in names(source_datasets)) {
  cat(sprintf(
    "- %s: %d rows, %d predictors, source prevalence = %.3f\n",
    nm,
    nrow(source_datasets[[nm]]),
    ncol(source_datasets[[nm]]) - 1,
    mean(source_datasets[[nm]]$y == "pos")
  ))
}
cat("\n")

# =============================================================================
# 3) METRIC FUNCTIONS
# =============================================================================

compute_auc <- function(true_y, pred_prob) {
  if (length(unique(na.omit(droplevels(as.factor(true_y))))) < 2) return(NA_real_)
  tryCatch({
    roc_obj <- pROC::roc(true_y, pred_prob, quiet = TRUE)
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) NA_real_)
}

brier_score <- function(true_y, pred_prob) {
  if (length(unique(na.omit(droplevels(as.factor(true_y))))) < 2) return(NA_real_)
  true_numeric <- as.numeric(true_y == levels(as.factor(true_y))[2])
  mean((true_numeric - pred_prob)^2, na.rm = TRUE)
}

precision_score <- function(true_y, pred_class) {
  if (length(unique(na.omit(droplevels(as.factor(true_y))))) < 2) return(NA_real_)
  tryCatch({
    cm <- caret::confusionMatrix(pred_class, true_y, positive = "pos")
    as.numeric(cm$byClass["Precision"])
  }, error = function(e) NA_real_)
}

confusion_metrics <- function(true_y, pred_class) {
  if (length(unique(na.omit(droplevels(as.factor(true_y))))) < 2) {
    return(list(
      Sensitivity = NA_real_,
      Specificity = NA_real_,
      F1          = NA_real_
    ))
  }
  
  cm <- caret::confusionMatrix(pred_class, true_y, positive = "pos")
  precision   <- as.numeric(cm$byClass["Precision"])
  sensitivity <- as.numeric(cm$byClass["Sensitivity"])
  specificity <- as.numeric(cm$byClass["Specificity"])
  
  f1 <- if (is.na(precision) || is.na(sensitivity) || (precision + sensitivity) == 0) {
    NA_real_
  } else {
    2 * precision * sensitivity / (precision + sensitivity)
  }
  
  list(
    Sensitivity = sensitivity,
    Specificity = specificity,
    F1          = f1
  )
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
  threshold_grid <- seq(0.01, 0.99, by = 0.01)
  
  threshold_results <- lapply(threshold_grid, function(th) {
    pred_class <- factor(ifelse(pred_prob >= th, "pos", "neg"), levels = c("neg", "pos"))
    cm_stats <- confusion_metrics(obs, pred_class)
    prec <- precision_score(obs, pred_class)
    
    data.frame(
      threshold   = th,
      Sensitivity = cm_stats$Sensitivity,
      Specificity = cm_stats$Specificity,
      Precision   = prec,
      F1          = cm_stats$F1
    )
  })
  
  threshold_results <- bind_rows(threshold_results)
  if (nrow(threshold_results) == 0) return(0.50)
  
  threshold_results$Youden <- threshold_results$Sensitivity + threshold_results$Specificity - 1
  
  threshold_results <- threshold_results[order(
    -threshold_results$Youden,
    -threshold_results$Sensitivity,
    -threshold_results$F1,
    threshold_results$threshold
  ), , drop = FALSE]
  
  as.numeric(threshold_results$threshold[1])
}

# =============================================================================
# 4) MISSINGNESS AND IMPUTATION HELPERS
# =============================================================================

induce_mar_ampute <- function(data_x, prop, seed = 123) {
  if (prop <= 0) return(as.data.frame(data_x))
  if (nrow(data_x) == 0 || ncol(data_x) == 0) return(as.data.frame(data_x))
  
  stopifnot(ncol(data_x) == 18)
  
  my_pattern <- matrix(c(
    1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 1)
  
  my_weights <- matrix(c(
    1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 1)
  
  set.seed(seed)
  
  amp_result <- mice::ampute(
    data     = as.data.frame(data_x),
    prop     = prop,
    bycases  = FALSE,
    mech     = "MAR",
    weights =   my_weights,
    patterns = my_pattern
  )
  
  out <- as.data.frame(amp_result$amp)
  names(out) <- names(data_x)
  out
}

impute_train_test_trainonly <- function(train_x, test_x, method_name, m = 5, seed = 123, maxit = 10) {
  train_x <- as.data.frame(train_x)
  test_x  <- as.data.frame(test_x)
  
  n_train <- nrow(train_x)
  n_test  <- nrow(test_x)
  
  if (n_train == 0 || n_test == 0) {
    stop("Both train_x and test_x must have positive number of rows.")
  }
  
  set.seed(seed)
  mids_train <- mice(
    data      = train_x,
    method    = method_name,
    m         = m,
    maxit     = maxit,
    printFlag = FALSE
  )
  
  train_list <- vector("list", m)
  test_list  <- vector("list", m)
  
  for (imp_idx in seq_len(m)) {
    train_completed <- complete(mids_train, action = imp_idx)
    
    combined_data <- bind_rows(train_completed, test_x)
    
    set.seed(seed + imp_idx)
    mids_test <- mice.mids(
      obj       = mids_train,
      newdata   = combined_data,
      maxit     = 1,
      printFlag = FALSE
    )
    
    completed_full <- complete(mids_test, action = imp_idx)
    completed_full[1:n_train, ] <- train_completed
    
    train_list[[imp_idx]] <- completed_full[1:n_train, , drop = FALSE]
    test_list[[imp_idx]]  <- completed_full[(n_train + 1):(n_train + n_test), , drop = FALSE]
  }
  
  list(train = train_list, test = test_list)
}

winsorize_complete_data <- function(train_x_complete) {
  numeric_cols <- names(train_x_complete)[sapply(train_x_complete, is.numeric)]
  
  limits_list <- list()
  for (col in numeric_cols) {
    stats <- boxplot.stats(train_x_complete[[col]])$stats
    limits_list[[col]] <- list(low = stats[1], high = stats[5])
  }
  
  return(limits_list)
}

apply_winsorization <- function(data_x, limits_list) {
  numeric_cols <- names(limits_list)
  for (col in numeric_cols) {
    if (col %in% names(data_x)) {
      limits <- limits_list[[col]]
      data_x[[col]] <- pmin(pmax(data_x[[col]], limits$low), limits$high)
    }
  }
  return(data_x)
}

summarize_missingness <- function(df_x, dataset_label, scenario_name, split_id, target_prop) {
  var_tbl <- data.frame(
    Scenario           = scenario_name,
    Split              = split_id,
    Dataset            = dataset_label,
    Missing_Proportion = target_prop,
    Variable           = names(df_x),
    Missing_Count      = sapply(df_x, function(x) sum(is.na(x))),
    Total_Count        = nrow(df_x),
    Missing_Percent    = sapply(df_x, function(x) mean(is.na(x))),
    stringsAsFactors   = FALSE
  )
  
  overall_tbl <- data.frame(
    Scenario                = scenario_name,
    Split                   = split_id,
    Dataset                 = dataset_label,
    Missing_Proportion      = target_prop,
    Total_Cells             = nrow(df_x) * ncol(df_x),
    Missing_Cells           = sum(is.na(as.matrix(df_x))),
    Overall_Missing_Percent = mean(is.na(as.matrix(df_x))),
    stringsAsFactors        = FALSE
  )
  
  list(
    variable_level = var_tbl,
    overall_level  = overall_tbl
  )
}

create_missingness_comparison_table <- function(overall_missingness_tbl) {
  overall_missingness_tbl %>%
    mutate(
      Target_Missing_Percent = Missing_Proportion,
      Gap_Overall            = Overall_Missing_Percent - Target_Missing_Percent
    ) %>%
    arrange(Scenario, Split, Missing_Proportion, Dataset)
}

# =============================================================================
# 5) MODELING HELPERS
# =============================================================================

get_knn_grid <- function() {
  fixed_knn_grid
}

get_dt_grid <- function() {
  expand.grid(
    cp = c(0.0005, 0.001, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.1)
  )
}

fit_model_with_cv_predictions <- function(train_df, test_df, model_type, seed, cv_folds = 5) {
  train_df <- as.data.frame(train_df)
  test_df  <- as.data.frame(test_df)
  
  train_control <- trainControl(
    method          = "cv",
    number          = cv_folds,
    classProbs      = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  
  set.seed(seed)
  
  model_fit <- tryCatch({
    if (model_type == "GLM") {
      suppressWarnings(caret::train(y ~ ., data = train_df, method = "glm",
                                    family = "binomial", trControl = train_control, metric = "ROC"))
    } else if (model_type == "KNN") {
      suppressWarnings(caret::train(y ~ ., data = train_df, method = "knn",
                                    trControl = train_control, metric = "ROC",
                                    preProcess = c("center", "scale"), tuneGrid = get_knn_grid()))
    } else if (model_type == "DT") {
      suppressWarnings(caret::train(y ~ ., data = train_df, method = "rpart",
                                    trControl = train_control, metric = "ROC", tuneGrid = get_dt_grid()))
    } else {
      stop("Unsupported model_type")
    }
  }, error = function(e) {
    message("  [fit_model_with_cv_predictions] caret failed (", model_type, "): ", e$message)
    NULL
  })
  
  if (is.null(model_fit)) {
    n_train <- nrow(train_df)
    base_rate <- mean(train_df$y == "pos", na.rm = TRUE)
    dummy_cv <- data.frame(rowIndex = seq_len(n_train),
                           obs      = train_df$y,
                           pos      = rep(base_rate, n_train))
    dummy_test <- rep(base_rate, nrow(test_df))
    return(list(model = NULL, cv_pred = dummy_cv, test_prob = dummy_test))
  }
  
  cv_pred <- model_fit$pred
  if (!is.null(model_fit$bestTune) && ncol(model_fit$bestTune) > 0) {
    for (nm in names(model_fit$bestTune)) {
      cv_pred <- cv_pred[cv_pred[[nm]] == model_fit$bestTune[[nm]], , drop = FALSE]
    }
  }
  
  cv_pred <- cv_pred %>%
    arrange(rowIndex) %>%
    select(rowIndex, obs, pos)
  
  test_prob <- predict(model_fit, newdata = test_df, type = "prob")$pos
  
  list(
    model     = model_fit,
    cv_pred   = cv_pred,
    test_prob = test_prob
  )
}

fit_classical_stacking_base <- function(train_df, test_df, model_type, seed, cv_folds = 5) {
  fit_obj <- fit_model_with_cv_predictions(
    train_df   = train_df,
    test_df    = test_df,
    model_type = model_type,
    seed       = seed,
    cv_folds   = cv_folds
  )
  
  oof_df <- fit_obj$cv_pred %>%
    arrange(rowIndex) %>%
    transmute(
      row_id = rowIndex,
      pred   = pos
    )
  
  list(
    oof_df    = oof_df,
    test_pred = as.numeric(fit_obj$test_prob)
  )
}

# =============================================================================
# 6) META-LEARNER FOR STACKING
# =============================================================================

fit_meta_learner_oof <- function(meta_train_df, meta_test_df, seed, cv_folds = 5) {
  meta_train_df <- as.data.frame(meta_train_df)
  meta_test_df  <- as.data.frame(meta_test_df)
  
  meta_cols <- setdiff(names(meta_train_df), "y")
  
  if (length(meta_cols) == 0 || nrow(meta_train_df) == 0) {
    return(list(
      train_oof_pred = numeric(0),
      test_pred      = rep(NA_real_, nrow(meta_test_df))
    ))
  }
  
  meta_train_df$y <- factor(meta_train_df$y, levels = c("neg", "pos"))
  
  if (length(unique(meta_train_df$y)) < 2) {
    base_rate <- mean(meta_train_df$y == "pos")
    return(list(
      train_oof_pred = rep(base_rate, nrow(meta_train_df)),
      test_pred      = rep(base_rate, nrow(meta_test_df))
    ))
  }
  
  set.seed(seed)
  fold_index <- caret::createFolds(meta_train_df$y, k = cv_folds, returnTrain = FALSE)
  train_oof_pred <- rep(NA_real_, nrow(meta_train_df))
  
  for (fold_id in seq_along(fold_index)) {
    valid_idx <- fold_index[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(meta_train_df)), valid_idx)
    
    fold_train <- meta_train_df[train_idx, , drop = FALSE]
    fold_valid <- meta_train_df[valid_idx, , drop = FALSE]
    
    if (length(unique(fold_train$y)) < 2) {
      train_oof_pred[valid_idx] <- mean(fold_train$y == "pos")
      next
    }
    
    fold_fit <- tryCatch(
      glm(y ~ ., data = fold_train, family = binomial()),
      error = function(e) NULL
    )
    
    if (is.null(fold_fit)) {
      fold_pred <- rowMeans(as.matrix(fold_valid[, meta_cols, drop = FALSE]), na.rm = TRUE)
      fold_pred[!is.finite(fold_pred)] <- mean(fold_train$y == "pos")
      train_oof_pred[valid_idx] <- fold_pred
    } else {
      fold_pred <- tryCatch(
        as.numeric(predict(fold_fit, newdata = fold_valid, type = "response")),
        error = function(e) rep(mean(fold_train$y == "pos"), nrow(fold_valid))
      )
      fold_pred[!is.finite(fold_pred)] <- mean(fold_train$y == "pos")
      train_oof_pred[valid_idx] <- fold_pred
    }
  }
  
  train_oof_pred[!is.finite(train_oof_pred)] <- mean(meta_train_df$y == "pos")
  
  final_fit <- tryCatch(
    glm(y ~ ., data = meta_train_df, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(final_fit)) {
    test_pred <- rowMeans(as.matrix(meta_test_df[, meta_cols, drop = FALSE]), na.rm = TRUE)
    test_pred[!is.finite(test_pred)] <- mean(meta_train_df$y == "pos")
  } else {
    test_pred <- tryCatch(
      as.numeric(predict(final_fit, newdata = meta_test_df, type = "response")),
      error = function(e) rowMeans(as.matrix(meta_test_df[, meta_cols, drop = FALSE]), na.rm = TRUE)
    )
    test_pred[!is.finite(test_pred)] <- mean(meta_train_df$y == "pos")
  }
  
  list(
    train_oof_pred = train_oof_pred,
    test_pred      = test_pred
  )
}

# =============================================================================
# 7) MAIN EXPERIMENT LOOP
# =============================================================================

all_results               <- list()
all_test_tables           <- list()
all_missingness_variable  <- list()
all_missingness_overall   <- list()

for (scenario_name in names(source_datasets)) {
  dat <- source_datasets[[scenario_name]]
  missing_props <- missing_props_by_scenario[[scenario_name]]
  
  # ---- CHECKPOINT: load previous progress for this scenario ----------------
  checkpoint_file <- file.path(output_dir, paste0("checkpoint_", scenario_name, "_raw_results.csv"))
  
  if (file.exists(checkpoint_file)) {
    cat("=== CHECKPOINT FOUND for", scenario_name, "— resuming ===\n")
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
    cat("No checkpoint found for", scenario_name, "— starting fresh.\n\n")
  }
  
  scenario_missingness_variable <- data.frame()
  scenario_missingness_overall  <- data.frame()
  
  cat("============================================================\n")
  cat("Running scenario:", scenario_name, "\n")
  cat("============================================================\n")
  
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
      
      cat("  Random split:", split_id, "of", n_random_splits, "\n")
      
      data_seed <- main_seed + (round(prop * 10) + 1) * 1000 + split_id * 10000
      n_before_split <- nrow(scenario_results)   # track new rows for checkpoint
      
      set.seed(data_seed + 1)
      index <- createDataPartition(dat$y, p = train_fraction, list = FALSE)
      train <- dat[index, , drop = FALSE]
      test  <- dat[-index, , drop = FALSE]
      
      train_y <- factor(train$y, levels = c("neg", "pos"))
      test_y  <- factor(test$y,  levels = c("neg", "pos"))
      train_x_complete <- train[, setdiff(names(train), "y"), drop = FALSE]
      test_x_complete  <- test[,  setdiff(names(test),  "y"), drop = FALSE]
      
      winsor_limits <- winsorize_complete_data(train_x_complete)
      train_x_winsorized <- apply_winsorization(train_x_complete, winsor_limits)
      test_x_winsorized <- apply_winsorization(test_x_complete, winsor_limits)
      
      if (prop > 0) {
        train_x_missing <- induce_mar_ampute(
          data_x = train_x_winsorized,
          prop   = prop,
          seed   = data_seed + 2
        )
        
        test_x_missing <- induce_mar_ampute(
          data_x = test_x_winsorized,
          prop   = prop,
          seed   = data_seed + 3
        )
        
        train_x <- train_x_missing
        test_x  <- test_x_missing
      } else {
        train_x <- train_x_winsorized
        test_x  <- test_x_winsorized
      }
      
      train <- data.frame(y = train_y, train_x, check.names = FALSE)
      test  <- data.frame(y = test_y,  test_x,  check.names = FALSE)
      
      train_y <- factor(train$y, levels = c("neg", "pos"))
      test_y  <- factor(test$y,  levels = c("neg", "pos"))
      train_x <- train[, setdiff(names(train), "y"), drop = FALSE]
      test_x  <- test[,  setdiff(names(test),  "y"), drop = FALSE]
      
      train_missing_summary <- summarize_missingness(
        df_x          = train_x,
        dataset_label = "Train",
        scenario_name = scenario_name,
        split_id      = split_id,
        target_prop   = prop
      )
      
      test_missing_summary <- summarize_missingness(
        df_x          = test_x,
        dataset_label = "Test",
        scenario_name = scenario_name,
        split_id      = split_id,
        target_prop   = prop
      )
      
      scenario_missingness_variable <- bind_rows(
        scenario_missingness_variable,
        train_missing_summary$variable_level,
        test_missing_summary$variable_level
      )
      
      scenario_missingness_overall <- bind_rows(
        scenario_missingness_overall,
        train_missing_summary$overall_level,
        test_missing_summary$overall_level
      )
      
      if (prop == 0) {
        complete_train <- data.frame(train_x, y = train_y)
        complete_test  <- test_x
        
        for (clf in classifier_models) {
          fit_complete <- fit_model_with_cv_predictions(
            train_df   = complete_train,
            test_df    = complete_test,
            model_type = clf,
            seed       = data_seed + 10 + match(clf, classifier_models),
            cv_folds   = cv_folds_standard
          )
          
          threshold <- find_optimal_threshold(
            obs       = fit_complete$cv_pred$obs,
            pred_prob = fit_complete$cv_pred$pos
          )
          
          train_prob  <- fit_complete$cv_pred$pos
          train_class <- factor(ifelse(train_prob >= threshold, "pos", "neg"), levels = c("neg", "pos"))
          train_cm_stats <- confusion_metrics(fit_complete$cv_pred$obs, train_class)
          
          test_prob   <- fit_complete$test_prob
          test_class  <- factor(ifelse(test_prob >= threshold, "pos", "neg"), levels = c("neg", "pos"))
          test_cm_stats <- confusion_metrics(test_y, test_class)
          
          scenario_results <- bind_rows(scenario_results, data.frame(
            Scenario           = scenario_name,
            Split              = split_id,
            Missing_Proportion = prop,
            Imputation_Method  = "Complete",
            Model              = clf,
            Train_AUC          = compute_auc(fit_complete$cv_pred$obs, train_prob),
            Test_AUC           = compute_auc(test_y, test_prob),
            Train_Accuracy     = mean(train_class == fit_complete$cv_pred$obs),
            Test_Accuracy      = mean(test_class == test_y),
            Train_Brier        = brier_score(fit_complete$cv_pred$obs, train_prob),
            Test_Brier         = brier_score(test_y, test_prob),
            Train_Precision    = precision_score(fit_complete$cv_pred$obs, train_class),
            Test_Precision     = precision_score(test_y, test_class),
            Train_Sensitivity  = train_cm_stats$Sensitivity,
            Test_Sensitivity   = test_cm_stats$Sensitivity,
            Train_Specificity  = train_cm_stats$Specificity,
            Test_Specificity   = test_cm_stats$Specificity,
            Train_F1           = train_cm_stats$F1,
            Test_F1            = test_cm_stats$F1,
            Threshold          = threshold
          ))
        }
        
        # Save checkpoint for prop==0 split
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
      
      imputed_objs <- list()
      
      for (method in methods) {
        imputed_objs[[method]] <- impute_train_test_trainonly(
          train_x     = train_x,
          test_x      = test_x,
          method_name = method,
          m           = m,
          seed        = data_seed + 100 + match(method, methods),
          maxit       = 10
        )
      }
      
      train_imputed <- list()
      test_imputed  <- list()
      
      for (method in methods) {
        train_imputed[[method]] <- vector("list", m)
        test_imputed[[method]]  <- vector("list", m)
        
        for (imp_idx in seq_len(m)) {
          train_imputed[[method]][[imp_idx]] <- data.frame(
            imputed_objs[[method]]$train[[imp_idx]],
            y = train_y
          )
          test_imputed[[method]][[imp_idx]] <- imputed_objs[[method]]$test[[imp_idx]]
        }
      }
      
      for (method in methods) {
        for (clf in classifier_models) {
          cv_pred_list   <- vector("list", m)
          test_prob_list <- vector("list", m)
          
          for (imp_idx in seq_len(m)) {
            fit_imp <- fit_model_with_cv_predictions(
              train_df   = train_imputed[[method]][[imp_idx]],
              test_df    = test_imputed[[method]][[imp_idx]],
              model_type = clf,
              seed       = data_seed + 1000 +
                match(method, methods) * 100 +
                match(clf, classifier_models) * 10 +
                imp_idx,
              cv_folds   = cv_folds_standard
            )
            
            cv_pred_list[[imp_idx]] <- fit_imp$cv_pred %>%
              rename(!!paste0("pos_imp", imp_idx) := pos)
            
            test_prob_list[[imp_idx]] <- fit_imp$test_prob
          }
          
          pooled_cv <- Reduce(
            function(x, y) full_join(x, y, by = c("rowIndex", "obs")),
            cv_pred_list
          ) %>%
            arrange(rowIndex)
          
          pooled_cv$pooled_prob <- rowMeans(
            as.matrix(pooled_cv[, grep("^pos_imp", names(pooled_cv)), drop = FALSE]),
            na.rm = TRUE
          )
          pooled_cv$pooled_prob[!is.finite(pooled_cv$pooled_prob)] <- mean(train_y == "pos")
          
          model_threshold <- find_optimal_threshold(
            obs       = pooled_cv$obs,
            pred_prob = pooled_cv$pooled_prob
          )
          
          train_prob  <- pooled_cv$pooled_prob
          train_class <- factor(ifelse(train_prob >= model_threshold, "pos", "neg"), levels = c("neg", "pos"))
          train_cm_stats <- confusion_metrics(pooled_cv$obs, train_class)
          
          pooled_test_prob <- rowMeans(do.call(cbind, test_prob_list), na.rm = TRUE)
          pooled_test_prob[!is.finite(pooled_test_prob)] <- mean(train_y == "pos")
          
          test_class <- factor(ifelse(pooled_test_prob >= model_threshold, "pos", "neg"), levels = c("neg", "pos"))
          test_cm_stats <- confusion_metrics(test_y, test_class)
          cal_mice <- compute_calibration(test_y, pooled_test_prob)
          
          scenario_results <- bind_rows(scenario_results, data.frame(
            Scenario           = scenario_name,
            Split              = split_id,
            Missing_Proportion = prop,
            Imputation_Method  = toupper(method),
            Model              = clf,
            Train_AUC          = compute_auc(pooled_cv$obs, train_prob),
            Test_AUC           = compute_auc(test_y, pooled_test_prob),
            Train_Accuracy     = mean(train_class == pooled_cv$obs),
            Test_Accuracy      = mean(test_class == test_y),
            Train_Brier        = brier_score(pooled_cv$obs, train_prob),
            Test_Brier         = brier_score(test_y, pooled_test_prob),
            Train_Precision    = precision_score(pooled_cv$obs, train_class),
            Test_Precision     = precision_score(test_y, test_class),
            Train_Sensitivity  = train_cm_stats$Sensitivity,
            Test_Sensitivity   = test_cm_stats$Sensitivity,
            Train_Specificity  = train_cm_stats$Specificity,
            Test_Specificity   = test_cm_stats$Specificity,
            Train_F1           = train_cm_stats$F1,
            Test_F1            = test_cm_stats$F1,
            Threshold          = model_threshold,
            Cal_Slope          = cal_mice$cal_slope,
            Cal_Intercept      = cal_mice$cal_intercept
          ))
        }
      }
      
      # ---- Simple baselines: CCA, Mean, MissInd, MissForest -------------------
      
      # Helper: store a single-imputation result row for COVID-19
      store_single_result <- function(method_name, clf, cv_p, test_prob_vec) {
        thresh  <- find_optimal_threshold(cv_p$obs, cv_p$pos)
        tr_cl   <- factor(ifelse(cv_p$pos      >= thresh, "pos", "neg"), levels = c("neg", "pos"))
        te_cl   <- factor(ifelse(test_prob_vec >= thresh, "pos", "neg"), levels = c("neg", "pos"))
        tr_cm   <- confusion_metrics(cv_p$obs, tr_cl)
        te_cm   <- confusion_metrics(test_y,   te_cl)
        cal     <- compute_calibration(test_y, test_prob_vec)
        data.frame(
          Scenario           = scenario_name,
          Split              = split_id,
          Missing_Proportion = prop,
          Imputation_Method  = method_name,
          Model              = clf,
          Train_AUC          = compute_auc(cv_p$obs, cv_p$pos),
          Test_AUC           = compute_auc(test_y, test_prob_vec),
          Train_Accuracy     = mean(tr_cl == cv_p$obs),
          Test_Accuracy      = mean(te_cl == test_y),
          Train_Brier        = brier_score(cv_p$obs, cv_p$pos),
          Test_Brier         = brier_score(test_y, test_prob_vec),
          Train_Precision    = precision_score(cv_p$obs, tr_cl),
          Test_Precision     = precision_score(test_y, te_cl),
          Train_Sensitivity  = tr_cm$Sensitivity,
          Test_Sensitivity   = te_cm$Sensitivity,
          Train_Specificity  = tr_cm$Specificity,
          Test_Specificity   = te_cm$Specificity,
          Train_F1           = tr_cm$F1,
          Test_F1            = te_cm$F1,
          Threshold          = thresh,
          Cal_Slope          = cal$cal_slope,
          Cal_Intercept      = cal$cal_intercept
        )
      }
      
      # 1) CCA - complete cases only in training; test missing filled with train means
      cat("    Method: CCA\n")
      cca_complete <- complete.cases(train_x)
      if (sum(cca_complete) >= 10) {
        train_x_cca     <- train_x[cca_complete, , drop = FALSE]
        train_y_cca     <- train_y[cca_complete]
        train_means_cca <- colMeans(train_x, na.rm = TRUE)
        test_x_cca      <- as.data.frame(lapply(names(test_x), function(col) {
          v <- test_x[[col]]; v[is.na(v)] <- train_means_cca[[col]]; v
        })); names(test_x_cca) <- names(test_x)
        train_df_cca <- data.frame(train_x_cca, y = train_y_cca)
        for (clf in classifier_models) {
          f <- fit_model_with_cv_predictions(train_df_cca, test_x_cca, clf,
                                             seed     = data_seed + 30000 + match(clf, classifier_models) * 10,
                                             cv_folds = cv_folds_standard)
          scenario_results <- bind_rows(scenario_results,
                                        store_single_result("CCA", clf, f$cv_pred, f$test_prob))
        }
      }
      
      # 2) Mean imputation - column mean from train, applied to train & test
      cat("    Method: Mean Imputation\n")
      train_means_m <- colMeans(train_x, na.rm = TRUE)
      train_x_mean  <- train_x; test_x_mean <- test_x
      for (col in names(train_x_mean)) {
        train_x_mean[[col]][is.na(train_x_mean[[col]])] <- train_means_m[[col]]
        test_x_mean[[col]][is.na(test_x_mean[[col]])]   <- train_means_m[[col]]
      }
      train_df_mean <- data.frame(train_x_mean, y = train_y)
      for (clf in classifier_models) {
        f <- fit_model_with_cv_predictions(train_df_mean, test_x_mean, clf,
                                           seed     = data_seed + 31000 + match(clf, classifier_models) * 10,
                                           cv_folds = cv_folds_standard)
        scenario_results <- bind_rows(scenario_results,
                                      store_single_result("Mean", clf, f$cv_pred, f$test_prob))
      }
      
      # 3) Missing Indicator - binary flag per missing predictor + mean-fill (Van Ness et al. 2023)
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
        f <- fit_model_with_cv_predictions(train_df_ind, test_x_ind, clf,
                                           seed     = data_seed + 32000 + match(clf, classifier_models) * 10,
                                           cv_folds = cv_folds_standard)
        scenario_results <- bind_rows(scenario_results,
                                      store_single_result("MissInd", clf, f$cv_pred, f$test_prob))
      }
      
      # 4) MissForest - single imputation via random forests (train+test combined)
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
          f <- fit_model_with_cv_predictions(train_df_mf, test_x_mf, clf,
                                             seed     = data_seed + 33000 + match(clf, classifier_models) * 10,
                                             cv_folds = cv_folds_standard)
          scenario_results <- bind_rows(scenario_results,
                                        store_single_result("MissForest", clf, f$cv_pred, f$test_prob))
        }
      }
      
      # ---- Imputation-based stacking ----------------------------------------
      for (clf in classifier_models) {
        stacking_train_list <- vector("list", m)
        stacking_test_mat   <- matrix(NA_real_, nrow = nrow(test_x), ncol = m)
        
        for (imp_idx in seq_len(m)) {
          oof_preds_list  <- list()
          test_preds_list <- list()
          
          for (method in methods) {
            fit_stack_base <- fit_classical_stacking_base(
              train_df   = train_imputed[[method]][[imp_idx]],
              test_df    = test_imputed[[method]][[imp_idx]],
              model_type = clf,
              seed       = data_seed + 5000 +
                match(method, methods) * 100 +
                match(clf, classifier_models) * 1000 +
                imp_idx,
              cv_folds   = cv_folds_stacking
            )
            
            oof_preds_list[[method]] <- fit_stack_base$oof_df
            names(oof_preds_list[[method]])[names(oof_preds_list[[method]]) == "pred"] <- paste0(method, "_pred")
            test_preds_list[[method]] <- fit_stack_base$test_pred
          }
          
          combined_oof <- Reduce(function(x, y) full_join(x, y, by = "row_id"), oof_preds_list) %>%
            arrange(row_id)
          
          meta_cols <- grep("_pred$", names(combined_oof), value = TRUE)
          
          meta_train <- data.frame(
            row_id = combined_oof$row_id,
            y      = train_y[combined_oof$row_id],
            combined_oof[, meta_cols, drop = FALSE]
          )
          
          meta_test_df <- as.data.frame(do.call(cbind, test_preds_list))
          colnames(meta_test_df) <- paste0(methods, "_pred")
          
          meta_train_df <- meta_train[, c("y", meta_cols), drop = FALSE]
          meta_test_df  <- meta_test_df[, meta_cols, drop = FALSE]
          
          meta_fit <- fit_meta_learner_oof(
            meta_train_df = meta_train_df,
            meta_test_df  = meta_test_df,
            seed          = data_seed + 9000 + match(clf, classifier_models) * 100 + imp_idx,
            cv_folds      = cv_folds_meta
          )
          
          stacking_train_list[[imp_idx]] <- data.frame(
            row_id = meta_train$row_id,
            obs    = meta_train$y,
            pred   = meta_fit$train_oof_pred
          )
          
          stacking_test_mat[, imp_idx] <- meta_fit$test_pred
        }
        
        pooled_stack_train <- Reduce(
          function(x, y) full_join(x, y, by = c("row_id", "obs")),
          lapply(seq_along(stacking_train_list), function(i) {
            out <- stacking_train_list[[i]]
            names(out)[names(out) == "pred"] <- paste0("pred_imp", i)
            out
          })
        ) %>%
          arrange(row_id)
        
        pooled_stack_train$pooled_pred <- rowMeans(
          as.matrix(pooled_stack_train[, grep("^pred_imp", names(pooled_stack_train)), drop = FALSE]),
          na.rm = TRUE
        )
        pooled_stack_train$pooled_pred[!is.finite(pooled_stack_train$pooled_pred)] <- mean(train_y == "pos")
        
        final_stack_threshold <- find_optimal_threshold(
          obs       = pooled_stack_train$obs,
          pred_prob = pooled_stack_train$pooled_pred
        )
        
        stack_train_prob  <- pooled_stack_train$pooled_pred
        stack_train_class <- factor(ifelse(stack_train_prob >= final_stack_threshold, "pos", "neg"), levels = c("neg", "pos"))
        stack_train_cm    <- confusion_metrics(pooled_stack_train$obs, stack_train_class)
        
        stacking_pred_avg <- rowMeans(stacking_test_mat, na.rm = TRUE)
        stacking_pred_avg[!is.finite(stacking_pred_avg)] <- mean(train_y == "pos")
        
        stacking_class <- factor(ifelse(stacking_pred_avg >= final_stack_threshold, "pos", "neg"), levels = c("neg", "pos"))
        stack_test_cm  <- confusion_metrics(test_y, stacking_class)
        cal_stack      <- compute_calibration(test_y, stacking_pred_avg)
        
        scenario_results <- bind_rows(scenario_results, data.frame(
          Scenario           = scenario_name,
          Split              = split_id,
          Missing_Proportion = prop,
          Imputation_Method  = "Stacking",
          Model              = clf,
          Train_AUC          = compute_auc(pooled_stack_train$obs, stack_train_prob),
          Test_AUC           = compute_auc(test_y, stacking_pred_avg),
          Train_Accuracy     = mean(stack_train_class == pooled_stack_train$obs),
          Test_Accuracy      = mean(stacking_class == test_y),
          Train_Brier        = brier_score(pooled_stack_train$obs, stack_train_prob),
          Test_Brier         = brier_score(test_y, stacking_pred_avg),
          Train_Precision    = precision_score(pooled_stack_train$obs, stack_train_class),
          Test_Precision     = precision_score(test_y, stacking_class),
          Train_Sensitivity  = stack_train_cm$Sensitivity,
          Test_Sensitivity   = stack_test_cm$Sensitivity,
          Train_Specificity  = stack_train_cm$Specificity,
          Test_Specificity   = stack_test_cm$Specificity,
          Train_F1           = stack_train_cm$F1,
          Test_F1            = stack_test_cm$F1,
          Threshold          = final_stack_threshold,
          Cal_Slope          = cal_stack$cal_slope,
          Cal_Intercept      = cal_stack$cal_intercept
        ))
      }
      
      # ---- MIE_SE: Multiple Imputation Ensemble Stacking (Aleryani et al.) ------
      # MICE_SE variant: PMM (m=5 draws) x GLM+KNN+DT = 15 OOF predictions as meta-features.
      # Meta-learner runs once per classifier → 3 rows: MIE_SE_GLM, MIE_SE_KNN, MIE_SE_DT.
      # Base model seeds: 20000 range. Meta-learner seeds: data_seed + 25000 range.
      cat("    Method: MIE_SE (PMM x GLM+KNN+DT, per-clf meta-learner)\n")
      
      mie_se_oof_list  <- list()
      mie_se_test_list <- list()
      mie_se_idx <- 1
      
      for (imp_idx in seq_len(m)) {
        for (clf in classifier_models) {
          fit_mie <- fit_model_with_cv_predictions(
            train_df   = train_imputed[["pmm"]][[imp_idx]],
            test_df    = test_imputed[["pmm"]][[imp_idx]],
            model_type = clf,
            seed       = data_seed + 20000 + match(clf, classifier_models) * 10 + imp_idx,
            cv_folds   = cv_folds_standard
          )
          col_name <- paste0("pmm_k", imp_idx, "_", clf, "_pred")
          mie_se_oof_list[[mie_se_idx]]  <- fit_mie$cv_pred %>%
            transmute(row_id = rowIndex, !!col_name := pos)
          mie_se_test_list[[mie_se_idx]] <- fit_mie$test_prob
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
        f_meta <- fit_model_with_cv_predictions(
          train_df   = meta_train_mie,
          test_df    = mie_se_test_df,
          model_type = clf_meta,
          seed       = data_seed + 25000 + match(clf_meta, classifier_models) * 100,
          cv_folds   = cv_folds_meta
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
          Scenario           = scenario_name,
          Split              = split_id,
          Missing_Proportion = prop,
          Imputation_Method  = mie_label,
          Model              = mie_label,
          Train_AUC          = compute_auc(mie_se_train_obs, mie_se_train_prob),
          Test_AUC           = compute_auc(test_y, mie_se_test_prob),
          Train_Accuracy     = mean(mie_se_train_class == mie_se_train_obs),
          Test_Accuracy      = mean(mie_se_test_class == test_y),
          Train_Brier        = brier_score(mie_se_train_obs, mie_se_train_prob),
          Test_Brier         = brier_score(test_y, mie_se_test_prob),
          Train_Precision    = precision_score(mie_se_train_obs, mie_se_train_class),
          Test_Precision     = precision_score(test_y, mie_se_test_class),
          Train_Sensitivity  = mie_se_train_cm$Sensitivity,
          Test_Sensitivity   = mie_se_test_cm$Sensitivity,
          Train_Specificity  = mie_se_train_cm$Specificity,
          Test_Specificity   = mie_se_test_cm$Specificity,
          Train_F1           = mie_se_train_cm$F1,
          Test_F1            = mie_se_test_cm$F1,
          Threshold          = thresh_mie_se,
          Cal_Slope          = cal_mie$cal_slope,
          Cal_Intercept      = cal_mie$cal_intercept
        ))
      }
      
      # ---- Save checkpoint after each split ----------------------------------
      new_rows <- if (nrow(scenario_results) > n_before_split)
        scenario_results[(n_before_split + 1):nrow(scenario_results), ] else data.frame()
      if (nrow(new_rows) > 0) {
        readr::write_csv(new_rows, checkpoint_file,
                         append = file.exists(checkpoint_file))
        completed_keys <- c(completed_keys, split_key)
        cat("  Checkpoint saved (prop=", prop, ", split=", split_id, ")\n", sep = "")
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
  
  scenario_results$Imputation_Method <- factor(
    scenario_results$Imputation_Method,
    levels = method_order
  )
  
  scenario_results <- scenario_results %>%
    arrange(Split, Missing_Proportion, Model, Imputation_Method)
  
  safe_write_csv(
    scenario_results,
    file.path(output_dir, "scenario_raw_results", paste0(scenario_name, "_raw_results.csv"))
  )
  
  scenario_test_table <- scenario_results %>%
    select(
      Scenario, Split, Missing_Proportion, Imputation_Method, Model,
      Test_AUC, Test_Accuracy, Test_Brier, Test_Precision,
      Test_Sensitivity, Test_Specificity, Test_F1, Threshold
    ) %>%
    arrange(Split, Missing_Proportion, Model, Imputation_Method)
  
  safe_write_csv(
    scenario_test_table,
    file.path(output_dir, "scenario_test_summaries", paste0(scenario_name, "_test_metrics.csv"))
  )
  
  scenario_overfit_table <- scenario_results %>%
    select(
      Scenario, Split, Missing_Proportion, Imputation_Method, Model,
      Overfit_Gap_AUC, Overfit_Gap_Accuracy, Overfit_Gap_Brier,
      Overfit_Gap_Precision, Overfit_Gap_Sensitivity,
      Overfit_Gap_Specificity, Overfit_Gap_F1
    ) %>%
    arrange(Split, Missing_Proportion, Model, Imputation_Method)
  
  safe_write_csv(
    scenario_overfit_table,
    file.path(output_dir, "scenario_overfitting_summaries", paste0(scenario_name, "_overfitting_summary.csv"))
  )
  
  safe_write_csv(
    scenario_missingness_variable,
    file.path(output_dir, paste0(scenario_name, "_actual_missingness_variable_level.csv"))
  )
  
  safe_write_csv(
    scenario_missingness_overall,
    file.path(output_dir, paste0(scenario_name, "_actual_missingness_overall.csv"))
  )
  
  scenario_missingness_comparison <- create_missingness_comparison_table(scenario_missingness_overall)
  
  safe_write_csv(
    scenario_missingness_comparison,
    file.path(output_dir, paste0(scenario_name, "_actual_missingness_comparison.csv"))
  )
  
  all_results[[scenario_name]]               <- scenario_results
  all_test_tables[[scenario_name]]           <- scenario_test_table
  all_missingness_variable[[scenario_name]]  <- scenario_missingness_variable
  all_missingness_overall[[scenario_name]]   <- scenario_missingness_overall
  
  cat("Saved outputs for scenario:", scenario_name, "\n\n")
}

raw_results_file        <- file.path(output_dir, "all_scenarios_raw_results.csv")
test_metrics_file       <- file.path(output_dir, "all_scenarios_test_metrics.csv")
missing_var_file        <- file.path(output_dir, "all_scenarios_actual_missingness_variable_level.csv")
missing_overall_file    <- file.path(output_dir, "all_scenarios_actual_missingness_overall.csv")
missing_comparison_file <- file.path(output_dir, "all_scenarios_actual_missingness_comparison.csv")

if (exists("all_results") &&
    exists("all_test_tables") &&
    exists("all_missingness_variable") &&
    exists("all_missingness_overall")) {
  
  combined_results      <- bind_rows(all_results)
  combined_test_metrics <- bind_rows(all_test_tables)
  
  combined_missingness_variable   <- bind_rows(all_missingness_variable)
  combined_missingness_overall    <- bind_rows(all_missingness_overall)
  combined_missingness_comparison <- create_missingness_comparison_table(combined_missingness_overall)
  
  safe_write_csv(combined_results, raw_results_file)
  safe_write_csv(combined_test_metrics, test_metrics_file)
  safe_write_csv(combined_missingness_variable, missing_var_file)
  safe_write_csv(combined_missingness_overall, missing_overall_file)
  safe_write_csv(combined_missingness_comparison, missing_comparison_file)
  
} else if (file.exists(raw_results_file) &&
           file.exists(test_metrics_file) &&
           file.exists(missing_var_file) &&
           file.exists(missing_overall_file) &&
           file.exists(missing_comparison_file)) {
  
  combined_results                <- read_csv(raw_results_file, show_col_types = FALSE)
  combined_test_metrics           <- read_csv(test_metrics_file, show_col_types = FALSE)
  combined_missingness_variable   <- read_csv(missing_var_file, show_col_types = FALSE)
  combined_missingness_overall    <- read_csv(missing_overall_file, show_col_types = FALSE)
  combined_missingness_comparison <- read_csv(missing_comparison_file, show_col_types = FALSE)
  
} else {
  stop(
    paste(
      "Neither in-memory objects nor saved combined CSV files were found.",
      "Run the full experiment loop first (Section 7), or make sure these files exist:",
      raw_results_file,
      test_metrics_file,
      missing_var_file,
      missing_overall_file,
      missing_comparison_file,
      sep = "\n"
    )
  )
}




















# =============================================================================
# POST-ANALYSIS / SUMMARY / RANKING / TABLES
# =============================================================================

# =============================================================================
# SAFETY FOR MISSING COLUMNS
# =============================================================================

if (!"Eligible_Missing_Percent" %in% names(combined_missingness_overall)) {
  combined_missingness_overall$Eligible_Missing_Percent <- NA_real_
}
if (!"Protected_Missing_Percent" %in% names(combined_missingness_overall)) {
  combined_missingness_overall$Protected_Missing_Percent <- NA_real_
}
if (!"Missing_Cells" %in% names(combined_missingness_overall)) {
  combined_missingness_overall$Missing_Cells <- NA_real_
}
if (!"Total_Cells" %in% names(combined_missingness_overall)) {
  combined_missingness_overall$Total_Cells <- NA_real_
}

# =============================================================================
# 8B) MISSINGNESS SUMMARY ACROSS RANDOM SPLITS
# =============================================================================

missingness_summary <- combined_missingness_overall %>%
  group_by(Scenario, Missing_Proportion, Dataset) %>%
  summarise(
    Mean_Overall_Missingness    = mean(Overall_Missing_Percent, na.rm = TRUE),
    SD_Overall_Missingness      = sd(Overall_Missing_Percent, na.rm = TRUE),
    Mean_Eligible_Missingness   = mean(Eligible_Missing_Percent, na.rm = TRUE),
    SD_Eligible_Missingness     = sd(Eligible_Missing_Percent, na.rm = TRUE),
    Mean_Protected_Missingness  = mean(Protected_Missing_Percent, na.rm = TRUE),
    SD_Protected_Missingness    = sd(Protected_Missing_Percent, na.rm = TRUE),
    Mean_Missing_Cells          = mean(Missing_Cells, na.rm = TRUE),
    SD_Missing_Cells            = sd(Missing_Cells, na.rm = TRUE),
    Mean_Total_Cells            = mean(Total_Cells, na.rm = TRUE),
    SD_Total_Cells              = sd(Total_Cells, na.rm = TRUE),
    n_splits                    = n(),
    .groups = "drop"
  ) %>%
  arrange(Scenario, Missing_Proportion, Dataset)

missingness_variable_summary <- combined_missingness_variable %>%
  group_by(Scenario, Missing_Proportion, Dataset, Variable) %>%
  summarise(
    Mean_Missing_Percent = mean(Missing_Percent, na.rm = TRUE),
    SD_Missing_Percent   = sd(Missing_Percent, na.rm = TRUE),
    Mean_Missing_Count   = mean(Missing_Count, na.rm = TRUE),
    SD_Missing_Count     = sd(Missing_Count, na.rm = TRUE),
    n_splits             = n(),
    .groups = "drop"
  ) %>%
  arrange(Scenario, Missing_Proportion, Dataset, Variable)

safe_write_csv(
  missingness_summary,
  file.path(output_dir, "post_analysis", "missingness_mean_sd_summary.csv")
)

safe_write_csv(
  missingness_variable_summary,
  file.path(output_dir, "post_analysis", "missingness_variable_mean_sd_summary.csv")
)

# =============================================================================
# 8C) TRAIN / TEST MISSINGNESS TABLES
# =============================================================================

missingness_train_test_table <- combined_missingness_overall %>%
  filter(Missing_Proportion > 0) %>%
  mutate(
    Missing_Label = paste0(as.integer(Missing_Proportion * 100), "%")
  ) %>%
  group_by(Scenario, Missing_Proportion, Missing_Label, Dataset) %>%
  summarise(
    Mean_Missingness = mean(Overall_Missing_Percent, na.rm = TRUE),
    SD_Missingness   = sd(Overall_Missing_Percent, na.rm = TRUE),
    n_splits         = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Missingness_Mean_SD = sprintf("%.3f ± %.3f", Mean_Missingness, SD_Missingness)
  ) %>%
  select(Scenario, Missing_Proportion, Missing_Label, Dataset, Missingness_Mean_SD) %>%
  pivot_wider(
    names_from  = Dataset,
    values_from = Missingness_Mean_SD
  ) %>%
  arrange(Scenario, Missing_Proportion) %>%
  select(
    Scenario,
    `Missingness Percentage` = Missing_Label,
    Train,
    Test
  )

safe_write_csv(
  missingness_train_test_table,
  file.path(output_dir, "post_analysis", "missingness_train_test_mean_sd_table.csv")
)

missingness_train_test_numeric <- combined_missingness_overall %>%
  filter(Missing_Proportion > 0) %>%
  mutate(
    Missing_Label = paste0(as.integer(Missing_Proportion * 100), "%")
  ) %>%
  group_by(Scenario, Missing_Proportion, Missing_Label, Dataset) %>%
  summarise(
    Mean_Missingness = mean(Overall_Missing_Percent, na.rm = TRUE),
    SD_Missingness   = sd(Overall_Missing_Percent, na.rm = TRUE),
    n_splits         = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Dataset,
    values_from = c(Mean_Missingness, SD_Missingness),
    names_sep = "_"
  ) %>%
  arrange(Scenario, Missing_Proportion)

safe_write_csv(
  missingness_train_test_numeric,
  file.path(output_dir, "post_analysis", "missingness_train_test_mean_sd_numeric.csv")
)

# =============================================================================
# 9) OVERFITTING SUMMARY
# =============================================================================

overfitting_summary <- combined_results %>%
  filter(!is.na(Missing_Proportion)) %>%
  mutate(
    Method = as.character(Imputation_Method),
    Model  = as.character(Model)
  ) %>%
  group_by(Scenario, Model, Missing_Proportion, Method) %>%
  summarise(
    Mean_Overfit_AUC         = mean(Overfit_Gap_AUC, na.rm = TRUE),
    SD_Overfit_AUC           = sd(Overfit_Gap_AUC, na.rm = TRUE),
    Mean_Overfit_Accuracy    = mean(Overfit_Gap_Accuracy, na.rm = TRUE),
    SD_Overfit_Accuracy      = sd(Overfit_Gap_Accuracy, na.rm = TRUE),
    Mean_Overfit_Brier       = mean(Overfit_Gap_Brier, na.rm = TRUE),
    SD_Overfit_Brier         = sd(Overfit_Gap_Brier, na.rm = TRUE),
    Mean_Overfit_Precision   = mean(Overfit_Gap_Precision, na.rm = TRUE),
    SD_Overfit_Precision     = sd(Overfit_Gap_Precision, na.rm = TRUE),
    Mean_Overfit_Sensitivity = mean(Overfit_Gap_Sensitivity, na.rm = TRUE),
    SD_Overfit_Sensitivity   = sd(Overfit_Gap_Sensitivity, na.rm = TRUE),
    Mean_Overfit_Specificity = mean(Overfit_Gap_Specificity, na.rm = TRUE),
    SD_Overfit_Specificity   = sd(Overfit_Gap_Specificity, na.rm = TRUE),
    Mean_Overfit_F1          = mean(Overfit_Gap_F1, na.rm = TRUE),
    SD_Overfit_F1            = sd(Overfit_Gap_F1, na.rm = TRUE),
    n_splits                 = n(),
    .groups = "drop"
  ) %>%
  arrange(Scenario, Model, Missing_Proportion, Method)

safe_write_csv(
  as.data.frame(overfitting_summary),
  file.path(output_dir, "post_analysis", "overfitting_mean_sd_summary.csv")
)

# =============================================================================
# 10) MAIN PERFORMANCE SUMMARY
# =============================================================================

missingness_join_tbl <- missingness_summary %>%
  select(
    Scenario,
    Missing_Proportion,
    Dataset,
    Mean_Overall_Missingness,
    SD_Overall_Missingness,
    Mean_Eligible_Missingness,
    SD_Eligible_Missingness,
    Mean_Protected_Missingness,
    SD_Protected_Missingness
  ) %>%
  pivot_wider(
    names_from = Dataset,
    values_from = c(
      Mean_Overall_Missingness,
      SD_Overall_Missingness,
      Mean_Eligible_Missingness,
      SD_Eligible_Missingness,
      Mean_Protected_Missingness,
      SD_Protected_Missingness
    ),
    names_sep = "_"
  )

final_results <- combined_results %>%
  filter(!is.na(Missing_Proportion)) %>%
  mutate(
    Method = as.character(Imputation_Method),
    Model  = as.character(Model)
  ) %>%
  group_by(Scenario, Model, Missing_Proportion, Method) %>%
  summarise(
    AUC            = mean(Test_AUC, na.rm = TRUE),
    SD_AUC         = sd(Test_AUC, na.rm = TRUE),
    Accuracy       = mean(Test_Accuracy, na.rm = TRUE),
    SD_Accuracy    = sd(Test_Accuracy, na.rm = TRUE),
    Brier          = mean(Test_Brier, na.rm = TRUE),
    SD_Brier       = sd(Test_Brier, na.rm = TRUE),
    Precision      = mean(Test_Precision, na.rm = TRUE),
    SD_Precision   = sd(Test_Precision, na.rm = TRUE),
    Sensitivity    = mean(Test_Sensitivity, na.rm = TRUE),
    SD_Sensitivity = sd(Test_Sensitivity, na.rm = TRUE),
    Specificity    = mean(Test_Specificity, na.rm = TRUE),
    SD_Specificity = sd(Test_Specificity, na.rm = TRUE),
    F1             = mean(Test_F1, na.rm = TRUE),
    SD_F1          = sd(Test_F1, na.rm = TRUE),
    Mean_Threshold = mean(Threshold, na.rm = TRUE),
    SD_Threshold   = sd(Threshold, na.rm = TRUE),
    n_splits       = n(),
    .groups = "drop"
  ) %>%
  left_join(
    overfitting_summary,
    by = c("Scenario", "Model", "Missing_Proportion", "Method")
  ) %>%
  left_join(
    missingness_join_tbl,
    by = c("Scenario", "Missing_Proportion")
  ) %>%
  mutate(
    AUC_Mean_SD         = sprintf("%.4f ± %.4f", AUC, SD_AUC),
    Accuracy_Mean_SD    = sprintf("%.4f ± %.4f", Accuracy, SD_Accuracy),
    Brier_Mean_SD       = sprintf("%.4f ± %.4f", Brier, SD_Brier),
    Precision_Mean_SD   = sprintf("%.4f ± %.4f", Precision, SD_Precision),
    Sensitivity_Mean_SD = sprintf("%.4f ± %.4f", Sensitivity, SD_Sensitivity),
    Specificity_Mean_SD = sprintf("%.4f ± %.4f", Specificity, SD_Specificity),
    F1_Mean_SD          = sprintf("%.4f ± %.4f", F1, SD_F1),
    Threshold_Mean_SD   = sprintf("%.4f ± %.4f", Mean_Threshold, SD_Threshold)
  )

final_results$Method <- factor(
  final_results$Method,
  levels = c("Complete", "PMM", "RF", "CART", "NORM", "MIDASTOUCH", "Stacking")
)

final_results$Model <- factor(
  final_results$Model,
  levels = c("GLM", "KNN", "DT")
)

safe_write_csv(
  as.data.frame(final_results),
  file.path(output_dir, "post_analysis", "final_results_by_model.csv")
)

glm_results <- final_results %>% filter(Model == "GLM")
knn_results <- final_results %>% filter(Model == "KNN")
dt_results  <- final_results %>% filter(Model == "DT")

safe_write_csv(as.data.frame(glm_results), file.path(output_dir, "post_analysis", "GLM_mean_sd_results.csv"))
safe_write_csv(as.data.frame(knn_results), file.path(output_dir, "post_analysis", "KNN_mean_sd_results.csv"))
safe_write_csv(as.data.frame(dt_results),  file.path(output_dir, "post_analysis", "DT_mean_sd_results.csv"))

# =============================================================================

# =============================================================================
# 12) RANKING FUNCTION
# =============================================================================

create_rank_table <- function(dt, metric, higher_better = TRUE) {
  dt_copy <- copy(dt)
  dt_copy[[metric]] <- as.numeric(dt_copy[[metric]])
  
  wide <- dcast(dt_copy, Method ~ Missing_Proportion, value.var = metric)
  prop_cols <- setdiff(names(wide), "Method")
  
  wide[, Mean_Value := rowMeans(.SD, na.rm = TRUE), .SDcols = prop_cols]
  
  if (higher_better) {
    wide[, Overall_Rank := rank(-Mean_Value, ties.method = "min")]
  } else {
    wide[, Overall_Rank := rank(Mean_Value, ties.method = "min")]
  }
  
  rank_cols <- paste0("Rank_", prop_cols)
  
  if (higher_better) {
    wide[, (rank_cols) := lapply(.SD, function(x) rank(-x, ties.method = "average")), .SDcols = prop_cols]
  } else {
    wide[, (rank_cols) := lapply(.SD, function(x) rank(x, ties.method = "average")), .SDcols = prop_cols]
  }
  
  rank_matrix <- as.matrix(wide[, ..rank_cols])
  kendall_res <- suppressWarnings(kendall(rank_matrix))
  
  list(
    rank_table = wide[, c("Method", rank_cols, "Overall_Rank"), with = FALSE],
    kendall    = kendall_res
  )
}

# =============================================================================
# 13) RANKING ANALYSIS BY MODEL
# =============================================================================

rank_results <- final_results %>%
  filter(Method != "Complete") %>%
  group_by(Model, Method, Missing_Proportion) %>%
  summarise(
    AUC         = mean(AUC, na.rm = TRUE),
    Accuracy    = mean(Accuracy, na.rm = TRUE),
    Brier       = mean(Brier, na.rm = TRUE),
    Precision   = mean(Precision, na.rm = TRUE),
    Sensitivity = mean(Sensitivity, na.rm = TRUE),
    Specificity = mean(Specificity, na.rm = TRUE),
    F1          = mean(F1, na.rm = TRUE),
    .groups = "drop"
  )

rank_glm_results <- rank_results %>% filter(Model == "GLM")
rank_knn_results <- rank_results %>% filter(Model == "KNN")
rank_dt_results  <- rank_results %>% filter(Model == "DT")

setDT(rank_glm_results)
setDT(rank_knn_results)
setDT(rank_dt_results)

glm_auc_res         <- create_rank_table(rank_glm_results, "AUC",         TRUE)
glm_acc_res         <- create_rank_table(rank_glm_results, "Accuracy",    TRUE)
glm_brier_res       <- create_rank_table(rank_glm_results, "Brier",       FALSE)
glm_precision_res   <- create_rank_table(rank_glm_results, "Precision",   TRUE)
glm_sensitivity_res <- create_rank_table(rank_glm_results, "Sensitivity", TRUE)
glm_specificity_res <- create_rank_table(rank_glm_results, "Specificity", TRUE)
glm_f1_res          <- create_rank_table(rank_glm_results, "F1",          TRUE)

knn_auc_res         <- create_rank_table(rank_knn_results, "AUC",         TRUE)
knn_acc_res         <- create_rank_table(rank_knn_results, "Accuracy",    TRUE)
knn_brier_res       <- create_rank_table(rank_knn_results, "Brier",       FALSE)
knn_precision_res   <- create_rank_table(rank_knn_results, "Precision",   TRUE)
knn_sensitivity_res <- create_rank_table(rank_knn_results, "Sensitivity", TRUE)
knn_specificity_res <- create_rank_table(rank_knn_results, "Specificity", TRUE)
knn_f1_res          <- create_rank_table(rank_knn_results, "F1",          TRUE)

dt_auc_res          <- create_rank_table(rank_dt_results, "AUC",         TRUE)
dt_acc_res          <- create_rank_table(rank_dt_results, "Accuracy",    TRUE)
dt_brier_res        <- create_rank_table(rank_dt_results, "Brier",       FALSE)
dt_precision_res    <- create_rank_table(rank_dt_results, "Precision",   TRUE)
dt_sensitivity_res  <- create_rank_table(rank_dt_results, "Sensitivity", TRUE)
dt_specificity_res  <- create_rank_table(rank_dt_results, "Specificity", TRUE)
dt_f1_res           <- create_rank_table(rank_dt_results, "F1",          TRUE)

# =============================================================================
# 14) HELPER TO COMBINE METRICS INTO ONE TABLE
# =============================================================================

extract_metric_table <- function(res_obj, metric_name) {
  dt <- as.data.frame(res_obj$rank_table)
  dt$Metric  <- metric_name
  dt$W       <- round(as.numeric(res_obj$kendall$value), 3)
  dt$p_value <- signif(as.numeric(res_obj$kendall$p.value), 4)
  
  dt[, c(
    "Metric", "Method",
    "Rank_0.1", "Rank_0.2", "Rank_0.3", "Rank_0.4", "Rank_0.5",
    "Overall_Rank", "W", "p_value"
  )]
}

# =============================================================================
# 15) COMBINED RANK TABLES
# =============================================================================

combined_rank_table_glm <- rbind(
  extract_metric_table(glm_auc_res, "AUC"),
  extract_metric_table(glm_acc_res, "Accuracy"),
  extract_metric_table(glm_brier_res, "Brier"),
  extract_metric_table(glm_precision_res, "Precision"),
  extract_metric_table(glm_sensitivity_res, "Sensitivity"),
  extract_metric_table(glm_specificity_res, "Specificity"),
  extract_metric_table(glm_f1_res, "F1")
)

combined_rank_table_knn <- rbind(
  extract_metric_table(knn_auc_res, "AUC"),
  extract_metric_table(knn_acc_res, "Accuracy"),
  extract_metric_table(knn_brier_res, "Brier"),
  extract_metric_table(knn_precision_res, "Precision"),
  extract_metric_table(knn_sensitivity_res, "Sensitivity"),
  extract_metric_table(knn_specificity_res, "Specificity"),
  extract_metric_table(knn_f1_res, "F1")
)

combined_rank_table_dt <- rbind(
  extract_metric_table(dt_auc_res, "AUC"),
  extract_metric_table(dt_acc_res, "Accuracy"),
  extract_metric_table(dt_brier_res, "Brier"),
  extract_metric_table(dt_precision_res, "Precision"),
  extract_metric_table(dt_sensitivity_res, "Sensitivity"),
  extract_metric_table(dt_specificity_res, "Specificity"),
  extract_metric_table(dt_f1_res, "F1")
)

safe_write_csv(combined_rank_table_glm, file.path(output_dir, "post_analysis", "GLM_combined_rank_table.csv"))
safe_write_csv(combined_rank_table_knn, file.path(output_dir, "post_analysis", "KNN_combined_rank_table.csv"))
safe_write_csv(combined_rank_table_dt,  file.path(output_dir, "post_analysis", "DT_combined_rank_table.csv"))

# =============================================================================
# 16) SAVE MEAN ± SD TABLES
# =============================================================================

fmt3 <- function(x) sprintf("%.3f", x)

prepare_main_results_table <- function(df_model_name) {
  final_results %>%
    filter(Model == df_model_name) %>%
    mutate(
      `Missingness Percentage` = paste0(as.integer(Missing_Proportion * 100), "%"),
      AUC         = paste0(fmt3(AUC), " ± ", fmt3(SD_AUC)),
      Accuracy    = paste0(fmt3(Accuracy), " ± ", fmt3(SD_Accuracy)),
      Brier       = paste0(fmt3(Brier), " ± ", fmt3(SD_Brier)),
      Precision   = paste0(fmt3(Precision), " ± ", fmt3(SD_Precision)),
      Sensitivity = paste0(fmt3(Sensitivity), " ± ", fmt3(SD_Sensitivity)),
      Specificity = paste0(fmt3(Specificity), " ± ", fmt3(SD_Specificity)),
      F1          = paste0(fmt3(F1), " ± ", fmt3(SD_F1))
    ) %>%
    select(
      Scenario,
      `Missingness Percentage`,
      Method,
      AUC,
      Accuracy,
      Brier,
      Precision,
      Sensitivity,
      Specificity,
      F1
    )
}

main_glm_word <- prepare_main_results_table("GLM")
main_knn_word <- prepare_main_results_table("KNN")
main_dt_word  <- prepare_main_results_table("DT")

safe_write_csv(main_glm_word, file.path(output_dir, "post_analysis", "GLM_mean_sd_table.csv"))
safe_write_csv(main_knn_word, file.path(output_dir, "post_analysis", "KNN_mean_sd_table.csv"))
safe_write_csv(main_dt_word,  file.path(output_dir, "post_analysis", "DT_mean_sd_table.csv"))

prepare_overfit_table <- function(df_model_name) {
  overfitting_summary %>%
    filter(Model == df_model_name) %>%
    mutate(
      `Missingness Percentage` = paste0(as.integer(Missing_Proportion * 100), "%"),
      AUC         = paste0(fmt3(Mean_Overfit_AUC), " ± ", fmt3(SD_Overfit_AUC)),
      Accuracy    = paste0(fmt3(Mean_Overfit_Accuracy), " ± ", fmt3(SD_Overfit_Accuracy)),
      Brier       = paste0(fmt3(Mean_Overfit_Brier), " ± ", fmt3(SD_Overfit_Brier)),
      Precision   = paste0(fmt3(Mean_Overfit_Precision), " ± ", fmt3(SD_Overfit_Precision)),
      Sensitivity = paste0(fmt3(Mean_Overfit_Sensitivity), " ± ", fmt3(SD_Overfit_Sensitivity)),
      Specificity = paste0(fmt3(Mean_Overfit_Specificity), " ± ", fmt3(SD_Overfit_Specificity)),
      F1          = paste0(fmt3(Mean_Overfit_F1), " ± ", fmt3(SD_Overfit_F1))
    ) %>%
    select(
      Scenario,
      `Missingness Percentage`,
      Method,
      AUC,
      Accuracy,
      Brier,
      Precision,
      Sensitivity,
      Specificity,
      F1
    )
}

overfit_glm_word <- prepare_overfit_table("GLM")
overfit_knn_word <- prepare_overfit_table("KNN")
overfit_dt_word  <- prepare_overfit_table("DT")

safe_write_csv(overfit_glm_word, file.path(output_dir, "post_analysis", "GLM_overfit_mean_sd_table.csv"))
safe_write_csv(overfit_knn_word, file.path(output_dir, "post_analysis", "KNN_overfit_mean_sd_table.csv"))
safe_write_csv(overfit_dt_word,  file.path(output_dir, "post_analysis", "DT_overfit_mean_sd_table.csv"))

# =============================================================================
# 17) SAVE BACKUP
# =============================================================================

save(
  final_results,
  combined_results,
  combined_missingness_overall,
  combined_missingness_variable,
  missingness_summary,
  missingness_variable_summary,
  missingness_train_test_table,
  missingness_train_test_numeric,
  overfitting_summary,
  combined_rank_table_glm,
  combined_rank_table_knn,
  combined_rank_table_dt,
  main_glm_word,
  main_knn_word,
  main_dt_word,
  overfit_glm_word,
  overfit_knn_word,
  overfit_dt_word,
  file = file.path(output_dir, "post_analysis", "results_backup.RData")
)

# =============================================================================
# 18) CONSOLE MESSAGE
# =============================================================================

cat("\n==================== OUTPUTS SAVED ====================\n")
cat("Output directory:", output_dir, "\n")
cat("\nGenerated files in post_analysis folder:\n")
cat("  - missingness_mean_sd_summary.csv\n")
cat("  - missingness_variable_mean_sd_summary.csv\n")
cat("  - missingness_train_test_mean_sd_table.csv\n")
cat("  - missingness_train_test_mean_sd_numeric.csv\n")
cat("  - final_results_by_model.csv\n")
cat("  - GLM_mean_sd_results.csv\n")
cat("  - KNN_mean_sd_results.csv\n")
cat("  - DT_mean_sd_results.csv\n")
cat("  - GLM_mean_sd_table.csv\n")
cat("  - KNN_mean_sd_table.csv\n")
cat("  - DT_mean_sd_table.csv\n")
cat("  - GLM_overfit_mean_sd_table.csv\n")
cat("  - KNN_overfit_mean_sd_table.csv\n")
cat("  - DT_overfit_mean_sd_table.csv\n")
cat("  - overfitting_mean_sd_summary.csv\n")
cat("  - GLM_combined_rank_table.csv\n")
cat("  - KNN_combined_rank_table.csv\n")
cat("  - DT_combined_rank_table.csv\n")
cat("
Run R/02_covid19_figures.R to generate manuscript figures.
")
cat("  - results_backup.RData\n")
