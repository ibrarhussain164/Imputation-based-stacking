# =============================================================================
# ABM analysis pipeline
# Dataset: ABM / HBiostat meningitis dataset
# Outcome: abm (0 = acute viral meningitis, 1 = acute bacterial meningitis)
# Purpose: compare individual imputation methods and imputation-based stacking
# under naturally occurring clinical missingness.
# =============================================================================

# =============================================================================
# 0) CLEAN SESSION AND LOAD PACKAGES
# =============================================================================

rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(caret)
  library(pROC)
  library(mice)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(haven)
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

# =============================================================================
# SAFE FILE WRITERS
# =============================================================================

safe_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  df <- as.data.frame(df)
  
  tryCatch(
    {
      readr::write_csv(df, path)
      message("Saved: ", path)
      return(invisible(path))
    },
    error = function(e1) {
      alt_path <- paste0(
        tools::file_path_sans_ext(path),
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".csv"
      )
      
      message("Could not save to: ", path)
      message("Reason: ", e1$message)
      message("Trying alternative path: ", alt_path)
      
      tryCatch(
        {
          utils::write.csv(df, alt_path, row.names = FALSE)
          message("Saved alternative: ", alt_path)
          return(invisible(alt_path))
        },
        error = function(e2) {
          message("Alternative save also failed: ", e2$message)
          return(invisible(NULL))
        }
      )
    }
  )
}

safe_write_text <- function(text, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    writeLines(text, con = path, useBytes = TRUE),
    error = function(e) {
      alt_path <- paste0(
        tools::file_path_sans_ext(path),
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".txt"
      )
      writeLines(text, con = alt_path, useBytes = TRUE)
    }
  )
}

safe_save_rdata <- function(object_names, path, env = parent.frame()) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  tryCatch(
    {
      save(list = object_names, file = path, envir = env)
      message("Saved RData: ", path)
      return(invisible(path))
    },
    error = function(e1) {
      alt_path <- paste0(
        tools::file_path_sans_ext(path),
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".RData"
      )
      message("Could not save RData to: ", path)
      message("Reason: ", e1$message)
      message("Trying alternative RData path: ", alt_path)
      
      tryCatch(
        {
          save(list = object_names, file = alt_path, envir = env)
          message("Saved alternative RData: ", alt_path)
          return(invisible(alt_path))
        },
        error = function(e2) {
          message("Alternative RData save also failed: ", e2$message)
          return(invisible(NULL))
        }
      )
    }
  )
}

safe_ggsave <- function(filename, plot, width = 14, height = 12, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  tryCatch(
    {
      ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi)
      message("Saved figure: ", filename)
      return(invisible(filename))
    },
    error = function(e1) {
      alt_path <- paste0(
        tools::file_path_sans_ext(filename),
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".png"
      )
      message("Could not save figure to: ", filename)
      message("Reason: ", e1$message)
      message("Trying alternative figure path: ", alt_path)
      ggplot2::ggsave(filename = alt_path, plot = plot, width = width, height = height, dpi = dpi)
      return(invisible(alt_path))
    }
  )
}

# =============================================================================
# 1) USER CONFIGURATION
# =============================================================================

base_dir <- Sys.getenv("ABM_OUTPUT_ROOT", unset = file.path(getwd(), "outputs"))
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
}

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(base_dir, paste0("ABM_NaturalMissingness_", run_stamp))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "scenario_raw_results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "scenario_test_summaries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "scenario_overfitting_summaries"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "post_analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "post_analysis", "latex_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "post_analysis", "figures"), recursive = TRUE, showWarnings = FALSE)

safe_write_text(output_dir, file.path(getwd(), "latest_abm_output_dir.txt"))

cat("Working directory:", getwd(), "\n")
cat("Output directory:", output_dir, "\n")

# -----------------------------------------------------------------------------
# DATA PATH
# -----------------------------------------------------------------------------
# Option 1: Excel file
abm_file_path <- Sys.getenv("ABM_DATA_FILE", unset = "~/Desktop/project 1/Abm.xlsx")

# Option 2: SPSS file
# abm_file_path <- "abm.sav"

command_args <- commandArgs(trailingOnly = TRUE)
if (length(command_args) >= 1 && nzchar(command_args[1])) {
  abm_file_path <- command_args[1]
}

# If you already created Abm in your R session, for example:
# Abm <- read_excel("Abm.xlsx")
# then the function below will use that object automatically.

outcome_name <- "abm"
id_cols <- c("casenum")

train_fraction   <- 0.80
m                <- 5      # Use m = 5 for publication; use m = 1 for quick testing
main_seed        <- 123
n_random_splits  <- 50

missing_setting_label <- "Natural"

methods <- c("pmm", "rf", "cart", "norm", "midastouch")
method_order <- c("CCA", "Mean", "MissInd", "MissForest",
                  "PMM", "RF", "CART", "NORM", "MIDASTOUCH",
                  "MIE_SE_GLM", "MIE_SE_KNN", "MIE_SE_DT", "Stacking")
classifier_models <- c("GLM", "KNN", "DT")

cv_folds_standard <- 5
cv_folds_stacking <- 5
cv_folds_meta     <- 5

fixed_knn_grid <- data.frame(k = seq(3, 17, by = 2))

# Pre-specified predictor panel:
# This panel avoids ID, demographic/categorical, administrative, treatment,
# follow-up, and direct diagnostic-result variables. It keeps numerical blood
# and cerebrospinal-fluid measurements that are clinically relevant for
# distinguishing bacterial vs viral meningitis.
selected_abm_numeric_predictors <- c(
  "age",
  "wbc",
  "pmn",
  "bands",
  "bloodgl",
  "gl",
  "pr",
  "reds",
  "whites",
  "polys",
  "lymphs",
  "monos",
  "others"
)

# Variables with fewer than this many observed values are excluded because MICE
# can become unstable when a predictor is extremely sparse.
min_observed_values_per_predictor <- 5

# =============================================================================
# 2) DATA PREPARATION: ABM NATURAL MISSINGNESS
# =============================================================================

read_abm_data <- function(file_path) {
  if (exists("Abm", envir = .GlobalEnv)) {
    message("Using existing object: Abm")
    return(as.data.frame(get("Abm", envir = .GlobalEnv)))
  }
  
  file_path_expanded <- path.expand(file_path)
  
  if (!file.exists(file_path_expanded)) {
    stop(paste0(
      "ABM file not found at: ", file_path,
      "\nPut Abm.xlsx in your working directory or change abm_file_path."
    ))
  }
  
  ext <- tolower(tools::file_ext(file_path_expanded))
  
  if (ext %in% c("xlsx", "xls")) {
    out <- readxl::read_excel(file_path_expanded)
  } else if (ext == "sav") {
    out <- haven::read_sav(file_path_expanded)
    out <- as.data.frame(out)
  } else if (ext == "csv") {
    out <- readr::read_csv(
      file_path_expanded,
      na = c("", "NA", "NaN", ".", "?", " "),
      show_col_types = FALSE
    )
  } else {
    stop("Unsupported file type. Use .xlsx, .xls, .sav, or .csv")
  }
  
  as.data.frame(out)
}

prepare_abm_data <- function(file_path) {
  
  abm_raw <- read_abm_data(file_path)
  
  # Convert labelled SPSS columns to numeric/character where needed.
  abm_raw <- abm_raw %>%
    mutate(across(everything(), ~ {
      if (inherits(.x, "haven_labelled")) {
        return(as.numeric(.x))
      } else {
        return(.x)
      }
    }))
  
  if (!outcome_name %in% names(abm_raw)) {
    stop(paste0("Outcome column not found: ", outcome_name))
  }
  
  # Remove rows where dependent variable is missing.
  n_before <- nrow(abm_raw)
  abm_raw <- abm_raw %>%
    filter(!is.na(.data[[outcome_name]]))
  n_after_drop_outcome <- nrow(abm_raw)
  
  # Keep only valid binary outcome values 0 and 1.
  abm_raw[[outcome_name]] <- suppressWarnings(as.numeric(abm_raw[[outcome_name]]))
  abm_raw <- abm_raw %>%
    filter(.data[[outcome_name]] %in% c(0, 1))
  n_after_binary_outcome <- nrow(abm_raw)
  
  # ---------------------------------------------------------------------------
  # Pre-specified predictor selection
  # ---------------------------------------------------------------------------
  # Do NOT automatically use all numeric variables. In ABM, several variables are
  # numerically coded but may represent diagnosis, treatment, follow-up, or
  # diagnostic test-result information. Using them can create leakage.
  #
  # Therefore, we use a pre-specified numerical blood/CSF clinical predictor
  # panel defined in selected_abm_numeric_predictors.
  numeric_cols <- names(abm_raw)[sapply(abm_raw, is.numeric)]
  
  missing_selected_predictors <- setdiff(selected_abm_numeric_predictors, names(abm_raw))
  if (length(missing_selected_predictors) > 0) {
    stop(paste0(
      "These pre-specified ABM predictors were not found in the data: ",
      paste(missing_selected_predictors, collapse = ", ")
    ))
  }
  
  non_numeric_selected_predictors <- setdiff(selected_abm_numeric_predictors, numeric_cols)
  if (length(non_numeric_selected_predictors) > 0) {
    stop(paste0(
      "These pre-specified ABM predictors are not numeric after import: ",
      paste(non_numeric_selected_predictors, collapse = ", ")
    ))
  }
  
  numeric_predictors <- selected_abm_numeric_predictors
  
  # Remove selected predictors with 100% missingness or too few observed values.
  observed_counts <- sapply(abm_raw[, numeric_predictors, drop = FALSE], function(x) sum(!is.na(x)))
  too_sparse_predictors <- names(observed_counts)[observed_counts < min_observed_values_per_predictor]
  numeric_predictors <- setdiff(numeric_predictors, too_sparse_predictors)
  
  # Remove selected predictors with zero variance after dropping missing outcome.
  zero_var_predictors <- numeric_predictors[sapply(abm_raw[, numeric_predictors, drop = FALSE], function(x) {
    ux <- unique(stats::na.omit(x))
    length(ux) < 2
  })]
  numeric_predictors <- setdiff(numeric_predictors, zero_var_predictors)
  
  # Variables excluded because they are not part of the pre-specified leakage-reduced panel.
  variables_not_in_prespecified_panel <- setdiff(
    names(abm_raw),
    c(outcome_name, id_cols, numeric_predictors)
  )
  
  # Final dataset with original variable names.
  dat_original_names <- abm_raw %>%
    select(all_of(c(outcome_name, numeric_predictors))) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
  
  dat <- dat_original_names
  dat$y <- factor(ifelse(dat[[outcome_name]] == 1, "pos", "neg"), levels = c("neg", "pos"))
  dat[[outcome_name]] <- NULL
  
  dat <- dat[, c("y", numeric_predictors), drop = FALSE]
  
  variable_map <- data.frame(
    Original_Name = numeric_predictors,
    Model_Name    = paste0("x", seq_along(numeric_predictors)),
    stringsAsFactors = FALSE
  )
  
  names(dat) <- c("y", variable_map$Model_Name)
  
  excluded_variables <- data.frame(
    Variable = c(
      id_cols[id_cols %in% names(abm_raw)],
      variables_not_in_prespecified_panel,
      too_sparse_predictors,
      zero_var_predictors
    ),
    Reason = c(
      rep("ID variable excluded", length(id_cols[id_cols %in% names(abm_raw)])),
      rep(
        "Excluded a priori: not part of leakage-reduced numerical blood/CSF predictor panel",
        length(variables_not_in_prespecified_panel)
      ),
      rep(paste0("Too sparse: fewer than ", min_observed_values_per_predictor, " observed values"), length(too_sparse_predictors)),
      rep("Zero variance after removing missing outcome", length(zero_var_predictors))
    ),
    stringsAsFactors = FALSE
  ) %>%
    distinct()
  
  predictor_missing <- mean(is.na(dat[, setdiff(names(dat), "y"), drop = FALSE]))
  
  cat("ABM dataset loaded successfully.\n")
  cat("Original rows:", n_before, "\n")
  cat("Rows after removing missing outcome:", n_after_drop_outcome, "\n")
  cat("Rows after keeping abm in {0,1}:", n_after_binary_outcome, "\n")
  cat("Pre-specified numeric predictors used:", ncol(dat) - 1, "\n")
  cat("Outcome counts:\n")
  print(table(dat$y))
  
  cat("\nOutcome prevalence (%):\n")
  print(round(prop.table(table(dat$y)) * 100, 2))
  
  cat("\nOverall natural missingness in selected numeric predictors (%):\n")
  print(round(predictor_missing * 100, 2))
  
  cat("\nPre-specified numeric predictors used:\n")
  print(variable_map)
  
  if (nrow(excluded_variables) > 0) {
    cat("\nExcluded variables:\n")
    print(excluded_variables)
  }
  
  attr(dat, "variable_map") <- variable_map
  attr(dat, "excluded_variables") <- excluded_variables
  attr(dat, "natural_missingness") <- predictor_missing
  attr(dat, "n_original_rows") <- n_before
  attr(dat, "n_after_missing_outcome_removed") <- n_after_drop_outcome
  attr(dat, "n_after_binary_outcome_filter") <- n_after_binary_outcome
  
  dat
}

source_datasets <- list(
  ABM = prepare_abm_data(abm_file_path)
)

abm_variable_map <- attr(source_datasets$ABM, "variable_map")
abm_excluded_variables <- attr(source_datasets$ABM, "excluded_variables")

safe_write_csv(
  abm_variable_map,
  file.path(output_dir, "abm_variable_map.csv")
)

safe_write_csv(
  abm_excluded_variables,
  file.path(output_dir, "abm_excluded_variables.csv")
)

cat("\nPrepared source datasets:\n")
for (nm in names(source_datasets)) {
  cat(sprintf(
    "- %s: %d rows, %d predictors, prevalence = %.3f, natural missingness = %.3f\n",
    nm,
    nrow(source_datasets[[nm]]),
    ncol(source_datasets[[nm]]) - 1,
    mean(source_datasets[[nm]]$y == "pos"),
    attr(source_datasets[[nm]], "natural_missingness")
  ))
}
cat("\n")

# =============================================================================
# 3) METRIC FUNCTIONS
# =============================================================================

compute_auc <- function(true_y, pred_prob) {
  true_y <- factor(true_y, levels = c("neg", "pos"))
  pred_prob <- as.numeric(pred_prob)
  if (length(unique(na.omit(droplevels(true_y)))) < 2) return(NA_real_)
  if (all(is.na(pred_prob))) return(NA_real_)
  
  roc_obj <- tryCatch(
    pROC::roc(true_y, pred_prob, levels = c("neg", "pos"), quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(roc_obj)) return(NA_real_)
  
  as.numeric(pROC::auc(roc_obj))
}

brier_score <- function(true_y, pred_prob) {
  true_y <- factor(true_y, levels = c("neg", "pos"))
  pred_prob <- as.numeric(pred_prob)
  if (length(unique(na.omit(droplevels(true_y)))) < 2) return(NA_real_)
  if (all(is.na(pred_prob))) return(NA_real_)
  
  true_numeric <- as.numeric(true_y == "pos")
  mean((true_numeric - pred_prob)^2, na.rm = TRUE)
}

precision_score <- function(true_y, pred_class) {
  true_y    <- factor(true_y,    levels = c("neg", "pos"))
  pred_class <- factor(pred_class, levels = c("neg", "pos"))
  
  if (length(unique(na.omit(droplevels(true_y)))) < 2) return(NA_real_)
  
  cm <- tryCatch(
    caret::confusionMatrix(pred_class, true_y, positive = "pos"),
    error = function(e) NULL
  )
  if (is.null(cm)) return(NA_real_)
  
  as.numeric(cm$byClass["Precision"])
}

confusion_metrics <- function(true_y, pred_class) {
  true_y <- factor(true_y, levels = c("neg", "pos"))
  pred_class <- factor(pred_class, levels = c("neg", "pos"))
  
  if (length(unique(na.omit(droplevels(true_y)))) < 2) {
    return(list(
      Sensitivity = NA_real_,
      Specificity = NA_real_,
      F1          = NA_real_
    ))
  }
  
  cm <- tryCatch(
    caret::confusionMatrix(pred_class, true_y, positive = "pos"),
    error = function(e) NULL
  )
  
  if (is.null(cm)) {
    return(list(
      Sensitivity = NA_real_,
      Specificity = NA_real_,
      F1          = NA_real_
    ))
  }
  
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
  if (length(unique(obs_y_num)) < 2 || length(pred_prob) < 5) {
    return(list(cal_slope = NA_real_, cal_intercept = NA_real_))
  }
  logit_pred <- log(pred_prob / (1 - pred_prob))
  cal_fit <- tryCatch(
    glm(obs_y_num ~ logit_pred, family = binomial()),
    error = function(e) NULL
  )
  list(
    cal_intercept = if (!is.null(cal_fit)) as.numeric(coef(cal_fit)[1]) else NA_real_,
    cal_slope     = if (!is.null(cal_fit)) as.numeric(coef(cal_fit)[2]) else NA_real_
  )
}

find_optimal_threshold <- function(obs, pred_prob) {
  obs <- factor(obs, levels = c("neg", "pos"))
  pred_prob <- as.numeric(pred_prob)
  
  if (length(unique(obs)) < 2 || all(is.na(pred_prob))) {
    return(0.50)
  }
  
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
  
  out <- as.numeric(threshold_results$threshold[1])
  if (!is.finite(out)) out <- 0.50
  out
}

# =============================================================================
# 4) MISSINGNESS, WINSORIZATION, AND IMPUTATION HELPERS
# =============================================================================

winsorize_train_data_with_missing <- function(train_x) {
  numeric_cols <- names(train_x)[sapply(train_x, is.numeric)]
  
  limits_list <- list()
  
  for (col in numeric_cols) {
    x <- train_x[[col]]
    x_nonmiss <- x[is.finite(x) & !is.na(x)]
    
    if (length(x_nonmiss) < 5) {
      limits_list[[col]] <- list(low = -Inf, high = Inf)
    } else {
      stats <- boxplot.stats(x_nonmiss)$stats
      limits_list[[col]] <- list(low = stats[1], high = stats[5])
    }
  }
  
  limits_list
}

apply_winsorization <- function(data_x, limits_list) {
  for (col in names(limits_list)) {
    if (col %in% names(data_x)) {
      limits <- limits_list[[col]]
      data_x[[col]] <- ifelse(
        is.na(data_x[[col]]),
        NA_real_,
        pmin(pmax(data_x[[col]], limits$low), limits$high)
      )
    }
  }
  data_x
}

summarize_missingness <- function(df_x, dataset_label, scenario_name, split_id, missing_setting) {
  var_tbl <- data.frame(
    Scenario                = scenario_name,
    Split                   = split_id,
    Dataset                 = dataset_label,
    Missing_Setting         = missing_setting,
    Variable                = names(df_x),
    Missing_Count           = sapply(df_x, function(x) sum(is.na(x))),
    Total_Count             = nrow(df_x),
    Missing_Percent         = sapply(df_x, function(x) mean(is.na(x))),
    stringsAsFactors        = FALSE
  )
  
  overall_tbl <- data.frame(
    Scenario                 = scenario_name,
    Split                    = split_id,
    Dataset                  = dataset_label,
    Missing_Setting          = missing_setting,
    Total_Cells              = nrow(df_x) * ncol(df_x),
    Missing_Cells            = sum(is.na(as.matrix(df_x))),
    Overall_Missing_Percent  = mean(is.na(as.matrix(df_x))),
    stringsAsFactors         = FALSE
  )
  
  list(
    variable_level = var_tbl,
    overall_level  = overall_tbl
  )
}

create_missingness_comparison_table <- function(overall_missingness_tbl) {
  overall_missingness_tbl %>%
    arrange(Scenario, Split, Missing_Setting, Dataset)
}

# Train-only MICE imputation matched to Code 2.
# Train set is imputed using MICE, and the test set is imputed using the
# fitted training-imputation object through mice.mids(newdata = ...).
# Outcome is never included in train_x or test_x.
impute_train_test_trainonly <- function(train_x, test_x, method_name, m = 5, seed = 123, maxit = 10) {
  train_x <- as.data.frame(train_x)
  test_x  <- as.data.frame(test_x)
  
  n_train <- nrow(train_x)
  n_test  <- nrow(test_x)
  
  if (n_train == 0 || n_test == 0) {
    stop("Both train_x and test_x must have positive number of rows.")
  }
  
  set.seed(seed)
  
  mids_train <- mice::mice(
    data      = train_x,
    method    = method_name,
    m         = m,
    maxit     = maxit,
    printFlag = FALSE
  )
  
  train_list <- vector("list", m)
  test_list  <- vector("list", m)
  
  for (imp_idx in seq_len(m)) {
    
    train_completed <- mice::complete(mids_train, action = imp_idx)
    train_completed <- as.data.frame(train_completed)
    train_completed <- train_completed[, names(train_x), drop = FALSE]
    
    combined_data <- dplyr::bind_rows(train_completed, test_x)
    
    set.seed(seed + imp_idx)
    
    mids_test <- mice::mice.mids(
      obj       = mids_train,
      newdata   = combined_data,
      maxit     = 1,
      printFlag = FALSE
    )
    
    completed_full <- mice::complete(mids_test, action = imp_idx)
    completed_full <- as.data.frame(completed_full)
    completed_full <- completed_full[, names(train_x), drop = FALSE]
    
    # Keep the already completed training data fixed.
    completed_full[seq_len(n_train), ] <- train_completed
    
    train_list[[imp_idx]] <- completed_full[seq_len(n_train), , drop = FALSE]
    test_list[[imp_idx]]  <- completed_full[(n_train + 1):(n_train + n_test), , drop = FALSE]
  }
  
  list(train = train_list, test = test_list)
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
  
  train_df$y <- factor(train_df$y, levels = c("neg", "pos"))
  
  # Safety: remove predictors with zero variance in this training split.
  predictor_cols <- setdiff(names(train_df), "y")
  zero_var_cols <- predictor_cols[sapply(train_df[, predictor_cols, drop = FALSE], function(x) {
    ux <- unique(stats::na.omit(x))
    length(ux) < 2
  })]
  
  if (length(zero_var_cols) > 0) {
    train_df <- train_df[, setdiff(names(train_df), zero_var_cols), drop = FALSE]
    test_df  <- test_df[,  setdiff(names(test_df),  zero_var_cols), drop = FALSE]
  }
  
  train_control <- caret::trainControl(
    method          = "cv",
    number          = cv_folds,
    classProbs      = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )
  
  set.seed(seed)
  
  model_fit <- tryCatch({
    if (model_type == "GLM") {
      suppressWarnings(
        caret::train(
          y ~ .,
          data      = train_df,
          method    = "glm",
          family    = "binomial",
          trControl = train_control,
          metric    = "ROC"
        )
      )
    } else if (model_type == "KNN") {
      suppressWarnings(
        caret::train(
          y ~ .,
          data       = train_df,
          method     = "knn",
          trControl  = train_control,
          metric     = "ROC",
          preProcess = c("center", "scale"),
          tuneGrid   = get_knn_grid()
        )
      )
    } else if (model_type == "DT") {
      suppressWarnings(
        caret::train(
          y ~ .,
          data      = train_df,
          method    = "rpart",
          trControl = train_control,
          metric    = "ROC",
          tuneGrid  = get_dt_grid()
        )
      )
    } else {
      stop("Unsupported model_type")
    }
  }, error = function(e) {
    message("Model failed: ", model_type, ". Error: ", e$message)
    NULL
  })
  
  if (is.null(model_fit)) {
    base_rate <- mean(train_df$y == "pos")
    cv_pred <- data.frame(
      rowIndex = seq_len(nrow(train_df)),
      obs      = train_df$y,
      pos      = base_rate
    )
    test_prob <- rep(base_rate, nrow(test_df))
    
    return(list(
      model     = NULL,
      cv_pred   = cv_pred,
      test_prob = test_prob
    ))
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
  
  test_prob <- tryCatch(
    predict(model_fit, newdata = test_df, type = "prob")$pos,
    error = function(e) rep(mean(train_df$y == "pos"), nrow(test_df))
  )
  
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
# 7) MAIN EXPERIMENT LOOP: NATURAL MISSINGNESS ONLY
# =============================================================================

all_results               <- list()
all_test_tables           <- list()
all_missingness_variable  <- list()
all_missingness_overall   <- list()
all_timing                <- list()    # <<< NEW: per-split timing

for (scenario_name in names(source_datasets)) {
  
  dat <- source_datasets[[scenario_name]]
  
  scenario_missingness_variable <- data.frame()
  scenario_missingness_overall  <- data.frame()
  scenario_timing <- data.frame()
  
  # ---- Checkpoint / resume --------------------------------------------------
  checkpoint_file <- file.path(output_dir, paste0("checkpoint_", scenario_name, "_raw_results.csv"))
  if (file.exists(checkpoint_file)) {
    scenario_results <- readr::read_csv(checkpoint_file, show_col_types = FALSE)
    scenario_results$Imputation_Method <- as.character(scenario_results$Imputation_Method)
    completed_keys <- unique(paste(scenario_results$Split, sep = "_"))
    cat("  Checkpoint loaded:", nrow(scenario_results), "rows,",
        length(completed_keys), "splits already done.\n")
  } else {
    scenario_results <- data.frame()
    completed_keys   <- character(0)
  }
  # ---------------------------------------------------------------------------
  
  cat("============================================================\n")
  cat("Running scenario:", scenario_name, "\n")
  cat("Missingness setting:", missing_setting_label, "\n")
  cat("============================================================\n")
  
  for (split_id in seq_len(n_random_splits)) {
    
    cat("  Random split:", split_id, "of", n_random_splits, "\n")
    
    # Skip if already completed
    split_key <- as.character(split_id)
    if (split_key %in% completed_keys) {
      cat("  Split", split_id, "- skipping (checkpoint)\n")
      next
    }
    n_before_split <- nrow(scenario_results)
    
    data_seed <- main_seed + split_id * 10000
    
    split_start_time <- proc.time()
    
    set.seed(data_seed + 1)
    index <- caret::createDataPartition(dat$y, p = train_fraction, list = FALSE)
    
    train <- dat[index, , drop = FALSE]
    test  <- dat[-index, , drop = FALSE]
    
    train_y <- factor(train$y, levels = c("neg", "pos"))
    test_y  <- factor(test$y,  levels = c("neg", "pos"))
    
    train_x_raw <- train[, setdiff(names(train), "y"), drop = FALSE]
    test_x_raw  <- test[,  setdiff(names(test),  "y"), drop = FALSE]
    
    # Winsorization is applied using non-missing training values only.
    # Missing values remain missing and are imputed later.
    winsor_limits <- winsorize_train_data_with_missing(train_x_raw)
    train_x <- apply_winsorization(train_x_raw, winsor_limits)
    test_x  <- apply_winsorization(test_x_raw,  winsor_limits)
    
    train_missing_summary <- summarize_missingness(
      df_x            = train_x,
      dataset_label   = "Train",
      scenario_name   = scenario_name,
      split_id        = split_id,
      missing_setting = missing_setting_label
    )
    
    test_missing_summary <- summarize_missingness(
      df_x            = test_x,
      dataset_label   = "Test",
      scenario_name   = scenario_name,
      split_id        = split_id,
      missing_setting = missing_setting_label
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
    
    # Impute using each individual method.
    imputed_objs <- list()
    
    for (method in methods) {
      cat("    Imputation method:", method, "\n")
      
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
    
    # Individual imputation methods.
    for (method in methods) {
      for (clf in classifier_models) {
        
        cat("    Model:", clf, "| Method:", toupper(method), "\n")
        
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
          Missing_Setting    = missing_setting_label,
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
    
    # ---- Simple baselines: CCA, Mean, MissInd, MissForest (R1-C6 / R2-C1) ---
    
    # Helper to build one result row (ABM version)
    store_abm_single <- function(method_name, clf, cv_pred_obj, test_prob_vec) {
      thresh    <- find_optimal_threshold(cv_pred_obj$obs, cv_pred_obj$pos)
      tr_prob   <- cv_pred_obj$pos
      tr_cl     <- factor(ifelse(tr_prob       >= thresh, "pos", "neg"), levels = c("neg","pos"))
      te_cl     <- factor(ifelse(test_prob_vec >= thresh, "pos", "neg"), levels = c("neg","pos"))
      tr_cm     <- confusion_metrics(cv_pred_obj$obs, tr_cl)
      te_cm     <- confusion_metrics(test_y,          te_cl)
      cal       <- compute_calibration(test_y, test_prob_vec)
      data.frame(
        Scenario = scenario_name, Split = split_id,
        Missing_Setting = missing_setting_label,
        Imputation_Method = method_name, Model = clf,
        Train_AUC      = compute_auc(cv_pred_obj$obs, tr_prob),
        Test_AUC       = compute_auc(test_y, test_prob_vec),
        Train_Accuracy = mean(tr_cl == cv_pred_obj$obs),
        Test_Accuracy  = mean(te_cl == test_y),
        Train_Brier    = brier_score(cv_pred_obj$obs, tr_prob),
        Test_Brier     = brier_score(test_y, test_prob_vec),
        Train_Precision   = precision_score(cv_pred_obj$obs, tr_cl),
        Test_Precision    = precision_score(test_y, te_cl),
        Train_Sensitivity = tr_cm$Sensitivity, Test_Sensitivity = te_cm$Sensitivity,
        Train_Specificity = tr_cm$Specificity, Test_Specificity = te_cm$Specificity,
        Train_F1 = tr_cm$F1, Test_F1 = te_cm$F1, Threshold = thresh,
        Cal_Slope = cal$cal_slope, Cal_Intercept = cal$cal_intercept
      )
    }
    
    # 1) CCA
    cat("    Method: CCA\n")
    cca_complete <- complete.cases(train_x)
    if (sum(cca_complete) >= 10) {
      train_x_cca <- train_x[cca_complete, , drop = FALSE]
      train_y_cca <- train_y[cca_complete]
      tmeans_cca  <- colMeans(train_x, na.rm = TRUE)
      test_x_cca  <- as.data.frame(lapply(names(test_x), function(col) {
        v <- test_x[[col]]; v[is.na(v)] <- tmeans_cca[[col]]; v
      })); names(test_x_cca) <- names(test_x)
      for (clf in classifier_models) {
        fit_cca <- fit_model_with_cv_predictions(
          data.frame(train_x_cca, y = train_y_cca), test_x_cca, clf,
          seed = data_seed + 30000 + match(clf, classifier_models) * 10,
          cv_folds = cv_folds_standard)
        scenario_results <- bind_rows(scenario_results,
                                      store_abm_single("CCA", clf, fit_cca$cv_pred, fit_cca$test_prob))
      }
    }
    
    # 2) Mean imputation
    cat("    Method: Mean Imputation\n")
    tmeans_m     <- colMeans(train_x, na.rm = TRUE)
    train_x_mean <- train_x; test_x_mean <- test_x
    for (col in names(train_x_mean)) {
      train_x_mean[[col]][is.na(train_x_mean[[col]])] <- tmeans_m[[col]]
      test_x_mean[[col]][is.na(test_x_mean[[col]])]   <- tmeans_m[[col]]
    }
    for (clf in classifier_models) {
      fit_m <- fit_model_with_cv_predictions(
        data.frame(train_x_mean, y = train_y), test_x_mean, clf,
        seed = data_seed + 31000 + match(clf, classifier_models) * 10,
        cv_folds = cv_folds_standard)
      scenario_results <- bind_rows(scenario_results,
                                    store_abm_single("Mean", clf, fit_m$cv_pred, fit_m$test_prob))
    }
    
    # 3) Missing indicator ─ binary flag per missing predictor + mean-fill (Sterne et al. 2009)
    cat("    Method: Missing Indicator\n")
    hmiss        <- names(train_x)[sapply(train_x, function(x) any(is.na(x)))]
    train_x_ind  <- train_x; test_x_ind <- test_x
    for (col in hmiss) {
      train_x_ind[[paste0(col, "_miss")]] <- as.numeric(is.na(train_x[[col]]))
      test_x_ind[[paste0(col,  "_miss")]] <- as.numeric(is.na(test_x[[col]]))
      train_x_ind[[col]][is.na(train_x_ind[[col]])] <- tmeans_m[[col]]
      test_x_ind[[col]][is.na(test_x_ind[[col]])]   <- tmeans_m[[col]]
    }
    for (clf in classifier_models) {
      fit_ind <- fit_model_with_cv_predictions(
        data.frame(train_x_ind, y = train_y), test_x_ind, clf,
        seed = data_seed + 32000 + match(clf, classifier_models) * 10,
        cv_folds = cv_folds_standard)
      scenario_results <- bind_rows(scenario_results,
                                    store_abm_single("MissInd", clf, fit_ind$cv_pred, fit_ind$test_prob))
    }
    
    # 4) MissForest
    cat("    Method: MissForest\n")
    n_tr_mf      <- nrow(train_x)
    comb_x_mf    <- rbind(as.data.frame(train_x), as.data.frame(test_x))
    set.seed(data_seed + 33000)
    mf_out_abm <- tryCatch(
      missForest::missForest(comb_x_mf, verbose = FALSE),
      error = function(e) { message("MissForest failed: ", e$message); NULL }
    )
    if (!is.null(mf_out_abm)) {
      cimp_abm    <- as.data.frame(mf_out_abm$ximp)
      train_x_mf  <- cimp_abm[seq_len(n_tr_mf), , drop = FALSE]
      test_x_mf   <- cimp_abm[(n_tr_mf + 1):nrow(cimp_abm), , drop = FALSE]
      for (clf in classifier_models) {
        fit_mf <- fit_model_with_cv_predictions(
          data.frame(train_x_mf, y = train_y), test_x_mf, clf,
          seed = data_seed + 33000 + match(clf, classifier_models) * 10,
          cv_folds = cv_folds_standard)
        scenario_results <- bind_rows(scenario_results,
                                      store_abm_single("MissForest", clf, fit_mf$cv_pred, fit_mf$test_prob))
      }
    }
    
    # Imputation-based stacking.
    for (clf in classifier_models) {
      
      cat("    Model:", clf, "| Method: Stacking\n")
      
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
        
        combined_oof <- Reduce(
          function(x, y) full_join(x, y, by = "row_id"),
          oof_preds_list
        ) %>%
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
        Missing_Setting    = missing_setting_label,
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
    cat("    Model: ALL | Method: MIE_SE (PMM x GLM+KNN+DT, per-clf meta-learner)\n")
    
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
        Missing_Setting    = missing_setting_label,
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
    
    # ---- Record per-split computational cost --------------------------------
    split_elapsed <- (proc.time() - split_start_time)[["elapsed"]]
    scenario_timing <- bind_rows(scenario_timing, data.frame(
      Scenario        = scenario_name,
      Split           = split_id,
      Missing_Setting = missing_setting_label,
      Elapsed_Seconds = split_elapsed
    ))
    
    # ---- Checkpoint: append this split's rows to disk -----------------------
    new_rows <- if (nrow(scenario_results) > n_before_split)
      scenario_results[(n_before_split + 1):nrow(scenario_results), ] else data.frame()
    if (nrow(new_rows) > 0) {
      readr::write_csv(new_rows, checkpoint_file, append = file.exists(checkpoint_file))
      completed_keys <- c(completed_keys, split_key)
      cat("  Split", split_id, "checkpointed (", nrow(new_rows), "rows )\n")
    }
    # -------------------------------------------------------------------------
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
    arrange(Split, Model, Imputation_Method)
  
  safe_write_csv(
    scenario_results,
    file.path(output_dir, "scenario_raw_results", paste0(scenario_name, "_raw_results.csv"))
  )
  
  # ---- Save per-split timing for this scenario ----------------------------
  timing_summary_scenario <- scenario_timing %>%
    summarise(
      Scenario             = scenario_name,
      Missing_Setting      = missing_setting_label,
      Mean_Elapsed_Seconds = mean(Elapsed_Seconds, na.rm = TRUE),
      SD_Elapsed_Seconds   = sd(Elapsed_Seconds, na.rm = TRUE),
      Min_Elapsed_Seconds  = min(Elapsed_Seconds, na.rm = TRUE),
      Max_Elapsed_Seconds  = max(Elapsed_Seconds, na.rm = TRUE),
      n_splits             = n()
    )
  safe_write_csv(
    scenario_timing,
    file.path(output_dir, paste0(scenario_name, "_per_split_timing.csv"))
  )
  safe_write_csv(
    timing_summary_scenario,
    file.path(output_dir, paste0(scenario_name, "_timing_summary.csv"))
  )
  
  all_timing[[scenario_name]] <- scenario_timing
  
  scenario_test_table <- scenario_results %>%
    select(
      Scenario, Split, Missing_Setting, Imputation_Method, Model,
      Test_AUC, Test_Accuracy, Test_Brier, Test_Precision,
      Test_Sensitivity, Test_Specificity, Test_F1, Threshold,
      Cal_Slope, Cal_Intercept
    ) %>%
    arrange(Split, Model, Imputation_Method)
  
  safe_write_csv(
    scenario_test_table,
    file.path(output_dir, "scenario_test_summaries", paste0(scenario_name, "_test_metrics.csv"))
  )
  
  scenario_overfit_table <- scenario_results %>%
    select(
      Scenario, Split, Missing_Setting, Imputation_Method, Model,
      Overfit_Gap_AUC, Overfit_Gap_Accuracy, Overfit_Gap_Brier,
      Overfit_Gap_Precision, Overfit_Gap_Sensitivity,
      Overfit_Gap_Specificity, Overfit_Gap_F1
    ) %>%
    arrange(Split, Model, Imputation_Method)
  
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

# =============================================================================
# 8) COMBINE OUTPUTS
# =============================================================================

raw_results_file        <- file.path(output_dir, "all_scenarios_raw_results.csv")
test_metrics_file       <- file.path(output_dir, "all_scenarios_test_metrics.csv")
missing_var_file        <- file.path(output_dir, "all_scenarios_actual_missingness_variable_level.csv")
missing_overall_file    <- file.path(output_dir, "all_scenarios_actual_missingness_overall.csv")
missing_comparison_file <- file.path(output_dir, "all_scenarios_actual_missingness_comparison.csv")

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

# ---- Save combined timing across all scenarios --------------------------
combined_timing <- bind_rows(all_timing)
safe_write_csv(combined_timing, file.path(output_dir, "all_scenarios_per_split_timing.csv"))

# =============================================================================
# POST-ANALYSIS / SUMMARY / RANKING / FIGURES
# =============================================================================

# =============================================================================
# 9) MISSINGNESS SUMMARY ACROSS RANDOM SPLITS
# =============================================================================

missingness_summary <- combined_missingness_overall %>%
  group_by(Scenario, Missing_Setting, Dataset) %>%
  summarise(
    Mean_Overall_Missingness = mean(Overall_Missing_Percent, na.rm = TRUE),
    SD_Overall_Missingness   = sd(Overall_Missing_Percent, na.rm = TRUE),
    Mean_Missing_Cells       = mean(Missing_Cells, na.rm = TRUE),
    SD_Missing_Cells         = sd(Missing_Cells, na.rm = TRUE),
    Mean_Total_Cells         = mean(Total_Cells, na.rm = TRUE),
    SD_Total_Cells           = sd(Total_Cells, na.rm = TRUE),
    n_splits                 = n(),
    .groups = "drop"
  ) %>%
  arrange(Scenario, Missing_Setting, Dataset)

missingness_variable_summary <- combined_missingness_variable %>%
  group_by(Scenario, Missing_Setting, Dataset, Variable) %>%
  summarise(
    Mean_Missing_Percent = mean(Missing_Percent, na.rm = TRUE),
    SD_Missing_Percent   = sd(Missing_Percent, na.rm = TRUE),
    Mean_Missing_Count   = mean(Missing_Count, na.rm = TRUE),
    SD_Missing_Count     = sd(Missing_Count, na.rm = TRUE),
    n_splits             = n(),
    .groups = "drop"
  ) %>%
  arrange(Scenario, Missing_Setting, Dataset, Variable)

safe_write_csv(
  missingness_summary,
  file.path(output_dir, "post_analysis", "missingness_mean_sd_summary.csv")
)

safe_write_csv(
  missingness_variable_summary,
  file.path(output_dir, "post_analysis", "missingness_variable_mean_sd_summary.csv")
)

missingness_train_test_table <- combined_missingness_overall %>%
  group_by(Scenario, Missing_Setting, Dataset) %>%
  summarise(
    Mean_Missingness = mean(Overall_Missing_Percent, na.rm = TRUE),
    SD_Missingness   = sd(Overall_Missing_Percent, na.rm = TRUE),
    n_splits         = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Missingness_Mean_SD = sprintf("%.3f ± %.3f", Mean_Missingness, SD_Missingness)
  ) %>%
  select(Scenario, Missing_Setting, Dataset, Missingness_Mean_SD) %>%
  pivot_wider(
    names_from  = Dataset,
    values_from = Missingness_Mean_SD
  ) %>%
  arrange(Scenario, Missing_Setting)

safe_write_csv(
  missingness_train_test_table,
  file.path(output_dir, "post_analysis", "missingness_train_test_mean_sd_table.csv")
)

missingness_train_test_numeric <- combined_missingness_overall %>%
  group_by(Scenario, Missing_Setting, Dataset) %>%
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
  arrange(Scenario, Missing_Setting)

safe_write_csv(
  missingness_train_test_numeric,
  file.path(output_dir, "post_analysis", "missingness_train_test_mean_sd_numeric.csv")
)

# =============================================================================
# 10) OVERFITTING SUMMARY
# =============================================================================

overfitting_summary <- combined_results %>%
  mutate(
    Method = as.character(Imputation_Method),
    Model  = as.character(Model)
  ) %>%
  group_by(Scenario, Model, Missing_Setting, Method) %>%
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
  arrange(Scenario, Model, Missing_Setting, Method)

safe_write_csv(
  as.data.frame(overfitting_summary),
  file.path(output_dir, "post_analysis", "overfitting_mean_sd_summary.csv")
)

# =============================================================================
# 11) MAIN PERFORMANCE SUMMARY
# =============================================================================

missingness_join_tbl <- missingness_summary %>%
  select(
    Scenario,
    Missing_Setting,
    Dataset,
    Mean_Overall_Missingness,
    SD_Overall_Missingness
  ) %>%
  pivot_wider(
    names_from = Dataset,
    values_from = c(
      Mean_Overall_Missingness,
      SD_Overall_Missingness
    ),
    names_sep = "_"
  )

final_results <- combined_results %>%
  mutate(
    Method = as.character(Imputation_Method),
    Model  = as.character(Model)
  ) %>%
  group_by(Scenario, Model, Missing_Setting, Method) %>%
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
    by = c("Scenario", "Model", "Missing_Setting", "Method")
  ) %>%
  left_join(
    missingness_join_tbl,
    by = c("Scenario", "Missing_Setting")
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
  levels = c("CCA", "Mean", "MissInd", "MissForest",
             "PMM", "RF", "CART", "NORM", "MIDASTOUCH",
             "MIE_SE_GLM", "MIE_SE_KNN", "MIE_SE_DT", "Stacking")
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
# 12) 95% CONFIDENCE INTERVALS FOR TEST METRICS
# =============================================================================

metric_info <- data.frame(
  Metric = c("AUC", "Accuracy", "Brier", "Precision", "Sensitivity", "Specificity", "F1"),
  Column = c("Test_AUC", "Test_Accuracy", "Test_Brier", "Test_Precision",
             "Test_Sensitivity", "Test_Specificity", "Test_F1"),
  stringsAsFactors = FALSE
)

make_ci_table <- function(results, metric_name, metric_col) {
  results %>%
    mutate(
      Method = as.character(Imputation_Method),
      Model  = as.character(Model)
    ) %>%
    group_by(Scenario, Model, Missing_Setting, Method) %>%
    summarise(
      Metric = metric_name,
      n      = sum(!is.na(.data[[metric_col]])),
      Mean   = mean(.data[[metric_col]], na.rm = TRUE),
      SD     = sd(.data[[metric_col]], na.rm = TRUE),
      SE     = SD / sqrt(n),
      CI_Lower = ifelse(n > 1, Mean - qt(0.975, df = n - 1) * SE, NA_real_),
      CI_Upper = ifelse(n > 1, Mean + qt(0.975, df = n - 1) * SE, NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      Mean_SD = sprintf("%.4f ± %.4f", Mean, SD),
      Mean_CI = sprintf("%.4f (95%% CI: %.4f, %.4f)", Mean, CI_Lower, CI_Upper)
    )
}

ci_95_results <- bind_rows(lapply(seq_len(nrow(metric_info)), function(i) {
  make_ci_table(
    results     = combined_results,
    metric_name = metric_info$Metric[i],
    metric_col  = metric_info$Column[i]
  )
}))

safe_write_csv(
  ci_95_results,
  file.path(output_dir, "post_analysis", "test_metrics_95CI.csv")
)

# =============================================================================
# 13) PAIRED WILCOXON TEST: STACKING VS EACH INDIVIDUAL METHOD
# =============================================================================

run_wilcox_vs_stacking <- function(results, metric_name, metric_col) {
  
  out_list <- list()
  counter <- 1
  
  for (scenario_name in unique(as.character(results$Scenario))) {
    for (model_name in unique(as.character(results$Model))) {
      
      dat_wide <- results %>%
        filter(
          as.character(Scenario) == scenario_name,
          as.character(Model) == model_name
        ) %>%
        mutate(Method = as.character(Imputation_Method)) %>%
        select(Split, Method, value = all_of(metric_col)) %>%
        pivot_wider(
          names_from  = Method,
          values_from = value
        )
      
      if (!"Stacking" %in% names(dat_wide)) next
      
      compare_methods <- setdiff(names(dat_wide), c("Split", "Stacking"))
      
      for (method_name in compare_methods) {
        
        paired_dat <- dat_wide %>%
          select(Stacking, all_of(method_name)) %>%
          filter(stats::complete.cases(.))
        
        if (nrow(paired_dat) < 3) {
          p_value <- NA_real_
          statistic <- NA_real_
        } else {
          test_obj <- tryCatch(
            stats::wilcox.test(
              paired_dat$Stacking,
              paired_dat[[method_name]],
              paired = TRUE,
              exact  = FALSE
            ),
            error = function(e) NULL
          )
          
          p_value <- if (is.null(test_obj)) NA_real_ else test_obj$p.value
          statistic <- if (is.null(test_obj)) NA_real_ else as.numeric(test_obj$statistic)
        }
        
        mean_stacking <- mean(paired_dat$Stacking, na.rm = TRUE)
        mean_method   <- mean(paired_dat[[method_name]], na.rm = TRUE)
        median_diff   <- median(paired_dat$Stacking - paired_dat[[method_name]], na.rm = TRUE)
        
        if (metric_name == "Brier") {
          better_by_mean <- ifelse(mean_stacking < mean_method, "Stacking", method_name)
          direction_note <- "For Brier, lower is better; negative median difference favors Stacking."
        } else {
          better_by_mean <- ifelse(mean_stacking > mean_method, "Stacking", method_name)
          direction_note <- "For this metric, higher is better; positive median difference favors Stacking."
        }
        
        out_list[[counter]] <- data.frame(
          Scenario = scenario_name,
          Model = model_name,
          Missing_Setting = missing_setting_label,
          Metric = metric_name,
          Comparison = paste0("Stacking vs ", method_name),
          Method_Compared = method_name,
          n_pairs = nrow(paired_dat),
          Mean_Stacking = mean_stacking,
          Mean_Method = mean_method,
          Median_Diff_Stacking_minus_Method = median_diff,
          Better_By_Mean = better_by_mean,
          Wilcox_Statistic = statistic,
          p_value = p_value,
          Direction_Note = direction_note,
          stringsAsFactors = FALSE
        )
        
        counter <- counter + 1
      }
    }
  }
  
  bind_rows(out_list)
}

wilcox_vs_stacking <- bind_rows(lapply(seq_len(nrow(metric_info)), function(i) {
  run_wilcox_vs_stacking(
    results     = combined_results,
    metric_name = metric_info$Metric[i],
    metric_col  = metric_info$Column[i]
  )
})) %>%
  group_by(Scenario, Model, Metric) %>%
  mutate(
    p_adjust_holm = p.adjust(p_value, method = "holm")
  ) %>%
  ungroup() %>%
  arrange(Scenario, Model, Metric, p_adjust_holm)

safe_write_csv(
  wilcox_vs_stacking,
  file.path(output_dir, "post_analysis", "wilcoxon_stacking_vs_methods.csv")
)

# =============================================================================
# 14) RANKING ANALYSIS BY MODEL
# =============================================================================

create_rank_table <- function(dt, metric, higher_better = TRUE) {
  dt_copy <- data.table::copy(dt)
  dt_copy[[metric]] <- as.numeric(dt_copy[[metric]])
  
  if (higher_better) {
    dt_copy[, Overall_Rank := rank(-get(metric), ties.method = "min")]
  } else {
    dt_copy[, Overall_Rank := rank(get(metric), ties.method = "min")]
  }
  
  out <- dt_copy[, .(Method, Overall_Rank)]
  out[, W := NA_real_]
  out[, p_value := NA_real_]
  out
}

rank_results <- final_results %>%
  group_by(Model, Method) %>%
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

data.table::setDT(rank_glm_results)
data.table::setDT(rank_knn_results)
data.table::setDT(rank_dt_results)

extract_metric_table <- function(rank_dt, metric_name, higher_better = TRUE) {
  res <- create_rank_table(rank_dt, metric_name, higher_better)
  out <- as.data.frame(res)
  out$Metric <- metric_name
  out[, c("Metric", "Method", "Overall_Rank", "W", "p_value")]
}

make_combined_rank_table <- function(rank_dt) {
  rbind(
    extract_metric_table(rank_dt, "AUC",         TRUE),
    extract_metric_table(rank_dt, "Accuracy",    TRUE),
    extract_metric_table(rank_dt, "Brier",       FALSE),
    extract_metric_table(rank_dt, "Precision",   TRUE),
    extract_metric_table(rank_dt, "Sensitivity", TRUE),
    extract_metric_table(rank_dt, "Specificity", TRUE),
    extract_metric_table(rank_dt, "F1",          TRUE)
  )
}

combined_rank_table_glm <- make_combined_rank_table(rank_glm_results)
combined_rank_table_knn <- make_combined_rank_table(rank_knn_results)
combined_rank_table_dt  <- make_combined_rank_table(rank_dt_results)

safe_write_csv(combined_rank_table_glm, file.path(output_dir, "post_analysis", "GLM_combined_rank_table.csv"))
safe_write_csv(combined_rank_table_knn, file.path(output_dir, "post_analysis", "KNN_combined_rank_table.csv"))
safe_write_csv(combined_rank_table_dt,  file.path(output_dir, "post_analysis", "DT_combined_rank_table.csv"))

make_paper_rank_table <- function(combined_df) {
  # Rename any Rank_0.x columns -> "x0%" (handles MCAR/MAR style; ABM has none)
  rank_col_pattern <- grep("^Rank_", names(combined_df), value = TRUE)
  if (length(rank_col_pattern) > 0) {
    new_names <- sapply(rank_col_pattern, function(m) {
      pct <- as.numeric(sub("Rank_", "", m)) * 100
      paste0(round(pct), "%")
    })
    names(combined_df)[match(rank_col_pattern, names(combined_df))] <- new_names
  }
  
  # Sort by Metric then Overall_Rank
  combined_df <- combined_df[order(combined_df$Metric, combined_df$Overall_Rank), ]
  
  # Move Metric to first column
  other_cols <- setdiff(names(combined_df), "Metric")
  combined_df <- combined_df[, c("Metric", other_cols)]
  
  # Show W and p_value only on first row of each metric group
  if ("p_value" %in% names(combined_df)) {
    combined_df$p_value <- ifelse(!is.na(combined_df$p_value) & combined_df$p_value < 0.001,
                                  "<0.001", as.character(combined_df$p_value))
    first_row <- !duplicated(combined_df$Metric)
    if ("W" %in% names(combined_df)) combined_df$W[!first_row] <- NA
    combined_df$p_value[!first_row] <- NA
  }
  
  combined_df
}

safe_write_csv(make_paper_rank_table(as.data.frame(combined_rank_table_glm)),
               file.path(output_dir, "post_analysis", "GLM_paper_rank_table.csv"))
safe_write_csv(make_paper_rank_table(as.data.frame(combined_rank_table_knn)),
               file.path(output_dir, "post_analysis", "KNN_paper_rank_table.csv"))
safe_write_csv(make_paper_rank_table(as.data.frame(combined_rank_table_dt)),
               file.path(output_dir, "post_analysis", "DT_paper_rank_table.csv"))

# =============================================================================
# 15) SAVE MEAN ± SD TABLES
# =============================================================================

fmt3 <- function(x) sprintf("%.3f", x)

prepare_main_results_table <- function(df_model_name) {
  final_results %>%
    filter(Model == df_model_name) %>%
    mutate(
      `Missingness Setting` = Missing_Setting,
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
      `Missingness Setting`,
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
      `Missingness Setting` = Missing_Setting,
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
      `Missingness Setting`,
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
# 16) SAVE BACKUP
# =============================================================================

objects_to_save <- c(
  "final_results",
  "combined_results",
  "combined_missingness_overall",
  "combined_missingness_variable",
  "missingness_summary",
  "missingness_variable_summary",
  "missingness_train_test_table",
  "missingness_train_test_numeric",
  "overfitting_summary",
  "ci_95_results",
  "wilcox_vs_stacking",
  "combined_rank_table_glm",
  "combined_rank_table_knn",
  "combined_rank_table_dt",
  "main_glm_word",
  "main_knn_word",
  "main_dt_word",
  "overfit_glm_word",
  "overfit_knn_word",
  "overfit_dt_word",
  "abm_variable_map",
  "abm_excluded_variables",
  "selected_abm_numeric_predictors",
  "combined_timing"
)

# =============================================================================
# PAIRED DIFFERENCES + WIN PROPORTIONS (R1-C5)
# =============================================================================

individual_methods_abm <- c("PMM", "RF", "CART", "NORM", "MIDASTOUCH",
                            "CCA", "Mean", "MissInd", "MissForest")

abm_paired_diff <- data.frame()

for (sc_name in unique(as.character(combined_results$Scenario))) {
  for (model_name in c("GLM", "KNN", "DT")) {
    sub <- combined_results %>%
      filter(as.character(Scenario) == sc_name,
             as.character(Model) == model_name,
             as.character(Imputation_Method) %in%
               c("Stacking", individual_methods_abm))
    if (nrow(sub) == 0) next
    wide <- sub %>%
      select(Split, Imputation_Method, Test_AUC) %>%
      mutate(Imputation_Method = as.character(Imputation_Method)) %>%
      pivot_wider(names_from = Imputation_Method, values_from = Test_AUC)
    if (!"Stacking" %in% names(wide)) next
    icols <- intersect(individual_methods_abm, names(wide))
    if (length(icols) == 0) next
    best_indiv <- apply(wide[, icols, drop = FALSE], 1, max, na.rm = TRUE)
    diff_vec   <- wide$Stacking - best_indiv
    n          <- sum(!is.na(diff_vec))
    if (n < 2) next
    mean_diff <- mean(diff_vec, na.rm = TRUE)
    sd_diff   <- sd(diff_vec, na.rm = TRUE)
    se_diff   <- sd_diff / sqrt(n)
    abm_paired_diff <- bind_rows(abm_paired_diff, data.frame(
      Scenario = sc_name, Model = model_name,
      Missing_Setting = missing_setting_label,
      Metric = "AUC", n_splits = n,
      Mean_Diff = mean_diff, SD_Diff = sd_diff,
      CI_Lower = mean_diff - qt(0.975, df = n-1) * se_diff,
      CI_Upper = mean_diff + qt(0.975, df = n-1) * se_diff,
      Win_Proportion = mean(diff_vec > 0, na.rm = TRUE)
    ))
  }
}

safe_write_csv(
  abm_paired_diff,
  file.path(output_dir, "post_analysis", "stacking_paired_diff_win_prop.csv")
)

# =============================================================================
# CALIBRATION ANALYSIS (R1-C8) — Cox slope/intercept summary
# (compute_calibration defined at top of script; called inline per split above)
# =============================================================================

# Summarise Cox calibration (slope + intercept) from per-split stored values.
# Cal_Slope and Cal_Intercept were computed inline for every split/method above.
# Here we aggregate across splits per Scenario × Model × Method.
cal_summary <- combined_results %>%
  group_by(Scenario, Missing_Setting, Imputation_Method, Model) %>%
  summarise(
    Mean_Cal_Slope      = mean(Cal_Slope,     na.rm = TRUE),
    SD_Cal_Slope        = sd(Cal_Slope,       na.rm = TRUE),
    Mean_Cal_Intercept  = mean(Cal_Intercept, na.rm = TRUE),
    SD_Cal_Intercept    = sd(Cal_Intercept,   na.rm = TRUE),
    n_splits_cal        = sum(!is.na(Cal_Slope)),
    Mean_Test_AUC       = mean(Test_AUC,      na.rm = TRUE),
    Mean_Test_Brier     = mean(Test_Brier,    na.rm = TRUE),
    Mean_Sensitivity    = mean(Test_Sensitivity, na.rm = TRUE),
    Mean_Specificity    = mean(Test_Specificity, na.rm = TRUE),
    n_splits            = n(),
    .groups = "drop"
  )

safe_write_csv(
  cal_summary,
  file.path(output_dir, "post_analysis", "calibration_summary.csv")
)

# Update objects to save
objects_to_save <- c(objects_to_save, "abm_paired_diff", "cal_summary")

safe_save_rdata(
  object_names = objects_to_save,
  path = file.path(output_dir, "post_analysis", "results_backup.RData"),
  env = environment()
)

# =============================================================================
# 17) CONSOLE MESSAGE
# =============================================================================

cat("\n==================== OUTPUTS SAVED ====================\n")
cat("Output directory:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("  - abm_variable_map.csv\n")
cat("  - abm_excluded_variables.csv\n")
cat("  - all_scenarios_raw_results.csv\n")
cat("  - all_scenarios_test_metrics.csv\n")
cat("  - all_scenarios_actual_missingness_variable_level.csv\n")
cat("  - all_scenarios_actual_missingness_overall.csv\n")
cat("  - all_scenarios_actual_missingness_comparison.csv\n")
cat("\nGenerated files in post_analysis folder:\n")
cat("  - missingness_mean_sd_summary.csv\n")
cat("  - missingness_variable_mean_sd_summary.csv\n")
cat("  - missingness_train_test_mean_sd_table.csv\n")
cat("  - missingness_train_test_mean_sd_numeric.csv\n")
cat("  - final_results_by_model.csv\n")
cat("  - test_metrics_95CI.csv\n")
cat("  - wilcoxon_stacking_vs_methods.csv\n")
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
cat("  - results_backup.RData\n")
cat("\nNew outputs (MIE_SE + timing):\n")
cat("  - ABM_per_split_timing.csv  (in output_dir root)\n")
cat("  - ABM_timing_summary.csv    (in output_dir root)\n")
cat("  - all_scenarios_per_split_timing.csv\n")
cat("  MIE_SE rows stored with Imputation_Method='MIE_SE', Model='MIE_SE' in all result files.\n")
cat("\nFigure files are generated separately by running R/02_abm_figures.R after this script.\n")
cat("=======================================================\n")
