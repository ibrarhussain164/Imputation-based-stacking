# =============================================================================
# ABM figure generation
# Purpose: create classifier-level summary figures from the completed ABM analysis.
# Run this script after R/01_abm_main_analysis.R.
# =============================================================================

rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

resolve_output_dir <- function() {
  env_dir <- Sys.getenv("ABM_OUTPUT_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(path.expand(env_dir))
  }
  pointer_file <- file.path(getwd(), "latest_abm_output_dir.txt")
  if (file.exists(pointer_file)) {
    dir_from_file <- trimws(readLines(pointer_file, warn = FALSE)[1])
    if (nzchar(dir_from_file)) {
      return(path.expand(dir_from_file))
    }
  }
  file.path(getwd(), "outputs")
}

output_dir <- resolve_output_dir()
results_file <- file.path(output_dir, "post_analysis", "final_results_by_model.csv")

if (!file.exists(results_file)) {
  stop(paste0(
    "Could not find final_results_by_model.csv at: ", results_file,
    "\nRun R/01_abm_main_analysis.R first, or set ABM_OUTPUT_DIR to the completed analysis output folder."
  ))
}

final_results <- readr::read_csv(results_file, show_col_types = FALSE)

# 14) FIGURES ONLY: 7 METRIC SUBPLOTS FOR EACH CLASSIFIER
# =============================================================================
# Output:
#   1) Logistic Regression: 7 metric subplots
#   2) K Nearest Neighbour: 7 metric subplots
#   3) Decision Tree: 7 metric subplots
#
# Saves:
#   - PNG at 900 dpi
#   - PDF vector format
#
# Font:
#   - Times New Roman for all text elements
#
# Important:
#   
# =============================================================================

# -----------------------------------------------------------------------------
# User-defined method colors
# -----------------------------------------------------------------------------

method_colors <- c(
  "Stacking"   = "blue",
  "PMM"        = "grey",
  "RF"         = "lightgreen",
  "CART"       = "lightpink",
  "NORM"       = "lightblue",
  "MIDASTOUCH" = "brown"
)

plot_method_order <- c(
  "Stacking",
  "PMM",
  "RF",
  "CART",
  "NORM",
  "MIDASTOUCH"
)

model_display_names <- c(
  "GLM" = "Logistic Regression",
  "KNN" = "K Nearest Neighbour",
  "DT"  = "Decision Tree"
)

times_font <- "Times New Roman"

# -----------------------------------------------------------------------------
# Create clean figure output folder
# -----------------------------------------------------------------------------

fig_root <- file.path(output_dir, "post_analysis", "figures")
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

# Remove old duplicate figure folders/files if they exist in the same run folder
old_fig_dirs <- c(
  file.path(fig_root, "separate_metric_figures"),
  file.path(fig_root, "separate_metric_figures_900dpi"),
  file.path(fig_root, "seven_combined_metric_figures_900dpi"),
  file.path(fig_root, "classifier_7_metric_subplots_900dpi")
)

for (old_dir in old_fig_dirs) {
  if (dir.exists(old_dir)) {
    unlink(old_dir, recursive = TRUE, force = TRUE)
    message("Removed old duplicate figure folder: ", old_dir)
  }
}

old_bar_files <- list.files(
  path = fig_root,
  pattern = "^(Bar_SD_|Separate_|Combined_|Updated_).*\\.(png|pdf)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(old_bar_files) > 0) {
  unlink(old_bar_files, force = TRUE)
  message("Removed old duplicate figure files from figure root.")
}

fig_dir <- file.path(
  fig_root,
  "classifier_7_metric_subplots_900dpi"
)

png_dir <- file.path(fig_dir, "PNG_900dpi")
pdf_dir <- file.path(fig_dir, "PDF")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Check required columns
# -----------------------------------------------------------------------------

required_figure_cols <- c(
  "Model", "Missing_Setting", "Method",
  "AUC", "SD_AUC",
  "Accuracy", "SD_Accuracy",
  "Brier", "SD_Brier",
  "Precision", "SD_Precision",
  "Sensitivity", "SD_Sensitivity",
  "Specificity", "SD_Specificity",
  "F1", "SD_F1"
)

missing_figure_cols <- setdiff(required_figure_cols, names(final_results))

if (length(missing_figure_cols) > 0) {
  stop(
    paste0(
      "These required columns are missing from final_results: ",
      paste(missing_figure_cols, collapse = ", ")
    )
  )
}

# -----------------------------------------------------------------------------
# Prepare data for plotting
# -----------------------------------------------------------------------------

metric_plot_long <- final_results %>%
  dplyr::mutate(
    Method = as.character(Method),
    Model  = as.character(Model)
  ) %>%
  dplyr::filter(
    Method %in% plot_method_order,
    Model %in% names(model_display_names)
  ) %>%
  dplyr::mutate(
    Model_Display = dplyr::recode(Model, !!!model_display_names),
    Method = factor(Method, levels = plot_method_order),
    Model_Display = factor(
      Model_Display,
      levels = c(
        "Logistic Regression",
        "K Nearest Neighbour",
        "Decision Tree"
      )
    )
  ) %>%
  dplyr::select(
    Model,
    Model_Display,
    Missing_Setting,
    Method,
    AUC,
    SD_AUC,
    Accuracy,
    SD_Accuracy,
    Brier,
    SD_Brier,
    Precision,
    SD_Precision,
    Sensitivity,
    SD_Sensitivity,
    Specificity,
    SD_Specificity,
    F1,
    SD_F1
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      AUC,
      Accuracy,
      Brier,
      Precision,
      Sensitivity,
      Specificity,
      F1
    ),
    names_to = "Metric",
    values_to = "Mean_Value"
  ) %>%
  dplyr::mutate(
    SD_Value = dplyr::case_when(
      Metric == "AUC"         ~ SD_AUC,
      Metric == "Accuracy"    ~ SD_Accuracy,
      Metric == "Brier"       ~ SD_Brier,
      Metric == "Precision"   ~ SD_Precision,
      Metric == "Sensitivity" ~ SD_Sensitivity,
      Metric == "Specificity" ~ SD_Specificity,
      Metric == "F1"          ~ SD_F1,
      TRUE                    ~ NA_real_
    ),
    Metric = factor(
      Metric,
      levels = c(
        "AUC",
        "Accuracy",
        "Brier",
        "Precision",
        "Sensitivity",
        "Specificity",
        "F1"
      )
    )
  )

# -----------------------------------------------------------------------------
# Safe save function: PNG 900 dpi + PDF
# -----------------------------------------------------------------------------

save_plot_png_pdf <- function(plot,
                              filename_base,
                              png_dir,
                              pdf_dir,
                              width = 14,
                              height = 12,
                              dpi = 900) {
  
  png_file <- file.path(png_dir, paste0(filename_base, ".png"))
  pdf_file <- file.path(pdf_dir, paste0(filename_base, ".pdf"))
  
  tryCatch(
    {
      ggplot2::ggsave(
        filename = png_file,
        plot     = plot,
        width    = width,
        height   = height,
        units    = "in",
        dpi      = dpi,
        bg       = "white"
      )
      message("Saved PNG: ", png_file)
    },
    error = function(e) {
      message("PNG save failed for ", filename_base, ": ", e$message)
    }
  )
  
  tryCatch(
    {
      ggplot2::ggsave(
        filename = pdf_file,
        plot     = plot,
        width    = width,
        height   = height,
        units    = "in",
        device   = grDevices::cairo_pdf,
        family   = times_font,
        bg       = "white"
      )
      message("Saved PDF: ", pdf_file)
    },
    error = function(e) {
      message("PDF save failed with cairo_pdf. Trying normal pdf device.")
      
      tryCatch(
        {
          ggplot2::ggsave(
            filename = pdf_file,
            plot     = plot,
            width    = width,
            height   = height,
            units    = "in",
            device   = "pdf",
            bg       = "white"
          )
          message("Saved PDF: ", pdf_file)
        },
        error = function(e2) {
          message("PDF save failed for ", filename_base, ": ", e2$message)
        }
      )
    }
  )
}

# -----------------------------------------------------------------------------
# Function: one combined figure containing 7 metric subplots for one classifier
# -----------------------------------------------------------------------------

plot_classifier_7_metrics <- function(df, model_name) {
  
  model_name <- as.character(model_name)
  
  df_sub <- df %>%
    dplyr::filter(Model == model_name) %>%
    dplyr::mutate(
      Method = factor(as.character(Method), levels = plot_method_order)
    )
  
  if (nrow(df_sub) == 0) {
    message("No data found for classifier: ", model_name)
    return(NULL)
  }
  
  classifier_title <- model_display_names[[model_name]]
  if (is.null(classifier_title) || is.na(classifier_title)) {
    classifier_title <- model_name
  }
  
  p <- ggplot2::ggplot(
    df_sub,
    ggplot2::aes(
      x = Method,
      y = Mean_Value,
      fill = Method
    )
  ) +
    ggplot2::geom_col(
      width = 0.75,
      color = "black",
      linewidth = 0.25
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = Mean_Value - SD_Value,
        ymax = Mean_Value + SD_Value
      ),
      width = 0.20,
      linewidth = 0.50
    ) +
    ggplot2::scale_fill_manual(
      values = method_colors,
      limits = plot_method_order,
      drop = FALSE
    ) +
    ggplot2::facet_wrap(
      ~ Metric,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::theme_minimal(
      base_size = 14,
      base_family = times_font
    ) +
    ggplot2::labs(
      title = classifier_title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = times_font,
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(
        family = times_font,
        angle = 45,
        hjust = 1,
        face = "bold",
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        family = times_font,
        color = "black"
      ),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        family = times_font,
        face = "bold",
        hjust = 0.5,
        color = "black",
        size = 20
      ),
      strip.text = ggplot2::element_text(
        family = times_font,
        face = "bold",
        color = "black",
        size = 14
      ),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )
  
  filename_base <- paste0("Figure_7_Metrics_", model_name)
  
  save_plot_png_pdf(
    plot          = p,
    filename_base = filename_base,
    png_dir       = png_dir,
    pdf_dir       = pdf_dir,
    width         = 14,
    height        = 12,
    dpi           = 900
  )
  
  return(p)
}

# -----------------------------------------------------------------------------
# Generate only the figures
# -----------------------------------------------------------------------------

models_to_plot <- c("GLM", "KNN", "DT")

classifier_7_metric_plots <- list()

for (model_name in models_to_plot) {
  classifier_7_metric_plots[[model_name]] <- plot_classifier_7_metrics(
    df = metric_plot_long,
    model_name = model_name
  )
}

cat("\n==================== FIGURES SAVED ====================\n")

cat("Figure folder:\n")
cat(fig_dir, "\n\n")

cat("PNG 900 dpi folder:\n")
cat(png_dir, "\n\n")

cat("PDF folder:\n")
cat(pdf_dir, "\n\n")

cat("Classifier figures created:\n")
cat(length(models_to_plot), "PNG files +",
    length(models_to_plot), "PDF files\n")

cat("Total figure files:\n")
cat(length(models_to_plot) * 2, "files\n")

cat("Each classifier figure contains 7 metric subplots.\n")
cat("===============================================================\n")


