# =============================================================================
# COVID-19 figure generation
# Purpose: create classifier-level summary figures from a completed analysis run.
# Run this script after R/01_covid19_mar_main_analysis.R.
# =============================================================================

rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
})

resolve_output_dir <- function() {
  env_dir <- Sys.getenv("COVID19_OUTPUT_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(path.expand(env_dir))
  }

  pointer_file <- file.path(getwd(), "latest_covid19_output_dir.txt")
  if (file.exists(pointer_file)) {
    dir_from_file <- trimws(readLines(pointer_file, warn = FALSE)[1])
    if (nzchar(dir_from_file)) {
      return(path.expand(dir_from_file))
    }
  }

  stop(paste0(
    "Could not resolve the analysis output folder. Run the main analysis first, ",
    "or set COVID19_OUTPUT_DIR to a completed output directory."
  ))
}

output_dir <- resolve_output_dir()
results_file <- file.path(output_dir, "post_analysis", "final_results_by_model.csv")

if (!file.exists(results_file)) {
  stop(paste0(
    "Could not find final_results_by_model.csv at: ", results_file, "\n",
    "Run R/01_covid19_mar_main_analysis.R first, or set COVID19_OUTPUT_DIR."
  ))
}

final_results <- readr::read_csv(results_file, show_col_types = FALSE)

fig_dir <- file.path(output_dir, "post_analysis", "figures", "classifier_7_metric_subplots")
png_dir <- file.path(fig_dir, "PNG_900dpi")
pdf_dir <- file.path(fig_dir, "PDF")
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

required_cols <- c(
  "Model", "Missing_Proportion", "Method",
  "AUC", "SD_AUC", "Accuracy", "SD_Accuracy", "Brier", "SD_Brier",
  "Precision", "SD_Precision", "Sensitivity", "SD_Sensitivity",
  "Specificity", "SD_Specificity", "F1", "SD_F1"
)

missing_cols <- setdiff(required_cols, names(final_results))
if (length(missing_cols) > 0) {
  stop(paste0(
    "These required columns are missing from final_results_by_model.csv: ",
    paste(missing_cols, collapse = ", ")
  ))
}

method_levels_plot <- c("Stacking", "PMM", "RF", "CART", "NORM", "MIDASTOUCH")

method_colors <- c(
  "Stacking"   = "blue",
  "PMM"        = "grey",
  "RF"         = "lightgreen",
  "CART"       = "lightpink",
  "NORM"       = "lightblue",
  "MIDASTOUCH" = "brown"
)

model_display_names <- c(
  "GLM" = "Logistic Regression",
  "KNN" = "K Nearest Neighbour",
  "DT"  = "Decision Tree"
)

metric_plot_long <- final_results %>%
  dplyr::filter(Method != "Complete") %>%
  dplyr::select(
    Model, Missing_Proportion, Method,
    AUC, SD_AUC,
    Accuracy, SD_Accuracy,
    Brier, SD_Brier,
    Precision, SD_Precision,
    Sensitivity, SD_Sensitivity,
    Specificity, SD_Specificity,
    F1, SD_F1
  ) %>%
  dplyr::mutate(
    Model = as.character(Model),
    Method = factor(as.character(Method), levels = method_levels_plot),
    Model_Display = dplyr::recode(Model, !!!model_display_names)
  ) %>%
  tidyr::pivot_longer(
    cols = c(AUC, Accuracy, Brier, Precision, Sensitivity, Specificity, F1),
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
      levels = c("AUC", "Accuracy", "Brier", "Precision", "Sensitivity", "Specificity", "F1")
    )
  )

plot_metric_panel <- function(data, metric_name, title_letter, model_name) {
  df_metric <- data %>% dplyr::filter(Metric == metric_name)
  dodge_width <- 0.05

  ggplot(
    df_metric,
    aes(
      x = Missing_Proportion,
      y = Mean_Value,
      color = Method
    )
  ) +
    geom_errorbar(
      aes(ymin = Mean_Value - SD_Value, ymax = Mean_Value + SD_Value),
      width = 0.015,
      linewidth = 0.6,
      position = position_dodge(width = dodge_width)
    ) +
    geom_point(
      size = 2.8,
      position = position_dodge(width = dodge_width)
    ) +
    scale_color_manual(
      values = method_colors[method_levels_plot],
      breaks = method_levels_plot,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = sort(unique(df_metric$Missing_Proportion)),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      title = paste0(title_letter, ". ", metric_name, " (", model_name, ")"),
      x = "Missing percentage",
      y = paste0("Mean ", metric_name, " (± SD)"),
      color = "Method"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )
}

make_metric_figure <- function(model_code, data_model) {
  display_name <- ifelse(model_code %in% names(model_display_names), model_display_names[[model_code]], model_code)

  metric_specs <- list(
    list(metric = "AUC",         letter = "A"),
    list(metric = "Accuracy",    letter = "B"),
    list(metric = "Brier",       letter = "C"),
    list(metric = "Precision",   letter = "D"),
    list(metric = "Sensitivity", letter = "E"),
    list(metric = "Specificity", letter = "F"),
    list(metric = "F1",          letter = "G")
  )

  plot_list <- lapply(metric_specs, function(m) {
    plot_metric_panel(
      data         = data_model,
      metric_name  = m$metric,
      title_letter = m$letter,
      model_name   = display_name
    )
  })

  fig_out <- ggpubr::ggarrange(
    plotlist = plot_list,
    ncol = 2,
    nrow = 4,
    common.legend = TRUE,
    legend = "bottom"
  )

  png_file <- file.path(png_dir, paste0("Fig_", model_code, "_7Metrics_MeanSD_NoComplete.png"))
  pdf_file <- file.path(pdf_dir, paste0("Fig_", model_code, "_7Metrics_MeanSD_NoComplete.pdf"))

  ggplot2::ggsave(filename = png_file, plot = fig_out, width = 14, height = 16, dpi = 900)
  ggplot2::ggsave(filename = pdf_file, plot = fig_out, width = 14, height = 16)

  message("Saved: ", png_file)
  message("Saved: ", pdf_file)

  fig_out
}

metric_data_list <- list(
  GLM = metric_plot_long %>% dplyr::filter(Model == "GLM"),
  KNN = metric_plot_long %>% dplyr::filter(Model == "KNN"),
  DT  = metric_plot_long %>% dplyr::filter(Model == "DT")
)

figure_list <- lapply(names(metric_data_list), function(model_code) {
  make_metric_figure(
    model_code = model_code,
    data_model = metric_data_list[[model_code]]
  )
})

names(figure_list) <- names(metric_data_list)

cat("\nFigure generation completed.\n")
cat("Figure directory:", fig_dir, "\n")
