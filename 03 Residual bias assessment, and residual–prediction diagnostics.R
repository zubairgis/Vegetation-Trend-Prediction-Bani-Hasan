# Reproducible workflow for residual distribution,
#Residual bias assessment, and residual–prediction diagnostics

# --------------------------------------------------
# 1. Required Packages
# --------------------------------------------------

packages <- c(
  "sf",
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr"
)

missing_packages <- packages[
  !sapply(packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# --------------------------------------------------
# 2. Input and Output Paths
# --------------------------------------------------

models_shp <- "./data/Models.shp"

out_dir_resid <- "./results/Plots_Q1/Residuals"

out_dir_rvf <- "./results/Plots_Q1/Residuals_vs_Fitted"

dir.create(
  out_dir_resid,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  out_dir_rvf,
  showWarnings = FALSE,
  recursive = TRUE
)

# --------------------------------------------------
# 3. Read Spatial Data
# --------------------------------------------------

sf_all <- st_read(
  models_shp,
  quiet = TRUE
)

required_fields <- c(
  "K_NDVI",
  "prd_ANN",
  "p_ANN_P",
  "p_ANN_G",
  "RF_MN",
  "XGB_MN"
)

missing_fields <- setdiff(
  required_fields,
  names(sf_all)
)

if (length(missing_fields) > 0) {
  
  stop(
    paste(
      "Missing required fields:",
      paste(missing_fields, collapse = ", ")
    )
  )
}

df <- sf_all %>%
  st_drop_geometry()

# --------------------------------------------------
# 4. Identify Test Dataset
# --------------------------------------------------

split_candidates <- c(
  "split",
  "Split",
  "dataset",
  "Dataset",
  "set",
  "Set",
  "data_split",
  "DATA_SPLIT",
  "foldset",
  "FOLDSET"
)

split_col <- intersect(
  split_candidates,
  names(df)
)

split_col <- if (
  length(split_col) > 0
) split_col[1] else NA_character_

if (!is.na(split_col)) {
  
  split_vals <- tolower(
    as.character(df[[split_col]])
  )
  
  is_test <- split_vals %in% c(
    "test",
    "testing",
    "valid",
    "validation",
    "holdout",
    "val"
  )
  
  if (sum(is_test, na.rm = TRUE) == 0) {
    
    warning(
      paste(
        "Split column found but no TEST rows detected.",
        "Using all rows."
      )
    )
    
    is_test <- rep(TRUE, nrow(df))
  }
  
} else {
  
  warning(
    paste(
      "No split column found.",
      "Using all rows for residual diagnostics."
    )
  )
  
  is_test <- rep(TRUE, nrow(df))
}

df_test <- df[
  is_test,
  ,
  drop = FALSE
]

# --------------------------------------------------
# 5. Residual Long Table
# --------------------------------------------------

res_long <- df_test %>%
  transmute(
    Observed = K_NDVI,
    
    ANN = Observed - prd_ANN,
    
    `ANN-PSO` = Observed - p_ANN_P,
    
    `ANN-GWO` = Observed - p_ANN_G,
    
    RF = Observed - RF_MN,
    
    XGB = Observed - XGB_MN
  ) %>%
  pivot_longer(
    cols = -Observed,
    names_to = "Model",
    values_to = "Residual"
  ) %>%
  filter(!is.na(Residual)) %>%
  mutate(
    Model = factor(
      Model,
      levels = c(
        "ANN",
        "ANN-PSO",
        "ANN-GWO",
        "RF",
        "XGB"
      )
    )
  )

# --------------------------------------------------
# 6. Residual Summary Statistics
# --------------------------------------------------

res_summary <- res_long %>%
  group_by(Model) %>%
  summarise(
    N = n(),
    
    Mean_Residual =
      mean(Residual),
    
    Median_Residual =
      median(Residual),
    
    SD_Residual =
      sd(Residual),
    
    IQR_Residual =
      IQR(Residual),
    
    P05 =
      quantile(Residual, 0.05),
    
    P95 =
      quantile(Residual, 0.95),
    
    .groups = "drop"
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 6)
    )
  )

write_csv(
  res_summary,
  file.path(
    out_dir_resid,
    "Residual_Summary_Test.csv"
  )
)

# --------------------------------------------------
# 7. Plot Theme
# --------------------------------------------------

q1_theme <- function() {
  
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# --------------------------------------------------
# 8. Residual Density Plot
# --------------------------------------------------

p_density <- ggplot(
  res_long,
  aes(
    x = Residual,
    group = Model
  )
) +
  geom_density(
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~Model,
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    title = "Residual density (test set)",
    x = "Residual (Observed − Predicted)",
    y = "Density"
  ) +
  q1_theme()

ggsave(
  file.path(
    out_dir_resid,
    "ResidualDensity_Test.png"
  ),
  p_density,
  width = 10.5,
  height = 6.0,
  dpi = 400
)

ggsave(
  file.path(
    out_dir_resid,
    "ResidualDensity_Test.pdf"
  ),
  p_density,
  width = 10.5,
  height = 6.0
)

# --------------------------------------------------
# 9. Residual Boxplot
# --------------------------------------------------

p_box <- ggplot(
  res_long,
  aes(
    x = Model,
    y = Residual
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_boxplot(
    outlier.alpha = 0.25,
    linewidth = 0.6
  ) +
  labs(
    title = "Residual distribution (test set)",
    x = NULL,
    y = "Residual (Observed − Predicted)"
  ) +
  q1_theme()

ggsave(
  file.path(
    out_dir_resid,
    "ResidualBoxplot_Test.png"
  ),
  p_box,
  width = 8.5,
  height = 5.5,
  dpi = 400
)

ggsave(
  file.path(
    out_dir_resid,
    "ResidualBoxplot_Test.pdf"
  ),
  p_box,
  width = 8.5,
  height = 5.5
)

# --------------------------------------------------
# 10. Residual–Prediction Relationships
# --------------------------------------------------

res_fit <- df_test %>%
  transmute(
    Observed = K_NDVI,
    
    ANN_fit = prd_ANN,
    
    `ANN-PSO_fit` = p_ANN_P,
    
    `ANN-GWO_fit` = p_ANN_G,
    
    RF_fit = RF_MN,
    
    XGB_fit = XGB_MN
  ) %>%
  pivot_longer(
    cols = ends_with("_fit"),
    names_to = "Model",
    values_to = "Fitted"
  ) %>%
  mutate(
    Model = gsub("_fit", "", Model),
    
    Residual = Observed - Fitted,
    
    Model = factor(
      Model,
      levels = c(
        "ANN",
        "ANN-PSO",
        "ANN-GWO",
        "RF",
        "XGB"
      )
    )
  ) %>%
  filter(
    !is.na(Fitted) &
      !is.na(Residual)
  )

# --------------------------------------------------
# 11. Residual vs Fitted Plot
# --------------------------------------------------

p_rvf <- ggplot(
  res_fit,
  aes(
    x = Fitted,
    y = Residual
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_point(
    alpha = 0.25,
    size = 0.6
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~Model,
    ncol = 3,
    scales = "free_x"
  ) +
  labs(
    title = "Residuals vs fitted values",
    x = "Fitted NDVI trend",
    y = "Residual (Observed − Predicted)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    out_dir_rvf,
    "Residual_vs_Fitted.png"
  ),
  p_rvf,
  width = 11,
  height = 6.2,
  dpi = 400
)

ggsave(
  file.path(
    out_dir_rvf,
    "Residual_vs_Fitted.pdf"
  ),
  p_rvf,
  width = 11,
  height = 6.2
)

# --------------------------------------------------
# 12. Binned Residual Spread
# --------------------------------------------------

bin_stats <- res_fit %>%
  group_by(Model) %>%
  mutate(
    Fitted_bin = ntile(Fitted, 20)
  ) %>%
  group_by(Model, Fitted_bin) %>%
  summarise(
    bin_n = n(),
    
    fitted_mean =
      mean(Fitted),
    
    resid_sd =
      sd(Residual),
    
    resid_median =
      median(Residual),
    
    .groups = "drop"
  )

hetero_summary <- bin_stats %>%
  group_by(Model) %>%
  summarise(
    bins = n(),
    
    mean_bin_resid_sd =
      mean(resid_sd, na.rm = TRUE),
    
    max_bin_resid_sd =
      max(resid_sd, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(mean_bin_resid_sd)

write_csv(
  hetero_summary,
  file.path(
    out_dir_rvf,
    "Heteroscedasticity_BinnedSD_Summary.csv"
  )
)

write_csv(
  bin_stats,
  file.path(
    out_dir_rvf,
    "Heteroscedasticity_BinnedSD_Long.csv"
  )
)

# --------------------------------------------------
# 13. Binned Residual SD Plot
# --------------------------------------------------

p_binsd <- ggplot(
  bin_stats,
  aes(
    x = fitted_mean,
    y = resid_sd
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  facet_wrap(
    ~Model,
    ncol = 3,
    scales = "free_x"
  ) +
  labs(
    title = "Residual spread across fitted-value bins",
    x = "Mean fitted NDVI trend (bin)",
    y = "SD of residuals (within bin)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    out_dir_rvf,
    "BinnedResidualSD_vs_Fitted.png"
  ),
  p_binsd,
  width = 11,
  height = 6.2,
  dpi = 400
)

ggsave(
  file.path(
    out_dir_rvf,
    "BinnedResidualSD_vs_Fitted.pdf"
  ),
  p_binsd,
  width = 11,
  height = 6.2
)
