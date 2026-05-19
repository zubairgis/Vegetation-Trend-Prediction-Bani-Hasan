# Reproducible workflow for spatial autocorrelation,
# Integrated model ranking and selection of the optimal predictive framework

# --------------------------------------------------
# 1. Required Packages
# --------------------------------------------------

packages <- c(
  "dplyr",
  "readr",
  "tidyr",
  "ggplot2",
  "tibble"
)

missing_packages <- packages[
  !sapply(packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    dependencies = TRUE
  )
}

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(tibble)

# --------------------------------------------------
# 2. Input and Output Paths
# --------------------------------------------------

stats_ann_path <-
  "./results/model_statistics.csv"

boot_sum_path <-
  "./results/RF_XGB_Bootstrap/bootstrap_stats_summary.csv"

resid_path <-
  "./results/Plots_Q1/Residuals/Residual_Summary_Test.csv"

morans_path <-
  "./results/Spatial_Autocorrelation/Global_MoransI_Table.csv"

lisa_path <-
  "./results/Spatial_Autocorrelation/LISA_Cluster_Counts_and_SignificantPercent.csv"

out_dir <-
  "./results/Model_Ranking"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# --------------------------------------------------
# 3. Min–Max Scaling Function
# --------------------------------------------------

minmax <- function(x) {
  
  rng <- range(
    x,
    na.rm = TRUE
  )
  
  if (
    is.infinite(rng[1]) ||
    is.infinite(rng[2]) ||
    diff(rng) == 0
  ) {
    
    return(
      rep(
        0.5,
        length(x)
      )
    )
  }
  
  (x - rng[1]) /
    (rng[2] - rng[1])
}

# --------------------------------------------------
# 4. Accuracy Metrics
# --------------------------------------------------

ann_stats <- read_csv(
  stats_ann_path,
  show_col_types = FALSE
) %>%
  filter(
    tolower(Split) == "test"
  ) %>%
  transmute(
    Model = Model,
    RMSE = RMSE,
    MAE = MAE,
    R2 = R2
  )

boot_stats <- read_csv(
  boot_sum_path,
  show_col_types = FALSE
) %>%
  transmute(
    Model = as.character(model),
    RMSE = RMSE_MN,
    MAE = MAE_MN,
    R2 = R2_MN
  )

acc <- bind_rows(
  ann_stats,
  boot_stats
) %>%
  filter(
    Model %in% c(
      "ANN",
      "ANN-PSO",
      "ANN-GWO",
      "RF",
      "XGB"
    )
  ) %>%
  group_by(Model) %>%
  summarise(
    RMSE =
      mean(
        RMSE,
        na.rm = TRUE
      ),
    
    MAE =
      mean(
        MAE,
        na.rm = TRUE
      ),
    
    R2 =
      mean(
        R2,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

# --------------------------------------------------
# 5. Residual Statistics
# --------------------------------------------------

resid <- read_csv(
  resid_path,
  show_col_types = FALSE
) %>%
  transmute(
    Model = Model,
    
    SD_Residual =
      SD_Residual,
    
    IQR_Residual =
      IQR_Residual
  ) %>%
  filter(
    Model %in% c(
      "ANN",
      "ANN-PSO",
      "ANN-GWO",
      "RF",
      "XGB"
    )
  )

# --------------------------------------------------
# 6. Global Moran's I
# --------------------------------------------------

mor <- read_csv(
  morans_path,
  show_col_types = FALSE
) %>%
  transmute(
    Variable = Variable,
    Moran_I = Moran_I
  )

mor_res <- tibble(
  
  Model = c(
    "ANN",
    "ANN-PSO",
    "ANN-GWO",
    "RF",
    "XGB"
  ),
  
  MoranI_Residual = c(
    
    mor$Moran_I[
      mor$Variable %in%
        c(
          "RES_ANN",
          "RES_ANN "
        )
    ][1],
    
    mor$Moran_I[
      mor$Variable %in%
        c(
          "RES_ANN_PSO",
          "RES_ANN-PSO"
        )
    ][1],
    
    mor$Moran_I[
      mor$Variable %in%
        c(
          "RES_ANN_GWO",
          "RES_ANN-GWO"
        )
    ][1],
    
    mor$Moran_I[
      mor$Variable %in%
        c("RES_RF")
    ][1],
    
    mor$Moran_I[
      mor$Variable %in%
        c("RES_XGB")
    ][1]
  )
)

# --------------------------------------------------
# 7. LISA Significant Percentage
# --------------------------------------------------

lisa <- read_csv(
  lisa_path,
  show_col_types = FALSE
) %>%
  transmute(
    Model = Model,
    SigPct = SigPct
  ) %>%
  filter(
    Model %in% c(
      "ANN",
      "ANN-PSO",
      "ANN-GWO",
      "RF",
      "XGB"
    )
  )

# --------------------------------------------------
# 8. Merge Evaluation Criteria
# --------------------------------------------------

allcrit <- acc %>%
  left_join(
    resid,
    by = "Model"
  ) %>%
  left_join(
    mor_res,
    by = "Model"
  ) %>%
  left_join(
    lisa,
    by = "Model"
  )

# --------------------------------------------------
# --------------------------------------------------
# 9. Ranking Weights
# --------------------------------------------------

w <- list(
  RMSE = 0.25,
  MAE = 0.20,
  R2 = 0.20,
  SD_Residual = 0.075,
  IQR_Residual = 0.075,
  MoranI_Residual = 0.10,
  SigPct = 0.10
)

# --------------------------------------------------
# 10. Composite Ranking
# --------------------------------------------------

rank_tbl <- allcrit %>%
  mutate(
    
    RMSE_s =
      1 - minmax(RMSE),
    
    MAE_s =
      1 - minmax(MAE),
    
    R2_s =
      minmax(R2),
    
    SDres_s =
      1 - minmax(SD_Residual),
    
    IQRres_s =
      1 - minmax(IQR_Residual),
    
    Moran_s =
      1 - minmax(MoranI_Residual),
    
    Sig_s =
      1 - minmax(SigPct),
    
    Score =
      w$RMSE * RMSE_s +
      w$MAE * MAE_s +
      w$R2 * R2_s +
      w$SD_Residual * SDres_s +
      w$IQR_Residual * IQRres_s +
      w$MoranI_Residual * Moran_s +
      w$SigPct * Sig_s
  ) %>%
  arrange(
    desc(Score)
  ) %>%
  mutate(
    Rank = row_number(),
    
    Score =
      round(
        Score,
        4
      )
  ) %>%
  select(
    Rank,
    Model,
    Score,
    RMSE,
    MAE,
    R2,
    SD_Residual,
    IQR_Residual,
    MoranI_Residual,
    SigPct
  )
# 11. Export Ranking Table
# --------------------------------------------------

write_csv(
  rank_tbl,
  file.path(
    out_dir,
    "Model_Ranking_Summary.csv"
  )
)

# --------------------------------------------------
# 12. Ranking Figure
# --------------------------------------------------

p_rank <- ggplot(
  rank_tbl,
  aes(
    x = reorder(
      Model,
      Score
    ),
    y = Score
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title =
      "Integrated model ranking",
    x = NULL,
    y = "Composite score"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      element_blank()
  )

ggsave(
  file.path(
    out_dir,
    "Figure_ModelRanking_Bars.png"
  ),
  p_rank,
  width = 7.5,
  height = 4.8,
  dpi = 400
)

ggsave(
  file.path(
    out_dir,
    "Figure_ModelRanking_Bars.pdf"
  ),
  p_rank,
  width = 7.5,
  height = 4.8
)
