# Reproducible workflow for 
#RF and XGBoost bootstrap modeling

# --------------------------------------------------
# 1. Required Packages
# --------------------------------------------------

packages <- c(
  "sf",
  "dplyr",
  "caret",
  "ranger",
  "xgboost",
  "Matrix"
)

missing_packages <- packages[
  !sapply(packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

library(sf)
library(dplyr)
library(caret)
library(ranger)
library(xgboost)
library(Matrix)

# --------------------------------------------------
# 2. Input and Output Paths
# --------------------------------------------------

shp_in <- "./data/Data.shp"

out_dir <- "./results/RF_XGB_Bootstrap"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------
# 3. Predictor and Response Variables
# --------------------------------------------------

selected_vars <- c(
  "SS_EVI", "K_NDWI", "SS_NDWI",
  "C_bio15", "C_bio10", "C_bio17",
  "C_bio8", "C_bio18", "C_bio12",
  "C_bio13", "K_LST", "C_bio7",
  "SS_LST", "C_bio6", "C_bio5",
  "C_bio1", "C_bio4", "C_bio19",
  "NL", "C_bio2", "C_bio3"
)

target_var <- "K_NDVI"

required_vars <- c(selected_vars, target_var)

# --------------------------------------------------
# 4. Data Preparation
# --------------------------------------------------

pts <- st_read(shp_in, quiet = TRUE) %>%
  mutate(.row_id = row_number())

missing_fields <- setdiff(required_vars, names(pts))

if (length(missing_fields) > 0) {
  stop(
    paste(
      "Missing required fields:",
      paste(missing_fields, collapse = ", ")
    )
  )
}

df0 <- pts %>%
  st_drop_geometry() %>%
  select(.row_id, all_of(required_vars)) %>%
  filter(if_all(all_of(required_vars), ~ !is.na(.)))

pts_model <- pts %>%
  inner_join(df0 %>% select(.row_id), by = ".row_id")

# --------------------------------------------------
# 5. Training and Testing Split
# --------------------------------------------------

set.seed(123)

train_idx <- createDataPartition(
  df0[[target_var]],
  p = 0.80,
  list = FALSE
)

train_df <- df0[train_idx, ]

test_df <- df0[-train_idx, ]

# --------------------------------------------------
# 6. Evaluation Metrics
# --------------------------------------------------

calc_metrics <- function(obs, pred) {
  
  obs <- as.numeric(obs)
  
  pred <- as.numeric(pred)
  
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  
  mae <- mean(abs(obs - pred), na.rm = TRUE)
  
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  
  ss_tot <- sum(
    (obs - mean(obs, na.rm = TRUE))^2,
    na.rm = TRUE
  )
  
  r2 <- 1 - (ss_res / ss_tot)
  
  data.frame(
    RMSE = rmse,
    MAE = mae,
    R2 = r2
  )
}

# --------------------------------------------------
# 7. XGBoost Matrix Functions
# --------------------------------------------------

make_dmat <- function(df) {
  
  x <- as.matrix(df[, selected_vars])
  
  x <- Matrix(x, sparse = TRUE)
  
  xgb.DMatrix(
    data = x,
    label = df[[target_var]]
  )
}

make_xmat <- function(df) {
  
  x <- as.matrix(df[, selected_vars])
  
  Matrix(x, sparse = TRUE)
}

# --------------------------------------------------
# 8. Bootstrap Settings
# --------------------------------------------------

B <- 30

set.seed(123)

n_all <- nrow(df0)

pred_rf_mat <- matrix(
  NA_real_,
  nrow = n_all,
  ncol = B
)

pred_xgb_mat <- matrix(
  NA_real_,
  nrow = n_all,
  ncol = B
)

stats_long <- data.frame(
  iter = integer(0),
  model = character(0),
  RMSE = numeric(0),
  MAE = numeric(0),
  R2 = numeric(0),
  stringsAsFactors = FALSE
)

all_feat <- df0[, c(".row_id", selected_vars, target_var)]

# --------------------------------------------------
# 9. Bootstrap Modeling
# --------------------------------------------------

for (b in seq_len(B)) {
  
  boot_ids <- sample(
    seq_len(nrow(train_df)),
    size = nrow(train_df),
    replace = TRUE
  )
  
  boot_df <- train_df[boot_ids, ]
  
  rf_fit <- ranger(
    formula = reformulate(selected_vars, response = target_var),
    data = boot_df[, c(selected_vars, target_var)],
    num.trees = 500,
    importance = "none",
    mtry = max(1, floor(sqrt(length(selected_vars)))),
    min.node.size = 5,
    seed = 123 + b
  )
  
  pred_rf_mat[, b] <- predict(
    rf_fit,
    data = all_feat[, selected_vars]
  )$predictions
  
  rf_te <- predict(
    rf_fit,
    data = test_df[, selected_vars]
  )$predictions
  
  m_rf <- calc_metrics(
    test_df[[target_var]],
    rf_te
  )
  
  stats_long <- rbind(
    stats_long,
    data.frame(
      iter = b,
      model = "RF",
      m_rf
    )
  )
  
  set.seed(1000 + b)
  
  idx_sub <- sample(
    seq_len(nrow(boot_df)),
    size = floor(0.80 * nrow(boot_df))
  )
  
  dtrain <- make_dmat(boot_df[idx_sub, ])
  
  dval <- make_dmat(boot_df[-idx_sub, ])
  
  params <- list(
    booster = "gbtree",
    objective = "reg:squarederror",
    eval_metric = "rmse",
    eta = 0.05,
    max_depth = 6,
    min_child_weight = 1,
    subsample = 0.80,
    colsample_bytree = 0.80,
    gamma = 0,
    lambda = 1
  )
  
  xgb_fit <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 800,
    watchlist = list(
      train = dtrain,
      val = dval
    ),
    early_stopping_rounds = 30,
    verbose = 0
  )
  
  all_xmat <- make_xmat(all_feat)
  
  pred_xgb_mat[, b] <- predict(
    xgb_fit,
    newdata = all_xmat
  )
  
  test_xmat <- make_xmat(test_df)
  
  xgb_te <- predict(
    xgb_fit,
    newdata = test_xmat
  )
  
  m_xgb <- calc_metrics(
    test_df[[target_var]],
    xgb_te
  )
  
  stats_long <- rbind(
    stats_long,
    data.frame(
      iter = b,
      model = "XGB",
      m_xgb
    )
  )
}

# --------------------------------------------------
# 10. Bootstrap Prediction Summary
# --------------------------------------------------

rf_mean <- rowMeans(
  pred_rf_mat,
  na.rm = TRUE
)

rf_sd <- apply(
  pred_rf_mat,
  1,
  sd,
  na.rm = TRUE
)

xgb_mean <- rowMeans(
  pred_xgb_mat,
  na.rm = TRUE
)

xgb_sd <- apply(
  pred_xgb_mat,
  1,
  sd,
  na.rm = TRUE
)

# --------------------------------------------------
# 11. Export Spatial Outputs
# --------------------------------------------------

out_sf <- pts_model %>%
  select(.row_id, geometry) %>%
  left_join(
    df0 %>% select(.row_id, all_of(target_var)),
    by = ".row_id"
  ) %>%
  mutate(
    RF_MN = as.numeric(rf_mean),
    RF_SD = as.numeric(rf_sd),
    XGB_MN = as.numeric(xgb_mean),
    XGB_SD = as.numeric(xgb_sd)
  )

out_shp <- file.path(
  out_dir,
  "RF_XGB_boot_pred.shp"
)

st_write(
  out_sf,
  out_shp,
  delete_layer = TRUE,
  quiet = TRUE
)

# --------------------------------------------------
# 12. Export Bootstrap Statistics
# --------------------------------------------------

long_csv <- file.path(
  out_dir,
  "bootstrap_stats_long.csv"
)

write.csv(
  stats_long,
  long_csv,
  row.names = FALSE
)

summary_stats <- stats_long %>%
  group_by(model) %>%
  summarise(
    B = n(),
    RMSE_MN = mean(RMSE, na.rm = TRUE),
    RMSE_SD = sd(RMSE, na.rm = TRUE),
    MAE_MN = mean(MAE, na.rm = TRUE),
    MAE_SD = sd(MAE, na.rm = TRUE),
    R2_MN = mean(R2, na.rm = TRUE),
    R2_SD = sd(R2, na.rm = TRUE),
    .groups = "drop"
  )

summary_csv <- file.path(
  out_dir,
  "bootstrap_stats_summary.csv"
)

write.csv(
  summary_stats,
  summary_csv,
  row.names = FALSE
)
