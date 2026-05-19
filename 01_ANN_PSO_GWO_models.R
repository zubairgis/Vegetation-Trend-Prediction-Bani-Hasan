# Reproducible workflow for vegetation trend prediction
#ANN, ANN-PSO, and ANN-GWO models

# --------------------------------------------------
# 1. Required Packages
# --------------------------------------------------

packages <- c(
  "sf",
  "dplyr",
  "caret",
  "nnet",
  "pso",
  "doParallel",
  "parallel"
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
library(nnet)
library(pso)
library(doParallel)
library(parallel)

# --------------------------------------------------
# 2. Input and Output Paths
# --------------------------------------------------

shp_in <- "./data/Data.shp"
out_root <- "./results"

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "ANN"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "ANN_GWO"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "ANN_PSO"), showWarnings = FALSE, recursive = TRUE)

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
# 5. Predictor Standardization
# --------------------------------------------------

preproc <- preProcess(
  df0[, selected_vars],
  method = c("center", "scale")
)

X_scaled <- predict(preproc, df0[, selected_vars])

dat_scaled <- bind_cols(
  df0 %>% select(.row_id, all_of(target_var)),
  X_scaled
)

# --------------------------------------------------
# 6. Training and Testing Split
# --------------------------------------------------

set.seed(123)

train_idx <- createDataPartition(
  dat_scaled[[target_var]],
  p = 0.80,
  list = FALSE
)

train_dat <- dat_scaled[train_idx, ]
test_dat  <- dat_scaled[-train_idx, ]

# --------------------------------------------------
# 7. Parallel Processing
# --------------------------------------------------

cores_to_use <- max(1, detectCores() - 1)

cl <- makePSOCKcluster(cores_to_use)

registerDoParallel(cl)

# --------------------------------------------------
# 8. Evaluation Metrics
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
# 9. Baseline ANN Model
# --------------------------------------------------

ctrl5 <- trainControl(
  method = "cv",
  number = 5,
  allowParallel = TRUE
)

set.seed(123)

ann_model <- caret::train(
  reformulate(selected_vars, response = target_var),
  data = train_dat %>%
    select(all_of(c(target_var, selected_vars))),
  method = "nnet",
  trControl = ctrl5,
  linout = TRUE,
  trace = FALSE,
  maxit = 140,
  tuneGrid = expand.grid(
    size = c(3, 5, 7),
    decay = c(0.00, 0.01, 0.03)
  )
)

pred_tr_ann <- predict(ann_model, newdata = train_dat)

pred_te_ann <- predict(ann_model, newdata = test_dat)

m_tr_ann <- calc_metrics(
  train_dat[[target_var]],
  pred_tr_ann
)

m_te_ann <- calc_metrics(
  test_dat[[target_var]],
  pred_te_ann
)

# --------------------------------------------------
# 10. Optimization Dataset and Objective Function
# --------------------------------------------------

set.seed(123)

opt_n <- min(4000, nrow(train_dat))

opt_ids <- sample(
  seq_len(nrow(train_dat)),
  size = opt_n,
  replace = FALSE
)

opt_dat <- train_dat[opt_ids, ]

ctrl2 <- trainControl(
  method = "cv",
  number = 2,
  allowParallel = TRUE
)

obj_cv_rmse <- function(par) {
  
  size <- as.integer(round(par[1]))
  
  decay <- par[2]
  
  size <- max(1L, min(20L, size))
  
  decay <- max(0, min(0.10, decay))
  
  set.seed(123)
  
  fit <- tryCatch(
    
    caret::train(
      reformulate(selected_vars, response = target_var),
      data = opt_dat %>%
        select(all_of(c(target_var, selected_vars))),
      method = "nnet",
      trControl = ctrl2,
      linout = TRUE,
      trace = FALSE,
      maxit = 140,
      tuneGrid = data.frame(
        size = size,
        decay = decay
      )
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(1e9)
  }
  
  fit$results$RMSE[1]
}

# --------------------------------------------------
# 11. ANN-PSO Optimization
# --------------------------------------------------

set.seed(123)

pso_res <- psoptim(
  par = c(7, 0.02),
  fn = obj_cv_rmse,
  lower = c(1, 0.00),
  upper = c(20, 0.10),
  control = list(
    maxit = 12,
    s = 22,
    trace = 1
  )
)

best_size_pso <- as.integer(round(pso_res$par[1]))

best_decay_pso <- pso_res$par[2]

best_size_pso <- max(1L, min(20L, best_size_pso))

best_decay_pso <- max(0, min(0.10, best_decay_pso))

set.seed(123)

ann_pso <- caret::train(
  reformulate(selected_vars, response = target_var),
  data = train_dat %>%
    select(all_of(c(target_var, selected_vars))),
  method = "nnet",
  trControl = trainControl(method = "none"),
  linout = TRUE,
  trace = FALSE,
  maxit = 140,
  tuneGrid = data.frame(
    size = best_size_pso,
    decay = best_decay_pso
  )
)

pred_tr_pso <- predict(ann_pso, newdata = train_dat)

pred_te_pso <- predict(ann_pso, newdata = test_dat)

m_tr_pso <- calc_metrics(
  train_dat[[target_var]],
  pred_tr_pso
)

m_te_pso <- calc_metrics(
  test_dat[[target_var]],
  pred_te_pso
)

# --------------------------------------------------
# 12. ANN-GWO Optimization
# --------------------------------------------------

gwo_optimize <- function(
    fn,
    lower,
    upper,
    n_wolves = 9,
    n_iter = 12,
    seed = 123
) {
  
  set.seed(seed)
  
  dimn <- length(lower)
  
  X <- matrix(
    runif(
      n_wolves * dimn,
      min = lower,
      max = upper
    ),
    nrow = n_wolves,
    byrow = TRUE
  )
  
  fit <- apply(X, 1, fn)
  
  ord <- order(fit)
  
  alpha <- X[ord[1], ]
  beta  <- X[ord[2], ]
  delta <- X[ord[3], ]
  
  alpha_fit <- fit[ord[1]]
  
  for (t in seq_len(n_iter)) {
    
    a <- 2 - t * (2 / n_iter)
    
    for (i in seq_len(n_wolves)) {
      
      for (d in seq_len(dimn)) {
        
        r1 <- runif(1)
        r2 <- runif(1)
        
        A1 <- 2 * a * r1 - a
        C1 <- 2 * r2
        
        D_alpha <- abs(C1 * alpha[d] - X[i, d])
        
        X1 <- alpha[d] - A1 * D_alpha
        
        r1 <- runif(1)
        r2 <- runif(1)
        
        A2 <- 2 * a * r1 - a
        C2 <- 2 * r2
        
        D_beta <- abs(C2 * beta[d] - X[i, d])
        
        X2 <- beta[d] - A2 * D_beta
        
        r1 <- runif(1)
        r2 <- runif(1)
        
        A3 <- 2 * a * r1 - a
        C3 <- 2 * r2
        
        D_delta <- abs(C3 * delta[d] - X[i, d])
        
        X3 <- delta[d] - A3 * D_delta
        
        X[i, d] <- (X1 + X2 + X3) / 3
      }
      
      X[i, ] <- pmax(
        lower,
        pmin(upper, X[i, ])
      )
    }
    
    fit <- apply(X, 1, fn)
    
    ord <- order(fit)
    
    alpha <- X[ord[1], ]
    beta  <- X[ord[2], ]
    delta <- X[ord[3], ]
    
    alpha_fit <- fit[ord[1]]
  }
  
  list(
    best_par = alpha,
    best_value = alpha_fit
  )
}

gwo_res <- gwo_optimize(
  fn = obj_cv_rmse,
  lower = c(1, 0.00),
  upper = c(20, 0.10),
  n_wolves = 9,
  n_iter = 12,
  seed = 123
)

best_size_gwo <- as.integer(
  round(gwo_res$best_par[1])
)

best_decay_gwo <- gwo_res$best_par[2]

best_size_gwo <- max(
  1L,
  min(20L, best_size_gwo)
)

best_decay_gwo <- max(
  0,
  min(0.10, best_decay_gwo)
)

set.seed(123)

ann_gwo <- caret::train(
  reformulate(selected_vars, response = target_var),
  data = train_dat %>%
    select(all_of(c(target_var, selected_vars))),
  method = "nnet",
  trControl = trainControl(method = "none"),
  linout = TRUE,
  trace = FALSE,
  maxit = 140,
  tuneGrid = data.frame(
    size = best_size_gwo,
    decay = best_decay_gwo
  )
)

pred_tr_gwo <- predict(
  ann_gwo,
  newdata = train_dat
)

pred_te_gwo <- predict(
  ann_gwo,
  newdata = test_dat
)

m_tr_gwo <- calc_metrics(
  train_dat[[target_var]],
  pred_tr_gwo
)

m_te_gwo <- calc_metrics(
  test_dat[[target_var]],
  pred_te_gwo
)

# --------------------------------------------------
# 13. Full-Area Prediction
# --------------------------------------------------

full_pred_df <- dat_scaled %>%
  select(.row_id, all_of(selected_vars))

pred_all_ann <- as.numeric(
  predict(ann_model, newdata = full_pred_df)
)

pred_all_pso <- as.numeric(
  predict(ann_pso, newdata = full_pred_df)
)

pred_all_gwo <- as.numeric(
  predict(ann_gwo, newdata = full_pred_df)
)

out_sf <- pts_model %>%
  select(.row_id, geometry) %>%
  left_join(
    dat_scaled %>%
      select(.row_id, all_of(target_var)),
    by = ".row_id"
  ) %>%
  mutate(
    pred_ANN = pred_all_ann,
    pred_ANN_PSO = pred_all_pso,
    pred_ANN_GWO = pred_all_gwo
  )

# --------------------------------------------------
# 14. Export Spatial Outputs
# --------------------------------------------------

st_write(
  out_sf,
  file.path(out_root, "ANN", "ANN_pred.shp"),
  delete_layer = TRUE,
  quiet = TRUE
)

st_write(
  out_sf,
  file.path(out_root, "ANN_PSO", "ANN_PSO_pred.shp"),
  delete_layer = TRUE,
  quiet = TRUE
)

st_write(
  out_sf,
  file.path(out_root, "ANN_GWO", "ANN_GWO_pred.shp"),
  delete_layer = TRUE,
  quiet = TRUE
)

# --------------------------------------------------
# 15. Export Model Statistics
# --------------------------------------------------

stats <- bind_rows(
  
  data.frame(
    Model = "ANN",
    Split = "Train",
    m_tr_ann,
    size = ann_model$bestTune$size,
    decay = ann_model$bestTune$decay
  ),
  
  data.frame(
    Model = "ANN",
    Split = "Test",
    m_te_ann,
    size = ann_model$bestTune$size,
    decay = ann_model$bestTune$decay
  ),
  
  data.frame(
    Model = "ANN-PSO",
    Split = "Train",
    m_tr_pso,
    size = best_size_pso,
    decay = best_decay_pso
  ),
  
  data.frame(
    Model = "ANN-PSO",
    Split = "Test",
    m_te_pso,
    size = best_size_pso,
    decay = best_decay_pso
  ),
  
  data.frame(
    Model = "ANN-GWO",
    Split = "Train",
    m_tr_gwo,
    size = best_size_gwo,
    decay = best_decay_gwo
  ),
  
  data.frame(
    Model = "ANN-GWO",
    Split = "Test",
    m_te_gwo,
    size = best_size_gwo,
    decay = best_decay_gwo
  )
)

write.csv(
  stats,
  file.path(out_root, "model_statistics.csv"),
  row.names = FALSE
)

# --------------------------------------------------
# 16. Close Parallel Backend
# --------------------------------------------------

stopCluster(cl)

registerDoSEQ()
