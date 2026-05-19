# Reproducible workflow for spatial autocorrelation,
#Global Moran's I, and LISA cluster analysis

# --------------------------------------------------
# 1. Required Packages
# --------------------------------------------------

packages <- c(
  "sf",
  "dplyr",
  "spdep",
  "tmap",
  "readr",
  "tibble"
)

missing_packages <- packages[
  !sapply(packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

library(sf)
library(dplyr)
library(spdep)
library(tmap)
library(readr)
library(tibble)

# --------------------------------------------------
# 2. Input and Output Paths
# --------------------------------------------------

models_shp <- "./data/Models.shp"

boundary_shp <- "./data/Boundary.shp"

out_dir <- "./results/Spatial_Autocorrelation"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# --------------------------------------------------
# 3. Parameters
# --------------------------------------------------

k_nn <- 8

alpha_lisa <- 0.05

use_row_standardized <- TRUE

# --------------------------------------------------
# 4. Read Spatial Data
# --------------------------------------------------

sf_all <- st_read(
  models_shp,
  quiet = TRUE
)

boundary <- NULL

if (file.exists(boundary_shp)) {
  
  boundary <- st_read(
    boundary_shp,
    quiet = TRUE
  )
}

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

sf_use <- sf_all %>%
  filter(
    if_all(
      all_of(required_fields),
      ~ !is.na(.)
    )
  )

# --------------------------------------------------
# 5. Coordinate Reference System Handling
# --------------------------------------------------

if (is.na(st_crs(sf_use))) {
  
  stop(
    "Models.shp has no CRS defined."
  )
}

if (st_is_longlat(sf_use)) {
  
  cen <- st_coordinates(
    st_centroid(
      st_union(
        st_geometry(sf_use)
      )
    )
  )
  
  lon <- cen[1]
  
  lat <- cen[2]
  
  utm_zone <- floor((lon + 180) / 6) + 1
  
  epsg <- ifelse(
    lat >= 0,
    32600 + utm_zone,
    32700 + utm_zone
  )
  
  sf_use <- st_transform(
    sf_use,
    epsg
  )
  
  if (
    !is.null(boundary) &&
    !is.na(st_crs(boundary)) &&
    st_crs(boundary) != st_crs(sf_use)
  ) {
    
    boundary <- st_transform(
      boundary,
      st_crs(sf_use)
    )
  }
}

# --------------------------------------------------
# 6. Spatial Weights Matrix
# --------------------------------------------------

coords <- st_coordinates(sf_use)

knn <- knearneigh(
  coords,
  k = k_nn
)

nb <- knn2nb(knn)

styleW <- ifelse(
  use_row_standardized,
  "W",
  "B"
)

lw <- nb2listw(
  nb,
  style = styleW,
  zero.policy = TRUE
)

# --------------------------------------------------
# 7. Residual Calculation
# --------------------------------------------------

sf_use <- sf_use %>%
  mutate(
    RES_ANN = K_NDVI - prd_ANN,
    
    RES_ANN_PSO = K_NDVI - p_ANN_P,
    
    RES_ANN_GWO = K_NDVI - p_ANN_G,
    
    RES_RF = K_NDVI - RF_MN,
    
    RES_XGB = K_NDVI - XGB_MN
  )

vars <- c(
  "Observed" = "K_NDVI",
  
  "RES_ANN" = "RES_ANN",
  
  "RES_ANN-PSO" = "RES_ANN_PSO",
  
  "RES_ANN-GWO" = "RES_ANN_GWO",
  
  "RES_RF" = "RES_RF",
  
  "RES_XGB" = "RES_XGB"
)

# --------------------------------------------------
# 8. Global Moran's I
# --------------------------------------------------

global_tbl <- lapply(
  names(vars),
  function(nm) {
    
    v <- sf_use[[vars[[nm]]]]
    
    mt <- moran.test(
      v,
      lw,
      randomisation = TRUE,
      zero.policy = TRUE
    )
    
    data.frame(
      Variable = nm,
      
      Moran_I =
        unname(
          mt$estimate[
            ["Moran I statistic"]
          ]
        ),
      
      Expected_I =
        unname(
          mt$estimate[
            ["Expectation"]
          ]
        ),
      
      Variance =
        unname(
          mt$estimate[
            ["Variance"]
          ]
        ),
      
      Z =
        unname(
          (
            mt$estimate[
              ["Moran I statistic"]
            ] -
              mt$estimate[
                ["Expectation"]
              ]
          ) /
            sqrt(
              mt$estimate[
                ["Variance"]
              ]
            )
        ),
      
      p_value = mt$p.value,
      
      kNN = k_nn,
      
      Weights = styleW,
      
      N = nrow(sf_use),
      
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    Moran_I =
      round(Moran_I, 6),
    
    Expected_I =
      round(Expected_I, 6),
    
    Variance =
      signif(Variance, 6),
    
    Z =
      round(Z, 3),
    
    p_value =
      signif(p_value, 3)
  )

write_csv(
  global_tbl,
  file.path(
    out_dir,
    "Global_MoransI_Table.csv"
  )
)

# --------------------------------------------------
# 9. Local Moran's I Function
# --------------------------------------------------

lisa_cluster <- function(
    x,
    lw,
    alpha = 0.05
) {
  
  x <- as.numeric(x)
  
  z <- as.numeric(scale(x))
  
  lagz <- lag.listw(
    lw,
    z,
    zero.policy = TRUE
  )
  
  lagz <- as.numeric(scale(lagz))
  
  li <- localmoran(
    z,
    lw,
    zero.policy = TRUE
  )
  
  p <- li[, 5]
  
  cl <- rep(
    "Not significant",
    length(z)
  )
  
  sig <- p <= alpha
  
  cl[
    sig & z >= 0 & lagz >= 0
  ] <- "High-High"
  
  cl[
    sig & z <= 0 & lagz <= 0
  ] <- "Low-Low"
  
  cl[
    sig & z >= 0 & lagz <= 0
  ] <- "High-Low"
  
  cl[
    sig & z <= 0 & lagz >= 0
  ] <- "Low-High"
  
  list(
    cluster = factor(
      cl,
      levels = c(
        "High-High",
        "Low-Low",
        "High-Low",
        "Low-High",
        "Not significant"
      )
    ),
    
    pvalue = p,
    
    Ii = li[, 1]
  )
}

# --------------------------------------------------
# 10. LISA Residual Clusters
# --------------------------------------------------

lisa_res <- list(
  
  ANN =
    lisa_cluster(
      sf_use$RES_ANN,
      lw,
      alpha_lisa
    ),
  
  `ANN-PSO` =
    lisa_cluster(
      sf_use$RES_ANN_PSO,
      lw,
      alpha_lisa
    ),
  
  `ANN-GWO` =
    lisa_cluster(
      sf_use$RES_ANN_GWO,
      lw,
      alpha_lisa
    ),
  
  RF =
    lisa_cluster(
      sf_use$RES_RF,
      lw,
      alpha_lisa
    ),
  
  XGB =
    lisa_cluster(
      sf_use$RES_XGB,
      lw,
      alpha_lisa
    )
)

sf_map <- sf_use %>%
  mutate(
    LISA_ANN =
      lisa_res$ANN$cluster,
    
    LISA_ANN_PSO =
      lisa_res$`ANN-PSO`$cluster,
    
    LISA_ANN_GWO =
      lisa_res$`ANN-GWO`$cluster,
    
    LISA_RF =
      lisa_res$RF$cluster,
    
    LISA_XGB =
      lisa_res$XGB$cluster
  )

# --------------------------------------------------
# 11. Export LISA Shapefile
# --------------------------------------------------

st_write(
  sf_map,
  file.path(
    out_dir,
    "Residuals_with_LISA_Clusters.shp"
  ),
  delete_layer = TRUE,
  quiet = TRUE
)

# --------------------------------------------------
# 12. LISA Mapping
# --------------------------------------------------

tmap_mode("plot")

pal_lisa <- c(
  "High-High" = "#b2182b",
  
  "Low-Low" = "#2166ac",
  
  "High-Low" = "#ef8a62",
  
  "Low-High" = "#67a9cf",
  
  "Not significant" = "#f0f0f0"
)

tm_lisa <- function(
    field,
    title_txt
) {
  
  base <- if (
    !is.null(boundary)
  ) {
    
    tm_shape(boundary) +
      tm_borders(
        lwd = 1,
        col = "black"
      )
    
  } else {
    
    tm_shape(sf_map[0, ])
  }
  
  base +
    tm_shape(sf_map) +
    tm_symbols(
      col = field,
      shape = 15,
      size = 0.06,
      palette = pal_lisa,
      border.col = NA,
      alpha = 1,
      title = "LISA cluster"
    ) +
    tm_layout(
      title = title_txt,
      title.size = 1.0,
      frame = FALSE,
      legend.outside = TRUE,
      legend.outside.position = "right"
    )
}

m1 <- tm_lisa(
  "LISA_ANN",
  "ANN residuals"
)

m2 <- tm_lisa(
  "LISA_ANN_PSO",
  "ANN-PSO residuals"
)

m3 <- tm_lisa(
  "LISA_ANN_GWO",
  "ANN-GWO residuals"
)

m4 <- tm_lisa(
  "LISA_RF",
  "RF residuals"
)

m5 <- tm_lisa(
  "LISA_XGB",
  "XGB residuals"
)

fig_lisa <- tmap_arrange(
  m1,
  m2,
  m3,
  m4,
  m5,
  ncol = 3
)

tmap_save(
  fig_lisa,
  file.path(
    out_dir,
    "LISA_Residual_Clusters.png"
  ),
  width = 14,
  height = 8,
  dpi = 600
)

tmap_save(
  fig_lisa,
  file.path(
    out_dir,
    "LISA_Residual_Clusters.pdf"
  ),
  width = 14,
  height = 8
)

# --------------------------------------------------
# 13. LISA Cluster Summary
# --------------------------------------------------

cluster_summary <- tibble(
  
  Model = c(
    "ANN",
    "ANN-PSO",
    "ANN-GWO",
    "RF",
    "XGB"
  ),
  
  HighHigh = c(
    sum(sf_map$LISA_ANN == "High-High"),
    sum(sf_map$LISA_ANN_PSO == "High-High"),
    sum(sf_map$LISA_ANN_GWO == "High-High"),
    sum(sf_map$LISA_RF == "High-High"),
    sum(sf_map$LISA_XGB == "High-High")
  ),
  
  LowLow = c(
    sum(sf_map$LISA_ANN == "Low-Low"),
    sum(sf_map$LISA_ANN_PSO == "Low-Low"),
    sum(sf_map$LISA_ANN_GWO == "Low-Low"),
    sum(sf_map$LISA_RF == "Low-Low"),
    sum(sf_map$LISA_XGB == "Low-Low")
  ),
  
  HighLow = c(
    sum(sf_map$LISA_ANN == "High-Low"),
    sum(sf_map$LISA_ANN_PSO == "High-Low"),
    sum(sf_map$LISA_ANN_GWO == "High-Low"),
    sum(sf_map$LISA_RF == "High-Low"),
    sum(sf_map$LISA_XGB == "High-Low")
  ),
  
  LowHigh = c(
    sum(sf_map$LISA_ANN == "Low-High"),
    sum(sf_map$LISA_ANN_PSO == "Low-High"),
    sum(sf_map$LISA_ANN_GWO == "Low-High"),
    sum(sf_map$LISA_RF == "Low-High"),
    sum(sf_map$LISA_XGB == "Low-High")
  ),
  
  NotSig = c(
    sum(sf_map$LISA_ANN == "Not significant"),
    sum(sf_map$LISA_ANN_PSO == "Not significant"),
    sum(sf_map$LISA_ANN_GWO == "Not significant"),
    sum(sf_map$LISA_RF == "Not significant"),
    sum(sf_map$LISA_XGB == "Not significant")
  )
) %>%
  mutate(
    N = rowSums(
      across(
        c(
          HighHigh,
          LowLow,
          HighLow,
          LowHigh,
          NotSig
        )
      )
    ),
    
    SigPct =
      round(
        100 * (N - NotSig) / N,
        2
      )
  )

write_csv(
  cluster_summary,
  file.path(
    out_dir,
    "LISA_Cluster_Counts_and_SignificantPercent.csv"
  )
)
