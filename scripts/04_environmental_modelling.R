## ============================================================================
## Supplementary R script
## Environmental modeling and temporal recurrence of canopy spectral diversity
## Sooretama Biological Reserve (REBIO Sooretama) and Vale Natural Reserve (RNV)


library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(mgcv)
library(ranger)
library(xgboost)
library(ggplot2)
library(patchwork)
library(spdep)

SEED <- 1
set.seed(SEED)

## ============================================================================
## 0.1) Input data, responses, predictors, and projected coordinates
## ============================================================================

DATA_FILE <- "df_RNV_July2026_clean.txt"

if (!file.exists(DATA_FILE)) {
  stop("Input file not found: ", DATA_FILE)
}

data_all <- read.table(
  DATA_FILE,
  header = TRUE, sep = "\t", dec = ".",
  stringsAsFactors = FALSE
)

# Projected coordinate columns are required for environmental fold records,
# spatial blocking, and Moran's I. If desired, replace NULL below with the exact
# two projected-coordinate column names, e.g. c("Easting", "Northing").
COORD_COLS <- NULL

# When COORD_COLS is left as NULL, common projected-coordinate names are detected.
coordinate_candidates <- list(
  c("x", "y"))
if (is.null(COORD_COLS)) {
  valid_coordinate_pairs <- Filter(
    function(z) all(z %in% names(data_all)),
    coordinate_candidates
  )
  if (length(valid_coordinate_pairs) == 0) {
    stop(
      "Projected coordinate columns could not be detected. ",
      "Set COORD_COLS <- c('your_x_column', 'your_y_column') near the top of the script."
    )
  }
  COORD_COLS <- valid_coordinate_pairs[[1]]
}
if (length(COORD_COLS) != 2 || !all(COORD_COLS %in% names(data_all))) {
  stop("COORD_COLS must contain exactly two existing projected-coordinate columns.")
}
message("Projected coordinate columns: ", paste(COORD_COLS, collapse = ", "))

indices_div <- c(
  "H_2024", "S_2024", "FRic_2024", 
  "FDiv_2024", "FDis_2024", "RaoQ_2024", "RaoRAW_2024"
)

predictors <- c(
  "aspectHorn_90m", "slopeHorn_90m",
  "TPI_90m", "TRIriley_90m", "roughness_90m",
  "flowdir_90m", "dist_to_river_m", "dist_to_edge_m", "DEM"
)


var_labels <- c(
  DEM = "Elevation (DEM)",
  aspectHorn_90m = "Aspect",
  slopeHorn_90m = "Slope",
  TPI_90m = "Topographic position (TPI)",
  TRIriley_90m = "Terrain ruggedness (TRI)",
  roughness_90m = "Roughness",
  flowdir_90m = "Flow direction",
  dist_to_edge_m = "Distance to edge",
  dist_to_river_m = "Distance to rivers"
)

var_labels <- c(
  DEM = "Elevation (DEM)",
  aspect = "Aspect",
  slope = "Slope",
  TPI = "Topographic position (TPI)",
  TRI = "Terrain ruggedness (TRI)",
  flowdir = "Flow direction",
  river_dist = "Distance to rivers",
  edge_dist = "Distance to edge",
  rough_abs = "Roughness",
)


res_labels <- c(
  H_2024 = "H'", S_2024 = "S", FRic_2024 = "FRic",
  FDiv_2024 = "FDiv", FDis_2024 = "FDis",
  RaoRAW_2024 = "RaoQ raw", RaoQ_2024 = "RaoQ sc"
)

## =========================================================
## 0.2) Collinearity among environmental predictors
## =========================================================

pred_mat <- data_all[, predictors, drop = FALSE]
cor_mat <- cor(pred_mat, use = "pairwise.complete.obs")

cor_df <- cor_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("var1") %>%
  tidyr::pivot_longer(cols = -var1, names_to = "var2", values_to = "r") %>%
  mutate(
    var1_lab = factor(var_labels[var1], levels = var_labels[predictors]),
    var2_lab = factor(var_labels[var2], levels = rev(var_labels[predictors]))
  )

fig_cor <- ggplot(cor_df, aes(x = var1_lab, y = var2_lab, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(
    name = "Pearson r", limits = c(-1, 1), midpoint = 0,
    low = "#2166ac", mid = "white", high = "#b2182b"
  ) +
  labs(x = "Predictors", y = "Predictors") +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

print(fig_cor)
ggsave("fig_correlation_predictors.png", fig_cor, width = 7, height = 6, dpi = 300)
write.csv(cor_mat, "correlation_matrix_predictors.csv")

calc_vif_values <- function(df, vars) {
  out <- vapply(vars, function(v) {
    others <- setdiff(vars, v)
    dat <- df[, c(v, others), drop = FALSE]
    dat <- dat[complete.cases(dat), , drop = FALSE]
    if (length(others) == 0 || nrow(dat) < length(others) + 3) return(NA_real_)
    fit <- lm(stats::reformulate(others, response = v), data = dat)
    r2 <- summary(fit)$r.squared
    if (!is.finite(r2) || r2 >= 1) return(Inf)
    1 / (1 - r2)
  }, numeric(1))
  tibble(Term = names(out), VIF = as.numeric(out))
}

vif_df <- calc_vif_values(data_all, predictors) %>%
  mutate(
    Term_plot = var_labels[Term],
    Term_plot = factor(Term_plot, levels = Term_plot[order(VIF)]),
    vif_class = cut(
      VIF, breaks = c(0, 5, 10, Inf),
      labels = c("< 5 (low)", "5-10 (moderate)", "> 10 (high)"),
      right = FALSE
    )
  )

fig_vif <- ggplot(vif_df, aes(x = VIF, y = Term_plot, fill = vif_class)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  scale_fill_manual(
    name = "VIF class",
    values = c(
      "< 5 (low)" = "#c7e9c0",
      "5-10 (moderate)" = "#fdd49e",
      "> 10 (high)" = "#fb6a4a"
    )
  ) +
  labs(x = "Variance inflation factor (VIF)", y = "Environmental predictors") +
  theme_bw()

print(fig_vif)
ggsave("fig_VIF_predictors.png", fig_vif, width = 7, height = 5, dpi = 300)
write.csv(vif_df, "VIF_predictors.csv", row.names = FALSE)

## =========================================================
## ENVIRONMENTALLY BLOCKED NESTED CROSS-VALIDATION (k-means) 
## -------------------------------
responses <- indices_div

preds_all <- c(
  "aspect",
  "slope",
  "TPI",
  "TRI",
  "flowdir",
  "river_dist",
  "edge_dist",
  "DEM",
  "rough_abs"
)



CENTERS  <- 50
K_OUTER  <- 5
K_INNER  <- 5

MIN_SUBSET <- 2

## 
TOPN      <- 8        # number of predictor subsets retained after Stage 1
STAGE1_N  <- 12000    # maximum number of observations used in Stage 1
K_STAGE1  <- 2        # environmentally blocked folds used in Stage 1

detected_cores <- parallel::detectCores()
NTHREADS <- if (is.na(detected_cores)) 1 else max(1, detected_cores - 2)

## Collinearity thresholds for candidate predictor subsets.
COR_MAX <- 0.75
VIF_MAX <- 5

## XGBoost configuration.
XGB_VALIDATION_FRACTION <- 0.10
XGB_MAX_NROUNDS <- 1200
XGB_EARLY_STOPPING <- 25

## Spatial diagnostics and spatial-block sensitivity settings.
MORAN_K <- 8
MORAN_NSIM <- 999
MORAN_MAX_N <- 5000
SPATIAL_BLOCK_SIZE <- NULL  # NULL: estimate block width to obtain approximately 50 blocks

CKPT_FILE <- "checkpoint_nested_env_blocking.rds"
LOG_DIR   <- "logs_nested_env_blocking_2024"
dir.create(LOG_DIR, showWarnings = FALSE)

UPDATE_EVERY <- 10     # logging interval during model selection


## -------------------------------
## 1) Helpers
## -------------------------------
log_msg <- function(msg, log_file){
  cat(msg, "\n")
  cat(paste0(msg, "\n"), file = log_file, append = TRUE)
  flush.console()
}

metrics <- function(obs, pred){
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]

  if(length(obs) < 3){
    return(tibble::tibble(
      R2_pred = NA_real_, nRMSE = NA_real_,
      Calibration_intercept = NA_real_, Calibration_slope = NA_real_,
      n = length(obs)
    ))
  }
  
  rng <- diff(range(obs))

  RMSE <- sqrt(mean((obs - pred)^2))
  nRMSE <- if(is.finite(rng) && rng > 0) 100 * RMSE / rng else NA_real_

  SSE <- sum((obs - pred)^2)
  SST <- sum((obs - mean(obs))^2)
  R2_pred <- if(is.finite(SST) && SST > 0) 1 - SSE/SST else NA_real_

  cal <- lm(obs ~ pred)
  cal_coef <- coef(cal)

  tibble::tibble(
    R2_pred = R2_pred,
    nRMSE = nRMSE,
    Calibration_intercept = unname(cal_coef[1]),
    Calibration_slope = unname(cal_coef[2]),
    n = length(obs)
  )
}

all_subsets <- function(preds, min_size = 2){
  subs <- unlist(
    lapply(seq_along(preds), function(m) combn(preds, m, simplify = FALSE)),
    recursive = FALSE
  )
  Filter(function(s) length(s) >= min_size, subs)
}

subset_is_valid <- function(df, vars, cor_max = COR_MAX, vif_max = VIF_MAX){
  dat <- df[, vars, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  if(nrow(dat) < length(vars) + 3) return(FALSE)

  cm <- suppressWarnings(cor(dat))
  upper <- abs(cm[upper.tri(cm)])
  if(length(upper) > 0 && any(upper >= cor_max, na.rm = TRUE)) return(FALSE)

  vifs <- calc_vif_values(dat, vars)$VIF
  if(any(!is.finite(vifs)) || any(vifs >= vif_max, na.rm = TRUE)) return(FALSE)
  TRUE
}

filter_subsets_collinearity <- function(df, subsets,
                                        cor_max = COR_MAX,
                                        vif_max = VIF_MAX){
  keep <- vapply(
    subsets,
    function(x) subset_is_valid(df, x, cor_max = cor_max, vif_max = vif_max),
    logical(1)
  )
  subsets[keep]
}

## -------------------------------
## 2) Folds por GRUPO (clusters) balanceados
## -------------------------------
make_group_folds <- function(group_id, K = 5, seed = 1){
  set.seed(seed)
  g <- as.character(group_id)
  tab <- sort(table(g), decreasing = TRUE)
  groups <- names(tab)
  sizes  <- as.integer(tab)

  fold_sizes <- rep(0, K)
  group_to_fold <- integer(length(groups))
  names(group_to_fold) <- groups

  ## Randomize tied group sizes using the fixed seed.
  size_levels <- unique(sizes)
  ordered_groups <- unlist(lapply(size_levels, function(sz){
    same <- groups[sizes == sz]
    sample(same, length(same))
  }), use.names = FALSE)

  for(i in seq_along(ordered_groups)){
    f <- which.min(fold_sizes)
    group_to_fold[ordered_groups[i]] <- f
    fold_sizes[f] <- fold_sizes[f] + as.integer(tab[ordered_groups[i]])
  }
  as.integer(group_to_fold[g])
}


make_spatial_blocks <- function(coords, K = 5, target_blocks = 50,
                                block_size = NULL, seed = 1){
  coords <- as.matrix(coords)
  if(ncol(coords) != 2) stop("coords must contain exactly two columns.")

  rx <- diff(range(coords[,1], na.rm = TRUE))
  ry <- diff(range(coords[,2], na.rm = TRUE))
  if(is.null(block_size)){
    area_bbox <- rx * ry
    if(!is.finite(area_bbox) || area_bbox <= 0) stop("Invalid projected coordinates.")
    block_size <- sqrt(area_bbox / target_blocks)
  }

  bx <- floor((coords[,1] - min(coords[,1], na.rm = TRUE)) / block_size)
  by <- floor((coords[,2] - min(coords[,2], na.rm = TRUE)) / block_size)
  block_id <- paste(bx, by, sep = "_")
  if(length(unique(block_id)) < K){
    stop("The number of spatial blocks is smaller than the number of folds.")
  }

  list(
    block_id = block_id,
    fold = make_group_folds(block_id, K = K, seed = seed),
    block_size = block_size,
    n_blocks = length(unique(block_id))
  )
}

## -------------------------------
## 2.1) Stage 1 sampling while preserving fold representation
## -------------------------------
sample_by_fold <- function(df, fold_id, n_total = 12000, seed = 1){
  if(is.null(n_total) || n_total <= 0 || n_total >= nrow(df)){
    return(list(df = df, idx = seq_len(nrow(df))))
  }
  set.seed(seed)
  K <- max(fold_id)
  per_fold <- max(1, floor(n_total / K))
  idx_keep <- integer(0)
  
  for(f in 1:K){
    idx_f <- which(fold_id == f)
    if(length(idx_f) == 0) next
    take <- min(per_fold, length(idx_f))
    idx_keep <- c(idx_keep, sample(idx_f, take))
  }
  list(df = df[idx_keep, , drop = FALSE], idx = idx_keep)
}

## -------------------------------
## 3) Model fits (train -> predict)
## -------------------------------
## 3.1) GAM via BAM 
fit_pred_gam <- function(train, test, resp, preds, k = 6, gamma = 1.0, nthreads = 1){
  form <- as.formula(
    paste(resp, "~", paste(paste0("s(", preds, ", bs='tp', k=", k, ")"), collapse = " + "))
  )
  m <- mgcv::bam(
    form, data = train,
    method = "fREML",
    select = TRUE,
    gamma = gamma,
    discrete = TRUE,
    nthreads = nthreads
  )
  as.numeric(predict(m, newdata = test))
}

fit_pred_rf <- function(train, test, resp, preds, mtry, min_node = 10, ntree = 400,
                        seed = 1, nthreads = 1){
  m <- ranger::ranger(
    formula = as.formula(paste(resp, "~ .")),
    data = train[, c(resp, preds), drop = FALSE],
    num.trees = ntree,
    mtry = mtry,
    min.node.size = min_node,
    importance = "permutation",
    seed = seed,
    num.threads = nthreads
  )
  as.numeric(predict(m, data = test[, preds, drop = FALSE])$predictions)
}

predict_xgb_best <- function(model, newdata){
  best_it <- model$best_iteration
  if(is.null(best_it) || !is.finite(best_it)){
    return(as.numeric(predict(model, newdata = newdata)))
  }

  tryCatch(
    as.numeric(predict(
      model,
      newdata = newdata,
      iterationrange = c(1L, as.integer(best_it) + 1L)
    )),
    error = function(e){
      
      as.numeric(predict(
        model,
        newdata = newdata,
        ntreelimit = as.integer(best_it)
      ))
    }
  )
}

fit_pred_xgb <- function(train, test, resp, preds, params,
                         nrounds = XGB_MAX_NROUNDS,
                         early = XGB_EARLY_STOPPING,
                         seed = 1,
                         validation_fraction = XGB_VALIDATION_FRACTION){
  set.seed(seed)
  n <- nrow(train)
  if(n < 20) stop("Insufficient sample size for the internal XGBoost validation split.")

  n_val <- max(10, floor(n * validation_fraction))
  n_val <- min(n_val, n - 10)
  idx_val <- sample(seq_len(n), n_val)

  train_fit <- train[-idx_val, , drop = FALSE]
  train_val <- train[idx_val, , drop = FALSE]

  dtr <- xgb.DMatrix(
    data = as.matrix(train_fit[, preds, drop = FALSE]),
    label = train_fit[[resp]]
  )
  dval <- xgb.DMatrix(
    data = as.matrix(train_val[, preds, drop = FALSE]),
    label = train_val[[resp]]
  )
  dte <- xgb.DMatrix(data = as.matrix(test[, preds, drop = FALSE]))

  m <- xgb.train(
    params = params,
    data = dtr,
    nrounds = nrounds,
    watchlist = list(train = dtr, validation = dval),
    early_stopping_rounds = early,
    maximize = FALSE,
    verbose = 0
  )

  pred <- predict_xgb_best(m, dte)
  attr(pred, "best_iteration") <- m$best_iteration
  pred
}

## -------------------------------
## 4) Inner CV: RMSE médio
## -------------------------------
cv_rmse_gam <- function(df, resp, preds, fold_id, k, gamma, nthreads = 1){
  K <- max(fold_id)
  rmses <- numeric(K)
  for(f in 1:K){
    tr <- df[fold_id != f, , drop = FALSE]
    te <- df[fold_id == f, , drop = FALSE]
    pr <- fit_pred_gam(tr, te, resp, preds, k = k, gamma = gamma, nthreads = nthreads)
    rmses[f] <- sqrt(mean((te[[resp]] - pr)^2))
  }
  mean(rmses)
}

cv_rmse_rf <- function(df, resp, preds, fold_id, mtry, min_node, ntree, seed, nthreads){
  K <- max(fold_id)
  rmses <- numeric(K)
  for(f in 1:K){
    tr <- df[fold_id != f, , drop = FALSE]
    te <- df[fold_id == f, , drop = FALSE]
    pr <- fit_pred_rf(tr, te, resp, preds, mtry = mtry, min_node = min_node,
                      ntree = ntree, seed = seed, nthreads = nthreads)
    rmses[f] <- sqrt(mean((te[[resp]] - pr)^2))
  }
  mean(rmses)
}

cv_rmse_xgb <- function(df, resp, preds, fold_id, params, seed,
                        nrounds = 1200, early = 25){
  K <- max(fold_id)
  rmses <- numeric(K)
  for(f in 1:K){
    tr <- df[fold_id != f, , drop = FALSE]
    te <- df[fold_id == f, , drop = FALSE]
    pr <- fit_pred_xgb(tr, te, resp, preds, params = params, seed = seed,
                       nrounds = nrounds, early = early)
    rmses[f] <- sqrt(mean((te[[resp]] - pr)^2))
  }
  mean(rmses)
}

## -------------------------------
## 5) Inner cross-validation: two-stage model selection
## -------------------------------
select_best_pipeline <- function(train_df, resp, preds_all, inner_fold,
                                 cluster_env,           # environmental-cluster labels for the training data
                                 nthreads = 1, seed = 1,
                                 min_subset = 2,
                                 log_file = NULL,
                                 update_every = 10,
                                 topN = 8,
                                 stage1_n = 12000,
                                 k_stage1 = 2){
  
  subs_all <- all_subsets(preds_all, min_size = min_subset)
  subs <- filter_subsets_collinearity(train_df, subs_all)
  if(length(subs) == 0) stop("No predictor subset passed the collinearity criteria.")
  if(length(subs) < topN) topN <- length(subs)

  if(!is.null(log_file)){
    log_msg(sprintf(
      "  [COLLINEARITY] %d/%d subsets retained | |r| < %.2f | VIF < %.1f",
      length(subs), length(subs_all), COR_MAX, VIF_MAX
    ), log_file)
  }
  
  ## ----------------
  ## STAGE 1: reduced sample and two environmental folds
  ## ----------------
  smp <- sample_by_fold(train_df, inner_fold, n_total = stage1_n, seed = seed)
  train_s1 <- smp$df
  cluster_s1 <- cluster_env[smp$idx]
  
  fold_s1 <- make_group_folds(cluster_s1, K = k_stage1, seed = seed)
  
  if(!is.null(log_file)){
    log_msg(sprintf("  [INNER] Stage1: n=%d/%d | K_stage1=%d | topN=%d",
                    nrow(train_s1), nrow(train_df), k_stage1, topN), log_file)
  }
  
  ## Stage 2 hyperparameter grids (full inner cross-validation)
  grid_gam2 <- expand.grid(k = c(5, 6, 8), gamma = c(1.0, 1.4))
  grid_rf2  <- expand.grid(min_node = c(5, 10), ntree = c(500))
  grid_xgb2 <- expand.grid(
    max_depth = c(3, 6),
    eta = c(0.05, 0.1),
    subsample = c(0.7, 1.0),
    colsample_bytree = c(0.7, 1.0)
  )
  
  best <- list(GAM = NULL, RF = NULL, XGB = NULL)
  
  ## =========================================================
  ## STAGE 1 (COARSE)
  ## =========================================================
  
  ## ---- GAM coarse screening
  if(!is.null(log_file)) log_msg(sprintf("  [INNER][STAGE1] GAM: %d subsets", length(subs)), log_file)
  gam_s1 <- vector("list", length(subs))
  for(i in seq_along(subs)){
    s <- subs[[i]]
    rmse <- cv_rmse_gam(train_s1, resp, s, fold_s1, k = 6, gamma = 1.0, nthreads = nthreads)
    gam_s1[[i]] <- list(preds = s, rmse = rmse)
    if(!is.null(log_file) && i %% update_every == 0){
      log_msg(sprintf("    STAGE1 GAM: %.1f%% (%d/%d)", 100*i/length(subs), i, length(subs)), log_file)
    }
  }
  ord <- order(vapply(gam_s1, `[[`, numeric(1), "rmse"))
  gam_top <- gam_s1[ord[1:topN]]
  if(!is.null(log_file)) log_msg(sprintf("  [STAGE1] GAM top RMSE=%.6f", gam_top[[1]]$rmse), log_file)
  
  ## ---- RF coarse
  if(!is.null(log_file)) log_msg(sprintf("  [INNER][STAGE1] RF: %d subsets", length(subs)), log_file)
  rf_s1 <- vector("list", length(subs))
  for(i in seq_along(subs)){
    s <- subs[[i]]
    p <- length(s)
    mtry0 <- max(1, floor(p/3))
    rmse <- cv_rmse_rf(train_s1, resp, s, fold_s1,
                       mtry = mtry0, min_node = 15, ntree = 200,
                       seed = seed, nthreads = nthreads)
    rf_s1[[i]] <- list(preds = s, rmse = rmse)
    if(!is.null(log_file) && i %% update_every == 0){
      log_msg(sprintf("    STAGE1 RF: %.1f%% (%d/%d)", 100*i/length(subs), i, length(subs)), log_file)
    }
  }
  ord <- order(vapply(rf_s1, `[[`, numeric(1), "rmse"))
  rf_top <- rf_s1[ord[1:topN]]
  if(!is.null(log_file)) log_msg(sprintf("  [STAGE1] RF top RMSE=%.6f", rf_top[[1]]$rmse), log_file)
  
  ## ---- XGB coarse
  if(!is.null(log_file)) log_msg(sprintf("  [INNER][STAGE1] XGB: %d subsets", length(subs)), log_file)
  xgb_s1 <- vector("list", length(subs))
  params0 <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    max_depth = 3,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8,
    nthread = nthreads
  )
  for(i in seq_along(subs)){
    s <- subs[[i]]
    rmse <- cv_rmse_xgb(train_s1, resp, s, fold_s1, params = params0, seed = seed,
                        nrounds = 250, early = 15)
    xgb_s1[[i]] <- list(preds = s, rmse = rmse)
    if(!is.null(log_file) && i %% update_every == 0){
      log_msg(sprintf("    STAGE1 XGB: %.1f%% (%d/%d)", 100*i/length(subs), i, length(subs)), log_file)
    }
  }
  ord <- order(vapply(xgb_s1, `[[`, numeric(1), "rmse"))
  xgb_top <- xgb_s1[ord[1:topN]]
  if(!is.null(log_file)) log_msg(sprintf("  [STAGE1] XGB top RMSE=%.6f", xgb_top[[1]]$rmse), log_file)
  
  ## =========================================================
  ## STAGE 2: refinement using the full outer-training data and K_INNER folds
  ## =========================================================
  
  if(!is.null(log_file)) log_msg(sprintf("  [INNER][STAGE2] GAM: %d subsets x %d configs",
                                         length(gam_top), nrow(grid_gam2)), log_file)
  best_rmse <- Inf
  done <- 0
  total <- length(gam_top) * nrow(grid_gam2)
  for(obj in gam_top){
    s <- obj$preds
    for(i in seq_len(nrow(grid_gam2))){
      rmse <- cv_rmse_gam(train_df, resp, s, inner_fold,
                          k = grid_gam2$k[i], gamma = grid_gam2$gamma[i], nthreads = nthreads)
      done <- done + 1
      if(!is.null(log_file) && done %% update_every == 0){
        log_msg(sprintf("    STAGE2 GAM: %.1f%% (%d/%d)", 100*done/total, done, total), log_file)
      }
      if(rmse < best_rmse){
        best_rmse <- rmse
        best$GAM <- list(preds = s, k = grid_gam2$k[i], gamma = grid_gam2$gamma[i], rmse = rmse, stage = "two-stage")
      }
    }
  }
  if(!is.null(log_file)) log_msg(sprintf("  [INNER] GAM best RMSE=%.6f | vars=%s | k=%d gamma=%.1f",
                                         best$GAM$rmse, paste(best$GAM$preds, collapse=", "),
                                         best$GAM$k, best$GAM$gamma), log_file)
  
  ## ---- RF refine
  if(!is.null(log_file)) log_msg(sprintf("  [INNER][STAGE2] RF: %d subsets x (mtry x %d configs)",
                                         length(rf_top), nrow(grid_rf2)), log_file)
  best_rmse <- Inf
  done <- 0
  total <- sum(vapply(rf_top, function(obj) length(obj$preds) * nrow(grid_rf2), numeric(1)))
  for(obj in rf_top){
    s <- obj$preds
    p <- length(s)
    for(i in seq_len(nrow(grid_rf2))){
      for(mtry in seq_len(p)){
        rmse <- cv_rmse_rf(train_df, resp, s, inner_fold,
                           mtry = mtry,
                           min_node = grid_rf2$min_node[i],
                           ntree = grid_rf2$ntree[i],
                           seed = seed,
                           nthreads = nthreads)
        done <- done + 1
        if(!is.null(log_file) && done %% update_every == 0){
          log_msg(sprintf("    STAGE2 RF: %.1f%% (%d/%d)", 100*done/total, done, total), log_file)
        }
        if(rmse < best_rmse){
          best_rmse <- rmse
          best$RF <- list(preds = s, mtry = mtry, min_node = grid_rf2$min_node[i],
                          ntree = grid_rf2$ntree[i], rmse = rmse, stage = "two-stage")
        }
      }
    }
  }
  if(!is.null(log_file)) log_msg(sprintf("  [INNER] RF best RMSE=%.6f | vars=%s | mtry=%d min_node=%d ntree=%d",
                                         best$RF$rmse, paste(best$RF$preds, collapse=", "),
                                         best$RF$mtry, best$RF$min_node, best$RF$ntree), log_file)
  
  ## ---- XGB refine
  if(!is.null(log_file)) log_msg(sprintf("  [INNER][STAGE2] XGB: %d subsets x %d configs",
                                         length(xgb_top), nrow(grid_xgb2)), log_file)
  best_rmse <- Inf
  done <- 0
  total <- length(xgb_top) * nrow(grid_xgb2)
  for(obj in xgb_top){
    s <- obj$preds
    for(i in seq_len(nrow(grid_xgb2))){
      params <- list(
        objective = "reg:squarederror",
        eval_metric = "rmse",
        max_depth = grid_xgb2$max_depth[i],
        eta = grid_xgb2$eta[i],
        subsample = grid_xgb2$subsample[i],
        colsample_bytree = grid_xgb2$colsample_bytree[i],
        nthread = nthreads
      )
      rmse <- cv_rmse_xgb(train_df, resp, s, inner_fold,
                          params = params, seed = seed,
                          nrounds = 1200, early = 25)
      done <- done + 1
      if(!is.null(log_file) && done %% update_every == 0){
        log_msg(sprintf("    STAGE2 XGB: %.1f%% (%d/%d)", 100*done/total, done, total), log_file)
      }
      if(rmse < best_rmse){
        best_rmse <- rmse
        best$XGB <- list(preds = s, params = params, nrounds = XGB_MAX_NROUNDS, early = XGB_EARLY_STOPPING, validation_fraction = XGB_VALIDATION_FRACTION, rmse = rmse, stage = "two-stage")
      }
    }
  }
  if(!is.null(log_file)) log_msg(sprintf("  [INNER] XGB best RMSE=%.6f | vars=%s | depth=%d eta=%.2f subs=%.1f col=%.1f",
                                         best$XGB$rmse, paste(best$XGB$preds, collapse=", "),
                                         best$XGB$params$max_depth, best$XGB$params$eta,
                                         best$XGB$params$subsample, best$XGB$params$colsample_bytree), log_file)
  
  best$stage1_top <- list(GAM = gam_top, RF = rf_top, XGB = xgb_top)
  best$stage1_n_used <- nrow(train_s1)
  best$k_stage1_used <- k_stage1
  
  best
}

## -------------------------------
## 6) Outer cross-validation runner with fold-level checkpoints
## -------------------------------
run_nested_env_blocking_all <- function(df, responses, preds_all, coord_cols,
                                        centers = 50, K_outer = 5, K_inner = 5,
                                        min_subset = 2,
                                        seed = 1, nthreads = 1,
                                        ckpt_file = "checkpoint_nested_env_blocking.rds",
                                        log_dir = "logs_nested_env_blocking",
                                        update_every = 10,
                                        topN = 8,
                                        stage1_n = 12000,
                                        k_stage1 = 2){
  
  if(file.exists(ckpt_file)){
    ckpt <- readRDS(ckpt_file)
    results_acc <- ckpt$results
    done_map <- ckpt$done_map
    best_acc <- ckpt$best_pipelines
    pred_acc <- ckpt$oof_predictions
    folds_acc <- ckpt$fold_assignments
  } else {
    results_acc <- list()
    done_map <- list()
    best_acc <- list()
    pred_acc <- list()
    folds_acc <- list()
  }
  
  total_tasks <- length(responses) * K_outer
  done_tasks <- sum(vapply(done_map, function(x) sum(unlist(x)), integer(1), USE.NAMES = FALSE), na.rm = TRUE)
  
  pb_total <- utils::txtProgressBar(min = 0, max = total_tasks, style = 3)
  utils::setTxtProgressBar(pb_total, done_tasks)
  
  for(r in seq_along(responses)){
    resp <- responses[r]
    log_file <- file.path(log_dir, paste0("NESTED_", resp, ".log"))
    
    log_msg("====================================================", log_file)
    log_msg(sprintf("[RESP %d/%d] %s", r, length(responses), resp), log_file)
    log_msg(sprintf("centers=%d | K_outer=%d | K_inner=%d | min_subset=%d | nthreads=%d | topN=%d | stage1_n=%d | K_stage1=%d",
                    centers, K_outer, K_inner, min_subset, nthreads, topN, stage1_n, k_stage1), log_file)
    log_msg(paste0("Start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), log_file)
    
    needed <- c(resp, preds_all, coord_cols)
    miss <- setdiff(needed, names(df))
    if(length(miss) > 0) stop("Missing columns for ", resp, ": ", paste(miss, collapse=", "))
    
    complete_idx <- which(stats::complete.cases(df[, needed, drop = FALSE]))
    df0 <- df[complete_idx, needed, drop = FALSE]
    df0 <- tibble::as_tibble(df0)
    df0$row_id_original <- complete_idx
    
    ## Environmental clustering
    X <- scale(df0[, preds_all])
    set.seed(seed)
    km <- kmeans(X, centers = centers, nstart = 10)
    df0$cluster_env <- km$cluster
    
    ## Assign complete environmental clusters to outer folds
    df0$fold_outer <- make_group_folds(df0$cluster_env, K = K_outer, seed = seed)

    folds_acc[[resp]] <- df0 %>%
      dplyr::select(row_id_original, dplyr::all_of(coord_cols), cluster_env, fold_outer)
    
    if(is.null(done_map[[resp]])) done_map[[resp]] <- as.list(rep(FALSE, K_outer))
    if(is.null(results_acc[[resp]])) results_acc[[resp]] <- vector("list", K_outer)
    if(is.null(best_acc[[resp]])) best_acc[[resp]] <- vector("list", K_outer)
    if(is.null(pred_acc[[resp]])) pred_acc[[resp]] <- vector("list", K_outer)
    
    for(f in 1:K_outer){
      
      if(isTRUE(done_map[[resp]][[f]])){
        log_msg(sprintf("[CHECKPOINT] Skipping %s | outer fold %d/%d (already completed)", resp, f, K_outer), log_file)
        next
      }
      
      log_msg(sprintf("---- OUTER fold %d/%d ----", f, K_outer), log_file)
      cat(sprintf("\n[%s] OUTER fold %d/%d\n", resp, f, K_OUTER)); flush.console()
      
      train_df <- dplyr::filter(df0, fold_outer != f)
      test_df  <- dplyr::filter(df0, fold_outer == f)
      
      ## Reassign training-set environmental clusters to inner folds
      train_df$fold_inner <- make_group_folds(train_df$cluster_env, K = K_inner, seed = seed + f)
      
      ## Two-stage predictor-subset and hyperparameter selection
      best <- select_best_pipeline(
        train_df, resp, preds_all,
        inner_fold = train_df$fold_inner,
        cluster_env = train_df$cluster_env,
        nthreads = nthreads,
        seed = seed + f,
        min_subset = min_subset,
        log_file = log_file,
        update_every = update_every,
        topN = topN,
        stage1_n = stage1_n,
        k_stage1 = k_stage1
      )
      
      ## Evaluate the selected configurations on the withheld outer fold
      pr_gam <- fit_pred_gam(train_df, test_df, resp, best$GAM$preds,
                             k = best$GAM$k, gamma = best$GAM$gamma, nthreads = nthreads)
      met_gam <- metrics(test_df[[resp]], pr_gam) %>%
        dplyr::mutate(Model="GAM", Response=resp, Fold=f, Validation="Environmental",
                      Vars=paste(best$GAM$preds, collapse=", "),
                      Params=paste0("k=",best$GAM$k,"; gamma=",best$GAM$gamma),
                      Inner_RMSE = best$GAM$rmse,
                      Stage = best$GAM$stage,
                      Stage1_n = best$stage1_n_used,
                      K_stage1 = best$k_stage1_used)
      
      pr_rf <- fit_pred_rf(train_df, test_df, resp, best$RF$preds,
                           mtry = best$RF$mtry, min_node = best$RF$min_node,
                           ntree = best$RF$ntree, seed = seed + f, nthreads = nthreads)
      met_rf <- metrics(test_df[[resp]], pr_rf) %>%
        dplyr::mutate(Model="RF", Response=resp, Fold=f, Validation="Environmental",
                      Vars=paste(best$RF$preds, collapse=", "),
                      Params=paste0("mtry=",best$RF$mtry,"; min_node=",best$RF$min_node,"; ntree=",best$RF$ntree),
                      Inner_RMSE = best$RF$rmse,
                      Stage = best$RF$stage,
                      Stage1_n = best$stage1_n_used,
                      K_stage1 = best$k_stage1_used)
      
      pr_xgb <- fit_pred_xgb(train_df, test_df, resp, best$XGB$preds,
                             params = best$XGB$params, seed = seed + f,
                             nrounds = best$XGB$nrounds,
                             early = best$XGB$early,
                             validation_fraction = best$XGB$validation_fraction)
      xgb_best_iteration <- attr(pr_xgb, "best_iteration")
      met_xgb <- metrics(test_df[[resp]], pr_xgb) %>%
        dplyr::mutate(Model="XGB", Response=resp, Fold=f, Validation="Environmental",
                      Vars=paste(best$XGB$preds, collapse=", "),
                      Params=paste0("depth=",best$XGB$params$max_depth,
                                    "; eta=",best$XGB$params$eta,
                                    "; subs=",best$XGB$params$subsample,
                                    "; col=",best$XGB$params$colsample_bytree,
                                    "; max_nrounds=", best$XGB$nrounds,
                                    "; early=", best$XGB$early,
                                    "; best_iteration=", xgb_best_iteration),
                      Inner_RMSE = best$XGB$rmse,
                      Stage = best$XGB$stage,
                      Stage1_n = best$stage1_n_used,
                      K_stage1 = best$k_stage1_used)
      
      results_acc[[resp]][[f]] <- dplyr::bind_rows(met_gam, met_rf, met_xgb)
      

      pred_acc[[resp]][[f]] <- dplyr::bind_rows(
        tibble::tibble(
          Response=resp, Model="GAM", Fold=f, Validation="Environmental",
          Observed=test_df[[resp]], Predicted=pr_gam,
          Residual=test_df[[resp]]-pr_gam,
          row_id_original=test_df$row_id_original,
          x=test_df[[coord_cols[1]]], y=test_df[[coord_cols[2]]]
        ),
        tibble::tibble(
          Response=resp, Model="RF", Fold=f, Validation="Environmental",
          Observed=test_df[[resp]], Predicted=pr_rf,
          Residual=test_df[[resp]]-pr_rf,
          row_id_original=test_df$row_id_original,
          x=test_df[[coord_cols[1]]], y=test_df[[coord_cols[2]]]
        ),
        tibble::tibble(
          Response=resp, Model="XGB", Fold=f, Validation="Environmental",
          Observed=test_df[[resp]], Predicted=as.numeric(pr_xgb),
          Residual=test_df[[resp]]-as.numeric(pr_xgb),
          row_id_original=test_df$row_id_original,
          x=test_df[[coord_cols[1]]], y=test_df[[coord_cols[2]]]
        )
      )
      
      ## Store the selected pipeline for later summaries and refitting
      best_acc[[resp]][[f]] <- list(
        resp = resp,
        fold = f,
        n_train = nrow(train_df),
        n_test  = nrow(test_df),
        best = best
      )
      
      done_map[[resp]][[f]] <- TRUE
      
      ## Save a checkpoint after each outer fold
      saveRDS(list(results = results_acc,
                   done_map = done_map,
                   best_pipelines = best_acc,
                   oof_predictions = pred_acc,
                   fold_assignments = folds_acc), ckpt_file)
      
      done_tasks <- done_tasks + 1
      utils::setTxtProgressBar(pb_total, done_tasks)
      pct_total <- 100 * done_tasks / total_tasks
      cat(sprintf("Overall progress: %.1f%% (%d/%d)\n", pct_total, done_tasks, total_tasks)); flush.console()
      log_msg(sprintf("Fold %d completed. Overall progress: %.1f%% (%d/%d)",
                      f, pct_total, done_tasks, total_tasks), log_file)
    }
    
    log_msg(paste0("End: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), log_file)
  }
  
  close(pb_total)
  
  res_nested <- dplyr::bind_rows(lapply(results_acc, function(lst) dplyr::bind_rows(lst)))
  
  list(
    res_nested = res_nested,
    best_pipelines = best_acc,
    oof_predictions = dplyr::bind_rows(lapply(pred_acc, dplyr::bind_rows)),
    fold_assignments = folds_acc
  )
}

## -------------------------------
## 7) Run environmentally blocked nested cross-validation
## -------------------------------
out <- run_nested_env_blocking_all(
  df = data_all,
  responses = responses,
  preds_all = preds_all,
  coord_cols = COORD_COLS,
  centers = CENTERS,
  K_outer = K_OUTER,
  K_inner = K_INNER,
  min_subset = MIN_SUBSET,
  seed = SEED,
  nthreads = NTHREADS,
  ckpt_file = CKPT_FILE,
  log_dir = LOG_DIR,
  update_every = UPDATE_EVERY,
  topN = TOPN,
  stage1_n = STAGE1_N,
  k_stage1 = K_STAGE1
)

res_nested <- out$res_nested
best_pipelines <- out$best_pipelines
oof_predictions_env <- out$oof_predictions

## -------------------------------
## 8) Performance summaries and best-performing models
## -------------------------------
summary_by_model <- res_nested %>%
  dplyr::group_by(Response, Model, Validation) %>%
  dplyr::summarise(
    R2_pred_mean = mean(R2_pred, na.rm=TRUE), R2_pred_sd = sd(R2_pred, na.rm=TRUE),
    nRMSE_mean = mean(nRMSE, na.rm=TRUE), nRMSE_sd = sd(nRMSE, na.rm=TRUE),
    Calibration_intercept_mean = mean(Calibration_intercept, na.rm=TRUE),
    Calibration_intercept_sd = sd(Calibration_intercept, na.rm=TRUE),
    Calibration_slope_mean = mean(Calibration_slope, na.rm=TRUE),
    Calibration_slope_sd = sd(Calibration_slope, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Response, RMSE_mean)

winner <- summary_by_model %>%
  dplyr::group_by(Response) %>%
  dplyr::arrange(nRMSE_mean, dplyr::desc(R2_pred_mean)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

print(summary_by_model)
print(winner)

# Complete performance diagnostics for supplementary reporting. Predictive R2
# and nRMSE are the primary manuscript metrics; the remaining statistics are
# retained as complementary diagnostics.
complete_performance_table_2024 <- summary_by_model %>%
  dplyr::filter(Validation == "Environmental") %>%
  dplyr::mutate(
    Index = unname(res_labels[Response]),
    Model = Model,
    `Predictive R2 (mean +/- SD)` = sprintf("%.3f +/- %.3f", R2_pred_mean, R2_pred_sd),
    `nRMSE (%) (mean +/- SD)` = sprintf("%.3f +/- %.3f", nRMSE_mean, nRMSE_sd),
    `Calibration intercept (mean +/- SD)` = sprintf(
      "%.3f +/- %.3f", Calibration_intercept_mean, Calibration_intercept_sd
    ),
    `Calibration slope (mean +/- SD)` = sprintf(
      "%.3f +/- %.3f", Calibration_slope_mean, Calibration_slope_sd
    )
  ) %>%
  dplyr::select(
    Index, Model,
    `Predictive R2 (mean +/- SD)`,
    `Squared Pearson r2 (mean +/- SD)`,
    `RMSE (mean +/- SD)`,
    `nRMSE (%) (mean +/- SD)`,
    `MAE (mean +/- SD)`,
    `Bias (mean +/- SD)`,
    `Calibration intercept (mean +/- SD)`,
    `Calibration slope (mean +/- SD)`
  )

write.csv(
  complete_performance_table_2024,
  "complete_model_performance_diagnostics_2024.csv",
  row.names = FALSE
)

write.csv(res_nested, "nested_env_blocking_fold_results_all_indices2024.csv", row.names = FALSE)
write.csv(summary_by_model, "nested_env_blocking_summary_all_indices2024.csv", row.names = FALSE)
write.csv(winner, "nested_env_blocking_winner_all_indices.csv", row.names = FALSE)

saveRDS(best_pipelines, "best_pipelines_by_response_and_fold2024.rds")
write.csv(oof_predictions_env, "nested_env_blocking_OOF_predictions_2024.csv", row.names = FALSE)
saveRDS(out$fold_assignments, "environmental_fold_assignments_2024.rds")

top_sets <- res_nested %>%
  dplyr::count(Response, Model, Vars, Params, sort = TRUE)

head(top_sets, 20)

## =========================================================
## Figures and tables derived from the nested cross-validation results
## =========================================================

res_labels <- c(
  H_2024 = "H'",
  S_2024 = "S",
  FRic_2024 = "FRic",
  FDiv_2024 = "FDiv",
  FDis_2024 = "FDis",
  RaoRAW_2024 = "RaoQ raw",
  RaoQ_2024 = "RaoQ sc"
)

plot_dat <- res_nested %>%
  filter(Response %in% names(res_labels)) %>%
  mutate(
    Response = as.character(Response),
    Index = unname(res_labels[Response]),
    Index = factor(Index, levels = unname(res_labels)),
    Model = factor(Model, levels = c("GAM", "RF", "XGB"))
  )

model_cols <- c(
  GAM = "#66c2a5",
  RF  = "#1b9e77",
  XGB = "#7570b3"
)

p_r2pred <- ggplot(plot_dat, aes(x = Index, y = R2_pred, fill = Model, color = Model)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.55,
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.75),
    size = 1.8,
    alpha = 0.7
  ) +
  scale_fill_manual(values = model_cols) +
  scale_color_manual(values = model_cols) +
  labs(
    x = "Spectral diversity index",
    y = expression(R[pred]^2),
    fill = "Model",
    color = "Model"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_nrmse <- ggplot(plot_dat, aes(x = Index, y = nRMSE, fill = Model, color = Model)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.55,
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.75),
    size = 1.8,
    alpha = 0.7
  ) +
  scale_fill_manual(values = model_cols) +
  scale_color_manual(values = model_cols) +
  labs(
    x = "Spectral diversity index",
    y = "Normalized RMSE (%)",
    fill = "Model",
    color = "Model"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

fig_2024 <- p_r2pred / p_nrmse +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "top")

print(fig_2024)

ggsave(
  "model_performance_predictive_R2_nRMSE_2024.png",
  fig_2024,
  width = 10,
  height = 6.5,
  dpi = 300
)


#######################################
################ model-performance table


# Labels
res_labels <- c(
  H_2024        = "H'",
  S_2024        = "S",
  FRic_2024     = "FRic",
  FDiv_2024     = "FDiv",
  FDis_2024     = "FDis",
  RaoRAW_2024 = "RaoQ raw",
  RaoQ_2024  = "RaoQ sc"
)

resp_order <- c("S", "RaoQ sc", "RaoQ raw", "H'", "FRic", "FDiv", "FDis")

predictor_labels <- c(
  DEM             = "Elevation (DEM)",
  aspectHorn_90m  = "Aspect",
  slopeHorn_90m   = "Slope",
  TPI_90m         = "Topographic position (TPI)",
  TRIriley_90m    = "Terrain ruggedness (TRI)",
  roughness_90m   = "Roughness",
  flowdir_90m     = "Flow direction",
  dist_to_edge_m  = "Distance to edge",
  dist_to_river_m = "Distance to rivers"
)

# -----------------------------
# 1) Extract selected predictors by fold and algorithm
# -----------------------------
best_2024 <- best_pipelines[grepl("_2024$", names(best_pipelines))]

extract_best_long <- function(resp_folds, resp_name) {
  purrr::map_dfr(resp_folds, function(fold_obj) {
    purrr::map_dfr(c("GAM", "RF", "XGB"), function(model_name) {
      if (!is.list(fold_obj) || is.null(fold_obj$best[[model_name]])) return(NULL)
      best_obj <- fold_obj$best[[model_name]]
      if (is.null(best_obj$preds)) return(NULL)
      
      tibble(
        Response = resp_name,
        Fold = fold_obj$fold,
        Model = model_name,
        Predictor = sort(best_obj$preds)
      )
    })
  })
}

best_predictors_long_2024 <- purrr::imap_dfr(best_2024, extract_best_long)

# -----------------------------
# 2) Count predictor-retention frequency
# -----------------------------
predictor_selection_frequency <- best_predictors_long_2024 %>%
  count(Response, Model, Predictor, name = "n_folds") %>%
  mutate(
    Predictor_label = predictor_labels[Predictor]
  )

# -----------------------------
# 3) Summarize frequently retained predictors (>=3/5 folds)
# -----------------------------
selected_predictors_text <- predictor_selection_frequency %>%
  filter(n_folds >= 3) %>%
  arrange(Response, Model, desc(n_folds), Predictor_label) %>%
  mutate(
    pred_txt = paste0(Predictor_label, " (", n_folds, "/5)")
  ) %>%
  group_by(Response, Model) %>%
  summarise(
    Frequently_retained_predictors = paste(pred_txt, collapse = "; "),
    .groups = "drop"
  )

# Fill with "None" when no predictor was retained in at least 3/5 folds
all_resp_model <- expand.grid(
  Response = unique(best_predictors_long_2024$Response),
  Model = c("GAM", "RF", "XGB"),
  stringsAsFactors = FALSE
)

selected_predictors_text <- all_resp_model %>%
  left_join(selected_predictors_text, by = c("Response", "Model")) %>%
  mutate(
    Frequently_retained_predictors = ifelse(is.na(Frequently_retained_predictors), "None", Frequently_retained_predictors)
  )

# -----------------------------
# 4) Performance table with frequently retained predictors
# -----------------------------
summary_2024 <- summary_by_model %>%
  filter(Response %in% names(res_labels)) %>%
  mutate(
    Response_label = factor(unname(res_labels[Response]), levels = resp_order),
    Model = factor(Model, levels = c("GAM", "RF", "XGB"))
  )

winner_2024 <- summary_2024 %>%
  group_by(Response, Response_label) %>%
  arrange(nRMSE_mean, desc(R2_pred_mean), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(Response, Response_label, Best_model = Model)

table_models_2024 <- summary_2024 %>%
  left_join(winner_2024, by = c("Response", "Response_label")) %>%
  left_join(selected_predictors_text, by = c("Response", "Model")) %>%
  mutate(
    Winner = ifelse(Model == Best_model, "Yes", ""),
    `Predictive R2 (mean +/- SD)` = sprintf("%.3f +/- %.3f", R2_pred_mean, R2_pred_sd),
    `nRMSE (%) (mean +/- SD)` = sprintf("%.3f +/- %.3f", nRMSE_mean, nRMSE_sd)
  ) %>%
  select(
    Index = Response_label,
    Model = Model,
    `Frequently retained predictors (>=3/5 outer folds)` = Frequently_retained_predictors,
    `Predictive R2 (mean +/- SD)`,
    `nRMSE (%) (mean +/- SD)`,
    Winner
  ) %>%
  arrange(Index, Model)

print(table_models_2024)

write.csv(
  table_models_2024,
  "model_performance_with_predictor_retention_2024.csv",
  row.names = FALSE
)
table_models_2024


########################### predictor-retention figure

# --------------------------------------------------
# 1) Order of spectral diversity indices and predictors
# --------------------------------------------------
resp_order <- c("S", "RaoQ sc", "RaoQ raw", "H'", "FRic", "FDiv", "FDis")

environmental_predictors <- c(
  "DEM",
  "dist_to_edge_m",
  "dist_to_river_m",
  "TPI_90m",
  "flowdir_90m",
  "aspectHorn_90m",
  "slopeHorn_90m",
  "TRIriley_90m",
  "roughness_90m"
)

predictor_labels <- c(
  DEM             = "Elevation (DEM)",
  dist_to_edge_m  = "Distance to edge",
  dist_to_river_m = "Distance to rivers",
  TPI_90m         = "Topographic position (TPI)",
  flowdir_90m     = "Flow direction",
  aspectHorn_90m  = "Aspect",
  slopeHorn_90m   = "Slope",
  TRIriley_90m    = "Terrain ruggedness (TRI)",
  roughness_90m   = "Roughness"
)

res_labels <- c(
  H_2024        = "H'",
  S_2024        = "S",
  FRic_2024     = "FRic",
  FDiv_2024     = "FDiv",
  FDis_2024     = "FDis",
  RaoRAW_2024 = "RaoQ raw",
  RaoQ_2024  = "RaoQ sc"
)

# --------------------------------------------------
# 2) Heatmap data
# --------------------------------------------------
predictor_selection_frequency <- best_predictors_long_2024 %>%
  count(Response, Model, Predictor, name = "n_folds") %>%
  complete(
    Response = names(res_labels),
    Model = c("GAM", "RF", "XGB"),
    Predictor = environmental_predictors,
    fill = list(n_folds = 0)
  ) %>%
  mutate(
    Index = factor(unname(res_labels[Response]), levels = resp_order),
    Predictor_label = factor(predictor_labels[Predictor], levels = predictor_labels[environmental_predictors]),
    Model = factor(Model, levels = c("GAM", "RF", "XGB")),
    text_col = ifelse(n_folds >= 4, "light", "dark")
  )

# --------------------------------------------------
# 3) Predictor-retention heatmap
# --------------------------------------------------
p_retention_by_model_refined <- ggplot(
  predictor_selection_frequency,
  aes(x = Predictor_label, y = Index, fill = n_folds)
) +
  geom_tile(color = "grey92", linewidth = 0.35) +
  geom_text(
    aes(label = n_folds, color = text_col),
    size = 3.2,
    show.legend = FALSE
  ) +
  facet_wrap(~ Model, nrow = 1) +
  scale_fill_gradient(
    low = "#f7f7f7",
    high = "#2f6db3",
    limits = c(0, 5),
    breaks = 0:5,
    name = "Outer folds"
  ) +
  scale_color_manual(
    values = c("dark" = "black", "light" = "white")
  ) +
  labs(
    x = "Environmental predictors",
    y = "Spectral diversity indices"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

print(p_retention_by_model_refined)

ggsave(
  "Heatmap_predictor_selection_frequency_2024_refined.png",
  p_retention_by_model_refined,
  width = 12.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)


######## observed versus predicted values
############################################################
## TEMPORAL RECURRENCE ASSESSMENT (2020 and 2022)
## Representative configurations selected in 2024 and refitted on the 2024 data

## ========================================================
## 1) Progress reporting
## ========================================================
t0_global <- Sys.time()

total_steps <- 8
current_step <- 0

show_progress <- function(step_name, step_number, total_steps, t0_global) {
  elapsed <- Sys.time() - t0_global
  pct_done <- round(100 * step_number / total_steps, 1)
  pct_left <- round(100 - pct_done, 1)
  
  cat(
    "\n----------------------------------------\n",
    "Step: ", step_number, " of ", total_steps, " - ", step_name, "\n",
    "Completed: ", pct_done, "% | Remaining: ", pct_left, "%\n",
    "Elapsed time: ", round(as.numeric(elapsed, units = "mins"), 2), " min\n",
    "----------------------------------------\n",
    sep = ""
  )
}

## ========================================================
## 2) General definitions
## ========================================================
responses_2024 <- c(
  "H_2024", "S_2024", "FRic_2024", 
  "FDiv_2024", "FDis_2024", "RaoRAW_2024", "RaoQ_2024"
)

res_labels <- c(
  H_2024        = "H'",
  S_2024        = "S",
  FRic_2024     = "FRic",
  FDiv_2024     = "FDiv",
  FDis_2024     = "FDis",
  RaoRAW_2024 = "RaoQ raw",
  RaoQ_2024  = "RaoQ sc"
)

resp_order <- c("S", "RaoQ sc", "RaoQ raw", "H'", "FRic", "FDiv", "FDis")

environmental_predictors <- c(
  "aspectHorn_90m", "slopeHorn_90m", "TPI_90m", "TRIriley_90m",
  "roughness_90m", "flowdir_90m", "dist_to_river_m",
  "dist_to_edge_m", "DEM"
)

# 2024 data
data_2024 <- data_all

best_2024 <- best_pipelines[grepl("_2024$", names(best_pipelines))]


## ========================================================
## 3) Step 1 - Select representative model configurations
## ========================================================
current_step <- 1
cat("\nStarting Step 1: selecting representative model configurations...\n")

extract_candidates <- function(resp_folds, resp_name, model_name) {
  map_dfr(resp_folds, function(fold_obj) {
    if (!is.list(fold_obj) || is.null(fold_obj$best[[model_name]])) return(NULL)
    
    obj <- fold_obj$best[[model_name]]
    if (is.null(obj$preds)) return(NULL)
    
    tibble(
      resp = resp_name,
      fold = fold_obj$fold,
      preds = list(sort(obj$preds)),
      preds_key = paste(sort(obj$preds), collapse = " + "),
      rmse = obj$rmse,
      params = list(obj)
    )
  })
}

# GAM
gam_candidates <- imap_dfr(best_2024, ~ extract_candidates(.x, .y, "GAM")) %>%
  mutate(
    k = map_dbl(params, ~ .x$k),
    gamma = map_dbl(params, ~ .x$gamma)
  )

gam_choice <- gam_candidates %>%
  group_by(resp, preds_key, k, gamma) %>%
  summarise(
    n_folds = n(),
    mean_rmse = mean(rmse, na.rm = TRUE),
    preds = list(preds[[1]]),
    .groups = "drop"
  ) %>%
  group_by(resp) %>%
  arrange(desc(n_folds), mean_rmse, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# RF
rf_candidates <- imap_dfr(best_2024, ~ extract_candidates(.x, .y, "RF")) %>%
  mutate(
    mtry = map_dbl(params, ~ .x$mtry),
    min_node = map_dbl(params, ~ .x$min_node),
    ntree = map_dbl(params, ~ .x$ntree)
  )

rf_choice <- rf_candidates %>%
  group_by(resp, preds_key, mtry, min_node, ntree) %>%
  summarise(
    n_folds = n(),
    mean_rmse = mean(rmse, na.rm = TRUE),
    preds = list(preds[[1]]),
    .groups = "drop"
  ) %>%
  group_by(resp) %>%
  arrange(desc(n_folds), mean_rmse, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# XGB
xgb_candidates <- imap_dfr(best_2024, ~ extract_candidates(.x, .y, "XGB")) %>%
  mutate(
    params_key = map_chr(params, ~ {
      if (is.null(.x$params)) return(NA_character_)
      paste(names(.x$params), .x$params, sep = "=", collapse = "; ")
    })
  )

xgb_choice <- xgb_candidates %>%
  group_by(resp, preds_key, params_key) %>%
  summarise(
    n_folds = n(),
    mean_rmse = mean(rmse, na.rm = TRUE),
    preds = list(preds[[1]]),
    params = list(params[[1]]$params),
    .groups = "drop"
  ) %>%
  group_by(resp) %>%
  arrange(desc(n_folds), mean_rmse, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Confirm that exactly one representative configuration is available for each
# response and algorithm before running the spatial and temporal assessments.
expected_responses <- sort(responses_2024)
if (!identical(sort(gam_choice$resp), expected_responses) ||
    !identical(sort(rf_choice$resp), expected_responses) ||
    !identical(sort(xgb_choice$resp), expected_responses)) {
  stop("Representative model configurations are incomplete or duplicated.")
}

show_progress("Representative model configuration selection", current_step, total_steps, t0_global)

## ========================================================
## 3B) Environmental versus spatial partitioning sensitivity analysis
## ========================================================

run_spatial_sensitivity <- function(df, coord_cols, responses,
                                    gam_choice, rf_choice, xgb_choice,
                                    K = 5, target_blocks = 50,
                                    block_size = NULL,
                                    seed = 1, nthreads = 1){
  results <- list()
  preds_out <- list()
  fold_info <- list()

  for(resp in responses){
    needed <- unique(c(resp, preds_all, coord_cols))
    df0 <- df[, needed, drop = FALSE]
    df0 <- df0[complete.cases(df0), , drop = FALSE]
    df0 <- tibble::as_tibble(df0)

    blocks <- make_spatial_blocks(
      coords = df0[, coord_cols, drop = FALSE],
      K = K,
      target_blocks = target_blocks,
      block_size = block_size,
      seed = seed
    )

    df0$spatial_block <- blocks$block_id
    df0$fold_spatial <- blocks$fold

    fold_info[[resp]] <- tibble::tibble(
      Response = resp,
      x = df0[[coord_cols[1]]],
      y = df0[[coord_cols[2]]],
      spatial_block = df0$spatial_block,
      fold_spatial = df0$fold_spatial,
      block_size = blocks$block_size,
      n_blocks = blocks$n_blocks
    )

    gam_row <- dplyr::filter(gam_choice, .data$resp == .env$resp)
    rf_row  <- dplyr::filter(rf_choice,  .data$resp == .env$resp)
    xgb_row <- dplyr::filter(xgb_choice, .data$resp == .env$resp)

    if(nrow(gam_row) != 1 || nrow(rf_row) != 1 || nrow(xgb_row) != 1){
      stop("Representative model configuration missing or duplicated for ", resp)
    }

    for(f in seq_len(K)){
      cat(
        "\nResponse:", resp,
        "| Fold:", f, "of", K,
        "| Progress:",
        round(
          100 * (
            (match(resp, responses) - 1) * K + f
          ) / (length(responses) * K),
          1
        ),
        "%\n"
      )
      train_df <- dplyr::filter(df0, fold_spatial != f)
      test_df  <- dplyr::filter(df0, fold_spatial == f)

      pr_gam <- fit_pred_gam(
        train_df, test_df, resp, gam_row$preds[[1]],
        k = gam_row$k,
        gamma = gam_row$gamma,
        nthreads = nthreads
      )

      pr_rf <- fit_pred_rf(
        train_df, test_df, resp, rf_row$preds[[1]],
        mtry = rf_row$mtry,
        min_node = rf_row$min_node,
        ntree = rf_row$ntree,
        seed = seed + f,
        nthreads = nthreads
      )

      pr_xgb <- fit_pred_xgb(
        train_df, test_df, resp, xgb_row$preds[[1]],
        params = xgb_row$params[[1]],
        nrounds = XGB_MAX_NROUNDS,
        early = XGB_EARLY_STOPPING,
        validation_fraction = XGB_VALIDATION_FRACTION,
        seed = seed + f
      )

      results[[length(results) + 1]] <- dplyr::bind_rows(
        metrics(test_df[[resp]], pr_gam) %>%
          dplyr::mutate(Response=resp, Model="GAM", Fold=f, Validation="Spatial"),
        metrics(test_df[[resp]], pr_rf) %>%
          dplyr::mutate(Response=resp, Model="RF", Fold=f, Validation="Spatial"),
        metrics(test_df[[resp]], pr_xgb) %>%
          dplyr::mutate(Response=resp, Model="XGB", Fold=f, Validation="Spatial")
      )

      preds_out[[length(preds_out) + 1]] <- dplyr::bind_rows(
        tibble::tibble(
          Response=resp, Model="GAM", Fold=f, Validation="Spatial",
          Observed=test_df[[resp]], Predicted=pr_gam,
          Residual=test_df[[resp]]-pr_gam,
          x=test_df[[coord_cols[1]]], y=test_df[[coord_cols[2]]]
        ),
        tibble::tibble(
          Response=resp, Model="RF", Fold=f, Validation="Spatial",
          Observed=test_df[[resp]], Predicted=pr_rf,
          Residual=test_df[[resp]]-pr_rf,
          x=test_df[[coord_cols[1]]], y=test_df[[coord_cols[2]]]
        ),
        tibble::tibble(
          Response=resp, Model="XGB", Fold=f, Validation="Spatial",
          Observed=test_df[[resp]], Predicted=as.numeric(pr_xgb),
          Residual=test_df[[resp]]-as.numeric(pr_xgb),
          x=test_df[[coord_cols[1]]], y=test_df[[coord_cols[2]]]
        )
      )
    }
  }

  list(
    metrics = dplyr::bind_rows(results),
    predictions = dplyr::bind_rows(preds_out),
    fold_info = dplyr::bind_rows(fold_info)
  )
}

spatial_out <- run_spatial_sensitivity(
  df = data_all,
  coord_cols = COORD_COLS,
  responses = responses_2024,
  gam_choice = gam_choice,
  rf_choice = rf_choice,
  xgb_choice = xgb_choice,
  K = K_OUTER,
  target_blocks = CENTERS,
  block_size = SPATIAL_BLOCK_SIZE,
  seed = SEED,
  nthreads = NTHREADS
)

res_spatial <- spatial_out$metrics
oof_predictions_spatial <- spatial_out$predictions

write.csv(res_spatial, "spatial_blocking_sensitivity_fold_results_2024.csv", row.names = FALSE)
write.csv(oof_predictions_spatial, "spatial_blocking_OOF_predictions_2024.csv", row.names = FALSE)
write.csv(spatial_out$fold_info, "spatial_block_assignments_2024.csv", row.names = FALSE)

validation_results <- dplyr::bind_rows(
  res_nested,
  res_spatial
)

validation_summary <- validation_results %>%
  dplyr::group_by(Response, Model, Validation) %>%
  dplyr::summarise(
    R2_pred_mean = mean(R2_pred, na.rm=TRUE),
    R2_pred_sd = sd(R2_pred, na.rm=TRUE),
    nRMSE_mean = mean(nRMSE, na.rm=TRUE),
    nRMSE_sd = sd(nRMSE, na.rm=TRUE),
    Calibration_intercept_mean = mean(Calibration_intercept, na.rm=TRUE),
    Calibration_intercept_sd = sd(Calibration_intercept, na.rm=TRUE),
    Calibration_slope_mean = mean(Calibration_slope, na.rm=TRUE),
    Calibration_slope_sd = sd(Calibration_slope, na.rm=TRUE),
    .groups = "drop"
  )

write.csv(
  validation_summary,
  "environmental_vs_spatial_validation_summary_2024.csv",
  row.names = FALSE
)

validation_plot_data <- validation_results %>%
  dplyr::mutate(
    Index = factor(unname(res_labels[Response]), levels = unname(res_labels)),
    Model = factor(Model, levels = c("GAM", "RF", "XGB")),
    Validation = factor(Validation, levels = c("Environmental", "Spatial"))
  )

p_validation_r2pred <- ggplot(
  validation_plot_data,
  aes(x = Index, y = R2_pred, fill = Validation)
) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(
    aes(group = Validation),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    size = 1.3, alpha = 0.6
  ) +
  facet_wrap(~ Model, nrow = 1) +
  labs(x = "Spectral diversity index", y = expression(R[pred]^2), fill = "Validation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_validation_nrmse <- ggplot(
  validation_plot_data,
  aes(x = Index, y = nRMSE, fill = Validation)
) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(
    aes(group = Validation),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    size = 1.3, alpha = 0.6
  ) +
  facet_wrap(~ Model, nrow = 1) +
  labs(x = "Spectral diversity index", y = "Normalized RMSE (%)", fill = "Validation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig_validation <- p_validation_r2pred / p_validation_nrmse +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

print(fig_validation)
ggsave(
  "figure_environmental_vs_spatial_validation_R2pred_nRMSE.png",
  fig_validation,
  width = 12, height = 9, dpi = 300
)

## Moran's I permutation test.
moran_permutation <- function(df, value_col, coord_cols = c("x", "y"),
                              k = 8, nsim = 999, max_n = 5000, seed = 1){
  dat <- df[, c(coord_cols, value_col), drop = FALSE]
  names(dat) <- c("x", "y", "value")
  dat <- dat[complete.cases(dat), , drop = FALSE]

  ## Average duplicate coordinates to avoid zero-distance neighbors.
  dat <- dat %>%
    dplyr::group_by(x, y) %>%
    dplyr::summarise(value = mean(value, na.rm=TRUE), .groups = "drop")

  if(nrow(dat) > max_n){
    set.seed(seed)
    dat <- dat[sample(seq_len(nrow(dat)), max_n), , drop = FALSE]
  }

  if(nrow(dat) <= k + 1){
    return(tibble::tibble(
      Moran_I=NA_real_, p_value=NA_real_, n_used=nrow(dat), k=k, nsim=nsim
    ))
  }

  coords <- as.matrix(dat[, c("x", "y")])
  knn <- spdep::knearneigh(coords, k=min(k, nrow(dat)-1), longlat=FALSE)
  nb <- spdep::knn2nb(knn, sym=TRUE)
  lw <- spdep::nb2listw(nb, style="W", zero.policy=TRUE)

  set.seed(seed)
  mor <- spdep::moran.mc(
    dat$value,
    listw=lw,
    nsim=nsim,
    zero.policy=TRUE
  )

  tibble::tibble(
    Moran_I=unname(as.numeric(mor$statistic)),
    p_value=mor$p.value,
    n_used=nrow(dat),
    k=min(k, nrow(dat)-1),
    nsim=nsim
  )
}

## Spatial autocorrelation of environmental predictors.
moran_predictors <- purrr::map_dfr(seq_along(preds_all), function(i){
  pred <- preds_all[i]
  moran_permutation(
    data_all,
    value_col=pred,
    coord_cols=COORD_COLS,
    k=MORAN_K,
    nsim=MORAN_NSIM,
    max_n=MORAN_MAX_N,
    seed=SEED+i
  ) %>%
    dplyr::mutate(Predictor=pred, Predictor_label=var_labels[pred])
})
write.csv(moran_predictors, "Moran_predictors_2024.csv", row.names=FALSE)

## Spatial autocorrelation of combined out-of-fold residuals.
oof_all <- dplyr::bind_rows(oof_predictions_env, oof_predictions_spatial)

moran_residuals <- oof_all %>%
  dplyr::group_by(Response, Model, Validation) %>%
  dplyr::group_modify(~ moran_permutation(
    .x,
    value_col="Residual",
    coord_cols=c("x", "y"),
    k=MORAN_K,
    nsim=MORAN_NSIM,
    max_n=MORAN_MAX_N,
    seed=SEED + sum(utf8ToInt(paste(.y, collapse="")))
  )) %>%
  dplyr::ungroup()

write.csv(
  moran_residuals,
  "Moran_OOF_model_residuals_2024.csv",
  row.names=FALSE
)


## ========================================================
## 4) Step 2 - Functions for model fitting and prediction
## ========================================================
## ========================================================
current_step <- 2
cat("\nStarting Step 2: defining model-fitting and prediction functions...\n")

fit_gam_2024 <- function(resp, choice_row, training_data) {
  vars_top <- choice_row$preds[[1]]
  vars_z <- paste0(vars_top, "_z")

  means <- sapply(training_data[, vars_top, drop = FALSE], mean, na.rm = TRUE)
  sds   <- sapply(training_data[, vars_top, drop = FALSE], sd, na.rm = TRUE)

  scaled_data <- training_data
  for (v in vars_top) {
    scaled_data[[paste0(v, "_z")]] <- (scaled_data[[v]] - means[v]) / sds[v]
  }

  df_fit <- scaled_data %>%
    dplyr::select(dplyr::all_of(c(resp, vars_z))) %>%
    dplyr::filter(complete.cases(.))

  form_gam <- as.formula(
    paste(
      resp, "~",
      paste(sprintf("s(%s, bs='tp', k=%s)", vars_z, choice_row$k), collapse = " + ")
    )
  )

  mod <- mgcv::bam(
    formula = form_gam,
    data = df_fit,
    family = gaussian(),
    method = "fREML",
    gamma = choice_row$gamma,
    select = TRUE,
    discrete = TRUE,
    nthreads = NTHREADS
  )

  list(model = mod, preds = vars_top, means = means, sds = sds)
}

predict_gam_transfer <- function(fit_obj, new_data) {
  vars_top <- fit_obj$preds
  vars_z <- paste0(vars_top, "_z")
  newdat <- new_data

  for (v in vars_top) {
    newdat[[paste0(v, "_z")]] <- (newdat[[v]] - fit_obj$means[v]) / fit_obj$sds[v]
  }

  as.numeric(predict(fit_obj$model, newdata = newdat[, vars_z, drop = FALSE]))
}

fit_rf_2024 <- function(resp, choice_row, training_data) {
  vars_top <- choice_row$preds[[1]]

  df_fit <- training_data %>%
    dplyr::select(dplyr::all_of(c(resp, vars_top))) %>%
    dplyr::filter(complete.cases(.))

  mod <- ranger::ranger(
    formula = as.formula(paste(resp, "~ .")),
    data = df_fit,
    num.trees = choice_row$ntree,
    mtry = choice_row$mtry,
    min.node.size = choice_row$min_node,
    importance = "permutation",
    seed = SEED,
    num.threads = NTHREADS
  )

  list(model = mod, preds = vars_top)
}

predict_rf_transfer <- function(fit_obj, new_data) {
  as.numeric(
    predict(
      fit_obj$model,
      data = new_data[, fit_obj$preds, drop = FALSE]
    )$predictions
  )
}

fit_xgb_2024 <- function(resp, choice_row, training_data) {
  vars_top <- choice_row$preds[[1]]
  params <- choice_row$params[[1]]

  df_fit <- training_data %>%
    dplyr::select(dplyr::all_of(c(resp, vars_top))) %>%
    dplyr::filter(complete.cases(.))

  set.seed(SEED)
  n_val <- max(10, floor(nrow(df_fit) * XGB_VALIDATION_FRACTION))
  n_val <- min(n_val, nrow(df_fit) - 10)
  idx_val <- sample(seq_len(nrow(df_fit)), n_val)

  train_fit <- df_fit[-idx_val, , drop = FALSE]
  train_val <- df_fit[idx_val, , drop = FALSE]

  dtr <- xgb.DMatrix(
    as.matrix(train_fit[, vars_top, drop = FALSE]),
    label = train_fit[[resp]]
  )
  dval <- xgb.DMatrix(
    as.matrix(train_val[, vars_top, drop = FALSE]),
    label = train_val[[resp]]
  )

  mod <- xgb.train(
    params = params,
    data = dtr,
    nrounds = XGB_MAX_NROUNDS,
    watchlist = list(train = dtr, validation = dval),
    early_stopping_rounds = XGB_EARLY_STOPPING,
    maximize = FALSE,
    verbose = 0
  )

  list(model = mod, preds = vars_top, best_iteration = mod$best_iteration)
}

predict_xgb_transfer <- function(fit_obj, new_data) {
  X <- xgb.DMatrix(as.matrix(new_data[, fit_obj$preds, drop = FALSE]))
  predict_xgb_best(fit_obj$model, X)
}

show_progress("Model-fitting and prediction functions", current_step, total_steps, t0_global)

## ========================================================
## 5) Step 3 - Refit representative configurations on the full 2024 data
## ========================================================
current_step <- 3
cat("\nStarting Step 3: refitting representative configurations on the 2024 data...\n")

gam_models_2024 <- setNames(vector("list", nrow(gam_choice)), gam_choice$resp)
for (i in seq_len(nrow(gam_choice))) {
  rr <- gam_choice$resp[i]
  cat("  GAM:", rr, "(", i, "/", nrow(gam_choice), ")\n")
  gam_models_2024[[rr]] <- fit_gam_2024(rr, gam_choice[i, ], data_2024)
}

rf_models_2024 <- setNames(vector("list", nrow(rf_choice)), rf_choice$resp)
for (i in seq_len(nrow(rf_choice))) {
  rr <- rf_choice$resp[i]
  cat("  RF:", rr, "(", i, "/", nrow(rf_choice), ")\n")
  rf_models_2024[[rr]] <- fit_rf_2024(rr, rf_choice[i, ], data_2024)
}

xgb_models_2024 <- setNames(vector("list", nrow(xgb_choice)), xgb_choice$resp)
for (i in seq_len(nrow(xgb_choice))) {
  rr <- xgb_choice$resp[i]
  cat("  XGB:", rr, "(", i, "/", nrow(xgb_choice), ")\n")
  xgb_models_2024[[rr]] <- fit_xgb_2024(rr, xgb_choice[i, ], data_2024)
}

show_progress("Refitting representative configurations on 2024 data", current_step, total_steps, t0_global)

## ========================================================
## 6) Step 4 - Performance metrics
## ========================================================
current_step <- 4
cat("\nStarting Step 4: defining performance metrics...\n")

calc_metrics <- function(obs, pred) {
  metrics(obs, pred)
}

show_progress("Performance metrics", current_step, total_steps, t0_global)

## ========================================================
## 7) Step 5 - Temporal predictions for 2020 and 2022
## ========================================================
current_step <- 5
cat("\nStarting Step 5: generating temporal predictions for 2020 and 2022...\n")

years_test <- c(2020, 2022)

pred_results <- list()
metrics_results <- list()

counter <- 0
total_jobs <- length(responses_2024) * length(years_test) * 3

for (resp_2024 in responses_2024) {
  for (yr in years_test) {
    resp_test <- sub("_2024$", paste0("_", yr), resp_2024)
    if (!resp_test %in% names(data_all)) next
    
    obs <- data_all[[resp_test]]
    
    counter <- counter + 1
    cat("  Job", counter, "/", total_jobs, "-", resp_2024, yr, "GAM\n")
    pred_gam <- predict_gam_transfer(gam_models_2024[[resp_2024]], data_all)
    metrics_results[[length(metrics_results) + 1]] <- calc_metrics(obs, pred_gam) %>%
      mutate(Response = resp_2024, Year = yr, Model = "GAM")
    pred_results[[length(pred_results) + 1]] <- tibble(
      Response = resp_2024, Year = yr, Model = "GAM",
      Observed = obs, Predicted = pred_gam
    )
    
    counter <- counter + 1
    cat("  Job", counter, "/", total_jobs, "-", resp_2024, yr, "RF\n")
    pred_rf <- predict_rf_transfer(rf_models_2024[[resp_2024]], data_all)
    metrics_results[[length(metrics_results) + 1]] <- calc_metrics(obs, pred_rf) %>%
      mutate(Response = resp_2024, Year = yr, Model = "RF")
    pred_results[[length(pred_results) + 1]] <- tibble(
      Response = resp_2024, Year = yr, Model = "RF",
      Observed = obs, Predicted = pred_rf
    )
    
    counter <- counter + 1
    cat("  Job", counter, "/", total_jobs, "-", resp_2024, yr, "XGB\n")
    pred_xgb <- predict_xgb_transfer(xgb_models_2024[[resp_2024]], data_all)
    metrics_results[[length(metrics_results) + 1]] <- calc_metrics(obs, pred_xgb) %>%
      mutate(Response = resp_2024, Year = yr, Model = "XGB")
    pred_results[[length(pred_results) + 1]] <- tibble(
      Response = resp_2024, Year = yr, Model = "XGB",
      Observed = obs, Predicted = pred_xgb
    )
  }
}

pred_results <- bind_rows(pred_results) %>%
  mutate(Index = res_labels[Response])

metrics_results <- bind_rows(metrics_results) %>%
  mutate(Index = res_labels[Response])

show_progress("Temporal predictions for 2020 and 2022", current_step, total_steps, t0_global)

## ========================================================
## 8) Step 6 - Save temporal performance tables
## ========================================================
current_step <- 6
cat("\nStarting Step 6: saving temporal performance tables...\n")

write.csv(metrics_results, "Temporal_recurrence_metrics_2020_2022.csv", row.names = FALSE)
write.csv(pred_results, "Temporal_recurrence_predictions_2020_2022.csv", row.names = FALSE)

metrics_table <- metrics_results %>%
  mutate(
    `Predictive R2` = sprintf("%.3f", R2_pred),
    nRMSE = sprintf("%.3f", nRMSE),
    `Calibration intercept` = sprintf("%.3f", Calibration_intercept),
    `Calibration slope` = sprintf("%.3f", Calibration_slope)
  ) %>%
  select(
    Index, Year, Model,
    `Predictive R2`, nRMSE,
    `Calibration intercept`, `Calibration slope`, n
  )

write.csv(metrics_table, "Temporal_recurrence_metrics_table_2020_2022.csv", row.names = FALSE)

show_progress("Saving temporal performance tables", current_step, total_steps, t0_global)

## ========================================================
## 9) Step 7 - Observed-versus-predicted plotting function
## ========================================================
current_step <- 7
cat("\nStarting Step 7: defining the observed-versus-predicted plotting function...\n")

plot_obs_pred_index <- function(index_name) {
  ggplot(
    pred_results %>% filter(Index == index_name),
    aes(x = Observed, y = Predicted)
  ) +
    geom_point(alpha = 0.35, size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    facet_grid(Year ~ Model, scales = "free") +
    labs(
      x = paste("Observed", index_name),
      y = paste("Predicted", index_name)
    ) +
    theme_bw(base_size = 11)
}

show_progress("Observed-versus-predicted plotting function", current_step, total_steps, t0_global)

## ========================================================
## 10) Step 8 - Generate temporal figures and heatmaps
## ========================================================
current_step <- 8
cat("\nStarting Step 8: generating temporal figures and heatmaps...\n")

# 8A) Observed-versus-predicted plots
idxs <- unique(pred_results$Index)

for (i in seq_along(idxs)) {
  idx <- idxs[i]
  cat("  Observed versus predicted plot", i, "/", length(idxs), "-", idx, "\n")
  p <- plot_obs_pred_index(idx) +
    theme_minimal(base_size = 16) +  # overall font size
    theme(
      plot.title = element_text(size = 19, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 16, face = "bold"),  # facet-strip labels
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  ggsave(
    filename = paste0("Obs_vs_Pred_", gsub("[^A-Za-z0-9]", "_", idx), "_2020_2022.png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# 8B) Predictive R2 heatmap
heat_r2pred <- metrics_results %>%
  mutate(
    Index = factor(Index, levels = resp_order),
    Model = factor(Model, levels = c("GAM", "RF", "XGB")),
    Year = factor(Year)
  )

p_heat_r2pred <- ggplot(heat_r2pred, aes(x = Model, y = Index, fill = R2_pred)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", R2_pred)), size = 3.2) +
  facet_wrap(~ Year, nrow = 1) +
  scale_fill_gradient(
    low = "#f7f7f7",
    high = "#2f6db3",
    name = expression(R[pred]^2)
  ) +
  labs(
    x = "Model",
    y = "Spectral diversity indices"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold")
  )

print(p_heat_r2pred)

ggsave(
  "Heatmap_R2pred_temporal_recurrence_2020_2022.png",
  p_heat_r2pred,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

# 8C) nRMSE heatmap
p_heat_nrmse <- ggplot(heat_r2pred, aes(x = Model, y = Index, fill = nRMSE)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", nRMSE)), size = 3.2) +
  facet_wrap(~ Year, nrow = 1) +
  scale_fill_gradient(
    low = "#2f6db3",
    high = "#f7f7f7",
    name = "nRMSE"
  ) +
  labs(
    x = "Model",
    y = "Spectral diversity indices"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold")
  )

print(p_heat_nrmse)

ggsave(
  "Heatmap_nRMSE_temporal_recurrence_2020_2022.png",
  p_heat_nrmse,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

show_progress("Temporal figures and heatmaps", current_step, total_steps, t0_global)

cat("\nANALYSIS COMPLETED.\n")


## ========================================================
## 11) Simplified workflow diagram
## ========================================================
workflow <- tibble::tibble(
  step = 1:8,
  label = c(
    "Data preparation\nand collinearity",
    "Environmental clusters\n(k-means, k = 50)",
    "Outer environmental\n5-fold CV",
    "Inner CV:\nsubset screening + tuning",
    "Out-of-fold predictions\nand performance",
    "Spatial-block\nsensitivity analysis",
    "Moran's I of\nOOF residuals",
    "Temporal recurrence:\n2020 and 2022"
  ),
  x = 1:8,
  y = 1
)

p_workflow <- ggplot(workflow, aes(x = x, y = y)) +
  geom_segment(
    data = workflow %>% dplyr::filter(step < max(step)),
    aes(x = x + 0.32, xend = x + 0.68, y = y, yend = y),
    arrow = grid::arrow(length = grid::unit(0.11, "inches")),
    inherit.aes = FALSE
  ) +
  geom_label(aes(label = label), size = 3.0, label.size = 0.25) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.04))) +
  coord_cartesian(clip = "off") +
  theme_void()

print(p_workflow)
ggsave("workflow_modeling.png", p_workflow, width = 16, height = 3.2, dpi = 300)

## ========================================================
## 12) Computational reproducibility
## ========================================================
analysis_config <- list(
  data_file = DATA_FILE,
  responses = responses,
  predictors = preds_all,
  coordinate_columns = COORD_COLS,
  environmental_clusters = CENTERS,
  outer_folds = K_OUTER,
  inner_folds = K_INNER,
  minimum_subset_size = MIN_SUBSET,
  stage1_max_n = STAGE1_N,
  stage1_folds = K_STAGE1,
  top_subsets_refined = TOPN,
  correlation_threshold = COR_MAX,
  vif_threshold = VIF_MAX,
  gam_basis = "thin plate regression spline, bs='tp'",
  gam_k_grid = c(5, 6, 8),
  gam_gamma_grid = c(1.0, 1.4),
  rf_min_node_grid = c(5, 10),
  rf_num_trees = 500,
  xgb_max_depth_grid = c(3, 6),
  xgb_eta_grid = c(0.05, 0.10),
  xgb_subsample_grid = c(0.7, 1.0),
  xgb_colsample_grid = c(0.7, 1.0),
  xgb_max_nrounds = XGB_MAX_NROUNDS,
  xgb_early_stopping_rounds = XGB_EARLY_STOPPING,
  xgb_internal_validation_fraction = XGB_VALIDATION_FRACTION,
  spatial_block_size = unique(spatial_out$fold_info$block_size),
  number_spatial_blocks = unique(spatial_out$fold_info$n_blocks),
  moran_k = MORAN_K,
  moran_permutations = MORAN_NSIM,
  moran_max_sample = MORAN_MAX_N,
  random_seed = SEED,
  threads = NTHREADS
)

saveRDS(analysis_config, "analysis_configuration.rds")
capture.output(sessionInfo(), file = "sessionInfo_modeling.txt")

package_names <- c(
  "dplyr", "purrr", "tidyr", "tibble", "mgcv",
  "ranger", "xgboost", "ggplot2", "patchwork", "spdep"
)
package_versions <- tibble::tibble(
  package = package_names,
  version = vapply(package_names, function(p) as.character(packageVersion(p)), character(1))
)
write.csv(package_versions, "package_versions_modeling.csv", row.names = FALSE)

cat("\nSCRIPT COMPLETED.\n")

