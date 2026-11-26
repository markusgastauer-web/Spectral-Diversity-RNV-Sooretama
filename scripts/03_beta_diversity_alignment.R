library(terra)
library(vegan)

# -------------------------------------------------------------------------
# Input directories
# Replace all "path/to/..." placeholders with the locations of your files.
#
# Required directories:
#   - beta_diversity_2020/ : beta-diversity iteration rasters for 2020
#   - beta_diversity_2022/ : beta-diversity iteration rasters for 2022
#   - beta_diversity_2024/ : beta-diversity iteration rasters for 2024
#
# Each directory should contain files matching:
#   BetaDiversity_BCdiss_PCO*
#
# Outputs:
#   - BetaDiversity_mean.tif
#   - BetaDiversity_sd.tif
# -------------------------------------------------------------------------

# --------------------------------
# Read all beta diversity rasters
# --------------------------------

path <- "path/to/beta_diversity_iterations"

files <- list.files(
  path,
  pattern = "BetaDiversity_BCdiss_PCO.*",
  full.names = TRUE
)

files <- files[!grepl("\\.hdr$", files)]
files <- files[!grepl("\\.aux.xml$", files)]

r_list <- lapply(files, rast)

# --------------------------------
# Reference raster (iteration 1)
# --------------------------------

ref_r <- r_list[[1]]

ref_vals <- values(ref_r)

valid_ref <- complete.cases(ref_vals)

ref_coords <- ref_vals[valid_ref, ]

# --------------------------------
# Align all iterations
# --------------------------------

aligned_rasters <- list()

aligned_rasters[[1]] <- ref_r

for(i in 2:length(r_list)) {
  
  cat("Aligning:", i, "\n")
  
  r <- r_list[[i]]
  
  vals <- values(r)
  
  valid <- complete.cases(vals)
  
  coords <- vals[valid, ]
  
  # Must use identical pixels
  common <- valid & valid_ref
  
  X <- ref_vals[common, ]
  
  Y <- vals[common, ]
  
  proc <- procrustes(
    X,
    Y,
    scale = TRUE
  )
  
  aligned <- predict(
    proc,
    newdata = coords
  )
  
  out <- matrix(
    NA,
    nrow = ncell(r),
    ncol = 3
  )
  
  out[valid, ] <- aligned
  
  r_aligned <- r
  
  values(r_aligned) <- out
  
  aligned_rasters[[i]] <- r_aligned
}

# --------------------------------
# Mean and SD
# --------------------------------

stacked <- rast(aligned_rasters)

groups <- rep(1:3, length(r_list))

mean_beta <- tapp(
  stacked,
  groups,
  mean,
  na.rm = TRUE
)

sd_beta <- tapp(
  stacked,
  groups,
  sd,
  na.rm = TRUE
)

names(mean_beta) <- c(
  "PCo1",
  "PCo2",
  "PCo3"
)

names(sd_beta) <- c(
  "PCo1_sd",
  "PCo2_sd",
  "PCo3_sd"
)

writeRaster(
  mean_beta,
  file.path(path, "BetaDiversity_mean.tif"),
  overwrite = TRUE
)

writeRaster(
  sd_beta,
  file.path(path, "BetaDiversity_sd.tif"),
  overwrite = TRUE
)

# --------------------------------
# Procrustes test among iterations
# --------------------------------

files <- list.files(
  path,
  pattern = "BetaDiversity_BCdiss_PCO.*",
  full.names = TRUE
)

files <- files[!grepl("\\.hdr$", files)]
files <- files[!grepl("\\.aux.xml$", files)]

r_list <- lapply(files, rast)

# Reference = first iteration

ref <- values(r_list[[1]])

valid_ref <- complete.cases(ref)

results <- data.frame(
  comparison = character(),
  correlation = numeric(),
  pvalue = numeric()
)

for(i in 2:length(r_list)) {
  
  vals <- values(r_list[[i]])
  
  common <- valid_ref &
    complete.cases(vals)
  
  X <- ref[common, ]
  
  Y <- vals[common, ]
  
  # Procrustes test
  ptest <- protest(
    X,
    Y,
    permutations = 999
  )
  
  results <- rbind(
    results,
    data.frame(
      comparison = paste0(
        "Iter001_vs_Iter",
        sprintf("%03d", i)
      ),
      correlation = ptest$t0,
      pvalue = ptest$signif
    )
  )
}

write.csv(
  results,
  file.path(
    path,
    "within_year_procrustes_statistics.csv"
  ),
  row.names = FALSE
)

# --------------------------------
# Comparison among years
# --------------------------------

years <- c("2020","2022","2024")

path_2020 <- "path/to/beta_diversity_2020"

path_2022 <- "path/to/beta_diversity_2022"

path_2024 <- "path/to/beta_diversity_2024"

# --------------------------------
# Temporal alignment of annual consensus ordinations
# --------------------------------

path_2020 <- "path/to/beta_diversity_2020"
path_2022 <- "path/to/beta_diversity_2022"
path_2024 <- "path/to/beta_diversity_2024"

beta_2020 <- rast(file.path(path_2020, "BetaDiversity_mean.tif"))
beta_2022 <- rast(file.path(path_2022, "BetaDiversity_mean.tif"))
beta_2024 <- rast(file.path(path_2024, "BetaDiversity_mean.tif"))

# Reference year = 2024

ref_vals <- values(beta_2024)

valid_ref <- complete.cases(ref_vals)

align_to_reference <- function(target_raster,
                               reference_vals,
                               reference_mask){
  
  target_vals <- values(target_raster)
  
  common <- reference_mask &
    complete.cases(target_vals)
  
  X <- reference_vals[common, ]
  Y <- target_vals[common, ]
  
  proc <- procrustes(
    X,
    Y,
    scale = TRUE
  )
  
  aligned <- predict(
    proc,
    newdata = Y
  )
  
  out <- matrix(
    NA,
    nrow = nrow(target_vals),
    ncol = ncol(target_vals)
  )
  
  out[common, ] <- aligned
  
  aligned_raster <- target_raster
  values(aligned_raster) <- out
  
  list(
    raster = aligned_raster,
    proc = proc
  )
}

# Align annual consensus ordinations to 2024

aligned_2020 <- align_to_reference(
  beta_2020,
  ref_vals,
  valid_ref
)

aligned_2022 <- align_to_reference(
  beta_2022,
  ref_vals,
  valid_ref
)

# 2024 remains unchanged because it is the reference

aligned_2024 <- beta_2024

# Export aligned rasters

writeRaster(
  aligned_2020$raster,
  file.path(
    path_2020,
    "BetaDiversity_2020_mean_aligned_to_2024.tif"
  ),
  overwrite = TRUE
)

writeRaster(
  aligned_2022$raster,
  file.path(
    path_2022,
    "BetaDiversity_2022_mean_aligned_to_2024.tif"
  ),
  overwrite = TRUE
)

writeRaster(
  aligned_2024,
  file.path(
    path_2024,
    "BetaDiversity_2024_mean_aligned_to_2024.tif"
  ),
  overwrite = TRUE
)

files <- list(
  file.path(path_2020,
            "BetaDiversity_2020_mean_aligned_to_2024.tif"),
  file.path(path_2022,
            "BetaDiversity_2022_mean_aligned_to_2024.tif"),
  file.path(path_2024,
            "BetaDiversity_2024_mean_aligned_to_2024.tif")
)

r_list <- lapply(files, rast)

# Reference year = 2024 consensus ordination
ref <- values(r_list[[3]])

valid_ref <- complete.cases(ref)

results <- data.frame(
  comparison = character(),
  correlation = numeric(),
  pvalue = numeric()
)

for(i in c(1,2)) {
  
  vals <- values(r_list[[i]])
  
  common <- valid_ref &
    complete.cases(vals)
  
  X <- ref[common, ]
  
  Y <- vals[common, ]
  
  # Procrustes test
  ptest <- protest(
    X,
    Y,
    permutations = 999
  )
  
  results <- rbind(
    results,
    data.frame(
      comparison = paste0(
        years[i],
        "_vs_2024"
      ),
      correlation = ptest$t0,
      pvalue = ptest$signif
    )
  )
}

results

write.csv(
  results,
  "between_year_procrustes_statistics.csv",
  row.names = FALSE
)
