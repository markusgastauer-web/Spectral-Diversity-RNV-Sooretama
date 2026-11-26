# -------------------------------------------------------------------------
# Libraries
# -------------------------------------------------------------------------

library(terra)
library(spatialEco)
library(sf)
library(corrplot)
library(lidR)

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

input_dir  <- file.path("data", "raw")
output_dir <- file.path("outputs", "predictors")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dem_file       <- file.path(input_dir, "DEM.tif")
diversity_file <- file.path(input_dir, "diversity_stack.tif")
streams_file   <- file.path(input_dir, "streams.gpkg")
landcover_file <- file.path(input_dir, "Coverage_RNV.tif")

# -------------------------------------------------------------------------
# Read input data
# -------------------------------------------------------------------------

dem_original <- rast(dem_file)
diversity    <- rast(diversity_file)
rivers       <- vect(streams_file)
landcover    <- rast(landcover_file)

# Confirm layer names
names(diversity)
names(landcover)

if (is.na(crs(diversity))) {
  crs(diversity) <- "EPSG:32724"
}

# -------------------------------------------------------------------------
# Check and harmonize coordinate systems
# -------------------------------------------------------------------------

if (is.na(crs(dem_original)) || crs(dem_original) == "") {
  stop("The DEM has no defined CRS.")
}

if (is.na(crs(diversity)) || crs(diversity) == "") {
  stop("The diversity raster has no defined CRS.")
}

if (is.na(crs(rivers)) || crs(rivers) == "") {
  stop("The river layer has no defined CRS.")
}

# Reproject instead of assigning CRS
if (!same.crs(rivers, diversity)) {
  rivers <- project(rivers, crs(diversity))
}

if (!same.crs(landcover, diversity)) {
  landcover <- project(
    landcover,
    diversity,
    method = "near"
  )
}

# -------------------------------------------------------------------------
# Terrain predictors at diversity-raster resolution
# -------------------------------------------------------------------------

dem <- project(
  dem_original,
  diversity,
  method = "bilinear"
)

names(dem) <- "DEM"

slope <- terrain(
  dem,
  v = "slope",
  unit = "degrees",
  neighbors = 8
)
names(slope) <- "slope"

aspect <- terrain(
  dem,
  v = "aspect",
  unit = "degrees",
  neighbors = 8
)
names(aspect) <- "aspect"

tri <- terrain(
  dem,
  v = "TRI",
  neighbors = 8
)
names(tri) <- "TRI"

tpi <- terrain(
  dem,
  v = "TPI",
  neighbors = 8
)
names(tpi) <- "TPI"

writeRaster(
  dem,
  file.path(output_dir, "resampled_dem.tif"),
  overwrite = TRUE
)

writeRaster(
  slope,
  file.path(output_dir, "slope.tif"),
  overwrite = TRUE
)

writeRaster(
  aspect,
  file.path(output_dir, "aspect.tif"),
  overwrite = TRUE
)

writeRaster(
  tri,
  file.path(output_dir, "tri.tif"),
  overwrite = TRUE
)

writeRaster(
  tpi,
  file.path(output_dir, "tpi.tif"),
  overwrite = TRUE
)

# -------------------------------------------------------------------------
# DEM at approximately 30 m for within-cell terrain heterogeneity
# -------------------------------------------------------------------------

dem_30m_template <- rast(
  ext(diversity),
  resolution = 90 / 3,
  crs = crs(diversity)
)

dem_30m <- project(
  dem_original,
  dem_30m_template,
  method = "bilinear"
)

names(dem_30m) <- "DEM_30m"

writeRaster(
  dem_30m,
  file.path(output_dir, "resampled_dem_30m.tif"),
  overwrite = TRUE
)

# -------------------------------------------------------------------------
# Vector ruggedness measure
# -------------------------------------------------------------------------

vrm_in <- spatialEco::vrm(dem_30m, s = 3)
names(vrm_in) <- "VRM_IN"

vrm_out <- spatialEco::vrm(dem_30m, s = 9)
names(vrm_out) <- "VRM_OUT"

writeRaster(
  vrm_in,
  file.path(output_dir, "vrm_in.tif"),
  overwrite = TRUE
)

writeRaster(
  vrm_out,
  file.path(output_dir, "vrm_out.tif"),
  overwrite = TRUE
)

# Align VRM rasters with the diversity raster for extraction
vrm_in_div <- resample(vrm_in, diversity[[1]], method = "bilinear")
vrm_out_div <- resample(vrm_out, diversity[[1]], method = "bilinear")

# -------------------------------------------------------------------------
# Create cell ID raster while retaining the original NA mask
# -------------------------------------------------------------------------

if (!"H" %in% names(diversity)) {
  stop("Layer 'H' is not present in the diversity raster.")
}

valid_cells <- !is.na(diversity[["H"]])

id_values <- seq_len(ncell(diversity[["H"]]))
id_values[!as.vector(values(valid_cells))] <- NA

id_raster <- diversity[["H"]]
values(id_raster) <- id_values
names(id_raster) <- "plot"

grid <- as.polygons(
  id_raster,
  dissolve = FALSE,
  na.rm = TRUE
)

names(grid) <- "plot"

writeVector(
  grid,
  file.path(output_dir, "grid_polygons.gpkg"),
  layer = "grid",
  overwrite = TRUE
)

grid_centroids <- centroids(grid)

# -------------------------------------------------------------------------
# Rumple and roughness for each diversity-raster cell
# -------------------------------------------------------------------------

terrain_results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  
  dem_cell <- crop(
    dem_30m,
    grid[i, ],
    mask = TRUE
  )
  
  # Coordinates and elevation of valid cells
  xyz <- as.data.frame(
    dem_cell,
    xy = TRUE,
    na.rm = TRUE
  )
  
  # Standardize elevation column name
  names(xyz)[3] <- "z"
  
  if (nrow(xyz) < 3) {
    
    terrain_results[[i]] <- data.frame(
      plot      = grid$plot[i],
      rumple    = NA_real_,
      rough_abs = NA_real_,
      rough_sd  = NA_real_
    )
    
  } else {
    
    rumple_value <- lidR::rumple_index(
      x = xyz$x,
      y = xyz$y,
      z = xyz$z
    )
    
    terrain_results[[i]] <- data.frame(
      plot      = grid$plot[i],
      rumple    = as.numeric(rumple_value),
      rough_abs = diff(range(xyz$z, na.rm = TRUE)),
      rough_sd  = sd(xyz$z, na.rm = TRUE)
    )
  }
}

terrain_df <- do.call(rbind, terrain_results)

write.table(
  terrain_df,
  file.path(output_dir, "RNV_rumple_roughness.txt"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# -------------------------------------------------------------------------
# Distance to rivers
# -------------------------------------------------------------------------

river_distance <- distance(grid_centroids, rivers)

if (inherits(river_distance, "dist")) {
  river_distance <- as.matrix(river_distance)
}

if (is.matrix(river_distance)) {
  river_distance <- apply(river_distance, 1, min, na.rm = TRUE)
}

river_df <- data.frame(
  plot = grid$plot,
  river_dist = as.numeric(river_distance)
)

# -------------------------------------------------------------------------
# Distance to forest edge
# Assumption: land-cover class 3 represents forest
# -------------------------------------------------------------------------

# MapBiomas Collection 9:
# class 3 = Forest Formation
forest <- landcover == 3
forest[forest == 0] <- NA

forest_polygons <- as.polygons(
  forest,
  dissolve = TRUE,
  na.rm = TRUE
)

# Dissolve forest polygons before creating boundaries
forest_polygons <- aggregate(forest_polygons)

forest_edges <- as.lines(forest_polygons)

edge_distance <- distance(grid_centroids, forest_edges)

if (inherits(edge_distance, "dist")) {
  edge_distance <- as.matrix(edge_distance)
}

if (is.matrix(edge_distance)) {
  edge_distance <- apply(edge_distance, 1, min, na.rm = TRUE)
}

edge_df <- data.frame(
  plot = grid$plot,
  edge_dist = as.numeric(edge_distance)
)

# -------------------------------------------------------------------------
# Extract raster predictors
# -------------------------------------------------------------------------

predictor_stack <- c(
  dem,
  slope,
  aspect,
  tri,
  tpi,
  vrm_in_div,
  vrm_out_div
)

predictor_values <- extract(
  predictor_stack,
  grid_centroids,
  ID = FALSE
)

predictor_df <- data.frame(
  plot = grid$plot,
  predictor_values
)

# -------------------------------------------------------------------------
# Extract diversity metrics
# -------------------------------------------------------------------------

required_metrics <- c(
  "Shade",
  "S",
  "H",
  "FRic",
  "FDiv",
  "FDis",
  "RaoQ"
)

missing_metrics <- setdiff(required_metrics, names(diversity))

if (length(missing_metrics) > 0) {
  stop(
    paste(
      "Missing diversity layers:",
      paste(missing_metrics, collapse = ", ")
    )
  )
}

diversity_values <- extract(
  diversity[[required_metrics]],
  grid_centroids,
  ID = FALSE
)

diversity_df <- data.frame(
  plot = grid$plot,
  diversity_values
)

# -------------------------------------------------------------------------
# Merge all information by plot ID
# -------------------------------------------------------------------------

merged_df <- Reduce(
  function(x, y) merge(x, y, by = "plot", all = TRUE),
  list(
    river_df,
    edge_df,
    predictor_df,
    terrain_df,
    diversity_df
  )
)

# Check whether merging generated duplicate rows
if (anyDuplicated(merged_df$plot)) {
  warning("Duplicate plot IDs were generated during merging.")
}

# -------------------------------------------------------------------------
# Remove cells based on shade threshold
# -------------------------------------------------------------------------

merged_df_cl <- subset(
  merged_df,
  !is.na(Shade) & Shade >= 77
)

write.table(
  merged_df_cl,
  file.path(
    output_dir,
    "df_RNV_July2026_clean.txt"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

