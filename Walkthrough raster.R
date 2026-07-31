# On the sensitivities to the modifiable areal unit problem
# Raster walk-through: PRISM precipitation and county zones
# YE, Xiang 叶翔; CHEN, Jiayi 陈佳怡
# 2026-07-01

# 13. Setup for raster data processing
library(terra)
library(sf)

# Set project directory.
# Put all raster-related input data under this folder.
project_dir <- "D:/Change_to_your_local_path"

# Source raster sensitivity functions.
source(file.path(project_dir, "MAUP sensitivity raster.R"))

# Parallel Monte Carlo settings.
# Set use_parallel <- FALSE if you want serial execution.
use_parallel <- TRUE
parallel_workers <- default_parallel_workers_raster(reserve_cores = 5L)

# Required input folders and files:
#
# (a) PRISM annual precipitation raster
#    Put this file here:
#    D:/MAUP_sen/data/raster_ppt_2023/prism_ppt_us_30s_2023.tif
#
# (b) 2023 U.S. county boundary shapefile
#    Put the full shapefile set here:
#    D:/MAUP_sen/data/US_county/
#
#    The folder should contain at least:
#    tl_2023_us_county.shp
#    tl_2023_us_county.dbf
#    tl_2023_us_county.shx
#    tl_2023_us_county.prj
#
# If your data are stored elsewhere, change the two paths below.
tif_file_path <- file.path(
  project_dir,
  "data",
  "raster_ppt_2023",
  "prism_ppt_us_30s_2023.tif"
)

county_shp <- file.path(
  project_dir,
  "data",
  "US_county",
  "tl_2023_us_county.shp"
)

# Create output folders.
output_dir <- file.path(project_dir, "outputs")
raster_dir <- file.path(output_dir, "raster")
table_dir <- file.path(output_dir, "tables")
rds_dir <- file.path(output_dir, "rds")
fig_dir <- file.path(output_dir, "figures")

dir.create(file.path(project_dir, "outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# 14. Check input files

if (!file.exists(tif_file_path)) {
  stop(
    "PRISM raster file was not found. Please put prism_ppt_us_30s_2023.tif in: ",
    file.path(project_dir, "data", "raster_ppt_2023")
  )
}

if (!file.exists(county_shp)) {
  stop(
    "County shapefile was not found. Please put tl_2023_us_county.* in: ",
    file.path(project_dir, "data", "US_county")
  )
}


# 15. Load raster and county boundary data

# PRISM 2023 annual precipitation raster for the conterminous United States.
# The original resolution is approximately 800 m.
r_raw <- rast(tif_file_path)

# Read the U.S. county boundary shapefile.
us_counties_sf <- st_read(county_shp, quiet = TRUE)

# Filter New York counties only. State FIPS code for New York is 36.
ny_sf <- us_counties_sf[us_counties_sf$STATEFP == "36", ]

# Convert the sf object to a terra SpatVector for raster operations.
ny_vect <- vect(ny_sf)


# 16. Align CRS, crop raster, and build county-zone raster

# Reproject the county vector to the CRS of the PRISM raster.
ny_vect <- project(ny_vect, crs(r_raw))

# Convert county GEOID to numeric values so it can be stored in a raster layer.
ny_vect$GEOID_NUM <- as.numeric(ny_vect$GEOID)

# Crop and mask the PRISM raster to New York.
r_ny <- crop(r_raw, ny_vect)
r_ny <- mask(r_ny, ny_vect)

# Rasterize county membership on the same grid as the precipitation raster.
# Each raster cell stores the GEOID of the county it belongs to.
r_zones <- rasterize(ny_vect, r_ny, field = "GEOID_NUM")

# Check whether the precipitation raster and the county-zone raster align.
compareGeom(r_ny, r_zones)

# Rename the two raster layers used in the following steps.
ppt_ny    <- terra::deepcopy(r_ny)
county_ny <- terra::deepcopy(r_zones)


# 17. Visual check

par(mfrow = c(1, 2))

plot(
  ppt_ny,
  main = "Layer 1: Precipitation (Value)",
  col = map.pal("viridis", 100),
  axes = FALSE
)

plot(
  county_ny,
  main = "Layer 2: County Affiliation (ID)",
  col = sample(rainbow(62)),
  axes = FALSE
)


# 18. Calculate county-level precipitation

calculate_county_precip <- function(r_data, r_zones) {
  
  message("Calculating county-level precipitation...")
  
  # Check whether the precipitation raster and zone raster align.
  if (!compareGeom(r_data, r_zones, stopOnError = FALSE)) {
    warning(
      "The precipitation raster and zone raster do not fully align. ",
      "Please check extent, resolution, origin, and CRS."
    )
  }
  
  # Calculate the physical area of each raster cell.
  # cellSize accounts for geographic CRS distortion when needed.
  pixel_area <- cellSize(r_data)
  
  # Weighted precipitation sum for each pixel.
  # If precipitation is in mm, this intermediate layer has units of mm * m^2.
  precip_vol <- r_data * pixel_area
  
  # Sum weighted precipitation by county zone.
  df_vol <- zonal(precip_vol, r_zones, fun = "sum", na.rm = TRUE)
  names(df_vol) <- c("GEOID", "Total_Volume")
  
  # Sum raster-cell area by county zone.
  df_area <- zonal(pixel_area, r_zones, fun = "sum", na.rm = TRUE)
  names(df_area) <- c("GEOID", "Total_Area")
  
  # Compute area-weighted mean precipitation for each county.
  county_stats <- merge(df_vol, df_area, by = "GEOID")
  county_stats$Precipitation <- county_stats$Total_Volume / county_stats$Total_Area
  
  result_table <- county_stats[, c("GEOID", "Precipitation")]
  
  message("County-level precipitation calculation finished.")
  
  return(result_table)
}

# Calculate the baseline county-level precipitation table.
county_ppt <- calculate_county_precip(ppt_ny, county_ny)

# Save the county-level precipitation table.
write.csv(
  county_ppt,
  file.path(table_dir, "ny_county_precipitation_2023.csv"),
  row.names = FALSE
)


# ============================================================
# 19. Raster merging sensitivity
# Manuscript Example set VI: 1st-3rd order adjacent county merging
# ============================================================

# This analysis randomly performs sequential adjacent county merges on the
# raster-zone layer. After each merge order k = 1, 2, and 3, the county-level
# precipitation values are recalculated and summarized.

res_raster_merging_k3 <- merging_sensitivity_raster(
  r_data           = ppt_ny,
  r_zones          = county_ny,
  summary_function = summary_90th_percentile_raster,
  k_order          = 3,
  n_iterations     = 100,
  merge_method     = "area_weighted",
  random_seed      = 123,
  directions       = 8,
  parallel         = use_parallel,
  n_workers        = parallel_workers,
  save_path        = file.path(rds_dir, "fig_raster_merging_k1_to_k3.rds")
)

plot_merging_sensitivity_raster(
  result_obj = res_raster_merging_k3,
  filename   = file.path(fig_dir, "fig_raster_merging_k1_to_k3.png"),
  bins       = 50,
  width      = 7,
  height     = 16.5,
  dpi        = 300
)


# ============================================================
# 20. Raster splitting sensitivity
# Manuscript Example set VI: 1st-3rd order county splitting
# ============================================================

# This analysis randomly splits eligible raster counties into two parts. After
# each split order k = 1, 2, and 3, the county-level precipitation values are
# recalculated and summarized.

res_raster_splitting_k3 <- splitting_sensitivity_raster(
  r_data           = ppt_ny,
  r_zones          = county_ny,
  summary_function = summary_90th_percentile_raster,
  k_order          = 3,
  n_iterations     = 100,
  split_method     = "area_weighted",
  random_seed      = 123,
  n_attempts       = 10,
  directions       = 8,
  parallel         = use_parallel,
  n_workers        = parallel_workers,
  save_path        = file.path(rds_dir, "fig_raster_splitting_k1_to_k3.rds")
)

plot_splitting_sensitivity_raster(
  result_obj = res_raster_splitting_k3,
  filename   = file.path(fig_dir, "fig_raster_splitting_k1_to_k3.png"),
  bins       = 50,
  width      = 7,
  height     = 16.5,
  dpi        = 300
)
# ============================================================
# 21. Raster pixel reassignment sensitivity
# Manuscript Example set VI: pixel reassignment sensitivity with alpha = 0.2s, 0.5s, and s
# ============================================================

# Here s is the number of pixels covering the average county size.
# For each alpha level, alpha_pixels = round(alpha_multiplier * s).
# The county-level precipitation value is recalculated after pixel
# reassignment, and the summary function is applied to the vector of
# county-level precipitation values.

pixel_alpha_plan <- list(
  fig_raster_pixel_020 = 0.2,
  fig_raster_pixel_050 = 0.5,
  fig_raster_pixel_100 = 1.0
)

pixel_reassignment_results <- vector("list", length(pixel_alpha_plan))
names(pixel_reassignment_results) <- names(pixel_alpha_plan)

for (fig_id in names(pixel_alpha_plan)) {
  
  alpha_multiplier <- pixel_alpha_plan[[fig_id]]
  
  message("Running ", fig_id, ": alpha = ", alpha_multiplier, "s")
  
  pixel_reassignment_results[[fig_id]] <- pixel_reassignment_sensitivity_raster(
    r_data           = ppt_ny,
    r_zones          = county_ny,
    alpha_multiplier = alpha_multiplier,
    summary_function = summary_90th_percentile_raster,
    n_iterations     = 100,
    random_seed      = 123,
    max_failure      = 1000,
    directions       = 4,
    parallel         = use_parallel,
    n_workers        = parallel_workers,
    save_path        = file.path(rds_dir, paste0(fig_id, ".rds"))
  )
}

plot_pixel_reassignment_grid_raster(
  result_list = pixel_reassignment_results,
  filename    = file.path(fig_dir, "fig_raster_pixel_reassignment_alpha_grid.png"),
  bins        = 50,
  width       = 7,
  height      = 16.5,
  dpi         = 300
)



