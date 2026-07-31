# On the sensitivities to the modifiable areal unit problem
# YE, Xiang 叶翔; CHEN, Jiayi 陈佳怡
# yexiang@nnu.edu.cn
# 2026-05-09

# 0. Setup
library(geojsonio)
library(readxl)
library(tidyverse)
library(sf)
library(lwgeom)
library(ggplot2)
library(spdep)

# Set project directory
project_dir <- "D:/Change_to_your_local_path"

# Source the core sensitivity functions and plotting function
source(file.path(project_dir, "2026-04-30 MAUP sensitivity.R"))

dir.create(file.path(project_dir, "outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "outputs", "rds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "outputs", "figures"), recursive = TRUE, showWarnings = FALSE)

s2701_path <- file.path(project_dir, "data", "ACSST5Y2023.S2701-Data.xlsx")
s0101_path <- file.path(project_dir, "data", "ACSST5Y2023.S0101-Data.xlsx")

# 1. Load spatial data
# Replace these file names with the actual shapefile names in your data folder
county_shp <- file.path(project_dir, "data", "US_county", "tl_2023_us_county.shp")
subcounty_shp <- file.path(project_dir, "data", "US_subcounty", "tl_2023_36_cousub.shp")

# State FIPS codes used in the manuscript
state_fips <- c(
  NY = "36",
  PA = "42",
  OH = "39",
  MI = "26"
)

# Read county shapefile first
county <- st_read(county_shp, quiet = TRUE)

# Filter to the four study states before geometry operations
county <- county %>%
  filter(STATEFP %in% state_fips) %>%
  st_make_valid() %>%
  st_transform(6350)

# Read New York county subdivision shapefile
subcounty <- st_read(subcounty_shp, quiet = TRUE)

# This file is already New York only, so transform it directly
subcounty <- subcounty %>%
  st_make_valid() %>%
  st_transform(6350)

# 2. Load ACS attribute tables
# Define ACS CSV paths
s2701_path <- file.path(project_dir, "data", "ACSST5Y2023.S2701-Data.csv")
s0101_path <- file.path(project_dir, "data", "ACSST5Y2023.S0101-Data.csv")

# Read ACS CSV tables
s2701_raw <- read_csv(
  s2701_path,
  col_types = cols(.default = col_character()),
  name_repair = "minimal",
  show_col_types = FALSE
)

s0101_raw <- read_csv(
  s0101_path,
  col_types = cols(.default = col_character()),
  name_repair = "minimal",
  show_col_types = FALSE
)

# Remove possible BOM characters from column names
names(s2701_raw) <- gsub("^\ufeff", "", names(s2701_raw))
names(s0101_raw) <- gsub("^\ufeff", "", names(s0101_raw))


# 3. Clean ACS S2701 variables
# Geography identifier column
col_geo <- "Geography"

# Original ACS field names in S2701
col_pop_civ <- 
  "Estimate!!Total!!Civilian noninstitutionalized population"

col_pop_insured <- 
  "Estimate!!Insured!!Civilian noninstitutionalized population"

col_ins_rate <- 
  "Estimate!!Percent Insured!!Civilian noninstitutionalized population"

col_hh_inc_100k_plus <- 
  "Estimate!!Total!!Civilian noninstitutionalized population!!HOUSEHOLD INCOME (IN 2023 INFLATION-ADJUSTED DOLLARS)!!Total household population!!$100,000 and over"

# Extract and rename selected S2701 variables
county_attr_s2701 <- s2701_raw %>%
  transmute(
    county_id = str_extract(.data[[col_geo]], "(?<=US)\\d{5}$"),
    pop_civ = parse_number(as.character(.data[[col_pop_civ]])),
    pop_insured = parse_number(as.character(.data[[col_pop_insured]])),
    ins_rate = parse_number(as.character(.data[[col_ins_rate]])) / 100,
    hh_inc_100k_plus = parse_number(as.character(.data[[col_hh_inc_100k_plus]]))
  )


# 4. Clean ACS S0101 total population variable

# Original ACS field name in S0101
col_pop_tol <- "Estimate!!Total!!Total population"

# Extract total population and aggregate subcounty data to the county level
county_attr_s0101 <- s0101_raw %>%
  mutate(
    # Remove the '$' anchor to extract the first 5 digits (State + County FIPS) after 'US'
    county_id = str_extract(.data[[col_geo]], "(?<=US)\\d{5}"),
    pop_tol = parse_number(as.character(.data[[col_pop_tol]]))
  ) %>%
  # Filter out invalid rows (e.g., the second header row)
  filter(!is.na(county_id)) %>%
  # Group by the 5-digit county ID
  group_by(county_id) %>%
  # Sum the population of all subcounties within the same county
  summarise(
    pop_tol = sum(pop_tol, na.rm = TRUE),
    .groups = "drop"
  )


# 5. Join ACS attributes to the county shapefile

# The shapefile usually contains a GEOID field.
# If your shapefile uses a different field name, replace GEOID below.
county <- county %>%
  mutate(
    county_id = as.character(GEOID)
  ) %>%
  left_join(county_attr_s2701, by = "county_id") %>%
  left_join(county_attr_s0101, by = "county_id") %>%
  mutate(
    ID = row_number()
  )

# Check whether all required variables exist
required_cols <- c(
  "ID",
  "county_id",
  "pop_tol",
  "pop_civ",
  "pop_insured",
  "ins_rate",
  "hh_inc_100k_plus"
)

missing_cols <- setdiff(required_cols, names(county))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}


# 6. Define value recalculation rules

# Extensive variables are summed when regions are merged.
# During splitting, they are allocated by area under the homogeneity assumption.
# Rate variables are recomputed after merging and kept unchanged after splitting.
field_rules <- list(
  pop_tol = list(
    merge = "sum",
    split = "area_weighted"
  ),
  pop_civ = list(
    merge = "sum",
    split = "area_weighted"
  ),
  pop_insured = list(
    merge = "sum",
    split = "area_weighted"
  ),
  hh_inc_100k_plus = list(
    merge = "sum",
    split = "area_weighted"
  ),
  ins_rate = list(
    merge = list(
      numerator = "pop_insured",
      denominator = "pop_civ"
    ),
    split = "keep"
  )
)

# 7. Define summary functions

# Mean
summary_mean <- function(x) {
  mean(x, na.rm = TRUE)
}

# Median
summary_median <- function(x) {
  median(x, na.rm = TRUE)
}

# Standard deviation
summary_sd <- function(x) {
  sd(x, na.rm = TRUE)
}

# Interquartile range
summary_iqr <- function(x) {
  IQR(x, na.rm = TRUE)
}

# Coefficient of variation
summary_cv <- function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}



# ============================================================
# 8. Merging sensitivity
# First-order merging sensitivity of insurance rate in New York
# ============================================================

# ------------------------------------------------------------
# 8.1 Prepare New York county data
# ------------------------------------------------------------

county_ny <- county %>%
  filter(STATEFP == "36") %>%
  filter(
    !is.na(pop_civ),
    !is.na(pop_insured),
    !is.na(ins_rate),
    !is.na(hh_inc_100k_plus)
  ) %>%
  mutate(
    ID = row_number(),
    ins_rate = if_else(ins_rate > 1, ins_rate / 100, ins_rate)
  ) %>%
  st_make_valid()

# Check the number of New York counties
message("Number of New York counties used for merging sensitivity: ", nrow(county_ny))


# ------------------------------------------------------------
# 8.2 Define output folders
# ------------------------------------------------------------

rds_dir <- file.path(project_dir, "outputs", "rds")
fig_dir <- file.path(project_dir, "outputs", "figures")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 8.3 Optional fallback for calculate_merge_complexity()
# ------------------------------------------------------------

# Some versions of the core function file call calculate_merge_complexity()
# inside merging_sensitivity(). If it is not defined in the sourced file,
# this fallback gives a conservative estimate based on rook-adjacent pairs.
if (!exists("calculate_merge_complexity")) {
  
  calculate_merge_complexity <- function(sf, n_simulations = 200) {
    
    nb <- poly2nb(sf, queen = FALSE)
    n_pairs <- sum(card(nb)) / 2
    
    data.frame(
      k_order = 1:10,
      estimated_total_combinations = n_pairs ^ (1:10)
    )
  }
}


# ------------------------------------------------------------
# 8.4 Define merging analysis plan
# ------------------------------------------------------------

merging_plan <- list(
  fig_merging_mean = list(
    summary_function = summary_mean,
    label = "mean"
  ),
  fig_merging_median = list(
    summary_function = summary_median,
    label = "median"
  ),
  fig_merging_sd = list(
    summary_function = summary_sd,
    label = "standard deviation"
  ),
  fig_merging_iqr = list(
    summary_function = summary_iqr,
    label = "interquartile range"
  )
)


# ------------------------------------------------------------
# 8.5 Run first-order merging sensitivity and save plots
# ------------------------------------------------------------

for (fig_id in names(merging_plan)) {
  
  message("Running ", fig_id, ": ", merging_plan[[fig_id]]$label)
  
  res <- merging_sensitivity(
    sf                   = county_ny,
    field_rules          = field_rules,
    summary_field        = "ins_rate",
    summary_function     = merging_plan[[fig_id]]$summary_function,
    k_order              = 1,
    exhaustive_threshold = Inf,
    n_iterations         = 100,
    random_seed          = 123,
    save_path            = file.path(rds_dir, paste0(fig_id, ".rds"))
  )
  
  plot_and_save(
    result_obj      = res,
    filename_prefix = file.path(fig_dir, fig_id)
  )
}


# ============================================================
# 9. Splitting sensitivity
# First-order splitting sensitivity of the standard deviation
# of insurance rate in NY, PA, OH, and MI
# ============================================================

# ------------------------------------------------------------
# 9.1 Define state FIPS codes
# ------------------------------------------------------------

state_fips <- c(
  NY = "36",
  PA = "42",
  OH = "39",
  MI = "26"
)

# ------------------------------------------------------------
# 9.2 Helper function to prepare a state-level county layer
# ------------------------------------------------------------

prepare_state_county <- function(sf_obj, statefp) {
  
  sf_obj %>%
    filter(STATEFP == statefp) %>%
    filter(
      !is.na(pop_civ),
      !is.na(pop_insured),
      !is.na(ins_rate)
    ) %>%
    mutate(
      ID = row_number(),
      ins_rate = if_else(ins_rate > 1, ins_rate / 100, ins_rate)
    ) %>%
    st_make_valid()
}

# ------------------------------------------------------------
# 9.3 Create county layers for the four states
# ------------------------------------------------------------

county_ny <- prepare_state_county(county, state_fips[["NY"]])
county_pa <- prepare_state_county(county, state_fips[["PA"]])
county_oh <- prepare_state_county(county, state_fips[["OH"]])
county_mi <- prepare_state_county(county, state_fips[["MI"]])

# ------------------------------------------------------------
# 9.4 Define the splitting analysis plan
# ------------------------------------------------------------

splitting_plan <- list(
  fig_splitting_ny = list(
    sf    = county_ny,
    label = "New York"
  ),
  fig_splitting_pa = list(
    sf    = county_pa,
    label = "Pennsylvania"
  ),
  fig_splitting_oh = list(
    sf    = county_oh,
    label = "Ohio"
  ),
  fig_splitting_mi = list(
    sf    = county_mi,
    label = "Michigan"
  )
)

# ------------------------------------------------------------
# 9.5 Run first-order splitting sensitivity and save plots
# ------------------------------------------------------------

for (fig_id in names(splitting_plan)) {
  
  message("Running ", fig_id, ": ", splitting_plan[[fig_id]]$label)
  
  sf_state <- splitting_plan[[fig_id]]$sf
  
  if (nrow(sf_state) == 0) {
    warning("No valid county records found for ", splitting_plan[[fig_id]]$label, ". Skipping.")
    next
  }
  
  res <- splitting_sensitivity(
    sf               = sf_state,
    field_rules      = field_rules,
    summary_field    = "ins_rate",
    summary_function = summary_sd,
    k_order          = 1,
    n_iterations     = 100,
    random_seed      = 123,
    save_path        = file.path(rds_dir, paste0(fig_id, ".rds"))
  )
  
  plot_and_save(
    result_obj      = res,
    filename_prefix = file.path(fig_dir, fig_id)
  )
}

# ============================================================
# 10. Continuous sensitivity
# Coefficient of variation of county population in New York
# ============================================================

# ------------------------------------------------------------
# 10.1 Define output folders
# ------------------------------------------------------------

rds_dir <- file.path(project_dir, "outputs", "rds")
fig_dir <- file.path(project_dir, "outputs", "figures")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 10.2 Prepare New York county layer
# ------------------------------------------------------------

county_ny_continuous <- county %>%
  filter(STATEFP == "36") %>%
  filter(
    !is.na(pop_tol),
    !is.na(pop_civ),
    !is.na(pop_insured),
    !is.na(ins_rate)
  ) %>%
  mutate(
    ID = row_number(),
    ins_rate = if_else(ins_rate > 1, ins_rate / 100, ins_rate)
  ) %>%
  st_make_valid()

if (nrow(county_ny_continuous) == 0) {
  stop("county_ny_continuous has 0 rows. Please check whether pop_tol was successfully joined.")
}

message("Number of New York counties used for continuous sensitivity: ", nrow(county_ny_continuous))


# ------------------------------------------------------------
# 10.3 Calculate average county area
# ------------------------------------------------------------

county_area_m2 <- as.numeric(st_area(county_ny_continuous))

avg_county_area_m2 <- mean(
  county_area_m2[is.finite(county_area_m2) & county_area_m2 > 0],
  na.rm = TRUE
)

if (!is.finite(avg_county_area_m2) || avg_county_area_m2 <= 0) {
  stop("Average county area is invalid. Please check geometry and CRS.")
}

avg_county_area_km2 <- avg_county_area_m2 / 1e6

message("Average New York county area: ", round(avg_county_area_km2, 2), " km^2")


# ------------------------------------------------------------
# 10.4 Define continuous sensitivity plan
# ------------------------------------------------------------

continuous_plan <- list(
  fig_continuous_020 = list(
    alpha = 0.2 * avg_county_area_m2,
    label = "20 percent of average county area"
  ),
  fig_continuous_050 = list(
    alpha = 0.5 * avg_county_area_m2,
    label = "50 percent of average county area"
  ),
  fig_continuous_100 = list(
    alpha = 1.0 * avg_county_area_m2,
    label = "100 percent of average county area"
  ),
  fig_continuous_300 = list(
    alpha = 3.0 * avg_county_area_m2,
    label = "300 percent of average county area"
  )
)


# ------------------------------------------------------------
# 10.5 Run continuous sensitivity and save plots
# ------------------------------------------------------------

for (fig_id in names(continuous_plan)) {
  
  message("Running ", fig_id, ": ", continuous_plan[[fig_id]]$label)
  
  res <- continuous_sensitivity(
    sf               = county_ny_continuous,
    field_rules      = field_rules,
    summary_field    = "pop_tol",
    summary_function = summary_cv,
    alpha            = continuous_plan[[fig_id]]$alpha,
    n_iterations     = 100,
    max_iter         = 1000,
    tol_ratio        = 0.1,
    random_seed      = 123,
    save_path        = file.path(rds_dir, paste0(fig_id, ".rds"))
  )
  
  if (nrow(res$distribution) == 0) {
    warning("No valid result collected for ", fig_id, ". Plot skipped.")
    next
  }
  
  plot_and_save(
    result_obj      = res,
    filename_prefix = file.path(fig_dir, fig_id)
  )
}


# ============================================================
# 11. Discrete sensitivity
# Moran's I of county population in the NY county-subdivision hierarchy
# ============================================================

# ------------------------------------------------------------
# 11.1 Prepare New York county layer as coarser_sf
# ------------------------------------------------------------

county_ny_discrete <- county %>%
  filter(STATEFP == "36") %>%
  filter(!is.na(pop_tol)) %>%
  mutate(
    ID = row_number()
  ) %>%
  st_make_valid()

# Lookup table: county_id -> coarse_ID used by discrete sensitivity functions
county_id_lookup <- county_ny_discrete %>%
  st_drop_geometry() %>%
  select(
    county_id,
    coarse_ID = ID
  )


# ------------------------------------------------------------
# 11.2 Prepare county subdivision population attributes from S0101
# ------------------------------------------------------------

subcounty_attr_s0101 <- s0101_raw %>%
  transmute(
    subcounty_id = str_extract(.data[[col_geo]], "(?<=US)\\d+"),
    pop_tol = parse_number(as.character(.data[[col_pop_tol]]))
  ) %>%
  filter(
    !is.na(subcounty_id),
    !is.na(pop_tol)
  )


# ------------------------------------------------------------
# 11.3 Prepare New York county subdivision layer as finer_sf
# ------------------------------------------------------------

subcounty_ny_discrete <- subcounty %>%
  mutate(
    subcounty_id = as.character(GEOID),
    county_id = str_sub(subcounty_id, 1, 5)
  ) %>%
  left_join(subcounty_attr_s0101, by = "subcounty_id") %>%
  left_join(county_id_lookup, by = "county_id") %>%
  filter(
    !is.na(pop_tol),
    !is.na(coarse_ID)
  ) %>%
  mutate(
    ID = row_number()
  ) %>%
  st_make_valid()

message("Number of NY counties: ", nrow(county_ny_discrete))
message("Number of NY county subdivisions: ", nrow(subcounty_ny_discrete))


# ------------------------------------------------------------
# 11.4 Define field rules for discrete reassignment
# ------------------------------------------------------------

# pop_tol is an extensive variable.
# After a fine unit is reassigned, county-level population is re-aggregated by sum.
field_rules_discrete <- list(
  pop_tol = list(
    merge = "sum"
  )
)


# ------------------------------------------------------------
# 11.5 Define Moran's I summary function
# ------------------------------------------------------------

# Construct rook contiguity weights for the original NY county layer.
# The discrete sensitivity functions keep the number of counties unchanged,
# so this fixed listw can be used with the county-level population vector.
neighbor_nb <- poly2nb(county_ny_discrete, queen = FALSE)

county_lw <- nb2listw(
  neighbor_nb,
  style       = "W",
  zero.policy = TRUE
)

summary_moran <- function(x) {
  moran.test(
    x,
    listw       = county_lw,
    zero.policy = TRUE,
    na.action   = na.omit
  )$estimate[["Moran I statistic"]]
}


# ------------------------------------------------------------
# 11.6 Calculate area budget
# ------------------------------------------------------------

county_area_m2 <- as.numeric(st_area(county_ny_discrete))

avg_county_area_m2 <- mean(
  county_area_m2[is.finite(county_area_m2) & county_area_m2 > 0],
  na.rm = TRUE
)

avg_county_area_km2 <- avg_county_area_m2 / 1e6

message("Average NY county area: ", round(avg_county_area_km2, 2), " km^2")

# Area budget for discrete area sensitivity.
# Use 1.0 * avg_county_area_m2 if you want one average county area.
discrete_alpha <- 1.0 * avg_county_area_m2


# ------------------------------------------------------------
# 11.7 Discrete sensitivity by reassigned area
# ------------------------------------------------------------

res_discrete_area <- discrete_area_sensitivity(
  finer_sf         = subcounty_ny_discrete,
  coarser_sf       = county_ny_discrete,
  field_rules      = field_rules_discrete,
  n_iterations     = 100,
  alpha            = discrete_alpha,
  summary_field    = "pop_tol",
  summary_function = summary_moran,
  tol_ratio        = 0.1,
  keep_maps        = FALSE,
  save_path        = file.path(rds_dir, "fig_discrete_area.rds")
)

plot_and_save(
  result_obj      = res_discrete_area,
  filename_prefix = file.path(fig_dir, "fig_discrete_area")
)


# ------------------------------------------------------------
# 11.8 Discrete sensitivity by number of reassigned fine regions
# ------------------------------------------------------------

res_discrete_region <- discrete_region_sensitivity(
  finer_sf         = subcounty_ny_discrete,
  coarser_sf       = county_ny_discrete,
  field_rules      = field_rules_discrete,
  k_regions        = 10,
  n_iterations     = 100,
  summary_function = summary_moran,
  summary_field    = "pop_tol",
  save_path        = file.path(rds_dir, "fig_discrete_region.rds")
)

plot_and_save(
  result_obj      = res_discrete_region,
  filename_prefix = file.path(fig_dir, "fig_discrete_region")
)
