# On the sensitivities to the modifiable areal unit problem
# YE, Xiang 叶翔; CHEN, Jiayi 陈佳怡
# yexiang@nnu.edu.cn
# 2026-07-02

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
source(file.path(project_dir, "MAUP sensitivity.R"))

dir.create(file.path(project_dir, "outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "outputs", "rds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "outputs", "figures"), recursive = TRUE, showWarnings = FALSE)

rds_dir <- file.path(project_dir, "outputs", "rds")
fig_dir <- file.path(project_dir, "outputs", "figures")

# Parallel Monte Carlo settings.
# Set use_parallel <- FALSE if you want the original serial behavior.
use_parallel <- TRUE
parallel_workers <- max(1L, parallel::detectCores(logical = TRUE) - 5L)

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

col_hh_income_total <- 
  "Estimate!!Total!!Civilian noninstitutionalized population!!HOUSEHOLD INCOME (IN 2023 INFLATION-ADJUSTED DOLLARS)!!Total household population"

# Extract and rename selected S2701 variables
county_attr_s2701 <- s2701_raw %>%
  transmute(
    county_id = str_extract(.data[[col_geo]], "(?<=US)\\d{5}$"),
    pop_civ = parse_number(as.character(.data[[col_pop_civ]])),
    pop_insured = parse_number(as.character(.data[[col_pop_insured]])),
    ins_rate = parse_number(as.character(.data[[col_ins_rate]])) / 100,
    hh_income_total = parse_number(as.character(.data[[col_hh_income_total]])),
    hh_inc_100k_plus = parse_number(as.character(.data[[col_hh_inc_100k_plus]])),
    hh_inc_100k_plus_rate = if_else(
      !is.na(hh_income_total) & hh_income_total > 0,
      hh_inc_100k_plus / hh_income_total,
      NA_real_
    )
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
  "hh_income_total",
  "hh_inc_100k_plus",
  "hh_inc_100k_plus_rate"
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
  hh_income_total = list(
    merge = "sum",
    split = "area_weighted"
  ),
  hh_inc_100k_plus_rate = list(
    merge = list(
      numerator = "hh_inc_100k_plus",
      denominator = "hh_income_total"
    ),
    split = "keep"
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
  # Compute the mean of a numeric vector while ignoring missing values.
  mean(x, na.rm = TRUE)
}

# Median
summary_median <- function(x) {
  # Compute the median of a numeric vector while ignoring missing values.
  median(x, na.rm = TRUE)
}

# Standard deviation
summary_sd <- function(x) {
  # Compute the standard deviation of a numeric vector while ignoring missing values.
  sd(x, na.rm = TRUE)
}

# Interquartile range
summary_iqr <- function(x) {
  # Compute the interquartile range of a numeric vector while ignoring missing values.
  IQR(x, na.rm = TRUE)
}

# Coefficient of variation
summary_cv <- function(x) {
  # Compute the coefficient of variation while ignoring missing values.
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

# Pearson correlation.
# Use with summary_field = c("field_x", "field_y").
summary_pearson <- function(x, y) {
  # Compute the Pearson correlation between two numeric vectors using complete cases.
  cor(x, y, use = "complete.obs", method = "pearson")
}



# ============================================================
# 8. Merging sensitivity
# Manuscript Example set I: first-order merging sensitivity of four summary statistics
# of percentage insured in New York
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
# 8.2 Define merging analysis plan
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
# 8.3 Run first-order merging sensitivity and save plots
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
    parallel             = use_parallel,
    n_workers            = parallel_workers,
    save_path            = file.path(rds_dir, paste0(fig_id, ".rds"))
  )
  
  plot_and_save(
    result_obj      = res,
    filename_prefix = file.path(fig_dir, fig_id)
  )
}


# ============================================================
# 9. Splitting sensitivity
# Manuscript Example set II: first-order splitting sensitivity of the standard deviation
# of percentage insured in NY, PA, OH, and MI
# ============================================================

# ------------------------------------------------------------
# 9.1 Helper function to prepare a state-level county layer
# ------------------------------------------------------------

prepare_state_county <- function(sf_obj, statefp) {
  # Prepare a state-level county layer with valid insurance fields and
  # sequential IDs for sensitivity analysis.
  
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
# 9.2 Create county layers for the four states
# ------------------------------------------------------------

county_ny <- prepare_state_county(county, state_fips[["NY"]])
county_pa <- prepare_state_county(county, state_fips[["PA"]])
county_oh <- prepare_state_county(county, state_fips[["OH"]])
county_mi <- prepare_state_county(county, state_fips[["MI"]])

# ------------------------------------------------------------
# 9.3 Define the splitting analysis plan
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
# 9.4 Run first-order splitting sensitivity and save plots
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
    parallel         = use_parallel,
    n_workers        = parallel_workers,
    save_path        = file.path(rds_dir, paste0(fig_id, ".rds"))
  )
  
  plot_and_save(
    result_obj      = res,
    filename_prefix = file.path(fig_dir, fig_id)
  )
}

# ============================================================
# 10. High-order sensitivity
# Manuscript Example set III: 1st-5th order merging and splitting sensitivities
# of Pearson correlation between percentage insured and percentage high income
# in New York
# ============================================================

# The S2701 table provides counts for total household population and
# household population in the $100,000-and-over category. The percentage
# variable used here is computed as:
# hh_inc_100k_plus_rate = hh_inc_100k_plus / hh_income_total

county_ny_high_order <- county %>%
  filter(STATEFP == "36") %>%
  filter(
    !is.na(pop_civ),
    !is.na(pop_insured),
    !is.na(ins_rate),
    !is.na(hh_income_total),
    !is.na(hh_inc_100k_plus),
    !is.na(hh_inc_100k_plus_rate)
  ) %>%
  mutate(
    ID = row_number(),
    ins_rate = if_else(ins_rate > 1, ins_rate / 100, ins_rate)
  ) %>%
  st_make_valid()

high_order_fields <- c("ins_rate", "hh_inc_100k_plus_rate")
high_order_k <- 1:5

high_order_merging_results <- vector("list", length(high_order_k))
names(high_order_merging_results) <- paste0("k", high_order_k)

for (k in high_order_k) {
  
  message("Running high-order merging sensitivity: k = ", k)
  
  high_order_merging_results[[paste0("k", k)]] <- merging_sensitivity(
    sf                   = county_ny_high_order,
    field_rules          = field_rules,
    summary_field        = high_order_fields,
    summary_function     = summary_pearson,
    k_order              = k,
    exhaustive_threshold = 0,
    n_iterations         = 100,
    random_seed          = 1000 + k,
    parallel             = use_parallel,
    n_workers            = parallel_workers,
    save_path            = file.path(rds_dir, paste0("fig_high_order_merging_k", k, ".rds"))
  )
}

saveRDS(
  high_order_merging_results,
  file.path(rds_dir, "fig_high_order_merging_all_k.rds")
)

high_order_merging_plots <- plot_and_save(
  result_obj      = high_order_merging_results,
  filename_prefix = file.path(fig_dir, "fig_high_order_merging")
)

combine_sensitivity_plots(
  plots    = high_order_merging_plots,
  filename = file.path(fig_dir, "fig_high_order_merging_combined.png"),
  ncol     = 1,
  width    = 7,
  height   = 26.25,
  dpi      = 300,
  device   = "png"
)


high_order_splitting_results <- vector("list", length(high_order_k))
names(high_order_splitting_results) <- paste0("k", high_order_k)

for (k in high_order_k) {
  
  message("Running high-order splitting sensitivity: k = ", k)
  
  high_order_splitting_results[[paste0("k", k)]] <- splitting_sensitivity(
    sf               = county_ny_high_order,
    field_rules      = field_rules,
    summary_field    = high_order_fields,
    summary_function = summary_pearson,
    k_order          = k,
    n_iterations     = 100,
    random_seed      = 2000 + k,
    parallel         = use_parallel,
    n_workers        = parallel_workers,
    save_path        = file.path(rds_dir, paste0("fig_high_order_splitting_k", k, ".rds"))
  )
}

saveRDS(
  high_order_splitting_results,
  file.path(rds_dir, "fig_high_order_splitting_all_k.rds")
)

high_order_splitting_plots <- plot_and_save(
  result_obj      = high_order_splitting_results,
  filename_prefix = file.path(fig_dir, "fig_high_order_splitting")
)

combine_sensitivity_plots(
  plots    = high_order_splitting_plots,
  filename = file.path(fig_dir, "fig_high_order_splitting_combined.png"),
  ncol     = 1,
  width    = 7,
  height   = 26.25,
  dpi      = 300,
  device   = "png"
)
# ============================================================
# 11. Continuous sensitivity
# Manuscript Example set IV: continuous reassignment sensitivity of the coefficient
# of variation of county population in New York
# ============================================================

# ------------------------------------------------------------
# 11.1 Prepare New York county layer
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
# 11.2 Calculate average county area
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
# 11.3 Define continuous sensitivity plan
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
# 11.4 Run continuous sensitivity and save plots
# ------------------------------------------------------------

for (fig_id in names(continuous_plan)) {
  
  message("Running ", fig_id, ": ", continuous_plan[[fig_id]]$label)
  
  res <- continuous_reassignment_sensitivity(
    sf               = county_ny_continuous,
    field_rules      = field_rules,
    summary_field    = "pop_tol",
    summary_function = summary_cv,
    alpha            = continuous_plan[[fig_id]]$alpha,
    n_iterations     = 100,
    max_iter         = 1000,
    tol_ratio        = 0.1,
    random_seed      = 123,
    parallel         = use_parallel,
    n_workers        = parallel_workers,
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
# 12. Discrete sensitivity
# Manuscript Example set V: discrete reassignment sensitivity of Moran's I
# of county population in the NY county-subdivision hierarchy
# ============================================================

# ------------------------------------------------------------
# 12.1 Prepare New York county layer as coarser_sf
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
# 12.2 Prepare county subdivision population attributes from S0101
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
# 12.3 Prepare New York county subdivision layer as finer_sf
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
# 12.4 Define field rules for discrete reassignment
# ------------------------------------------------------------

# pop_tol is an extensive variable.
# After a fine unit is reassigned, county-level population is re-aggregated by sum.
field_rules_discrete <- list(
  pop_tol = list(
    merge = "sum"
  )
)


# ------------------------------------------------------------
# 12.5 Define Moran's I summary function
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
  # Compute Moran's I for the county-level population vector using the fixed
  # neighborhood weights defined above.
  moran.test(
    x,
    listw       = county_lw,
    zero.policy = TRUE,
    na.action   = na.omit
  )$estimate[["Moran I statistic"]]
}


# ------------------------------------------------------------
# 12.6 Calculate area budget
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
# 12.7 Discrete sensitivity by reassigned area
# ------------------------------------------------------------

res_discrete_area <- discrete_reassignment_sensitivity_area(
  finer_sf         = subcounty_ny_discrete,
  coarser_sf       = county_ny_discrete,
  field_rules      = field_rules_discrete,
  n_iterations     = 100,
  alpha            = discrete_alpha,
  summary_field    = "pop_tol",
  summary_function = summary_moran,
  random_seed      = 123,
  parallel         = use_parallel,
  n_workers        = parallel_workers,
  tol_ratio        = 0.1,
  keep_maps        = FALSE,
  save_path        = file.path(rds_dir, "fig_discrete_area.rds")
)

plot_and_save(
  result_obj      = res_discrete_area,
  filename_prefix = file.path(fig_dir, "fig_discrete_area")
)


# ------------------------------------------------------------
# 12.8 Discrete sensitivity by number of reassigned fine regions
# ------------------------------------------------------------

res_discrete_region <- discrete_reassignment_sensitivity_regions(
  finer_sf         = subcounty_ny_discrete,
  coarser_sf       = county_ny_discrete,
  field_rules      = field_rules_discrete,
  k_regions        = 10,
  n_iterations     = 100,
  summary_function = summary_moran,
  summary_field    = "pop_tol",
  random_seed      = 123,
  parallel         = use_parallel,
  n_workers        = parallel_workers,
  save_path        = file.path(rds_dir, "fig_discrete_region.rds")
)

plot_and_save(
  result_obj      = res_discrete_region,
  filename_prefix = file.path(fig_dir, "fig_discrete_region")
)



