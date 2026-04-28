# =============================================================================
# Forest Plot Tree Geolocation Script
# =============================================================================
# Purpose:
#   Reads forest plot survey data from an Excel workbook and produces a
#   combined sf object with precise WGS84 → UTM Zone 17N (EPSG:26917) coords
#   for both plot centres and individual surveyed trees.
#
# Input sheets used:
#   - "Plot Information"  : plot centre lat/lon (WGS84), plot-level attributes
#   - "Tree records"      : per-tree azimuth, distance-to-edge, DBH, and
#                           all other tree attributes
#
# Tree centre geometry:
#   The field measurement gives azimuth (degrees from N, clockwise) and
#   distance to the *edge* of the tree. The tree centre offset is therefore:
#     distance_to_centre = distance_to_edge + (DBH_m / 2)
#   Cartesian offsets in a local UTM metre frame:
#     delta_easting  = distance_to_centre * sin(azimuth_rad)
#     delta_northing = distance_to_centre * cos(azimuth_rad)
#
# Output:
#   forest_plot_spatial.gpkg  — GeoPackage with two layers:
#     "plot_centres"  : one point per plot, all plot-level attributes
#     "trees"         : one point per tree, all tree attributes + plot_id join
#
# Notes on data quality:
#   - One tree (NL_01, tree 62) has azimuth = 1258.6°, which is taken modulo
#     360 (= 178.6°) and flagged in the `azimuth_flag` column.
#   - Two trees have non-numeric DBH values ("10 x <9", "3 x <9"); DBH is set
#     to NA and flagged in `dbh_flag`. Distance offset uses edge distance only
#     (i.e., no DBH/2 correction) for these records.
#   - The "Region" column in Tree records contains a spreadsheet #NAME? error
#     imported as a string; it is retained as-is for the analyst to resolve.
# =============================================================================

library(readxl)
library(dplyr)
library(sf)
library(janitor) # clean_names() for tidy column names

# -----------------------------------------------------------------------------
# 0. File path — adjust if needed
# -----------------------------------------------------------------------------
xlsx_path <- here("data/raw/ES_NRCAN_forest_plots.xlsx")

# =============================================================================
# 1. PLOT CENTRES
# =============================================================================

plots_raw <- read_excel(xlsx_path, sheet = "Plot Information")

plots_sf <- plots_raw |>
  # Standardise column names to snake_case
  clean_names() |>
  # Ensure lat/lon are numeric
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  ) |>
  # Build sf object in WGS84
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |>
  # Reproject to UTM Zone 17N
  st_transform(crs = 26917) |>
  # Extract projected coordinates for later use in tree offset calculations
  mutate(
    plot_easting = st_coordinates(geometry)[, "X"],
    plot_northing = st_coordinates(geometry)[, "Y"]
  )

plot(plots_sf)

# =============================================================================
# 2. TREE RECORDS
# =============================================================================

trees_raw <- read_excel(xlsx_path, sheet = "Tree records")

# Clean and flag data quality issues before computing geometry
trees_clean <- trees_raw |>
  clean_names() |>
  rename(
    azimuth_deg = azimuth_angle_from_n_to_tree_from_plot_centre,
    dist_to_edge = distance_from_plot_centre_m,
    dbh_cm = dbh_cm,
    height_m = height_m
  ) |>
  mutate(
    # ── DBH: coerce to numeric; flag non-numeric entries ──────────────────────
    dbh_flag = if_else(
      is.na(suppressWarnings(as.numeric(dbh_cm))),
      paste0("Non-numeric DBH recorded as: ", dbh_cm),
      NA_character_
    ),
    dbh_cm = suppressWarnings(as.numeric(dbh_cm)),

    # ── Azimuth: coerce to numeric; apply modulo 360 for values > 360 ─────────
    azimuth_deg = as.numeric(azimuth_deg),
    azimuth_flag = if_else(
      azimuth_deg > 360,
      paste0(
        "Original azimuth ",
        azimuth_deg,
        "° wrapped to ",
        azimuth_deg %% 360,
        "°"
      ),
      NA_character_
    ),
    azimuth_deg = azimuth_deg %% 360,

    # ── Distance and offsets ──────────────────────────────────────────────────
    # Convert DBH from cm → m; use 0 if DBH is NA (distance to edge only)
    dbh_m = if_else(is.na(dbh_cm), 0, dbh_cm / 100),
    dist_to_centre = dist_to_edge + (dbh_m / 2),

    # Convert azimuth to radians for trig
    azimuth_rad = azimuth_deg * (pi / 180),

    # UTM offsets: north = 0°, east = 90°, etc.
    delta_easting = dist_to_centre * sin(azimuth_rad),
    delta_northing = dist_to_centre * cos(azimuth_rad)
  )

# =============================================================================
# 3. JOIN PLOT COORDINATES TO TREES
# =============================================================================

# Pull only the columns needed for the spatial join (keep it lean)
plot_coords <- plots_sf |>
  st_drop_geometry() |>
  select(plot, plot_easting, plot_northing)

trees_with_coords <- trees_clean |>
  # Join on plot ID — left join retains all trees, flags unmatched plots
  left_join(plot_coords, by = "plot") |>
  mutate(
    tree_easting = plot_easting + delta_easting,
    tree_northing = plot_northing + delta_northing
  )

# =============================================================================
# 4. BUILD TREE sf OBJECT
# =============================================================================

# Drop working columns not needed in the final output
drop_cols <- c(
  "dbh_m",
  "dist_to_centre",
  "azimuth_rad",
  "delta_easting",
  "delta_northing",
  "plot_easting",
  "plot_northing"
)

trees_sf <- trees_with_coords |>
  # Build sf object directly from the computed UTM coordinates
  st_as_sf(
    coords = c("tree_easting", "tree_northing"),
    crs = 26917,
    remove = FALSE
  ) |>
  # Keep the coordinate columns (useful for QA) but drop the working ones
  select(-any_of(drop_cols))

# =============================================================================
# 5. CLEAN UP PLOT CENTRES (drop working columns)
# =============================================================================

plots_sf <- plots_sf |>
  select(-plot_easting, -plot_northing)

plot_area_sf <- plots_sf |>
  st_buffer(16) # creates 30 m circular plot

# =============================================================================
# 6. WRITE OUTPUT GEOPACKAGE
# =============================================================================

output_path <- here("data/processed/ES_forest_plot_spatial.gpkg")

st_write(
  plots_sf,
  output_path,
  layer = "plot_centres",
  driver = "GPKG",
  delete_dsn = TRUE, # overwrite existing file on re-runs
  quiet = FALSE
)

st_write(
  trees_sf,
  output_path,
  layer = "trees",
  driver = "GPKG",
  append = TRUE, # add second layer to same file
  quiet = FALSE
)

st_write(
  plot_area_sf,
  output_path,
  layer = "areas",
  driver = "GPKG",
  append = TRUE, # add third layer to same file
  quiet = FALSE
)

# =============================================================================
# 7. QUICK SUMMARY
# =============================================================================

message("\n── Output written to: ", output_path)
message("   Layer 'plot_centres' : ", nrow(plots_sf), " plots")
message("   Layer 'trees'        : ", nrow(trees_sf), " trees")
message("   CRS (both layers)    : EPSG:26917 (UTM Zone 17N)")

# Report data quality flags
n_dbh_flag <- sum(!is.na(trees_sf$dbh_flag))
n_az_flag <- sum(!is.na(trees_sf$azimuth_flag))

if (n_dbh_flag > 0) {
  message(
    "\n⚠ DBH flags (",
    n_dbh_flag,
    " trees — DBH set to NA, edge distance used):"
  )
  trees_sf |>
    st_drop_geometry() |>
    filter(!is.na(dbh_flag)) |>
    select(plot, tree, dbh_flag) |>
    print()
}

if (n_az_flag > 0) {
  message("\n⚠ Azimuth flags (", n_az_flag, " trees — wrapped to 0–360°):")
  trees_sf |>
    st_drop_geometry() |>
    filter(!is.na(azimuth_flag)) |>
    select(plot, tree, azimuth_flag) |>
    print()
}
