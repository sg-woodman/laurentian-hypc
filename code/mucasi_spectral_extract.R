# =============================================================================
# µCASI Hyperspectral Point Cloud -- Spectral Profile Extractor
# =============================================================================
# Purpose:  Memory-efficient extraction of full spectral profiles (all bands)
#           for points within a polygon or within a buffer around a point.
#           Uses a two-pass strategy:
#             Pass 1 -- load only X, Y to identify candidate points (cheap)
#             Pass 2 -- reload only matched rows with all 238 bands (targeted)
#
#           This avoids ever loading all 238 bands for the full point cloud,
#           making it practical even on 16 GB RAM systems with 28+ GB datasets.
#
# Outputs:  1. <stem>_spectra_raw.csv    -- one row per point, all band values
#           2. <stem>_spectra_mean.csv   -- mean spectrum across all matched points
#           3. <stem>_spectra_plot.png   -- spectral profile plot
#           4. <stem>_spectra_points.gpkg -- spatial locations of extracted points
#
# Requires: pts data.table and wave_lookup tibble from process_mucasi_cloudcompare.R
#           OR run standalone by setting input_dir / file_stem below.
# =============================================================================

library(data.table)
library(tidyverse)
library(sf)

# =============================================================================
# -- CONFIGURATION -------------------------------------------------------------
# =============================================================================

# Directory containing your µCASI .txt files
input_dir <- "/Users/sam/Documents/cfs/laurentian-recovery/data/raw/hspc/"

# Output directory
output_dir <- "/Users/sam/Documents/cfs/laurentian-recovery/data/processed/"

# File stem -- matches a single flight line or site prefix
# Use the same value as in process_mucasi_cloudcompare.R
file_stem <- "DaisyLake"
output_stem <- "DaisyLake_spectra"

# Short identifier to distinguish outputs (e.g. site, catchment, date, run)
# Prepended to all output filenames: <id>_<output_stem>_raw.csv etc.
# Set to "" to omit (no prefix added)
extract_id <- "catch_j_"

# EPSG code for your point cloud CRS
pointcloud_epsg <- 32617 # UTM Zone 17N WGS84 -- typical for Sudbury

# --- Extraction mode ----------------------------------------------------------
# "polygon" -- extract all points within a polygon (shapefile, gpkg, geojson)
# "point"   -- extract points within a buffer radius around one or more points
extraction_mode <- "point"

catch_i <- vect("/Users/sam/Downloads/Daisy Subwatersheds/Subcatchment_I.shp")
catch_j <- vect("/Users/sam/Downloads/Daisy Subwatersheds/Subcatchment J.shp")

catch_ij <- bind_spat_rows(catch_i, catch_j)
writeVector(catch_ij, "/Users/sam/Downloads/catchment_IJ.gpkg")

# For "polygon" mode: path to polygon file
polygon_path <- "/Users/sam/Downloads/catchment_IJ.gpkg"
polygon_reproject <- TRUE # set TRUE if polygon is in a different CRS

# For "point" mode: coordinates (in the same CRS as point cloud) and buffer
# Can be a single point or a data frame of multiple points with an id column
extract_points <- data.frame(
  id = c("catch_k", "catch_j", "catch_i", "catch_d", "catch_c"),
  X  = c(509654.51, 509255.52, 509101.14, 507917.39, 507618.04), # replace with your actual coordinates
  Y  = c(5145079.74, 5145079.74, 5144847.58, 5144032.87, 5143731.14)
)


buffer_radius <- 3 # metres -- radius around each point to extract

# --- Spectral profile plot options -------------------------------------------
# Which points to include in the mean spectrum plot:
# "all"     -- single mean spectrum across all extracted points
# "by_id"   -- one mean spectrum per polygon feature or point id (if multiple)
plot_grouping <- "by_id"

# Smooth the spectral profile with a rolling mean (reduces sensor noise)
# Set to 1 to disable smoothing, or an odd integer e.g. 5, 7, 9
smooth_window <- 1

# Wavelength range to plot (nm) -- set to NULL to use full range
plot_wavelength_min <- NULL # e.g. 500
plot_wavelength_max <- NULL # e.g. 850

# =============================================================================
# -- FUNCTIONS -----------------------------------------------------------------
# =============================================================================

#' Parse µCASI header (same as main script -- reproduced here for standalone use)
parse_mucasi_header <- function(filepath) {
  lines <- readLines(filepath, n = 2)
  wavelengths <- as.numeric(
    strsplit(trimws(sub("wavelengths_nm:", "", lines[1])), "\\s+")[[1]]
  )
  wavelengths <- wavelengths[!is.na(wavelengths)]
  col_names <- strsplit(trimws(sub("#\\s*columns:", "", lines[2])), "\\s+")[[1]]
  band_names <- col_names[4:length(col_names)]
  wave_lookup <- tibble(band = band_names, wavelength_nm = wavelengths)
  list(col_names = col_names, wavelengths = wavelengths, wave_lookup = wave_lookup)
}


#' Detect flight line parts matching a file stem
detect_parts <- function(input_dir, file_stem) {
  all_txt <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)
  matched <- sort(all_txt[grepl(file_stem, basename(all_txt))])
  if (length(matched) == 0) stop("No files found for stem: ", file_stem)
  matched
}


#' Pass 1: load only X and Y from a single part file
#' Returns a data.table with columns: X, Y, .part (file index), .row (row in file)
#' .part and .row together uniquely identify each point for the second pass
load_xy_only <- function(filepath, col_names, part_index) {
  # Read only X and Y columns by position (col 1 and 2 after skip)
  # Using select= with column positions is much faster than reading all cols
  xy <- fread(filepath,
    header       = FALSE,
    skip         = 2,
    select       = c(1L, 2L), # X is col 1, Y is col 2
    col.names    = c("X", "Y"),
    showProgress = FALSE
  )

  xy[, `:=`(
    .part = part_index,
    .row = .I
  )] # .I = row number within this file

  xy
}


#' Build spatial filter from extraction mode settings
#' Returns a list with: bbox (named vector), polygon_sf (sf object or NULL)
build_spatial_filter <- function(extraction_mode, pointcloud_epsg) {
  pc_crs <- st_crs(pointcloud_epsg)

  if (extraction_mode == "polygon") {
    if (!file.exists(polygon_path)) stop("Polygon not found: ", polygon_path)
    poly <- st_read(polygon_path, quiet = TRUE)

    if (nrow(poly) > 1) {
      cat(sprintf(
        "  %d polygon features detected -- keeping separate for by_id grouping\n",
        nrow(poly)
      ))
    }

    if (polygon_reproject) {
      poly <- st_transform(poly, crs = pc_crs)
    } else if (is.na(st_crs(poly))) {
      poly <- st_set_crs(poly, pc_crs)
    }

    # Add an id column if not present (used for by_id grouping)
    if (!"id" %in% names(poly)) poly$id <- as.character(seq_len(nrow(poly)))

    # Use union for bbox (covers all features), keep original for intersection
    bbox <- st_bbox(st_union(poly))
    return(list(bbox = bbox, filter_sf = poly))
  } else if (extraction_mode == "point") {
    # Convert extract_points to sf with buffer
    pts_sf <- st_as_sf(extract_points, coords = c("X", "Y"), crs = pc_crs)
    buffered <- st_buffer(pts_sf, dist = buffer_radius)

    bbox <- st_bbox(st_union(buffered))
    return(list(bbox = bbox, filter_sf = buffered))
  } else {
    stop("extraction_mode must be 'polygon' or 'point'")
  }
}


#' Rolling mean smoother for spectral profiles
#' @param x   Numeric vector (reflectance values across bands)
#' @param w   Window size (must be odd)
#' @return Smoothed numeric vector, same length as x
smooth_spectrum <- function(x, w = smooth_window) {
  if (w <= 1) {
    return(x)
  }
  stats::filter(x, rep(1 / w, w), sides = 2) %>%
    as.numeric() %>%
    # filter() returns NA at edges -- fill with original values
    {
      ifelse(is.na(.), x, .)
    }
}


# =============================================================================
# -- MAIN EXTRACTION PIPELINE --------------------------------------------------
# =============================================================================

cat("\n=== µCASI Spectral Profile Extractor ===\n\n")
cat(sprintf("Mode: %s | Stem: %s\n\n", extraction_mode, file_stem))

# Build output filename prefix from extract_id and output_stem.
# If extract_id is empty, use output_stem alone (no leading underscore).
file_prefix <- if (nzchar(extract_id)) {
  paste0(extract_id, "_", output_stem)
} else {
  output_stem
}
cat(sprintf("Output prefix: %s\n\n", file_prefix))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -- 1. Find input files and parse header -------------------------------------
part_files <- detect_parts(input_dir, file_stem)
n_parts <- length(part_files)

cat(sprintf("Found %d part file(s)\n", n_parts))

cat("Parsing header...\n")
header <- parse_mucasi_header(part_files[1])
col_names <- header$col_names
wave_lookup <- header$wave_lookup
band_cols <- wave_lookup$band # all 238 band column names
cat(sprintf(
  "  %d bands, %.1f-%.1f nm\n\n",
  nrow(wave_lookup), min(wave_lookup$wavelength_nm),
  max(wave_lookup$wavelength_nm)
))

# -- 2. Build spatial filter --------------------------------------------------
cat("Building spatial filter...\n")
spatial_filter <- build_spatial_filter(extraction_mode, pointcloud_epsg)
bbox <- spatial_filter$bbox
filter_sf <- spatial_filter$filter_sf

cat(sprintf(
  "  Filter extent: X [%.0f, %.0f], Y [%.0f, %.0f]\n\n",
  bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]
))

# -- 3. PASS 1: Load X, Y only and find candidate rows -------------------------
# This is very fast -- we only load 2 of 241 columns across the full dataset.
# We store (.part, .row) pairs that pass both the bbox and polygon tests.
cat("Pass 1: scanning X, Y coordinates across all files...\n")

candidate_list <- vector("list", n_parts)
pc_crs <- st_crs(pointcloud_epsg)

for (i in seq_along(part_files)) {
  cat(sprintf("  [%d/%d] %s\n", i, n_parts, basename(part_files[i])))

  # Load X, Y only
  xy <- load_xy_only(part_files[i], col_names, part_index = i)

  # Fast bbox pre-filter
  xy_bbox <- xy[X >= bbox["xmin"] & X <= bbox["xmax"] &
    Y >= bbox["ymin"] & Y <= bbox["ymax"]]

  if (nrow(xy_bbox) == 0) {
    cat(sprintf("    -> 0 points in bounding box\n"))
    candidate_list[[i]] <- data.table()
    rm(xy, xy_bbox)
    gc(verbose = FALSE)
    next
  }

  # Precise polygon/buffer intersection
  xy_sf <- st_as_sf(xy_bbox, coords = c("X", "Y"), crs = pc_crs, remove = FALSE)
  xy_in <- st_join(xy_sf, filter_sf, join = st_within, left = FALSE)

  if (nrow(xy_in) == 0) {
    cat(sprintf("    -> 0 points within polygon/buffer\n"))
    candidate_list[[i]] <- data.table()
  } else {
    # Keep (.part, .row, X, Y) plus the id from the polygon/buffer for grouping
    result <- xy_in %>%
      st_drop_geometry() %>%
      as.data.table() %>%
      select(any_of(c(".part", ".row", "X", "Y", "id")))

    cat(sprintf("    -> %s candidate points\n", format(nrow(result), big.mark = ",")))
    candidate_list[[i]] <- result
  }

  rm(xy, xy_bbox, xy_sf, xy_in)
  gc(verbose = FALSE)
}

candidates <- rbindlist(candidate_list, fill = TRUE)
rm(candidate_list)
gc(verbose = FALSE)

n_candidates <- nrow(candidates)
cat(sprintf(
  "\nPass 1 complete: %s candidate points found\n\n",
  format(n_candidates, big.mark = ",")
))

if (n_candidates == 0) {
  stop(
    "No points found within the specified polygon/buffer. ",
    "Check CRS alignment and coordinate values."
  )
}

# -- 4. PASS 2: Reload only candidate rows with all 238 bands -----------------
# We split candidates by part file, then use the stored .row indices to
# load only those rows. fread does not support arbitrary row selection,
# so we read the full file once more but immediately filter to candidate rows.
# This is still efficient because we only do it for files that had candidates.
cat("Pass 2: loading full spectra for candidate points...\n")
cat(sprintf(
  "  (Reading all 238 bands for %s points only)\n\n",
  format(n_candidates, big.mark = ",")
))

all_cols_to_load <- c("X", "Y", "Z", band_cols)
spectra_list <- vector("list", n_parts)

for (i in seq_along(part_files)) {
  part_candidates <- candidates[.part == i]

  if (nrow(part_candidates) == 0) next

  cat(sprintf(
    "  [%d/%d] %s -- extracting %s rows\n",
    i, n_parts, basename(part_files[i]),
    format(nrow(part_candidates), big.mark = ",")
  ))

  # Read full file with all columns
  full_part <- fread(part_files[i],
    header       = FALSE,
    skip         = 2,
    col.names    = col_names,
    showProgress = FALSE
  )

  # Filter to candidate rows using the stored .row indices
  matched <- full_part[part_candidates$.row]

  # Coerce to numeric
  matched[, (all_cols_to_load) := lapply(.SD, as.numeric),
    .SDcols = all_cols_to_load
  ]

  # Attach group id for later aggregation
  matched[, id := part_candidates$id]

  spectra_list[[i]] <- matched
  rm(full_part)
  gc(verbose = FALSE)
}

spectra <- rbindlist(spectra_list, fill = TRUE)
rm(spectra_list)
gc(verbose = FALSE)

# Add a unique point identifier for plotting
spectra[, .point_id := .I]

cat(sprintf(
  "\nExtracted %s points with full spectra (%d bands)\n\n",
  format(nrow(spectra), big.mark = ","), length(band_cols)
))

# -- 5. Export raw spectra CSV ------------------------------------------------
cat("Exporting raw spectra...\n")

raw_path <- file.path(output_dir, paste0(file_prefix, "_raw.csv"))
fwrite(spectra, raw_path)
size_mb <- round(file.size(raw_path) / 1e6, 1)
cat(sprintf("  Raw spectra: %s (%s MB)\n", basename(raw_path), size_mb))

# -- 6. Compute mean spectrum without pivot_longer ----------------------------
# pivot_longer on 1.1M points x 238 bands creates ~262M rows and blows
# past 16 GB. Instead we compute mean and SD directly in wide format using
# data.table column operations, then pivot only the small summary table
# (238 rows) into long format for plotting.
cat("Computing mean spectrum (memory-efficient wide-format approach)...\n")

# Helper: apply rolling mean smoother across band columns for each row
# Operates in-place on the data.table to avoid copying
if (smooth_window > 1) {
  cat(sprintf("  Smoothing band columns (window = %d)...\n", smooth_window))
  # Apply filter() across each row using a transposed approach:
  # work column-wise with a shifted mean rather than row-wise pivot
  for (b in seq_along(band_cols)) {
    # For each band, average with its neighbours within the window
    half_w <- floor(smooth_window / 2)
    idx_min <- max(1, b - half_w)
    idx_max <- min(length(band_cols), b + half_w)
    neighbour_cols <- band_cols[idx_min:idx_max]
    # Replace band value with mean of neighbouring bands
    spectra[, (band_cols[b]) := rowMeans(.SD, na.rm = TRUE),
      .SDcols = neighbour_cols
    ]
  }
}

# Compute mean and SD per band column -- returns a 1-row summary per group
if (plot_grouping == "by_id" && "id" %in% names(spectra)) {
  # Reshape per-group summary into long format (small: n_groups x 238 rows)
  # Note: id is provided by data.table's by = id grouping -- do not repeat it
  # inside the list() or left_join() will complain about duplicate id columns
  mean_spectra <- spectra[,
    {
      lapply(band_cols, function(b) {
        list(
          band = b,
          mean_reflectance = mean(get(b), na.rm = TRUE),
          sd_reflectance = sd(get(b), na.rm = TRUE),
          n_points = .N
        )
      }) %>% rbindlist()
    },
    by = id
  ] %>%
    left_join(wave_lookup, by = "band") %>%
    arrange(id, wavelength_nm)
} else {
  # Single mean spectrum across all points -- compute per band column
  cat(sprintf(
    "  Summarising %d band columns across %s points...\n",
    length(band_cols), format(nrow(spectra), big.mark = ",")
  ))

  mean_vals <- spectra[, lapply(.SD, mean, na.rm = TRUE), .SDcols = band_cols]
  sd_vals <- spectra[, lapply(.SD, sd, na.rm = TRUE), .SDcols = band_cols]
  n_pts <- nrow(spectra)

  # Build long-format summary from the single-row results (238 rows only)
  mean_spectra <- tibble(
    band             = band_cols,
    mean_reflectance = as.numeric(mean_vals),
    sd_reflectance   = as.numeric(sd_vals),
    n_points         = n_pts
  ) %>%
    left_join(wave_lookup, by = "band") %>%
    arrange(wavelength_nm)
}

mean_path <- file.path(output_dir, paste0(file_prefix, "_mean.csv"))
write_csv(mean_spectra, mean_path)
cat(sprintf("  Mean spectrum saved: %s\n", basename(mean_path)))

# -- 7. Plot spectral profile -------------------------------------------------
# Plot directly from the small mean_spectra summary (238 rows) --
# no need for the full long-format spectra object
cat("Plotting spectral profile...\n")

plot_path <- file.path(output_dir, paste0(file_prefix, "_plot.png"))

plot_title <- sprintf(
  "Spectral Profile -- %s (%s)", file_prefix,
  ifelse(extraction_mode == "polygon",
    basename(polygon_path),
    sprintf(
      "%d point(s), %.0fm buffer",
      nrow(extract_points), buffer_radius
    )
  )
)

# Wavelength range filter
plot_data <- mean_spectra
if (!is.null(plot_wavelength_min)) {
  plot_data <- plot_data %>% filter(wavelength_nm >= plot_wavelength_min)
}
if (!is.null(plot_wavelength_max)) {
  plot_data <- plot_data %>% filter(wavelength_nm <= plot_wavelength_max)
}

# Spectral region annotations
regions <- tribble(
  ~xmin, ~xmax, ~label, ~fill,
  400, 500, "Blue", "#AED6F1",
  500, 600, "Green", "#A9DFBF",
  600, 700, "Red", "#F1948A",
  700, 750, "Red Edge", "#D7BDE2",
  750, 900, "NIR", "#FAD7A0"
)

p <- ggplot() +
  geom_rect(
    data = regions,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = label),
    alpha = 0.12, inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = setNames(regions$fill, regions$label),
    name = "Spectral region"
  )

if (plot_grouping == "by_id" && "id" %in% names(plot_data)) {
  p <- p +
    geom_ribbon(
      data = plot_data,
      aes(
        x = wavelength_nm,
        ymin = mean_reflectance - sd_reflectance,
        ymax = mean_reflectance + sd_reflectance,
        fill = id
      ),
      alpha = 0.15, show.legend = FALSE
    ) +
    geom_line(
      data = plot_data,
      aes(x = wavelength_nm, y = mean_reflectance, colour = id),
      linewidth = 0.9
    ) +
    scale_colour_brewer(palette = "Set2", name = "Feature")
} else {
  p <- p +
    geom_ribbon(
      data = plot_data,
      aes(
        x = wavelength_nm,
        ymin = mean_reflectance - sd_reflectance,
        ymax = mean_reflectance + sd_reflectance
      ),
      fill = "steelblue", alpha = 0.25
    ) +
    geom_line(
      data = plot_data,
      aes(x = wavelength_nm, y = mean_reflectance),
      colour = "steelblue", linewidth = 1.0
    )
}

p <- p +
  labs(
    title = plot_title,
    subtitle = sprintf(
      "n = %s points | smoothing window = %d bands",
      format(nrow(spectra), big.mark = ","), smooth_window
    ),
    x = "Wavelength (nm)",
    y = "Reflectance (DN)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(plot_path, p, width = 10, height = 5, dpi = 200)
cat(sprintf("  Spectral plot saved: %s\n", basename(plot_path)))

# -- 8. Export extraction locations as geopackage ----------------------------
cat("Saving extraction point locations...\n")

pts_out <- spectra %>%
  select(.point_id, id, X, Y, Z) %>%
  st_as_sf(coords = c("X", "Y"), crs = pc_crs)

gpkg_path <- file.path(output_dir, paste0(file_prefix, "_points.gpkg"))
st_write(pts_out, gpkg_path, delete_dsn = TRUE, quiet = TRUE)
cat(sprintf("  Point locations: %s\n", basename(gpkg_path)))

cat("\n=== Spectral extraction complete ===\n")
cat(sprintf("Output directory: %s\n\n", output_dir))

# =============================================================================
# -- QUICK INTERACTIVE PLOT (run manually after extraction) --------------------
# =============================================================================
# To view the spectral profile interactively in R after running the script:
#
#   library(plotly)
#
#   # Mean spectrum with SD ribbon -- hover to see wavelength/reflectance
#   mean_spectra %>%
#     plot_ly(x = ~wavelength_nm, y = ~mean_reflectance,
#             type = "scatter", mode = "lines",
#             line = list(color = "steelblue")) %>%
#     add_ribbons(ymin = ~mean_reflectance - sd_reflectance,
#                 ymax = ~mean_reflectance + sd_reflectance,
#                 fillcolor = "rgba(70, 130, 180, 0.2)",
#                 line = list(color = "transparent")) %>%
#     layout(xaxis = list(title = "Wavelength (nm)"),
#            yaxis = list(title = "Reflectance (DN)"))
# =============================================================================

# =============================================================================
# -- SECTION 9: FOREST HEALTH METRICS -----------------------------------------
# =============================================================================
# Computes a suite of spectral indices relevant to boreal forest condition
# in Ontario, using band combinations appropriate for µCASI (407-897 nm).
#
# Indices are computed per point (appended to spectra) and summarised
# (mean +/- SD per group) in a dedicated health metrics CSV and plot.
#
# References:
#   Zarco-Tejada et al. (2001) -- chlorophyll fluorescence / red edge indices
#   Gitelson et al. (2006)     -- anthocyanin reflectance index
#   Chappelle et al. (1992)    -- ratio analysis of reflectance spectra (RARS)
#   Blackburn (1998)           -- carotenoid and chlorophyll indices
#   Ceccato et al. (2001)      -- moisture / fire vulnerability indices
#   Merton & Huntington (1999) -- RVSI red edge vegetation stress
#   Daughtry et al. (2000)     -- MCARI chlorophyll content
#   Peñuelas et al. (1994)     -- PRI photochemical reflectance (xanthophyll)
#   Gamon et al. (1992)        -- PRI original formulation
#
# Boreal context:
#   Black spruce, jack pine, trembling aspen, and white birch dominate Ontario
#   boreal stands. Key stress signatures include:
#     - Reduced red edge slope (conifer defoliation, moisture stress)
#     - Elevated anthocyanin (early senescence in aspen/birch)
#     - Reduced carotenoid:chlorophyll ratio (chronic stress)
#     - Low moisture indices (fire vulnerability in jack pine stands)
#     - PRI as a proxy for light-use efficiency (photosynthetic activity)
# =============================================================================

cat("\n=== Forest Health Metrics ===\n\n")

# -- Band lookup for all health metric wavelengths ----------------------------
# Using get_band() to find the closest available µCASI band for each target.
# Target wavelengths follow published index formulations where possible,
# adjusted to the nearest available band in the 407-897 nm µCASI range.

hb <- list(
  # Chlorophyll-sensitive bands
  b445  = get_band(445, wave_lookup), # blue -- chlorophyll a absorption
  b500  = get_band(500, wave_lookup),
  b531  = get_band(531, wave_lookup), # PRI reference band (xanthophyll)
  b550  = get_band(550, wave_lookup), # green reflectance peak
  b570  = get_band(570, wave_lookup), # PRI measurement band
  b620  = get_band(620, wave_lookup), # orange -- carotenoid absorption
  b670  = get_band(670, wave_lookup), # red -- chlorophyll a absorption
  b680  = get_band(680, wave_lookup), # red -- chlorophyll absorption peak
  b695  = get_band(695, wave_lookup), # red edge base
  b700  = get_band(700, wave_lookup), # red edge onset
  b705  = get_band(705, wave_lookup), # red edge / Chl fluorescence
  b710  = get_band(710, wave_lookup), # red edge
  b720  = get_band(720, wave_lookup), # red edge inflection
  b740  = get_band(740, wave_lookup), # red edge shoulder
  b750  = get_band(750, wave_lookup), # red edge / NIR transition
  b760  = get_band(760, wave_lookup), # NIR -- oxygen absorption region
  b780  = get_band(780, wave_lookup)[1], # NIR plateau
  b800  = get_band(800, wave_lookup), # NIR reference
  b860  = get_band(860, wave_lookup) # NIR2 -- moisture / structure
)

# Print band assignments so the user can verify closest matches
cat("Band assignments (target -> actual µCASI band):\n")
target_nms <- c(445, 500, 531, 550, 570, 620, 670, 680, 695, 700, 705, 710, 720, 740, 750, 760, 780, 800, 860)
for (j in seq_along(hb)) {
  actual_nm <- wave_lookup$wavelength_nm[wave_lookup$band == hb[[j]]]
  cat(sprintf(
    "  %3.0fnm -> %s (%.1fnm)\n",
    target_nms[j], hb[[j]], actual_nm
  ))
}
cat("\n")

# -- Compute per-point health metrics -----------------------------------------
cat("Computing per-point health indices...\n")

spectra[, `:=`(
  # ── Chlorophyll content ─────────────────────────────────────────────────

  # Chlorophyll Index Green (CIg) -- Gitelson et al. 2003
  # NIR/Green - 1; sensitive to total chlorophyll in broadleaf and conifer
  # High values = high chlorophyll; decline indicates stress or senescence
  CIg = (get(hb$b800) / get(hb$b550)) - 1,

  # Red Edge Chlorophyll Index (CIre) -- Gitelson et al. 2003
  # NIR/RedEdge - 1; less soil background influence than CIg
  # Particularly sensitive in dense boreal canopies (black spruce, balsam fir)
  CIre = (get(hb$b800) / get(hb$b720)) - 1,

  # Transformed Chlorophyll Absorption Reflectance Index (TCARI)
  # Daughtry et al. 2000 -- minimises soil and non-photosynthetic vegetation
  # TCARI > 0.2 suggests stressed or sparse canopy in boreal context
  TCARI = 3 * ((get(hb$b700) - get(hb$b670)) -
    0.2 * (get(hb$b700) - get(hb$b550)) *
      (get(hb$b700) / get(hb$b670))),

  # TCARI/OSAVI ratio -- decouples chlorophyll from canopy structure
  # More robust than TCARI alone for heterogeneous boreal stands
  TCARI_OSAVI = (3 * ((get(hb$b700) - get(hb$b670)) -
    0.2 * (get(hb$b700) - get(hb$b550)) *
      (get(hb$b700) / get(hb$b670)))) /
    ((1 + 0.16) * (get(hb$b800) - get(hb$b670)) /
      (get(hb$b800) + get(hb$b670) + 0.16)),

  # Red Edge Position proxy (REP) -- Guyot & Baret 1988 simplified
  # Uses linear interpolation between 700 and 740 nm bands
  # REP shifts toward longer wavelengths with higher chlorophyll content;
  # blueshift (toward 700) indicates chlorophyll loss / defoliation stress
  REP = 700 + 40 * ((((get(hb$b670) + get(hb$b780)) / 2) - get(hb$b700)) /
    (get(hb$b740) - get(hb$b700))),

  # ── Carotenoids ─────────────────────────────────────────────────────────

  # Carotenoid Reflectance Index 1 (CRI1) -- Gitelson et al. 2002
  # 1/Blue - 1/Green; carotenoids absorb strongly in blue
  # Increases with carotenoid:chlorophyll ratio -- early stress signal in
  # trembling aspen and white birch before visible yellowing
  CRI1 = (1 / get(hb$b445)) - (1 / get(hb$b550)),

  # Carotenoid Reflectance Index 2 (CRI2) -- Gitelson et al. 2002
  # Adds red band for improved carotenoid sensitivity in conifers
  CRI2 = (1 / get(hb$b445)) - (1 / get(hb$b700)),

  # ── Anthocyanins ────────────────────────────────────────────────────────

  # Anthocyanin Reflectance Index 1 (ARI1) -- Gitelson et al. 2001
  # 1/Green - 1/RedEdge; anthocyanins mask chlorophyll in the green band
  # Elevated values indicate stress-induced anthocyanin production --
  # common in autumn senescence, drought stress, and insect defoliation
  # in boreal broadleafs (aspen, birch)
  ARI1 = (1 / get(hb$b550)) - (1 / get(hb$b700)),

  # Anthocyanin Reflectance Index 2 (ARI2) -- Gitelson et al. 2001
  # NIR-weighted version; more robust for varying canopy densities
  ARI2 = get(hb$b800) * ((1 / get(hb$b550)) - (1 / get(hb$b700))),

  # ── Photochemical efficiency / light-use efficiency ─────────────────────

  # Photochemical Reflectance Index (PRI) -- Gamon et al. 1992
  # (531 - 570) / (531 + 570); tracks xanthophyll cycle activity
  # Negative PRI = active photoprotection (stressed / high light)
  # Positive PRI = efficient photosynthesis (healthy boreal canopy)
  # Particularly useful for jack pine and black spruce light-use monitoring
  PRI = (get(hb$b531) - get(hb$b570)) / (get(hb$b531) + get(hb$b570)),

  # ── Moisture content and fire vulnerability ──────────────────────────────

  # Normalised Difference Water Index proxy (NDWIcanopy)
  # Uses 860 nm as a moisture-sensitive NIR band (true NDWI uses SWIR at
  # 1240 nm, outside µCASI range). 780/860 ratio responds to canopy water
  # content in Ontario boreal species -- lower values = drier canopy
  # Ceccato et al. 2001 approach adapted to visible-NIR sensor range
  NDWIc = (get(hb$b780) - get(hb$b860)) / (get(hb$b780) + get(hb$b860)),

  # Moisture Stress Index (MSI) -- Hunt & Rock 1989 adapted
  # 860/780; inversely related to leaf water content
  # MSI > 1.1 suggests moisture stress in boreal conifers
  # Higher MSI = drier = higher fire vulnerability
  MSI = get(hb$b860) / get(hb$b780),

  # Plant Senescence Reflectance Index (PSRI) -- Merzlyak et al. 1999
  # (Red - Green) / RedEdge; rises sharply during senescence and drought
  # Elevated PSRI in jack pine indicates pre-fire moisture stress;
  # high PSRI across a stand = elevated fire vulnerability
  PSRI = (get(hb$b680) - get(hb$b500)) / get(hb$b750),

  # ── Canopy structure / stress integration ────────────────────────────────

  # Red Edge Vegetation Stress Index (RVSI) -- Merton & Huntington 1999
  # ((740 + 700) / 2) - 720; measures concavity of red edge
  # Positive RVSI = stressed/senescing vegetation (flattened red edge)
  # Negative RVSI = healthy actively growing vegetation (steep red edge)
  # Sensitive to defoliation in boreal conifers (spruce budworm, bark beetle)
  RVSI = ((get(hb$b740) + get(hb$b700)) / 2) - get(hb$b720),

  # Ratio Analysis of Reflectance Spectra -- chlorophyll a (RARS_a)
  # Chappelle et al. 1992: R675 / R700
  # Sensitive to chlorophyll a degradation; useful for conifer stress detection
  RARS_a = get(hb$b670) / get(hb$b700),

  # Vogelmann Red Edge Index 1 (VOG1) -- Vogelmann et al. 1993
  # R740 / R720; sensitive to chlorophyll content and red edge slope
  # Robust to illumination variation -- useful for airborne µCASI data
  VOG1 = get(hb$b740) / get(hb$b720),

  # Vogelmann Red Edge Index 2 (VOG2) -- Vogelmann et al. 1993
  # (R734 - R747) / (R715 + R726) -- uses red edge curvature
  # More sensitive to chlorophyll at high biomass (dense boreal stands)
  VOG2 = (get(hb$b740) - get(hb$b750)) / (get(hb$b710) + get(hb$b720))
)]

# Fix b500 reference (used in PSRI) -- not in hb list, add separately
b500 <- get_band(500, wave_lookup)
spectra[, PSRI := (get(hb$b680) - get(b500)) / get(hb$b750)]

health_metric_cols <- c(
  "CIg", "CIre", "TCARI", "TCARI_OSAVI", # chlorophyll
  "REP", # red edge position
  "CRI1", "CRI2", # carotenoids
  "ARI1", "ARI2", # anthocyanins
  "PRI", # photochemical efficiency
  "NDWIc", "MSI", "PSRI", # moisture / fire vulnerability
  "RVSI", "RARS_a", "VOG1", "VOG2" # stress / structure
)

cat(sprintf("  Computed %d health indices per point\n\n", length(health_metric_cols)))

# -- Export per-point health metrics ------------------------------------------
cat("Exporting per-point health metrics...\n")

health_pts_path <- file.path(output_dir, paste0(file_prefix, "_health_points.csv"))
fwrite(
  spectra[, c(".point_id", "id", "X", "Y", "Z", health_metric_cols),
    with = FALSE
  ],
  health_pts_path
)
size_mb <- round(file.size(health_pts_path) / 1e6, 1)
cat(sprintf(
  "  Per-point metrics: %s (%s MB)\n",
  basename(health_pts_path), size_mb
))

# -- Summary table: mean +/- SD per group or overall --------------------------
cat("Summarising health metrics...\n")

if (plot_grouping == "by_id" && "id" %in% names(spectra)) {
  health_summary <- spectra[, lapply(.SD, function(x) {
    list(
      mean = round(mean(x, na.rm = TRUE), 4),
      sd = round(sd(x, na.rm = TRUE), 4)
    )
  }),
  .SDcols = health_metric_cols, by = id
  ]

  # Reshape to tidy long format: id, metric, mean, sd
  health_long <- spectra[,
    {
      lapply(health_metric_cols, function(m) {
        list(
          metric = m,
          mean = round(mean(get(m), na.rm = TRUE), 4),
          sd = round(sd(get(m), na.rm = TRUE), 4),
          n = sum(!is.na(get(m)))
        )
      }) %>% rbindlist()
    },
    by = id
  ] %>%
    arrange(id, metric)
} else {
  health_long <- data.table(
    metric = health_metric_cols,
    mean = round(sapply(health_metric_cols, function(m) {
      mean(spectra[[m]], na.rm = TRUE)
    }), 4),
    sd = round(sapply(health_metric_cols, function(m) {
      sd(spectra[[m]], na.rm = TRUE)
    }), 4),
    n = sapply(health_metric_cols, function(m) {
      sum(!is.na(spectra[[m]]))
    })
  ) %>% arrange(metric)
}

health_summary_path <- file.path(
  output_dir,
  paste0(file_prefix, "_health_summary.csv")
)
fwrite(health_long, health_summary_path)
cat(sprintf("  Summary table: %s\n\n", basename(health_summary_path)))

# -- Forest health dashboard plot ---------------------------------------------
cat("Plotting forest health dashboard...\n")

# Interpretation thresholds for healthy boreal forest
# Based on published ranges for Ontario boreal species
# (black spruce, jack pine, trembling aspen, white birch)
# Sources: Zarco-Tejada et al. 2001, Gitelson et al. 2006, Peñuelas et al. 1994
health_thresholds <- tribble(
  ~metric, ~healthy_min, ~healthy_max, ~label, ~category,
  "CIg", 2.0, 8.0, "Chl Index Green", "Chlorophyll",
  "CIre", 1.5, 5.0, "Chl Index Red Edge", "Chlorophyll",
  "TCARI", 0.0, 0.2, "TCARI (low=healthy)", "Chlorophyll",
  "TCARI_OSAVI", 0.0, 0.4, "TCARI/OSAVI (low=healthy)", "Chlorophyll",
  "REP", 715.0, 730.0, "Red Edge Position (nm)", "Chlorophyll",
  "CRI1", 0.0, 0.03, "Carotenoid Index 1", "Carotenoids",
  "CRI2", 0.0, 0.05, "Carotenoid Index 2", "Carotenoids",
  "ARI1", 0.0, 0.02, "Anthocyanin Index 1", "Anthocyanins",
  "ARI2", 0.0, 0.10, "Anthocyanin Index 2", "Anthocyanins",
  "PRI", 0.0, 0.05, "Photochem. Reflectance", "Photochemistry",
  "NDWIc", 0.0, 0.15, "Canopy Water Index", "Moisture/Fire",
  "MSI", 0.8, 1.1, "Moisture Stress (low=wet)", "Moisture/Fire",
  "PSRI", -0.1, 0.05, "Senescence Index (low=green)", "Moisture/Fire",
  "RVSI", -0.5, 0.0, "Red Edge Stress (neg=healthy)", "Stress",
  "RARS_a", 0.3, 0.6, "Chl-a Ratio", "Stress",
  "VOG1", 2.0, 6.0, "Vogelmann Index 1", "Stress",
  "VOG2", -0.05, 0.0, "Vogelmann Index 2", "Stress"
)

# Join computed means with thresholds for plotting
if (plot_grouping == "by_id" && "id" %in% names(health_long)) {
  plot_health <- health_long %>%
    left_join(health_thresholds, by = "metric")
} else {
  plot_health <- health_long %>%
    left_join(health_thresholds, by = "metric") %>%
    mutate(id = "All points")
}

# Flag each metric as healthy, stressed, or outside published range
plot_health <- plot_health %>%
  mutate(
    status = case_when(
      mean >= healthy_min & mean <= healthy_max ~ "Healthy",
      TRUE ~ "Outside healthy range"
    ),
    # Normalise mean to 0-1 within the healthy range for a common axis
    mean_norm = (mean - healthy_min) / (healthy_max - healthy_min)
  )

# Dashboard: one panel per category, bars coloured by health status
p_health <- ggplot(
  plot_health,
  aes(
    x = reorder(label, mean_norm),
    y = mean,
    fill = status
  )
) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.25, colour = "grey30", linewidth = 0.4
  ) +
  # Shade the healthy reference range
  geom_rect(
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = healthy_min, ymax = healthy_max
    ),
    fill = "green", alpha = 0.06, inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Healthy" = "#2E7D32",
      "Outside healthy range" = "#C62828"
    ),
    name = "Status"
  ) +
  facet_wrap(~category, scales = "free", ncol = 2) +
  coord_flip() +
  labs(
    title = sprintf("Forest Health Dashboard -- %s", file_prefix),
    subtitle = sprintf(
      "n = %s points | Healthy range = published boreal thresholds",
      format(nrow(spectra), big.mark = ",")
    ),
    x = NULL,
    y = "Index value"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

health_plot_path <- file.path(
  output_dir,
  paste0(file_prefix, "_health_dashboard.png")
)
ggsave(health_plot_path, p_health, width = 14, height = 10, dpi = 200)
cat(sprintf("  Health dashboard: %s\n\n", basename(health_plot_path)))

# -- Print console summary ----------------------------------------------------
cat("── Health Metric Summary ────────────────────────────────────────────────\n")
cat(sprintf("  %-20s  %8s  %8s  %s\n", "Metric", "Mean", "SD", "Status"))
cat(sprintf("  %-20s  %8s  %8s  %s\n", "------", "----", "--", "------"))

summary_print <- if ("id" %in% names(health_long)) {
  health_long %>%
    group_by(metric) %>%
    summarise(
      mean = mean(mean, na.rm = TRUE),
      sd = mean(sd, na.rm = TRUE), .groups = "drop"
    )
} else {
  health_long
}

summary_print %>%
  left_join(health_thresholds %>% select(metric, healthy_min, healthy_max, label),
    by = "metric"
  ) %>%
  mutate(status = ifelse(mean >= healthy_min & mean <= healthy_max,
    "OK", "CHECK"
  )) %>%
  arrange(category) %>%
  {
    for (i in seq_len(nrow(.))) {
      cat(sprintf(
        "  %-20s  %8.4f  %8.4f  %s\n",
        .$metric[i], .$mean[i], .$sd[i], .$status[i]
      ))
    }
    .
  } %>%
  invisible()

cat("─────────────────────────────────────────────────────────────────────────\n")
cat("\nNote: 'CHECK' flags values outside published boreal reference ranges.\n")
cat("Thresholds are indicative -- interpret in context of stand type,\n")
cat("phenological stage, and local site conditions.\n\n")

cat("=== Forest health metrics complete ===\n")
cat(sprintf("Output directory: %s\n\n", output_dir))
