# =============================================================================
# Daisy North — Tree Crown Segmentation & Validation
# =============================================================================
# Workflow:
#   1. Filter LAS catalog to Daisy North tiles
#   2. Segment trees using li2012 (point cloud, no raster CHM required)
#   3. Extract crown polygons
#   4. Validate against field plot data (daisy_plot sf object)
# =============================================================================

library(lidR)
library(sf)
library(dplyr)
library(purrr)
library(future)
library(stringr)
library(ggplot2)

# ── 0. Configuration ──────────────────────────────────────────────────────────

LAZ_DIR <- "data/fri_hag" # directory of downloaded HAG LAZ tiles
OUTPUT_DIR <- "data/fri_seg" # output directory for segmented tiles
CROWNS_OUT <- "data/crowns_daisy_north.gpkg" # merged crown polygons output
MATCH_DIST_M <- 3 # max metres for field tree → crown match
N_WORKERS <- 40 # parallel workers for HPC

# li2012 parameters (based on tuning test — adjust as needed)
SEG_PARAMS <- li2012(dt1 = 1.0, dt2 = 1.5, R = 8, Zu = 10, hmin = 3)

# Daisy North tile names
DAISY_NORTH_SITES <- c(
  "1kmZ175070514302022L",
  "1kmZ175080514402022L",
  "1kmZ175070514402022L",
  "1kmZ175090514502022L",
  "1kmZ175080514502022L",
  "1kmZ175090514402022L"
)

# =============================================================================
# STAGE 1: Filter LAS catalog to Daisy North tiles
# =============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Build full catalog then filter to target tiles
all_laz_files <- list.files(LAZ_DIR, pattern = "\\.laz$", full.names = TRUE)

filtered_files <- keep(
  all_laz_files,
  ~ str_detect(.x, paste(DAISY_NORTH_SITES, collapse = "|"))
)

cat(sprintf(
  "Tiles matched: %d of %d\n",
  length(filtered_files),
  length(all_laz_files)
))
stopifnot(
  "No matching tiles found — check LAZ_DIR and tile names" = length(
    filtered_files
  ) >
    0
)

# Build catalog from filtered tiles only
ctg <- readLAScatalog(filtered_files)

cat("Catalog summary:\n")
print(ctg)

# =============================================================================
# STAGE 2: Tree segmentation via LAS catalog engine
# =============================================================================

# Catalog processing options
opt_filter(ctg) <- "-keep_class 4 5" # medium + high vegetation only
opt_chunk_buffer(ctg) <- 20 # 20m buffer for edge handling
opt_chunk_size(ctg) <- 0 # use original tile size
opt_output_files(ctg) <- file.path(OUTPUT_DIR, "{ORIGINALFILENAME}_seg")
opt_progress(ctg) <- TRUE # per-tile timing output

# Parallel processing
plan(multisession, workers = N_WORKERS)

cat(sprintf(
  "\nRunning tree segmentation on %d tiles with %d workers...\n",
  length(filtered_files),
  N_WORKERS
))

seg_time <- system.time(
  ctg_seg <- segment_trees(ctg, SEG_PARAMS)
)

plan(sequential)

cat(sprintf(
  "\nSegmentation complete in %.1f minutes\n",
  seg_time["elapsed"] / 60
))

# =============================================================================
# STAGE 3: Extract and merge crown polygons
# =============================================================================

cat("\nExtracting crown polygons...\n")

seg_files <- list.files(OUTPUT_DIR, pattern = "_seg\\.laz$", full.names = TRUE)
n <- length(seg_files)

crowns_list <- vector("list", n)

for (i in seq_along(seg_files)) {
  message(sprintf(
    "[%d/%d] Extracting crowns: %s",
    i,
    n,
    basename(seg_files[i])
  ))

  tryCatch(
    {
      las_seg <- readLAS(seg_files[i])

      crowns_list[[i]] <- crown_metrics(
        las = las_seg,
        func = ~ list(
          Z_max = max(Z),
          Z_mean = mean(Z),
          Z_sd = sd(Z),
          npoints = .N
        ),
        geom = "convex"
      )
    },
    error = \(e) message(sprintf("  FAILED: %s", conditionMessage(e)))
  )
}

# Merge and assign globally unique treeIDs
crowns_all <- bind_rows(crowns_list) |>
  mutate(treeID = row_number())

cat(sprintf("\nTotal crowns segmented: %d\n", nrow(crowns_all)))

# Save crown polygons
st_write(crowns_all, CROWNS_OUT, delete_if_exists = TRUE)
cat(sprintf("Crown polygons saved to: %s\n", CROWNS_OUT))

# =============================================================================
# STAGE 4: Validation against field plot data
# =============================================================================

cat("\nRunning validation against field plot data...\n")

# ── 4a. Align CRS ─────────────────────────────────────────────────────────────

# daisy_plot is assumed to be loaded in the environment as an sf object
stopifnot("daisy_plot not found — load field data first" = exists("daisy_plot"))

if (st_crs(daisy_plot) != st_crs(crowns_all)) {
  daisy_plot <- st_transform(daisy_plot, st_crs(crowns_all))
}

# ── 4b. Match each field tree to nearest crown centroid ───────────────────────

nearest <- st_nearest_feature(daisy_plot, crowns_all)
distances <- st_distance(daisy_plot, crowns_all[nearest, ], by_element = TRUE)

trees_matched <- daisy_plot |>
  mutate(
    matched_treeID = crowns_all$treeID[nearest],
    match_dist_m = as.numeric(distances),
    matched = match_dist_m <= MATCH_DIST_M,
    crown_class = factor(
      crown_class,
      levels = c("Dominant", "Co-dominant", "Intermediate", "Suppressed")
    )
  )

# ── 4c. Overall detection metrics ─────────────────────────────────────────────

plot_extent <- st_convex_hull(st_union(daisy_plot)) |> st_buffer(15)
crowns_in_aoi <- crowns_all[
  st_intersects(crowns_all, plot_extent, sparse = FALSE),
]

n_field <- nrow(trees_matched)
n_detected <- sum(trees_matched$matched)
n_crowns <- nrow(crowns_in_aoi)
n_false_pos <- max(0, n_crowns - n_detected)

precision <- n_detected / n_crowns
recall <- n_detected / n_field

overall_metrics <- tibble(
  n_field_trees = n_field,
  n_detected = n_detected,
  n_crowns_lidar = n_crowns,
  recall = round(recall, 3),
  precision = round(precision, 3),
  f_score = round(2 * precision * recall / (precision + recall), 3),
  omission_rate = round(1 - recall, 3),
  commission_rate = round(n_false_pos / n_crowns, 3)
)

cat(
  "\n── Overall Detection Metrics ────────────────────────────────────────────\n"
)
print(overall_metrics)

# ── 4d. Detection by crown class ──────────────────────────────────────────────

detection_by_class <- trees_matched |>
  st_drop_geometry() |>
  group_by(crown_class) |>
  summarise(
    n_field = n(),
    n_detected = sum(matched),
    recall = round(n_detected / n_field, 3),
    omission = round(1 - n_detected / n_field, 3),
    mean_match_dist_m = round(mean(match_dist_m[matched], na.rm = TRUE), 2),
    .groups = "drop"
  ) |>
  arrange(crown_class)

cat(
  "\n── Detection by Crown Class ─────────────────────────────────────────────\n"
)
print(detection_by_class)

# ── 4e. Detection by crown class and species ──────────────────────────────────

detection_by_class_species <- trees_matched |>
  st_drop_geometry() |>
  group_by(crown_class, species) |>
  summarise(
    n_field = n(),
    n_detected = sum(matched),
    recall = round(n_detected / n_field, 3),
    .groups = "drop"
  ) |>
  filter(n_field >= 3) |>
  arrange(crown_class, desc(recall))

cat(
  "\n── Detection by Crown Class and Species ─────────────────────────────────\n"
)
print(detection_by_class_species, n = Inf)

# ── 4f. Detection by plot ──────────────────────────────────────────────────────

detection_by_plot <- trees_matched |>
  st_drop_geometry() |>
  group_by(plot) |>
  summarise(
    n_field = n(),
    n_detected = sum(matched),
    recall = round(n_detected / n_field, 3),
    .groups = "drop"
  ) |>
  arrange(recall)

cat(
  "\n── Detection by Plot ────────────────────────────────────────────────────\n"
)
print(detection_by_plot)

# ── 4g. DBH comparison: matched vs missed trees ────────────────────────────────

dbh_comparison <- trees_matched |>
  st_drop_geometry() |>
  group_by(matched) |>
  summarise(
    n = n(),
    mean_dbh = round(mean(dbh_cm, na.rm = TRUE), 1),
    median_dbh = round(median(dbh_cm, na.rm = TRUE), 1),
    sd_dbh = round(sd(dbh_cm, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  mutate(matched = if_else(matched, "Detected", "Missed"))

cat(
  "\n── DBH Comparison: Detected vs Missed Trees ─────────────────────────────\n"
)
print(dbh_comparison)

# =============================================================================
# STAGE 5: Validation plots
# =============================================================================

cat("\nGenerating validation plots...\n")

# ── 5a. Detection rate by crown class (bar chart) ─────────────────────────────

p_class <- ggplot(
  detection_by_class,
  aes(x = crown_class, y = recall, fill = crown_class)
) +
  geom_col(show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.0f%%\n(n=%d)", recall * 100, n_field)),
    vjust = -0.3,
    size = 3.5
  ) +
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.8) +
  scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
  labs(
    title = "Tree Detection Rate by Crown Class",
    subtitle = sprintf(
      "Overall recall = %.1f%%,  F-score = %.3f",
      recall * 100,
      overall_metrics$f_score
    ),
    x = "Crown Class",
    y = "Recall (Detection Rate)"
  ) +
  theme_minimal(base_size = 12)

# ── 5b. DBH distribution: detected vs missed ─────────────────────────────────

p_dbh <- trees_matched |>
  st_drop_geometry() |>
  mutate(status = if_else(matched, "Detected", "Missed")) |>
  ggplot(aes(x = dbh_cm, fill = status)) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("Detected" = "#2196F3", "Missed" = "#F44336")) +
  labs(
    title = "DBH Distribution: Detected vs Missed Trees",
    x = "DBH (cm)",
    y = "Count",
    fill = NULL
  ) +
  theme_minimal(base_size = 12)

# ── 5c. Spatial map: matched field trees over crown polygons ──────────────────

p_spatial <- ggplot() +
  geom_sf(
    data = crowns_in_aoi,
    fill = "lightblue",
    colour = "steelblue",
    linewidth = 0.2,
    alpha = 0.5
  ) +
  geom_sf(
    data = trees_matched,
    aes(colour = matched, shape = crown_class),
    size = 1.8
  ) +
  scale_colour_manual(
    values = c("TRUE" = "#2196F3", "FALSE" = "#F44336"),
    labels = c("TRUE" = "Detected", "FALSE" = "Missed")
  ) +
  labs(
    title = "Field Trees vs Segmented Crowns",
    colour = "Match Status",
    shape = "Crown Class"
  ) +
  theme_minimal(base_size = 11)

print(p_class)
print(p_dbh)
print(p_spatial)

cat("\nDone.\n")
