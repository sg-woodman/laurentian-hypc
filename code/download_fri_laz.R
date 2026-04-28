library(here)
library(tidyverse)
library(sf)


# ── 1. Load tile index and AOI ────────────────────────────────────────────────

tile_index <- st_read(
  here("data/raw/FRI_Leaf_On_Tile_Index_GeoPackage.gpkg"),
  layer = "FRI_Tile_Index"
)

aoi <- st_read(here("data/processed/uCASI_flight_union.gpkg")) |>
  filter(site == "daisy")

# Ensure matching CRS
tile_index <- st_transform(tile_index, st_crs(aoi))

# ── 2. Intersect AOI with tile index ─────────────────────────────────────────

tiles_needed <- st_filter(tile_index, aoi, .predicate = st_intersects)
cat("Tiles to download:", nrow(tiles_needed), "\n")
print(tiles_needed$Tilename)

# ── 3. Batch download a given product column ──────────────────────────────────

download_tiles <- function(tiles, url_col, download_dir) {
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

  urls <- tiles[[url_col]]
  filenames <- basename(urls)
  destfiles <- file.path(download_dir, filenames)
  n <- length(urls)

  for (i in seq_along(urls)) {
    if (file.exists(destfiles[i])) {
      message(sprintf("[%d/%d] Skipping (exists): %s", i, n, filenames[i]))
      next
    }

    message(sprintf("[%d/%d] Downloading: %s", i, n, filenames[i]))
    tryCatch(
      download.file(urls[i], destfiles[i], mode = "wb", quiet = TRUE),
      error = \(e) {
        message(sprintf(
          "[%d/%d] FAILED: %s — %s",
          i,
          n,
          filenames[i],
          conditionMessage(e)
        ))
      }
    )
  }
}

# ── 4. Download HAG LAZ and/or CHM tiles ─────────────────────────────────────
download_dir <- "data/raw/fri_laz"
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

download_tiles(tiles_needed, "Download_HAG", "data/raw/fri_hag")
download_tiles(tiles_needed, "Download_CHM", "data/fri_chm")
