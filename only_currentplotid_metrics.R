library(terra)
library(sf)
library(landscapemetrics)
library(dplyr)
library(data.table)

# ====================================================
# PATHS
# ====================================================

scenario_name <- "holistic_hotspot"
suffix        <- "hh"

raster_dir   <- "R:/Chapter_3_fragmentation/2026_NEE_R2/binary_forest_scenarios/holistic_hotspot"
fishnet_dir  <- "R:/Chapter_3_fragmentation/2026_NEE_R2/fishnets_projected"
current_dir  <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/frag_metrics_current_10m"
output_dir   <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FutureScenario_metrics/holistic_hotspot_binary_metrics"
lookup_csv   <- "R:/Chapter_3_fragmentation/2026_NEE_R2/gfw_pnr_lookup.csv"
gfw_folder   <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/current_forest_10m"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ====================================================
# IDENTIFY NO-PNR-OVERLAP TILES — same logic as Python script
# These are GFW tiles on disk that never appear in the lookup
# for this scenario. They were direct-copied from GFW (CopyRaster)
# rather than overlaid, so their "scenario" raster = current raster
# exactly. No new forest = no FFI change possible, so they add
# nothing to the comparison and can be skipped here.
# ====================================================

lookup <- fread(lookup_csv)
pnr_tiles_hh <- unique(lookup[scenario == "holistic_hotspot", gfw_name])

all_gfw_files <- list.files(gfw_folder, pattern = "\\.tif$", full.names = FALSE)
all_gfw_names <- gsub("_10m\\.tif$", "", all_gfw_files)

no_pnr_tiles <- setdiff(all_gfw_names, pnr_tiles_hh)

cat("Total GFW tiles on disk:        ", length(all_gfw_names), "\n")
cat("Tiles WITH PNR overlap (hh):    ", length(pnr_tiles_hh), "\n")
cat("Tiles with NO PNR overlap (hh): ", length(no_pnr_tiles), "\n")
cat("No-PNR tile names:\n")
print(no_pnr_tiles)
cat("\n")

# ====================================================
# FILE LISTS — all tiles (no trial restriction)
# ====================================================

raster_files <- list.files(
  raster_dir,
  pattern = paste0("_", suffix, "_b\\.tif$"),
  full.names = TRUE
)

tile_names_all <- gsub(paste0("_", suffix, "_b\\.tif$"), "", basename(raster_files))

cat("Total raster tiles found:", length(tile_names_all), "\n")

# Remove no-PNR tiles — their scenario raster is identical to current
# forest (direct copy), so recalculating metrics would just reproduce
# the current forest metrics exactly. No point spending compute time.
tile_names_all <- setdiff(tile_names_all, no_pnr_tiles)
raster_files   <- raster_files[
  gsub(paste0("_", suffix, "_b\\.tif$"), "", basename(raster_files)) %in% tile_names_all
]

cat("Tiles after removing no-PNR overlap:", length(tile_names_all), "\n")

# ====================================================
# CHECKPOINTING — skip already completed tiles
# ====================================================

already_done <- gsub(
  paste0("_", suffix, "_b\\.csv$"), "",
  list.files(output_dir, pattern = paste0("_", suffix, "_b\\.csv$"))
)

tiles_to_run <- setdiff(tile_names_all, already_done)

cat("Already completed:    ", length(already_done), "\n")
cat("Remaining to process: ", length(tiles_to_run), "\n\n")

raster_files <- raster_files[tile_names_all %in% tiles_to_run]

overall_start <- Sys.time()

# ====================================================
# TILE LOOP
# ====================================================

for (r_file in raster_files) {
  
  tile_start  <- Sys.time()
  tile_name   <- gsub(paste0("_", suffix, "_b\\.tif$"), "", basename(r_file))
  output_file <- file.path(output_dir, paste0(tile_name, "_", suffix, "_b.csv"))
  
  cat("\n====================================================\n")
  cat("Tile:", tile_name, "| Scenario:", scenario_name, "\n")
  cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("====================================================\n")
  
  # ------------------------------------------------
  # LOAD CURRENT FOREST METRICS — get valid plot_ids
  # Only calculate scenario metrics for plots that
  # already have forest in the current baseline.
  # Plots with new PNR forest on previously non-forest
  # land have no FFI reference and are excluded here.
  # ------------------------------------------------
  
  current_csv <- file.path(current_dir, paste0(tile_name, "_metrics.csv"))
  
  if (!file.exists(current_csv)) {
    cat(" -> No current forest metrics found for", tile_name, "— skipping\n")
    next
  }
  
  current_metrics <- fread(current_csv)
  valid_plot_ids  <- unique(current_metrics[class == 1, plot_id])
  
  cat(" -> Current forest metrics loaded |",
      length(valid_plot_ids), "plots with existing forest\n")
  
  rm(current_metrics); gc()
  
  if (length(valid_plot_ids) == 0) {
    cat(" -> No forested plots in current metrics for this tile — skipping\n")
    next
  }
  
  # ------------------------------------------------
  # MATCH FISHNET
  # ------------------------------------------------
  
  fishnet_path <- file.path(fishnet_dir, paste0(tile_name, "_fishnet_projected.gpkg"))
  
  if (!file.exists(fishnet_path)) {
    cat(" -> No matching fishnet found, skipping\n")
    next
  }
  
  # ------------------------------------------------
  # LOAD RASTER
  # ------------------------------------------------
  
  r <- tryCatch(rast(r_file), error = function(e) {
    cat(" -> ERROR reading raster:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(r)) next
  cat(" -> Raster loaded |", nrow(r), "x", ncol(r),
      "| Cells:", format(ncell(r), big.mark = ","), "\n")
  
  # ------------------------------------------------
  # LOAD, ALIGN, AND FILTER FISHNET
  # Only keep cells that had forest in current metrics
  # ------------------------------------------------
  
  fishnet_sf <- st_read(fishnet_path, quiet = TRUE)
  fishnet_sf <- st_transform(fishnet_sf, crs = crs(r))
  
  n_before   <- nrow(fishnet_sf)
  fishnet_sf <- fishnet_sf[fishnet_sf$plot_id %in% valid_plot_ids, ]
  n_after    <- nrow(fishnet_sf)
  
  cat(" -> Fishnet filtered:", n_before, "->", n_after,
      "cells (skipped", n_before - n_after, "non-forest plots)\n")
  
  if (n_after == 0) {
    cat(" -> No matching forested plots in fishnet — skipping\n")
    rm(r, fishnet_sf); gc()
    next
  }
  
  fishnet_vect <- vect(fishnet_sf)
  has_row_col  <- all(c("row", "col") %in% names(fishnet_sf))
  
  # ====================================================
  # CELL LOOP
  # ====================================================
  
  metrics_list <- list()
  empty_count  <- 0
  valid_count  <- 0
  
  for (i in seq_len(nrow(fishnet_sf))) {
    
    gc()
    
    if (i %% 500 == 0) {
      elapsed <- round(as.numeric(difftime(Sys.time(), tile_start, units = "mins")), 1)
      cat("  Progress:", i, "/", nrow(fishnet_sf), "| Elapsed:", elapsed, "mins\n")
    }
    
    poly_sf   <- fishnet_sf[i, ]
    poly_vect <- fishnet_vect[i, ]
    
    r_crop <- tryCatch(crop(r, poly_vect),  error = function(e) NULL)
    if (is.null(r_crop)) { rm(poly_sf, poly_vect); gc(); next }
    
    r_mask <- tryCatch(mask(r_crop, poly_vect), error = function(e) NULL)
    if (is.null(r_mask)) { rm(poly_sf, poly_vect, r_crop); gc(); next }
    
    valid_cells <- global(!is.na(r_mask), "sum")[1, 1]
    if (is.na(valid_cells) || valid_cells == 0) {
      empty_count <- empty_count + 1
      rm(poly_sf, poly_vect, r_crop, r_mask); gc()
      next
    }
    
    valid_count <- valid_count + 1
    
    metrics <- tryCatch({
      calculate_lsm(
        r_mask,
        what = c("lsm_c_ed",
                 "lsm_c_pd",
                 "lsm_c_area_mn",
                 "lsm_c_pland")
      )
    }, error = function(e) NULL)
    
    if (is.null(metrics)) {
      rm(poly_sf, poly_vect, r_crop, r_mask); gc()
      next
    }
    
    centroid <- st_coordinates(st_centroid(poly_sf))
    
    metrics$plot_id  <- poly_sf$plot_id
    metrics$center_x <- centroid[1]
    metrics$center_y <- centroid[2]
    metrics$tile_id  <- tile_name
    metrics$scenario <- scenario_name
    
    if (has_row_col) {
      metrics$row <- poly_sf$row
      metrics$col <- poly_sf$col
    }
    
    metrics_list[[length(metrics_list) + 1]] <- metrics
    
    rm(poly_sf, poly_vect, r_crop, r_mask, metrics, centroid)
    gc()
  }
  
  # ====================================================
  # SAVE
  # ====================================================
  
  if (length(metrics_list) == 0) {
    cat(" -> No valid metrics for this tile, skipping save\n")
    rm(r, fishnet_sf, fishnet_vect, metrics_list); gc()
    next
  }
  
  final_metrics <- bind_rows(metrics_list)
  
  col_order <- c("plot_id", "center_x", "center_y",
                 "layer", "level", "class", "id",
                 "metric", "value", "tile_id", "scenario")
  if (has_row_col) col_order <- c(col_order, "row", "col")
  
  final_metrics <- final_metrics %>% select(any_of(col_order))
  
  write.csv(final_metrics, output_file, row.names = FALSE)
  
  cat("\nValid grids:", valid_count,
      "| Empty grids:", empty_count,
      "| Rows saved:", nrow(final_metrics), "\n")
  cat("Tile time:",
      round(as.numeric(difftime(Sys.time(), tile_start, units = "mins")), 2),
      "mins\n")
  
  rm(r, fishnet_sf, fishnet_vect, metrics_list, final_metrics)
  gc()
}

overall_time <- round(
  as.numeric(difftime(Sys.time(), overall_start, units = "hours")), 2
)
done_count <- length(list.files(
  output_dir, pattern = paste0("_", suffix, "_b\\.csv$")
))

cat("\n====================================================\n")
cat("ALL HOLISTIC_HOTSPOT TILES COMPLETE\n")
cat("Total runtime:  ", overall_time, "hours\n")
cat("CSVs in output: ", done_count, "\n")
cat("====================================================\n")


#############################

library(terra)
library(sf)
library(landscapemetrics)
library(dplyr)
library(data.table)

# ====================================================
# PATHS
# ====================================================

scenario_name <- "all_pnr"
suffix        <- "ap"

raster_dir   <- "R:/Chapter_3_fragmentation/2026_NEE_R2/binary_forest_scenarios/all_pnr"
fishnet_dir  <- "R:/Chapter_3_fragmentation/2026_NEE_R2/fishnets_projected"
current_dir  <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/frag_metrics_current_10m"
output_dir   <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FutureScenario_metrics/all_pnr_metrics"
lookup_csv   <- "R:/Chapter_3_fragmentation/2026_NEE_R2/gfw_pnr_lookup.csv"
gfw_folder   <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/current_forest_10m"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ====================================================
# IDENTIFY NO-PNR-OVERLAP TILES
# ====================================================

lookup <- fread(lookup_csv)
pnr_tiles_ap <- unique(lookup[scenario == "all_pnr", gfw_name])

all_gfw_files <- list.files(gfw_folder, pattern = "\\.tif$", full.names = FALSE)
all_gfw_names <- gsub("_10m\\.tif$", "", all_gfw_files)

no_pnr_tiles <- setdiff(all_gfw_names, pnr_tiles_ap)

cat("Total GFW tiles on disk:         ", length(all_gfw_names), "\n")
cat("Tiles WITH PNR overlap (all_pnr):", length(pnr_tiles_ap), "\n")
cat("Tiles with NO PNR overlap:       ", length(no_pnr_tiles), "\n")
cat("No-PNR tile names:\n")
print(no_pnr_tiles)
cat("\n")

# ====================================================
# FILE LISTS
# ====================================================

raster_files <- list.files(
  raster_dir,
  pattern = paste0("_", suffix, "_b\\.tif$"),
  full.names = TRUE
)

tile_names_all <- gsub(paste0("_", suffix, "_b\\.tif$"), "", basename(raster_files))

cat("Total raster tiles found:", length(tile_names_all), "\n")

tile_names_all <- setdiff(tile_names_all, no_pnr_tiles)
raster_files   <- raster_files[
  gsub(paste0("_", suffix, "_b\\.tif$"), "", basename(raster_files)) %in% tile_names_all
]

cat("Tiles after removing no-PNR overlap:", length(tile_names_all), "\n")

# ====================================================
# CHECKPOINTING
# ====================================================

already_done <- gsub(
  paste0("_", suffix, "_b\\.csv$"), "",
  list.files(output_dir, pattern = paste0("_", suffix, "_b\\.csv$"))
)

tiles_to_run <- setdiff(tile_names_all, already_done)

cat("Already completed:    ", length(already_done), "\n")
cat("Remaining to process: ", length(tiles_to_run), "\n\n")

raster_files <- raster_files[tile_names_all %in% tiles_to_run]

overall_start <- Sys.time()

# ====================================================
# TILE LOOP
# ====================================================

for (r_file in raster_files) {
  
  tile_start  <- Sys.time()
  tile_name   <- gsub(paste0("_", suffix, "_b\\.tif$"), "", basename(r_file))
  output_file <- file.path(output_dir, paste0(tile_name, "_", suffix, "_b.csv"))
  
  cat("\n====================================================\n")
  cat("Tile:", tile_name, "| Scenario:", scenario_name, "\n")
  cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("====================================================\n")
  
  # Load current forest metrics to get valid plot_ids
  current_csv <- file.path(current_dir, paste0(tile_name, "_metrics.csv"))
  
  if (!file.exists(current_csv)) {
    cat(" -> No current forest metrics found for", tile_name, "skipping\n")
    next
  }
  
  current_metrics <- fread(current_csv)
  valid_plot_ids  <- unique(current_metrics[class == 1, plot_id])
  
  cat(" -> Current forest metrics loaded |",
      length(valid_plot_ids), "plots with existing forest\n")
  
  rm(current_metrics); gc()
  
  if (length(valid_plot_ids) == 0) {
    cat(" -> No forested plots in current metrics for this tile skipping\n")
    next
  }
  
  # Match fishnet
  fishnet_path <- file.path(fishnet_dir,
                            paste0(tile_name, "_fishnet_projected.gpkg"))
  
  if (!file.exists(fishnet_path)) {
    cat(" -> No matching fishnet found, skipping\n")
    next
  }
  
  # Load raster
  r <- tryCatch(rast(r_file), error = function(e) {
    cat(" -> ERROR reading raster:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(r)) next
  cat(" -> Raster loaded |", nrow(r), "x", ncol(r),
      "| Cells:", format(ncell(r), big.mark = ","), "\n")
  
  # Load, align, and filter fishnet
  fishnet_sf <- st_read(fishnet_path, quiet = TRUE)
  fishnet_sf <- st_transform(fishnet_sf, crs = crs(r))
  
  n_before   <- nrow(fishnet_sf)
  fishnet_sf <- fishnet_sf[fishnet_sf$plot_id %in% valid_plot_ids, ]
  n_after    <- nrow(fishnet_sf)
  
  cat(" -> Fishnet filtered:", n_before, "->", n_after,
      "cells (skipped", n_before - n_after, "non-forest plots)\n")
  
  if (n_after == 0) {
    cat(" -> No matching forested plots in fishnet skipping\n")
    rm(r, fishnet_sf); gc()
    next
  }
  
  fishnet_vect <- vect(fishnet_sf)
  has_row_col  <- all(c("row", "col") %in% names(fishnet_sf))
  
  # Cell loop
  metrics_list <- list()
  empty_count  <- 0
  valid_count  <- 0
  
  for (i in seq_len(nrow(fishnet_sf))) {
    
    gc()
    
    if (i %% 500 == 0) {
      elapsed <- round(as.numeric(difftime(Sys.time(), tile_start,
                                           units = "mins")), 1)
      cat("  Progress:", i, "/", nrow(fishnet_sf),
          "| Elapsed:", elapsed, "mins\n")
    }
    
    poly_sf   <- fishnet_sf[i, ]
    poly_vect <- fishnet_vect[i, ]
    
    r_crop <- tryCatch(crop(r, poly_vect),  error = function(e) NULL)
    if (is.null(r_crop)) { rm(poly_sf, poly_vect); gc(); next }
    
    r_mask <- tryCatch(mask(r_crop, poly_vect), error = function(e) NULL)
    if (is.null(r_mask)) { rm(poly_sf, poly_vect, r_crop); gc(); next }
    
    valid_cells <- global(!is.na(r_mask), "sum")[1, 1]
    if (is.na(valid_cells) || valid_cells == 0) {
      empty_count <- empty_count + 1
      rm(poly_sf, poly_vect, r_crop, r_mask); gc()
      next
    }
    
    valid_count <- valid_count + 1
    
    metrics <- tryCatch({
      calculate_lsm(
        r_mask,
        what = c("lsm_c_ed",
                 "lsm_c_pd",
                 "lsm_c_area_mn",
                 "lsm_c_pland")
      )
    }, error = function(e) NULL)
    
    if (is.null(metrics)) {
      rm(poly_sf, poly_vect, r_crop, r_mask); gc()
      next
    }
    
    centroid <- st_coordinates(st_centroid(poly_sf))
    
    metrics$plot_id  <- poly_sf$plot_id
    metrics$center_x <- centroid[1]
    metrics$center_y <- centroid[2]
    metrics$tile_id  <- tile_name
    metrics$scenario <- scenario_name
    
    if (has_row_col) {
      metrics$row <- poly_sf$row
      metrics$col <- poly_sf$col
    }
    
    metrics_list[[length(metrics_list) + 1]] <- metrics
    
    rm(poly_sf, poly_vect, r_crop, r_mask, metrics, centroid)
    gc()
  }
  
  # Save
  if (length(metrics_list) == 0) {
    cat(" -> No valid metrics for this tile, skipping save\n")
    rm(r, fishnet_sf, fishnet_vect, metrics_list); gc()
    next
  }
  
  final_metrics <- bind_rows(metrics_list)
  
  col_order <- c("plot_id", "center_x", "center_y",
                 "layer", "level", "class", "id",
                 "metric", "value", "tile_id", "scenario")
  if (has_row_col) col_order <- c(col_order, "row", "col")
  
  final_metrics <- final_metrics %>% select(any_of(col_order))
  
  write.csv(final_metrics, output_file, row.names = FALSE)
  
  cat("\nValid grids:", valid_count,
      "| Empty grids:", empty_count,
      "| Rows saved:", nrow(final_metrics), "\n")
  cat("Tile time:",
      round(as.numeric(difftime(Sys.time(), tile_start, units = "mins")), 2),
      "mins\n")
  
  rm(r, fishnet_sf, fishnet_vect, metrics_list, final_metrics)
  gc()
}

overall_time <- round(
  as.numeric(difftime(Sys.time(), overall_start, units = "hours")), 2
)
done_count <- length(list.files(
  output_dir, pattern = paste0("_", suffix, "_b\\.csv$")
))

cat("\n====================================================\n")
cat("ALL ALL_PNR TILES COMPLETE\n")
cat("Total runtime:  ", overall_time, "hours\n")
cat("CSVs in output: ", done_count, "\n")
cat("====================================================\n")
