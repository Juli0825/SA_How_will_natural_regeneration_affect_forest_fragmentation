##############################################################
# Remake ONE tile's all_pnr metrics using the CURRENT FOREST
# fishnet, so the scenario grid matches current forest exactly.
# Same tile and same fishnet fix as the holistic hotspot run.
##############################################################

library(terra)
library(sf)
library(landscapemetrics)
library(dplyr)
library(data.table)

# ================= CONFIG =================
tile_name <- "20S_050W"

# --- all_pnr ---
scenario_name <- "all_pnr"
suffix        <- "ap"
raster_dir    <- "R:/Chapter_3_fragmentation/2026_NEE_R2/binary_forest_scenarios/all_pnr"
output_dir    <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FutureScenario_metrics/all_pnr_metrics"
# NOTE: output_dir must be the SAME folder your FFI script reads all_pnr metrics from.

current_fishnet_dir   <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/fishnet_from10m"
projected_fishnet_dir <- "R:/Chapter_3_fragmentation/2026_NEE_R2/fishnets_projected"
current_metrics_dir   <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/frag_metrics_current_10m"

trim_to_current <- TRUE   # FALSE to process every cell in the fishnet

fishnet_path <- file.path(current_fishnet_dir, paste0(tile_name, "_fishnet.gpkg"))
output_file  <- file.path(output_dir, paste0(tile_name, "_", suffix, "_b.csv"))
r_file       <- file.path(raster_dir, paste0(tile_name, "_", suffix, "_b.tif"))

tile_start <- Sys.time()
cat("====================================================\n")
cat("Tile:", tile_name, "| Scenario:", scenario_name, "\n")
cat("Using CURRENT FOREST fishnet\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("====================================================\n")

# ---------- 1. OFFSET CHECK (confirm the diagnosis) ----------
proj_path <- file.path(projected_fishnet_dir,
                       paste0(tile_name, "_fishnet_projected.gpkg"))
if (file.exists(proj_path)) {
  cf <- st_read(fishnet_path, quiet = TRUE)
  pf <- st_read(proj_path,    quiet = TRUE)
  if (!is.na(st_crs(pf)) && !is.na(st_crs(cf)) && st_crs(pf) != st_crs(cf))
    pf <- st_transform(pf, st_crs(cf))
  cc <- st_coordinates(st_centroid(st_geometry(cf)))
  cp <- st_coordinates(st_centroid(st_geometry(pf)))
  grid <- 5000
  wrap <- function(d) { d <- d %% grid; ifelse(d > grid / 2, d - grid, d) }
  dx <- wrap((min(cp[, 1]) %% grid) - (min(cc[, 1]) %% grid))
  dy <- wrap((min(cp[, 2]) %% grid) - (min(cc[, 2]) %% grid))
  cat(sprintf(" -> Sub cell fishnet offset (projected minus current): dx = %.1f m, dy = %.1f m\n",
              dx, dy))
  if (max(abs(c(dx, dy))) > 1) {
    cat("    Offset confirmed, this rerun is the right fix.\n")
  } else {
    cat("    NOTE: fishnets look aligned, so the cause may be something else.\n")
  }
  rm(cf, pf, cc, cp); gc()
} else {
  cat(" -> No projected fishnet found to compare, skipping offset check.\n")
}

# ---------- 2. BACKUP THE EXISTING METRICS ----------
if (file.exists(output_file)) {
  bak <- paste0(output_file, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  file.copy(output_file, bak)
  cat(" -> Backed up existing metrics to:\n    ", bak, "\n")
} else {
  cat(" -> No existing metrics file, writing fresh.\n")
}

# ---------- 3. LOAD RASTER AND FISHNET ----------
r <- rast(r_file)
cat(" -> Raster loaded |", nrow(r), "x", ncol(r), "\n")

fishnet_sf <- st_read(fishnet_path, quiet = TRUE)
fishnet_sf <- st_transform(fishnet_sf, crs = crs(r))
cat(" -> Current forest fishnet loaded:", nrow(fishnet_sf), "cells\n")

extent_poly <- st_as_sf(as.polygons(ext(r), crs = crs(r)))
fishnet_sf  <- fishnet_sf[
  st_intersects(fishnet_sf, extent_poly, sparse = FALSE)[, 1], ]
cat(" -> After extent filter:", nrow(fishnet_sf), "cells\n")

if (trim_to_current) {
  cm_files <- list.files(current_metrics_dir,
                         pattern = paste0(tile_name, ".*\\.csv$"),
                         full.names = TRUE)
  if (length(cm_files) > 0) {
    cm  <- fread(cm_files[1])
    keep <- unique(cm$plot_id)
    before <- nrow(fishnet_sf)
    fishnet_sf <- fishnet_sf[fishnet_sf$plot_id %in% keep, ]
    cat(" -> Trimmed to current forest plot_ids:", before, "->",
        nrow(fishnet_sf), "cells\n")
    rm(cm); gc()
  } else {
    cat(" -> No current metrics file matched, skipping trim.\n")
  }
}
cat("\n")

fishnet_vect <- vect(fishnet_sf)
has_row_col  <- all(c("row", "col") %in% names(fishnet_sf))

# ---------- 4. CELL LOOP ----------
metrics_list <- list()
empty_count  <- 0
valid_count  <- 0

for (i in seq_len(nrow(fishnet_sf))) {
  
  if (i %% 500 == 0) {
    elapsed <- round(as.numeric(difftime(Sys.time(), tile_start, units = "mins")), 1)
    cat("  Progress:", i, "/", nrow(fishnet_sf),
        "| Valid:", valid_count, "| Elapsed:", elapsed, "mins\n")
    gc()
  }
  
  poly_sf   <- fishnet_sf[i, ]
  poly_vect <- fishnet_vect[i, ]
  
  r_crop <- tryCatch(crop(r, poly_vect), error = function(e) NULL)
  if (is.null(r_crop)) next
  
  r_mask <- tryCatch(mask(r_crop, poly_vect), error = function(e) NULL)
  if (is.null(r_mask)) next
  
  valid_cells <- global(!is.na(r_mask), "sum")[1, 1]
  if (is.na(valid_cells) || valid_cells == 0) {
    empty_count <- empty_count + 1
    next
  }
  valid_count <- valid_count + 1
  
  metrics <- tryCatch(
    calculate_lsm(r_mask,
                  what = c("lsm_c_ed", "lsm_c_pd",
                           "lsm_c_area_mn", "lsm_c_pland")),
    error = function(e) NULL)
  if (is.null(metrics)) next
  
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
}

cat("\nValid grids:", valid_count, "| Empty grids:", empty_count, "\n")

if (length(metrics_list) == 0) {
  rm(r, fishnet_sf, fishnet_vect, metrics_list); gc()
  stop("No valid metrics produced, nothing written. Original file untouched.")
}

final_metrics <- bind_rows(metrics_list)

col_order <- c("plot_id", "center_x", "center_y",
               "layer", "level", "class", "id",
               "metric", "value", "tile_id", "scenario")
if (has_row_col) col_order <- c(col_order, "row", "col")
final_metrics <- final_metrics %>% select(any_of(col_order))

fwrite(final_metrics, output_file)
cat("Saved:", nrow(final_metrics), "rows to\n", output_file, "\n")

# ---------- 5. QUICK CHECK AGAINST CURRENT FOREST ----------
cm_files <- list.files(current_metrics_dir,
                       pattern = paste0(tile_name, ".*\\.csv$"), full.names = TRUE)
if (length(cm_files) > 0) {
  cm <- fread(cm_files[1])
  n_cur <- uniqueN(cm[class == 1 & metric == "ed", plot_id])
  n_new <- uniqueN(final_metrics[final_metrics$class == 1 &
                                   final_metrics$metric == "ed", "plot_id"])
  cat(sprintf("\nCells with class 1 ED: current %d, new scenario %d (%.1f%%)\n",
              n_cur, n_new, 100 * n_new / n_cur))
  cat("Close to or above 100% means the grids now line up.\n")
}

cat("Tile time:",
    round(as.numeric(difftime(Sys.time(), tile_start, units = "mins")), 2), "mins\n")
cat("\nNEXT: rerun the all_pnr FFI script so the new metrics flow through\n")
cat("to the FFI CSV, delta CSV and rasters.\n")

rm(r, fishnet_sf, fishnet_vect, metrics_list, final_metrics); gc()