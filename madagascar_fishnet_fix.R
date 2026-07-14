library(terra)
library(sf)
library(landscapemetrics)
library(dplyr)
library(data.table)

scenario_name <- "holistic_hotspot"
suffix        <- "hh"
tile_name     <- "20S_040E"

raster_dir  <- "R:/Chapter_3_fragmentation/2026_NEE_R2/binary_forest_scenarios/holistic_hotspot"
output_dir  <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FutureScenario_metrics/holistic_hotspot_binary_metrics"

# KEY FIX: use the current forest fishnet, not the projected one
fishnet_path <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/fishnet_from10m/20S_040E_fishnet.gpkg"
output_file  <- file.path(output_dir, paste0(tile_name, "_", suffix, "_b.csv"))

tile_start <- Sys.time()

cat("====================================================\n")
cat("Tile:", tile_name, "| Using CURRENT FOREST fishnet\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("====================================================\n")

# Load raster
r_file <- file.path(raster_dir, paste0(tile_name, "_", suffix, "_b.tif"))
r      <- rast(r_file)
cat(" -> Raster loaded |", nrow(r), "x", ncol(r), "\n")

# Load CURRENT FOREST fishnet — same grid as current metrics
fishnet_sf   <- st_read(fishnet_path, quiet = TRUE)
fishnet_sf   <- st_transform(fishnet_sf, crs = crs(r))
cat(" -> Current forest fishnet loaded:", nrow(fishnet_sf), "cells\n")

# Filter to hh extent only
hh_extent_poly <- st_as_sf(as.polygons(ext(r), crs = crs(r)))
fishnet_sf     <- fishnet_sf[
  st_intersects(fishnet_sf, hh_extent_poly, sparse = FALSE)[, 1], ]
cat(" -> After extent filter:", nrow(fishnet_sf), "cells\n\n")

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
        "| Valid:", valid_count,
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
    calculate_lsm(r_mask,
                  what = c("lsm_c_ed", "lsm_c_pd",
                           "lsm_c_area_mn", "lsm_c_pland"))
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

cat("\nValid grids:", valid_count,
    "| Empty grids:", empty_count, "\n")

if (length(metrics_list) == 0) {
  cat("No valid metrics produced\n")
  rm(r, fishnet_sf, fishnet_vect, metrics_list); gc()
  stop("Stopping.")
}

final_metrics <- bind_rows(metrics_list)

col_order <- c("plot_id", "center_x", "center_y",
               "layer", "level", "class", "id",
               "metric", "value", "tile_id", "scenario")
if (has_row_col) col_order <- c(col_order, "row", "col")

final_metrics <- final_metrics %>% select(any_of(col_order))

fwrite(final_metrics, output_file)
cat("Saved:", nrow(final_metrics), "rows to\n", output_file, "\n")
cat("Tile time:",
    round(as.numeric(difftime(Sys.time(), tile_start, units = "mins")), 2),
    "mins\n")

rm(r, fishnet_sf, fishnet_vect, metrics_list, final_metrics)
gc()