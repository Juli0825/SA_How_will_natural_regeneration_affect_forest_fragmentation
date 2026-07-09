library(terra)
library(sf)
library(landscapemetrics)
library(dplyr)
library(data.table)

scenario_name <- "holistic_hotspot"
suffix        <- "hh"

raster_dir  <- "R:/Chapter_3_fragmentation/2026_NEE_R2/binary_forest_scenarios/holistic_hotspot"
fishnet_dir <- "R:/Chapter_3_fragmentation/2026_NEE_R2/fishnets_projected"
output_dir  <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FutureScenario_metrics/holistic_hotspot_binary_metrics"

# The 5 suspicious tiles
target_tiles <- c("30N_110W")

overall_start <- Sys.time()

for (tile_name in target_tiles) {
  
  tile_start  <- Sys.time()
  output_file <- file.path(output_dir,
                           paste0(tile_name, "_", suffix, "_b.csv"))
  
  cat("\n====================================================\n")
  cat("Tile:", tile_name, "\n")
  cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("====================================================\n")
  
  # Overwrite any previously copied current forest CSV
  if (file.exists(output_file)) {
    file.remove(output_file)
    cat(" -> Removed previous copied file\n")
  }
  
  # Load scenario raster
  r_file <- file.path(raster_dir, paste0(tile_name, "_", suffix, "_b.tif"))
  if (!file.exists(r_file)) {
    cat(" -> Raster not found, skipping\n")
    next
  }
  
  r <- tryCatch(rast(r_file), error = function(e) {
    cat(" -> ERROR reading raster:", e$message, "\n")
    return(NULL)
  })
  if (is.null(r)) next
  
  cat(" -> Raster loaded |", nrow(r), "x", ncol(r), "\n")
  cat("    hh_b extent:", as.vector(ext(r)), "\n")
  
  # Load full fishnet
  fishnet_path <- file.path(fishnet_dir,
                            paste0(tile_name, "_fishnet_projected.gpkg"))
  if (!file.exists(fishnet_path)) {
    cat(" -> Fishnet not found, skipping\n")
    rm(r); gc(); next
  }
  
  fishnet_sf <- st_read(fishnet_path, quiet = TRUE)
  fishnet_sf <- st_transform(fishnet_sf, crs = crs(r))
  cat(" -> Full fishnet loaded:", nrow(fishnet_sf), "cells\n")
  
  # Filter fishnet to cells intersecting hh_b extent only
  # No current forest plot_id filter this time
  hh_extent_poly <- st_as_sf(as.polygons(ext(r), crs = crs(r)))
  
  fishnet_sf <- fishnet_sf[
    st_intersects(fishnet_sf, hh_extent_poly, sparse = FALSE)[, 1], ]
  
  cat(" -> After extent filter:", nrow(fishnet_sf),
      "cells intersect hh_b extent\n\n")
  
  if (nrow(fishnet_sf) == 0) {
    cat(" -> No cells inside hh extent, skipping\n")
    rm(r); gc(); next
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
  
  cat("\nValid grids:", valid_count,
      "| Empty grids:", empty_count, "\n")
  
  if (length(metrics_list) == 0) {
    cat(" -> No valid metrics for this tile\n")
    rm(r, fishnet_sf, fishnet_vect, metrics_list); gc()
    next
  }
  
  final_metrics <- bind_rows(metrics_list)
  
  col_order <- c("plot_id", "center_x", "center_y",
                 "layer", "level", "class", "id",
                 "metric", "value", "tile_id", "scenario")
  if (has_row_col) col_order <- c(col_order, "row", "col")
  
  final_metrics <- final_metrics %>% select(any_of(col_order))
  
  fwrite(final_metrics, output_file)
  
  cat("Saved:", nrow(final_metrics), "rows\n")
  cat("Unique plot_ids:", length(unique(final_metrics$plot_id)), "\n")
  cat("Tile time:",
      round(as.numeric(difftime(Sys.time(), tile_start, units = "mins")), 2),
      "mins\n")
  
  rm(r, fishnet_sf, fishnet_vect, metrics_list, final_metrics)
  gc()
}

overall_time <- round(
  as.numeric(difftime(Sys.time(), overall_start, units = "hours")), 2
)
cat("\n====================================================\n")
cat("ALL 5 TILES DONE\n")
cat("Total runtime:", overall_time, "hours\n")
cat("====================================================\n")
cat("\nNext step: check what plot_ids came back and see\n")
cat("if they match current forest plot_ids or are new ones\n")

library(data.table)

current_metrics_dir <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/frag_metrics_current_10m"
output_dir          <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FutureScenario_metrics/holistic_hotspot_binary_metrics"

tile_name <- "30N_100W"

hh <- fread(file.path(output_dir,
                      paste0(tile_name, "_hh_b.csv")))
cf <- fread(file.path(current_metrics_dir,
                      paste0(tile_name, "_metrics.csv")))

cf_forest <- cf[class == 1]

cat("hh_b plot_id range:\n")
cat("  min:", min(hh$plot_id), "\n")
cat("  max:", max(hh$plot_id), "\n")
cat("  first 10:", paste(sort(unique(hh$plot_id))[1:10], collapse = ", "), "\n\n")

cat("Current forest plot_id range:\n")
cat("  min:", min(cf_forest$plot_id), "\n")
cat("  max:", max(cf_forest$plot_id), "\n")
cat("  first 10:", paste(sort(unique(cf_forest$plot_id))[1:10], collapse = ", "), "\n\n")

# Check center_x and center_y ranges too
cat("hh_b center_x range:", min(hh$center_x), "to", max(hh$center_x), "\n")
cat("cf  center_x range: ", min(cf$center_x),  "to", max(cf$center_x),  "\n\n")

cat("hh_b center_y range:", min(hh$center_y), "to", max(hh$center_y), "\n")
cat("cf  center_y range: ", min(cf$center_y),  "to", max(cf$center_y),  "\n\n")

# The real test — do coordinates match even if plot_ids dont?
cat("Are center_x values overlapping?\n")
hh_x <- unique(hh$center_x)
cf_x <- unique(cf$center_x)
cat("  Matching center_x values:", length(intersect(hh_x, cf_x)), "\n")

rm(hh, cf, cf_forest); gc()







