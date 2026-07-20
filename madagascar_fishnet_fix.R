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


# ============================================================
# Diagnose striping in scenario FFI rasters (hh and ap)
# ------------------------------------------------------------
# Current forest FFI is clean, so the fault is in something the
# two scenarios SHARE. Prime suspect: the projected fishnet grid
# drifted off the 5km origin (same failure mode as 20S_040E,
# which had a 936m offset).
#
# PART A checks the FFI CSVs: is the data complete for the tile,
#   and are its cell coordinates on the 5km grid?
# PART B compares the two fishnets and reports any origin offset,
#   plus the true lon/lat of the tile so we confirm we are even
#   looking at the right one.
# ============================================================

library(data.table)
library(sf)

# ---------------- SETTINGS ----------------
tile <- "10N_000E"     # suspect tile — switch once we confirm the real one
grid <- 5000           # 5km grid

ap_csv  <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/all_pnr_results/all_pnr_ap_b_FFI_10m.csv"
hh_csv  <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/holistic_hotspot_results/holistic_hotspot_hh_b_FFI_10m.csv"
cur_csv <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/current_forest_FFI_10m.csv"  # confirm path

proj_fn_path <- "R:/Chapter_3_fragmentation/2026_NEE_R2/fishnets_projected/10N_000E_fishnet_projected.gpkg"
cur_fn_path  <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/fishnet_from10m/10N_000E_fishnet.gpkg"

# Coordinate column names in the FFI CSVs.
# If the printed column list shows different names, change these two.
xcol <- "center_x"
ycol <- "center_y"

# ============================================================
# PART A — INSPECT THE FFI CSVs FOR THE SUSPECT TILE
# ============================================================
inspect_csv <- function(path, label) {
  cat("\n=================", label, "=================\n")
  if (!file.exists(path)) { cat("FILE NOT FOUND:", path, "\n"); return(invisible()) }
  dt <- fread(path)
  cat("Columns:", paste(names(dt), collapse = ", "), "\n")
  
  tcol <- if ("tile_id" %in% names(dt)) "tile_id" else
    names(dt)[grep("tile", names(dt), ignore.case = TRUE)][1]
  cat("Tile column:", tcol, " | unique tiles:", length(unique(dt[[tcol]])),
      " | suspect present:", tile %in% dt[[tcol]], "\n")
  
  sub <- dt[dt[[tcol]] == tile, ]
  cat("Rows for", tile, ":", nrow(sub), "\n")
  if (nrow(sub) == 0) return(invisible())
  
  fficol <- names(sub)[grep("FFI", names(sub), ignore.case = TRUE)][1]
  cat("FFI column:", fficol, " | NA FFI values:", sum(is.na(sub[[fficol]])), "\n")
  
  if (all(c(xcol, ycol) %in% names(sub))) {
    ux <- sort(unique(round(sub[[xcol]])))
    uy <- sort(unique(round(sub[[ycol]])))
    dx <- diff(ux)
    cat("X columns:", length(ux), " | X range", min(ux), "to", max(ux), "\n")
    cat("Origin residual to", grid, "grid:  x =", ux[1] %% grid,
        " y =", uy[1] %% grid, "\n")
    cat("Missing columns inside the CSV (gaps > 1.5 cells):",
        sum(dx > grid * 1.5), "\n")
  } else {
    cat("Coord columns", xcol, "/", ycol,
        "not found — set xcol/ycol from the column list above.\n")
  }
  invisible(sub)
}

inspect_csv(ap_csv,  "ALL_PNR (ap)")
inspect_csv(hh_csv,  "HOLISTIC HOTSPOT (hh)")
inspect_csv(cur_csv, "CURRENT FOREST")

# ============================================================
# PART B — COMPARE THE FISHNETS
# ============================================================
cat("\n\n================= FISHNET COMPARISON =================\n")

cur_fn  <- st_read(cur_fn_path,  quiet = TRUE)
proj_fn <- st_read(proj_fn_path, quiet = TRUE)

cat("Current fishnet:   cells =", nrow(cur_fn),
    " | CRS =", st_crs(cur_fn)$epsg,  "\n")
cat("Projected fishnet: cells =", nrow(proj_fn),
    " | CRS =", st_crs(proj_fn)$epsg, "\n")

# Where is this tile really?  (settles the 10N_000E vs image question)
ll <- st_bbox(st_transform(st_geometry(cur_fn), 4326))
cat("\nReal-world extent of", tile, "(lon/lat):\n")
cat("  lon", round(ll["xmin"], 2), "to", round(ll["xmax"], 2),
    " |  lat", round(ll["ymin"], 2), "to", round(ll["ymax"], 2), "\n")

# Align CRS before comparing grids
if (st_crs(proj_fn) != st_crs(cur_fn)) {
  cat("\nCRS differ — reprojecting projected fishnet to current CRS for comparison.\n")
  proj_fn <- st_transform(proj_fn, st_crs(cur_fn))
}

cc <- st_coordinates(st_centroid(st_geometry(cur_fn)))
cp <- st_coordinates(st_centroid(st_geometry(proj_fn)))

ux_c <- sort(unique(round(cc[, 1], 1))); uy_c <- sort(unique(round(cc[, 2], 1)))
ux_p <- sort(unique(round(cp[, 1], 1))); uy_p <- sort(unique(round(cp[, 2], 1)))

cat("\nCell size (median centroid spacing):\n")
cat("  current   x =", median(diff(ux_c)), " y =", median(diff(uy_c)), "\n")
cat("  projected x =", median(diff(ux_p)), " y =", median(diff(uy_p)), "\n")

cat("\nGrid origin (first centroid):\n")
cat("  current   x =", ux_c[1], " y =", uy_c[1], "\n")
cat("  projected x =", ux_p[1], " y =", uy_p[1], "\n")

off_x <- ux_p[1] - ux_c[1]
off_y <- uy_p[1] - uy_c[1]
cat("\nOrigin offset (projected minus current):\n")
cat("  dx =", off_x, "  dy =", off_y, "\n")
cat("  dx mod", grid, "=", round(off_x %% grid, 1),
    "  dy mod", grid, "=", round(off_y %% grid, 1), "\n")

mx <- off_x %% grid; my <- off_y %% grid
off_in_x <- mx > 1 && mx < grid - 1
off_in_y <- my > 1 && my < grid - 1

if (off_in_x || off_in_y) {
  cat("\n>>> OFFSET DETECTED: the projected fishnet is misaligned to the 5km grid.\n")
  cat("    Same failure mode as 20S_040E. Fix: replace this projected fishnet\n")
  cat("    with a copy of the current forest fishnet, then regenerate the\n")
  cat("    metrics and FFI for", tile, "in both scenarios.\n")
} else {
  cat("\nGrids appear aligned — the offset is NOT in the fishnet origin.\n")
  cat("If the raster still stripes, the fault is in the FFI raster write\n")
  cat("template (snap raster / cell size / extent), which we check next.\n")
}


##############
############## There are other tiles have the issue, but for now
############## Let's only fix the striped one that is obvious on map
# ============================================================
# Fix the striped tile 10N_000E in the FFI CSVs (ap and hh)
# ------------------------------------------------------------
# The per-cell FFI values are correct; only center_x / center_y
# were collapsed by the round(x/5000)*5000 snap. plot_id is
# untouched, so we recover the true coordinates by joining the
# tile's rows to the fishnet centroids on plot_id, then write
# the corrected CSV back in place (original backed up first).
#
# Only 10N_000E rows change. Every other tile is passed through
# unchanged, so you keep a full set to rebuild the FFI rasters.
# ============================================================

library(data.table)
library(sf)

tile <- "10N_000E"
grid <- 5000

fishnet_path <- "R:/Chapter_3_fragmentation/2026_NEE_R2/fishnets_projected/10N_000E_fishnet_projected.gpkg"

files <- list(
  ap = "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/all_pnr_results/all_pnr_ap_b_FFI_10m.csv",
  hh = "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/holistic_hotspot_results/holistic_hotspot_hh_b_FFI_10m.csv"
)

# ---------------- TRUE COORDINATES FROM THE FISHNET ----------------
fn <- st_read(fishnet_path, quiet = TRUE)

pid_col <- intersect(c("plot_id", "PLOT_ID", "plotid", "FID", "fid", "id"), names(fn))[1]
if (is.na(pid_col)) {
  stop(paste("No plot_id column found in fishnet. Columns:",
             paste(names(fn), collapse = ", ")))
}

cent <- st_coordinates(st_centroid(st_geometry(fn)))
cent_dt <- data.table(pkey = as.character(fn[[pid_col]]),
                      cx = cent[, 1], cy = cent[, 2])
cent_dt <- cent_dt[!duplicated(pkey)]          # one centroid per cell
cat("Fishnet centroids loaded:", nrow(cent_dt),
    "cells | plot_id column:", pid_col, "\n\n")

# ---------------- FIX ONE FILE ----------------
fix_tile <- function(csv_path, label) {
  cat("=======================", label, "=======================\n")
  dt  <- fread(csv_path)
  idx <- dt$tile_id == tile
  cat("Rows for", tile, ":", sum(idx),
      "| distinct plot_id:", uniqueN(dt$plot_id[idx]), "\n")
  
  bx <- sort(unique(round(dt$center_x[idx])))
  cat("BEFORE: columns", length(bx), "| x residual", bx[1] %% grid, "\n")
  
  sub <- dt[idx]
  sub[, pkey := as.character(plot_id)]
  sub2 <- merge(sub, cent_dt, by = "pkey", all.x = TRUE, sort = FALSE)
  
  n_na <- sum(is.na(sub2$cx))
  if (n_na > 0) {
    cat("WARNING:", n_na,
        "rows did not match a fishnet cell — file left unchanged.\n\n")
    return(invisible())
  }
  
  sub2[, center_x := cx][, center_y := cy]
  sub2[, c("pkey", "cx", "cy") := NULL]
  
  out <- rbind(dt[!idx], sub2, use.names = TRUE)
  setcolorder(out, names(dt))
  setorder(out, tile_id, plot_id)
  
  ax <- sort(unique(round(out$center_x[out$tile_id == tile])))
  cat("AFTER:  columns", length(ax), "| x residual", ax[1] %% grid, "\n")
  
  # back up once, then replace in place
  bak <- paste0(csv_path, ".bak_before_10N000E_fix")
  if (!file.exists(bak)) file.copy(csv_path, bak)
  fwrite(out, csv_path)
  cat("Backup :", bak, "\n")
  cat("Replaced:", csv_path, "\n")
  
  # standalone corrected tile, handy for a quick single-tile raster check
  tile_out <- file.path(dirname(csv_path),
                        paste0(tile, "_", label, "_FFI_fixed.csv"))
  fwrite(sub2, tile_out)
  cat("Tile-only file:", tile_out, "\n\n")
}

fix_tile(files$ap, "ap")
fix_tile(files$hh, "hh")

cat("Done. Re-run make_raster on the corrected CSVs and the stripes\n")
cat("over", tile, "should be gone. Originals are kept as .bak files.\n")











