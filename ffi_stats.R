library(sf)
library(terra)

tile_name <- "20S_040E"

fishnet_current <- st_read(
  paste0("R:/Chapter_3_fragmentation/frag_2026_exct_median/fishnet_from10m/",
         tile_name, "_fishnet.gpkg"),
  quiet = TRUE)

fishnet_scenario <- st_read(
  paste0("R:/Chapter_3_fragmentation/2026_NEE_R2/fishnets_projected/",
         tile_name, "_fishnet_projected.gpkg"),
  quiet = TRUE)

cat("===== Current forest fishnet =====\n")
cat("CRS:       ", st_crs(fishnet_current)$proj4string, "\n")
cat("Rows:      ", nrow(fishnet_current), "\n")
cat("Bbox:      ", as.vector(st_bbox(fishnet_current)), "\n")
cat("Columns:   ", paste(names(fishnet_current), collapse = ", "), "\n")

cat("\n===== Scenario fishnet projected =====\n")
cat("CRS:       ", st_crs(fishnet_scenario)$proj4string, "\n")
cat("Rows:      ", nrow(fishnet_scenario), "\n")
cat("Bbox:      ", as.vector(st_bbox(fishnet_scenario)), "\n")
cat("Columns:   ", paste(names(fishnet_scenario), collapse = ", "), "\n")

# Check centroid coordinates of first 5 cells in each
cat("\n===== First 5 cell centroids in current forest fishnet =====\n")
coords_cf <- st_coordinates(st_centroid(fishnet_current[1:5, ]))
print(round(coords_cf, 2))

cat("\n===== First 5 cell centroids in scenario fishnet =====\n")
coords_sc <- st_coordinates(st_centroid(fishnet_scenario[1:5, ]))
print(round(coords_sc, 2))

# Check if plot_ids match
cat("\n===== plot_id comparison =====\n")
cat("Current forest plot_id range:",
    min(fishnet_current$plot_id), "to", max(fishnet_current$plot_id), "\n")
cat("Scenario    plot_id range:",
    min(fishnet_scenario$plot_id), "to", max(fishnet_scenario$plot_id), "\n")

# Check coordinate offset between matching cells
cat("\n===== Coordinate offset check =====\n")
n_check <- min(10, nrow(fishnet_current), nrow(fishnet_scenario))
cf_centroids <- st_coordinates(st_centroid(fishnet_current[1:n_check, ]))
sc_centroids <- st_coordinates(st_centroid(fishnet_scenario[1:n_check, ]))

for (i in 1:n_check) {
  cat(sprintf("  Cell %d: cf=(%.1f, %.1f)  sc=(%.1f, %.1f)  offset=(%.1f, %.1f)\n",
              i,
              cf_centroids[i, "X"], cf_centroids[i, "Y"],
              sc_centroids[i, "X"], sc_centroids[i, "Y"],
              sc_centroids[i, "X"] - cf_centroids[i, "X"],
              sc_centroids[i, "Y"] - cf_centroids[i, "Y"]))
}

library(sf)
library(terra)

tile_name <- "20S_040E"

fishnet_current <- st_read(
  paste0("R:/Chapter_3_fragmentation/frag_2026_exct_median/fishnet_from10m/",
         tile_name, "_fishnet.gpkg"),
  quiet = TRUE)

fishnet_scenario <- st_read(
  paste0("R:/Chapter_3_fragmentation/2026_NEE_R2/fishnets_projected/",
         tile_name, "_fishnet_projected.gpkg"),
  quiet = TRUE)

# Plot both fishnets overlaid
# Sample a small subset so it is not too slow
set.seed(42)
cf_sample <- fishnet_current[sample(nrow(fishnet_current), 200), ]
sc_sample <- fishnet_scenario[sample(nrow(fishnet_scenario), 200), ]

plot(st_geometry(cf_sample),
     border = "blue",
     lwd    = 0.5,
     main   = paste(tile_name, "fishnet comparison"),
     sub    = "Blue = current forest fishnet | Red = scenario fishnet projected")

plot(st_geometry(sc_sample),
     border = "red",
     lwd    = 0.5,
     add    = TRUE)

legend("topleft",
       legend = c("Current forest fishnet", "Scenario fishnet projected"),
       col    = c("blue", "red"),
       lty    = 1,
       lwd    = 1.5,
       cex    = 0.8)

#######################
#######################
library(sf)
library(terra)

tile_name <- "20S_040E"

r <- rast(paste0("R:/Chapter_3_fragmentation/2026_NEE_R2/",
                 "binary_forest_scenarios/holistic_hotspot/",
                 tile_name, "_hh_b.tif"))

fishnet_cf <- st_read(
  paste0("R:/Chapter_3_fragmentation/frag_2026_exct_median/",
         "fishnet_from10m/", tile_name, "_fishnet.gpkg"),
  quiet = TRUE)

fishnet_cf <- st_transform(fishnet_cf, crs = crs(r))

cat("hh_b raster extent:\n")
cat("  xmin:", xmin(ext(r)), "xmax:", xmax(ext(r)), "\n")
cat("  ymin:", ymin(ext(r)), "ymax:", ymax(ext(r)), "\n\n")

cat("Current forest fishnet bbox:\n")
bb <- st_bbox(fishnet_cf)
cat("  xmin:", bb["xmin"], "xmax:", bb["xmax"], "\n")
cat("  ymin:", bb["ymin"], "ymax:", bb["ymax"], "\n\n")

# How many fishnet cells intersect the hh_b extent?
hh_poly    <- st_as_sf(as.polygons(ext(r), crs = crs(r)))
n_inside   <- sum(st_intersects(fishnet_cf, hh_poly, sparse = FALSE)[, 1])

cat("Fishnet cells intersecting hh_b extent:", n_inside, "\n")
cat("Total fishnet cells:                   ", nrow(fishnet_cf), "\n")
cat("Expected (from previous hh_b run):      ~5,427\n\n")

cat("Does current forest fishnet FULLY cover hh_b extent?\n")
cat("  hh_b ymax:", ymax(ext(r)), "\n")
cat("  fishnet ymax:", bb["ymax"], "\n")
cat("  Covered:", bb["ymax"] >= ymax(ext(r)), "\n")
cat("  hh_b ymin:", ymin(ext(r)), "\n")
cat("  fishnet ymin:", bb["ymin"], "\n")
cat("  Covered:", bb["ymin"] <= ymin(ext(r)), "\n")


library(terra)

tile_name <- "20S_040E"

r_hh <- rast(paste0("R:/Chapter_3_fragmentation/2026_NEE_R2/",
                    "binary_forest_scenarios/holistic_hotspot/",
                    tile_name, "_hh_b.tif"))

r_ap <- rast(paste0("R:/Chapter_3_fragmentation/2026_NEE_R2/",
                    "binary_forest_scenarios/all_pnr/",
                    tile_name, "_ap_b.tif"))

cat("hh_b extent:\n")
cat("  xmin:", xmin(ext(r_hh)), "xmax:", xmax(ext(r_hh)), "\n")
cat("  ymin:", ymin(ext(r_hh)), "ymax:", ymax(ext(r_hh)), "\n\n")

cat("all_pnr extent:\n")
cat("  xmin:", xmin(ext(r_ap)), "xmax:", xmax(ext(r_ap)), "\n")
cat("  ymin:", ymin(ext(r_ap)), "ymax:", ymax(ext(r_ap)), "\n\n")

cat("Extents identical?\n")
cat("  xmin match:", xmin(ext(r_hh)) == xmin(ext(r_ap)), "\n")
cat("  xmax match:", xmax(ext(r_hh)) == xmax(ext(r_ap)), "\n")
cat("  ymin match:", ymin(ext(r_hh)) == ymin(ext(r_ap)), "\n")
cat("  ymax match:", ymax(ext(r_hh)) == ymax(ext(r_ap)), "\n")











##### stats summary
library(dplyr)

output_dir <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results"
scenario_label <- "holistic_hotspot_hh_b"
cell_area      <- 25

# Load from saved CSVs
cat("Loading saved results...\n")
delta_ffi <- read.csv(file.path(output_dir,
                                "delta_FFI_holistic_hotspot_hh_b_10m.csv"))
cat("Rows loaded:", format(nrow(delta_ffi), big.mark = ","), "\n\n")

# Classify
delta_ffi <- delta_ffi %>%
  mutate(
    delta_class = case_when(
      delta_FFI < -0.001 ~ -1L,
      delta_FFI >  0.001 ~  1L,
      TRUE               ~  0L
    )
  )

n_total <- nrow(delta_ffi)
n_dec   <- sum(delta_ffi$delta_class == -1, na.rm = TRUE)
n_unch  <- sum(delta_ffi$delta_class ==  0, na.rm = TRUE)
n_inc   <- sum(delta_ffi$delta_class ==  1, na.rm = TRUE)

# ====================================================
# TABLE 1: OVERALL FFI SUMMARY
# ====================================================

ffi_summary <- data.frame(
  scenario          = scenario_label,
  mean_FFI_current  = round(mean(delta_ffi$FFI_current,  na.rm = TRUE), 4),
  mean_FFI_scenario = round(mean(delta_ffi$FFI_scenario, na.rm = TRUE), 4),
  mean_delta_FFI    = round(mean(delta_ffi$delta_FFI,    na.rm = TRUE), 4),
  pct_change_FFI    = round(mean(delta_ffi$delta_FFI,    na.rm = TRUE) /
                              mean(delta_ffi$FFI_current,  na.rm = TRUE) * 100, 2),
  median_delta_FFI  = round(median(delta_ffi$delta_FFI,  na.rm = TRUE), 4),
  min_delta_FFI     = round(min(delta_ffi$delta_FFI,     na.rm = TRUE), 4),
  max_delta_FFI     = round(max(delta_ffi$delta_FFI,     na.rm = TRUE), 4),
  n_total_plots     = n_total,
  total_area_km2    = n_total * cell_area
)

write.csv(ffi_summary,
          file.path(output_dir,
                    paste0("summary_FFI_overall_", scenario_label, ".csv")),
          row.names = FALSE)
cat("Saved: summary_FFI_overall\n")

# ====================================================
# TABLE 2: SPATIAL BREAKDOWN
# ====================================================

spatial_summary <- data.frame(
  scenario     = scenario_label,
  category     = c("FFI decreased", "unchanged", "FFI increased", "total"),
  n_cells      = c(n_dec, n_unch, n_inc, n_total),
  area_km2     = c(n_dec, n_unch, n_inc, n_total) * cell_area,
  pct_of_total = round(c(n_dec, n_unch, n_inc, n_total) / n_total * 100, 1)
)

write.csv(spatial_summary,
          file.path(output_dir,
                    paste0("summary_spatial_", scenario_label, ".csv")),
          row.names = FALSE)
cat("Saved: summary_spatial\n")

# ====================================================
# TABLE 3: INDIVIDUAL METRIC CHANGES
# ====================================================

metrics_info <- list(
  list(name = "FFI",     cur = "FFI_current",  del = "delta_FFI",
       improve = "negative"),
  list(name = "ED",      cur = "ED_current",   del = "delta_ED",
       improve = "negative"),
  list(name = "PD",      cur = "PD_current",   del = "delta_PD",
       improve = "negative"),
  list(name = "MPA",     cur = "MPA_current",  del = "delta_MPA",
       improve = "positive")
)

metric_summary <- do.call(rbind, lapply(metrics_info, function(m) {
  cur_val <- mean(delta_ffi[[m$cur]], na.rm = TRUE)
  del_val <- mean(delta_ffi[[m$del]], na.rm = TRUE)
  pct     <- del_val / cur_val * 100
  direction <- if (m$improve == "negative") {
    ifelse(del_val < 0, "improved", "worsened")
  } else {
    ifelse(del_val > 0, "improved", "worsened")
  }
  data.frame(
    scenario          = scenario_label,
    metric            = m$name,
    mean_current      = round(cur_val, 4),
    mean_scenario     = round(cur_val + del_val, 4),
    mean_delta        = round(del_val, 4),
    pct_change        = round(pct, 2),
    direction         = direction
  )
}))

write.csv(metric_summary,
          file.path(output_dir,
                    paste0("summary_metrics_", scenario_label, ".csv")),
          row.names = FALSE)
cat("Saved: summary_metrics\n")

# ====================================================
# TABLE 4: DISTRIBUTION BUCKETS
# ====================================================

breaks <- c(-Inf, -0.1, -0.05, -0.01, -0.001,
            0.001, 0.01, 0.05, 0.1, Inf)
labels <- c("< -0.10", "-0.10 to -0.05", "-0.05 to -0.01",
            "-0.01 to -0.001", "-0.001 to +0.001",
            "+0.001 to +0.01", "+0.01 to +0.05",
            "+0.05 to +0.10", "> +0.10")

delta_ffi$bucket <- cut(delta_ffi$delta_FFI,
                        breaks = breaks, labels = labels)

distribution_summary <- delta_ffi %>%
  group_by(bucket) %>%
  summarise(n_cells  = n(), .groups = "drop") %>%
  mutate(
    scenario     = scenario_label,
    area_km2     = n_cells * cell_area,
    pct_of_total = round(n_cells / n_total * 100, 1)
  ) %>%
  select(scenario, bucket, n_cells, area_km2, pct_of_total)

write.csv(distribution_summary,
          file.path(output_dir,
                    paste0("summary_distribution_", scenario_label, ".csv")),
          row.names = FALSE)
cat("Saved: summary_distribution\n")

cat("\n===== ALL SUMMARY CSVs SAVED =====\n")
cat("Location:", output_dir, "\n")
cat("1. summary_FFI_overall_",    scenario_label, ".csv\n", sep = "")
cat("2. summary_spatial_",        scenario_label, ".csv\n", sep = "")
cat("3. summary_metrics_",        scenario_label, ".csv\n", sep = "")
cat("4. summary_distribution_",   scenario_label, ".csv\n", sep = "")



############ stats in grids show a change
library(dplyr)

# Split into changed vs unchanged
decreased <- delta_ffi %>% filter(delta_class == -1)
increased <- delta_ffi %>% filter(delta_class ==  1)
changed   <- delta_ffi %>% filter(delta_class != 0)

cat("============================================================\n")
cat("  WITHIN CHANGED GRIDS ONLY\n")
cat("============================================================\n\n")

cat("--- GRIDS WHERE FFI DECREASED (holistic hotspot improved landscape) ---\n")
cat(sprintf("  n grids:             %s\n",
            format(nrow(decreased), big.mark = ",")))
cat(sprintf("  Mean FFI current:    %.4f\n",
            mean(decreased$FFI_current,  na.rm = TRUE)))
cat(sprintf("  Mean FFI scenario:   %.4f\n",
            mean(decreased$FFI_scenario, na.rm = TRUE)))
cat(sprintf("  Mean delta FFI:      %+.4f\n",
            mean(decreased$delta_FFI,    na.rm = TRUE)))
cat(sprintf("  Pct change in FFI:   %+.2f%%\n",
            mean(decreased$delta_FFI,    na.rm = TRUE) /
              mean(decreased$FFI_current,  na.rm = TRUE) * 100))
cat(sprintf("  Range delta FFI:     %+.4f  to  %+.4f\n",
            min(decreased$delta_FFI, na.rm = TRUE),
            max(decreased$delta_FFI, na.rm = TRUE)))
cat(sprintf("  Median delta FFI:    %+.4f\n",
            median(decreased$delta_FFI, na.rm = TRUE)))

cat("\n  Individual metrics in decreased grids:\n")
cat(sprintf("    ED  change: %+.4f  (%+.2f%%)\n",
            mean(decreased$delta_ED,  na.rm = TRUE),
            mean(decreased$delta_ED,  na.rm = TRUE) /
              mean(decreased$ED_current, na.rm = TRUE) * 100))
cat(sprintf("    PD  change: %+.4f  (%+.2f%%)\n",
            mean(decreased$delta_PD,  na.rm = TRUE),
            mean(decreased$delta_PD,  na.rm = TRUE) /
              mean(decreased$PD_current, na.rm = TRUE) * 100))
cat(sprintf("    MPA change: %+.4f  (%+.2f%%)\n",
            mean(decreased$delta_MPA, na.rm = TRUE),
            mean(decreased$delta_MPA, na.rm = TRUE) /
              mean(decreased$MPA_current, na.rm = TRUE) * 100))

cat("\n--- GRIDS WHERE FFI INCREASED ---\n")
cat(sprintf("  n grids:             %s\n",
            format(nrow(increased), big.mark = ",")))
cat(sprintf("  Mean delta FFI:      %+.4f\n",
            mean(increased$delta_FFI, na.rm = TRUE)))
cat(sprintf("  Pct change in FFI:   %+.2f%%\n",
            mean(increased$delta_FFI,   na.rm = TRUE) /
              mean(increased$FFI_current, na.rm = TRUE) * 100))

cat("\n--- ALL CHANGED GRIDS COMBINED ---\n")
cat(sprintf("  n grids:             %s  (%.1f%% of total)\n",
            format(nrow(changed), big.mark = ","),
            nrow(changed) / nrow(delta_ffi) * 100))
cat(sprintf("  Mean delta FFI:      %+.4f\n",
            mean(changed$delta_FFI, na.rm = TRUE)))
cat(sprintf("  Pct change in FFI:   %+.2f%%\n",
            mean(changed$delta_FFI,   na.rm = TRUE) /
              mean(changed$FFI_current, na.rm = TRUE) * 100))

cat("\n============================================================\n")
cat("SUGGESTED REPORTING:\n")
cat("  The holistic hotspot scenario affected X% of tropical forest\n")
cat("  grid cells. Among those cells where FFI changed, Y% showed\n")
cat("  a decrease in fragmentation, with a mean FFI reduction of\n")
cat("  Z% relative to the current baseline.\n")
cat("============================================================\n")

# Save
changed_summary <- data.frame(
  scenario         = "holistic_hotspot_hh_b",
  group            = c("FFI decreased", "FFI increased", "all changed"),
  n_grids          = c(nrow(decreased), nrow(increased), nrow(changed)),
  pct_of_total     = round(c(nrow(decreased), nrow(increased),
                             nrow(changed)) / nrow(delta_ffi) * 100, 1),
  mean_FFI_current = round(c(mean(decreased$FFI_current, na.rm = TRUE),
                             mean(increased$FFI_current, na.rm = TRUE),
                             mean(changed$FFI_current,   na.rm = TRUE)), 4),
  mean_delta_FFI   = round(c(mean(decreased$delta_FFI, na.rm = TRUE),
                             mean(increased$delta_FFI, na.rm = TRUE),
                             mean(changed$delta_FFI,   na.rm = TRUE)), 4),
  pct_change_FFI   = round(c(
    mean(decreased$delta_FFI, na.rm = TRUE) /
      mean(decreased$FFI_current, na.rm = TRUE) * 100,
    mean(increased$delta_FFI, na.rm = TRUE) /
      mean(increased$FFI_current, na.rm = TRUE) * 100,
    mean(changed$delta_FFI, na.rm = TRUE) /
      mean(changed$FFI_current, na.rm = TRUE) * 100), 2),
  mean_delta_ED    = round(c(mean(decreased$delta_ED,  na.rm = TRUE),
                             mean(increased$delta_ED,  na.rm = TRUE),
                             mean(changed$delta_ED,    na.rm = TRUE)), 4),
  mean_delta_PD    = round(c(mean(decreased$delta_PD,  na.rm = TRUE),
                             mean(increased$delta_PD,  na.rm = TRUE),
                             mean(changed$delta_PD,    na.rm = TRUE)), 4),
  mean_delta_MPA   = round(c(mean(decreased$delta_MPA, na.rm = TRUE),
                             mean(increased$delta_MPA, na.rm = TRUE),
                             mean(changed$delta_MPA,   na.rm = TRUE)), 4)
)

write.csv(changed_summary,
          file.path(output_dir,
                    "summary_changed_grids_holistic_hotspot_hh_b.csv"),
          row.names = FALSE)
cat("\nSaved: summary_changed_grids_holistic_hotspot_hh_b.csv\n")


####### A clean tidy up version of getting the stats
####### Test before commit
# ============================================================
# Holistic hotspot FFI results summary
# ------------------------------------------------------------
# Produces ONE tidy CSV: one row per metric (FFI, ED, PD, MPA)
# reporting, side by side:
#   - all grids           : the global average (diluted by the
#                           many cells the scenario never touched)
#   - changed footprint   : grids where FFI actually changed,
#                           i.e. where holistic hotspot did something
#   - improved only       : within the footprint, only the cells
#                           moving in the beneficial direction
#                           (so you can omit cells moving the wrong way)
#
# Improvement direction per metric:
#   FFI, ED, PD improve when they DECREASE
#   MPA (mean patch area) improves when it INCREASES
# ============================================================

library(dplyr)

# ---------------- PATHS ----------------
read_dir       <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results"
output_dir     <- file.path(read_dir, "holistic_hotspot_results")
scenario_label <- "holistic_hotspot_hh_b"
cell_area      <- 25          # 5km x 5km grid = 25 km2

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------- LOAD ----------------
cat("Loading delta FFI results...\n")
delta_ffi <- read.csv(file.path(read_dir,
                                "delta_FFI_holistic_hotspot_hh_b_10m.csv"))
cat("Rows loaded:", format(nrow(delta_ffi), big.mark = ","), "\n\n")

# ---------------- FOOTPRINT (grids the scenario changed) ----------------
delta_ffi <- delta_ffi %>%
  mutate(
    delta_class = case_when(
      delta_FFI < -0.001 ~ -1L,
      delta_FFI >  0.001 ~  1L,
      TRUE               ~  0L
    )
  )

n_total  <- nrow(delta_ffi)
footprint <- delta_ffi$delta_class != 0    # logical: cell was altered
n_changed <- sum(footprint, na.rm = TRUE)

# ---------------- ONE ROW PER METRIC ----------------
# improve = "decrease" -> beneficial when delta < 0
# improve = "increase" -> beneficial when delta > 0
metrics_info <- list(
  list(name = "FFI", cur = "FFI_current", del = "delta_FFI", improve = "decrease"),
  list(name = "ED",  cur = "ED_current",  del = "delta_ED",  improve = "decrease"),
  list(name = "PD",  cur = "PD_current",  del = "delta_PD",  improve = "decrease"),
  list(name = "MPA", cur = "MPA_current", del = "delta_MPA", improve = "increase")
)

pct <- function(delta_mean, current_mean) {
  if (is.na(current_mean) || current_mean == 0) return(NA_real_)
  round(delta_mean / current_mean * 100, 2)
}

build_row <- function(m) {
  cur_all <- delta_ffi[[m$cur]]
  del_all <- delta_ffi[[m$del]]
  
  # within the FFI-changed footprint
  cur_chg <- cur_all[footprint]
  del_chg <- del_all[footprint]
  
  # improved subset within the footprint
  imp <- if (m$improve == "decrease") del_chg < 0 else del_chg > 0
  imp[is.na(imp)] <- FALSE
  cur_imp <- cur_chg[imp]
  del_imp <- del_chg[imp]
  
  data.frame(
    scenario              = scenario_label,
    metric                = m$name,
    improve_direction     = m$improve,
    
    # footprint counts
    n_grids_total         = n_total,
    n_grids_changed       = n_changed,
    pct_grids_changed     = round(n_changed / n_total * 100, 1),
    changed_area_km2      = n_changed * cell_area,
    
    # scope 1: all grids
    mean_current_all      = round(mean(cur_all, na.rm = TRUE), 4),
    mean_scenario_all     = round(mean(cur_all + del_all, na.rm = TRUE), 4),
    mean_delta_all        = round(mean(del_all, na.rm = TRUE), 4),
    pct_change_all        = pct(mean(del_all, na.rm = TRUE),
                                mean(cur_all, na.rm = TRUE)),
    
    # scope 2: changed footprint
    mean_current_changed  = round(mean(cur_chg, na.rm = TRUE), 4),
    mean_scenario_changed = round(mean(cur_chg + del_chg, na.rm = TRUE), 4),
    mean_delta_changed    = round(mean(del_chg, na.rm = TRUE), 4),
    pct_change_changed    = pct(mean(del_chg, na.rm = TRUE),
                                mean(cur_chg, na.rm = TRUE)),
    
    # scope 3: improved only (beneficial direction within footprint)
    n_grids_improved      = sum(imp),
    pct_grids_improved    = round(sum(imp) / n_changed * 100, 1),
    mean_delta_improved   = round(mean(del_imp, na.rm = TRUE), 4),
    pct_change_improved   = pct(mean(del_imp, na.rm = TRUE),
                                mean(cur_imp, na.rm = TRUE)),
    
    stringsAsFactors = FALSE
  )
}

summary_tbl <- do.call(rbind, lapply(metrics_info, build_row))

out_csv <- file.path(output_dir,
                     paste0("summary_", scenario_label, ".csv"))
write.csv(summary_tbl, out_csv, row.names = FALSE)

# ---------------- CONSOLE READOUT ----------------
ffi <- summary_tbl[summary_tbl$metric == "FFI", ]

cat("============================================================\n")
cat("  HOLISTIC HOTSPOT SUMMARY\n")
cat("============================================================\n")
cat(sprintf("  Total grids:            %s\n", format(n_total, big.mark = ",")))
cat(sprintf("  Grids changed:          %s  (%.1f%% of total, %s km2)\n",
            format(n_changed, big.mark = ","),
            ffi$pct_grids_changed,
            format(ffi$changed_area_km2, big.mark = ",")))
cat("\n  Within changed grids:\n")
cat(sprintf("    Mean FFI current:     %.4f\n", ffi$mean_current_changed))
cat(sprintf("    Mean FFI scenario:    %.4f\n", ffi$mean_scenario_changed))
cat(sprintf("    Mean delta FFI:       %+.4f  (%+.2f%% vs baseline)\n",
            ffi$mean_delta_changed, ffi$pct_change_changed))
cat(sprintf("    Grids improved:       %.1f%%\n", ffi$pct_grids_improved))

cat("\n  Per metric within changed grids (pct vs baseline):\n")
for (i in seq_len(nrow(summary_tbl))) {
  r <- summary_tbl[i, ]
  cat(sprintf("    %-4s  delta %+.4f  (%+.2f%%)   improved cells %+.2f%%\n",
              r$metric, r$mean_delta_changed, r$pct_change_changed,
              r$pct_change_improved))
}

cat("\n------------------------------------------------------------\n")
cat("SUGGESTED REPORTING (fill straight from the numbers above):\n")
cat(sprintf(
  "  The holistic hotspot scenario altered %.1f%% of tropical forest\n  grid cells (%s of %s, %s km2). Within these affected cells,\n  mean FFI fell from %.4f to %.4f, a reduction of %.1f%% relative\n  to the current baseline, and %.1f%% of affected cells showed\n  reduced fragmentation.\n",
  ffi$pct_grids_changed,
  format(n_changed, big.mark = ","),
  format(n_total,   big.mark = ","),
  format(ffi$changed_area_km2, big.mark = ","),
  ffi$mean_current_changed,
  ffi$mean_scenario_changed,
  abs(ffi$pct_change_changed),
  ffi$pct_grids_improved))
cat("------------------------------------------------------------\n")
cat("\nSaved:", out_csv, "\n")



























