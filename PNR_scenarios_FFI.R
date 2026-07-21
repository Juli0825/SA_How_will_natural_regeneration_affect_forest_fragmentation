library(data.table)
library(dplyr)
library(tidyr)
library(terra)

# ====================================================
# CHANGE THESE THREE LINES FOR EACH SCENARIO
# ====================================================

scenario_name        <- "all_pnr"
scenario_metrics_dir <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FutureScenario_metrics/all_pnr_metrics"
current_ffi_path     <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/current_forest_FFI_10m.csv"

# ====================================================
# FIXED PATHS
# ====================================================

boundaries_path <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/global_boundaries_10m_current_forest.csv"
output_dir      <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/ap_fix_imcomplete"
mollweide_crs   <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
round_digits    <- 10

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#snap_5km <- function(x) round(x / 5000) * 5000 # this will cause stripes at West Africa and some other tiles
                                                # When those residual 0 points are written to a raster that expects centres at 2500, 
                                                # every second column falls into a gap. That is your barcode.
#snap_5km <- function(x) round(x) # this line fix the stripe problem
                                  # nope!hh is fine, but delta ffi ap is domed by half of plots missing

#snap_5km <- function(x) floor(x / 5000) * 5000 + 2500 # misting parts, still not working

snap_5km <- function(x) floor(x / 5000 + 0.5) * 5000
# ====================================================
# STEP 1: LOAD GLOBAL BOUNDARIES
# ====================================================

cat("===== STEP 1: Loading global boundaries =====\n")

boundaries <- read.csv(boundaries_path)
print(boundaries)

get_b <- function(m, type) {
  boundaries$value[boundaries$metric == m &
                     boundaries$boundary_type == type]
}

ed_min  <- get_b("ed",      "lower");  ed_max  <- get_b("ed",      "upper")
pd_min  <- get_b("pd",      "lower");  pd_max  <- get_b("pd",      "upper")
mpa_min <- get_b("area_mn", "lower");  mpa_max <- get_b("area_mn", "upper")

cat("ED  boundaries:", ed_min,  "to", ed_max,  "\n")
cat("PD  boundaries:", pd_min,  "to", pd_max,  "\n")
cat("MPA boundaries:", mpa_min, "to", mpa_max, "\n")

# ====================================================
# STEP 2: LOAD SCENARIO METRICS
# ====================================================

cat("\n===== STEP 2: Loading", scenario_name, "metrics =====\n")

scenario_files <- list.files(scenario_metrics_dir,
                             pattern = "_ap_b\\.csv$",
                             full.names = TRUE)
cat("CSV files found:", length(scenario_files), "\n")

if (length(scenario_files) == 0) {
  stop("No metrics files found — check scenario_metrics_dir")
}

scenario_all <- rbindlist(lapply(scenario_files, fread), fill = TRUE)
cat("Total rows loaded:", format(nrow(scenario_all), big.mark = ","), "\n")
cat("Tiles present:    ", length(unique(scenario_all$tile_id)), "\n")

scenario_forest <- scenario_all[class == 1 &
                                  metric %in% c("ed", "pd", "area_mn")]
cat("Rows after filter:", format(nrow(scenario_forest), big.mark = ","), "\n")

rm(scenario_all)
gc()

# ====================================================
# STEP 3: CALCULATE SCENARIO FFI
# ====================================================

cat("\n===== STEP 3: Calculating", scenario_name, "FFI =====\n")

scenario_ffi <- scenario_forest %>%
  as_tibble() %>%
  select(plot_id, tile_id, center_x, center_y, metric, value) %>%
  mutate(
    center_x = snap_5km(center_x),
    center_y = snap_5km(center_y)
  ) %>%
  pivot_wider(names_from  = metric,
              values_from = value) %>%
  rename(ED = ed, PD = pd, MPA = area_mn) %>%
  mutate(
    ED_capped  = pmax(ed_min,  pmin(ed_max,  ED)),
    PD_capped  = pmax(pd_min,  pmin(pd_max,  PD)),
    MPA_capped = pmax(mpa_min, pmin(mpa_max, MPA))
  ) %>%
  mutate(
    ED_norm  = (ED_capped  - ed_min)  / (ed_max  - ed_min),
    PD_norm  = (PD_capped  - pd_min)  / (pd_max  - pd_min),
    MPA_norm = (MPA_capped - mpa_min) / (mpa_max - mpa_min)
  ) %>%
  mutate(FFI = (ED_norm + PD_norm + (1 - MPA_norm)) / 3) %>%
  mutate(
    ED_norm  = round(ED_norm,  round_digits),
    PD_norm  = round(PD_norm,  round_digits),
    MPA_norm = round(MPA_norm, round_digits),
    FFI      = round(FFI,      round_digits)
  )

cat("FFI range:", round(range(scenario_ffi$FFI, na.rm = TRUE), 4), "\n")
cat("FFI mean: ", round(mean(scenario_ffi$FFI,  na.rm = TRUE), 4), "\n")
cat("Total plots:", format(nrow(scenario_ffi), big.mark = ","), "\n")

scenario_ffi_path <- file.path(output_dir,
                               paste0(scenario_name, "_ap_b_FFI_10m.csv"))
write.csv(scenario_ffi, scenario_ffi_path, row.names = FALSE)
cat("Saved:", scenario_ffi_path, "\n")

# ====================================================
# STEP 4: DELTA FFI
# ====================================================

cat("\n===== STEP 4: Calculating delta FFI =====\n")

current_ffi <- read.csv(current_ffi_path) %>%
  mutate(
    center_x = snap_5km(center_x),
    center_y = snap_5km(center_y)
  )

cat("Current forest plots:", format(nrow(current_ffi), big.mark = ","), "\n")

delta_ffi <- current_ffi %>%
  select(plot_id, tile_id, center_x, center_y,
         FFI_current = FFI,
         ED_current  = ED,
         PD_current  = PD,
         MPA_current = MPA) %>%
  inner_join(
    scenario_ffi %>%
      select(center_x, center_y, tile_id,
             FFI_scenario = FFI,
             ED_scenario  = ED,
             PD_scenario  = PD,
             MPA_scenario = MPA),
    by = c("center_x", "center_y", "tile_id")
  ) %>%
  mutate(
    delta_FFI = round(FFI_scenario - FFI_current, round_digits),
    delta_ED  = round(ED_scenario  - ED_current,  round_digits),
    delta_PD  = round(PD_scenario  - PD_current,  round_digits),
    delta_MPA = round(MPA_scenario - MPA_current, round_digits)
  )

n_exact_zero <- sum(delta_ffi$delta_FFI == 0, na.rm = TRUE)
n_tiny       <- sum(abs(delta_ffi$delta_FFI) > 0 &
                      abs(delta_ffi$delta_FFI) < 1e-9, na.rm = TRUE)

cat("\nFloating point check:\n")
cat("  Exactly zero:       ",
    format(n_exact_zero, big.mark = ","), "\n")
cat("  Tiny (<1e-9):       ",
    format(n_tiny, big.mark = ","), "(should be 0)\n\n")

cat("Plots in current FFI:  ", format(nrow(current_ffi),  big.mark = ","), "\n")
cat("Plots in scenario FFI: ", format(nrow(scenario_ffi), big.mark = ","), "\n")
cat("Matched plots:         ", format(nrow(delta_ffi),    big.mark = ","), "\n")
cat("Unmatched:             ",
    format(nrow(current_ffi) - nrow(delta_ffi), big.mark = ","), "\n")

cat("\ndelta FFI summary:\n")
cat("Mean:   ", round(mean(delta_ffi$delta_FFI,  na.rm = TRUE), 4), "\n")
cat("Range:  ", round(range(delta_ffi$delta_FFI, na.rm = TRUE), 4), "\n")

delta_path <- file.path(output_dir,
                        paste0("delta_FFI_", scenario_name, "_ap_b_10m.csv"))
write.csv(delta_ffi, delta_path, row.names = FALSE)
cat("\nSaved:", delta_path, "\n")

# ====================================================
# STEP 5: CLASSIFY DELTA FFI
# No threshold — use strict zero for unchanged
# round_digits = 10 already eliminates floating point
# noise so any non-zero value is a genuine signal
# ====================================================

cat("\n===== STEP 5: Classifying delta FFI =====\n")

delta_ffi <- delta_ffi %>%
  mutate(
    delta_class = case_when(
      delta_FFI < 0 ~ -1L,
      delta_FFI > 0 ~  1L,
      TRUE          ~  0L
    ),
    delta_label = case_when(
      delta_class == -1 ~ "FFI decreased",
      delta_class ==  0 ~ "unchanged",
      delta_class ==  1 ~ "FFI increased"
    )
  )

cat("FFI decreased: ",
    format(sum(delta_ffi$delta_class == -1), big.mark = ","), "\n")
cat("Unchanged:     ",
    format(sum(delta_ffi$delta_class ==  0), big.mark = ","), "\n")
cat("FFI increased: ",
    format(sum(delta_ffi$delta_class ==  1), big.mark = ","), "\n")

# ====================================================
# STEP 6: RASTERIZE
# ====================================================

cat("\n===== STEP 6: Rasterizing =====\n")

make_raster <- function(df, field,
                        crs = mollweide_crs, res = 5000) {
  pts <- vect(df, geom = c("center_x", "center_y"), crs = crs)
  tpl <- rast(ext = ext(pts), resolution = res, crs = crs)
  rasterize(pts, tpl, field = field, fun = "mean")
}

cat("1. Current forest FFI...\n")
r_current_ffi  <- make_raster(current_ffi,  "FFI")

cat("2. Scenario FFI...\n")
r_scenario_ffi <- make_raster(scenario_ffi, "FFI")

cat("3. Delta FFI continuous...\n")
r_delta_cont   <- make_raster(delta_ffi, "delta_FFI")

cat("4. Delta FFI classified...\n")
r_delta_class  <- make_raster(delta_ffi, "delta_class")

cat("5. Delta ED...\n")
r_delta_ed     <- make_raster(delta_ffi, "delta_ED")

cat("6. Delta PD...\n")
r_delta_pd     <- make_raster(delta_ffi, "delta_PD")

cat("7. Delta MPA...\n")
r_delta_mpa    <- make_raster(delta_ffi, "delta_MPA")

# ====================================================
# STEP 7: SAVE RASTERS
# ====================================================

cat("\n===== STEP 7: Saving rasters =====\n")

writeRaster(r_current_ffi,
            file.path(output_dir, "current_forest_FFI_raster.tif"),
            overwrite = TRUE)
writeRaster(r_scenario_ffi,
            file.path(output_dir,
                      paste0(scenario_name, "_ap_b_FFI_raster.tif")),
            overwrite = TRUE)
writeRaster(r_delta_cont,
            file.path(output_dir,
                      paste0("delta_FFI_", scenario_name,
                             "_ap_b_continuous.tif")),
            overwrite = TRUE)
writeRaster(r_delta_class,
            file.path(output_dir,
                      paste0("delta_FFI_", scenario_name,
                             "_ap_b_classified.tif")),
            overwrite = TRUE)
writeRaster(r_delta_ed,
            file.path(output_dir,
                      paste0("delta_ED_", scenario_name, "_ap_b.tif")),
            overwrite = TRUE)
writeRaster(r_delta_pd,
            file.path(output_dir,
                      paste0("delta_PD_", scenario_name, "_ap_b.tif")),
            overwrite = TRUE)
writeRaster(r_delta_mpa,
            file.path(output_dir,
                      paste0("delta_MPA_", scenario_name, "_ap_b.tif")),
            overwrite = TRUE)

cat("All 7 rasters saved\n")

# ====================================================
# STEP 8: METRIC CHANGE SUMMARY
# ====================================================

cat("\n===== STEP 8: Metric change summary =====\n")

mean_FFI_current <- mean(delta_ffi$FFI_current, na.rm = TRUE)
mean_ED_current  <- mean(delta_ffi$ED_current,  na.rm = TRUE)
mean_PD_current  <- mean(delta_ffi$PD_current,  na.rm = TRUE)
mean_MPA_current <- mean(delta_ffi$MPA_current, na.rm = TRUE)

mean_delta_FFI   <- mean(delta_ffi$delta_FFI, na.rm = TRUE)
mean_delta_ED    <- mean(delta_ffi$delta_ED,  na.rm = TRUE)
mean_delta_PD    <- mean(delta_ffi$delta_PD,  na.rm = TRUE)
mean_delta_MPA   <- mean(delta_ffi$delta_MPA, na.rm = TRUE)

cat(sprintf("\n  %-22s  current: %8.4f  delta: %+.4f  (%+.2f%%)\n",
            "FFI",
            mean_FFI_current, mean_delta_FFI,
            mean_delta_FFI / mean_FFI_current * 100))
cat(sprintf("  %-22s  current: %8.4f  delta: %+.4f  (%+.2f%%)\n",
            "Edge density (ED)",
            mean_ED_current, mean_delta_ED,
            mean_delta_ED / mean_ED_current * 100))
cat(sprintf("  %-22s  current: %8.4f  delta: %+.4f  (%+.2f%%)\n",
            "Patch density (PD)",
            mean_PD_current, mean_delta_PD,
            mean_delta_PD / mean_PD_current * 100))
cat(sprintf("  %-22s  current: %8.2f  delta: %+.4f  (%+.2f%%)\n",
            "Mean patch area (MPA)",
            mean_MPA_current, mean_delta_MPA,
            mean_delta_MPA / mean_MPA_current * 100))

cat("\n  ED:  negative delta = less edge = less fragmented\n")
cat("  PD:  negative delta = fewer patches = less fragmented\n")
cat("  MPA: positive delta = larger patches = less fragmented\n")

cat("\n===== ALL DONE =====\n")
cat("Output files saved to:", output_dir, "\n")