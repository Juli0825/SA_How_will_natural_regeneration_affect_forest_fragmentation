### Each PNR future forest scenario FFI
### Using global boundaries created from current_forest_10m

library(data.table)
library(dplyr)
library(tidyr)

# ====================================================
# Just for thesis/ holistic hotspot + all pnr
# ====================================================

scenario_name        <- "holistic_hotspot"
scenario_metrics_dir <- "2026_NEE_R2/FutureScenario_metrics/holistic_hotspot"
current_ffi_path     <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/current_forest_FFI_10m.csv"

# ====================================================
# FIXED PATHS - do not change these
# ====================================================

boundaries_path <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/global_boundaries_10m_current_forest.csv"
output_dir      <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ====================================================
# STEP 1: LOAD GLOBAL BOUNDARIES FROM LOCAL DRIVE
# Calculated once from current forest - never changes
# ====================================================

cat("===== STEP 1: Loading global boundaries =====\n")

boundaries <- read.csv(boundaries_path)

cat("Boundaries loaded from:\n", boundaries_path, "\n\n")
print(boundaries)

get_b <- function(m, type) {
  boundaries$value[boundaries$metric == m &
                     boundaries$boundary_type == type]
}

ed_min  <- get_b("ed",      "lower");  ed_max  <- get_b("ed",      "upper")
pd_min  <- get_b("pd",      "lower");  pd_max  <- get_b("pd",      "upper")
mpa_min <- get_b("area_mn", "lower");  mpa_max <- get_b("area_mn", "upper")

cat("\nED  boundaries:", ed_min,  "to", ed_max,  "\n")
cat("PD  boundaries:", pd_min,  "to", pd_max,  "\n")
cat("MPA boundaries:", mpa_min, "to", mpa_max, "\n")

# ====================================================
# STEP 2: LOAD SCENARIO METRICS
# ====================================================

cat("\n===== STEP 2: Loading", scenario_name, "metrics =====\n")

scenario_files <- list.files(scenario_metrics_dir,
                             pattern = paste0("_", scenario_name, "_metrics\\.csv$"),
                             full.names = TRUE)
cat("CSV files found:", length(scenario_files), "\n")

if (length(scenario_files) == 0) {
  stop("No metrics files found - check scenario_metrics_dir and scenario_name")
}

scenario_all <- rbindlist(lapply(scenario_files, fread), fill = TRUE)
cat("Total rows loaded:", format(nrow(scenario_all), big.mark = ","), "\n")
cat("Tiles present:    ", length(unique(scenario_all$tile_id)), "\n")

scenario_forest <- scenario_all[class == 1 & metric %in% c("ed", "pd", "area_mn")]
cat("Rows after filter:", format(nrow(scenario_forest), big.mark = ","), "\n")

rm(scenario_all)
gc()

# ====================================================
# STEP 3: CALCULATE SCENARIO FFI
# Same boundaries as current forest - directly comparable
# ====================================================

cat("\n===== STEP 3: Calculating", scenario_name, "FFI =====\n")

scenario_ffi <- scenario_forest %>%
  as_tibble() %>%
  select(plot_id, tile_id, center_x, center_y, metric, value) %>%
  pivot_wider(names_from  = metric,
              values_from = value) %>%
  rename(ED = ed, PD = pd, MPA = area_mn) %>%
  
  # Cap at boundaries BEFORE normalizing
  mutate(
    ED_capped  = pmax(ed_min,  pmin(ed_max,  ED)),
    PD_capped  = pmax(pd_min,  pmin(pd_max,  PD)),
    MPA_capped = pmax(mpa_min, pmin(mpa_max, MPA))
  ) %>%
  
  # Normalize using fixed boundary range
  mutate(
    ED_norm  = (ED_capped  - ed_min)  / (ed_max  - ed_min),
    PD_norm  = (PD_capped  - pd_min)  / (pd_max  - pd_min),
    MPA_norm = (MPA_capped - mpa_min) / (mpa_max - mpa_min)
  ) %>%
  
  # FFI = (ED_norm + PD_norm + (1 - MPA_norm)) / 3
  mutate(FFI = (ED_norm + PD_norm + (1 - MPA_norm)) / 3)

cat("FFI range:", round(range(scenario_ffi$FFI, na.rm = TRUE), 4), "\n")
cat("FFI mean: ", round(mean(scenario_ffi$FFI,  na.rm = TRUE), 4), "\n")
cat("Total plots:", format(nrow(scenario_ffi), big.mark = ","), "\n")

scenario_ffi_path <- file.path(output_dir,
                               paste0(scenario_name, "_FFI_10m.csv"))
write.csv(scenario_ffi, scenario_ffi_path, row.names = FALSE)
cat("Saved:", scenario_ffi_path, "\n")

# ====================================================
# STEP 4: LOAD CURRENT FOREST FFI AND CALCULATE
#         DELTA FFI
#         negative = fragmentation improved
#         positive = fragmentation worsened
# ====================================================

cat("\n===== STEP 4: Calculating delta FFI =====\n")

cat("Loading current forest FFI from:\n", current_ffi_path, "\n")
current_ffi <- read.csv(current_ffi_path)
cat("Current forest plots:", format(nrow(current_ffi), big.mark = ","), "\n")

delta_ffi <- current_ffi %>%
  select(plot_id, tile_id, center_x, center_y,
         FFI_current = FFI,
         ED_current  = ED,
         PD_current  = PD,
         MPA_current = MPA) %>%
  inner_join(
    scenario_ffi %>%
      select(plot_id, tile_id,
             FFI_scenario = FFI,
             ED_scenario  = ED,
             PD_scenario  = PD,
             MPA_scenario = MPA),
    by = c("plot_id", "tile_id")
  ) %>%
  mutate(
    delta_FFI = FFI_scenario - FFI_current,
    delta_ED  = ED_scenario  - ED_current,
    delta_PD  = PD_scenario  - PD_current,
    delta_MPA = MPA_scenario - MPA_current
  )

cat("\nPlots in current FFI:  ", format(nrow(current_ffi),  big.mark = ","), "\n")
cat("Plots in scenario FFI: ", format(nrow(scenario_ffi), big.mark = ","), "\n")
cat("Matched plots:         ", format(nrow(delta_ffi),    big.mark = ","), "\n")
cat("Unmatched (no-PNR):   ",
    format(nrow(current_ffi) - nrow(delta_ffi), big.mark = ","), "\n")

cat("\ndelta FFI summary:\n")
cat("Mean:                  ", round(mean(delta_ffi$delta_FFI,  na.rm = TRUE), 4), "\n")
cat("Range:                 ", round(range(delta_ffi$delta_FFI, na.rm = TRUE), 4), "\n")
cat("Plots improved (< 0): ", format(sum(delta_ffi$delta_FFI < 0,  na.rm = TRUE), big.mark = ","), "\n")
cat("Plots worsened (> 0): ", format(sum(delta_ffi$delta_FFI > 0,  na.rm = TRUE), big.mark = ","), "\n")
cat("Plots unchanged (= 0):", format(sum(delta_ffi$delta_FFI == 0, na.rm = TRUE), big.mark = ","), "\n")

delta_path <- file.path(output_dir,
                        paste0("delta_FFI_", scenario_name, "_10m.csv"))
write.csv(delta_ffi, delta_path, row.names = FALSE)
cat("\nSaved:", delta_path, "\n")

cat("\n===== ALL DONE =====\n")
cat("Output files:\n")
cat("1.", scenario_ffi_path, "\n")
cat("2.", delta_path, "\n")

#####
##### Stats
library(dplyr)

# ====================================================
# EXPLORE DELTA FFI DISTRIBUTION
# ====================================================

cat("===== Delta FFI distribution =====\n\n")

cat("--- Basic statistics ---\n")
cat(sprintf("  Mean:    %+.4f\n", mean(delta_ffi$delta_FFI, na.rm = TRUE)))
cat(sprintf("  Median:  %+.4f\n", median(delta_ffi$delta_FFI, na.rm = TRUE)))
cat(sprintf("  SD:      %.4f\n",  sd(delta_ffi$delta_FFI, na.rm = TRUE)))
cat(sprintf("  Min:     %+.4f\n", min(delta_ffi$delta_FFI, na.rm = TRUE)))
cat(sprintf("  Max:     %+.4f\n", max(delta_ffi$delta_FFI, na.rm = TRUE)))

cat("\n--- Quantiles ---\n")
q <- quantile(delta_ffi$delta_FFI, 
              probs = c(0.01, 0.05, 0.10, 0.25, 
                        0.50, 0.75, 0.90, 0.95, 0.99),
              na.rm = TRUE)
for (i in seq_along(q)) {
  cat(sprintf("  %3.0f%%:  %+.4f\n", as.numeric(names(q)[i]) * 100, q[i]))
}

cat("\n--- Distribution buckets ---\n")
breaks <- c(-Inf, -0.1, -0.05, -0.01, -0.001,
            0.001, 0.01, 0.05, 0.1, Inf)
labels <- c("< -0.10", "-0.10 to -0.05", "-0.05 to -0.01",
            "-0.01 to -0.001", "-0.001 to +0.001",
            "+0.001 to +0.01", "+0.01 to +0.05",
            "+0.05 to +0.10", "> +0.10")

delta_ffi$bucket <- cut(delta_ffi$delta_FFI,
                        breaks = breaks,
                        labels = labels)

bucket_summary <- delta_ffi %>%
  group_by(bucket) %>%
  summarise(
    n_cells  = n(),
    area_km2 = n() * 25,
    pct      = n() / nrow(delta_ffi) * 100,
    .groups  = "drop"
  )

for (i in 1:nrow(bucket_summary)) {
  cat(sprintf("  %-22s  %7s cells  %8.0f km²  %5.1f%%\n",
              bucket_summary$bucket[i],
              format(bucket_summary$n_cells[i], big.mark = ","),
              bucket_summary$area_km2[i],
              bucket_summary$pct[i]))
}

cat("\n--- Zero and near-zero ---\n")
cat(sprintf("  Exactly zero:       %s cells\n",
            format(sum(delta_ffi$delta_FFI == 0, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  Within ±0.001:      %s cells\n",
            format(sum(abs(delta_ffi$delta_FFI) < 0.001, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  Within ±0.01:       %s cells\n",
            format(sum(abs(delta_ffi$delta_FFI) < 0.01, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  Within ±0.05:       %s cells\n",
            format(sum(abs(delta_ffi$delta_FFI) < 0.05, na.rm = TRUE), big.mark = ",")))

#####


# Rasterise delta FFI and get stats
#####
library(terra)
library(dplyr)

scenario_name <- "holistic_hotspot"
output_dir    <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results"
cell_area_km2 <- 25

# ====================================================
# STEP 1: CLASSIFY DELTA FFI
# < -0.001 = FFI decreased
# -0.001 to +0.001 = unchanged (floating point buffer)
# > +0.001 = FFI increased
# ====================================================

cat("===== STEP 1: Classifying delta FFI =====\n")

delta_ffi <- delta_ffi %>%
  mutate(
    delta_class = case_when(
      delta_FFI < -0.001 ~ -1,
      delta_FFI >  0.001 ~  1,
      TRUE               ~  0
    ),
    delta_label = case_when(
      delta_class == -1 ~ "FFI decreased",
      delta_class ==  0 ~ "unchanged",
      delta_class ==  1 ~ "FFI increased"
    )
  )

cat("Classification counts:\n")
cat(sprintf("  FFI decreased: %s cells\n",
            format(sum(delta_ffi$delta_class == -1), big.mark = ",")))
cat(sprintf("  Unchanged:     %s cells\n",
            format(sum(delta_ffi$delta_class ==  0), big.mark = ",")))
cat(sprintf("  FFI increased: %s cells\n",
            format(sum(delta_ffi$delta_class ==  1), big.mark = ",")))

# ====================================================
# STEP 2: RASTERIZE
# ====================================================

cat("\n===== STEP 2: Rasterizing =====\n")

mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

delta_points <- vect(
  delta_ffi,
  geom = c("center_x", "center_y"),
  crs  = mollweide_crs
)

template <- rast(
  ext        = ext(delta_points),
  resolution = 5000,
  crs        = mollweide_crs
)

# Continuous delta FFI raster
delta_rast <- rasterize(delta_points, template,
                        field = "delta_FFI",   fun = "mean")

# Classified raster (-1 / 0 / 1)
class_rast <- rasterize(delta_points, template,
                        field = "delta_class", fun = "mean")

# Save
delta_rast_path <- file.path(output_dir,
                             paste0("delta_FFI_", scenario_name, "_continuous.tif"))
class_rast_path <- file.path(output_dir,
                             paste0("delta_FFI_", scenario_name, "_classified.tif"))

writeRaster(delta_rast, delta_rast_path, overwrite = TRUE)
writeRaster(class_rast, class_rast_path, overwrite = TRUE)

cat("Continuous raster saved to:\n ", delta_rast_path, "\n")
cat("Classified raster saved to:\n ", class_rast_path, "\n")
cat("  -1 = FFI decreased\n")
cat("   0 = unchanged\n")
cat("  +1 = FFI increased\n")

# ====================================================
# STEP 3: SUMMARY STATISTICS
# ====================================================

cat("\n===== STEP 3: Summary statistics =====\n")

n_total_current <- nrow(current_ffi)
n_matched       <- nrow(delta_ffi)
n_unmatched     <- n_total_current - n_matched

n_decreased     <- sum(delta_ffi$delta_class == -1, na.rm = TRUE)
n_unchanged     <- sum(delta_ffi$delta_class ==  0, na.rm = TRUE)
n_increased     <- sum(delta_ffi$delta_class ==  1, na.rm = TRUE)

# FFI change
mean_ffi_current  <- mean(delta_ffi$FFI_current,  na.rm = TRUE)
mean_ffi_scenario <- mean(delta_ffi$FFI_scenario, na.rm = TRUE)
mean_delta_ffi    <- mean(delta_ffi$delta_FFI,    na.rm = TRUE)
pct_change_ffi    <- (mean_delta_ffi / mean_ffi_current) * 100

# Individual metric changes
mean_ED_current  <- mean(delta_ffi$ED_current,  na.rm = TRUE)
mean_PD_current  <- mean(delta_ffi$PD_current,  na.rm = TRUE)
mean_MPA_current <- mean(delta_ffi$MPA_current, na.rm = TRUE)

mean_delta_ED    <- mean(delta_ffi$delta_ED,  na.rm = TRUE)
mean_delta_PD    <- mean(delta_ffi$delta_PD,  na.rm = TRUE)
mean_delta_MPA   <- mean(delta_ffi$delta_MPA, na.rm = TRUE)

pct_change_ED    <- (mean_delta_ED  / mean_ED_current)  * 100
pct_change_PD    <- (mean_delta_PD  / mean_PD_current)  * 100
pct_change_MPA   <- (mean_delta_MPA / mean_MPA_current) * 100

# ====================================================
# STEP 4: PRINT RESULTS
# ====================================================

cat("\n")
cat("============================================================\n")
cat(sprintf("  SCENARIO: %s\n", toupper(gsub("_", " ", scenario_name))))
cat("============================================================\n\n")

cat("--- OVERALL FFI CHANGE ---\n")
cat(sprintf("  Mean FFI current:    %.4f\n",   mean_ffi_current))
cat(sprintf("  Mean FFI scenario:   %.4f\n",   mean_ffi_scenario))
cat(sprintf("  Mean delta FFI:      %+.4f  (%+.2f%%)\n",
            mean_delta_ffi, pct_change_ffi))
cat(sprintf("  Median delta FFI:    %+.4f\n",
            median(delta_ffi$delta_FFI, na.rm = TRUE)))
cat(sprintf("  Range delta FFI:     %+.4f  to  %+.4f\n",
            min(delta_ffi$delta_FFI, na.rm = TRUE),
            max(delta_ffi$delta_FFI, na.rm = TRUE)))

cat("\n--- SPATIAL BREAKDOWN ---\n")
cat(sprintf("  %-16s  %8s cells  %10.0f km²  %5.1f%%\n",
            "FFI decreased",
            format(n_decreased, big.mark = ","),
            n_decreased * cell_area_km2,
            n_decreased / n_total_current * 100))
cat(sprintf("  %-16s  %8s cells  %10.0f km²  %5.1f%%\n",
            "Unchanged",
            format(n_unchanged + n_unmatched, big.mark = ","),
            (n_unchanged + n_unmatched) * cell_area_km2,
            (n_unchanged + n_unmatched) / n_total_current * 100))
cat(sprintf("  %-16s  %8s cells  %10.0f km²  %5.1f%%\n",
            "FFI increased",
            format(n_increased, big.mark = ","),
            n_increased * cell_area_km2,
            n_increased / n_total_current * 100))
cat(sprintf("  %-16s  %8s cells  %10.0f km²  %5.1f%%\n",
            "TOTAL",
            format(n_total_current, big.mark = ","),
            n_total_current * cell_area_km2,
            100.0))
cat(sprintf("\n  Note: unchanged includes %s no-PNR tiles\n",
            format(n_unmatched, big.mark = ",")))
cat("        identical to current forest\n")

cat("\n--- INDIVIDUAL METRIC CHANGES ---\n")
cat("  (mean across matched plots)\n\n")
cat(sprintf("  %-22s  current: %8.4f  delta: %+.4f  (%+.2f%%)\n",
            "Edge density (ED)",
            mean_ED_current, mean_delta_ED, pct_change_ED))
cat(sprintf("  %-22s  current: %8.4f  delta: %+.4f  (%+.2f%%)\n",
            "Patch density (PD)",
            mean_PD_current, mean_delta_PD, pct_change_PD))
cat(sprintf("  %-22s  current: %8.2f  delta: %+.4f  (%+.2f%%)\n",
            "Mean patch area (MPA)",
            mean_MPA_current, mean_delta_MPA, pct_change_MPA))
cat("\n  ED:  negative delta = less edge = less fragmented\n")
cat("  PD:  negative delta = fewer patches = less fragmented\n")
cat("  MPA: positive delta = larger patches = less fragmented\n")

cat("\n============================================================\n")

# ====================================================
# STEP 5: SAVE SUMMARY CSVs
# ====================================================

metric_summary <- data.frame(
  metric        = c("FFI", "ED", "PD", "MPA"),
  mean_current  = round(c(mean_ffi_current,  mean_ED_current,
                          mean_PD_current,   mean_MPA_current), 4),
  mean_scenario = round(c(mean_ffi_scenario,
                          mean_ED_current  + mean_delta_ED,
                          mean_PD_current  + mean_delta_PD,
                          mean_MPA_current + mean_delta_MPA), 4),
  mean_delta    = round(c(mean_delta_ffi, mean_delta_ED,
                          mean_delta_PD,  mean_delta_MPA), 4),
  pct_change    = round(c(pct_change_ffi, pct_change_ED,
                          pct_change_PD,  pct_change_MPA), 2)
)

spatial_summary <- data.frame(
  category     = c("FFI decreased", "unchanged",
                   "FFI increased", "total"),
  n_cells      = c(n_decreased,
                   n_unchanged + n_unmatched,
                   n_increased,
                   n_total_current),
  area_km2     = c(n_decreased,
                   n_unchanged + n_unmatched,
                   n_increased,
                   n_total_current) * cell_area_km2,
  pct_of_total = round(
    c(n_decreased,
      n_unchanged + n_unmatched,
      n_increased,
      n_total_current) / n_total_current * 100, 1)
)

metric_path  <- file.path(output_dir,
                          paste0("metric_summary_",  scenario_name, ".csv"))
spatial_path <- file.path(output_dir,
                          paste0("spatial_summary_", scenario_name, ".csv"))

write.csv(metric_summary,  metric_path,  row.names = FALSE)
write.csv(spatial_summary, spatial_path, row.names = FALSE)

cat("\nSummary CSVs saved:\n")
cat("1.", metric_path,  "\n")
cat("2.", spatial_path, "\n")
cat("\n===== ALL DONE =====\n")

#####

## rasterize current FFI and scenario FFI for comparison
## a quick Madagascar diagnostic:

mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
output_dir    <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results"

# ====================================================
# RASTERIZE CURRENT FOREST FFI
# ====================================================

cat("Rasterizing current forest FFI...\n")

current_points <- vect(
  current_ffi,
  geom = c("center_x", "center_y"),
  crs  = mollweide_crs
)

template_current <- rast(
  ext        = ext(current_points),
  resolution = 5000,
  crs        = mollweide_crs
)

current_rast <- rasterize(current_points, template_current,
                          field = "FFI", fun = "mean")

current_rast_path <- file.path(output_dir, "current_forest_FFI_raster.tif")
writeRaster(current_rast, current_rast_path, overwrite = TRUE)
cat("Saved:", current_rast_path, "\n")

# ====================================================
# RASTERIZE HOLISTIC HOTSPOT FFI
# ====================================================

cat("\nRasterizing holistic hotspot FFI...\n")

scenario_points <- vect(
  scenario_ffi,
  geom = c("center_x", "center_y"),
  crs  = mollweide_crs
)

template_scenario <- rast(
  ext        = ext(scenario_points),
  resolution = 5000,
  crs        = mollweide_crs
)

scenario_rast <- rasterize(scenario_points, template_scenario,
                           field = "FFI", fun = "mean")

scenario_rast_path <- file.path(output_dir,
                                "holistic_hotspot_FFI_raster.tif")
writeRaster(scenario_rast, scenario_rast_path, overwrite = TRUE)
cat("Saved:", scenario_rast_path, "\n")

# ====================================================
# MADAGASCAR DIAGNOSTIC
# Madagascar approx bbox in Mollweide:
# x: 4,400,000 to 5,200,000
# y: -2,800,000 to -1,200,000
# ====================================================

cat("\n===== Madagascar diagnostic =====\n")

madag_delta <- delta_ffi %>%
  filter(center_x > 4400000 & center_x < 5200000 &
           center_y > -2800000 & center_y < -1200000)

cat("Plots in Madagascar region:", format(nrow(madag_delta), big.mark = ","), "\n")

if (nrow(madag_delta) > 0) {
  
  cat("\nFFI change summary:\n")
  cat(sprintf("  Mean delta FFI:  %+.4f\n",
              mean(madag_delta$delta_FFI, na.rm = TRUE)))
  cat(sprintf("  Mean FFI current:   %.4f\n",
              mean(madag_delta$FFI_current, na.rm = TRUE)))
  cat(sprintf("  Mean FFI scenario:  %.4f\n",
              mean(madag_delta$FFI_scenario, na.rm = TRUE)))
  
  cat("\nDirection breakdown:\n")
  cat(sprintf("  FFI decreased: %s cells  (%.1f%%)\n",
              format(sum(madag_delta$delta_class == -1), big.mark = ","),
              mean(madag_delta$delta_class == -1) * 100))
  cat(sprintf("  Unchanged:     %s cells  (%.1f%%)\n",
              format(sum(madag_delta$delta_class ==  0), big.mark = ","),
              mean(madag_delta$delta_class ==  0) * 100))
  cat(sprintf("  FFI increased: %s cells  (%.1f%%)\n",
              format(sum(madag_delta$delta_class ==  1), big.mark = ","),
              mean(madag_delta$delta_class ==  1) * 100))
  
  cat("\nIndividual metric changes:\n")
  cat(sprintf("  Mean delta ED:   %+.4f  (%.1f%%)\n",
              mean(madag_delta$delta_ED,  na.rm = TRUE),
              mean(madag_delta$delta_ED,  na.rm = TRUE) /
                mean(madag_delta$ED_current, na.rm = TRUE) * 100))
  cat(sprintf("  Mean delta PD:   %+.4f  (%.1f%%)\n",
              mean(madag_delta$delta_PD,  na.rm = TRUE),
              mean(madag_delta$delta_PD,  na.rm = TRUE) /
                mean(madag_delta$PD_current, na.rm = TRUE) * 100))
  cat(sprintf("  Mean delta MPA:  %+.4f  (%.1f%%)\n",
              mean(madag_delta$delta_MPA, na.rm = TRUE),
              mean(madag_delta$delta_MPA, na.rm = TRUE) /
                mean(madag_delta$MPA_current, na.rm = TRUE) * 100))
  
  cat("\nIf PD increased strongly and MPA decreased, it confirms\n")
  cat("scattered new patches are driving the FFI increase there.\n")
}

cat("\n===== DONE =====\n")
cat("Three rasters ready to load in ArcGIS:\n")
cat("1.", current_rast_path, "\n")
cat("2.", scenario_rast_path, "\n")
cat("3.", file.path(output_dir,
                    paste0("delta_FFI_holistic_hotspot_continuous.tif")), "\n")