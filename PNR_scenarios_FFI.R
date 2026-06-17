### Each PNR future forest scenario FFI
### Using global boundaries created from current_forest_10m
### Fixed with round_digits
### So delta_ffi will be 0 if ffi_current and ffi_scenario is the same 

library(data.table)
library(dplyr)
library(tidyr)

# ====================================================
# CHANGE THESE THREE LINES FOR EACH SCENARIO
# ====================================================

scenario_name        <- "holistic_hotspot"
scenario_metrics_dir <- "2026_NEE_R2/FutureScenario_metrics/holistic_hotspot"
current_ffi_path     <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/current_forest_FFI_10m.csv"

# ====================================================
# FIXED PATHS — do not change these
# ====================================================

boundaries_path <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results/global_boundaries_10m_current_forest.csv"
output_dir      <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results"
round_digits    <- 10   # eliminates floating point noise, preserves real differences

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ====================================================
# STEP 1: LOAD GLOBAL BOUNDARIES FROM LOCAL DRIVE
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
  stop("No metrics files found — check scenario_metrics_dir and scenario_name")
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
# Round to consistent digits before saving —
# prevents floating point noise in delta FFI
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
  
  # FFI formula from Ma et al 2023
  mutate(FFI = (ED_norm + PD_norm + (1 - MPA_norm)) / 3) %>%
  
  # Round to consistent digits — CRITICAL to prevent
  # floating point noise creating spurious delta FFI values
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
    # Round delta values — collapses any residual floating
    # point noise to exactly 0 for genuinely unchanged cells
    delta_FFI = round(FFI_scenario - FFI_current, round_digits),
    delta_ED  = round(ED_scenario  - ED_current,  round_digits),
    delta_PD  = round(PD_scenario  - PD_current,  round_digits),
    delta_MPA = round(MPA_scenario - MPA_current, round_digits)
  )

# Verify floating point fix worked
n_exact_zero <- sum(delta_ffi$delta_FFI == 0, na.rm = TRUE)
n_tiny       <- sum(abs(delta_ffi$delta_FFI) > 0 &
                      abs(delta_ffi$delta_FFI) < 1e-9, na.rm = TRUE)

cat(sprintf("\nFloating point check:\n"))
cat(sprintf("  Exactly zero delta FFI:    %s cells\n",
            format(n_exact_zero, big.mark = ",")))
cat(sprintf("  Tiny non-zero (<1e-9):     %s cells\n",
            format(n_tiny,       big.mark = ",")))
cat(sprintf("  (tiny count should be 0)\n\n"))

cat("Plots in current FFI:  ", format(nrow(current_ffi),  big.mark = ","), "\n")
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