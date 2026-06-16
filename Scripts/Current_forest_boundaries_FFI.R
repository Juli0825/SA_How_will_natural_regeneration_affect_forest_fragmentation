#### Create global boundaris
#### Get current FFI
#### Lower limit of global bounderies is 0, FFI from 10m - 5km


library(data.table)
library(dplyr)
library(tidyr)

# ====================================================
# PATHS
# ====================================================

current_metrics_dir <- "R:/Chapter_3_fragmentation/frag_2026_exct_median/frag_metrics_current_10m"
output_dir          <- "R:/Chapter_3_fragmentation/2026_NEE_R2/FFI_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ====================================================
# STEP 1: LOAD ALL CURRENT FOREST 10M METRICS
# ====================================================

cat("===== STEP 1: Loading current forest 10m metrics =====\n")

current_files <- list.files(current_metrics_dir,
                            pattern = "_metrics\\.csv$",
                            full.names = TRUE)
cat("CSV files found:", length(current_files), "\n")

all_metrics <- rbindlist(lapply(current_files, fread), fill = TRUE)
cat("Total rows loaded:", format(nrow(all_metrics), big.mark = ","), "\n")
cat("Tiles present:    ", length(unique(all_metrics$tile_id)), "\n")

# Filter to forest class and FFI metrics only
forest_metrics <- all_metrics[class == 1 & metric %in% c("ed", "pd", "area_mn")]
cat("Rows after filter:", format(nrow(forest_metrics), big.mark = ","), "\n")

rm(all_metrics)
gc()

# ====================================================
# STEP 2: CALCULATE AND SAVE GLOBAL BOUNDARIES
# IQR method from Ma et al 2023
# Negative lower bounds floored at 0 — metrics are
# physically impossible to be negative, and flooring
# ensures a cell with ED/PD/area_mn = 0 normalizes
# correctly to 0 rather than a spurious positive value
# ====================================================

cat("\n===== STEP 2: Calculating global IQR boundaries =====\n")

boundaries <- forest_metrics %>%
  as_tibble() %>%
  group_by(metric) %>%
  summarise(
    q1  = quantile(value, 0.25, na.rm = TRUE),
    q3  = quantile(value, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr,
    .groups = "drop"
  ) %>%
  select(metric, lower, upper) %>%
  pivot_longer(cols      = c(lower, upper),
               names_to  = "boundary_type",
               values_to = "value") %>%
  
  # Floor negative lower bounds at zero
  # IQR can produce negative lowers for right-skewed
  # distributions even when no real values are negative
  mutate(value = ifelse(boundary_type == "lower" & value < 0, 0, value))

cat("\nGlobal boundaries (with negative lowers floored at 0):\n")
print(boundaries)

# Sanity checks
area_upper <- boundaries$value[boundaries$metric == "area_mn" &
                                 boundaries$boundary_type == "upper"]
area_lower <- boundaries$value[boundaries$metric == "area_mn" &
                                 boundaries$boundary_type == "lower"]

if (area_upper < 200) {
  warning("area_mn upper looks too low (", round(area_upper, 1),
          ") — check input data is 10m metrics")
} else {
  cat("\nSanity check passed: area_mn upper =", round(area_upper, 1),
      "ha (correct for 10m data)\n")
}

if (any(boundaries$value[boundaries$boundary_type == "lower"] < 0)) {
  warning("Some lower bounds still negative — check data")
} else {
  cat("Sanity check passed: all lower bounds >= 0\n")
}

# Save — this file gets copied to all other servers
boundaries_path <- file.path(output_dir,
                             "global_boundaries_10m_current_forest.csv")
write.csv(boundaries, boundaries_path, row.names = FALSE)
cat("\nBoundaries saved to:\n", boundaries_path, "\n")

# ====================================================
# STEP 3: CALCULATE CURRENT FOREST FFI
# ====================================================

cat("\n===== STEP 3: Calculating current forest FFI =====\n")

# Extract boundary values
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

current_ffi <- forest_metrics %>%
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
  
  # Normalize using boundary values as fixed range
  mutate(
    ED_norm  = (ED_capped  - ed_min)  / (ed_max  - ed_min),
    PD_norm  = (PD_capped  - pd_min)  / (pd_max  - pd_min),
    MPA_norm = (MPA_capped - mpa_min) / (mpa_max - mpa_min)
  ) %>%
  
  # FFI = (ED_norm + PD_norm + (1 - MPA_norm)) / 3
  mutate(FFI = (ED_norm + PD_norm + (1 - MPA_norm)) / 3)

cat("\nFFI range:", round(range(current_ffi$FFI, na.rm = TRUE), 4), "\n")
cat("FFI mean: ", round(mean(current_ffi$FFI,  na.rm = TRUE), 4), "\n")
cat("Total plots:", format(nrow(current_ffi), big.mark = ","), "\n")

# Quick check: fully forested plots should have FFI near 0
fully_forested <- current_ffi %>%
  filter(ED == 0 & PD <= 0.01)
if (nrow(fully_forested) > 0) {
  cat("\nSpot check — near-intact forest plots (ED=0):\n")
  cat("Mean FFI:", round(mean(fully_forested$FFI, na.rm = TRUE), 4),
      "(should be near 0)\n")
}

ffi_path <- file.path(output_dir, "current_forest_FFI_10m.csv")
write.csv(current_ffi, ffi_path, row.names = FALSE)
cat("\nCurrent forest FFI saved to:\n", ffi_path, "\n")

cat("\n===== DONE =====\n")
cat("Files ready to copy to other servers:\n")
cat("1.", boundaries_path, "\n")
cat("2.", ffi_path, "\n")
