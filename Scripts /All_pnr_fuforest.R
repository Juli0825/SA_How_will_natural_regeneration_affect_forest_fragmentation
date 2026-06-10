# Take out 30m pnr cost-benefit cells from 1km cost-benefit pnr layer
# Resample to 10m
# Overlay with current forest
# Generate a current forest + 9M km2 FUTURE FOREST LAYER
# Load the saved RData file
load("Data/df_natRege.RData")

# Check dataframe
print(head(pnr_df_updated))

# Convert all 30m pnr to dataframe first
library(terra)
library(dplyr)
library(purrr)

a <- rast("Data/pnr_30m/pnv_pct_30m_tile_0_10_-5_5.tif")
print(a)
plot(a)

b<- rast("Data/pnv_bin_30m/pnv_bin_30m_tile_1.tif")
print(b)
plot(b)

c <- rast("Data/mo_bi_pnr/mo_pnv_bin_30m_tile_1.tif")
print(c)
plot(c)

s_time <-Sys.time()
s_time
# Function to process and save each tile as RData
process_and_save_tile <- function(tile_path, output_dir = "processed_pnr_tiles") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  cat("Processing:", basename(tile_path), "\n")
  
  tile <- rast(tile_path)
  pnr_points <- tile %>%  # Give it a meaningful name
    as.data.frame(xy = TRUE) %>%
    filter(.[[3]] == 1) %>%
    select(x, y)
  
  # Save as RData
  tile_name <- tools::file_path_sans_ext(basename(tile_path))
  output_file <- file.path(output_dir, paste0(tile_name, "_pnr_points.RData"))
  save(pnr_points, file = output_file)
  
  cat("  PNR cells:", nrow(pnr_points), "- saved to:", output_file, "\n")
  return(pnr_points)
}

projected_tiles <- list.files("Data/mo_bi_pnr", pattern = "*.tif$", full.names = TRUE)  # Adjust path

# Process all tiles
for (tile_path in projected_tiles) {
  process_and_save_tile(tile_path)
}

#### This doesn't work
#### Too many values that hit integer overflow, exceeds R's memory limits
#### Gonna keep it in seperate tiles and .RData
# Combine all converted pnr binary mollweide 30m dataframe into one
# # Use base R approach with error checking
# rdata_files <- list.files("processed_pnr_tiles", pattern = "*.RData$", full.names = TRUE)
# 
# cat("Found", length(rdata_files), "RData files\n")
# 
# # Load and check each file individually
# all_pnr_dfs <- list()
# for(i in seq_along(rdata_files)) {
#   cat("Loading file", i, ":", basename(rdata_files[i]), "\n")
#   
#   tryCatch({
#     load(rdata_files[i])  # Loads 'pnr_points'
#     
#     # Check if pnr_points exists and has data
#     if(exists("pnr_points") && nrow(pnr_points) > 0) {
#       all_pnr_dfs[[i]] <- pnr_points
#       cat("  -> Added", nrow(pnr_points), "rows\n")
#     } else {
#       cat("  -> Empty or missing data, skipping\n")
#     }
#     
#     rm(pnr_points)  # Clean up
#   }, error = function(e) {
#     cat("  -> Error loading file:", e$message, "\n")
#   })
# }
# 
# # Remove NULL entries
# all_pnr_dfs <- all_pnr_dfs[!sapply(all_pnr_dfs, is.null)]
# 
# cat("Successfully loaded", length(all_pnr_dfs), "dataframes\n")
# 
# # Combine all dataframes
# pnr_combined_df <- do.call(rbind, all_pnr_dfs)
# 
# # Check result
# cat("Total PNR cells:", nrow(pnr_combined_df), "\n")
# head(pnr_combined_df)

e_time <- Sys.time()
total_time <- e_time - s_time
total_time
