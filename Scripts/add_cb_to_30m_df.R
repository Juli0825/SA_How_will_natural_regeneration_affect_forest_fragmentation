# Step 2 spatailly join with cost and benefit data to 30m dataframe
# This script is to filter out 30m cells in 1km cells
# process each 30m tile separately and spatially join with cost-benefit data

# Load the saved RData file
load("Data/df_natRege.RData")

# Check dataframe
print(head(pnr_df_updated))

load("Data/processed_pnr_tiles/mo_pnv_bin_30m_tile_1_pnr_points.RData")
print(head(pnr_points))

# Load required libraries
library(dplyr)
library(data.table)
library(RANN)

s_time <-Sys.time()
s_time

# Get list of RData files
rdata_files <- list.files("Data/processed_pnr_tiles", pattern = "*.RData$", full.names = TRUE)
print(paste("Found", length(rdata_files), "tile files to process"))

# Convert pnr_df_updated to data.table for efficient processing
dt_costs <- as.data.table(pnr_df_updated)
print(paste("Total 1km pixels:", nrow(dt_costs)))

# Create output directory
output_dir <- "pnr_30m_with_cb"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Set distance threshold (1000m worked perfectly in our test)
distance_threshold <- 1000

# Initialize counters
total_processed <- 0
total_kept <- 0
tiles_with_data <- 0

# Process each tile
for (i in seq_along(rdata_files)) {
  start_time <- Sys.time()
  file_path <- rdata_files[i]
  file_name <- basename(file_path)
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("Processing tile", i, "of", length(rdata_files), ":", file_name, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  # Load the tile data
  load(file_path)
  
  # Convert to data.table
  tile_data <- as.data.table(pnr_points)
  
  print(paste("Loaded tile with", nrow(tile_data), "30m cells"))
  
  # Get tile coordinate ranges
  tile_x_range <- range(tile_data$x)
  tile_y_range <- range(tile_data$y)
  
  cat("Tile coordinate ranges:\n")
  cat("  X:", round(tile_x_range[1]), "to", round(tile_x_range[2]), "\n")
  cat("  Y:", round(tile_y_range[1]), "to", round(tile_y_range[2]), "\n")
  
  # Pre-filter CB data to relevant area (with buffer)
  buffer <- 10000  # 10km buffer
  relevant_cb <- dt_costs[
    Longitude >= (tile_x_range[1] - buffer) & 
      Longitude <= (tile_x_range[2] + buffer) &
      Latitude >= (tile_y_range[1] - buffer) & 
      Latitude <= (tile_y_range[2] + buffer)
  ]
  
  print(paste("Relevant CB points in this area:", nrow(relevant_cb)))
  
  if (nrow(relevant_cb) == 0) {
    cat("No CB points found in this area. Skipping tile...\n")
    rm(pnr_points, tile_data)
    next
  }
  
  # Perform nearest neighbor search
  coords_30m <- as.matrix(tile_data[, .(x, y)])
  coords_cb <- as.matrix(relevant_cb[, .(Longitude, Latitude)])
  
  cat("Performing nearest neighbor search...\n")
  nn_result <- nn2(coords_cb, coords_30m, k = 1)
  nearest_indices <- nn_result$nn.idx[,1]
  nearest_distances <- nn_result$nn.dists[,1]
  
  # Apply distance threshold
  valid_cells <- nearest_distances <= distance_threshold
  n_valid <- sum(valid_cells)
  
  cat("Distance statistics:\n")
  cat("  Min:", round(min(nearest_distances), 1), "m\n")
  cat("  Max:", round(max(nearest_distances), 1), "m\n")
  cat("  Median:", round(median(nearest_distances), 1), "m\n")
  cat("  Cells within", distance_threshold, "m:", n_valid, "\n")
  
  if (n_valid == 0) {
    cat("No cells within distance threshold. Skipping tile...\n")
    rm(pnr_points, tile_data)
    next
  }
  
  # Create final dataframe with cost-benefit values
  cat("Assigning cost-benefit values...\n")
  tile_with_costs <- tile_data[valid_cells]
  
  # Add cost-benefit values from nearest CB pixels
  valid_cb_indices <- nearest_indices[valid_cells]
  tile_with_costs[, carbon := relevant_cb$carbon[valid_cb_indices]]
  tile_with_costs[, bio := relevant_cb$bio[valid_cb_indices]]
  tile_with_costs[, establishment := relevant_cb$establishment[valid_cb_indices]]
  tile_with_costs[, landcosts := relevant_cb$landcosts[valid_cb_indices]]
  
  # Keep only required columns: x, y, carbon, bio, establishment, landcosts
  tile_with_costs <- tile_with_costs[, .(x, y, carbon, bio, establishment, landcosts)]
  
  # Save the processed tile
  output_file <- file.path(output_dir, gsub("\\.RData$", "_with_cb.RData", file_name))
  save(tile_with_costs, file = output_file)
  
  # Update counters
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  total_processed <- total_processed + nrow(tile_data)
  total_kept <- total_kept + nrow(tile_with_costs)
  tiles_with_data <- tiles_with_data + 1
  
  # Calculate retention rate for this tile
  retention_rate <- round(100 * nrow(tile_with_costs) / nrow(tile_data), 1)
  
  cat("Results:\n")
  cat("- Original 30m cells:", nrow(tile_data), "\n")
  cat("- Cells with CB data:", nrow(tile_with_costs), "\n")
  cat("- Retention rate:", retention_rate, "%\n")
  cat("- Processing time:", round(processing_time, 1), "seconds\n")
  cat("- Output saved to:", output_file, "\n")
  
  # Show sample of results
  cat("\nSample of processed data:\n")
  print(head(tile_with_costs, 3))
  
  # Clean up memory
  rm(pnr_points, tile_data, tile_with_costs, coords_30m, coords_cb, nn_result)
  gc()
}

# Final summary
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PROCESSING COMPLETE!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Total tiles processed:", length(rdata_files), "\n")
cat("Tiles with CB data:", tiles_with_data, "\n")
cat("Tiles skipped (no CB coverage):", length(rdata_files) - tiles_with_data, "\n")
cat("Total 30m cells processed:", format(total_processed, big.mark = ","), "\n")
cat("Total cells with CB data:", format(total_kept, big.mark = ","), "\n")
if (total_processed > 0) {
  overall_retention <- round(100 * total_kept / total_processed, 1)
  cat("Overall retention rate:", overall_retention, "%\n")
}

cat("\nOutput files saved in:", output_dir, "\n")
cat("Ready for scenario analysis!\n")

# Show what to do next
cat("\n", paste(rep("-", 40), collapse = ""), "\n")
cat("NEXT STEPS:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("1. Your processed tiles are in:", output_dir, "\n")
cat("2. Each tile has: x, y, carbon, bio, establishment, landcosts\n")
cat("3. Ready for 6 scenario filtering:\n")
cat("   - Cost-based: sort by (establishment + landcosts)\n")
cat("   - Carbon-based: sort by carbon (descending)\n")
cat("   - Bio-based: sort by bio (descending)\n")
cat("   - Holistic: combine metrics as needed\n")
cat("   - Fragmentation-minimizing: spatial optimization\n")

e_time <- Sys.time()
total_time <- e_time - s_time
total_time
 


#### debugging
# Let's understand the geographic coverage of both datasets

cat("=== GEOGRAPHIC COVERAGE ANALYSIS ===\n")

# 1. Look at the FIRST few points in the cost-benefit dataframe
cat("\n1. FIRST few rows of cost-benefit data:\n")
print(head(pnr_df_updated[, c("Longitude", "Latitude", "carbon", "bio")], 10))

# 2. Look at RANDOM sample from cost-benefit data
cat("\n2. RANDOM sample from cost-benefit data:\n")
set.seed(123)
random_indices <- sample(nrow(pnr_df_updated), 10)
print(pnr_df_updated[random_indices, c("Longitude", "Latitude", "carbon", "bio")])

# 3. Overall ranges of cost-benefit data
cat("\n3. OVERALL coordinate ranges of cost-benefit data:\n")
cat("X (Longitude) range:", range(pnr_df_updated$Longitude), "\n")
cat("Y (Latitude) range:", range(pnr_df_updated$Latitude), "\n")

# 4. Check a few different tiles to see their coordinate ranges
cat("\n4. COORDINATE RANGES of different tiles:\n")

# Get first few tile files
tile_files <- list.files("Data/processed_pnr_tiles", pattern = "*.RData$", full.names = TRUE)[1:5]

for (i in 1:min(5, length(tile_files))) {
  file_path <- tile_files[i]
  file_name <- basename(file_path)
  
  # Load tile
  load(file_path)
  tile_data <- pnr_points  # Use base R instead of data.table
  
  cat("\nTile", i, "(", file_name, "):\n")
  cat("  X range:", range(tile_data$x), "\n")
  cat("  Y range:", range(tile_data$y), "\n")
  cat("  Sample coordinates:", tile_data$x[1], ",", tile_data$y[1], "\n")
  
  # Clean up
  rm(pnr_points, tile_data)
}

# 5. Try to find overlap
cat("\n5. LOOKING FOR OVERLAPS:\n")
cat("Checking if any tile coordinates fall within CB data ranges...\n")

cb_x_range <- range(pnr_df_updated$Longitude)
cb_y_range <- range(pnr_df_updated$Latitude)

cat("CB data covers:\n")
cat("  X:", cb_x_range[1], "to", cb_x_range[2], "\n")
cat("  Y:", cb_y_range[1], "to", cb_y_range[2], "\n")

# Check each tile for overlap
for (i in 1:min(10, length(tile_files))) {
  file_path <- tile_files[i]
  file_name <- basename(file_path)
  
  load(file_path)
  tile_data <- pnr_points
  
  tile_x_range <- range(tile_data$x)
  tile_y_range <- range(tile_data$y)
  
  # Check for overlap
  x_overlap <- !(tile_x_range[2] < cb_x_range[1] || tile_x_range[1] > cb_x_range[2])
  y_overlap <- !(tile_y_range[2] < cb_y_range[1] || tile_y_range[1] > cb_y_range[2])
  
  overlap <- x_overlap && y_overlap
  
  cat("Tile", i, "overlap:", overlap, "\n")
  
  if (overlap) {
    cat("  *** POTENTIAL MATCH FOUND! ***\n")
    cat("  Tile X:", tile_x_range, "\n")
    cat("  Tile Y:", tile_y_range, "\n")
    
    # Test a few points from this tile
    buffer <- 5000
    relevant_cb <- pnr_df_updated[
      pnr_df_updated$Longitude >= (tile_x_range[1] - buffer) & 
        pnr_df_updated$Longitude <= (tile_x_range[2] + buffer) &
        pnr_df_updated$Latitude >= (tile_y_range[1] - buffer) & 
        pnr_df_updated$Latitude <= (tile_y_range[2] + buffer), ]
    
    cat("  CB points in this area:", nrow(relevant_cb), "\n")
    
    if (nrow(relevant_cb) > 0) {
      # Calculate actual distances to verify
      library(RANN)
      coords_30m <- as.matrix(tile_data[1:100, c("x", "y")])  # Just test first 100 points
      coords_cb <- as.matrix(relevant_cb[, c("Longitude", "Latitude")])
      
      nn_result <- nn2(coords_cb, coords_30m, k = 1)
      min_dist <- min(nn_result$nn.dists)
      
      cat("  Minimum distance between tile and CB points:", round(min_dist), "meters\n")
      
      if (min_dist <= 1000) {
        cat("  *** SUCCESS! This tile has matching CB data! ***\n")
        break
      }
    }
  }
  
  rm(pnr_points, tile_data)
}

cat("\n=== CONCLUSION ===\n")
cat("This analysis will help us understand:\n")
cat("1. Which geographic regions have CB data\n")
cat("2. Which tiles (if any) overlap with CB data\n")
cat("3. Whether we need to process all tiles or focus on specific ones\n")


#####
# Test on just tile 3
library(dplyr)
library(data.table)
library(RANN)

# Load tile 3
load("Data/processed_pnr_tiles/mo_pnv_bin_30m_tile_11_pnr_points.RData")
tile_data <- as.data.table(pnr_points)

# Convert CB data
dt_costs <- as.data.table(pnr_df_updated)

# Test the spatial join on this tile
# [rest of the spatial joining code from the main script]

# Test spatial join on tile 3
print(paste("Tile 3 - 30m cells:", nrow(tile_data)))
print(paste("Tile 3 coordinate range:"))
print(paste("X:", range(tile_data$x)))
print(paste("Y:", range(tile_data$y)))

# Get tile coordinate ranges for filtering CB data
tile_x_range <- range(tile_data$x)
tile_y_range <- range(tile_data$y)

# Pre-filter CB data to relevant area
buffer <- 10000  # 10km buffer
relevant_cb <- dt_costs[
  Longitude >= (tile_x_range[1] - buffer) & 
    Longitude <= (tile_x_range[2] + buffer) &
    Latitude >= (tile_y_range[1] - buffer) & 
    Latitude <= (tile_y_range[2] + buffer)
]

print(paste("Relevant CB points:", nrow(relevant_cb)))

# Perform nearest neighbor search
coords_30m <- as.matrix(tile_data[, .(x, y)])
coords_cb <- as.matrix(relevant_cb[, .(Longitude, Latitude)])

nn_result <- nn2(coords_cb, coords_30m, k = 1)
nearest_indices <- nn_result$nn.idx[,1]
nearest_distances <- nn_result$nn.dists[,1]

# Check distance distribution
print(paste("Distance statistics:"))
print(summary(nearest_distances))

# Apply distance threshold
distance_threshold <- 1000  # 1km
valid_cells <- nearest_distances <= distance_threshold
print(paste("Cells within", distance_threshold, "m:", sum(valid_cells)))

if (sum(valid_cells) > 0) {
  # Create result with cost-benefit values
  result <- tile_data[valid_cells]
  
  # Add cost-benefit values
  valid_cb_indices <- nearest_indices[valid_cells]
  result[, carbon := relevant_cb$carbon[valid_cb_indices]]
  result[, bio := relevant_cb$bio[valid_cb_indices]]
  result[, establishment := relevant_cb$establishment[valid_cb_indices]]
  result[, landcosts := relevant_cb$landcosts[valid_cb_indices]]
  
  # Keep only required columns
  result <- result[, .(x, y, carbon, bio, establishment, landcosts)]
  
  print("SUCCESS! Sample results:")
  print(head(result))
  print(paste("Final result:", nrow(result), "cells with cost-benefit data"))
} else {
  print("No cells within distance threshold")
}




