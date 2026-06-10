library(terra)
library(sf)
library(sp)
library(landscapemetrics)
library(raster)


tile_1 <- rast("Data/moll_binary/moll_binary_00N_000E.tif")


start_time <- Sys.time()
cat("Starting analysis at:", as.character(start_time), "\n")

# Function to process a single tile
process_tile <- function(tile_path, output_folder) {
  # Extract tile name from path
  tile_name <- tools::file_path_sans_ext(basename(tile_path))
  cat("Processing tile:", tile_name, "\n")
  
  # Load the tile (already in Mollweide projection and binary)
  tile_binary <- tryCatch({
    rast(tile_path)
  }, error = function(e) {
    cat("Error loading tile", tile_name, ":", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(tile_binary)) return(NULL)
  
  # Get the extent of the raster
  ext_moll <- ext(tile_binary)
  
  # Create the grid
  grid_res <- 5000
  cat("Creating grid with resolution:", grid_res, "m\n")
  grid <- rast(ext=ext_moll, resolution=grid_res, crs=crs(tile_binary))
  
  # Fill with cell numbers
  values(grid) <- 1:ncell(grid)
  
  # Convert to polygons using terra
  cat("Converting grid to polygons...\n")
  grid_poly <- as.polygons(grid)
  cat("Created", nrow(grid_poly), "grid polygons\n")
  
  # Clean up grid to free memory
  rm(grid)
  gc()
  
  # Add ID column to the SpatVector
  grid_poly$FID <- 1:nrow(grid_poly)
  
  # Calculate metrics using terra objects
  cat("Calculating landscape metrics...\n")
  metrics <- tryCatch({
    # Try with terra SpatVector first
    result <- sample_lsm(tile_binary, grid_poly, what = c("lsm_c_ed", "lsm_c_pd", "lsm_c_area_mn"))
    
    # Immediately clean up large objects
    rm(tile_binary, grid_poly)
    gc()
    result
  }, error = function(e) {
    cat("Error calculating metrics with terra objects:", e$message, "\n")
    
    # If there was an error, try converting to sf
    if(exists("tile_binary") && exists("grid_poly")) {
      cat("Trying with sf objects instead...\n")
      
      tryCatch({
        # Convert to sf
        grid_poly_sf <- sf::st_as_sf(grid_poly)
        
        # Calculate metrics with sf
        result <- sample_lsm(tile_binary, grid_poly_sf, what = c("lsm_c_ed", "lsm_c_pd", "lsm_c_area_mn"))
        
        # Cleanup
        rm(tile_binary, grid_poly, grid_poly_sf)
        gc()
        
        return(result)
      }, error = function(e2) {
        cat("Error with sf objects too:", e2$message, "\n")
        
        # Clean up any objects that might still exist
        if(exists("tile_binary")) rm(tile_binary)
        if(exists("grid_poly")) rm(grid_poly)
        if(exists("grid_poly_sf")) rm(grid_poly_sf)
        gc()
        return(NULL)
      })
    } else {
      # Clean up any objects that might still exist
      if(exists("tile_binary")) rm(tile_binary)
      if(exists("grid_poly")) rm(grid_poly)
      gc()
      return(NULL)
    }
  })
  
  if(is.null(metrics)) return(NULL)
  
  # Add tile ID
  metrics$tile_id <- tile_name
  
  # Proper NA filtering
  metrics_clean <- metrics[!is.na(metrics$value), ]
  
  # Save metrics
  write.csv(metrics_clean, file.path(output_folder, paste0(tile_name, "_metrics.csv")), 
            row.names = FALSE)
  
  cat("Processed tile:", tile_name, "- Saved", nrow(metrics_clean), "rows\n")
  
  # Final cleanup
  rm(metrics, metrics_clean)
  gc()
  
  return(tile_name)
}

# Create output folder
output_folder <- "Data/metrics_before_pnr"
dir.create(output_folder, showWarnings = FALSE)

# Get list of all tile files
tiles_folder <- "Data/moll_binary/"
tile_files <- list.files(tiles_folder, pattern = "\\.tif$", full.names = TRUE)

# Get list of already processed files
processed_files <- list.files(output_folder, pattern = "_metrics.csv")
processed_tiles <- gsub("_metrics.csv", "", processed_files)

cat("Found", length(tile_files), "total tile files\n")
cat("Already processed", length(processed_tiles), "tiles\n")

# Process all tiles, skipping the already processed ones
for(tile_path in tile_files) {
  tile_name <- tools::file_path_sans_ext(basename(tile_path))
  
  # Skip if already processed
  if(tile_name %in% processed_tiles) {
    cat("Skipping already processed tile:", tile_name, "\n")
    next
  }
  
  tryCatch({
    process_tile(tile_path, output_folder)
    # Force garbage collection after each tile
    gc()
  }, error = function(e) {
    cat("Error processing tile", tile_name, ":", e$message, "\n")
    # Force cleanup even after errors
    gc()
  })
  
  # Display memory usage after each tile
  cat("Memory in use after tile:", format(utils::object.size(globalenv()), units = "auto"), "\n")
}

cat("Processing complete!\n")

end_time <- Sys.time()
cat("Finished at:", as.character(end_time), "\n")
cat("Total time:", format(end_time - start_time), "\n")


