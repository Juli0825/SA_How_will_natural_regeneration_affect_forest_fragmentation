# This script is creating pnr scenarios
# Took forever in R, and R keeps aborted
# So do it in Arcgis Pro Python window
# It's loading GFW tiles and match with pnr scenario tiles
# Then adding new pnr scneario forests to GFW current forest
# Save in raster format

import arcpy
import os
import time
import pandas as pd
from arcpy.sa import *

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True
arcpy.env.compression     = "LZW"

# ===== PATHS =====
lookup_csv  = r"R:\Chapter_3_fragmentation\2026_NEE_R2\gfw_pnr_lookup.csv"
output_base = r"R:\Chapter_3_fragmentation\2026_NEE_R2\future_forests"
log_file    = os.path.join(output_base, "processing_log.txt")

scenarios = ["all_pnr", "low_cost", "high_carbon", "high_bio", "holistic_hotspot"]

for s in scenarios:
    os.makedirs(os.path.join(output_base, s), exist_ok=True)

def log(msg):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    full_msg  = f"[{timestamp}] {msg}"
    print(full_msg)
    with open(log_file, "a") as f:
        f.write(full_msg + "\n")

def reset_env():
    arcpy.env.extent     = None
    arcpy.env.cellSize   = None
    arcpy.env.snapRaster = None

# ===== READ LOOKUP =====
lookup_df = pd.read_csv(lookup_csv)
log(f"Lookup loaded: {len(lookup_df)} rows")

gfw_folder = r"R:\Chapter_3_fragmentation\frag_2026_exct_median\current_forest_10m"
all_gfw    = {os.path.basename(f).replace("_10m.tif", ""):
              os.path.join(gfw_folder, f)
              for f in os.listdir(gfw_folder) if f.endswith("_10m.tif")}
log(f"Total GFW tiles: {len(all_gfw)}")

# ===== PROCESS =====
total_start = time.time()

for scenario in scenarios:

    scenario_start = time.time()
    log(f"\n{'='*60}")
    log(f"Processing scenario: {scenario}")
    log(f"{'='*60}")

    scenario_df = lookup_df[lookup_df["scenario"] == scenario]
    success = 0
    failed  = 0
    skipped = 0

    for gfw_name, gfw_file in all_gfw.items():

        tile_start  = time.time()
        output_file = os.path.join(output_base, scenario,
                                   f"{gfw_name}_{scenario}_10m.tif")

        if arcpy.Exists(output_file):
            log(f"  {gfw_name}: already exists, skipping")
            skipped += 1
            continue

        temp_mosaic   = os.path.join(output_base, f"tmp_mosaic_{gfw_name}.tif")
        temp_resample = os.path.join(output_base, f"tmp_resample_{gfw_name}.tif")

        pnr_files = list(scenario_df[
            scenario_df["gfw_name"] == gfw_name]["pnr_file"])

        try:
            arcpy.env.snapRaster             = gfw_file
            arcpy.env.cellSize               = gfw_file
            arcpy.env.extent                 = gfw_file
            arcpy.env.outputCoordinateSystem = arcpy.Describe(gfw_file).spatialReference

            gfw_rast = Raster(gfw_file)

            if len(pnr_files) == 0:
                log(f"  {gfw_name}: no PNR overlap, copying current forest")
                arcpy.management.CopyRaster(gfw_file, output_file,
                                            pixel_type="8_BIT_UNSIGNED")
                tile_time = round((time.time() - tile_start) / 60, 2)
                log(f"  {gfw_name}: done in {tile_time} mins")
                success += 1
                continue

            if len(pnr_files) == 1:
                pnr_raw = Raster(pnr_files[0])
            else:
                arcpy.management.MosaicToNewRaster(
                    input_rasters = pnr_files,
                    output_location = output_base,
                    raster_dataset_name_with_extension = os.path.basename(temp_mosaic),
                    pixel_type      = "8_BIT_UNSIGNED",
                    mosaic_method   = "MAXIMUM",
                    number_of_bands = 1
                )
                pnr_raw = Raster(temp_mosaic)

            arcpy.management.Resample(
                in_raster       = pnr_raw,
                out_raster      = temp_resample,
                resampling_type = "NEAREST"
            )
            pnr_10m = Raster(temp_resample)

            gfw_0    = Con(IsNull(gfw_rast), 0, gfw_rast)
            pnr_0    = Con(IsNull(pnr_10m),  0, pnr_10m)
            combined = Con((gfw_0 == 1) | (pnr_0 == 1), 1, None)
            combined.save(output_file)

            tile_time = round((time.time() - tile_start) / 60, 2)
            log(f"  {gfw_name}: done in {tile_time} mins")
            success += 1

        except Exception as e:
            log(f"  {gfw_name}: FAILED - {str(e)}")
            failed += 1

        finally:
            for tmp in [temp_mosaic, temp_resample]:
                if os.path.exists(tmp):
                    try:
                        arcpy.management.Delete(tmp)
                    except:
                        pass
            reset_env()

    scenario_time = round((time.time() - scenario_start) / 60, 1)
    log(f"Scenario {scenario}: success={success} failed={failed} "
        f"skipped={skipped} time={scenario_time} mins")

total_time = round((time.time() - total_start) / 60, 1)
log(f"\nALL COMPLETE - Total time: {total_time} mins")
log(f"Output: {output_base}")

arcpy.CheckInExtension("Spatial")
