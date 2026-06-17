import arcpy
import os
import time
import pandas as pd
from arcpy.sa import *

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

# ===== PATHS =====
lookup_csv  = r"R:\Chapter_3_fragmentation\2026_NEE_R2\gfw_pnr_lookup.csv"
gfw_folder  = r"R:\Chapter_3_fragmentation\frag_2026_exct_median\current_forest_10m"
output_base = r"R:\Chapter_3_fragmentation\2026_NEE_R2\binary_forest_scenarios"
log_file    = os.path.join(output_base, f"processing_log_{time.strftime('%Y%m%d_%H%M%S')}.txt")

# Priority order: holistic_hotspot and all_pnr first, then the rest
scenarios = ["holistic_hotspot", "all_pnr", "high_bio", "low_cost", "high_carbon"]

# Create output subfolders
os.makedirs(output_base, exist_ok=True)
for s in scenarios:
    os.makedirs(os.path.join(output_base, s), exist_ok=True)

# ===== LOGGING =====
def log(msg):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    full_msg  = f"[{timestamp}] {msg}"
    print(full_msg)
    with open(log_file, "a") as f:
        f.write(full_msg + "\n")

log("=" * 60)
log("BINARY FOREST SCENARIO RASTER CREATION")
log("Fix 1: Con(..., 1, 0) instead of Con(..., 1, None)")
log("Fix 2: in_memory paths have no .tif extension")
log("Scenarios: holistic_hotspot + all_pnr first, then high_bio, low_cost, high_carbon")
log("=" * 60)

# ===== READ LOOKUP =====
lookup_df = pd.read_csv(lookup_csv)
log(f"Lookup loaded: {len(lookup_df)} rows")
log(f"GFW tiles in lookup: {lookup_df['gfw_name'].nunique()}")
log(f"Scenarios in lookup: {list(lookup_df['scenario'].unique())}")

# Build dict of all GFW tiles on disk
all_gfw = {
    os.path.basename(f).replace("_10m.tif", ""): os.path.join(gfw_folder, f)
    for f in os.listdir(gfw_folder) if f.endswith("_10m.tif")
}
log(f"Total GFW tiles on disk: {len(all_gfw)}")

# ===== MAIN PROCESSING LOOP =====
total_start   = time.time()
grand_success = 0
grand_failed  = 0
grand_copied  = 0

for scenario in scenarios:

    scenario_start = time.time()
    log(f"\n{'=' * 60}")
    log(f"SCENARIO: {scenario}")
    log(f"{'=' * 60}")

    scenario_df   = lookup_df[lookup_df["scenario"] == scenario]
    output_folder = os.path.join(output_base, scenario)
    pnr_tiles     = set(scenario_df["gfw_name"].unique())

    log(f"GFW tiles with PNR overlap: {len(pnr_tiles)}")
    log(f"GFW tiles with no PNR overlap (copy GFW direct): {len(all_gfw) - len(pnr_tiles)}")

    success = 0
    failed  = 0
    copied  = 0

    for i, (gfw_name, gfw_path) in enumerate(sorted(all_gfw.items()), 1):

        output_file = os.path.join(output_folder, f"{gfw_name}_binary.tif")
        tile_start  = time.time()

        try:
            # ----------------------------------------------------------
            # CASE A: No PNR overlap — future forest == current forest
            # Copy GFW raster directly (already a true binary 0/1)
            # ----------------------------------------------------------
            if gfw_name not in pnr_tiles:
                log(f"  [{i}/{len(all_gfw)}] {gfw_name}: no PNR overlap — copying GFW raster")
                arcpy.management.CopyRaster(
                    gfw_path,
                    output_file,
                    format="TIFF",
                    pixel_type="8_BIT_UNSIGNED"
                )
                tile_time = round((time.time() - tile_start) / 60, 2)
                log(f"  [{i}/{len(all_gfw)}] {gfw_name}: copied in {tile_time} mins")
                copied += 1
                grand_copied += 1
                continue

            # ----------------------------------------------------------
            # CASE B: PNR overlap — overlay GFW + PNR scenario tiles
            # ----------------------------------------------------------
            pnr_files = scenario_df[scenario_df["gfw_name"] == gfw_name]["pnr_file"].tolist()

            valid_pnr = [p for p in pnr_files if os.path.exists(p)]
            missing   = [p for p in pnr_files if not os.path.exists(p)]
            for m in missing:
                log(f"  WARNING: PNR file not found: {m}")

            if not valid_pnr:
                log(f"  [{i}/{len(all_gfw)}] {gfw_name}: ERROR — no valid PNR files, skipping")
                failed += 1
                grand_failed += 1
                continue

            log(f"  [{i}/{len(all_gfw)}] {gfw_name}: {len(valid_pnr)} PNR tile(s)")

            # Snap environment to GFW raster
            gfw_rast = Raster(gfw_path)
            arcpy.env.snapRaster             = gfw_path
            arcpy.env.outputCoordinateSystem = arcpy.Describe(gfw_path).spatialReference
            arcpy.env.extent                 = gfw_rast.extent
            arcpy.env.cellSize               = gfw_rast.meanCellWidth

            # Mosaic PNR tiles if more than one overlaps this GFW tile
            # NOTE: in_memory paths must NOT have a file extension in ArcGIS Pro
            if len(valid_pnr) > 1:
                log(f"  Mosaicking {len(valid_pnr)} PNR tiles...")
                arcpy.management.MosaicToNewRaster(
                    valid_pnr,
                    "in_memory",
                    "pnr_mosaic",          # no .tif extension
                    arcpy.Describe(gfw_path).spatialReference,
                    "32_BIT_FLOAT",
                    None,
                    1,
                    "MAXIMUM",
                    "FIRST"
                )
                pnr_mosaic_path = "in_memory/pnr_mosaic"
            else:
                pnr_mosaic_path = valid_pnr[0]

            # Resample PNR from 30m to 10m, snapped to GFW grid
            # NOTE: output also in in_memory with no extension
            log(f"  Resampling PNR 30m -> 10m (nearest neighbour)...")
            arcpy.management.Resample(
                pnr_mosaic_path,
                "in_memory/pnr_resampled",     # no .tif extension
                f"{gfw_rast.meanCellWidth} {gfw_rast.meanCellHeight}",
                "NEAREST"
            )
            pnr_10m = Raster("in_memory/pnr_resampled")

            # Convert NoData to 0 in both rasters
            gfw_0 = Con(IsNull(gfw_rast), 0, gfw_rast)
            pnr_0 = Con(IsNull(pnr_10m),  0, pnr_10m)

            # KEY FIX: false value = 0 (not None) so non-forest pixels are
            # written as 0 rather than NoData — gives landscapemetrics a
            # true binary raster with class 0 AND class 1
            log(f"  Overlaying GFW + PNR...")
            combined = Con((gfw_0 == 1) | (pnr_0 == 1), 1, 0)

            log(f"  Saving: {output_file}")
            combined.save(output_file)
            arcpy.management.BuildPyramids(output_file)

            # Clean up in_memory
            try:
                arcpy.management.Delete("in_memory/pnr_mosaic")
                arcpy.management.Delete("in_memory/pnr_resampled")
            except Exception:
                pass

            tile_time = round((time.time() - tile_start) / 60, 2)
            log(f"  [{i}/{len(all_gfw)}] {gfw_name}: SUCCESS in {tile_time} mins")
            success += 1
            grand_success += 1

        except Exception as e:
            log(f"  [{i}/{len(all_gfw)}] {gfw_name}: FAILED — {str(e)}")
            failed += 1
            grand_failed += 1
            try:
                arcpy.management.Delete("in_memory")
            except Exception:
                pass

    scenario_time = round((time.time() - scenario_start) / 60, 1)
    log(f"\nScenario {scenario} complete:")
    log(f"  Overlaid (PNR + GFW): {success}")
    log(f"  Copied (no PNR):      {copied}")
    log(f"  Failed:               {failed}")
    log(f"  Time: {scenario_time} mins")

# ===== FINAL SUMMARY =====
total_time = round((time.time() - total_start) / 60, 1)
log(f"\n{'=' * 60}")
log(f"ALL SCENARIOS COMPLETE")
log(f"Total time: {total_time} mins")
log(f"Output folder: {output_base}")
log(f"Grand total — overlaid: {grand_success} | copied: {grand_copied} | failed: {grand_failed}")
log(f"{'=' * 60}")

log("\nOutput raster counts per scenario:")
for s in scenarios:
    folder    = os.path.join(output_base, s)
    n         = len([f for f in os.listdir(folder) if f.endswith(".tif")])
    expected  = len(all_gfw)
    status    = "OK" if n == expected else f"WARNING — expected {expected}, got {n}"
    log(f"  {s}: {n} rasters [{status}]")

arcpy.CheckInExtension("Spatial")
log("Done.")
