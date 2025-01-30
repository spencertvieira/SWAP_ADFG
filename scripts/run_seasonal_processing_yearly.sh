#!/bin/bash

# Set variables
NETCDF_DIR="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw/monthly"
SHAPEFILE="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/resources/BiogeographicRegions_AKSWAP_ABRedited_edited_clipped_again_again_again_again_again_again.shp"
OUTPUT_DIR="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/outputs"
LOG_DIR="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/logs"

# Ensure necessary directories exist
if [[ ! -d "$NETCDF_DIR" || ! -d "$OUTPUT_DIR" || ! -d "$LOG_DIR" ]]; then
    echo "Error: One or more directories do not exist."
    exit 1
fi
if [[ ! -f "$SHAPEFILE" ]]; then
    echo "Error: Shapefile does not exist."
    exit 1
fi

# Variables and seasons to process
VARIABLES=("2m_temperature" "sea_surface_temperature" "total_precipitation" "snowfall" "evaporation" "10m_u_component_of_wind" "10m_v_component_of_wind")
SEASONS=("Annual" "DJF" "MAM" "JJA" "SON" "MJJAS" "ONDJFMAM")

# Submit jobs for each variable
for VARIABLE in "${VARIABLES[@]}"; do
    OUTPUT_CSV="${OUTPUT_DIR}/${VARIABLE}_seasonal_summaries.csv"
    # Submit a job for the variable with all seasons
    sbatch --job-name="${VARIABLE}_processing" \
        --output="${LOG_DIR}/${VARIABLE}_%j.out" \
        --error="${LOG_DIR}/${VARIABLE}_%j.err" \
        --export=NETCDF_DIR="$NETCDF_DIR",VARIABLE="$VARIABLE",SHAPEFILE="$SHAPEFILE",SEASONS="${SEASONS[*]}",OUTPUT_CSV="$OUTPUT_CSV" \
        run_seasonal_processing_yearly.slurm
done
