#!/bin/bash

# Set variables
NETCDF_DIR="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/raw"
SHAPEFILE="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/resources/shapefiles/BiogeographicRegions_AKSWAP_ABRedited/BiogeographicRegions_AKSWAP_ABRedited.shp"
OUTPUT_DIR="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/results"

VARIABLES=("2m_temperature" "total_precipitation" "snowfall" "evaporation" "10m_u_component_of_wind" "10m_v_component_of_wind")
SEASONS=("DJF" "MAM" "JJA" "SON" "Annual" "MJJAS" "ONDJFMAM")

# Submit jobs for each variable and season
for VARIABLE in "${VARIABLES[@]}"; do
    for SEASON in "${SEASONS[@]}"; do
        OUTPUT_CSV="${OUTPUT_DIR}/${VARIABLE}_seasonal_means.csv"
        sbatch run_seasonal_processing_yearly.slurm \
            $NETCDF_DIR \
            $VARIABLE \
            $SEASON \
            $SHAPEFILE \
            $OUTPUT_CSV
    done
done

