#!/bin/bash

# Variables and parameters
VARIABLES=("2m_temperature" "total_precipitation" "snowfall" "evaporation" "10m_u_component_of_wind" "10m_v_component_of_wind")
PERIODS=("1950-1999" "2005-2024")
VARIABLE_TYPES=("mean" "sum")
SHAPEFILE="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/resources/shapefiles/BiogeographicRegions_AKSWAP_ABRedited/BiogeographicRegions_AKSWAP_ABRedited.shp"
OUTPUT_DIR="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/results"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Submit jobs for each variable and period
for VARIABLE in "${VARIABLES[@]}"; do
  for PERIOD in "${PERIODS[@]}"; do
    START_YEAR=$(echo $PERIOD | cut -d'-' -f1)
    END_YEAR=$(echo $PERIOD | cut -d'-' -f2)
    VARIABLE_TYPE=$( [[ $VARIABLE == "2m_temperature" ]] && echo "mean" || echo "sum" )
    OUTPUT_CSV="${OUTPUT_DIR}/${VARIABLE}_${PERIOD}.csv"
    sbatch run_seasonal_processing.slurm $VARIABLE $START_YEAR $END_YEAR $SHAPEFILE $OUTPUT_CSV $VARIABLE_TYPE
  done
done

