#!/bin/bash

# Exit immediately if any command fails
set -e

# Define variables and year ranges
VARIABLES=("2m_temperature" "sea_surface_temperature" "total_precipitation" "snowfall" "evaporation" "10m_u_component_of_wind" "10m_v_component_of_wind")
START_YEARS=(1950 1975 2000)
END_YEARS=(1974 1999 2024)

# Log file for job submissions
LOG_FILE="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/logs/era5_submission.log"
mkdir -p $(dirname "$LOG_FILE")  # Ensure the logs directory exists
echo "Job submissions log - $(date)" > "$LOG_FILE"

# Submit a job for each variable and year range
for VARIABLE in "${VARIABLES[@]}"; do
  for i in "${!START_YEARS[@]}"; do
    START_YEAR=${START_YEARS[$i]}
    END_YEAR=${END_YEARS[$i]}

    JOB_ID=$(sbatch /beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/scripts/run_era5_download.slurm \
      "$VARIABLE" "$START_YEAR" "$END_YEAR" | awk '{print $4}')

    echo "Submitted job for $VARIABLE ($START_YEAR - $END_YEAR) with Job ID: $JOB_ID" | tee -a "$LOG_FILE"
  done
done
