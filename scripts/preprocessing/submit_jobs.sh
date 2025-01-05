#!/bin/bash

# Define variables and year ranges
VARIABLES=("2m_temperature" "total_precipitation" "snowfall" "evaporation" "10m_u_component_of_wind" "10m_v_component_of_wind")
START_YEARS=(1950 1975 2000)
END_YEARS=(1974 1999 2024)

# Submit a job for each variable and year range
for VARIABLE in "${VARIABLES[@]}"; do
  for i in "${!START_YEARS[@]}"; do
    sbatch run_era5_download.slurm $VARIABLE ${START_YEARS[$i]} ${END_YEARS[$i]}
  done
done


