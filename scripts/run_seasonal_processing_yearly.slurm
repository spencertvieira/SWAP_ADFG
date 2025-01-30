#!/bin/bash
#SBATCH --job-name=seasonal_processing
#SBATCH --output=/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/logs/%x_%j.log
#SBATCH --error=/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --partition=t2small
#SBATCH --time=4:00:00  # Adjust time based on expected run duration

# Ensure the log directory exists
mkdir -p /beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/logs

# Load modules (if needed)
module spider GCC/12.3.0 &> /dev/null
if [ $? -eq 0 ]; then
  module load GCC/12.3.0
  export LD_LIBRARY_PATH=/usr/local/pkg/Core/GCCcore/12.3.0/lib64:$LD_LIBRARY_PATH
else
  echo "Warning: GCC/12.3.0 not found, skipping module load."
fi

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate spatial

# Validate input arguments
if [ -z "$NETCDF_DIR" ] || [ -z "$VARIABLE" ] || [ -z "$SHAPEFILE" ] || [ -z "$SEASONS" ] || [ -z "$OUTPUT_CSV" ]; then
  echo "Error: Missing required environment variables."
  exit 1
fi

# Convert SEASONS string to an array
IFS=' ' read -r -a SEASON_ARRAY <<< "$SEASONS"

# Loop through seasons and run processing for each
for SEASON in "${SEASON_ARRAY[@]}"; do
  echo "Processing variable: $VARIABLE for season: $SEASON"
  
  python /beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/scripts/calculate_seasonal_means_yearly.py \
    --netcdf_dir "$NETCDF_DIR" --variable "$VARIABLE" --shapefile "$SHAPEFILE" \
    --season "$SEASON" --output_csv "$OUTPUT_CSV"
done

echo "Processing complete for $VARIABLE."
