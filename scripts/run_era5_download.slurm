#!/bin/bash
#SBATCH --job-name=era5_download
#SBATCH --output=/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/logs/%x_%j.log
#SBATCH --error=/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --partition=t2small
#SBATCH --time=2:00:00  # Set a reasonable max runtime

# Ensure SLURM logs directory exists
mkdir -p /beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/logs

# Load necessary modules (Check if GCC 12.3.0 exists first)
module spider GCC/12.3.0 &> /dev/null
if [ $? -eq 0 ]; then
  module load GCC/12.3.0
  export LD_LIBRARY_PATH=/usr/local/pkg/Core/GCCcore/12.3.0/lib64:$LD_LIBRARY_PATH
else
  echo "Warning: GCC/12.3.0 not found, skipping module load."
fi

# Activate the conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate spatial

# Check if the correct number of arguments is provided
if [ $# -ne 3 ]; then
  echo "Usage: sbatch run_era5_download.slurm <VARIABLE> <START_YEAR> <END_YEAR>"
  exit 1
fi

# Assign variables from arguments
VARIABLE=$1
START_YEAR=$2
END_YEAR=$3
SHAPEFILE_PATH="/beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/data/resources/bounding_box.shp"

# Verify the shapefile exists
if [ ! -f "$SHAPEFILE_PATH" ]; then
  echo "Error: Shapefile $SHAPEFILE_PATH not found!"
  exit 1
fi

# Print job info for debugging
echo "Running ERA5 download for:"
echo "  Variable: $VARIABLE"
echo "  Start Year: $START_YEAR"
echo "  End Year: $END_YEAR"
echo "  Shapefile: $SHAPEFILE_PATH"

# Run the Python script
python /beegfs/CMIP6/stvieira/projects/2025/SWAP_ADFG/scripts/preprocessing/download_era5.py \
    --variable "$VARIABLE" --start_year "$START_YEAR" --end_year "$END_YEAR" \
    --shapefile "$SHAPEFILE_PATH"
