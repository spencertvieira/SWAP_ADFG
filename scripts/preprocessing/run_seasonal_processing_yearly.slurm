#!/bin/bash
#SBATCH --job-name=calculate_means
#SBATCH --output=../logs/%x_%j.log
#SBATCH --error=../logs/%x_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=32GB
#SBATCH --partition=t2small

# Activate the environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate SWAP_ADFG

# Variables passed as arguments
NETCDF_DIR=$1
VARIABLE=$2
SEASON=$3
SHAPEFILE=$4
OUTPUT_CSV=$5

# Run the Python script
python ../scripts/calculate_seasonal_means_yearly.py \
    --netcdf_dir $NETCDF_DIR \
    --variable $VARIABLE \
    --season $SEASON \
    --shapefile $SHAPEFILE \
    --output_csv $OUTPUT_CSV
