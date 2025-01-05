#!/bin/bash
#SBATCH --job-name=era5_download
#SBATCH --output=../logs/%x_%j.log
#SBATCH --error=../logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=16GB
#SBATCH --partition=t2small

source ~/miniconda3/etc/profile.d/conda.sh
conda activate SWAP_ADFG

# Variables passed as arguments
VARIABLE=$1
START_YEAR=$2
END_YEAR=$3

# Run the Python script
python ../scripts/download_era5.py --variable $VARIABLE --start_year $START_YEAR --end_year $END_YEAR

