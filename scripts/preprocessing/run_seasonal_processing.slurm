#!/bin/bash
#SBATCH --job-name=seasonal_processing
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

python ../scripts/calculate_seasonal_means.py \
    --variable $1 \
    --start_year $2 \
    --end_year $3 \
    --shapefile $4 \
    --output_csv $5 \
    --variable_type $6

