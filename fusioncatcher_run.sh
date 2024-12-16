#!/bin/bash

################################################################################
# Script: fusioncatcher_run.sh
# Description: This script processes a list of directories containing FASTQ files,
# runs the FusionCatcher program to detect fusion genes, and saves the results in
# corresponding output directories. 
# Usage: sbatch fusioncatcher_run.sh
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=220GB
#SBATCH --time=10-12
#SBATCH --export=ALL

echo "SLURM job started on $(date '+%Y-%m-%d %H:%M:%S')"

# Set locale to avoid Perl locale warnings
export LC_ALL=C
export LANG=C

# Activate the conda environment
echo "Activating conda environment 'fusioncatcher'..."
source /zfs/omics/personal/cmoses/miniconda3/etc/profile.d/conda.sh
conda activate fusioncatcher
echo "Conda environment 'fusioncatcher' activated."

# Define the function to run FusionCatcher
run_fusioncatcher() {
    local input_dir=$1
    local dir_name=$(basename "$input_dir")
    local output_dir="${input_dir}${dir_name}_output"
    
    echo "Starting FusionCatcher for directory: $input_dir at $(date '+%Y-%m-%d %H:%M:%S')"
    fusioncatcher -d /zfs/jacobs/Colette/fusion_transcripts/fusion_detection/fusioncatcher_build/human_v102/ -i "$input_dir" -o "$output_dir" -p 8
    echo "Finished FusionCatcher for directory: $input_dir at $(date '+%Y-%m-%d %H:%M:%S')"
}

# Export the function to be used by parallel
export -f run_fusioncatcher

# Generate a list of directories containing .fastq.gz files
dir_list=()
for dir in */; do
    if ls "$dir"*.fastq.gz 1> /dev/null 2>&1; then
        dir_list+=("$dir")
    fi
done

# Run the function in parallel using GNU Parallel
parallel --jobs 4 run_fusioncatcher ::: "${dir_list[@]}"

echo "FusionCatcher processing complete."

# Deactivate the conda environment
echo "Deactivating conda environment 'fusioncatcher'..."
conda deactivate
echo "Conda environment 'fusioncatcher' deactivated."

echo "SLURM job finished on $(date '+%Y-%m-%d %H:%M:%S')"
