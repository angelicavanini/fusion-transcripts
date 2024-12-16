#!/usr/bin/env python

################################################################################
# Script: isoquant.py
# Description: This script takes fasta files containing Pacbio reads supporting 
# fusion transcripts and constructs transcript models. Results are saved in 
# subdirectories with the same name as the input file. Run "conda activate 
# isoquant" before use.
# Usage: sbatch isoquant.py
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=300GB
#SBATCH --time=1-12
#SBATCH --export=ALL

import os
from multiprocessing import Pool
import subprocess

# Function to run isoquant.py on a single fasta file
def run_isoquant(fasta_file):
    base_name = os.path.basename(fasta_file).replace('.fasta', '')
    output_dir = f"{base_name}.output"
    command = [
        "isoquant.py",
        "--reference", "/zfs/jacobs/genomes/human/hg38/hg38.fa",
        "--fastq", fasta_file,
        "--data_type", "pacbio_ccs",
        "-o", output_dir
    ]
    subprocess.run(command)

if __name__ == "__main__":
    fasta_dir = "/zfs/jacobs/Colette/fusion_transcripts/pacbio_variable_transcript_structure/fasta_pacbio"
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith('.fasta')]
    
    print("starting isoquant")

    # Use multiprocessing to run isoquant on multiple files in parallel
    with Pool(15) as p:
        p.map(run_isoquant, fasta_files)
    
    print("all done")
