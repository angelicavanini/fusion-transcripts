#!/usr/bin/env python

################################################################################
# Script: bowtie2_alignment_mayo.py
# Description: This script processes a list of fastq.gz files, aligns the reads 
# to the indexed transcripts using Bowtie2, converts SAM to BAM, sorts and indexes 
# the BAM files, and extracts reads for each transcript. The results are saved 
# in separate directories for SAM/BAM files and FASTA files for each fusion 
# transcript and sample combination. This version uses parallel processing for
# each sample.
# Use conda activate bowtie2_py38 before running
# Usage: sbatch bowtie2_alignment_mayo.py
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100GB
#SBATCH --time=7-12
#SBATCH --export=ALL

import os
import subprocess
import pysam
import gzip
from Bio import SeqIO
from multiprocessing import Pool
import math

# Paths to files and tools
transcripts_fasta = "FT_sequences_86bp_extraKANSL1.fasta"
fastq_list_file = "fastq_files.txt"  # Text file with list of fastq.gz files
fastq_dir = "fastq"  # Subdirectory containing fastq.gz files
bowtie2_index_base = "FT_sequences_86bp_extraKANSL1_index"
output_dir = "bowtie2_aligned_reads"
fasta_dir = os.path.join(output_dir, "fasta")
sam_bam_dir = os.path.join(output_dir, "sam_bam")

# Create output directories if they don't exist
os.makedirs(fasta_dir, exist_ok=True)
os.makedirs(sam_bam_dir, exist_ok=True)

print("Starting the script...")

# Step 1: Index the transcripts using Bowtie2
print("Indexing transcripts with Bowtie2...")
subprocess.run(["bowtie2-build", transcripts_fasta, bowtie2_index_base], check=True)
print("Indexing completed.")

# Read the list of fastq.gz files
print("Reading the list of fastq.gz files...")
with open(fastq_list_file, "r") as f:
    fastq_files = [os.path.join(fastq_dir, line.strip()) for line in f]
print(f"Found {len(fastq_files)} fastq.gz files to process.")

# Load the transcripts
print("Loading transcripts from the FASTA file...")
transcripts = {record.id: str(record.seq) for record in SeqIO.parse(transcripts_fasta, "fasta")}
print(f"Loaded {len(transcripts)} transcripts.")

# Create subdirectories for each transcript in the fasta directory
print("Creating subdirectories for each transcript in the FASTA directory...")
for transcript in transcripts:
    os.makedirs(os.path.join(fasta_dir, transcript), exist_ok=True)
print("Subdirectories created.")

def process_sample(rna_seq_fastq):
    # Modify sample_name extraction to exclude file extension
    sample_name = os.path.basename(rna_seq_fastq).replace('.fastq.gz', '')
    print(f"Processing sample: {sample_name}...")

    # Generate paths for SAM, BAM, and sorted BAM files
    aligned_sam = os.path.join(sam_bam_dir, f"{sample_name}_aligned_reads.sam")
    aligned_bam = os.path.join(sam_bam_dir, f"{sample_name}_aligned_reads.bam")
    sorted_bam = os.path.join(sam_bam_dir, f"{sample_name}_aligned_reads.sorted.bam")

    # Step 2: Align the rna-seq reads to the indexed transcripts
    print(f"Aligning reads for {rna_seq_fastq} with Bowtie2...")
    subprocess.run(["bowtie2", "--local", "--no-unal", "-N", "1", "--score-min", "G,32,25", "-x", bowtie2_index_base, "-U", rna_seq_fastq, "-S", aligned_sam], check=True) # Use G,32,25 or G,34,20
    print(f"Alignment completed for {rna_seq_fastq}.")

    # Step 3: Convert SAM to BAM, sort and index BAM
    print(f"Converting SAM to BAM for {sample_name}...")
    subprocess.run(["samtools", "view", "-bS", aligned_sam, "-o", aligned_bam], check=True)
    print(f"SAM to BAM conversion completed for {sample_name}.")

    print(f"Sorting BAM file for {sample_name}...")
    subprocess.run(["samtools", "sort", aligned_bam, "-o", sorted_bam], check=True)
    print(f"BAM sorting completed for {sample_name}.")

    print(f"Indexing sorted BAM file for {sample_name}...")
    subprocess.run(["samtools", "index", sorted_bam], check=True)
    print(f"BAM indexing completed for {sample_name}.")

    # Open the sorted BAM file
    print(f"Opening sorted BAM file for {sample_name}...")
    bamfile = pysam.AlignmentFile(sorted_bam, "rb")

    # Dictionary to store reads for each transcript
    transcript_reads = {transcript: [] for transcript in transcripts}

    # Step 4: Extract reads for each transcript
    print(f"Extracting reads for each transcript for {sample_name}...")
    for read in bamfile.fetch():
        if not read.is_unmapped:
            transcript = bamfile.get_reference_name(read.reference_id)
            transcript_reads[transcript].append(read)
    print(f"Read extraction completed for {sample_name}.")

    # Step 5: Write each transcript and its reads to a separate FASTA file
    print(f"Writing FASTA files for {sample_name}...")
    for transcript, reads in transcript_reads.items():
        transcript_fasta_path = os.path.join(fasta_dir, transcript, f"{transcript}_{sample_name}.fasta")
        with open(transcript_fasta_path, "w") as f:
            # Write the reads only if there are any
            if reads:
                for read in reads:
                    read_seq = read.query_sequence
                    read_id = read.query_name
                    f.write(f">{transcript}_{read_id}\n{read_seq}\n")
    print(f"FASTA files written for {sample_name}.")

    bamfile.close()
    print(f"Processing completed successfully for {rna_seq_fastq}.")

# Parallel processing of samples
print("Starting parallel processing of samples...")
with Pool(32) as pool:
    pool.map(process_sample, fastq_files)

print("All samples processed successfully.")
