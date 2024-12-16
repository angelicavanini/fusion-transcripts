#!/usr/bin/env python3

################################################################################
# Script: H1H2_SNP_genotyping.py
# Description: This script performs SNP genotyping for H1 and H2 cell lines
#              using SAMtools and BCFtools.
# Input bam files are generated from coverage_analysis.sh.
# Usage: sbatch H1H2_SNP_genotyping.py
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=4:00:00
#SBATCH --export=ALL

# Import necessary libraries
import os
import vcf
import subprocess
import csv
import glob

# Parse cell_lines.txt to extract unique cell lines
print("Parsing cell_lines.txt")
cell_lines = set()
with open("cell_lines.txt", "r") as f:
    for line in f:
        cell_lines.add(line.strip())
print("Completed")

# Iterate through sorted.bam files and extract sample information from the header
print("Extracting sample information from BAM headers")
bam_files = [f for f in os.listdir() if f.endswith(".sorted.bam")]
cell_line_samples = {}
for bam_file in bam_files:
    cell_line = "_".join(bam_file.split("_")[:2])  # Include the numbers from the filename
    if cell_line in cell_lines:
        # Use samtools to extract sample information from the header
        command = ["samtools", "view", "-H", bam_file]
        header_info = subprocess.check_output(command, universal_newlines=True)
        sample_info = None
        for line in header_info.split("\n"):
            if line.startswith("@RG"):
                parts = line.split("\t")
                for part in parts:
                    if part.startswith("SM:"):
                        sample_info = part[3:]
                        break
        if sample_info:
            if cell_line not in cell_line_samples:  # Avoid duplicates
                cell_line_samples[cell_line] = sample_info
print("Completed")

# Write results to a CSV file
print("Writing results to CSV file")
output_file_path = "cell_line_sample_mapping.csv"
with open(output_file_path, "w", newline="") as output_file:
    writer = csv.writer(output_file)
    writer.writerow(["Cell Line", "Sample ID"])
    for cell_line, sample_info in cell_line_samples.items():
        writer.writerow([cell_line, sample_info])
print("Completed")

# Generate list-of-BAMs.txt
print("Generating list-of-BAMs.txt")
sorted_bam_files = sorted(bam_files)  # Sort the list alphabetically
with open("list-of-BAMs.txt", "w") as f:
    f.write("\n".join(sorted_bam_files))

# Define the target batch size
target_batch_size = 100

# Initialize variables to track current batch and current cell line
current_batch = []
current_cell_line = None

# Function to write current batch to a text file
def write_batch(batch, batch_number):
    with open(f"batch_{batch_number}.txt", "w") as batch_file:
        batch_file.write("\n".join(batch))

# Loop through the sorted BAM files
batch_number = 1
for bam_file in sorted_bam_files:
    # Extract cell line name from the file name
    cell_line = bam_file.split("_")[0]

    # Check if the current batch is full and the cell line changes
    if len(current_batch) >= target_batch_size and current_cell_line != cell_line:
        # Write the current batch to a file and start a new batch
        write_batch(current_batch, batch_number)
        current_batch = [bam_file]
        batch_number += 1
    else:
        # Append the current file to the current batch
        current_batch.append(bam_file)

    current_cell_line = cell_line

# Write the last batch to a file
if current_batch:
    write_batch(current_batch, batch_number)
    
# Loop through the text files and perform mpileup for each batch
for txt_file in os.listdir():
    if txt_file.startswith("batch_") and txt_file.endswith(".txt"):
        with open(txt_file) as f:
            batch_files = f.read().splitlines()
        batch_file_list = " ".join(batch_files)
        output_prefix = f"H1H2_SNP_{txt_file[:-4]}"  # Remove the ".txt" extension from the file name
        pileup_command = f"bcftools mpileup -f /zfs/jacobs/genomes/human/hg19/hg19.fasta -r chr17:44081064-44081064 -Q 13 -q 0 -b {txt_file} -Ou | bcftools call -mv -Ov -o {output_prefix}.vcf"
        os.system(pileup_command)
        
# Print status: finished pileup for all cell lines
print("All pileup operations completed.")

# Open VCF files, parse them, and add cell line information
for txt_file in os.listdir():
    if txt_file.startswith("batch_") and txt_file.endswith(".txt"):
        with open(txt_file) as f:
            batch_files = f.read().splitlines()
        batch_file_list = " ".join(batch_files)
        # Remove the suffix ".txt" from the txt_file name to get the correct prefix
        batch_output_prefix = f"H1H2_SNP_{txt_file[:-4]}"  # Remove the ".txt" extension from the file name
        pileup_command = f"bcftools mpileup -f /zfs/jacobs/genomes/human/hg19/hg19.fasta -r chr17:44081064-44081064 -Q 13 -q 0 -b {txt_file} -Ou | bcftools call -mv -Ov -o {batch_output_prefix}.vcf"
        os.system(pileup_command)
        
        vcf_file = f"{batch_output_prefix}.vcf"
        output_vcf_file = f"{batch_output_prefix}_variants.vcf"
        with open(vcf_file, "r") as f, open(output_vcf_file, "w") as out_file:
            vcf_reader = vcf.Reader(f)
            
            # Write header with additional "CELL_LINE" column
            out_file.write("\t".join(["CHROM", "POS", "REF", "ALT", "SAMPLE", "CELL_LINE", "GENOTYPE", "DP", "AD", "PL"]) + "\n")
            
            # Iterate over records
            for record in vcf_reader:
                chrom = record.CHROM
                pos = record.POS
                ref = record.REF
                alt = record.ALT[0] if record.ALT else "."

                # Iterate over samples
                for sample in record.samples:
                    sample_name = sample.sample
                    cell_line = ""  # Default value if not found
                    for key, value in cell_line_samples.items():
                        if value == sample_name:
                            cell_line = key
                            break
                    
                    genotype = sample["GT"]
                    dp = sample.data.DP if hasattr(sample.data, "DP") else "."
                    ad = ",".join(map(str, sample.data.AD)) if hasattr(sample.data, "AD") else "."
                    pl = ",".join(map(str, sample.data.PL)) if hasattr(sample.data, "PL") else "."

                    out_file.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{sample_name}\t{cell_line}\t{genotype}\t{dp}\t{ad}\t{pl}\n")

# Get a list of all variant files
variant_files = glob.glob("H1H2_SNP_batch_*_variants.vcf")

# Combine all variant files into one
with open("H1H2_SNP_all_variants.vcf", "w") as outfile:
    # Write header from the first variant file
    with open(variant_files[0], "r") as first_file:
        for line in first_file:
            outfile.write(line)
    
    # Write variant records from all other files
    for variant_file in variant_files[1:]:
        with open(variant_file, "r") as infile:
            # Skip the header line
            next(infile)
            for line in infile:
                outfile.write(line)

print("All variant files combined into H1H2_SNP_all_variants.vcf")
                    
print("All done")