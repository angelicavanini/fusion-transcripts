#!/usr/bin/env python3

################################################################################
# Script: calculate_ratios.py
# Description: This script combines read counts from multiple VCF files for
#              different cell lines, calculates ratios of read counts for specific
#              regions, and writes the results to a file, grouped by region type.
# Input vcf files are generated from coverage_analysis.sh.
#   Each line should contain cell line name, region name, and coverage/read count
#   separated by tabs.
# Usage: sbatch calculate_ratios.py
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=1:00:00
#SBATCH --export=ALL

import os
import sys
from collections import defaultdict

# Define input and output file names
cell_lines_file = "cell_lines.txt"
output_file = "ratios.txt"

def parse_vcf(file_path):
    data = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            region_name = parts[4]
            read_count = int(parts[5])
            if region_name not in data:
                data[region_name] = read_count
            else:
                data[region_name] += read_count
    return data

# Function to calculate ratios
def calculate_ratios(input_files):
    # Dictionary to store summed read counts for each region for each cell line
    cell_line_sum_read_counts = defaultdict(lambda: defaultdict(int))

    # Iterate over input files
    for input_file in input_files:
        # Extract cell line name from the input file name
        cell_line = input_file.split('_')[0]

        # Print message indicating start of processing for the current cell line
        print(f"Processing cell line: {cell_line}", file=sys.stderr)

        # Read the input file
        with open(input_file, 'r') as f:
            for line in f:
                # Split the line
                parts = line.strip().split('\t')
                region = parts[1]  # Region is in the second column
                read_count = int(parts[2])  # Read count is in the third column

                # Add read count to the summed total for the corresponding region and cell line
                cell_line_sum_read_counts[cell_line][region] += read_count

        # Print message indicating completion of processing for the current cell line
        print(f"Finished processing cell line: {cell_line}", file=sys.stderr)

    # Open the output file for writing
    with open(output_file, 'w') as out_file:
        # Calculate ratios for each cell line
        for cell_line, regions in cell_line_sum_read_counts.items():
            # Write the cell line name to the output file
            out_file.write(f"\n{cell_line}\n")
            # Iterate over test regions
            for region, read_count in regions.items():
                # Skip control region
                if region.endswith('_control_region'):
                    continue
                control_region = region.split('_test_region')[0] + '_control_region'
                # Check if control region exists for the test region
                if control_region in regions:
                    ratio = read_count / regions[control_region]
                    out_file.write(f"{region}: {ratio:.10f}\n")
                else:
                    out_file.write(f"Warning: Control region not found for {region}\n")


def main():
    output_files = []  # List to store output filenames
    
    with open(cell_lines_file, 'r') as f:
        cell_lines = list(set(f.read().splitlines()))  # Read unique cell lines
    
    for cell_line in cell_lines:
        print(f"Processing cell line: {cell_line}")  # Print start message
        vcf_files = [file for file in os.listdir() if file.startswith(cell_line) and file.endswith('.vcf')]
        combined_data = {}
        
        for vcf_file in vcf_files:
            file_data = parse_vcf(vcf_file)
            for region_name, read_count in file_data.items():
                if region_name not in combined_data:
                    combined_data[region_name] = read_count
                else:
                    combined_data[region_name] += read_count
        
        output_file = f"{cell_line}_summed.vcf"  # Change output file name
        output_files.append(output_file)  # Append output filename to the list
        with open(output_file, 'w') as f:
            # Commenting out the line below to remove the unnecessary line
            # f.write("# Combined VCF File\n")
            for region_name, read_count in combined_data.items():
                f.write(f"{cell_line}\t{region_name}\t{read_count}\n")
        
        print(f"Finished processing cell line: {cell_line}")  # Print finish message

    # Call the function to calculate ratios
    calculate_ratios(output_files)

if __name__ == "__main__":
    main()
