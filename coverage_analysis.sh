#!/usr/bin/bash

################################################################################
# Script: coverage_analysis.sh
# Description: This script performs coverage analysis for hipsci WGS BAM files 
#              based on regions specified in a BED file.
# Files required:
#   1. CSV file (named: hipsci-files_wgs_cram.csv) containing information about
#      cell lines (cell_line, accession, file_name), downloaded from hipsci.org.
#   2. (Optional) Text file (named: cell_lines.txt) listing particular cell lines
#      for analysis.
#   3. BED file (named: copy_regions.bed) specifying regions of interest for
#      coverage analysis.
# Usage: sbatch coverage_analysis.sh
# Notes: 
# Manually adjust the samtools view regions in the script (line 134) as needed for your analysis.
# Use "conda activate genotyping" to ACTIVATE CORRECT ENVIRONMENT before running.
# Duplicate the first cell line in the cell_lines.txt file.
# Follow with calculate_ratios.py.
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50GB
#SBATCH --time=1-12
#SBATCH --export=ALL

# Function to convert chromosome names from Ensembl to UCSC format
convert_chromosome_names() {
    LC_ALL=C perl -lpe 's/SN:([0-9]+|[XY]|MT)\b/SN:chr$1/'
}

# Function to fetch the cram file from FTP link and extract specific regions to a single BAM file
fetch_and_extract() {
    cell_line=$1
    accession=$2
    file_name=$3
    ftp_link=$4  # Use ftp_link directly from CSV
    regions=("${@:5}")

    # Generate output BAM file name
    output_bam="${cell_line}_${file_name%.cram}.bam"

    # Check if the output BAM file already exists
    if [ -f "$output_bam" ]; then
        echo "Output BAM file $output_bam already exists. Skipping extraction."
        return
    fi

    # Trim leading and trailing spaces from file_name
    file_name=$(echo "$file_name" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')

    # Trim leading and trailing spaces from ftp_link and remove line breaks
    ftp_link=$(echo "$ftp_link" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e 's/\r//g')

    # Log the output BAM file name
    echo "Output BAM file: $output_bam"

    # Log the command being executed
    echo "Executing: samtools view -bo \"$output_bam\" -h \"$ftp_link\" ${regions[@]}"

    # Fetch and extract specific regions
    if samtools view -bo "$output_bam" -h "$ftp_link" "${regions[@]}"; then
        echo "Extraction successful: $output_bam"
    else
        echo "Error: Failed to extract regions from $ftp_link"
    fi
}

# Function to reheader BAM file to UCSC format
reheader_to_ucsc() {
    input_bam=$1
    output_bam="${input_bam%.bam}.ucsc.bam"

    # Reheader BAM file to UCSC format
    samtools view -H "$input_bam" | convert_chromosome_names > "${input_bam%.bam}.header.ucsc.sam"
    samtools reheader "${input_bam%.bam}.header.ucsc.sam" "$input_bam" > "$output_bam"
}

# Function to sort BAM file
sort_bam() {
    input_bam=$1
    sorted_bam="${input_bam%.ucsc.bam}.sorted.bam"

    # Sort BAM file
    samtools sort -o "$sorted_bam" "$input_bam"
}

# Function to index BAM file
index_bam() {
    input_bam=$1

    # Index BAM file
    samtools index "$input_bam"
}

# Function to perform coverage analysis using bedtools multicov
perform_coverage_analysis() {
    bed_file=$1

    # Perform coverage analysis
    for bam_file in *.sorted.bam; do
        if [ -f "$bam_file" ]; then
            bam_name=$(basename "$bam_file" .sorted.bam)  # Extract the base name without extension
            output_prefix="${bam_name}.vcf"  # Update output file name format

            # Perform coverage analysis for current BAM file
            bedtools multicov -bams "$bam_file" -bed "$bed_file" | sed "s/^/$bam_name\t/" > "$output_prefix"
        fi
    done
}

# Function to filter cell lines based on a list in a text file
filter_cell_lines() {
    cell_lines_file=$1
    csv_file=$2
    filtered_csv="filtered_$csv_file"

    # Remove existing filtered CSV file if it exists
    rm -f "$filtered_csv"

    # Filter CSV file based on cell lines list
    while IFS= read -r line; do
        grep "^$line," "$csv_file" >> "$filtered_csv"
    done < "$cell_lines_file"

    echo "$filtered_csv"
}

# Memory profiling function
memory_profile() {
    echo "Memory usage (RSS): $(ps -p $$ --no-headers -o rss) KB"
}

# Parse CSV file to get necessary information
csv_file="hipsci-files_wgs_cram.csv"

# Check if a cell lines text file is provided
if [ -f "cell_lines.txt" ]; then
    # Filter cell lines based on the provided text file
    cell_lines_file="cell_lines.txt"
    filtered_csv=$(filter_cell_lines "$cell_lines_file" "$csv_file")
else
    # Use the original CSV file without filtering
    filtered_csv="$csv_file"
fi

# Start memory profiling
echo "Memory profiling started..."
memory_profile

# Skip header and read cell line, accession, file name, and ftp_link columns
tail -n +2 "$filtered_csv" | while IFS=',' read -r cell_line _ _ _ _ _ _ _ _ file_name _ ftp_link; do
    # Manually specify multiple regions
    regions=("17:44037000-45000000" "5:68000000-71800000")

    # Fetch and extract regions from cram file
    fetch_and_extract "$cell_line" "$accession" "$file_name" "$ftp_link" "${regions[@]}"

    # Reheader BAM file to UCSC format
    reheader_to_ucsc "${cell_line}_${file_name%.cram}.bam"
    
    # Sort BAM file
    sorted_bam="${cell_line}_${file_name%.cram}.ucsc.sorted.bam"  # Define sorted BAM file name
    sort_bam "${cell_line}_${file_name%.cram}.ucsc.bam" "$sorted_bam"
    
    # Index sorted BAM file
    index_bam "$sorted_bam"
done

# Perform coverage analysis
bed_file="copy_regions.bed"  # Update this to your bed file name if needed

# Memory profiling after processing
echo "Memory profiling after processing..."
memory_profile

perform_coverage_analysis "$bed_file"

# List generated VCF files and save to coverage_vcf_files.txt
ls *.vcf > coverage_vcf_files.txt

echo All done
