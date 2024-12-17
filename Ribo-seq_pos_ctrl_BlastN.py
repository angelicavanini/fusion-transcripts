#BlastN analysis for Positive control filtering
################################################################################
# # Script: Ribo-seq_pos_ctrl_BlastN.py
# Description: 
# The analysis for positive control are identical, exept for the last step of 
# filtering out the BlastN results of the RS reads. 
# This script processes BlastN results of the RS reads in the positive control.
# RS reads matching other known transcripts are discarded. 
################################################################################

# Step 4: Process BlastN results and filter out reads mapping to multiple genes
print("Processing BlastN results...")

import os

def parse_blastn_file(blastn_file):
    """Parse the BlastN results file to identify reads that map to multiple genes with ≥100% coverage."""
    reads_dict = {}

    with open(blastn_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            if len(columns) < 7:
                continue
            read_name = columns[0]
            second_col = columns[1]
            coverage = float(columns[4])  # Coverage percentage
            ensg_id = second_col.split('|')[1][4:16]  # Extract ENSG ID

            if coverage >= 100:
                if read_name not in reads_dict:
                    reads_dict[read_name] = set()
                reads_dict[read_name].add(ensg_id)

    return reads_dict

def identify_reads_to_discard(reads_dict):
    """Identify reads that map to multiple genes."""
    reads_to_discard = {read_name for read_name, ensg_ids in reads_dict.items() if len(ensg_ids) > 1}
    return reads_to_discard

def filter_fasta_file(fasta_file, reads_to_discard, output_file):
    """Filter out reads that map to multiple genes from the FASTA file."""
    with open(fasta_file, 'r') as original_file, open(output_file, 'w') as filtered_file:
        first_line = original_file.readline()
        filtered_file.write(first_line)  # Write the >ft header
        
        write_sequence = True
        for line in original_file:
            if line.startswith(">"):
                read_name = line[1:].strip()  # Remove the '>'
                write_sequence = read_name not in reads_to_discard
                if write_sequence:
                    filtered_file.write(line)
            else:
                if write_sequence:
                    filtered_file.write(line)

def process_files(fasta_dir, blastn_results_dir, output_dir):
    """Process each file in the FASTA directory with corresponding BlastN results."""
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(fasta_dir):
        if filename.endswith("_reads.fasta"):
            base_filename = os.path.splitext(filename)[0]
            blastn_filename = f"{base_filename}_blastN.txt"
            fasta_file_path = os.path.join(fasta_dir, filename)
            blastn_file_path = os.path.join(blastn_results_dir, blastn_filename)
            output_file_path = os.path.join(output_dir, filename)

            if os.path.exists(blastn_file_path):
                reads_dict = parse_blastn_file(blastn_file_path)
                reads_to_discard = identify_reads_to_discard(reads_dict)
                filter_fasta_file(fasta_file_path, reads_to_discard, output_file_path)
                print(f"Processed {filename}, saved to {output_file_path}")
            else:
                print(f"BlastN results not found for {filename}")

# Define the directories
fasta_dir = "/zfs/jacobs/AngelicaVanini/posctrl_db3/set1/changed_headers_only_5bp"
blastn_results_dir = "/zfs/jacobs/AngelicaVanini/posctrl_db3/set1/blastn_results"
output_dir = "/zfs/jacobs/AngelicaVanini/posctrl_db3/set1/after_blastN_100"

# Process the files
process_files(fasta_dir, blastn_results_dir, output_dir)
print("All steps completed successfully.")
