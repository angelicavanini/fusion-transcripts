#Investigating Coding Potential of Fusion Transcripts through Ribosomal Profiling Data
#database3_SRR5047771

################################################################################
# # Script: Ribo-seq_FT.py
# Description: 
# Step 1. This script processes a .fastq file containing Ribo-Seq (RS) reads,
# aligns the reads to the indexed transcripts using Bowtie2, converts SAM to BAM, 
# sorts and indexes the BAM files, and extracts reads for each transcript. 
# The results are saved in separate directories for SAM/BAM files and FASTA 
# files for each fusion transcript and sample combination.
# Header names are assigned based on the gene name the RS reads matched to. 
# Step 2. This script filter out RS reads that are not spanning the fusion 
# junction of the FT they matched. 
# Step 3. All the positive selected RS reads, are screened with BlastN, to find
# whether they match some known transcripts as well or they're uniquely matching
# the novel FT. 
# Step 4. This script finally filter out all the RS reads that have matched some
# other known transcripts in the genome. 
################################################################################

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100GB
#SBATCH --time=7-12
#SBATCH --export=ALL

import os
import glob
import subprocess
import pysam
from Bio import SeqIO

# Pathnames to important files 
transcripts_fasta = "/zfs/jacobs/AngelicaVanini/db3/FT_sequences_86bp.fasta"
ribo_seq_fastq = "/zfs/jacobs/AngelicaVanini/db3/SRR5047771_1_trimmed.fq"
bowtie2_index_base = "FT_sequences_86bp_index"
aligned_sam = "aligned_reads_bowtie2.sam"
aligned_bam = "aligned_reads_bowtie2.bam"
sorted_bam = "aligned_reads_bowtie2.sorted.bam"

# Step 1.0: BowTie2 Index of the FT
print("Indexing transcripts...")
subprocess.run(["bowtie2-build", transcripts_fasta, bowtie2_index_base], check=True)

# Step 1.1: Alignment of the ribo-seq reads to the indexed transcripts
print("Aligning reads...")
subprocess.run([
    "bowtie2", 
    "-x", bowtie2_index_base, 
    "-U", ribo_seq_fastq, 
    "-S", aligned_sam, 
    "-L", "16",  # Seed length = 16
    "-N", "0", # Allowing no mismatches in the seed
    "-a"  # all hits 
], check=True)

# Step 1.2: Conversion of SAM files into BAM files, then sorting and indexing of the BAM
print("Processing SAM/BAM files...")
subprocess.run(["samtools", "view", "-bS", aligned_sam, "-o", aligned_bam], check=True)
subprocess.run(["samtools", "sort", aligned_bam, "-o", sorted_bam], check=True)
subprocess.run(["samtools", "index", sorted_bam], check=True)

# To load the transcripts with Seq.io
print("Loading transcripts...")
transcripts = {record.id: str(record.seq) for record in SeqIO.parse(transcripts_fasta, "fasta")}

# The last sorted BAM file neeeds to be open 
bamfile = pysam.AlignmentFile(sorted_bam, "rb")

# Dictionary to store reads for each transcript
transcript_reads = {transcript: [] for transcript in transcripts}

# Step 1.3: Extract reads for each transcript
print("Extracting reads for each transcript...")
for read in bamfile.fetch():
    if not read.is_unmapped:
        transcript = bamfile.get_reference_name(read.reference_id)
        transcript_reads[transcript].append(read)

# Step 1.4: Write each transcript and its reads to a separate FASTA file
print("Writing FASTA files...")
for transcript, reads in transcript_reads.items():
    with open(f"{transcript}_reads.fasta", "w") as f:
        # Write the transcript
        f.write(f">{transcript}\n{transcripts[transcript]}\n")
        
        # Write the reads
        for read in reads:
            read_seq = read.query_sequence
            read_id = read.query_name
            f.write(f">{read_id}\n{read_seq}\n")

bamfile.close()
print("Process completed successfully.")

# Step 1.5: Change headers
def change_headers(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        first_line = f_in.readline().strip()
        if first_line.startswith(">"):
            new_first_header = ">ft_" + first_line[1:]  # Add "ft_" after ">"
            f_out.write(new_first_header + "\n")

            # Extract the base name from the first header
            base_name = new_first_header[4:]  # Remove ">ft_" to get the base name
            count = 1
            
            for line in f_in:
                if line.startswith(">"):
                    new_header = f">rs_{base_name}_{count}"
                    f_out.write(new_header + "\n")
                    count += 1
                else:
                    f_out.write(line)

input_directory = "/zfs/jacobs/AngelicaVanini/db3"
output_directory = "/zfs/jacobs/AngelicaVanini/db3/changed_headers"

os.makedirs(output_directory, exist_ok=True)
input_files = glob.glob(os.path.join(input_directory, "*_reads.fasta"))

for input_file in input_files:
    output_file = os.path.join(output_directory, os.path.basename(input_file))
    change_headers(input_file, output_file)
    print(f"Processed {input_file} and saved to {output_file}")

# Step 2: Filter reads to keep only those crossing 5bp upstream and downstream of the midpoint
def filter_reads(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        lines = f_in.readlines()
        ref_id = lines[0].strip()
        ref_seq = lines[1].strip()
        
        midpoint = len(ref_seq) // 2
        upstream = midpoint - 5
        downstream = midpoint + 5

        f_out.write(f"{ref_id}\n{ref_seq}\n")
        
        for i in range(2, len(lines), 2):
            read_id = lines[i].strip()
            read_seq = lines[i+1].strip()
            
            read_start = lines[1].find(read_seq)
            read_end = read_start + len(read_seq)

            if read_start <= upstream and read_end >= downstream:
                f_out.write(f"{read_id}\n{read_seq}\n")

input_directory = "/zfs/jacobs/AngelicaVanini/db3/changed_headers"
output_directory = "/zfs/jacobs/AngelicaVanini/db3/changed_headers_only_5bp"

os.makedirs(output_directory, exist_ok=True)
input_files = glob.glob(os.path.join(input_directory, "*_reads.fasta"))

for input_file in input_files:
    base_name = os.path.basename(input_file)
    output_file = os.path.join(output_directory, base_name)
    filter_reads(input_file, output_file)
    print(f"Processed {input_file} and saved to {output_file}")

print("Step 6: All files processed successfully.")

# Step 3: Filter for only rs for blastn
def filter_fasta(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        skip_next = False
        for line in f_in:
            if skip_next:
                skip_next = False
                continue
            if line.startswith(">ft_"):
                skip_next = True
                continue
            f_out.write(line)

input_directory = "/zfs/jacobs/AngelicaVanini/db3/changed_headers_only_5bp"
output_directory = "/zfs/jacobs/AngelicaVanini/db3/noFT_justRS_5bp"

os.makedirs(output_directory, exist_ok=True)
input_files = glob.glob(os.path.join(input_directory, "*_reads.fasta"))

for input_file in input_files:
    output_file = os.path.join(output_directory, os.path.basename(input_file))
    filter_fasta(input_file, output_file)
    print(f"Processed {input_file} and saved to {output_file}")

# Step 8: Perform blastN_FT
def run_blastn(input_file, output_file, db_path):
    blastn_command = [
        "blastn",
        "-query", input_file,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid sstart send qcovs pident stitle"
    ]
    subprocess.run(blastn_command)

input_directory = "/zfs/jacobs/AngelicaVanini/db3/noFT_justRS_5bp"
output_directory = "/zfs/jacobs/AngelicaVanini/db3/blastn_results"
db_path = "/zfs/jacobs/genomes/human/Genecode_v44/gencode_v44_transcripts_filtered"

os.makedirs(output_directory, exist_ok=True)
input_files = glob.glob(os.path.join(input_directory, "*_reads.fasta"))

for input_file in input_files:
    base_name = os.path.basename(input_file)
    output_file = os.path.join(output_directory, f"{os.path.splitext(base_name)[0]}_blastN.txt")
    run_blastn(input_file, output_file, db_path)
    print(f"Processed {input_file} and saved to {output_file}")

# Step 4: Filter blastn results
def get_reads_to_discard(blastn_file):
    reads_to_discard = set()
    read_pcov = {}

    with open(blastn_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            if len(columns) < 6:
                continue
            read_name = columns[0]
            pcov = float(columns[4])

            if read_name not in read_pcov:
                read_pcov[read_name] = []

            read_pcov[read_name].append(pcov)

    for read_name, pcovs in read_pcov.items():
        if any(pcov >= 100 for pcov in pcovs):
            reads_to_discard.add(read_name)

    return reads_to_discard

def filter_reads(original_file, reads_to_discard, output_file):
    with open(original_file, 'r') as f_in, open(output_file, 'w') as f_out:
        write_sequence = True
        for line in f_in:
            if line.startswith(">"):
                read_name = line[1:].strip()  # Remove the '>'
                write_sequence = read_name not in reads_to_discard
                if write_sequence:
                    f_out.write(line)
            else:
                if write_sequence:
                    f_out.write(line)

def process_files(changed_headers_dir, blastn_results_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(changed_headers_dir):
        if filename.endswith("_reads.fasta"):
            base_filename = os.path.splitext(filename)[0]
            blastn_filename = f"{base_filename}_blastN.txt"
            original_file_path = os.path.join(changed_headers_dir, filename)
            blastn_file_path = os.path.join(blastn_results_dir, blastn_filename)
            output_file_path = os.path.join(output_dir, filename)

            if os.path.exists(blastn_file_path):
                reads_to_discard = get_reads_to_discard(blastn_file_path)
                filter_reads(original_file_path, reads_to_discard, output_file_path)
                print(f"Processed {filename}, saved to {output_file_path}")
            else:
                print(f"BlastN results not found for {filename}")

changed_headers_dir = "/zfs/jacobs/AngelicaVanini/db3/changed_headers_only_5bp"
blastn_results_dir = "/zfs/jacobs/AngelicaVanini/db3/blastn_results"
output_dir = "/zfs/jacobs/AngelicaVanini/db3/after_blastN"

process_files(changed_headers_dir, blastn_results_dir, output_dir)
