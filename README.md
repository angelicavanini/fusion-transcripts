## Workflow Overview

1. Run FusionCatcher Analysis (fusioncatcher_run.sh):  
   - Processes a list of directories containing FASTQ files.  
   - Executes the FusionCatcher program to detect fusion genes.  
   - Saves the results in corresponding output directories.  

2. Align Reads to Fusion Transcript Junction Sequences:  
   - Takes fusion transcripts identified by FusionCatcher and maps supporting reads to the junction sequences using Bowtie2 as a custom reference.  
   - This step uses three different scripts to accommodate the naming conventions of specific datasets:  
     - bowtie2_alignment_mayo.py: Designed for the Mayo RNA sequencing study.  
     - bowtie2_alignment_hipsci.py: Designed for the HipSci RNA sequencing dataset.  
     - bowtie2_alignment_hipsci_controls.py: For the HipSci control dataset.  
   - Steps performed by these scripts:
     - Align reads to indexed fusion transcripts using Bowtie2.  
     - Convert SAM to BAM files, sort, and index them.  
     - Extract reads for each fusion transcript.  
   - Outputs:
     - SAM/BAM files stored in structured directories.  
     - FASTA files for each fusion transcript and sample combination.  

3. Count Reads Mapped to Fusion Transcripts (bowtie2_counts.py):  
   - Counts the number of reads mapped to each fusion transcript for each sample using the FASTA files generated from Bowtie2 alignments.  
   - Outputs results as Counts_pivot_[date+time].txt.

4. Perform Coverage Analysis for WGS BAM Files (coverage_analysis.sh):  
   - Analyzes coverage in specified regions of HipSci WGS BAM files.  
   - Requires the following files:  
     - CSV file (hipsci-files_wgs_cram.csv) with metadata for cell lines.  
     - Optional text file (cell_lines.txt) listing specific cell lines to analyze.  
     - BED file (copy_regions.bed) specifying regions of interest.  

5. Calculate Read Count Ratios (calculate_ratios.py):  
   - Combines read counts from multiple VCF files generated by coverage_analysis.sh.  
   - Calculates ratios of read counts for specific regions and groups the results by region type.  

6. Perform SNP Genotyping (H1H2_SNP_genotyping.py):  
   - Performs SNP genotyping for H1 and H2 cell lines using SAMtools and BCFtools.  
   - Requires BAM files generated from coverage_analysis.sh.  

7. Construct Transcript Models (isoquant.py):  
   - Uses FASTA files containing PacBio reads supporting fusion transcripts to construct transcript models.  
   - Saves results in subdirectories corresponding to the input file names.

8. FusionCatcher Processing (FusionCatcher_processing.R):  
   - Takes outputs from all previous steps and combines them for visualization.  

9. Analyze Inter-Individual Variation in Fusion Transcript Expression (variable_fusion_expression.R):  
   - Focuses on fusion transcripts with inter-individual variation.  
   - Combines outputs from previous steps for targeted visualization.

---

## Repository Contents

### Fusion Gene Detection
- fusioncatcher_run.sh:  
  Runs FusionCatcher on a list of directories containing FASTQ files.  
  Outputs FusionCatcher results in structured directories.

### Bowtie2 Alignment and Read Processing
- bowtie2_alignment_mayo.py:  
  Aligns reads to fusion transcript junction sequences for the Mayo RNA sequencing dataset.

- bowtie2_alignment_hipsci.py:  
  Aligns reads to fusion transcript junction sequences for the HipSci RNA sequencing dataset.

- bowtie2_alignment_hipsci_controls.py:  
  Aligns reads to fusion transcript junction sequences for HipSci control samples.

- bowtie2_counts.py:  
  Counts the number of reads mapped to each fusion transcript for each sample.  
  Outputs the results in a pivot table format.

### Coverage Analysis and Genotyping
- coverage_analysis.sh:  
  Performs coverage analysis for WGS BAM files in specified regions.  
  Requires metadata and BED files.

- calculate_ratios.py:  
  Combines read counts from multiple VCF files to calculate ratios of read counts for specific regions.

- H1H2_SNP_genotyping.py:  
  Performs SNP genotyping for H1 and H2 cell lines using BAM files.

### Transcript Modeling and Visualization
- isoquant.py:  
  Constructs transcript models from PacBio reads supporting fusion transcripts.  

- FusionCatcher_processing.R:  
  Combines outputs from various scripts for comprehensive visualization.  

- variable_fusion_expression.R:  
  Focuses on visualizing fusion transcripts with inter-individual variation.


# Investigating Coding Potential of Fusion Transcripts through Ribosomal Profiling Data

## Ribo-seq analysis: 
- Step 1. This script processes a .fastq file containing Ribo-Seq (RS) reads, aligns the reads to the indexed transcripts using Bowtie2, converts SAM to BAM, sorts and indexes the BAM files, and extracts reads for each transcript. The results are saved in separate directories for SAM/BAM files and FASTA files for each fusion transcript and sample combination. Header names are assigned based on the gene name the RS reads matched to.
- Step 2. This script filter out RS reads that are not spanning the fusion junction of the FT they matched.
- Step 3. All the positive selected RS reads, are screened with BlastN, to find whether they match some known transcripts as well or they're uniquely matching the novel FT.
- Step 4. This script finally filter out all the RS reads that have matched some other known transcripts in the genome.

scripts used: 
- Ribo-seq_FT.py : whole workflow 
- Ribo-seq_pos_ctrl_BlastN.py : BlastN filetering step just for positive control sets 

## Generation of 10 random sets of 717 control transcripts spanning exon-exon junctions. 
- script 0: final_posctrl_creation_0.sh & final_posctrl_creation_python_to_run_with_0.py
#script 0 (bash+python): retrieve exons and introns coordinates of all genes from GENCODE_v44 file 

- script 1: final_posctrl_creation_1.sh
#script 1 (bash): shuffle and select 717 random exon junctions (43bp) and get fasta sequences
#after script 1 , to get ENSG  
- https://www.ensembl.org/biomart/martview/e3d502e53e5b353cf21726a897d51ac4 (tsv)
- NB: manually add names ENSG when missing

- script 2: final_posctrl_creation_2.py
#script 2 (python): match the ENSG name to the specific ENST name

- script 3: final_posctrl_creation_3.sh
#script 3 (bash): create fasta file with ENSG + sequence 

- script 4: final_posctrl_creation_4.py
#script 4 (python): generate final fasta file with the merged end (43bp) + start sequences (43bp) = 86bp



# scripts used to generate Figures 3E and 3F: 
- fig_3E_bubblechart.py
- fig_3E_piechart.py
- fig_3F_funnelchart.py






