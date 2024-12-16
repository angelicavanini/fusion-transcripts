# Clean up the R environment
rm(list=ls())

# Load necessary libraries
library(readxl)
library(writexl)
library(GenomicRanges)
#library(dplyr)
library(biomaRt)
library(karyoploteR)
library(tidyverse)
library(gplots)
library(ggplot2)
library(tibble)

# Set the working directory to the location where your files (of the type *final-list_candidate-fusion-genes.txt) are stored
setwd("/Users/colettemoses/Desktop/Good_scripts/fusioncatcher_new/final_output")

# Create an empty list to store the data frames
data_frames <- list()

# Get the list of file names from all files in the directory
file_names <- list.files()

# Read each file, add the additional column with the modified file name, and store it in the list
for (file in file_names) {
  df <- read.delim(file, header = TRUE, fill = TRUE)
  
  # Modify the file name
  modified_name <- sub("_final-list_candidate-fusion-genes.txt", "", file)
  
  # Add the modified file name as a new column
  df$FileName <- modified_name
  
  # Create FusionPair column
  df$FusionPair <- paste(df[,1], df[,2], sep = " ")
  
  # Create FusionID column
  df$FusionID <- paste(df$Gene_1_symbol.5end_fusion_partner., df$Gene_2_symbol.3end_fusion_partner., 
                       df$Fusion_point_for_gene_1.5end_fusion_partner., df$Fusion_point_for_gene_2.3end_fusion_partner., sep = "_")
  
  # Extract Strand information
  df$Strand.5end_partner <- substr(df$Fusion_point_for_gene_1.5end_fusion_partner., nchar(df$Fusion_point_for_gene_1.5end_fusion_partner.), nchar(df$Fusion_point_for_gene_1.5end_fusion_partner.))
  df$Strand.3end_partner <- substr(df$Fusion_point_for_gene_2.3end_fusion_partner., nchar(df$Fusion_point_for_gene_2.3end_fusion_partner.), nchar(df$Fusion_point_for_gene_2.3end_fusion_partner.))
  
  # Store the data frame in the list with the modified file name as the key
  data_frames[[modified_name]] <- df
}

# Combine all the data frames in the list into one
fusioncatcher_raw <- do.call(rbind, data_frames)

# Extract chromosome and position information
chrom <- as.data.frame(strsplit(fusioncatcher_raw$Fusion_point_for_gene_1.5end_fusion_partner., ":"))
chrom <- as.data.frame(t(chrom))
fusioncatcher_raw$Chromosome_1.5end_partner <- chrom[,1]
fusioncatcher_raw$Position_on_chromosome_1.5end_partner <- chrom[,2]

chrom <- as.data.frame(strsplit(fusioncatcher_raw$Fusion_point_for_gene_2.3end_fusion_partner., ":"))
chrom <- as.data.frame(t(chrom))
fusioncatcher_raw$Chromosome_2.3end_partner <- chrom[,1]
fusioncatcher_raw$Position_on_chromosome_2.3endpartner <- chrom[,2]

# Calculate the number of unique FusionIDs
n_fusion_ids <- length(unique(fusioncatcher_raw$FusionID))

# Print the number of unique FusionIDs
cat("Number of unique FusionIDs:", n_fusion_ids, "\n")

# Initialize an empty data frame to store descriptive information for each FusionID
descriptive_info_by_FusionID <- data.frame()

# Loop through each FusionID
for(i in 1:n_fusion_ids) {
  # Get the FusionID
  fusion_id <- unique(fusioncatcher_raw$FusionID)[i]
  
  # Find the indices of rows corresponding to the current FusionID
  ind <- which(fusioncatcher_raw$FusionID == fusion_id)
  
  # Extract descriptive information from the first row of the FusionID
  descriptive <- fusioncatcher_raw[ind[1], ]
  
  # Append the descriptive information to the new data frame
  descriptive_info_by_FusionID <- rbind(descriptive_info_by_FusionID, descriptive)
}

# Remove duplicate rows to get unique FusionIDs
collapsed_by_FusionID <- unique(descriptive_info_by_FusionID)

# Initialize vectors to store calculated values
present_in_n_samples <- vector()
file_names <- vector()
counts_sum <- vector()
counts_median <- vector()
counts_mean <- vector()
counts_sd <- vector()
pairs_sum <- vector()
pairs_median <- vector()
pairs_mean <- vector()
pairs_sd <- vector()
unique_reads_sum <- vector()
unique_reads_median <- vector()
unique_reads_mean <- vector()
unique_reads_sd <- vector()

# Loop through each unique FusionID
for(i in 1:nrow(collapsed_by_FusionID)){
  # Get the FusionID
  fusion_id <- collapsed_by_FusionID[i, "FusionID"]
  
  # Find the indices of rows corresponding to the current FusionID
  ind <- which(fusioncatcher_raw$FusionID == fusion_id)
  
  # Calculate the number of samples where the FusionID is present
  present_in_n_samples[i] <- length(unique(fusioncatcher_raw[ind, "FileName"]))
  
  # Get the names of files where the FusionID is present
  x <- unique(fusioncatcher_raw[ind, "FileName"])
  file_names[i] <- paste(x, collapse = ", ")
  
  # Calculate statistics for common mapping reads, spanning pairs, and spanning unique reads
  counts_sum[i] <- sum(fusioncatcher_raw[ind, "Counts_of_common_mapping_reads"])
  counts_median[i] <- median(fusioncatcher_raw[ind, "Counts_of_common_mapping_reads"])
  counts_mean[i] <- mean(fusioncatcher_raw[ind, "Counts_of_common_mapping_reads"])
  counts_sd[i] <- sd(fusioncatcher_raw[ind, "Counts_of_common_mapping_reads"], na.rm=TRUE)
  
  pairs_sum[i] <- sum(fusioncatcher_raw[ind, "Spanning_pairs"])
  pairs_median[i] <- median(fusioncatcher_raw[ind, "Spanning_pairs"])
  pairs_mean[i] <- mean(fusioncatcher_raw[ind, "Spanning_pairs"])
  pairs_sd[i] <- sd(fusioncatcher_raw[ind, "Spanning_pairs"],  na.rm=TRUE)
  
  unique_reads_sum[i] <- sum(fusioncatcher_raw[ind, "Spanning_unique_reads"])
  unique_reads_median[i] <- median(fusioncatcher_raw[ind, "Spanning_unique_reads"])
  unique_reads_mean[i] <- mean(fusioncatcher_raw[ind, "Spanning_unique_reads"])
  unique_reads_sd[i] <- sd(fusioncatcher_raw[ind, "Spanning_unique_reads"],  na.rm=TRUE)
}

# Calculate the total number of samples
n_samples <- length(unique(fusioncatcher_raw$FileName))

# Calculate the percentage of samples where each FusionID is present
percentage <- (present_in_n_samples / n_samples) * 100

# Combine calculated values with the collapsed data frame
collapsed_by_FusionID <- cbind(collapsed_by_FusionID, 
                               present_in_n_samples,
                               percentage,
                               file_names ,
                               round(counts_sum, 1),
                               round(counts_median, 1),
                               round(counts_mean, 1),
                               round(counts_sd, 1),
                               pairs_sum,
                               round(pairs_median, 1),
                               round(pairs_mean, 1),
                               round(pairs_sd, 1),
                               unique_reads_sum,
                               round(unique_reads_median, 1),
                               round(unique_reads_mean, 1),
                               round(unique_reads_sd, 1))

# Convert columns containing counts and reads to numeric
ind <- grep("counts|reads", names(collapsed_by_FusionID))
collapsed_by_FusionID[ind] <- lapply(collapsed_by_FusionID[ind], as.numeric)

# Define new column names
new_colnames <- c("Present in n samples",
                  "Present in % samples",
                  "File names",
                  "Common mapping reads sum",
                  "Common mapping reads median",
                  "Common mapping reads mean", 
                  "Common mapping reads sd", 
                  "Spanning pairs sum",
                  "Spanning pairs median", 
                  "Spanning pairs mean",
                  "Spanning pairs sd",
                  "Spanning unique reads sum",
                  "Spanning unique reads median", 
                  "Spanning unique reads mean",
                  "Spanning unique reads sd")

# Rename columns
colnames(collapsed_by_FusionID)[26:length(collapsed_by_FusionID)] <- new_colnames

# Delete the "FileName" column
collapsed_by_FusionID <- collapsed_by_FusionID[, !(names(collapsed_by_FusionID) %in% c("FileName"))]

# Add a column containing the fusion sequence with no asterisk (to be used for BLAT)
collapsed_by_FusionID$Fusion_sequence_no_asterisk <- gsub("\\*", "", collapsed_by_FusionID$Fusion_sequence)

# Properly format chromosomes and coordinates in fusion transcripts dataframe
collapsed_by_FusionID$Chromosome_1.5end_partner <- paste0("chr", collapsed_by_FusionID$Chromosome_1.5end_partner)
collapsed_by_FusionID$Chromosome_2.3end_partner <- paste0("chr", collapsed_by_FusionID$Chromosome_2.3end_partner)

collapsed_by_FusionID$Position_on_chromosome_1.5end_partner <- as.numeric(collapsed_by_FusionID$Position_on_chromosome_1.5end_partner)
collapsed_by_FusionID$Position_on_chromosome_2.3endpartner <- as.numeric(collapsed_by_FusionID$Position_on_chromosome_2.3endpartner)

# Sort the data frame by FusionID
collapsed_by_FusionID <- collapsed_by_FusionID[order(collapsed_by_FusionID$FusionID), ]
head(collapsed_by_FusionID)

# Add new columns indicating whether each fusion is on the same strand, intrachromosomal and occurs in expected direction of transcription
collapsed_by_FusionID$Intrachromosomal <- collapsed_by_FusionID$Chromosome_1.5end_partner == collapsed_by_FusionID$Chromosome_2.3end_partner

collapsed_by_FusionID$Same_strand <- ifelse(collapsed_by_FusionID$Intrachromosomal, 
                                            collapsed_by_FusionID$Strand.5end_partner == collapsed_by_FusionID$Strand.3end_partner,
                                            NA)

collapsed_by_FusionID$Expected_direction <- ifelse(collapsed_by_FusionID$Same_strand,
                                                   ifelse(collapsed_by_FusionID$Strand.5end_partner == "+",
                                                          collapsed_by_FusionID$Position_on_chromosome_1.5end_partner < collapsed_by_FusionID$Position_on_chromosome_2.3endpartner,
                                                          collapsed_by_FusionID$Position_on_chromosome_1.5end_partner > collapsed_by_FusionID$Position_on_chromosome_2.3endpartner),
                                                   NA)

## Filtering

# Filter out rows containing any of the 'red-flag' terms in the fusion_description column
# terms_to_remove <- c("similar_reads", "duplicates", "short_repeats", "long_repeats", "pair_pseudo_genes", "paralogs", "ambiguous")
collapsed_by_FusionID_filtered <- collapsed_by_FusionID[!grepl("similar_reads", collapsed_by_FusionID$Fusion_description, ignore.case = TRUE), ]
collapsed_by_FusionID_filtered <- collapsed_by_FusionID_filtered[!grepl("duplicates", collapsed_by_FusionID_filtered$Fusion_description, ignore.case = TRUE), ]
collapsed_by_FusionID_filtered <- collapsed_by_FusionID_filtered[!grepl("short_repeats", collapsed_by_FusionID_filtered$Fusion_description, ignore.case = TRUE), ]
collapsed_by_FusionID_filtered <- collapsed_by_FusionID_filtered[!grepl("long_repeats", collapsed_by_FusionID_filtered$Fusion_description, ignore.case = TRUE), ]
collapsed_by_FusionID_filtered <- collapsed_by_FusionID_filtered[!grepl("pair_pseudo_genes", collapsed_by_FusionID_filtered$Fusion_description, ignore.case = TRUE), ]
collapsed_by_FusionID_filtered <- collapsed_by_FusionID_filtered[!grepl("paralogs", collapsed_by_FusionID_filtered$Fusion_description, ignore.case = TRUE), ]
collapsed_by_FusionID_filtered <- collapsed_by_FusionID_filtered[!grepl("ambiguous", collapsed_by_FusionID_filtered$Fusion_description, ignore.case = TRUE), ]

# Filter out any rows with high number of common mapping reads
common_mapping <- collapsed_by_FusionID_filtered$`Common mapping reads sum` <= 10
collapsed_by_FusionID_filtered <- collapsed_by_FusionID_filtered[common_mapping, ]

## Filter out fusion transcripts intersecting with RepeatMasker

# Read the repeatmasker bed file and convert coordinates to GRanges object
repeatmasker_bed <- read.table("/Users/colettemoses/Desktop/Good_scripts/FusionCatcher/hg38_RepeatMasker.bed", header = FALSE)
colnames(repeatmasker_bed) <- c("chrom", "start", "end")

repeatmasker_ranges <- GRanges(seqnames = repeatmasker_bed$chrom,
                               ranges = IRanges(start = repeatmasker_bed$start,
                                                end = repeatmasker_bed$end))

# Create GRanges objects for fusion breakpoints with ±1 bp extension
fusion_breakpoints_5end <- GRanges(
  seqnames = collapsed_by_FusionID_filtered$Chromosome_1.5end_partner,
  ranges = IRanges(
    start = collapsed_by_FusionID_filtered$Position_on_chromosome_1.5end_partner - 1,
    end = collapsed_by_FusionID_filtered$Position_on_chromosome_1.5end_partner + 1
  )
)

fusion_breakpoints_3end <- GRanges(
  seqnames = collapsed_by_FusionID_filtered$Chromosome_2.3end_partner,
  ranges = IRanges(
    start = collapsed_by_FusionID_filtered$Position_on_chromosome_2.3endpartner - 1,
    end = collapsed_by_FusionID_filtered$Position_on_chromosome_2.3endpartner + 1
  )
)

# Intersect with repeatmasker track
overlap_5end_rmsk <- findOverlaps(fusion_breakpoints_5end, repeatmasker_ranges)
overlap_3end_rmsk <- findOverlaps(fusion_breakpoints_3end, repeatmasker_ranges)

# Get the indices of fusion transcripts where either end overlaps with repeatmasker
intersecting_indices_rmsk <- unique(c(queryHits(overlap_5end_rmsk), queryHits(overlap_3end_rmsk)))

# Filter out fusion transcripts with overlapping ends
collapsed_by_FusionID_filtered_rmsk <- collapsed_by_FusionID_filtered[-intersecting_indices_rmsk, ]

## Compare the percent of intra- and inter-chromosomal fusions before and after filtering

# Calculate the number of rows where Chromosome_1.5end_partner equals Chromosome_2.3end_partner
intrachromosomal_count_before <- sum(collapsed_by_FusionID$Intrachromosomal)
intrachromosomal_count_after <- sum(collapsed_by_FusionID_filtered_rmsk$Intrachromosomal)

# Calculate the number of rows where Chromosome_1.5end_partner does not equal Chromosome_2.3end_partner
interchromosomal_count_before <- nrow(collapsed_by_FusionID) - intrachromosomal_count_before
interchromosomal_count_after <- nrow(collapsed_by_FusionID_filtered_rmsk) - intrachromosomal_count_after

# Calculate the percentage of rows where Chromosome_1.5end_partner equals Chromosome_2.3end_partner
intrachromosomal_percent_before <- (intrachromosomal_count_before / nrow(collapsed_by_FusionID)) * 100
intrachromosomal_percent_after <- (intrachromosomal_count_after / nrow(collapsed_by_FusionID_filtered_rmsk)) * 100

# Calculate the percentage of rows where Chromosome_1.5end_partner does not equal Chromosome_2.3end_partner
interchromosomal_percent_before <- (interchromosomal_count_before / nrow(collapsed_by_FusionID)) * 100
interchromosomal_percent_after <- (interchromosomal_count_after / nrow(collapsed_by_FusionID_filtered_rmsk)) * 100

# Print the results
cat("Before filtering:\n")
cat("Intrachromosomal fusions:", intrachromosomal_count_before, "(", intrachromosomal_percent_before, "%)\n")
cat("Interchromosomal fusions:", interchromosomal_count_before, "(", interchromosomal_percent_before, "%)\n\n")

cat("After filtering:\n")
cat("Intrachromosomal fusions:", intrachromosomal_count_after, "(", intrachromosomal_percent_after, "%)\n")
cat("Interchromosomal fusions:", interchromosomal_count_after, "(", interchromosomal_percent_after, "%)\n")

## Perform intersection of fusion transcripts and segmental duplications

# Read the seg dups bed file and convert coordinates to GRanges object
segdups_bed <- read.table("/Users/colettemoses/Desktop/Good_scripts/FusionCatcher/hg38_genomicSuperDups.bed", header = FALSE)
colnames(segdups_bed) <- c("chrom", "start", "end")

segdups_ranges <- GRanges(seqnames = segdups_bed$chrom,
                          ranges = IRanges(start = segdups_bed$start,
                                           end = segdups_bed$end))

# Create GRanges objects for fusion breakpoints with ±1 bp extension
fusion_breakpoints_5end <- GRanges(
  seqnames = collapsed_by_FusionID_filtered_rmsk$Chromosome_1.5end_partner,
  ranges = IRanges(
    start = collapsed_by_FusionID_filtered_rmsk$Position_on_chromosome_1.5end_partner - 1,
    end = collapsed_by_FusionID_filtered_rmsk$Position_on_chromosome_1.5end_partner + 1
  )
)

fusion_breakpoints_3end <- GRanges(
  seqnames = collapsed_by_FusionID_filtered_rmsk$Chromosome_2.3end_partner,
  ranges = IRanges(
    start = collapsed_by_FusionID_filtered_rmsk$Position_on_chromosome_2.3endpartner - 1,
    end = collapsed_by_FusionID_filtered_rmsk$Position_on_chromosome_2.3endpartner + 1
  )
)

# Intersect with seg dups track
overlap_5end_dups <- findOverlaps(fusion_breakpoints_5end, segdups_ranges)
overlap_3end_dups <- findOverlaps(fusion_breakpoints_3end, segdups_ranges)

# Get the indices of fusion transcripts where either end overlaps with seg dups
intersecting_indices_dups <- unique(c(queryHits(overlap_5end_dups), queryHits(overlap_3end_dups)))

# Initialize a new column with FALSE for all entries
collapsed_by_FusionID_filtered_rmsk$Overlap_segdups <- FALSE

# Mark rows with intersecting indices as TRUE in the new column
collapsed_by_FusionID_filtered_rmsk$Overlap_segdups[intersecting_indices_dups] <- TRUE

## Filter results to find entries that could not be originating from readthrough according to reference genome structure

# Remove fusions flagged as adjacent and oriented in the expected direction of transcription
subset_to_remove <- collapsed_by_FusionID_filtered_rmsk[
  grepl("adjacent", collapsed_by_FusionID_filtered_rmsk$Fusion_description) & collapsed_by_FusionID_filtered_rmsk$Expected_direction,
]
collapsed_by_FusionID_nonreadthrough <- collapsed_by_FusionID_filtered_rmsk[!rownames(collapsed_by_FusionID_filtered_rmsk) %in% rownames(subset_to_remove), ]

# Select fusions on the opposite strand, not oriented in the expected direction of transcription, or with a fusion distance > 500 kb
nonRT_condition <- collapsed_by_FusionID_nonreadthrough$Same_strand == FALSE |
  (collapsed_by_FusionID_nonreadthrough$Same_strand == TRUE & collapsed_by_FusionID_nonreadthrough$Expected_direction == FALSE) |
  abs(collapsed_by_FusionID_nonreadthrough$Position_on_chromosome_1.5end_partner - collapsed_by_FusionID_nonreadthrough$Position_on_chromosome_2.3endpartner) > 500000

# Select only intrachromosomal fusions
nonRT_condition <- nonRT_condition & collapsed_by_FusionID_nonreadthrough$Intrachromosomal == TRUE

collapsed_by_FusionID_nonreadthrough <- collapsed_by_FusionID_nonreadthrough[nonRT_condition, ]

############################ RUN FROM HERE TO REPEAT WITH NEW BOWTIE SEARCH RESULTS ############################ 

## Remove fusion transcripts detected in only one individual, and collapse by FusionPair (optional)
collapsed_by_FusionID_filtered_rmsk_2ormore <- collapsed_by_FusionID_filtered_rmsk[collapsed_by_FusionID_filtered_rmsk$`Present in n samples` > 1, ]
collapsed_by_FusionID_nonreadthrough_2ormore <- collapsed_by_FusionID_nonreadthrough[collapsed_by_FusionID_nonreadthrough$`Present in n samples` > 1, ]

collapsed_by_FusionPair_2ormore <- collapsed_by_FusionID_filtered_rmsk_2ormore %>%
  distinct(FusionPair, .keep_all = TRUE)
collapsed_by_FusionPair_nonreadthrough_2ormore <- collapsed_by_FusionID_nonreadthrough_2ormore %>%
  distinct(FusionPair, .keep_all = TRUE)

## Adding in bowtie2 results for all 660 fusions

bowtie2_660_results_human <- read.table("/Users/colettemoses/Desktop/Good_scripts/bowtie2/Counts_pivot_2024-06-28_10-54-50_mayo717.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
bowtie2_660_results_hipsci <- read.table("/Users/colettemoses/Desktop/Good_scripts/bowtie2/Combined_Counts_pivot_2024-07-02_07-55-58.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# bowtie2_660_results_comparative1 <- read.table("/Users/colettemoses/Desktop/Good_scripts/bowtie2/Counts_pivot_2024-07-02_07-31-48_comparative_717_newersamples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# bowtie2_660_results_comparative2 <- read.table("/Users/colettemoses/Desktop/Good_scripts/bowtie2/Counts_pivot_2024-07-02_11-06-48_comparative_717_oldersamples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# bowtie2_660_results_comparative <- full_join(bowtie2_660_results_comparative1, bowtie2_660_results_comparative2, by = "Fusion_Transcript")
bowtie2_660_results_comparative <- read.table("/Users/colettemoses/Desktop/Good_scripts/bowtie2/Counts_pivot_2024-07-04_15-53-21_comparative_171_multiscore.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Choose (uncomment) one of the 3 options below:

# 1
# Merge the datasets, first checking all fusions are identical
bowtie2_660_results_comparative <- subset(bowtie2_660_results_comparative, Fusion_Transcript != "KANSL1_LRRC37A3_ex5")
identical(bowtie2_660_results_human$Fusion_Transcript, bowtie2_660_results_comparative$Fusion_Transcript)
# # Columns present in bowtie2_660_results_human but not in bowtie2_660_results_comparative
# missing_in_comparative <- setdiff(colnames(bowtie2_660_results_human), colnames(bowtie2_660_results_comparative))
# print(missing_in_comparative)
# # Columns present in bowtie2_660_results_comparative but not in bowtie2_660_results_human
# missing_in_human <- setdiff(colnames(bowtie2_660_results_comparative), colnames(bowtie2_660_results_human))
# print(missing_in_human)
# bowtie2_660_results_comparative <- bowtie2_660_results_comparative[, -c(662, 663)]
# identical(colnames(bowtie2_660_results_human), colnames(bowtie2_660_results_comparative))
# Perform join
bowtie2_660_results <- full_join(bowtie2_660_results_human, bowtie2_660_results_comparative, by = "Fusion_Transcript")

# # 2
bowtie2_660_results <- (bowtie2_660_results_comparative)

# 3
bowtie2_660_results <- (bowtie2_660_results_human)
# OR
bowtie2_660_results <- (bowtie2_660_results_hipsci)

colnames(bowtie2_660_results) <- sub("X", "", colnames(bowtie2_660_results)) # Mayo only
colnames(bowtie2_660_results)[1] <- "Fusion_Transcript_Name"
# To remove "KANSL1_ex5" if necessary
# bowtie2_660_results <- bowtie2_660_results[-334,]

# Convert columns to numeric where necessary (except Fusion_Transcript_Name)
bowtie2_660_results[, -1] <- lapply(bowtie2_660_results[, -1], as.numeric)

# Calculate counts > 0, percentages and summary statistics
bowtie2_660_results_summarised <- bowtie2_660_results %>%
  rowwise() %>%
  mutate(
    bowtie2_present_in_n_samples = sum(c_across(-Fusion_Transcript_Name) > 0),
    bowtie2_present_in_percent_samples = (bowtie2_present_in_n_samples / (ncol(.) - 1)) * 100,
    bowtie2_sum = rowSums(across(-c(Fusion_Transcript_Name, starts_with("bowtie2")))),
    bowtie2_median = apply(across(-c(Fusion_Transcript_Name, starts_with("bowtie2"))), 1, function(x) ifelse(sum(x > 0) > 0, median(x[x > 0]), NA)),
    bowtie2_mean = apply(across(-c(Fusion_Transcript_Name, starts_with("bowtie2"))), 1, function(x) ifelse(sum(x > 0) > 0, mean(x[x > 0]), NA)),
    bowtie2_SD = apply(across(-c(Fusion_Transcript_Name, starts_with("bowtie2"))), 1, function(x) ifelse(sum(x > 0) > 0, sd(x[x > 0]), NA))
  ) %>%
  ungroup() %>%
  dplyr::select(Fusion_Transcript_Name, bowtie2_present_in_n_samples, bowtie2_present_in_percent_samples, bowtie2_sum, bowtie2_median, bowtie2_mean, bowtie2_SD)

# Calculate the number of rows where the mean is 0
row_means <- rowMeans(bowtie2_660_results[, -1], na.rm = TRUE)
num_rows_with_zero_mean <- sum(row_means == 0)
print(num_rows_with_zero_mean)

# Join the bowtie2 results onto original FusionCatcher results
collapsed_by_FusionID_filtered_rmsk_2ormore <- left_join(collapsed_by_FusionID_filtered_rmsk_2ormore, bowtie2_660_results_summarised,
                    by = c("FusionID" = "Fusion_Transcript_Name"))

collapsed_by_FusionID_nonreadthrough_2ormore <- left_join(collapsed_by_FusionID_nonreadthrough_2ormore, bowtie2_660_results_summarised,
                                                         by = c("FusionID" = "Fusion_Transcript_Name"))

# Make sure to discount the summary statistics in collapsed_by_FusionID_filtered_rmsk_2ormore and collapsed_by_FusionID_nonreadthrough_2ormore
# if you have merged both the human and cross-species data

## Write the data frames to Excel files

original_file_path <- "/Users/colettemoses/Desktop/Good_scripts/fusioncatcher_new/Mayo_FusionCatcher_compiled_2ormore_v2.xlsx"
write_xlsx(list(collapsed_by_FusionID_filtered_rmsk_2ormore), original_file_path, col_names = TRUE)

# nonRT_file_path <- "/Users/colettemoses/Desktop/Good_scripts/FusionCatcher/Mayo_FusionCatcher_nonreadthrough.xlsx"
# write_xlsx(list(collapsed_by_FusionID_nonreadthrough), nonRT_file_path, col_names = TRUE)

nonRT_2ormore_file_path <- "/Users/colettemoses/Desktop/Good_scripts/FusionCatcher/Mayo_FusionCatcher_nonreadthrough_2ormore_plusbowtie2.xlsx"
write_xlsx(list(collapsed_by_FusionID_nonreadthrough_2ormore), nonRT_2ormore_file_path, col_names = TRUE)

# # Generate FASTA file of fusion sequences to be used for BLAT
# 
# writeFASTA <- function(file, fusion_data, append = FALSE) {
#   mode <- ifelse(append, "a", "w")  # Determine file writing mode
#   
#   # Open the file in the specified mode
#   con <- file(file, open = mode)
#   
#   # Write FusionID and sequence for each row
#   for (i in 1:nrow(fusion_data)) {
#     cat(">", fusion_data$FusionID[i], "\n", fusion_data$Fusion_sequence_no_asterisk[i], "\n", file = con, append = TRUE, sep = "")
#   }
#   
#   # Close the file connection
#   close(con)
# }
# 
# writeFASTA("/Users/colettemoses/Desktop/Luca_Mayo_FusionCatcher/fusion_junctions_BLAT.fasta", collapsed_by_FusionID_nonreadthrough_2ormore, append = TRUE)

########################################################################################################################

## Figure 1 karyoplot

# Preparing data for karyoplot

intrachromosomal_only <- collapsed_by_FusionID_filtered_rmsk_2ormore[collapsed_by_FusionID_filtered_rmsk_2ormore$Intrachromosomal == TRUE, ]
intrachromosomal_only <- collapsed_by_FusionID_filtered_rmsk_2ormore

fusionlinks1 <- intrachromosomal_only[, c(
  "Chromosome_1.5end_partner",
  "Position_on_chromosome_1.5end_partner",
  "Strand.5end_partner",
  "bowtie2_present_in_n_samples",
  "bowtie2_mean",
  "FusionID"
)]

# Add a new column for Position_on_chromosome_1.5end_partner_plus_1
fusionlinks1$Position_on_chromosome_1.5end_partner_plus_1 <- fusionlinks1$Position_on_chromosome_1.5end_partner + 1

# Create fusionlinks2 dataframe
fusionlinks2 <- intrachromosomal_only[, c(
  "Chromosome_2.3end_partner",
  "Position_on_chromosome_2.3endpartner",
  "Strand.3end_partner",
  "bowtie2_present_in_n_samples",
  "bowtie2_mean",
  "FusionID"
)]

# Add a new column for Position_on_chromosome_2.3end_partner_plus_1
fusionlinks2$Position_on_chromosome_2.3endpartner_plus_1 <- fusionlinks2$Position_on_chromosome_2.3endpartner + 1

# Rename columns
fusionlinks1 <- setNames(fusionlinks1, c("chr", "start", "strand", "n_individuals", "bowtie2_reads_mean", "FusionID", "end"))
fusionlinks2 <- setNames(fusionlinks2, c("chr", "start", "strand", "n_individuals", "bowtie2_reads_mean", "FusionID", "end"))

# Reorder columns
fusionlinks1 <- fusionlinks1 %>%
  dplyr::select(chr, start, end, strand, n_individuals, bowtie2_reads_mean, FusionID)

fusionlinks2 <- fusionlinks2 %>%
  dplyr::select(chr, start, end, strand, n_individuals, bowtie2_reads_mean, FusionID)

# Convert data to GRanges objects
links1 <- toGRanges(fusionlinks1)
links2 <- toGRanges(fusionlinks2)

# Visualize links
links.col <- colByValue(links1$n_individuals, colors=c("#2E72AD20", "#025F99"))
pp <- getDefaultPlotParams(plot.type=1)
kp <- plotKaryotype(genome="hg38", plot.type=1, plot.params=pp, ideogram.plotter=NULL)
kpAddCytobandsAsLine(kp)
kpPlotLinks(kp, data=links1, data2=links2, col=links.col, border=links.col, r0=-0.1, r1=1.1)

# Export to pdf
pdf("/Users/colettemoses/Desktop/karyoplot_bowtie2_mayo_717_all.pdf")
kp <- plotKaryotype(genome="hg38", plot.type=1, plot.params=pp, ideogram.plotter=NULL)
kpAddCytobandsAsLine(kp)
kpPlotLinks(kp, data=links1, data2=links2, col=links.col, border=links.col, r0=-0.1, r1=1.1)
dev.off()

# pdf("/Users/colettemoses/Desktop/karyoplot_chr17.pdf")
# kp <- plotKaryotype(genome="hg38", plot.type=1, plot.params=pp, ideogram.plotter=NULL, chromosomes="chr17")
# kpAddCytobandsAsLine(kp)
# kpPlotLinks(kp, data=links1, data2=links2, col=links.col, border=links.col, r0=-0.1, r1=1.1)
# dev.off()

# ## Fusion locations (to import into Illustrator)
# 
# # Import coordinates from Excel
# genestocoords <- as.data.frame(read_excel("/Users/colettemoses/Desktop/Work/EMBO-ENW_project/Dry_lab/fusion_analysis/R_scripts_fusioncatcher_mayo/genestocoord_consolidated.xlsx"))
# 
# # Extract chromosome, position, and fusion pairs
# chrs <- genestocoords$Chromosome_1.5end_partner
# pos <- genestocoords$Position_on_chromosome_1.5end_partner
# fusion <- genestocoords$FusionPair
# markers <- data.frame(chr=chrs, pos=pos, labels=fusion)
# 
# # Plotting karyotypes for specific chromosomes
# chromosomes_to_plot <- c("chr5", "chr16", "chr17", "chr22")
# for (chromosome in chromosomes_to_plot) {
#   kp <- plotKaryotype(genome="hg38", chromosomes=chromosome, plot.type=1, plot.params=pp)
#   kpPlotMarkers(kp, chr=chrs, x=pos, labels=fusion, adjust.label.position=FALSE)
# }

########################################################################################################################

# ## Figure 1 heatmaps
# 
# ## All fusions: FusionCatcher results
# 
# # Filter fusioncatcher_raw based on FusionID
# fusioncatcher_raw_filtered <- fusioncatcher_raw[fusioncatcher_raw$FusionID %in% collapsed_by_FusionID_filtered_rmsk_2ormore$FusionID, ]
# 
# # Convert values to numeric
# fusioncatcher_raw_filtered$Position_on_chromosome_1.5end_partner <- as.numeric(fusioncatcher_raw_filtered$Position_on_chromosome_1.5end_partner)
# fusioncatcher_raw_filtered$Spanning_unique_reads <- as.numeric(fusioncatcher_raw_filtered$Spanning_unique_reads)
# fusioncatcher_raw_filtered$Chromosome_1.5end_partner_numeric <- as.numeric(fusioncatcher_raw_filtered$Chromosome_1.5end_partner)
# 
# # Reformat FusionID to be more readable
# fusioncatcher_raw_filtered$FusionID_reformatted <- paste(fusioncatcher_raw_filtered$Gene_1_symbol.5end_fusion_partner,
#                                                "-", 
#                                                fusioncatcher_raw_filtered$Gene_2_symbol.3end_fusion_partner,
#                                                " chr",
#                                                fusioncatcher_raw_filtered$Chromosome_1.5end_partner,
#                                                ":",
#                                                fusioncatcher_raw_filtered$Position_on_chromosome_1.5end_partner,
#                                                "-",
#                                                fusioncatcher_raw_filtered$Position_on_chromosome_2.3endpartner,
#                                                sep = "")
# 
# # Function to create heatmap
# create_heatmap <- function(data, name){
#   # Order data by chromosome position
#   data <- data[order(data$Chromosome_1.5end_partner_numeric, data$Position_on_chromosome_1.5end_partner), ]
#   
#   # Select necessary columns
#   short <- data %>%
#     select(FileName, FusionID_reformatted, Spanning_unique_reads)
#   
#   # Create pivot table
#   pivot_data <- pivot_wider(short, names_from = FusionID_reformatted, values_from = Spanning_unique_reads, values_fill = 0)
#   
#   # Print the pivot table
#   print(pivot_data)
#   
#   # Convert pivot data to matrix and transpose
#   heatmap_matrix <- t(as.matrix(pivot_data[, -1]))  # Exclude FileName column
#   
#   # Set the color palette for the heatmap
#   heatmap_colors <- colorRampPalette(c("#C8C8C880", "#025F99"))(100)
#   
#   # Extract column names from transposed matrix (files)
#   labCol <- colnames(heatmap_matrix)
#   
#   # Extract row names from transposed matrix (fusion transcripts)
#   labRow <- rownames(heatmap_matrix)
#   
#   # Plot the heatmap using heatmap.2 with the updated labCol and labRow
#   heatmap.2(heatmap_matrix,
#             col = heatmap_colors,
#             scale = "none",
#             cexCol = 0.1,
#             cexRow = 0.2, 
#             Rowv = FALSE,
#             Colv = TRUE,
#             labCol = rep(" ", ncol(heatmap_matrix)),
#             labRow = labRow,
#             key = FALSE, 
#             trace = "none", 
#             dendrogram = "none", 
#             main = name,
#             cex = 3)
# }
# 
# # Loop over each chromosome
# for (chr in c(1:22, "X")) {
#   # Filter fusion data for the current chromosome
#   filtered_data <- fusioncatcher_raw_filtered %>%
#     filter(Chromosome_1.5end_partner == chr & Chromosome_2.3end_partner == chr)
#   
#   # Create heatmap for the current chromosome
#   create_heatmap(filtered_data, paste("Fusion transcripts\nchromosome", chr))
#   
#   # Output heatmap to PDF
#   pdf(paste0("/Users/colettemoses/Desktop/heatmaps/heatmap_chr_", chr, ".pdf"))
#   create_heatmap(filtered_data, paste("Fusion transcripts\nchromosome", chr))
#   dev.off()
# }
# 
# # Create heatmap for transcripts on all chromosomes
# create_heatmap(fusioncatcher_raw_filtered, paste("Fusion transcripts\nall chromosomes"))
# 
# # Output heatmap to PDF
# pdf("/Users/colettemoses/Desktop/heatmaps/heatmap_all.pdf", width = 8, height = 40)
# create_heatmap(fusioncatcher_raw_filtered, paste("Fusion transcripts\nall chromosomes"))
# dev.off()
# 
# ## 'Non-readthrough' fusions
# 
# # Filter fusioncatcher_raw based on FusionID from non-readthrough
# fusioncatcher_raw_filtered_nonreadthrough <- fusioncatcher_raw[fusioncatcher_raw$FusionID %in% collapsed_by_FusionID_nonreadthrough_2ormore$FusionID, ]
# 
# # Convert values to numeric
# fusioncatcher_raw_filtered_nonreadthrough$Position_on_chromosome_1.5end_partner <- as.numeric(fusioncatcher_raw_filtered_nonreadthrough$Position_on_chromosome_1.5end_partner)
# fusioncatcher_raw_filtered_nonreadthrough$Spanning_unique_reads <- as.numeric(fusioncatcher_raw_filtered_nonreadthrough$Spanning_unique_reads)
# fusioncatcher_raw_filtered_nonreadthrough$Chromosome_1.5end_partner_numeric <- as.numeric(fusioncatcher_raw_filtered_nonreadthrough$Chromosome_1.5end_partner)
# 
# # Reformat FusionID to be more readable
# fusioncatcher_raw_filtered_nonreadthrough$FusionID_reformatted <- paste(fusioncatcher_raw_filtered_nonreadthrough$Gene_1_symbol.5end_fusion_partner,
#                                                          "-", 
#                                                          fusioncatcher_raw_filtered_nonreadthrough$Gene_2_symbol.3end_fusion_partner,
#                                                          " chr",
#                                                          fusioncatcher_raw_filtered_nonreadthrough$Chromosome_1.5end_partner,
#                                                          ":",
#                                                          fusioncatcher_raw_filtered_nonreadthrough$Position_on_chromosome_1.5end_partner,
#                                                          "-",
#                                                          fusioncatcher_raw_filtered_nonreadthrough$Position_on_chromosome_2.3endpartner,
#                                                          sep = "")
# 
# # Create heatmap for non-readthrough transcripts on all chromosomes
# create_heatmap(fusioncatcher_raw_filtered_nonreadthrough, paste("Fusion transcripts\nnon-readthrough\nall chromosomes"))
# 
# # Output heatmap to PDF
# pdf("/Users/colettemoses/Desktop/heatmaps/heatmap_all_nonreadthrough.pdf")
# create_heatmap(fusioncatcher_raw_filtered_nonreadthrough, paste("Fusion transcripts\nnon-readthrough\nall chromosomes"))
# dev.off()

#################################################

## All fusions: bowtie2 results

# Reformat FusionID to be more readable
collapsed_by_FusionID_filtered_rmsk_2ormore$FusionID_reformatted <- paste(collapsed_by_FusionID_filtered_rmsk_2ormore$Gene_1_symbol.5end_fusion_partner,
                                                         "-", 
                                                         collapsed_by_FusionID_filtered_rmsk_2ormore$Gene_2_symbol.3end_fusion_partner,
                                                         " ",
                                                         collapsed_by_FusionID_filtered_rmsk_2ormore$Chromosome_1.5end_partner,
                                                         ":",
                                                         collapsed_by_FusionID_filtered_rmsk_2ormore$Position_on_chromosome_1.5end_partner,
                                                         "-",
                                                         collapsed_by_FusionID_filtered_rmsk_2ormore$Position_on_chromosome_2.3endpartner,
                                                         sep = "")

# Join the individual bowtie2 results with chromosomes and coordinates
bowtie2_heatmap_data <- bowtie2_660_results %>%
  left_join(collapsed_by_FusionID_filtered_rmsk_2ormore, by = c("Fusion_Transcript_Name" = "FusionID"))
bowtie2_heatmap_data <- bowtie2_heatmap_data %>%
  mutate(Fusion_Transcript_Name = FusionID_reformatted) 

# For Mayo only:
bowtie2_heatmap_data <- bowtie2_heatmap_data %>%
  dplyr::select(-c(278:296, 301:327)) # Change if necessary
bowtie2_heatmap_data <- bowtie2_heatmap_data[, c(1, 278:281, 2:277)]

# For hipsci only:
bowtie2_heatmap_data <- bowtie2_heatmap_data %>%
  dplyr::select(-c(193:211, 216:242)) # Change if necessary
bowtie2_heatmap_data <- bowtie2_heatmap_data[, c(1, 193:196, 2:192)]

# For comparative only:
bowtie2_heatmap_data <- bowtie2_heatmap_data %>%
  dplyr::select(-c(42:60, 65:91))
bowtie2_heatmap_data <- bowtie2_heatmap_data[, c(1, 42:45, 2:41)]

# For combined Mayo and comparative:
bowtie2_heatmap_data <- bowtie2_heatmap_data %>%
  dplyr::select(-c(311:329, 334:360))
bowtie2_heatmap_data <- bowtie2_heatmap_data[, c(1, 311:314, 2:310)]

# Define the color palette function with a logarithmic scale
heatmap_colors <- colorRampPalette(c("#C8C8C830", "#025F99"))

# Vector of chromosome names
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                 "chr20", "chr21", "chr22", "chrX", "chrY")

# # Function to create heatmap for a specific chromosome or all chromosomes with logarithmic color scaling
# create_heatmap_bowtie2 <- function(data, chromosome = NULL, name) {
#   if (!is.null(chromosome)) {
#     # Filter data for the specified chromosome
#     filtered_data <- data %>%
#       filter(Chromosome_1.5end_partner == chromosome & Chromosome_2.3end_partner == chromosome)
#   } else {
#     # Use all data if chromosome is not specified
#     filtered_data <- data
#   }
#   
#   if (nrow(filtered_data) > 0) {
#     # Order filtered data by Fusion_Transcript_Name
#     filtered_data <- filtered_data[order(filtered_data$Position_on_chromosome_1.5end_partner), ]
#     
#     # Convert filtered data to matrix
#     heatmap_matrix <- as.matrix(filtered_data[, -c(1:5)])  # Exclude metadata columns
#     
#     # Apply a logarithmic transformation to the data
#     heatmap_matrix_log <- log(heatmap_matrix + 1)  # Add 1 to avoid log(0)
#     
#     # Set color palette for the heatmap with logarithmic scale
#     heatmap_colors_log <- heatmap_colors(100)
#     
#     # Extract column names from transposed matrix (samples)
#     labCol <- colnames(heatmap_matrix)
#     
#     # Extract row names from transposed matrix (fusion transcripts)
#     labRow <- filtered_data$Fusion_Transcript_Name
#     
#     # Plot the heatmap using heatmap.2 with updated labCol and labRow
#     heatmap.2(heatmap_matrix_log,
#               col = heatmap_colors_log,
#               scale = "none",
#               cexCol = 0.1,
#               cexRow = 0.2, 
#               Rowv = FALSE,
#               Colv = TRUE,
#               labCol = rep(" ", ncol(heatmap_matrix)),
#               labRow = labRow,
#               key = FALSE, 
#               trace = "none", 
#               dendrogram = "none", 
#               main = name,
#               cex = 3)
#   } else {
#     cat("No data found for chromosome", chromosome, "\n")
#   }
# }

# To disable automatic clustering of samples (for comparative analysis)
# Function to create heatmap for a specific chromosome with logarithmic color scaling
create_heatmap_bowtie2 <- function(data, chromosome = NULL, name, column_order = NULL, column_labels = NULL) {
  if (!is.null(chromosome)) {
    # Filter data for the specified chromosome
    filtered_data <- data %>%
      filter(Chromosome_1.5end_partner == chromosome & Chromosome_2.3end_partner == chromosome)
  } else {
    # Use all data if chromosome is not specified
    filtered_data <- data
  }

  if (nrow(filtered_data) > 0) {
    # Order filtered data by Position_on_chromosome_1.5end_partner
    filtered_data <- filtered_data[order(filtered_data$Position_on_chromosome_1.5end_partner), ]

    # Convert filtered data to matrix
    heatmap_matrix <- as.matrix(filtered_data[, -c(1:5)])  # Exclude metadata columns

    # Apply a logarithmic transformation to the data
    heatmap_matrix_log <- log(heatmap_matrix + 1)  # Add 1 to avoid log(0)

    # Reorder the columns if a specific order is provided
    if (!is.null(column_order)) {
      heatmap_matrix_log <- heatmap_matrix_log[, column_order, drop = FALSE]
    }

    # Set color palette for the heatmap with logarithmic scale
    heatmap_colors_log <- heatmap_colors(100)

    # Extract row names from filtered data (fusion transcripts)
    labRow <- filtered_data$Fusion_Transcript_Name

    # Ensure column_labels length matches number of columns
    if (!is.null(column_labels) && length(column_labels) == ncol(heatmap_matrix_log)) {
      labCol <- column_labels
    } else {
      labCol <- colnames(heatmap_matrix_log)  # Default to column names if labels mismatch
    }

    # Plot the heatmap using heatmap.2 with updated labCol and labRow
    heatmap.2(heatmap_matrix_log,
              col = heatmap_colors_log,
              scale = "none",
              cexCol = 0.7,  # Adjust for better visibility
              cexRow = 0.4,
              Rowv = FALSE,
              Colv = FALSE,  # Disable column clustering
              labCol = labCol,
              labRow = labRow,
              key = FALSE,
              trace = "none",
              dendrogram = "none",
              main = name,
              cex = 1.5)
  } else {
    cat("No data found for chromosome", chromosome, "\n")
  }
}

# Specify the order of the columns
column_order <- colnames(bowtie2_heatmap_data)
column_order <- column_order[-(1:5)]
# desired_order <- c("SRR1758917", "SRR5804451", "SRR5804456", "SRR5804464", "SRR5804471", "SRR5804478", "SRR5804486", "SRR6334913", 
#                    "SRR6334974", "SRR6334986", "SRR6335020", "SRR6335047", "SRR2040584", "SRR2040585", "ERR3474012", "ERR3474022", 
#                    "SRR17380399", "SRR17380400", "SRR17380401", "SRR17380402", "SRR6334947", "SRR6334962", "SRR6334989", "SRR6334996", 
#                    "SRR6335010", "SRR2040592", "ERR3474062", "ERR3474072", "ERR3474082", "SRR17380395", "SRR17380396", "SRR17380397", 
#                    "SRR17380398", "SRR2040596", "SRR2040597", "SRR2040598", "SRR27109800", "SRR27109801", "SRR27109804", "SRR27109805")
# column_order <- desired_order
# first_part <- column_order[1:276]
# last_part <- column_order[277:316]
# if (all(sort(last_part) == sort(desired_order))) {
#   # Reorder last_part to match desired_order
#   reordered_last_part <- desired_order[match(desired_order, last_part)]
#   
#   # Combine the first 276 elements with the reordered last 40 elements
#   column_order <- c(first_part, reordered_last_part)
#   
#   # Print or return the new column order
#   print(column_order)
# } else {
#   print("The last 40 elements of column_order do not match the elements in desired_order.")
# }
column_labels <- column_order

# Loop through each chromosome and create heatmap
for (chromosome in chromosomes) {
  pdf(paste0("/Users/colettemoses/Desktop/heatmaps/heatmap_comparative_bowtie2_", chromosome, ".pdf"))
  create_heatmap_bowtie2(bowtie2_heatmap_data, chromosome, paste("bowtie2", chromosome), column_order, column_labels)
  dev.off()
}

# for (chromosome in chromosomes) {
#   pdf(paste0("/Users/colettemoses/Desktop/heatmaps/heatmap_", chromosome, ".pdf"))
#   create_heatmap_bowtie2(bowtie2_heatmap_data, chromosome, paste("bowtie2", chromosome))
#   dev.off()
# }

# # Generate heatmap for all chromosomes combined
# pdf("/Users/colettemoses/Desktop/heatmaps/heatmap_all_chromosomes.pdf")
# create_heatmap_bowtie2(bowtie2_heatmap_data, name = "All Chromosomes")
# dev.off()

# #### For linear colouring
# 
# # Function to create heatmap for a specific chromosome or all chromosomes without logarithmic color scaling
# create_heatmap_bowtie2_linear <- function(data, chromosome = NULL, name) {
#   if (!is.null(chromosome)) {
#     # Filter data for the specified chromosome
#     filtered_data <- data %>%
#       filter(Chromosome_1.5end_partner == chromosome)
#   } else {
#     # Use all data if chromosome is not specified
#     filtered_data <- data
#   }
#   
#   if (nrow(filtered_data) > 0) {
#     # Order filtered_data by Fusion_Transcript_Name
#     filtered_data <- filtered_data[order(filtered_data$Position_on_chromosome_1.5end_partner), ]
#     
#     # Convert filtered_data to matrix
#     heatmap_matrix <- as.matrix(filtered_data[, -c(1:5)])  # Exclude metadata columns
#     
#     # Set color palette for the heatmap without logarithmic scale
#     heatmap_colors_linear <- heatmap_colors(100)
#     
#     # Extract column names from transposed matrix (samples)
#     labCol <- colnames(heatmap_matrix)
#     
#     # Extract row names from transposed matrix (fusion transcripts)
#     labRow <- filtered_data$Fusion_Transcript_Name
#     
#     # Plot the heatmap using heatmap.2 with updated labCol and labRow
#     heatmap.2(heatmap_matrix,
#               col = heatmap_colors_linear,
#               scale = "none",
#               cexCol = 0.1,
#               cexRow = 0.2, 
#               Rowv = FALSE,
#               Colv = TRUE,
#               labCol = rep(" ", ncol(heatmap_matrix)),
#               labRow = labRow,
#               key = FALSE, 
#               trace = "none", 
#               dendrogram = "none", 
#               main = name,
#               cex = 3)
#   } else {
#     cat("No data found for chromosome", chromosome, "\n")
#   }
# }
# 
# # Vector of chromosome names
# chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
#                  "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
#                  "chr20", "chr21", "chr22", "chrX", "chrY")
# 
# for (chromosome in chromosomes) {
#   pdf(paste0("/Users/colettemoses/Desktop/heatmaps/heatmap_", chromosome, ".pdf"))
#   create_heatmap_bowtie2_linear(bowtie2_heatmap_data, chromosome, paste("bowtie2", chromosome))
#   dev.off()
# }

#########################################################

# To generate a heatmap of fusions where there was no expression seen in chimps and high expression in humans:

# Calculating the frequency of transcript expression so that selecting 20 chimp individuals with no expression p < 0.05
root_value <- 0.05^(1/20)
# OR
root_value <- 0.05^(1/33)
root_value <- 0.01^(1/33)
p <- 1 - root_value
print(p) # p = 0.1391083 OR 0.08678119 OR 0.130251
# p = 0.1391083 means that if the frequency of transcript X expression in the population is 13.9108% (or 8.678119%) or higher,
# the probability of selecting 20 individuals with no expression is less than or equal to 0.05.
# Any fusions that are expressed in more than 13.9108% of human individuals and zero chimp individuals are further candidates for human specificity.

# Read the list of fusions exceeding this cutoff into a vector
txt_file_path <- "/Users/colettemoses/Desktop/Good_scripts/Present_Mayo_greaterthan13pc.txt"
txt_file_path_p.05 <- "/Users/colettemoses/Desktop/Good_scripts/Present_Mayo_greaterthan8.678119pc.txt"
txt_file_path_p.01 <- "/Users/colettemoses/Desktop/Good_scripts/Present_Mayo_greaterthan13.0251pc.txt" # chimp and rhesus, p<0.01
transcript_names_p.05 <- readLines(txt_file_path_p.05)
transcript_names_p.01 <- readLines(txt_file_path_p.01)

# # Re-merge bowtie2_660_results_t with collapsed_by_FusionID_filtered_rmsk_2ormore using Fusion_Transcript_Name
# bowtie2_heatmap_data <- bowtie2_660_results_t %>%
#   left_join(collapsed_by_FusionID_filtered_rmsk_2ormore, by = c("Fusion_Transcript_Name" = "FusionID_new"))

# Filter bowtie2_heatmap_data to include only rows where Fusion_Transcript_Name is in the transcript_names vector
bowtie2_heatmap_data_filtered_p.05 <- bowtie2_heatmap_data[bowtie2_heatmap_data$Fusion_Transcript_Name %in% transcript_names_p.05, ]
bowtie2_heatmap_data_filtered_p.01 <- bowtie2_heatmap_data[bowtie2_heatmap_data$Fusion_Transcript_Name %in% transcript_names_p.01, ]

# # Clean up the dataframe to be in proper format for heatmap
# bowtie2_heatmap_data_filtered <- bowtie2_heatmap_data_filtered %>%
#   mutate(Fusion_Transcript_Name = FusionID_reformatted)
# bowtie2_heatmap_data_filtered <- bowtie2_heatmap_data_filtered %>%
#   dplyr::select(-c(289:308, 313:339))
# bowtie2_heatmap_data_filtered <- bowtie2_heatmap_data_filtered[, c(1, 289:292, 2:288)]

# Filtering for all chimp samples = 0
# columns_to_check <- grep("chimp", colnames(bowtie2_heatmap_data_filtered), value = TRUE, ignore.case = TRUE)
# columns_to_check <- c("SRR1758917", "SRR5804451", "SRR5804456", "SRR5804464", 
#                       "SRR5804471", "SRR5804478", "SRR5804486", "SRR6334913", 
#                       "SRR6334974", "SRR6334986", "SRR6335020", "SRR6335047", 
#                       "SRR2040584", "SRR2040585", "ERR3474012", "ERR3474022", 
#                       "SRR17380399", "SRR17380400", "SRR17380401", "SRR17380402")
# bowtie2_heatmap_data_chimpzero <- bowtie2_heatmap_data_filtered[apply(bowtie2_heatmap_data_filtered[columns_to_check], 1, function(row) all(row == 0)), ]

# Filtering for all chimp and rhesus samples = 0
columns_to_check <- grep("chimp|rhesus", colnames(bowtie2_heatmap_data_filtered_p.01), value = TRUE, ignore.case = TRUE)
# columns_to_check <- c("SRR1758917", "SRR5804451", "SRR5804456", "SRR5804464", 
#                       "SRR5804471", "SRR5804478", "SRR5804486", "SRR6334913", 
#                       "SRR6334974", "SRR6334986", "SRR6335020", "SRR6335047", 
#                       "SRR2040584", "SRR2040585", "ERR3474012", "ERR3474022", 
#                       "SRR17380399", "SRR17380400", "SRR17380401", "SRR17380402", 
#                       "SRR6334947", "SRR6334962", "SRR6334989", "SRR6334996", 
#                       "SRR6335010", "SRR2040592", "ERR3474062", "ERR3474072", 
#                       "ERR3474082", "SRR17380395", "SRR17380396", "SRR17380397", 
#                       "SRR17380398")
bowtie2_heatmap_data_chimprhesuszero_p.01 <- bowtie2_heatmap_data_filtered_p.01[apply(bowtie2_heatmap_data_filtered_p.01[columns_to_check], 1, function(row) all(row == 0)), ]

# Vector of chromosome names
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                 "chr20", "chr21", "chr22", "chrX", "chrY")
# Vector of chromosome names
# chromosomes <- c("chr18")

# Loop through each chromosome and create heatmap
for (chromosome in chromosomes) {
  # Create a PDF for each chromosome heatmap
  pdf(paste0("/Users/colettemoses/Desktop/heatmaps/heatmap_", chromosome, ".pdf"))
  
  # Create heatmap for the current chromosome
  create_heatmap_bowtie2(bowtie2_heatmap_data_chimprhesuszero, chromosome, paste("bowtie2", chromosome), column_order, column_labels)
  
  # Close the PDF file
  dev.off()
}

# Filter master list of fusions based on those that are 'human-specific'
humanspecificp.05 <- bowtie2_heatmap_data_chimprhesuszero_p.05$Fusion_Transcript_Name
humanspecificp.01 <- bowtie2_heatmap_data_chimprhesuszero_p.01$Fusion_Transcript_Name
collapsed_by_FusionID_filtered_rmsk_2ormore_HSp.05 <- collapsed_by_FusionID_filtered_rmsk_2ormore[collapsed_by_FusionID_filtered_rmsk_2ormore$FusionID_reformatted %in% humanspecificp.05, ]
collapsed_by_FusionID_filtered_rmsk_2ormore_HSp.01 <- collapsed_by_FusionID_filtered_rmsk_2ormore[collapsed_by_FusionID_filtered_rmsk_2ormore$FusionID_reformatted %in% humanspecificp.01, ]

## Write the data frames to Excel files
HS_file_pathp0.05 <- "/Users/colettemoses/Desktop/Good_scripts/comparative/Mayo_FusionCatcher_compiled_2ormore_v2_HSp0.05.xlsx"
HS_file_pathp0.01 <- "/Users/colettemoses/Desktop/Good_scripts/comparative/Mayo_FusionCatcher_compiled_2ormore_v2_HSp0.01.xlsx"
write_xlsx(list(collapsed_by_FusionID_filtered_rmsk_2ormore_HSp.05), HS_file_pathp0.05, col_names = TRUE)
write_xlsx(list(collapsed_by_FusionID_filtered_rmsk_2ormore_HSp.01), HS_file_pathp0.01, col_names = TRUE)
