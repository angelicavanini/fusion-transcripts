# Clean up the R environment
rm(list=ls())

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(car)
library(FSA)

setwd("/Users/colettemoses/Desktop")

#################### Copy number information ########################################
# Read the file
data <- readLines("/Users/colettemoses/Desktop/Good_scripts/genotyping/ratios_all.txt")

# Initialize variables to store data
cell_lines <- c()
regions <- c()
ratios <- c()

# Loop through each line in the file
for (line in data) {
  # Skip empty lines
  if (nchar(trimws(line)) == 0) {
    next
  }
  
  # Check if the line contains a cell line identifier
  if (grepl("^HPSI", line)) {
    # Extract cell line identifier
    cell_line <- line
  } else {
    # Split the line into region and ratio
    parts <- strsplit(line, ": ")[[1]]
    # Extract region and ratio
    region <- parts[1]
    ratio <- as.numeric(parts[2])
    # Store the data
    cell_lines <- c(cell_lines, cell_line)
    regions <- c(regions, region)
    ratios <- c(ratios, ratio)
  }
}

# Create a data frame from the extracted data
df_genotyperatios <- data.frame(Cell_Line = cell_lines, Region = regions, Ratio = ratios)

# Pivot the data to have regions as columns
df_pivoted <- df_genotyperatios %>%
  pivot_wider(names_from = Region, values_from = Ratio)

# Perform operations on the specified columns
df_normalised <- df_pivoted %>%
  mutate(`17q2131_test_region_3` = `17q2131_test_region_3` * 30000 / 197363,
         `17q2131_test_region_4` = `17q2131_test_region_4` * 30000 / 197499,
         `5q132_test_region_1` = `5q132_test_region_1` * 100000 / 18968)
df_normalised <- df_normalised %>%
  mutate(`17q2131_test_region_3_4_sum` = `17q2131_test_region_3` + `17q2131_test_region_4`)

# Add the new summed region to df_long
df_long <- tidyr::pivot_longer(df_normalised, cols = -Cell_Line, names_to = "Region", values_to = "Ratio")

# Plot the data with jittered x-coordinates and violin plot
plot <- ggplot(df_long, aes(x = Region, y = Ratio)) +
  geom_violin(trim = FALSE, fill = "lightgrey", color = "darkgrey") +  # Add violin plot with specified fill and outline color
  geom_jitter(width = 0.3, height = 0, alpha = 0.7, color = "black") +  # Add jitter to x-coordinates
  labs(title = "Coverage Ratios by Region",
       y = "Coverage Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),  # Remove x-axis title
        panel.grid.major.x = element_blank(),  # Remove vertical gridlines
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  ylim(0, NA)  # Set y-axis to start at zero

# Save the plot
ggsave("coverage_ratios_plot.pdf", plot = plot)
print(plot)

# Loop through each region
for (region in unique(df_long$Region)) {
  # Subset the data for the current region
  df_region <- subset(df_long, Region == region)
  
  # Plot the data with jittered x-coordinates and violin plot
  p <- ggplot(df_region, aes(x = Region, y = Ratio)) +
    geom_violin(trim = FALSE, fill = "lightgrey", color = "darkgrey") +  # Add violin plot with specified fill and outline color
    geom_jitter(width = 0.3, height = 0, alpha = 0.7, color = "black") +  # Add jitter to x-coordinates
    labs(title = paste("Coverage Ratios for", region),
         y = "Coverage Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),  # Remove x-axis title
          panel.grid.major.x = element_blank(),  # Remove vertical gridlines
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(size = 0.5)) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    ylim(0, NA)  # Set y-axis to start at zero
  
  # Save the plot
  ggsave(paste("coverage_ratios_", gsub(":", "_", region), ".pdf", sep = ""), plot = p)
}

#################### 17q21.31 SNP information ########################################

# Read the VCF file into a dataframe
pileup_results <- read.table("/Users/colettemoses/Desktop/Good_scripts/genotyping/H1H2_SNP_all_variants.vcf", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create a new dataframe to store the H1H2 status
H1H2_status <- pileup_results %>%
  # Extract the prefix part of the Cell_Line column
  mutate(Cell_Line = gsub("_.*", "", CELL_LINE)) %>%
  # Group by the Cell_Line
  group_by(Cell_Line) %>%
  # If there are multiple rows for a Cell_Line, combine the Genotype values
  summarize(Genotype = if(n() > 1) {
    # If all Genotype values are the same, return that value, otherwise return "Genotype not matching"
    if(all(GENOTYPE == GENOTYPE[1])) {
      GENOTYPE[1]
    } else {
      "Genotype not matching"
    }
  } else {
    # If there is only one row, return the Genotype value
    GENOTYPE[1]
  })

#################### 17q21.31 H1D/H2D genotyping based on CN and SNPs ########################################

# Define a function to categorize the values
categorize <- function(value) {
  ifelse(value < 1.25, 1, ifelse(value < 1.75, 1.5, 2))
}

# Apply the categorization function to the columns and create new columns
df_categorized <- df_normalised %>%
  mutate(CN_H1D_region = categorize(`17q2131_test_region_1`),
         CN_H2D_region = categorize(`17q2131_test_region_2`))
df_categorized <- df_categorized %>%
  mutate(H1Ddups = CN_H1D_region - 1)
df_categorized <- df_categorized %>%
  mutate(H2Ddups = CN_H2D_region - CN_H1D_region)

# Add the H1H2_genotype column to df_categorized using left_join
df_categorized <- left_join(df_categorized, H1H2_status %>% dplyr::select(Cell_Line, Genotype), by = "Cell_Line")
df_categorized <- df_categorized %>%
  mutate(Full_genotype = case_when(
    Genotype == "0/0" & H1Ddups == 0.0 & H2Ddups == 0.0 ~ "H1/H1",
    Genotype == "0/0" & H1Ddups == 0.5 & H2Ddups == 0.0 ~ "H1/H1D",
    Genotype == "0/0" & H1Ddups == 1.0 & H2Ddups == 0.0 ~ "H1D/H1D",
    Genotype == "0/1" & H1Ddups == 0.0 & H2Ddups == 0.0 ~ "H1/H2",
    Genotype == "0/1" & H1Ddups == 0.5 & H2Ddups == 0.0 ~ "H1D/H2",
    Genotype == "0/1" & H1Ddups == 0.0 & H2Ddups == 0.5 ~ "H1/H2D",
    Genotype == "0/1" & H1Ddups == 0.5 & H2Ddups == 0.5 ~ "H1D/H2D",
    Genotype == "1/1" & H1Ddups == 0.0 & H2Ddups == 0.0 ~ "H2/H2",
    Genotype == "1/1" & H1Ddups == 0.0 & H2Ddups == 0.5 ~ "H2/H2D",
    Genotype == "1/1" & H1Ddups == 0.0 & H2Ddups == 1.0 ~ "H2D/H2D",
    TRUE ~ "Unexpected"
  ))

# Loop through each region
for (region in unique(df_long$Region)) {
  # Subset the data for the current region
  df_region <- subset(df_long, Region == region)
  
  # Merge df_region with df_categorized based on Cell_Line
  df_merged <- merge(df_region, df_categorized, by = "Cell_Line", all.x = TRUE)
  
  # Plot the data with jittered x-coordinates and violin plot, grouped by Full_genotype
  p <- ggplot(df_merged, aes(x = Full_genotype, y = Ratio)) +
    geom_violin(trim = FALSE, fill = "lightgrey", color = "darkgrey") +  # Add violin plot with specified fill and outline color
    geom_jitter(width = 0.3, height = 0, alpha = 0.7, color = "black") +  # Add jitter to x-coordinates
    labs(title = paste("Coverage Ratios for", region),
         y = "Coverage Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),  # Remove x-axis title
          panel.grid.major.x = element_blank(),  # Remove vertical gridlines
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(size = 0.5)) +
    ylim(0, NA)  # Set y-axis to start at zero
  
  # Save the plot
  ggsave(paste("coverage_ratios_grouped", gsub(":", "_", region), ".pdf", sep = ""), plot = p)
}

#################### 17q21.31 H1D (ONLY) genotyping based on CN ########################################

# Define the cell lines to highlight for each category
blue_cell_lines <- c("HPSI1013i-wopl", "HPSI1014i-nosn", "HPSI1113i-hayt", "HPSI0913i-eika", "HPSI0314i-bubh")
red_cell_lines <- c("HPSI0114i-zoxy", "HPSI0414i-oikd", "HPSI0114i-zapk")
orange_cell_line <- "HPSI1113i-uofv"

# Create a list of region names
region_names <- c("17q2131_test_region_1")

# Define a function to categorize the values
categorize_H1D <- function(value) {
  ifelse(value < 1.25, "CN<1.25",
         ifelse(value < 1.75, "CN<1.75", "CN>1.75"))
}

# Apply the categorization function to the column and create a new column
df_categorized <- df_categorized %>%
  mutate(CN_H1D_region = categorize_H1D(`17q2131_test_region_1`))

# Subset the data for the specific region
df_region <- subset(df_long, Region == "17q2131_test_region_1")

# Merge df_region with df_categorized based on Cell_Line
df_merged <- merge(df_region, df_categorized, by = "Cell_Line", all.x = TRUE)

# Plot the data with jittered x-coordinates and violin plot, grouped by CN_H1D_region
p <- ggplot(df_merged, aes(x = CN_H1D_region, y = Ratio)) +
  geom_violin(trim = FALSE, fill = "lightgrey", color = "darkgrey") +  # Add violin plot with specified fill and outline color
  geom_jitter(width = 0.3, height = 0, alpha = 0.7) +  # Add jitter to x-coordinates, removed color aesthetic
  labs(title = "Coverage Ratios for 17q2131_test_region_1",
       y = "Coverage Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),  # Remove x-axis title
        panel.grid.major.x = element_blank(),  # Remove vertical gridlines
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  ylim(0, NA) +  # Set y-axis to start at zero
  # Highlight specific cell lines with slightly transparent points
  geom_point(data = subset(df_merged, Cell_Line %in% red_cell_lines), aes(color = "red"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_merged, Cell_Line %in% red_cell_lines), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "red") +  # Add text labels for red points with adjusted parameters
  geom_point(data = subset(df_merged, Cell_Line %in% blue_cell_lines), aes(color = "blue"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_merged, Cell_Line %in% blue_cell_lines), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "blue") +  # Add text labels for blue points with adjusted parameters
  geom_point(data = subset(df_merged, Cell_Line == orange_cell_line), aes(color = "orange"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_merged, Cell_Line == orange_cell_line), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "orange") +  # Add text labels for orange point with adjusted parameters
  # Add legend
  scale_color_manual(name = "Expression Status",
                     values = c("red" = "red", 
                                "blue" = "blue", 
                                "orange" = "orange"),
                     labels = c("legend text", 
                                "legend text", 
                                "legend text"))

# Save the plot
ggsave("H1D_CN_17q2131_test_region_1_grouped.pdf", plot = p)

#################### 17q21.31 H2D (ONLY) genotyping based on CN ########################################

# Define the cell lines to highlight for each category
blue_cell_lines <- c("HPSI1013i-wopl", "HPSI1014i-nosn", "HPSI1113i-hayt", "HPSI0913i-eika", "HPSI0314i-bubh")
red_cell_lines <- c("HPSI0114i-zoxy", "HPSI0414i-oikd", "HPSI0114i-zapk")
orange_cell_line <- "HPSI1113i-uofv"

# Create a list of region names
region_names <- c("17q2131_test_region_2_1_difference")

# Add an extra column to find H2D region - H1D region
df_categorized <- df_categorized %>%
  mutate(`17q2131_test_region_2_1_difference` = `17q2131_test_region_2` - `17q2131_test_region_1` + 1)

# Define a function to categorize the values
categorize_H2D <- function(value) {
  ifelse(value < 1.25, "CN<1.25",
         ifelse(value < 1.75, "CN<1.75", "CN>1.75"))
}

# Apply the categorization function to the column and create a new column
df_categorized <- df_categorized %>%
  mutate(CN_H2D_region = categorize_H2D(`17q2131_test_region_2_1_difference`))

# # Subset the data for the specific region
# df_region <- subset(df_long, Region == "17q2131_test_region_2_1_difference")
# 
# # Merge df_region with df_categorized based on Cell_Line
# df_merged <- merge(df_region, df_categorized, by = "Cell_Line", all.x = TRUE)

# Plot the data with jittered x-coordinates and violin plot, grouped by CN_H2D_region
p <- ggplot(df_categorized, aes(x = CN_H2D_region, y = `17q2131_test_region_2_1_difference`)) +
  geom_violin(trim = FALSE, fill = "lightgrey", color = "darkgrey") +  # Add violin plot with specified fill and outline color
  geom_jitter(width = 0.3, height = 0, alpha = 0.7) +  # Add jitter to x-coordinates, removed color aesthetic
  labs(title = "Coverage Ratios for 17q2131_test_region_2-1",
       y = "Coverage Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),  # Remove x-axis title
        panel.grid.major.x = element_blank(),  # Remove vertical gridlines
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  ylim(0, NA) +  # Set y-axis to start at zero
  # Highlight specific cell lines with slightly transparent points
  geom_point(data = subset(df_merged, Cell_Line %in% red_cell_lines), aes(color = "red"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_merged, Cell_Line %in% red_cell_lines), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "red") +  # Add text labels for red points with adjusted parameters
  geom_point(data = subset(df_merged, Cell_Line %in% blue_cell_lines), aes(color = "blue"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_merged, Cell_Line %in% blue_cell_lines), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "blue") +  # Add text labels for blue points with adjusted parameters
  geom_point(data = subset(df_merged, Cell_Line == orange_cell_line), aes(color = "orange"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_merged, Cell_Line == orange_cell_line), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "orange") +  # Add text labels for orange point with adjusted parameters
  # Add legend
  scale_color_manual(name = "Expression Status",
                     values = c("red" = "red", 
                                "blue" = "blue", 
                                "orange" = "orange"),
                     labels = c("legend text", 
                                "legend text", 
                                "legend text"))

# Save the plot
ggsave("H2D_CN_17q2131_test_region_2-1_grouped.pdf", plot = p)

#################### OCLNP1 genotyping based on CN ########################################

###### Plot showing expression in our cell lines from PCR ###### 

# Define the cell lines to highlight for each category
blue_cell_lines <- c("HPSI1013i-wopl", "HPSI1014i-nosn", "HPSI0114i-zapk")
red_cell_lines <- c("HPSI1113i-uofv", "HPSI0114i-zoxy", "HPSI1113i-hayt", "HPSI0913i-eika", "HPSI0314i-bubh")
orange_cell_line <- "HPSI0414i-oikd"

# Subset the data for the region of interest
df_region <- subset(df_long, Region == "5q132_test_region_1")

# Plot the data with jittered x-coordinates and violin plot
p <- ggplot(df_region, aes(x = Region, y = Ratio)) +
  geom_violin(trim = FALSE, fill = "lightgrey", color = "lightgrey") +  # Add violin plot with specified fill and outline color
  geom_jitter(width = 0.3, height = 0, alpha = 0.7, color = "darkgrey") +  # Add jitter to x-coordinates
  labs(title = "Coverage Ratios for 5q132_test_region_1",
       y = "Coverage Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),  # Remove x-axis title
        panel.grid.major.x = element_blank(),  # Remove vertical gridlines
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  ylim(0, NA) +  # Set y-axis to start at zero
  # Highlight specific cell lines with slightly transparent points
  geom_point(data = subset(df_region, Cell_Line %in% red_cell_lines), aes(color = "red"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_region, Cell_Line %in% red_cell_lines), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "red") +  # Add text labels for red points with adjusted parameters
  geom_point(data = subset(df_region, Cell_Line %in% blue_cell_lines), aes(color = "blue"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_region, Cell_Line %in% blue_cell_lines), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "blue") +  # Add text labels for blue points with adjusted parameters
  geom_point(data = subset(df_region, Cell_Line == orange_cell_line), aes(color = "orange"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_region, Cell_Line == orange_cell_line), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "orange") +  # Add text labels for orange point with adjusted parameters
  # Add legend
  scale_color_manual(name = "Expression Status",
                     values = c("red" = "red", 
                                "blue" = "blue", 
                                "orange" = "orange"),
                     labels = c("No NAIP-OCLN expression, n=3", 
                                "NAIP-OCLN expression in PCR only, n=1", 
                                "NAIP-OCLN expression in FusionCatcher and PCR, n=5"))

# Save the plot
ggsave("OCLNP1_CN_5q132_test_region_1.pdf", plot = p)

#################### Same plot but grouped according to OCLNP1 copy number

# Define a function to categorize the values
categorize_OCLNP1 <- function(value) {
  ifelse(value < 0.5, "CN<0.5",
         ifelse(value < 0.9, "CN<0.9",
                ifelse(value < 1.3, "CN<1.3", "CN>1.3")))
}

# Apply the categorization function to the column and create a new column
df_categorized <- df_categorized %>%
  mutate(CN_OCLNP1_region = categorize_OCLNP1(`5q132_test_region_1`))

# Subset the data for the specific region
df_region <- subset(df_long, Region == "5q132_test_region_1")

# Merge df_region with df_categorized based on Cell_Line
df_merged <- merge(df_region, df_categorized, by = "Cell_Line", all.x = TRUE)

# Plot the data with jittered x-coordinates and violin plot, grouped by CN_OCLNP1_region
p <- ggplot(df_merged, aes(x = CN_OCLNP1_region, y = Ratio)) +
  geom_violin(trim = FALSE, fill = "lightgrey", color = "darkgrey") +  # Add violin plot with specified fill and outline color
  geom_jitter(width = 0.3, height = 0, alpha = 0.7) +  # Add jitter to x-coordinates, removed color aesthetic
  labs(title = "Coverage Ratios for 5q132_test_region_1",
       y = "Coverage Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),  # Remove x-axis title
        panel.grid.major.x = element_blank(),  # Remove vertical gridlines
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  ylim(0, NA) +  # Set y-axis to start at zero
  # Highlight specific cell lines with slightly transparent points
  geom_point(data = subset(df_merged, Cell_Line %in% red_cell_lines), aes(color = "red"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_merged, Cell_Line %in% red_cell_lines), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "red") +  # Add text labels for red points with adjusted parameters
  geom_point(data = subset(df_merged, Cell_Line %in% blue_cell_lines), aes(color = "blue"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_merged, Cell_Line %in% blue_cell_lines), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "blue") +  # Add text labels for blue points with adjusted parameters
  geom_point(data = subset(df_merged, Cell_Line == orange_cell_line), aes(color = "orange"), size = 3, alpha = 0.5) +
  geom_text_repel(data = subset(df_merged, Cell_Line == orange_cell_line), aes(label = Cell_Line), 
                  hjust = -0.2, vjust = 0,
                  box.padding = 0.5, point.padding = 0.5, force = 2, color = "orange") +  # Add text labels for orange point with adjusted parameters
  # Add legend
  scale_color_manual(name = "Expression Status",
                     values = c("red" = "red", 
                                "blue" = "blue", 
                                "orange" = "orange"),
                     labels = c("No NAIP-OCLN expression, n=3", 
                                "NAIP-OCLN expression in PCR only, n=1", 
                                "NAIP-OCLN expression in FusionCatcher and PCR, n=5"))

# Save the plot
ggsave("OCLNP1_CN_5q132_test_region_1_grouped.pdf", plot = p)

#################### 17q21.31 NSF genotyping based on CN ########################################

# Define the cell lines to highlight for each category
blue_cell_lines <- c("HPSI1013i-wopl", "HPSI1014i-nosn", "HPSI1113i-hayt", "HPSI0913i-eika", "HPSI0314i-bubh")
red_cell_lines <- c("HPSI0114i-zoxy", "HPSI0414i-oikd", "HPSI0114i-zapk")
orange_cell_line <- "HPSI1113i-uofv"

# Create a list of region names
region_names <- c("17q2131_test_region_3", "17q2131_test_region_4", "17q2131_test_region_3_4_sum")

for (region_name in region_names) {
  # Subset the data for the region of interest
  df_region <- subset(df_long, Region == region_name)
  
  # Plot the data with jittered x-coordinates and violin plot
  p <- ggplot(df_region, aes(x = Region, y = Ratio)) +
    geom_violin(trim = FALSE, fill = "lightgrey", color = "lightgrey") +  # Add violin plot with specified fill and outline color
    geom_jitter(width = 0.3, height = 0, alpha = 0.7, color = "darkgrey") +  # Add jitter to x-coordinates
    labs(title = paste("Coverage Ratios for", region_name),
         y = "Coverage Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),  # Remove x-axis title
          panel.grid.major.x = element_blank(),  # Remove vertical gridlines
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(size = 0.5)) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    ylim(0, NA) +  # Set y-axis to start at zero
    # Highlight specific cell lines with slightly transparent points
    geom_point(data = subset(df_region, Cell_Line %in% red_cell_lines), aes(color = "red"), size = 3, alpha = 0.5) +
    geom_text_repel(data = subset(df_region, Cell_Line %in% red_cell_lines), aes(label = Cell_Line), 
                    hjust = -0.2, vjust = 0,
                    box.padding = 0.5, point.padding = 0.5, force = 2, color = "red") +  # Add text labels for red points with adjusted parameters
    geom_point(data = subset(df_region, Cell_Line %in% blue_cell_lines), aes(color = "blue"), size = 3, alpha = 0.5) +
    geom_text_repel(data = subset(df_region, Cell_Line %in% blue_cell_lines), aes(label = Cell_Line), 
                    hjust = -0.2, vjust = 0,
                    box.padding = 0.5, point.padding = 0.5, force = 2, color = "blue") +  # Add text labels for blue points with adjusted parameters
    geom_point(data = subset(df_region, Cell_Line == orange_cell_line), aes(color = "orange"), size = 3, alpha = 0.5) +
    geom_text_repel(data = subset(df_region, Cell_Line == orange_cell_line), aes(label = Cell_Line), 
                    hjust = -0.2, vjust = 0,
                    box.padding = 0.5, point.padding = 0.5, force = 2, color = "orange") +  # Add text labels for orange point with adjusted parameters
    # Add legend
    scale_color_manual(name = "Expression Status",
                       values = c("red" = "red", 
                                  "blue" = "blue", 
                                  "orange" = "orange"),
                       labels = c("No NSF-LRRC37A3 expression, n=5", 
                                  "Weak NSF-LRRC37A3 expression, n=1", 
                                  "Strong NSF-LRRC37A3 expression, n=3"))
  
  # Save the plot
  ggsave(paste("NSF_CN_", region_name, ".pdf", sep = ""), plot = p)
}

#################### Same plot but grouped according to NSF copy number

# Define a function to categorize the values for region 17q2131_test_region_3
# categorize_NSF_region3 <- function(value) {
#   ifelse(value < 0.6, "Low", "High")
# }

categorize_NSF_region3 <- function(value) {
  ifelse(value < 0.6, "CN<0.60",
         ifelse(value < 0.85, "CN<0.85",
                ifelse(value < 1.1, "CN<1.10",
                       ifelse(value < 1.35, "CN<1.35", "CN>1.35")
                )
         )
  )
}

# Define a function to categorize the values for region 17q2131_test_region_4
categorize_NSF_region4 <- function(value) {
  ifelse(value < 0.875, "Low", "High")
}

# Define a function to categorize the values for the summed region
categorize_NSF_region_sum <- function(value) {
  ifelse(value < 1.75, "Low", "High")
}

# Apply the categorization function to the respective columns and create new columns
df_categorized <- df_categorized %>%
  mutate(
    CN_NSF_17q2131_test_region_3 = categorize_NSF_region3(`17q2131_test_region_3`),
    CN_NSF_17q2131_test_region_4 = categorize_NSF_region4(`17q2131_test_region_4`),
    CN_NSF_17q2131_test_region_3_4_sum = categorize_NSF_region_sum(`17q2131_test_region_3_4_sum`)
  )

# Create a list of region names
region_names <- c("17q2131_test_region_3", "17q2131_test_region_4", "17q2131_test_region_3_4_sum")

# Loop through each region to generate plots
for (region_name in region_names) {
  # Subset the data for the specific region
  df_region <- subset(df_long, Region == region_name)
  
  # Merge df_region with df_categorized based on Cell_Line
  df_merged <- merge(df_region, df_categorized, by = "Cell_Line", all.x = TRUE)
  
  # Plot the data with jittered x-coordinates and violin plot, grouped by CN_OCLNP1_region
  p <- ggplot(df_merged, aes(x = get(paste0("CN_NSF_", region_name)), y = Ratio)) +
    geom_violin(trim = FALSE, fill = "lightgrey", color = "darkgrey") +  # Add violin plot with specified fill and outline color
    geom_jitter(width = 0.3, height = 0, alpha = 0.7) +  # Add jitter to x-coordinates, removed color aesthetic
    labs(title = paste("Coverage Ratios for", region_name),
         y = "Coverage Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),  # Remove x-axis title
          panel.grid.major.x = element_blank(),  # Remove vertical gridlines
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(size = 0.5)) +
    ylim(0, NA) +  # Set y-axis to start at zero
    # Highlight specific cell lines with slightly transparent points
    geom_point(data = subset(df_merged, Cell_Line %in% red_cell_lines), aes(color = "red"), size = 3, alpha = 0.5) +
    geom_text_repel(data = subset(df_merged, Cell_Line %in% red_cell_lines), aes(label = Cell_Line), 
                    hjust = -0.2, vjust = 0,
                    box.padding = 0.5, point.padding = 0.5, force = 2, color = "red") +  # Add text labels for red points with adjusted parameters
    geom_point(data = subset(df_merged, Cell_Line %in% blue_cell_lines), aes(color = "blue"), size = 3, alpha = 0.5) +
    geom_text_repel(data = subset(df_merged, Cell_Line %in% blue_cell_lines), aes(label = Cell_Line), 
                    hjust = -0.2, vjust = 0,
                    box.padding = 0.5, point.padding = 0.5, force = 2, color = "blue") +  # Add text labels for blue points with adjusted parameters
    geom_point(data = subset(df_merged, Cell_Line == orange_cell_line), aes(color = "orange"), size = 3, alpha = 0.5) +
    geom_text_repel(data = subset(df_merged, Cell_Line == orange_cell_line), aes(label = Cell_Line), 
                    hjust = -0.2, vjust = 0,
                    box.padding = 0.5, point.padding = 0.5, force = 2, color = "orange") +  # Add text labels for orange point with adjusted parameters
    # Add legend
    scale_color_manual(name = "Expression Status",
                       values = c("red" = "red", 
                                  "blue" = "blue", 
                                  "orange" = "orange"),
                       labels = c("No NSF-LRRC37A3 expression, n=5", 
                                  "Weak NSF-LRRC37A3 expression, n=1", 
                                  "Strong NSF-LRRC37A3 expression, n=3"))
  
  # Save the plot
  ggsave(paste("NSF_CN_", region_name, "_grouped.pdf", sep = ""), plot = p)
}

#################### Results of bowtie2 search ########################################

bowtie2_counts <- readLines("/Users/colettemoses/Desktop/Good_scripts/bowtie2/Counts_pivot_2024-07-08_16-38-45_hipsci_fulllength_G,32,21.txt")

# Split each line by "\t" delimiter and convert to a dataframe
data_list <- lapply(bowtie2_counts, function(x) unlist(strsplit(x, "\t")))
bowtie2_counts_df <- as.data.frame(do.call(rbind, data_list), stringsAsFactors = FALSE)

bowtie2_counts_df <- as.data.frame(t(bowtie2_counts_df), stringsAsFactors = FALSE)

# Set the first row as column names
colnames(bowtie2_counts_df) <- bowtie2_counts_df[1, ]
bowtie2_counts_df <- bowtie2_counts_df[-1, ]
colnames(bowtie2_counts_df)[1] <- "Prefix"

# Convert all columns except the first one to numeric
bowtie2_counts_df[, -1] <- lapply(bowtie2_counts_df[, -1], function(x) as.numeric(as.character(x)))

# Add a new column for the prefix
bowtie2_counts_df <- bowtie2_counts_df %>%
  mutate(Prefix = sub("_.*", "", Prefix))

# Group rows by prefix, calculate mean for each group, and create a new dataframe
bowtie2_counts_df <- bowtie2_counts_df %>%
  group_by(Prefix) %>%
  summarise_all(mean)

# Calculate the mean of each column
column_means <- colMeans(bowtie2_counts_df[, -1], na.rm = TRUE)

# Calculate the column means needed for normalization
kansl1_exon3_mean <- column_means["KANSL1_KANSL1_17:46094560"]
kansl1_exon2_mean <- column_means["KANSL1_KANSL1_17:46170855"]
kansl1_novelexon_mean <- column_means["KANSL1_KANSL1_17:46152904"]
ocln_mean <- column_means["OCLN_OCLN_5:69534694"]
nsf_mean <- column_means["NSF_NSF_17:46704854"]
# nsfp1_mean <- column_means["NSFP1_17:46704854"]
naip_mean <- column_means["NAIP_NAIP_5:70974129"] + column_means["NAIP_NAIP_5:70979869"] + column_means["NAIP_NAIP_5:70983775"] + column_means["NAIP_NAIP_5:70998724"]

# Create the bowtie2_counts_normalised dataframe
bowtie2_counts_normalised <- data.frame(
  Prefix = bowtie2_counts_df$Prefix,
  KANSL1_KANSL1_exon3 = bowtie2_counts_df$"KANSL1_KANSL1_17:46094560" / kansl1_exon3_mean,
  KANSL1_ARL17 = (bowtie2_counts_df$"KANSL1_ARL17A_17:46094560_17:46570869" + bowtie2_counts_df$"KANSL1_ARL17A_17:46094560_17:46517233") / kansl1_exon3_mean,
  KANSL1_KANSL1_exon2 = bowtie2_counts_df$"KANSL1_KANSL1_17:46170855" / kansl1_exon2_mean,
  KANSL1_KANSL1_novelexon = bowtie2_counts_df$"KANSL1_KANSL1_17:46152904" / kansl1_novelexon_mean,
  KANSL1_LRRC37A3_out_of_frame = (bowtie2_counts_df$"KANSL1_LRRC37A3_17:46152904_17:64869166" + bowtie2_counts_df$"KANSL1_LRRC37A3_17:46152904_17:64892600" + bowtie2_counts_df$"KANSL1_LRRC37A3_17:46152904_17:64868536") / kansl1_exon2_mean,
  KANSL1_LRRC37A3_in_frame = (bowtie2_counts_df$"KANSL1_LRRC37A3_17:46170855_17:64868536" + bowtie2_counts_df$"KANSL1_LRRC37A3_17:46170855_17:64892600" + bowtie2_counts_df$"KANSL1_LRRC37A3_ex5") / kansl1_exon2_mean,
#  LRRC37A3_LRRC37A3 = bowtie2_counts_df$"LRRC37A3_LRRC37A3_17:64868536" / kansl1_exon2_mean,
  NSF_NSF = bowtie2_counts_df$"NSF_NSF_17:46704854" / nsf_mean,
#  NSFP1 = bowtie2_counts_df$"NSFP1_17:46704854" / nsfp1_mean,
  NSF_LRRC37A3 = (bowtie2_counts_df$"NSF_LRRC37A3_17:46704854_17:64869166" + bowtie2_counts_df$"NSF_LRRC37A3_17:46704854_17:64892600" + bowtie2_counts_df$"NSF_LRRC37A3_17:46704854_17:64860973" + bowtie2_counts_df$"NSF_LRRC37A3_17:46704854_17:64897654") / nsf_mean,
#  NSF_LRRC37A3_p1 = (bowtie2_counts_df$"NSF_LRRC37A3_17:46704854_17:64869166" + bowtie2_counts_df$"NSF_LRRC37A3_17:46704854_17:64892600" + bowtie2_counts_df$"NSF_LRRC37A3_17:46704854_17:64860973" + bowtie2_counts_df$"NSF_LRRC37A3_17:46704854_17:64897654") / nsfp1_mean,
  OCLN_OCLN = bowtie2_counts_df$"OCLN_OCLN_5:69534694" / ocln_mean,
  NAIP_OCLN = (bowtie2_counts_df$"NAIP_OCLN_5:70974129_5:69534694" + bowtie2_counts_df$"NAIP_OCLN_5:70979869_5:69534694" + bowtie2_counts_df$"NAIP_OCLN_5:70983775_5:69534694" + bowtie2_counts_df$"NAIP_OCLN_5:70998724_5:69534694") / ocln_mean,
#  NAIP_OCLN_naip = (bowtie2_counts_df$"NAIP_OCLN_5:70974129_5:69534694" + bowtie2_counts_df$"NAIP_OCLN_5:70979869_5:69534694" + bowtie2_counts_df$"NAIP_OCLN_5:70983775_5:69534694" + bowtie2_counts_df$"NAIP_OCLN_5:70998724_5:69534694") / naip_mean,
  NAIP_NAIP = (bowtie2_counts_df$"NAIP_NAIP_5:70974129" + bowtie2_counts_df$"NAIP_NAIP_5:70979869" + bowtie2_counts_df$"NAIP_NAIP_5:70983775" + bowtie2_counts_df$"NAIP_NAIP_5:70998724") / naip_mean
)

# For normalisation, use...
# NAIPP4-OCLNP1: OCLN_OCLN, OCLN_OCLN_5:69534694
# NSFP1-LRRC37A2: NSF_NSF, NSF_NSF_17:46704854
# KANSL1-LRRC37A3 in-frame: KANSL1_KANSL1_exon2, KANSL1_KANSL1_17:46170855
# KANSL1-LRRC37A3 out-of-frame: KANSL1_KANSL1_novelexon, KANSL1_KANSL1_17:46152904 (KANSL1 transcript ENST00000639099.1, ENSE00003807679-ENSE00002635959... ENSE00002635959 is the one fusing to LRRC, but it is the last exon in the transcript)
# KANSL1-ARL17: KANSL1_KANSL1_exon3, KANSL1_KANSL1_17:46094560

##### Plotting #####

# H1D

# Select the columns you want to add from df_categorized
columns_to_add <- c("H1Ddups")

# Left join df_categorized to bowtie2_counts_df based on Prefix and Cell_Line, and select specific columns to add
bowtie2_counts_genotypes <- left_join(bowtie2_counts_normalised, dplyr::select(df_categorized, Cell_Line, all_of(columns_to_add)), by = c("Prefix" = "Cell_Line"))

# Convert H1Ddups to factor
bowtie2_counts_genotypes$H1Ddups <- factor(bowtie2_counts_genotypes$H1Ddups, levels = c(0, 0.5, 1))

# Select the columns you want to use
selected_columns <- c("Prefix", "KANSL1_KANSL1_exon3", "KANSL1_ARL17", "H1Ddups")

# Subset the dataframe
bowtie2_counts_genotypes_subset <- bowtie2_counts_genotypes[selected_columns]

# Rearrange the filtered dataframe to long format
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_subset %>%
  pivot_longer(cols = -c(Prefix, H1Ddups), 
               names_to = "Category", 
               values_to = "Value")

# Filter to include only "KANSL1_KANSL1_exon3" and "KANSL1_ARL17" categories
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_long %>%
  filter(Category %in% c("KANSL1_KANSL1_exon3", "KANSL1_ARL17"))

# Calculate mean and SEM for each group
mean_se_values <- bowtie2_counts_genotypes_long %>%
  group_by(H1Ddups, Category) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
            sem = sd(Value, na.rm = TRUE) / sqrt(n()))

# Plotting with bar to mark mean and error bars for SEM
plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = H1Ddups, y = Value)) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
  geom_bar(data = mean_se_values, aes(x = H1Ddups, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
  geom_errorbar(data = mean_se_values, aes(x = H1Ddups, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
  facet_grid(Category ~ ., scales = "free_y") +
  labs(title = "Normalised Counts by Category",
       y = "Normalised Counts",
       x = "H1Ddups") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  guides(color = FALSE)  # Remove the legend for the color aesthetic

# Print the plot
print(plot)
ggsave("KANSL1-ARL17_fusion_bowtie2_G,32,21.pdf", plot = plot)

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(KANSL1_ARL17 ~ H1Ddups, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result

kruskal_test_result_control <- kruskal.test(KANSL1_KANSL1_exon3 ~ H1Ddups, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result_control

# Perform Dunn's test with Bonferroni correction
dunn_test_result <- dunnTest(KANSL1_ARL17 ~ H1Ddups, 
                             data = bowtie2_counts_genotypes, 
                             method = "bonferroni")
dunn_test_result

# # Keeping y axis consistent:
# 
# # Calculate the maximum value across both categories
# max_value <- max(max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "KANSL1_KANSL1_exon3"]),
#                  max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "KANSL1_ARL17"]))
# 
# # Plotting with bar to mark mean and error bars for SEM, setting y-axis limits
# plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = H1Ddups, y = Value)) +
#   geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
#   geom_bar(data = mean_se_values, aes(x = H1Ddups, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
#   geom_errorbar(data = mean_se_values, aes(x = H1Ddups, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
#   facet_grid(Category ~ ., scales = "free_y") +
#   labs(title = "Normalised Counts by Category",
#        y = "Normalised Counts",
#        x = "H1Ddups") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
#         axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_line(size = 0.5)) +
#   guides(color = FALSE) +  # Remove the legend for the color aesthetic
#   ylim(0, max_value)  # Set y-axis limits
# 
# # Print the plot
# print(plot)
# ggsave("KANSL1-ARL17_fusion_bowtie2_G,32,21_ymax.pdf", plot = plot)

# # Plotting with line to mark mean
# plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = H1Ddups, y = Value)) +
#   geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
#   geom_point(data = mean_values, aes(x = H1Ddups, y = mean_value, color = "red"), shape = "-", size = 10) +  # Add mean values as points
#   facet_grid(Category ~ ., scales = "free_y") +
#   labs(title = "Normalised Counts by Category",
#        y = "Normalised Counts",
#        x = "H1Ddups") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
#         axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_line(size = 0.5)) +
# #  ylim(0, max_value) +
#   guides(color = FALSE)  # Remove the legend for the color aesthetic
# 
# # Print the plot
# print(plot)

# H2D LRRC out-of-frame

# Select the columns you want to add from df_categorized
columns_to_add <- c("H2Ddups")

# Left join df_categorized to bowtie2_counts_df based on Prefix and Cell_Line, and select specific columns to add
bowtie2_counts_genotypes <- left_join(bowtie2_counts_normalised, dplyr::select(df_categorized, Cell_Line, all_of(columns_to_add)), by = c("Prefix" = "Cell_Line"))

# Convert H2Ddups to factor
bowtie2_counts_genotypes$H2Ddups <- factor(bowtie2_counts_genotypes$H2Ddups, levels = c(0, 0.5, 1))

# Select the columns you want to use
selected_columns <- c("Prefix", "KANSL1_KANSL1_exon2", "KANSL1_LRRC37A3_out_of_frame", "H2Ddups")

# Subset the dataframe
bowtie2_counts_genotypes_subset <- bowtie2_counts_genotypes[selected_columns]

# Rearrange the filtered dataframe to long format
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_subset %>%
  pivot_longer(cols = -c(Prefix, H2Ddups), 
               names_to = "Category", 
               values_to = "Value")

# Filter to include only "KANSL1_KANSL1_exon2" and "KANSL1_LRRC37A3_out_of_frame" categories
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_long %>%
  filter(Category %in% c("KANSL1_KANSL1_exon2", "KANSL1_LRRC37A3_out_of_frame"))

# Calculate mean and SEM for each group
mean_se_values <- bowtie2_counts_genotypes_long %>%
  group_by(H2Ddups, Category) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
            sem = sd(Value, na.rm = TRUE) / sqrt(n()))

# Plotting with bar to mark mean and error bars for SEM
plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = H2Ddups, y = Value)) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
  geom_bar(data = mean_se_values, aes(x = H2Ddups, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
  geom_errorbar(data = mean_se_values, aes(x = H2Ddups, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
  facet_grid(Category ~ ., scales = "free_y") +
  labs(title = "Normalised Counts by Category",
       y = "Normalised Counts",
       x = "H2Ddups") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  guides(color = FALSE)  # Remove the legend for the color aesthetic

# Print the plot
print(plot)
ggsave("KANSL1-LRRC_outofframe_fusion_bowtie2_G,32,21.pdf", plot = plot)

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(KANSL1_LRRC37A3_out_of_frame ~ H2Ddups, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result

kruskal_test_result_control <- kruskal.test(KANSL1_KANSL1_exon2 ~ H2Ddups, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result_control

# Perform Dunn's test with Bonferroni correction
dunn_test_result <- dunnTest(KANSL1_LRRC37A3_out_of_frame ~ H2Ddups, 
                             data = bowtie2_counts_genotypes, 
                             method = "bonferroni")
dunn_test_result

# # Y max
# 
# # Calculate the maximum value across both categories
# max_value <- max(max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "KANSL1_KANSL1_exon2"]),
#                  max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "KANSL1_LRRC37A3_out_of_frame"]))
# 
# # Plotting with bar to mark mean and error bars for SEM, setting y-axis limits
# plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = H2Ddups, y = Value)) +
#   geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
#   geom_bar(data = mean_se_values, aes(x = H2Ddups, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
#   geom_errorbar(data = mean_se_values, aes(x = H2Ddups, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
#   facet_grid(Category ~ ., scales = "free_y") +
#   labs(title = "Normalised Counts by Category",
#        y = "Normalised Counts",
#        x = "H2Ddups") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
#         axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_line(size = 0.5)) +
#   guides(color = FALSE) +  # Remove the legend for the color aesthetic
#   ylim(0, max_value)  # Set y-axis limits
# 
# # Print the plot
# print(plot)
# ggsave("KANSL1-LRRC_outofframe_fusion_bowtie2_G,32,21_ymax.pdf", plot = plot)

# H2D LRRC in-frame

# Select the columns you want to add from df_categorized
columns_to_add <- c("H2Ddups")

# Left join df_categorized to bowtie2_counts_df based on Prefix and Cell_Line, and select specific columns to add
bowtie2_counts_genotypes <- left_join(bowtie2_counts_normalised, dplyr::select(df_categorized, Cell_Line, all_of(columns_to_add)), by = c("Prefix" = "Cell_Line"))

# Convert H2Ddups to factor
bowtie2_counts_genotypes$H2Ddups <- factor(bowtie2_counts_genotypes$H2Ddups, levels = c(0, 0.5, 1))

# Select the columns you want to use
selected_columns <- c("Prefix", "KANSL1_KANSL1_exon2", "KANSL1_LRRC37A3_in_frame", "H2Ddups")

# Subset the dataframe
bowtie2_counts_genotypes_subset <- bowtie2_counts_genotypes[selected_columns]

# Rearrange the filtered dataframe to long format
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_subset %>%
  pivot_longer(cols = -c(Prefix, H2Ddups), 
               names_to = "Category", 
               values_to = "Value")

# Filter to include only "KANSL1_KANSL1_exon2" and "KANSL1_LRRC37A3_in_frame" categories
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_long %>%
  filter(Category %in% c("KANSL1_KANSL1_exon2", "KANSL1_LRRC37A3_in_frame"))

# Calculate mean and SEM for each group
mean_se_values <- bowtie2_counts_genotypes_long %>%
  group_by(H2Ddups, Category) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
            sem = sd(Value, na.rm = TRUE) / sqrt(n()))

# Plotting with bar to mark mean and error bars for SEM
plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = H2Ddups, y = Value)) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
  geom_bar(data = mean_se_values, aes(x = H2Ddups, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
  geom_errorbar(data = mean_se_values, aes(x = H2Ddups, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
  facet_grid(Category ~ ., scales = "free_y") +
  labs(title = "Normalised Counts by Category",
       y = "Normalised Counts",
       x = "H2Ddups") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  guides(color = FALSE)  # Remove the legend for the color aesthetic

# Print the plot
print(plot)
ggsave("KANSL1-LRRC_inframe_fusion_bowtie2_G,32,21.pdf", plot = plot)

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(KANSL1_LRRC37A3_in_frame ~ H2Ddups, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result

# Perform Dunn's test with Bonferroni correction
dunn_test_result <- dunnTest(KANSL1_LRRC37A3_in_frame ~ H2Ddups, 
                             data = bowtie2_counts_genotypes, 
                             method = "bonferroni")
dunn_test_result

# control result is above

# # Y max
# 
# # Calculate the maximum value across both categories
# max_value <- max(max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "KANSL1_KANSL1_exon2"]),
#                  max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "KANSL1_LRRC37A3_in_frame"]))
# 
# # Plotting with bar to mark mean and error bars for SEM, setting y-axis limits
# plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = H2Ddups, y = Value)) +
#   geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
#   geom_bar(data = mean_se_values, aes(x = H2Ddups, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
#   geom_errorbar(data = mean_se_values, aes(x = H2Ddups, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
#   facet_grid(Category ~ ., scales = "free_y") +
#   labs(title = "Normalised Counts by Category",
#        y = "Normalised Counts",
#        x = "H2Ddups") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
#         axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_line(size = 0.5)) +
#   guides(color = FALSE) +  # Remove the legend for the color aesthetic
#   ylim(0, max_value)  # Set y-axis limits
# 
# # Print the plot
# print(plot)
# ggsave("KANSL1-LRRC_inframe_fusion_bowtie2_G,32,21_Ymax.pdf", plot = plot)

# NSF-LRRC37A3

# Select the columns you want to add from df_categorized
columns_to_add <- c("CN_NSF_17q2131_test_region_3")

# Left join df_categorized to bowtie2_counts_df based on Prefix and Cell_Line, and select specific columns to add
bowtie2_counts_genotypes <- left_join(bowtie2_counts_normalised, dplyr::select(df_categorized, Cell_Line, all_of(columns_to_add)), by = c("Prefix" = "Cell_Line"))

# Select the columns you want to use
selected_columns <- c("Prefix", "NSF_NSF", "NSF_LRRC37A3", "CN_NSF_17q2131_test_region_3")

# Subset the dataframe
bowtie2_counts_genotypes_subset <- bowtie2_counts_genotypes[selected_columns]

# Rearrange the filtered dataframe to long format
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_subset %>%
  pivot_longer(cols = -c(Prefix, CN_NSF_17q2131_test_region_3), 
               names_to = "Category", 
               values_to = "Value")

# Filter to include only "NSF_NSF" and "NSF_LRRC37A3" categories
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_long %>%
  filter(Category %in% c("NSF_NSF", "NSF_LRRC37A3"))

# Calculate mean and SEM for each group
mean_se_values <- bowtie2_counts_genotypes_long %>%
  group_by(CN_NSF_17q2131_test_region_3, Category) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
            sem = sd(Value, na.rm = TRUE) / sqrt(n()))

# Plotting with bar to mark mean and error bars for SEM
plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = CN_NSF_17q2131_test_region_3, y = Value)) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
  geom_bar(data = mean_se_values, aes(x = CN_NSF_17q2131_test_region_3, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
  geom_errorbar(data = mean_se_values, aes(x = CN_NSF_17q2131_test_region_3, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
  facet_grid(Category ~ ., scales = "free_y") +
  labs(title = "Normalised Counts by Category",
       y = "Normalised Counts",
       x = "CN_NSF_17q2131_test_region_3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  guides(color = FALSE)  # Remove the legend for the color aesthetic

# Print the plot
print(plot)
ggsave("NSF-LRRC_fusion_nsfnorm_bowtie2_G,32,21.pdf", plot = plot)

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(NSF_LRRC37A3 ~ CN_NSF_17q2131_test_region_3, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result

kruskal_test_result_control <- kruskal.test(NSF_NSF ~ CN_NSF_17q2131_test_region_3, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result_control

# Perform Dunn's test with Bonferroni correction
dunn_test_result <- dunnTest(NSF_LRRC37A3 ~ CN_NSF_17q2131_test_region_3, 
                             data = bowtie2_counts_genotypes, 
                             method = "bonferroni")
dunn_test_result

# # Y max
# 
# # Calculate the maximum value across both categories
# max_value <- max(max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "NSF_NSF"]),
#                  max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "NSF_LRRC37A3"]))
# 
# # Plotting with bar to mark mean and error bars for SEM, setting y-axis limits
# plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = CN_NSF_17q2131_test_region_3, y = Value)) +
#   geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
#   geom_bar(data = mean_se_values, aes(x = CN_NSF_17q2131_test_region_3, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
#   geom_errorbar(data = mean_se_values, aes(x = CN_NSF_17q2131_test_region_3, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
#   facet_grid(Category ~ ., scales = "free_y") +
#   labs(title = "Normalised Counts by Category",
#        y = "Normalised Counts",
#        x = "CN_NSF_17q2131_test_region_3") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
#         axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_line(size = 0.5)) +
#   guides(color = FALSE) +  # Remove the legend for the color aesthetic
#   ylim(0, max_value)  # Set y-axis limits
# 
# # Print the plot
# print(plot)
# ggsave("NSF-LRRC_fusion_nsfnorm_bowtie2_G,32,21_ymax.pdf", plot = plot)

# # NSF-LRRC37A3 with NSFP1 normalisation
# 
# # Select the columns you want to add from df_categorized
# columns_to_add <- c("CN_NSF_17q2131_test_region_3")
# 
# # Left join df_categorized to bowtie2_counts_df based on Prefix and Cell_Line, and select specific columns to add
# bowtie2_counts_genotypes <- left_join(bowtie2_counts_normalised, dplyr::select(df_categorized, Cell_Line, all_of(columns_to_add)), by = c("Prefix" = "Cell_Line"))
# 
# # Select the columns you want to use
# selected_columns <- c("Prefix", "NSFP1", "NSF_LRRC37A3_p1", "CN_NSF_17q2131_test_region_3")
# 
# # Subset the dataframe
# bowtie2_counts_genotypes_subset <- bowtie2_counts_genotypes[selected_columns]
# 
# # Rearrange the filtered dataframe to long format
# bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_subset %>%
#   pivot_longer(cols = -c(Prefix, CN_NSF_17q2131_test_region_3),
#                names_to = "Category",
#                values_to = "Value")
# 
# # Filter to include only "NSFP1" and "NSF_LRRC37A3" categories
# bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_long %>%
#   filter(Category %in% c("NSFP1", "NSF_LRRC37A3_p1"))
# 
# # Find the maximum value across both categories
# max_value <- max(max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "NSFP1"]),
#                  max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "NSF_LRRC37A3_p1"]))
# 
# # Calculate mean for each group
# mean_values <- bowtie2_counts_genotypes_long %>%
#   group_by(CN_NSF_17q2131_test_region_3, Category) %>%
#   summarise(mean_value = mean(Value, na.rm = TRUE))
# 
# # # Plotting
# # plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = CN_NSF_17q2131_test_region_3, y = Value)) +
# #   geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
# #   geom_point(data = mean_values, aes(x = CN_NSF_17q2131_test_region_3, y = mean_value, color = "red"), shape = "-", size = 10) +  # Add mean values as points
# #   facet_grid(Category ~ ., scales = "free_y") +
# #   labs(title = "Normalised Counts by Category",
# #        y = "Normalised Counts",
# #        x = "CN_NSF_17q2131_test_region_3") +
# #   theme_minimal() +
# #   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
# #         axis.title.x = element_blank(),
# #         panel.grid.major.x = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         panel.grid.major.y = element_line(size = 0.5)) +
# #   guides(color = FALSE)  # Remove the legend for the color aesthetic
# # 
# # # Print the plot
# # print(plot)
# 
# # Plotting with bar to mark mean
# plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = CN_NSF_17q2131_test_region_3, y = Value)) +
#   geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
#   geom_bar(data = mean_values, aes(x = CN_NSF_17q2131_test_region_3, y = mean_value), stat = "identity", fill = "red", alpha = 0.5) +  # Add bars for mean values
#   facet_grid(Category ~ ., scales = "free_y") +
#   labs(title = "Normalised Counts by Category",
#        y = "Normalised Counts",
#        x = "CN_NSF_17q2131_test_region_3") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
#         axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_line(size = 0.5)) +
#   guides(color = FALSE)  # Remove the legend for the color aesthetic
# 
# # Print/save the plot
# print(plot)
# ggsave("NSF-LRRC_fusion_nsfp1norm_bowtie2.pdf", plot = plot)
# 
# # Perform the Kruskal-Wallis test
# kruskal_test_result <- kruskal.test(NSF_LRRC37A3 ~ CN_NSF_17q2131_test_region_3, 
#                                     data = bowtie2_counts_genotypes)
# kruskal_test_result
# 
# kruskal_test_result_control <- kruskal.test(NSFP1 ~ CN_NSF_17q2131_test_region_3, 
#                                     data = bowtie2_counts_genotypes)
# kruskal_test_result_control
# 
# # Perform Dunn's test with Bonferroni correction
# dunn_test_result <- dunnTest(NSF_LRRC37A3 ~ CN_NSF_17q2131_test_region_3, 
#                              data = bowtie2_counts_genotypes, 
#                              method = "bonferroni")
# dunn_test_result
# 
# dunn_test_result_control <- dunnTest(NSFP1 ~ CN_NSF_17q2131_test_region_3, 
#                              data = bowtie2_counts_genotypes, 
#                              method = "bonferroni")
# dunn_test_result_control

# NAIP-OCLN normalised to NAIP

# Select the columns you want to add from df_categorized
columns_to_add <- c("CN_OCLNP1_region")

# Left join df_categorized to bowtie2_counts_df based on Prefix and Cell_Line, and select specific columns to add
bowtie2_counts_genotypes <- left_join(bowtie2_counts_normalised, dplyr::select(df_categorized, Cell_Line, all_of(columns_to_add)), by = c("Prefix" = "Cell_Line"))

# Select the columns you want to use
selected_columns <- c("Prefix", "NAIP_NAIP", "NAIP_OCLN", "CN_OCLNP1_region")

# Subset the dataframe
bowtie2_counts_genotypes_subset <- bowtie2_counts_genotypes[selected_columns]

# Rearrange the filtered dataframe to long format
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_subset %>%
  pivot_longer(cols = -c(Prefix, CN_OCLNP1_region), 
               names_to = "Category", 
               values_to = "Value")

# Filter to include only "NAIP_NAIP" and "NAIP_OCLN" categories
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_long %>%
  filter(Category %in% c("NAIP_NAIP", "NAIP_OCLN"))

# Calculate mean and SEM for each group
mean_se_values <- bowtie2_counts_genotypes_long %>%
  group_by(CN_OCLNP1_region, Category) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
            sem = sd(Value, na.rm = TRUE) / sqrt(n()))

# Plotting with bar to mark mean and error bars for SEM
plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = CN_OCLNP1_region, y = Value)) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
  geom_bar(data = mean_se_values, aes(x = CN_OCLNP1_region, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
  geom_errorbar(data = mean_se_values, aes(x = CN_OCLNP1_region, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
  facet_grid(Category ~ ., scales = "free_y") +
  labs(title = "Normalised Counts by Category",
       y = "Normalised Counts",
       x = "CN_OCLNP1_region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  guides(color = FALSE)  # Remove the legend for the color aesthetic

# Print the plot
print(plot)
ggsave("NAIP-OCLN_fusion_NAIPnorm_bowtie2_G,32,21.pdf", plot = plot)

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(NAIP_OCLN ~ CN_OCLNP1_region, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result

kruskal_test_result_control <- kruskal.test(NAIP_NAIP ~ CN_OCLNP1_region, 
                                            data = bowtie2_counts_genotypes)
kruskal_test_result_control

# Perform Dunn's test with Bonferroni correction
dunn_test_result <- dunnTest(NAIP_OCLN ~ CN_OCLNP1_region, 
                             data = bowtie2_counts_genotypes, 
                             method = "bonferroni")
dunn_test_result

dunn_test_result_control <- dunnTest(NAIP_NAIP ~ CN_OCLNP1_region, 
                             data = bowtie2_counts_genotypes, 
                             method = "bonferroni")
dunn_test_result_control

# NAIP-OCLN normalised to OCLN

# Select the columns you want to add from df_categorized
columns_to_add <- c("CN_OCLNP1_region")

# Left join df_categorized to bowtie2_counts_df based on Prefix and Cell_Line, and select specific columns to add
bowtie2_counts_genotypes <- left_join(bowtie2_counts_normalised, dplyr::select(df_categorized, Cell_Line, all_of(columns_to_add)), by = c("Prefix" = "Cell_Line"))

# Select the columns you want to use
selected_columns <- c("Prefix", "OCLN_OCLN", "NAIP_OCLN", "CN_OCLNP1_region")

# Subset the dataframe
bowtie2_counts_genotypes_subset <- bowtie2_counts_genotypes[selected_columns]

# Rearrange the filtered dataframe to long format
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_subset %>%
  pivot_longer(cols = -c(Prefix, CN_OCLNP1_region), 
               names_to = "Category", 
               values_to = "Value")

# Filter to include only "OCLN_OCLN" and "NAIP_OCLN" categories
bowtie2_counts_genotypes_long <- bowtie2_counts_genotypes_long %>%
  filter(Category %in% c("OCLN_OCLN", "NAIP_OCLN"))

# Calculate mean and SEM for each group
mean_se_values <- bowtie2_counts_genotypes_long %>%
  group_by(CN_OCLNP1_region, Category) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
            sem = sd(Value, na.rm = TRUE) / sqrt(n()))

# Plotting with bar to mark mean and error bars for SEM
plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = CN_OCLNP1_region, y = Value)) +
  geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
  geom_bar(data = mean_se_values, aes(x = CN_OCLNP1_region, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
  geom_errorbar(data = mean_se_values, aes(x = CN_OCLNP1_region, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
  facet_grid(Category ~ ., scales = "free_y") +
  labs(title = "Normalised Counts by Category",
       y = "Normalised Counts",
       x = "CN_OCLNP1_region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.5)) +
  guides(color = FALSE)  # Remove the legend for the color aesthetic

# Print the plot
print(plot)
ggsave("NAIP-OCLN_fusion_OCLNnorm_bowtie2_G,32,21.pdf", plot = plot)

# Perform the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(NAIP_OCLN ~ CN_OCLNP1_region, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result

kruskal_test_result_control <- kruskal.test(OCLN_OCLN ~ CN_OCLNP1_region, 
                                    data = bowtie2_counts_genotypes)
kruskal_test_result_control

# Perform Dunn's test with Bonferroni correction
dunn_test_result <- dunnTest(NAIP_OCLN ~ CN_OCLNP1_region, 
                             data = bowtie2_counts_genotypes, 
                             method = "bonferroni")
dunn_test_result

# # Y max
# 
# # Calculate the maximum value across both categories
# max_value <- max(max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "NAIP_NAIP"]),
#                  max(bowtie2_counts_genotypes_long$Value[bowtie2_counts_genotypes_long$Category == "NAIP_OCLN_naip"]))
# 
# # Plotting with bar to mark mean and error bars for SEM, setting y-axis limits
# plot <- ggplot(bowtie2_counts_genotypes_long, aes(x = CN_OCLNP1_region, y = Value)) +
#   geom_jitter(width = 0.3, height = 0, alpha = 0.7, shape = 1) +
#   geom_bar(data = mean_se_values, aes(x = CN_OCLNP1_region, y = mean_value), stat = "identity", fill = "red", alpha = 0.5, inherit.aes = FALSE) +  # Add bars for mean values
#   geom_errorbar(data = mean_se_values, aes(x = CN_OCLNP1_region, ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2, inherit.aes = FALSE) +  # Add error bars for SEM
#   facet_grid(Category ~ ., scales = "free_y") +
#   labs(title = "Normalised Counts by Category",
#        y = "Normalised Counts",
#        x = "CN_OCLNP1_region") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
#         axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_line(size = 0.5)) +
#   guides(color = FALSE) +  # Remove the legend for the color aesthetic
#   ylim(0, max_value)  # Set y-axis limits
# 
# # Print the plot
# print(plot)
# ggsave("NAIP-OCLN_fusion_naipnorm_bowtie2_G,32,21_ymax.pdf", plot = plot)


