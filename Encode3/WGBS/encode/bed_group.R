# Documentation -----------------------------------------------------------
##
##
##  Author: Daniel Batyrev (HUJI 777634015)
##
# Set up Work Environment --------------------------------------------------

#Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/HumanEvo/HumanEvo/"
  picuture_file_extension <- "pdf"
} else{
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picuture_file_extension <- "png"
}
setwd(this.dir)
detachAllPackages <- function() {
  basic.packages <-
    c(
      "package:stats",
      "package:graphics",
      "package:grDevices",
      "package:utils",
      "package:datasets",
      "package:methods",
      "package:base"
    )
  package.list <-
    search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]
  package.list <- setdiff(package.list, basic.packages)
  if (length(package.list) > 0)
    for (package in package.list)
      detach(package, character.only = TRUE)
}
detachAllPackages()
# MOVE files  --------------------------------------------------
# Load required library
library(data.table)

# Read the metadata.tsv file
metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Print the first few rows of the metadata to check if it's read correctly
head(metadata)

# Create a list of unique Biosamples
biosamples <- unique(metadata$Biosample.term.name)

# Create folders for each Biosample
for (biosample in biosamples) {
  dir.create(biosample, showWarnings = FALSE)
}

# Loop through each file and move it to the corresponding folder
for (file in list.files(pattern = ".bed")) {
  print(file)
  # Extract Biosample from metadata.tsv based on filename
  file_accession <- sub(".bed", "", file)
  biosample <- metadata$Biosample.term.name[metadata$File.accession == file_accession]
  
  # Move the file to the corresponding folder
  if (!is.na(biosample)) {
    file.rename(file, file.path(biosample, file))
  } else {
    print(paste("Biosample not found for file:", file))
  }
}

#save RDS ############################################################

bedColnames <- c(
"chr",
"start",
"end",
"name",
"score",
"strand",
"StartThick",
"EndThick",
"color",
"n_read",
"f_read",
"Reference_genotype",
"Sample_genotype",
"Quality_score")

for(folder in biosamples) {
  for (file in list.files(path = folder, pattern = ".bed")) {
    file_name  <- file.path(folder, file)
    print(file_name)
    bed <-
      read.table(
        file_name,
        header = FALSE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
    colnames(bed) <- bedColnames
    bed$m_read <- bed$n_read * (1 - bed$f_read / 100)
    saveRDS(
      object = bed[, c("chr",
                       "start",
                       "end",
                       "name",
                       "score",
                       "strand",
                       "n_read",
                       "f_read",
                       "m_read")],
      file = gsub(
        pattern = ".bed",
        replacement = ".rds",
        x = file_name
      )
    )
    print(dim(bed))
  }
}

# merge by chr and sample  ################################################################
library(dplyr)

# Initialize an empty data frame to store the aggregated data
aggregated_data <- data.frame(chr = character(), start = integer(), end = integer(), 
                              name = integer(), score = integer(), strand = character(),
                              nRead_HepG2 = integer(), nRead_K562 = integer(), 
                              nRead_A549 = integer(), nRead_GM12878 = integer(),
                              fRead_HepG2 = numeric(), fRead_K562 = numeric(),
                              fRead_A549 = numeric(), fRead_GM12878 = numeric(),
                              m_read = integer(), stringsAsFactors = FALSE)

# Loop through each biosample folder
for (folder in biosamples) {
  print(folder)
  # Initialize counters for each biosample
  nRead_biosample <- 0
  fRead_biosample <- 0
  
  # Initialize an empty data frame for the current biosample
  biosample_data <- data.frame()
  
  # Loop through each file in the biosample folder
  for (file in list.files(path = folder, pattern = ".rds")) {
    print(file)
    # Read the data from the file
    data <- readRDS(file.path(folder, file))
    
    # Add nRead for the current file to the biosample total
    nRead_biosample <- nRead_biosample + sum(data$n_read)
    
    # Add fRead for the current file to the biosample total (weighted sum)
    fRead_biosample <- fRead_biosample + sum(data$n_read * data$f_read) / sum(data$n_read)
    
    # Merge data with biosample_data using left join
    biosample_data <- left_join(biosample_data, data, by = c("chr", "start", "end", "strand"))
  }
  
  # Calculate fRead for the biosample (percentage of methylated reads)
  fRead_biosample <- fRead_biosample / nRead_biosample
  
  # Add biosample data to aggregated_data
  aggregated_data <- bind_rows(aggregated_data, biosample_data)
  
  # Add total nRead and fRead for the biosample to aggregated_data
  aggregated_data[nrow(aggregated_data), paste0("nRead_", basename(folder))] <- nRead_biosample
  aggregated_data[nrow(aggregated_data), paste0("fRead_", basename(folder))] <- fRead_biosample
}

# Remove duplicated rows
aggregated_data <- distinct(aggregated_data)

# View the aggregated data
head(aggregated_data)




##############################################################################################
# rm(bed)
# gc()

# Get unique chromosomes
CHR_NAMES <-
  c(
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22"
  )

# Loop through each chromosome
for (chr in unique(data$chr)) {
  # Subset data for the current chromosome
  data_chr <- data[data$chr == chr, ]
  
  # Ensure consistent column names between chr_merged_data and data_chr
  colnames(data_chr) <- colnames(chr_merged_data)
  
  # Merge the data frames by common columns: chr, start, end, and strand
  chr_merged_data <-
    merge(
      chr_merged_data,
      data_chr,
      by = c("chr", "start", "end", "strand"),
      all = TRUE
    )
}

library(dplyr)

# Loop through each file
for(s in 1:length(biosamples) {
  sample <- biosamples[s]
  file_list <- list.files(path = sample, pattern = ".rds")
  
  for(n_file in 1:length(file_list) {
    file <- file_list[n_file]
    file_name <- file.path(sample, file)
    print(file_name)
    # Read the data from the file
    data <- readRDS(file_name)

    if(n_file == 1){
      merged_data <- data
    }else{
      
    }
   
    merged_data <- merged_data %>%
      mutate(
        nRead = ifelse(is.na(n_read.y), n_read.x, n_read.x + n_read.y),
        fRead = ifelse(
          is.na(n_read.y),
          round(n_read.x * f_read.x / 100),
          (round(n_read.x * f_read.x / 100) + round(n_read.y * f_read.y / 100)) / nRead
        ),
        score = ifelse(is.na(score.y), score.x, score.x + score.y)
      ) %>%
      select(chr, start, end, name, score, strand, nRead, fRead, m_read) %>%
      distinct() # Remove duplicated rows
    
    # Save the merged data frame for the current sample
    saveRDS(object = chr_merged_data, file = file.path(sample, paste0(sample, "_merged.rds")))
  }
}

# Loop through each chromosome
for (chr in CHR_NAMES) {
  print(chr)
  # Subset data for the current chromosome
  chr_data <- merged_data %>% filter(chr == .data$chr)
  
  # Group by chromosome and start position, summarize reads according to biosample
  chr_data <- chr_data %>%
    group_by(chr, start) %>%
    summarise(
      across(starts_with("n_read_"), ~ sum(.x)),
      across(starts_with("f_read_"), ~ sum(.x * n_read) / sum(n_read)),
      across(starts_with("m_read_"), ~ sum(.x * n_read) / sum(n_read)),
      .groups = "drop"
    )
  
  # Recalculate f_read column
  chr_data <- chr_data %>% 
    mutate(across(starts_with("f_read_"), ~ round(. * 100, digits = 2)))
  
  # Save the result to a file
  saveRDS(object = chr_data,file = paste0(chr, "_merged_WGBS.rds"))
  cat("Chromosome", chr, "processed and saved as", file_name, "\n")
}
