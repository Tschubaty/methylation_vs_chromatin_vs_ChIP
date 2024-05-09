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


##########################################################################

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

s <- 1
n_file <- 1
chr <- CHR_NAMES[22]

# Loop through each file
for(s in 1:length(biosamples)) {
  sample <- biosamples[s]
  file_list <- list.files(path = sample, pattern = "\\.rds$")
  file_list <- file_list[!grepl("chr", files)]
  
  
  for(n_file in 1:length(file_list)) {
    file <- file_list[n_file]
    file_name <- file.path(sample, file)
    print(file_name)
    # Read the data from the file
    data <- readRDS(file_name)
    for(chr in CHR_NAMES) {
      chr_data <- data[data$chr == chr,]
      chr_data$m_read <-
        round(chr_data$n_read * chr_data$f_read / 100)
      colnames(chr_data)[7:9] <-
        paste(
          colnames(chr_data)[7:9],
          sample,
          gsub(
            pattern = ".rds",
            replacement = "",
            x = file
          ),
          sep = "_"
        )
      saveRDS(object = chr_data, file = file.path(sample, paste(chr, sample, file, sep = ".")))
    }
  }
}   
######################################################
for(chr in CHR_NAMES) {
  chr_data <- readRDS(file = paste0(chr, ".WGBS.rds"))
  print(chr)
  # Define the columns you want to replace NA values with 0
  columns_to_replace_na <-
    colnames(chr_data)[grepl("^n_read_|^m_read_", colnames(chr_data))]
  
  # Replace NA values with 0
  chr_data <- chr_data %>%
    mutate_at(vars(columns_to_replace_na), ~ ifelse(is.na(.), 0, .))
  
  # # Define the biosamples you want to sum up
  # biosamples <- c("HepG2", "K562", "A549", "GM12878")
  
  # Create new columns by summing up the values for each biosample
  for (sample in biosamples) {
    n_read_cols <-
      grep(paste0("^n_read_", sample), colnames(chr_data), value = TRUE)
    m_read_cols <-
      grep(paste0("^m_read_", sample), colnames(chr_data), value = TRUE)
    
    chr_data <- chr_data %>%
      mutate(
        !!paste0("nRead_", sample) := rowSums(select(., all_of(n_read_cols)), na.rm = TRUE),!!paste0("mRead_", sample) := rowSums(select(., all_of(m_read_cols)), na.rm = TRUE)
      )
    chr_data <- chr_data %>%
      mutate(!!paste0("fRead_", sample) := get(paste0("mRead_", sample)) / get(paste0("nRead_", sample)))
  }
  # Select columns
  selected_data <- chr_data %>%
    select(chr, start, end, name, score, strand, starts_with("nRead_"), starts_with("fRead_"), starts_with("mRead_"))
  
  # Save the selected data to an RDS file
  saveRDS(selected_data, file = paste(chr,"merged","WGBS","rds",sep = "."))
}
