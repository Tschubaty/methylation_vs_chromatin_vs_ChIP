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
###########################################################################

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

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Function to read and combine all chromosome files for a specific biosample
combine_chromosome_files <- function(biosample, input_dir) {
  file_paths <- file.path(input_dir,paste0(CHR_NAMES,"_", biosample, "_chromatin_histones.bed"))

  combined_data <- do.call(rbind, lapply(file_paths, function(file) {
    if (!file.exists(file)) {
      stop("File not found: ", file)
    }else{
      print(file)
    }
    
    data <- read.table(file, header = TRUE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
    data
  }))
  
  return(combined_data)
}

# Function to create a heatmap
create_heatmap <- function(data, output_file, biosample) {
  # Melt the data for ggplot
  melted_data <- melt(data, id.vars = c("Chromosome", "Start", "End", "CpG_Island_Status", "Strand", paste0("Chromatin_State_", biosample)),
                      variable.name = "Histone_Mark", value.name = "Presence")
  
  # Summarize the data to count the presence of histone marks in chromatin states
  summary_data <- melted_data %>%
    group_by_at(vars(paste0("Chromatin_State_", biosample), "Histone_Mark")) %>%
    summarize(Presence_Count = sum(Presence)) %>%
    ungroup() %>%
    rename(Chromatin_State = paste0("Chromatin_State_", biosample))
  
  # Create the heatmap
  heatmap_plot <- ggplot(summary_data, aes(x = Histone_Mark, y = Chromatin_State, fill = Presence_Count)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("Histone Mark Presence in Chromatin States for", biosample), fill = "Presence Count")
  
  # Save the plot
  ggsave(output_file, plot = heatmap_plot, width = 12, height = 8)
}

# Main script
biosample <- "K562"  # Change this to the desired biosample
input_dir <- "WGBS/byChr/histone_annotated"
output_file <- paste0(biosample, "_chromatin_histone_heatmap.png")

# Combine the data from all chromosome files
combined_data <- combine_chromosome_files(biosample, input_dir)

# Create and save the heatmap
create_heatmap(combined_data, output_file, biosample)

cat("Heatmap created and saved to", output_file, "\n")

