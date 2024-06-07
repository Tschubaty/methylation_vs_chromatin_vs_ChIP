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

  combined_data <- do.call(rbind, lapply(CHR_NAMES, function(chr) {
    # chr <- CHR_NAMES[22]
    file <- file.path(input_dir, paste0(chr, "_cpgs_island_", biosample, "_chromatin_histones.bed"))
    if (file.exists(file)) {
      print(file)
      # Read the header separately
      header <- readLines(file, n = 1)
      col_names <- strsplit(header, " ")[[1]]
      col_names[[1]] <- "chr"
      
      data <- read.table(file, header = FALSE, skip = 0, sep = "\t", stringsAsFactors = FALSE)
      colnames(data) <- col_names
      return(data)
    } else {
      stop("File not found: ", file)
      return(NULL)
    }
  }))
  
  if (is.null(combined_data)) {
    stop("No files found for the biosample: ", biosample)
  }
  
  return(combined_data)
}



# Main script
biosample <- "K562"  # Change this to the desired biosample
input_dir <- "WGBS/byChr/histone_annotated"
output_file <- paste0(biosample, "_chromatin_histone_heatmap.png")

# Combine the data from all chromosome files
combined_data <- combine_chromosome_files(biosample, input_dir)



colnames_state_summery <- c("biosample",names(combined_data)[6:length(names(combined_data))],"n_CpG")
state_summery <- data.frame(matrix(ncol = length(colnames_state_summery ), nrow = max(combined_data$Chromatin_State)))
colnames(state_summery) <- colnames_state_summery
state_summery$biosample <- biosample


for(state in 1:18) {
  print(state)
  state_data <- combined_data[combined_data$Chromatin_State == state, ]
  state_summery$Chromatin_State[state] <- state
  state_summery$n_CpG[state] <- nrow(state_data)
  
  for (hist in names(combined_data)[7:length(names(combined_data))]) {
    print(hist)
    state_summery[state, hist] <- mean(state_data[, hist])
    
  }
}

state_summery_long$Histone_Mark <-
  factor(state_summery_long$Histone_Mark, levels = sort(as.character(
    unique(state_summery_long$Histone_Mark)
  )))

state_summery_long <- melt(state_summery, id.vars = c("biosample", "Chromatin_State","n_CpG"),
     variable.name = "Histone_Mark", value.name = "Occurrence_Probability")

saveRDS(file = file.path("Histones",paste0(biosample,"_histone_summery.rds")),object = state_summery)

heatmap_plot <- ggplot(state_summery_long, aes(x = Histone_Mark, y = Chromatin_State, fill = Occurrence_Probability)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  labs(title = paste("Histone Mark Presence Probability in Chromatin States for", biosample), fill = "Probability")