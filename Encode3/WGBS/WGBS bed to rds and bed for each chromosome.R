#################################################################
##  Parse WGBS bed files in small chromosome rds files
##
##  input:  bed file of whole genomes
##  output: rds files for each chromosome + QC plots
##  v04 - 19.06.2023
##  Author: Daniel Batyrev 777634015
#################################################################

#Clear R working environment
rm(list = ls())

this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)

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

########################### input ##########################

#libs:
library(ggplot2)

# Input file should be in the follwing format:
# ENCODE3 methylBED
#
# Reference chromosome or scaffold
# Start position in chromosome
# End position in chromosome
# Name of item
# Score from 0-1000. Capped number of reads
# Strandedness, plus (+), minus (-), or unknown (.)
# Start of where display should be thick (start codon)
# End of where display should be thick (stop codon)
# Color value (RGB)
# Coverage, or number of reads
# Percentage of reads that show methylation at this position in the genome

output_folder <- "out"

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
#######################################################################################
# WGBS files
input_files <-
  c(
    file.path(
      "WGBS K562 GEOGSE86747 ENCODE4 v1.1.6 GRCh38 (ENCAN334EDE) processed data",
      "ENCFF660IHA.bed"
    ),
    file.path(
      "WGBS K562 GEOGSE86747 ENCODE4 v1.1.6 GRCh38 (ENCAN334EDE) processed data",
      "ENCFF328NMN.bed"
    ),
    file.path(
      "WGBS HepG2 GEOGSE127318 ENCODE4 v1.1.6 GRCh38 (ENCAN831UDN) processed data",
      "ENCFF820ATI.bed"
    ),
    file.path(
      "WGBS HepG2 GEOGSE127318 ENCODE4 v1.1.6 GRCh38 (ENCAN831UDN) processed data",
      "ENCFF690FNR.bed"
    ),
    file.path(
      "WGBS HepG2 GEOGSE86764 ENCODE4 v1.1.6 GRCh38 (ENCAN888RDF) processed data",
      "ENCFF817LMT.bed"
    ),
    file.path(
      "WGBS HepG2 GEOGSE86764 ENCODE4 v1.1.6 GRCh38 (ENCAN888RDF) processed data",
      "ENCFF453UDK.bed"
    )
  )
#
data_names <-
  gsub(pattern = ".bed",
       replacement = "",
       x = input_files)
data_names <-
  gsub(pattern = "GEO.*.data/",
       replacement = "",
       x = data_names)
data_names <- gsub(pattern = " ",
                   replacement = "_",
                   x = data_names)

############################ run script ####################
#f <- 1
for (f in 1:length(input_files)) {
  
  input_file <- input_files[f]
  
  WGBS <- read.table(file = input_file)
  
  # Each column represents the following:
  #
  # 1- Reference chromosome or scaffold
  # 2- Start position in chromosome
  # 3- End position in chromosome
  # 4- Name of item
  # 5- Score from 0-1000. Capped number of reads
  # 6- Strandedness, plus (+), minus (-), or unknown (?)
  # 7 -Start of where display should be thick
  # 8 -End of where display should be thick
  # 9 -Color value (RGB)
  # 10 Coverage, or number of reads
  # 11 - Percentage of reads that show methylation at this position in the genome
  # 12 - Reference genotype
  # 13 - Sample genotype
  # 14 - Quality score for genotype call
  
  
  colnames(WGBS) <- c(
    "chr",
    "start",
    "end",
    "name",
    "K1cappedReads",
    "strand",
    "startDisplay",
    "endDisplay",
    "RGB",
    "readCoverage",
    "percentMethylated",
    "RefGenotype",
    "SampleGenotype",
    "qScore"
  )
  
  
  plot_whole_Genome <- TRUE
  if(plot_whole_Genome){
    myplot <- ggplot(data = WGBS,aes(x=K1cappedReads)) +
      geom_histogram(breaks = c(1:100))
    ggsave(filename = file.path(output_folder,paste(input_file, "_read_distribution.png",sep = "")))
    
    myplot <-ggplot(data = WGBS,aes(x=percentMethylated)) +
      geom_histogram(binwidth=1)
    ggsave(filename = file.path(output_folder,paste(input_file, "_methyl_distribution.png",sep = "")))
}
  
  # debug chr <- CHR_NAMES[1]
  for (chr in CHR_NAMES) {
    print(chr)
    chr_df <-  WGBS[WGBS$chr == chr, ]
    rds_file_name <- paste(chr, data_names[f], "rds", sep = ".")
    saveRDS(object = chr_df[, c(
      "chr",
      "start",
      "end",
      "name",
      "K1cappedReads",
      "strand",
      "readCoverage",
      "percentMethylated"
    )],
    file = file.path(output_folder, rds_file_name))
    print(rds_file_name)
  }
}

############################## merge #################################################

# debug chr <- CHR_NAMES[1]
for (chr in CHR_NAMES) {
  # chr <- CHR_NAMES[length(CHR_NAMES)]
  print(chr)
  rds_files <-
    list.files(path = "out", pattern = paste(chr, "WGBS", sep = "."))
  
  for (rds_file_name in rds_files) {
    # rds_file_name <- rds_files[1]
    # rds_file_name <- rds_files[2]
    print(rds_file_name)
    df_load <- readRDS(file = file.path("out", rds_file_name))
    colnames(df_load)[which(colnames(df_load) == 'K1cappedReads')] <-
      'score'
    df_load[, paste(
      "Methylated",
      gsub(
        pattern = ".*.WGBS_|.rds",
        replacement = "",
        x = rds_file_name
      ),
      sep = "_"
    )] <- round(df_load$readCoverage * df_load$percentMethylated / 100)
    df_load[, paste(
      "UnMethylated",
      gsub(
        pattern = ".*.WGBS_|.rds",
        replacement = "",
        x = rds_file_name
      ),
      sep = "_"
    )] <- round(df_load$readCoverage * (1 - df_load$percentMethylated / 100))
    df_load <-
      df_load[, !colnames(df_load) %in% c("name", "score", "readCoverage", "percentMethylated")]
    
    if (rds_file_name == rds_files[1]) {
      df_merge <- head(df_load)
    } else{
      df_merge <-
        merge(
          x = df_merge,
          y = head(df_load),
          by = c("chr", "start", "end", "strand"),
          all.x = TRUE,
          all.y = TRUE
        )
    }
    
    
    # rename new columns
    # for(r in 4:ncol(df_load)){
    #   # r <- 4
    #   new_col_name <- paste(colnames(df_load)[r],gsub(pattern = ".*.WGBS_|.rds",replacement = "",x = rds_file_name),sep = "_")
    #   colnames(df_load)[which(colnames(df_load) == colnames(df_load)[r])] <- new_col_name
    # }
    #
  }
  
 merge_file_name <- paste(chr, "full", "rds", sep = ".")
  saveRDS(object = df_merge,
  file = file.path(output_folder, merge_file_name))
  print(merge_file_name)
}


