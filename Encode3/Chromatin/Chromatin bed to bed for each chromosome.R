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

sample_names <- c("K562", "HepG2")

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
# file
# E123 = K562
#E118 = HepG2

sample_names <- c("K562", "HepG2")

input_files <- c(K562 = "E123_18_core_K27ac_hg38lift_dense.bed",HepG2 = "E118_18_core_K27ac_hg38lift_dense.bed")
#
output_folder <- "by_chr"

# data_names <-
#   gsub(pattern = ".bed",
#        replacement = "",
#        x = input_files)
# data_names <-
#   gsub(pattern = "GEO.*.data/",
#        replacement = "",
#        x = data_names)
# data_names <- gsub(pattern = " ",
#                    replacement = "_",
#                    x = data_names)

############################ run script ####################
#f <- 1
for (f in 1:length(input_files)) {
  
  input_file <- input_files[f]
  
  bed <- readRDS("combined_annaotation.rds")
  
  #bed <- read.table(file = input_file,row.names = FALSE,col.names = FALSE,skip = 5)
  
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
  
  
  colnames(bed) <- c(
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
    
    names(df_load)[names(df_load) == 'percentMethylated'] <- paste(
      "percentMethylated",
      gsub(
        pattern = ".*.WGBS_|.rds",
        replacement = "",
        x = rds_file_name
      ),
      sep = "_"
    )
    
    
    df_load <-
      df_load[, !colnames(df_load) %in% c("name", "score", "readCoverage", "percentMethylated")]
    
    if (rds_file_name == rds_files[1]) {
      df_merge <- df_load
    } else{
      df_merge <-
        merge(
          x = df_merge,
          y = df_load,
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
  
 full_file_name <- paste(chr, "full", "rds", sep = ".")
  saveRDS(object = df_merge,
  file = file.path(output_folder, full_file_name))
  print(full_file_name)
}

############################## merge #################################################


# debug chr <- CHR_NAMES[1]
for (chr in CHR_NAMES) {
  # chr <- CHR_NAMES[length(CHR_NAMES)]
  print(chr)
  
  full_file_name <- paste(chr, "full", "rds", sep = ".")
  df_full <-  readRDS(file = file.path(output_folder, full_file_name))
  df_merge <- df_full[,c(1:4)]
  
  for(s in sample_names){
    # s <- sample_names[1]
    UnMethylated <- rowSums(x = df_full[,grep(pattern = paste("UnMethylated",s,sep = "_"),x = colnames(df_full))],na.rm = TRUE)
    methylated <- rowSums(x = df_full[,grep(pattern = paste("Methylated",s,sep = "_"),x = colnames(df_full))],na.rm = TRUE)
    df_merge[,paste("fRead",s,sep = "_")] <- methylated/(methylated + UnMethylated)
    df_merge[,paste("nRead",s,sep = "_")] <- methylated + UnMethylated
    
  }
  full_file_name <- paste(chr, "merge", "rds", sep = ".")
  saveRDS(object = df_merge,
          file = file.path(output_folder, full_file_name))
  print(full_file_name)
}

##############################plot full #######################

# debug chr <- CHR_NAMES[1]
for (chr in CHR_NAMES) {
  # chr <- CHR_NAMES[length(CHR_NAMES)]


  full_file_name <- paste(chr, "full", "rds", sep = ".")
  print(chr)

  df <- readRDS(file = file.path(output_folder, full_file_name))
  print(full_file_name)

 # long fomrat for plotting 
  col_index <- grepl(pattern = "chr|start|end|strand|percentMethylated",x = colnames(df))
  df_long <-reshape2::melt(data = df[,col_index], id.vars = c("chr","start","end","strand"), variable.name = "sample",value.name = "fRead")

p <- ggplot(data = df_long,mapping = aes(x = sample,y = fRead,fill = sample))+geom_boxplot()+ggtitle(chr)
ggsave(filename = file.path("plot",paste(chr,"boxplots","png",sep = ".")),plot = p,width = 10,height = 10)
}


# 
# names(df)[5:ncol(df)]


     
     
     
     
     
     
     