#################################################################
##  Chip-seq
##
##  input: Encode3\ChIP-seq\2. filter samples
##  output: Encode3\ChIP-seq\3. merge samples
##  v 03 - 26.07.2023
##  Author: Daniel Batyrev 777634015
#################################################################

#Clear R working environment
rm(list = ls())

this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"

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
library(gridExtra)
library(ggVennDiagram)

# set up parallel processing
n_cores <- detectCores() - 1
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


# function definitions

merge_sorted_bed <- function(bed, cut_off) {
  #for each row look at dist and merge with previous if dist <= cut off
  bed$state <- NA
  bed$hits <- NA
  bed$merge_peak <- NA
  rows <- c()
  states <- c()
  peak_coor <- c()
  for (r in 1:nrow(bed)) {
    row_pos <- bed[r,]
    # print(row_pos)
    # print(rows)
    # print(states)
    # print(peak_coor)
    if (abs(row_pos$distance) <= OVERLAP_CUTTOF) {
      # memorize rows
      rows <-  c(rows, r)
      states <- c(states, bed$sample[r])
      peak_coor <- c(peak_coor, bed$peak_coor[r])
      # in case of last row write it donw
      if (r == nrow(bed)) {
        # merge all previous
        if (length(unique(states)) > 1) {
          bed$state[rows] <- "both"
        } else{
          bed$state[rows] <- states
        }
        bed$merge_peak[rows] <- mean(peak_coor)
        bed$hits[rows] <-  length(rows)
      }
    }else{
      # merge all previous
      if (length(unique(states)) > 1) {
        bed$state[rows] <- "both"
      } else{
        bed$state[rows] <- states
      }
      bed$merge_peak[rows] <- mean(peak_coor)
      bed$hits[rows] <-  length(rows)
      
      # start new potential merging
      if (r < nrow(bed)) {
        rows <- c(r)
        states <- c(bed$sample[r])
        peak_coor <- c(bed$peak_coor[r])
      } else{
        # if not last row
        bed$state[r] <- bed$sample[r]
        bed$merge_peak[r] <- bed$peak_coor[r]
        bed$hits[r] <- 1
      }
    }
  }
  return(bed)
}

################################################ INPUT ###################################################
sample_names <- c("HepG2", "GM12878") # sample_names <- c("K562", "HepG2", "A549", "GM12878")
OVERLAP_CUTTOF <- 100
input_folder <- "2. filter samples" 
output_folder <- "3. merge samples"
################################################ INPUT ###################################################
dir.create(file.path(output_folder, paste(sample_names,collapse = "_")), showWarnings = FALSE)
output_folder  <- file.path(output_folder,paste(sample_names,collapse = "_"))
################################################ CODE ###################################################
#Read Input 

proteins <- list.dirs(path = input_folder,full.names = FALSE)


#bed_files <- list.files(path = fil input_folder, pattern = ".bed")

for (target in proteins) {
  # target <- proteins[which(proteins == "CTCF")]
  print(target)
  bed_all <- data.frame()
  protein_files <- list()
  for (s in sample_names) {
    # s <- sample_names[1]
    print(s)
    protein_files[[s]] <-
      list.files(path = file.path(input_folder,target), pattern = paste(s, target, sep = "_"))
    for (file_name in protein_files[[s]]) {
      # debug file_name <- protein_files[[s]][1]
      print(file_name)
      df_temp <-
        read.delim(
          file = file.path(input_folder,target, file_name),
          header = FALSE,
          col.names = c(
            "chr",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "signalValue",
            "pValue",
            "qValue",
            "peak"
          )
        )
      df_temp$sample <- s
      df_temp$file <-
        gsub(pattern = ".*_|.bed",
             replacement = "",
             x = file_name)
      bed_all <- rbind(bed_all, df_temp)
    }
  }
  
  if(length(unique(bed_all$sample)) > 1){
  
  # sort bed
  bed_all$chr <-  as.character(bed_all$chr)
  bed_all$peak_coor <- bed_all$start + bed_all$peak
  bed_all <-
    bed_all[order(bed_all$chr, bed_all$peak_coor),]
  
  # merge all bed files
  final_bed <- data.frame()
  for (c in 1:length(CHR_NAMES)) {
    print(CHR_NAMES[c])
    df_chr <- bed_all[bed_all$chr == CHR_NAMES[c], ]
    if (nrow(df_chr) > 0) {
      # compute peak distances 
      df_chr$distance <-
        c(Inf, df_chr$peak_coor[2:nrow(df_chr)] - df_chr$peak_coor[1:nrow(df_chr) -1])
      df_chr <- merge_sorted_bed(df_chr, OVERLAP_CUTTOF)
      final_bed <- rbind(final_bed, df_chr)
    }
  }
  saveRDS(object = final_bed, file = file.path(output_folder, paste(target,paste(sample_names,collapse = "_"), "RDS", sep = ".")))
  }
}
