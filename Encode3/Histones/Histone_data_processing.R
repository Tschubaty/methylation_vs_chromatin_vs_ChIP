#################################################################
##  Chip-seq
##
##  input:   https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/
##  output: Encode3\
##  v 01 - 14.06.2023
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


################################################ INPUT ###################################################

sample_names <- c("K562", "HepG2")
encode_meta_file_name <- "Histone-data.xlsx"
input_folder <- file.path("download", "original")
output_folder <- "processed"
WGBS_folder <- "C:/Users/Daniel Batyrev/Documents/GitHub/WGBS_vs_ChIP_ENCODE3/Encode3/WGBS/merged"
################################################ INPUT ###################################################
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

###################################################################################
#
# https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/
#
#   E123 = K562
# E118 = HepG2

build_meta <- FALSE
if(build_meta){


file_names <-
  list.files(path = input_folder, pattern = ".narrowPeak")
# rename experiment names
experiment_names <-
  gsub(pattern = ".narrowPeak",
       replacement = "",
       x = file_names)
experiment_names <-
  gsub(pattern = "-",
       replacement = "_",
       x = experiment_names)
experiment_names <-
  gsub(pattern = "E118",
       replacement = "HepG2",
       x = experiment_names)
experiment_names <-
  gsub(pattern = "E123",
       replacement = "K562",
       x = experiment_names)

names(file_names) <- experiment_names

meta <-
  data.frame(
    sample = gsub("\\_.*", "", experiment_names),
    target = gsub(pattern = ".*_", replacement =  "", experiment_names),
    file_name =  file_names
  )

meta$n_peaks <- NA
meta$covergage <- NA



for(f in 1:nrow(meta)){
  file_name <- file.path(input_folder,meta$file_name[f])
  print(file_name)
  histon_file <-read.table(file = file_name,col.names = c("chr","start","end" ,"name","score","strand","signalValue" ,"pValue", "qValue",  "peak"))
  df_histon <- histon_file[histon_file$chr %in% CHR_NAMES,]
  colName <- paste(meta$sample[f],meta$target[f],sep = "_")
  meta$n_peaks[f] <- nrow(df_histon)
  meta$covergage[f] <- sum(df_histon$end - df_histon$start)
}
saveRDS(object = meta,file = file.path(output_folder,"meta.rds"))

}else{
  meta <- readRDS(file = file.path(output_folder,"meta.rds"))
}

for(chr in CHR_NAMES){
  start_time <- Sys.time()


  # chr <- CHR_NAMES[1]
  print(chr)
  df_meth <- readRDS(file = file.path(WGBS_folder,paste(chr,"K562_vs_HepG2.merged_CpG_methylation.RDS",sep = ".")))
  df_meth[,rownames(meta) ] <- NA
  
  for(f in 1:nrow(meta)){
    start_hist <- Sys.time()
    file_name <- file.path(input_folder,meta$file_name[f])
    print(file_name)
    histon_file <-read.table(file = file_name,col.names = c("chr","start","end" ,"name","score","strand","signalValue" ,"pValue", "qValue",  "peak"))
    histon_file_chr <- histon_file[histon_file$chr == chr, c("chr","start","end" ,"score")]
    
    #pb = txtProgressBar(min = 1, max = nrow(histon_file_chr), initial = 1) 
    #histon_file_chr <- head(histon_file_chr) ###### debug
    for(r in 1:nrow(histon_file_chr)){
      #setTxtProgressBar(pb,r)
      peak <- histon_file_chr[r,]
      df_meth[peak$start <= df_meth$start & df_meth$end <= peak$end,rownames(meta)[f]] <- peak$score
      #print(r/nrow(histon_file_chr))
    }
    end_hist <- Sys.time()
    print("hist time:")
    print(end_hist - start_hist)
  }
  saveRDS(object = df_meth,file = file.path(output_folder,paste(chr,"WGBS","histones","rds",sep = ".")))
  end_time <- Sys.time()
  print("chromosome time:")
  print(end_time - start_time)
}

