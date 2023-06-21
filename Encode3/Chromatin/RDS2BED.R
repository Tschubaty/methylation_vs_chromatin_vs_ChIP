#################################################################
##  #RDS  to bed
##
##  input: RDS
##
##
##  output: BED
##  v_01 01.09.2021
##  Author: Daniel Batyrev 777634015
#################################################################
#Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"
} else{
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
}

setwd(this.dir)

output_folder <- "by_chr"

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

df <- readRDS(file = "combined_annaotation.RDS")

for (chr in CHR_NAMES) {
  write.table(
    x = df[df$chr == chr, ],
    file = file.path(output_folder, paste(chr, "chromatin", "bed", sep =
                                            ".")),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
}