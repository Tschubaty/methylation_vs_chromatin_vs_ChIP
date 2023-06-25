#################################################################
##  samppe methylaTION DIFFERENCE
##
##  input: everything in HG38
##
##          Encode3\
##
##
##  output:   Encode3\
##
##  v_1
##  Author: Daniel Batyrev 777634015
#################################################################
#Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"
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
#################################################################
# E123 = K562
# E118 = HepG2

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
sample_names <- c("K562", "HepG2")
##########################   library   #######################################

# library("ggplot2")
library(foreach)
library(doParallel)
n_cores <- detectCores() - 1
########################## code ############################

input_folder_data <- "Chromatin/out"
plot_folder <- "Chromatin/plot"
output_folder <- plot_folder <- "Chromatin/meta"
  

file_names <- paste(CHR_NAMES,"WGBS.histones.chromatin.RDS",sep = ".")

# field	example	SQL type	info	description
# bin	585	smallint(6)	range	Indexing field to speed chromosome range queries.
# chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
# chromStart	28735	int(10) unsigned	range	Start position in chromosome
# chromEnd	29737	int(10) unsigned	range	End position in chromosome
# name	CpG: 111	varchar(255)	values	CpG Island
# length	1002	int(10) unsigned	range	Island Length
# cpgNum	111	int(10) unsigned	range	Number of CpGs in island
# gcNum	731	int(10) unsigned	range	Number of C and G in island
# perCpg	22.2	float	range	Percentage of island that is CpG
# perGc	73	float	range	Percentage of island that is C or G
# obsExp	0.85	float	range	Ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island
CpG_bed_colum_names <-
  c(
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "length",
    "cpgNum",
    "gcNum",
    "perCpg",
    "perGc",
    "obsExp"
  )
CpG_file_name <- "WGBS/islands/cpgIslandExt.hg38.bed"
island_bed <- read.csv(
  file = CpG_file_name,
  header = FALSE,
  sep = "\t",
  col.names = CpG_bed_colum_names
)


# loop over every chrom to annonated mehtylation value to island 
for (chr in CHR_NAMES) {
  print(chr)
  # debug chr <- CHR_NAMES[21]
  chrom_islands <- island_bed[island_bed$chrom == chr,]
  chrom_df <-
    readRDS(file = file.path(
      input_folder_data,
      paste(chr, "WGBS.histones.chromatin.RDS", sep = ".")
    ))
  colnames(chrom_df)[colnames(chrom_df) == "chr.1"] <-
    "cpgNum_Island"
  chrom_df$cpgNum_Island <- 0
  
  
  # vecorize wgbs
  start <- chrom_df$start
  end <- chrom_df$end
  cpgNum_Island <- chrom_df$cpgNum_Island
  # vevorzie island bed
  chromStart <- chrom_islands$chromStart
  chromEnd <- chrom_islands$chromEnd
  cpgNum <- chrom_islands$cpgNum
  
  
  # start loop over each island
  for (i in 1:length(chromStart)) {
    if (i %% 100 == 0) {
      print(c("i=", i, "completion =", i / length(chromStart)))
    }
    cpgNum_Island[chromStart[i] <= start &
                    end <= chromEnd[i]] <- cpgNum[i]
  }
  # back to df
  chrom_df$cpgNum_Island <- cpgNum_Island
  
  
  saveRDS(object = chrom_df,
          file = file.path("WGBS",
            "complete",
            paste(chr, "WGBS.histones.chromatin.islands.RDS", sep = ".")
          ))
}