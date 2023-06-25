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

df <- paste(CHR_NAMES,"WGBS.histones.chromatin.RDS",sep = ".")

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
  # debug chr <- CHR_NAMES[1]
  chrom_islands <- island_bed[island_bed$chrom == chr, ]
  chrom_df <-
    readRDS(file = file.path(
      input_folder,
      paste(chr, "K562_vs_HepG2.merged_CpG_methylation.RDS", sep = ".")
    ))
  colnames(chrom_df)[colnames(chrom_df) == "chr.1"] <- "CpG_island"
  chrom_islands$mean_meth_HepG2 <- NA
  chrom_islands$mean_meth_K562 <- NA
  
  # start loop over each island
  for (i in 1:nrow(chrom_islands)) {
    if(i %% 100 == 0){
      print(c("i=",i,"completion =",i/nrow(chrom_islands)))
    }
    for (s in sample_names) {
      CpG_index <- chrom_islands$chromStart[i] <= chrom_df$start &
        chrom_df$end <= chrom_islands$chromEnd[i]
      
      mean_meth <-
        mean(x = chrom_df[CpG_index, paste("fRead", s, sep = "_")], na.rm = TRUE)
      chrom_islands[i, paste("mean_meth", s, sep = "_")] <-
        mean_meth
      
    }
  }
  saveRDS(object = chrom_islands,
          file = file.path(input_folder,paste(chr,"islands","rds",sep = ".")))
}


ggplot2::ggplot(data = chrom_islands,
                mapping = ggplot2::aes(x = chrom,y = mean_meth_HepG2 - mean_meth_K562))+
  ggplot2::geom_boxplot()+ggplot2::geom_point(position = "jitter")



i_bed <- data.frame()
for (c in CHR_NAMES) {
  #deug c <- CHR_NAMES[1]
  print(c)
  chrom_islands <- island_bed[island_bed$chrom == c, ]
  chrom_df <- df[df$chr == c, ]
  #debug
  chrom_islands <- chrom_islands[1:10, ]
  
  registerDoParallel(n_cores)  # use multicore
  processed_results_summery <-
    foreach (i = 1:1, .combine = rbind) %dopar% {
      # for(i in 1:nrow(chrom_islands)){
      result <-
        data.frame(numeric(nrow(chrom_islands)), numeric(nrow(chrom_islands)))
      colnames(result) <- paste("mean_meth", sample_names, sep =  "_")
      for (s in sample_names) {
        CpG_index <- chrom_islands$chrom[i] == chrom_df$chr &
          chrom_islands$chromStart[i] <= chrom_df$start &
          chrom_df$end <= chrom_islands$chromEnd[i]
        
        mean_meth <-
          mean(x = chrom_df[CpG_index, paste("fRead", s, sep = "_")], na.rm = TRUE)
        result[i, paste("mean_meth", s, sep = "_")] <- mean_meth
      }
      return(result)
    }
  stopImplicitCluster()
  i_bed <- rbind(i_bed, processed_results_summery)
}


# registerDoParallel(n_cores)  # use multicore
# processed_results_summery <- foreach (c = 1:length(CHR_NAMES), .combine = rbind) %dopar% {
#
#   chrom_islands <- island_bed[island_bed$chrom == c,]
#   chrom_df <- df[df$chr == c,]
#   #debug
#   chrom_islands <- chrom_islands[1:10,]
#
#
#
#       result <- data.frame(row.names = sample_names)
#       for(s in sample_names){
#
#         CpG_index <- chrom_islands$chrom[i] == chrom_df$chr &
#           chrom_islands$chromStart[i] <= chrom_df$start &
#           chrom_df$end <= chrom_islands$chromEnd[i]
#
#         mean_meth <- mean(x = chrom_df[CpG_index,paste("fRead",s,sep = "_")],na.rm = TRUE)
#         result[s,paste("mean_meth",s,sep = "_")] <- mean_meth
#       }
#       return(result)
#     }
#   i_bed <- rbind(i_bed,processed_results_summery)
# }



#debug chromatin_state_names <- chromatin_state_names[1]


for (i in 1:nrow(island_bed)) {
  if (i %% 1000 == 1) {
    print(i)
  }
  
  for (s in sample_names) {
    CpG_index <- island_bed$chrom[i] == df$chr &
      island_bed$chromStart[i] <= df$start &
      df$end <= island_bed$chromEnd[i]
    
    mean_meth <-
      mean(x = df[CpG_index, paste("fRead", s, sep = "_")], na.rm = TRUE)
    island_bed[i, paste("mean_meth", s, sep = "_")] <- mean_meth
  }
  
}
