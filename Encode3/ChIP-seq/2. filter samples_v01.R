#################################################################
##  Chip-seq
##
##  input:  Encode3\ChIP-seq\1. download from Encode3
##  output: Encode3\ChIP-seq\2. filter samples
##  v 01 - 14.09.2021
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
encode_meta_file_name <- "metadata.tsv"
input_folder <- "1. download from Encode3"
output_folder <- "2. filter samples"
################################################ INPUT ###################################################


################################################ CODE ###################################################
#### take meta files from #Sample*Sept2021 and copy files to output_folder   #######
for (s in sample_names) {
  input_folder_data <- paste(s, "Sept2021", sep = "")
  print(input_folder_data)
  
  meta <-
    read.csv(
      file = file.path(input_folder,input_folder_data, encode_meta_file_name),
      header = TRUE,
      sep = "\t"
    )
  
  #meta$Biosample.term.name
  #meta$File.assembly
  #conservative_ChiP_sept21
  
  uniform_columns <- c()
  for (cn in colnames(meta)) {
    print(cn)
    different <- length(unique(meta[, cn]))
    print(different)
    if (different <= 1) {
      uniform_columns <- c(uniform_columns, cn)
    }
  }
  
  my_meta <- meta[, setdiff(colnames(meta), uniform_columns)]
  
  my_meta <-
    meta[meta$File.analysis.status == "released", c(
      "Experiment.target",
      "Biosample.term.name",
      "File.assembly",
      "File.accession",
      "Experiment.accession",
      "Biosample.genetic.modifications.methods",
      "Biosample.genetic.modifications.categories",
      "Biosample.genetic.modifications.targets",
      "Biosample.genetic.modifications.site.coordinates",
      "Library.lysis.method",
      "Library.crosslinking.method",
      "Library.fragmentation.method",
      "Library.size.range"
    )]
  
  for (r in 1:nrow(my_meta)) {
    f <- my_meta$File.accession[r]
    file_name  <- paste(f, "bed", sep = ".")
    out_file_name <- paste(my_meta$Biosample.term.name[r],
                           my_meta$Experiment.target[r],
                           file_name,
                           sep = "_")
    file.copy(
      from = file.path(input_folder,input_folder_data, file_name),
      to = file.path(output_folder, out_file_name)
    )
  }
  
  saveRDS(object = my_meta,
          file = file.path(
            output_folder,
            paste(s, "meta_summery", "RDS", sep = ".")
          ))
}