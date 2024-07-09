#################################################################
##  Chip-seq
##
##  input: Encode3\ChIP-seq\3. merge samples
##  output: Encode3\ChIP-seq\4. plot and summery
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

########################### input ##########################

#libs:
library(ggplot2)
library(gridExtra)
library(ggVennDiagram)

################################################ INPUT ###################################################
sample_names <- c("K562", "HepG2")
OVERLAP_CUTTOF <- 100
input_folder <- "3. merge samples"
output_folder <- "4. plot and summery"
################################################ INPUT ###################################################


################################################ CODE ###################################################
#Read Input 
rds_files <- list.files(path = input_folder, pattern = "*.RDS")
proteins <- gsub(pattern = ".RDS",
       replacement = "",
       x = unique(rds_files))

# meta_summery <- data.frame(matrix(nrow = 0,ncol = length(proteins)))
# colnames(meta_summery) <- proteins
# 
# for (target in proteins) {
#   # target <- proteins[which(proteins == "AGO1")]
#   meta_summeryp[target] <- length(unique(df$file))
#   }

for (target in proteins) {
  # target <- proteins[which(proteins == "AGO1")]
  print(target)
  
  df <- readRDS(file = file.path(input_folder,paste(target,"RDS",sep = ".")))

  p1 <- ggplot2::ggplot(data = df,mapping = ggplot2::aes(x = file,fill = state)) + 
    ggplot2::geom_bar()+coord_flip()
  
  
  p2 <- ggplot2::ggplot()
  
  
  p3 <- ggplot2::ggplot(data = df,mapping = ggplot2::aes(x = hits,fill = file)) + ggplot2::geom_bar()
  
  
  p4 <- ggplot2::ggplot(data = df,mapping = ggplot2::aes(x = sample,fill = state)) + ggplot2::geom_bar()

  
  plot_list = list(p1,p2,p3,p4)
  
  p_all <- gridExtra::grid.arrange(grobs= plot_list, nrow = 2,shared_legend,top = target)
  
  ggplot2::ggsave(filename = file.path(output_folder,
                                       paste(target,"png",sep = ".")),
                  width = 16, height = 8,
                  plot = p_all)
}


###########

meta_summary <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/ChIP-seq/2.filter_samples/meta_summery.RDS")

# Extract protein names without the "-human" suffix
protein_names <- gsub("-human", "", meta_summary$Experiment.target)

# Write to a text file
writeLines(text = protein_names,con = file.path(this.dir,"protein_names.txt"))

# Verify the output
cat(protein_names, sep = "\n")
