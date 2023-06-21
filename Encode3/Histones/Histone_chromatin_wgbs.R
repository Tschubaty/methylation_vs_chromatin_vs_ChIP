#################################################################
##  Chip-seq
##
##  input:  C:\Users\Daniel Batyrev\Documents\GitHub\WGBS_vs_ChIP_ENCODE3\Encode3\Histones\processed
##
##   E123 = K562
##    E118 = HepG2
##
##  v 01 - 15.06.2023
##  Author: Daniel Batyrev 777634015
#################################################################

#Clear R working environment
rm(list = ls())

this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"

setwd(this.dir)

detachAllPackages <- function() {
  basic.packages <-
    c("package:stats",
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
################################################ libs #################################################
library(foreach)
library(doParallel)
# library(ggplot2)
# library(gridExtra)
# library(ggVennDiagram)

################################################ INPUT ###################################################

sample_names <- c("K562", "HepG2")
encode_meta_file_name <- "Histone-data.xlsx"
input_folder <- "processed"
output_folder <- "processed"
input_folder_chromatin <- "C:/Users/Daniel Batyrev/Documents/GitHub/WGBS_vs_ChIP_ENCODE3/Encode3/Chromatine/states"
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

# set up parallel processing
n_cores <- detectCores() - 1
###################################################################################

################################################ INPUT ###################################################


################################################ CODE ###################################################
#Read Input 
bed_files <- paste(CHR_NAMES,"WGBS.histones.rds",sep = ".")

chromatine_files <- list()
for (s in sample_names) {
  file_name <- paste(s,"chromatine_18state.RDS",sep = ".")
  chromatine_files[[s]] <- readRDS(file = file.path(input_folder_chromatin,file_name))
}


#debug

registerDoParallel(n_cores)  # use multicore
processed_results_summery <-
  foreach (chr_n = 1:length(CHR_NAMES ), .combine = rbind) %dopar% {
    # chr_n <- 22
    
    
    file_name <- bed_files[chr_n]
    print(file_name)
    df <- readRDS(file.path(input_folder , file_name))
    #df$chromatin_state <- NA
    df$chromatin_state_K562 <- NA
    df$chromatin_state_HepG2 <- NA

    # get chrn specific locations
    K562_state <- chromatine_files[["K562"]][chromatine_files[["K562"]]$chr == CHR_NAMES[chr_n],] 
    HepG2_state <- chromatine_files[["HepG2"]][chromatine_files[["HepG2"]]$chr == CHR_NAMES[chr_n],] 
    
    print(nrow(df))
    # initiate k562 counter and hepg2 counter
    k <- 1
    h <- 1
    
    current_K562 <- NA
    current_HepG2 <- NA
    
    # debug
    #df <- rbind(head(df,20000),tail(df,20000))
    
    for(r in 1:nrow(df)){
      #r <- 1
      if(r%% 1000 == 0){
        print(paste("CpG row number:",r,"K562 k:",k,"Hepg2 h:" ,h,sep=" "))
      }
      ########################### K562 ##########################
      # check if need to go to next segment
      while (df$end[r] > K562_state$end[k] & k < nrow(K562_state)) {
        k <- k+1
        current_K562 <- K562_state$name[k]
      }
      
      # if there is annotation in segement
      if(K562_state$start[k] <= df$start[r]){
          #DO ANNOATTE
        df$chromatin_state_K562[r] <- current_K562
      }
      ######################### HepG2 ###########################
      
      # check if need to go to next segment
      while (df$end[r] > HepG2_state$end[h] & h < nrow(HepG2_state)) {
        h <- h+1
        current_HepG2 <- HepG2_state$name[h]
      }
      
      # if there is annotation in segement
      if(HepG2_state$start[h] <= df$start[r]){
        #DO ANNOATTE
        df$chromatin_state_HepG2[r] <- current_HepG2
      }
      
      
    }
    
    saveRDS(object = df,file = file.path(output_folder,sub(pattern = ".rds",replacement = ".chromatin.RDS",x =  file_name)))
    
    # p_island <- ggplot2::ggplot(data = ADNP.ChiP_vs_WGBS_vs_chromatin ,mapping = ggplot2::aes(x = chromatin_state,fill = fRead_K562 < 0.25))+
    #   ggplot2::geom_bar(position="dodge") + ggplot2::coord_flip()
    # p_ChiP <- ggplot2::ggplot(data = ADNP.ChiP_vs_WGBS_vs_chromatin ,mapping = ggplot2::aes(x = chromatin_state,fill = ChiP_state))+
    #   ggplot2::geom_bar(group = "ChiP_state",position="dodge") + ggplot2::coord_flip()
    # 
    # plot_list <-  list(p_island ,p_ChiP)
    # lay <- rbind(
    #   c(1),
    #   c(2)
    # )
    # 
    # final_plot <-
    #   gridExtra::grid.arrange(
    #     grobs = plot_list,
    #     shared_legend,
    #     top = protein_name,
    #     layout_matrix = lay
    #   )
    return(nrow(df))
  }    
