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
#encode_meta_file_name <- "Histone-data.xlsx"
input_folder_WGBS_HISTONE <- "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/Histones/processed"
output_folder <- "out"
input_folder_chromatin <- "by_chr"
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

# input histones
meth_files <- paste(CHR_NAMES,"WGBS.histones.rds",sep = ".")

chromatin_files <- paste(CHR_NAMES,"chromatin.bed",sep=".")

#debug

# registerDoParallel(n_cores)  # use multicore
# processed_results_summery <-
#   foreach (chr_n = 1:length(CHR_NAMES ), .combine = rbind) %dopar% {
    chr_n <- 21
    
    file_name <- meth_files[chr_n]
    print(file_name)
    df <- readRDS(file.path(input_folder_WGBS_HISTONE, file_name))
    #df$chromatin_state <- NA
    df$chromatin_state_K562 <- NA
    df$chromatin_state_HepG2 <- NA

    # get chrn specific locations
    df_chromatin <-
      read.table(
        file = file.path(input_folder_chromatin,chromatin_files[chr_n]),
        header = FALSE,
        col.names = c(
          "chr" ,
          "start",
          "end",
          "chromatin_state_HepG2",
          "chromatin_state_K562"
        )
      )
    
    
    
 
    
    # debug
    #df <- rbind(head(df,20000),tail(df,20000))
    #df <- df[200000:300000,]
    print(nrow(df))
    i <- 1

    ######################pb = txtProgressBar(min = 0, max = nrow(df), initial = 1) 
    start_time <- Sys.time()
    
      chromatin_state_K562 <- df$chromatin_state_K562 
      chromatin_state_HepG2 <- df$chromatin_state_HepG2 
      start <- df$start
      end <- df$end
      
      state_HepG2 <- df_chromatin$chromatin_state_HepG2
      state_K562 <- df_chromatin$chromatin_state_HepG2
      state_start <- df_chromatin$start
      state_end <- df_chromatin$end
      
      for(r in 1:length(start)){
        #  for(r in 200000:210000){
        ######################setTxtProgressBar(pb,i)
        #r <- 1
        if(r%% 10000 == 0){
          print(paste("CpG row number:",r,"chromatin row number" ,i,"percent = ",r/length(start)*100,sep=" "))
        }
        ########################### K562 ##########################
        # check if need to go to next segment
        while (end[r] > state_end[i] & i <= length(state_end)) {
          i <- i+1
        }
        
        # if there is annotation in segement
        if(state_start[i] <= start[r]){
          #DO ANNOATTE
          chromatin_state_K562[r] <- state_K562[i]
          chromatin_state_HepG2[r] <- state_HepG2[i]
        }
        
      }
      
      df$chromatin_state_K562 <- chromatin_state_K562
      df$chromatin_state_HepG2 <- chromatin_state_HepG2

    end_time <- Sys.time()
    print(end_time - start_time) 
    
    ######################close(pb)
    saveRDS(object = df,file = file.path(output_folder,sub(pattern = ".rds",replacement = ".chromatin.RDS",x =  file_name)))
    
}    
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
  #   return(nrow(df))
  # }    
