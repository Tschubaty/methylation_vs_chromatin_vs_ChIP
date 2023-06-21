#################################################################
##  Chip-seq
##
##  input:  Encode3\
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
input_folder <- "download"
output_folder <- "processed"
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
# load meta files
meta <- openxlsx::read.xlsx(
  xlsxFile = encode_meta_file_name,colNames = FALSE
)
# remove duplicates
meta <- meta[-which(duplicated(meta)), ]

colnames(meta) <- c("file_type", "peak_type", "type" ,"pipeline" , "experiment","Target" , "Assay" ,"ChIP-seq" ,"Biosample" ,"cell_line", "status" ,"file_code")

meta <- meta[,c("peak_type", "Target" , "cell_line", "file_code")]

meta <- meta[order(meta$Target),]

###############################################
df_cpg <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/WGBS_vs_ChIP_ENCODE3/Encode3/WGBS/merged/all.merged_CpG_meth.RDS")

meta$n_peaks <- NA
meta$covergage <- NA

for(f in 1:nrow(meta)){
  
  file_name <- file.path(input_folder,meta$cell_line[f],paste(meta$file_code[f],"bed",sep = "."))
  print(file_name)
  histon_file <-read.table(file = file_name,col.names = c("chr","start","end" ,"name","score","strand","signalValue" ,"pValue", "qValue",  "peak"))
  df_histon <- histon_file[histon_file$chr %in% CHR_NAMES,]
  colName <- paste(meta$cell_line[f],meta$Target[f],sep = "_")
  meta$n_peaks[f] <- nrow(df_histon)
  meta$covergage[f] <- sum(df_histon$end - df_histon$start)
}















#################################################################################
j <- 1
s <- sample_names[j]
files <- meta$file_code[meta$cell_line == s]





i <- 1
file_name <- file.path(input_folder,s,paste(files[i],"bed",sep = "."))
colName <- paste(s,meta$Target[meta$file_code == files[i]],sep = "_")
df_cpg[,colName ] <- NA

histon_file <- read.table(file = file_name,col.names = c("chr","start","end" ,"name","score","strand","signalValue" ,"pValue", "qValue",  "peak")) 

for(chr in CHR_NAMES){
  print(chr)
  df_cpg_chr <- df_cpg[peak$chr == df_cpg$chr,] 
  histon_file_chr <- histon_file[histon_file$chr == chr,]
  for(r in 1:nrow(histon_file_chr)){
    peak <- histon_file_chr[r,]
    df_cpg_chr[peak$start <= df_cpg_chr$start & df_cpg_chr$end <= peak$end,colName] <- peak$score
  }
  df_cpg[peak$chr == df_cpg$chr,] <- df_cpg_chr
}


# registerDoParallel(n_cores)  # use multicore
# processed_results_summery <-
#   foreach (f = 1:length(protein_files), .combine = rbind) %dopar% {
#     #for (f in protein_files){
#     
#     #debug f <- which(grepl(x = protein_files,pattern = "HNRNPL"))[1]
#     # f <- 2
#     print(f)
#     file_name <- protein_files[f]
#     protein_name <-
#       sub(pattern = ".ChiP_vs_WGBS.RDS",
#           replacement = "",
#           x = file_name)
#     print(protein_name)
#     df <- readRDS(file.path(input_folder_protein , file_name))
#     #df$chromatin_state <- NA
#     df$chromatin_state_HepG2 <- ""
#     df$chromatin_state_K562 <- ""
#     
#     print(nrow(df))
#     for(r in 1:nrow(df)){
#       #r <- 1
#       if(r%% 1000 == 0){
#         #print(r)
#       }
#       
#       s <- "HepG2"
#       
#       temp <- as.character(chromatine_files[[s]][chromatine_files[[s]]$chr == df$chr[r] & 
#                                                    chromatine_files[[s]]$start <= df$merge_peak[r] & 
#                                                    df$merge_peak[r] < chromatine_files[[s]]$end ,"name"])
#       if(!length(temp) == 0){df$chromatin_state_HepG2[r] <- temp}
#       
#       # r <- 1579
#       # chromatine_files[[s]][31436,]
#       
#       s <- "K562"
#       temp <- as.character(chromatine_files[[s]][chromatine_files[[s]]$chr == df$chr[r] & 
#                                                    chromatine_files[[s]]$start <= df$merge_peak[r] & 
#                                                    df$merge_peak[r] < chromatine_files[[s]]$end ,"name"])
#       
#       if(!length(temp) == 0){df$chromatin_state_K562[r] <- temp}
#       
#       # if(length(chromatin_state) <= 0){
#       #   df[r,"chromatin_state"] <- NA
#       # }else{
#       #   df[r,"chromatin_state"] <- chromatin_state
#       # }
#       
#     }
#     
#     saveRDS(object = df,file = file.path(output_folder,sub(pattern = ".RDS",replacement = "_vs_chromatin.RDS",x =  file_name)))
#     
#     # p_island <- ggplot2::ggplot(data = ADNP.ChiP_vs_WGBS_vs_chromatin ,mapping = ggplot2::aes(x = chromatin_state,fill = fRead_K562 < 0.25))+
#     #   ggplot2::geom_bar(position="dodge") + ggplot2::coord_flip()
#     # p_ChiP <- ggplot2::ggplot(data = ADNP.ChiP_vs_WGBS_vs_chromatin ,mapping = ggplot2::aes(x = chromatin_state,fill = ChiP_state))+
#     #   ggplot2::geom_bar(group = "ChiP_state",position="dodge") + ggplot2::coord_flip()
#     # 
#     # plot_list <-  list(p_island ,p_ChiP)
#     # lay <- rbind(
#     #   c(1),
#     #   c(2)
#     # )
#     # 
#     # final_plot <-
#     #   gridExtra::grid.arrange(
#     #     grobs = plot_list,
#     #     shared_legend,
#     #     top = protein_name,
#     #     layout_matrix = lay
#     #   )
#     return(nrow(df))
#   }
