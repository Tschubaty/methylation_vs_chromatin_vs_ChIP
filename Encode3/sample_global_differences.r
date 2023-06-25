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

# library(foreach)
# library(doParallel)
library(ggplot2)
# library(hrbrthemes)
n_cores <- detectCores() - 1
########################## code ############################

input_folder_data <- "WGBS/complete"
plot_folder <- "Chromatin/plot"
output_folder <- "Chromatin/meta"
  

##########################################################################################################
#################################### plot and meta ######################################################
##########################################################################################################

chromatin_state_names <-
  c(
    "1_TssA",
    "2_TssFlnk",
    "3_TssFlnkU",
    "4_TssFlnkD",
    "5_Tx",
    "6_TxWk",
    "7_EnhG1",
    "8_EnhG2",
    "9_EnhA1",
    "10_EnhA2",
    "11_EnhWk",
    "12_ZNF_Rpts",
    "13_Het",
    "14_TssBiv",
    "15_EnhBiv",
    "16_ReprPC",
    "17_ReprPCWk",
    "18_Quies"
  )

group.colors <-
  c(HepG2 = "#F8766D",
    K562 = "#00BFC4",
    both = "#e4d00a")


############################ merge all data ###################################

for (chr in CHR_NAMES) {
  #deug chr <- CHR_NAMES[1]
  print(chr)
  
  df_chr <- readRDS(file = file.path(
      input_folder_data,
      paste(chr, "WGBS.histones.chromatin.islands.RDS", sep = ".")
    ))
  
  if(chr == CHR_NAMES[1]){
    df <- df_chr
  }else{
    df <- rbind(df,df_chr)
  }
  
}

saveRDS(object = df,file = file.path("WGBS","all.WGBS.histones.chromatin.islands.RDS", sep = "."))
######################################################################################
meta <- data.frame()
#colnames(meta) <- c("chromatin_state","sample","property","value")


r <- 1
for(state in 1:length(chromatin_state_names)){
  print(state)
  for(s in sample_names){
    targets <- colnames(df)[grep(pattern = paste(s,"_",sep = ""),x = colnames(df),ignore.case = TRUE)]
    print(targets)
    sample_chromatin_state_index <- df[,paste("chromatin_state",s,sep = "_")] == state
    state_cpg <- sum(sample_chromatin_state_index ,na.rm = TRUE)
    for(target in targets){
      print(target)
      meta[r,"chromatin_state"] <- state
      meta[r,"sample"] <- s
      meta[r,"property"] <- gsub(pattern = paste(s,"_",sep = ""),replacement = "",x = target)
      meta[r,"state_cpg"] <- state_cpg
      meta[r,"value"] <- sum(!is.na(df[sample_chromatin_state_index,target]))
      meta[r,"target_CpGs"] <- sum(!is.na(df[,target]))
      r <- r+1
    }
  }
}

# calc enrichemnt
max_cpg <- nrow(df)
meta$enrichment <- meta$value/(meta$target_CpGs * (meta$state_cpg/max_cpg))
meta$log_enrichment <- log(meta$value/(meta$target_CpGs * (meta$state_cpg/max_cpg)))
meta$normalized_enrichment <- meta$enrichment 


for(target in unique(meta$property)){
  print(target)
  for(s in sample_names){
    print(s)
    meta$normalized_enrichment[s == meta$sample & meta$property == target] <- meta$enrichment[s == meta$sample & meta$property == target]/max(meta$enrichment[s == meta$sample & meta$property == target])
  }
  
}

#save file  
saveRDS(object = meta,file = "meta.WGBS.histones.chromatin.islands.txt")


plot_list <- list()
for (s in sample_names) {
  plot_list[[s]] <- ggplot(
    data = meta[meta$sample == s, ],
    mapping = aes(x = property , chromatin_state, fill = normalized_enrichment)
  ) +
    geom_tile() +
    ggtitle(s) +
    scale_fill_gradient(low = "white", high = "blue") +
    theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
}

final_plot <-
  gridExtra::grid.arrange(
    grobs = plot_list,
    shared_legend,
    top = "Figure 2C",
    ncol = 2
  )

ggplot2::ggsave(filename = file.path(plot_folder,"chromatin_state_vs_histones.png"),plot = final_plot,width = 10,height = 6)




#   
#   chrom_islands <- island_bed[island_bed$chrom == c, ]
#   chrom_df <- df[df$chr == c, ]
#   #debug
#   chrom_islands <- chrom_islands[1:10, ]
#   
#   registerDoParallel(n_cores)  # use multicore
#   processed_results_summery <-
#     foreach (i = 1:1, .combine = rbind) %dopar% {
#       # for(i in 1:nrow(chrom_islands)){
#       result <-
#         data.frame(numeric(nrow(chrom_islands)), numeric(nrow(chrom_islands)))
#       colnames(result) <- paste("mean_meth", sample_names, sep =  "_")
#       for (s in sample_names) {
#         CpG_index <- chrom_islands$chrom[i] == chrom_df$chr &
#           chrom_islands$chromStart[i] <= chrom_df$start &
#           chrom_df$end <= chrom_islands$chromEnd[i]
#         
#         mean_meth <-
#           mean(x = chrom_df[CpG_index, paste("fRead", s, sep = "_")], na.rm = TRUE)
#         result[i, paste("mean_meth", s, sep = "_")] <- mean_meth
#       }
#       return(result)
#     }
#   stopImplicitCluster()
#   i_bed <- rbind(i_bed, processed_results_summery)
# }


# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)

# Give extreme colors:
ggplot(data, aes(X, Y, fill= Z)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  theme_ipsum()

# Color Brewer palette
ggplot(data, aes(X, Y, fill= Z)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdPu") +
  theme_ipsum()

# Color Brewer palette
library(viridis)
ggplot(data, aes(X, Y, fill= Z)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme_ipsum()


file_names <- paste(CHR_NAMES,"WGBS.histones.chromatin.RDS",sep = ".")


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
