#################################################################
##  Chip-seq
##
##  input:  Encode3\ChIP-seq\1. download from Encode3
##  output: Encode3\ChIP-seq\2. filter samples
##  v 02 - 24.07.2023
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

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}


################################################ INPUT ###################################################
sample_names <- c("K562", "HepG2","A549","GM12878")
#encode_meta_file_name <- "metadata.tsv"
input_folder <- "1. download from Encode3"
output_folder <- "2. filter samples"
picuture_file_extension = "png"
################################################ INPUT ###################################################


################################################ CODE ###################################################
#### take meta files from #Sample*Sept2021 and copy files to output_folder   #######

df_meta <- data.frame()

for (s in sample_names) {
  #s <- sample_names[1]
  meta_file_name <- paste("experiment_report",s,"conservative IDR thresholded peaks.tsv",sep = "_")
  print(meta_file_name)
  meta <-
    read.csv(
      file = file.path(input_folder, paste(s,"july2023",sep = ""),meta_file_name),
      header = TRUE,
      sep = "\t",
      skip = 1
    )
  
  df_meta <- rbind(df_meta,meta)
  
  
  
  
  # my_meta <- meta[, setdiff(colnames(meta), uniform_columns)]
  # 
  # print(head(meta$Files)[1])
  #

}

df_overview <- data.frame(target = unique(df_meta$Target.of.assay))
for(s in sample_names){
  df_overview[,s] <- df_overview$target %in% df_meta$Target.of.assay[df_meta$Biosample.term.name == s]
}
df_overview$sum <- apply(X = df_overview[,-1],MARGIN = 1,FUN = sum)

# ven diagram
ven_values <-  c(
  # single
  "K562" = sum(df_overview$K562 & df_overview$sum == 1),
  "HepG2" = sum(df_overview$HepG2 & df_overview$sum == 1),
  "A549" = sum(df_overview$A549 & df_overview$sum == 1),
  "GM12878" =  sum(df_overview$GM12878 & df_overview$sum == 1),
  # two
  "K562&HepG2" = sum(df_overview$K562 & df_overview$HepG2 & df_overview$sum == 2),
  "K562&A549" = sum(df_overview$K562 & df_overview$A549 & df_overview$sum == 2),
  "K562&GM12878" = sum(df_overview$K562 & df_overview$GM12878 & df_overview$sum == 2),
  "HepG2&A549" = sum(df_overview$HepG2 & df_overview$A549 & df_overview$sum == 2),
  "HepG2&GM12878" = sum(df_overview$HepG2 & df_overview$GM12878 & df_overview$sum == 2),
  "A549&GM12878" = sum(df_overview$A549 & df_overview$GM12878 & df_overview$sum == 2),
  # three
  "HepG2&GM12878&A549" = sum(!df_overview$K562 & df_overview$sum == 3),
  "K562&GM12878&A549" = sum(!df_overview$HepG2 & df_overview$sum == 3),
  "K562&HepG2&GM12878" = sum(!df_overview$A549 & df_overview$sum == 3),
  "K562&HepG2&A549" =  sum(!df_overview$GM12878 & df_overview$sum == 3),
  # all
  "K562&HepG2&A549&GM12878" = sum(df_overview$sum == 4)
)


group.colors <-
  c(K562 = "#F8766D",
    HepG2 = "#00BFC4",
    A549 = "#F0A000",
    GM12878 = "#50BF6D")

plot_ven <- plot(
  eulerr::euler(ven_values),
  fills = list(fill = group.colors),
  legend = list(side = "right"),
  quantities = list(cex = 2),
  labels = c()
)


ggplot2::ggsave(
  filename = paste(
    "Venn_Diagramm_ChiP_data",
    picuture_file_extension,
    sep = "."
  ),
  plot = plot_ven,
  device = picuture_file_extension,
  path =  output_folder,
  width = 1920,
  height = 1080,
  units = "px"
)


for (r in 1:nrow(df_overview)) {
  # check if corresponding chip data is available
  # r <- 4
  target <- df_overview$target[r]
  print(target)
  if (df_overview$sum[r] >= 2) {
    # create folder
    dir.create(file.path(output_folder, target), showWarnings = FALSE)
    for (s in sample_names[as.logical(df_overview[r,c(-1,-6)])]) {
      #s <- sample_names[3]
      # make search pattern
      print(s)
      search_patter <-
        gsub(pattern = "/|files",
             replacement = "",
             x = df_meta$Files[df_meta$Target.of.assay == target &
                                 df_meta$Biosample.term.name == s])
      search_patter <- unlist(strsplit(x = search_patter, split = ","))
      #search_patter <- paste(search_patter ,"bed",collapse = ".")
      search_patter <-
        paste(sapply(
          search_patter,
          FUN = function(x)
            paste(x, "bed", sep = ".")
        ), collapse = "|")
      
      found_files <-
        list.files(path = file.path(input_folder, paste(s, "july2023", sep = "")),
                   pattern = search_patter)
      
      #stopifnot( exprs = length(found_files) > 0, local = TRUE)
      for (f in found_files) {
        #f <- f
        out_file_name <- paste(s, target, f, sep = "_")
        file.copy(
          from = file.path(input_folder, paste(s, "july2023", sep = ""), f),
          to = file.path(output_folder, target , out_file_name)
        )
      }
    }
  }
}

targets_with_data <- list.dirs(path = output_folder,full.names = FALSE)
df_overview2 <- data.frame()

for (t in targets_with_data[-1]) {
  # t <- targets_with_data[1]
  #print(t)
  files <- list.files(path = file.path(output_folder, t))
  df_overview2 <- rbind(df_overview2 ,
                        data.frame(
                          target = t,
                          cell_line = gsub(paste("_",t,"-human",".*",sep = ""), "", files),
                          file_name = files)
                        )
}

#   saveRDS(object = my_meta,
#           file = file.path(
#             output_folder,
#             paste(s, "meta_summery", "RDS", sep = ".")
#           ))
# }