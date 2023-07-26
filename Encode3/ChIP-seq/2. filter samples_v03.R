#################################################################
##  Chip-seq
##
##  input:  Encode3\ChIP-seq\1. download from Encode3
##  output: Encode3\ChIP-seq\2. filter samples
##  v 03 - 26.07.2023
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
sample_names <- c("K562", "HepG2", "A549", "GM12878")
encode_meta_file_name <- "metadata.tsv"
input_folder <- "1. download from Encode3"
output_folder <- "2. filter samples"
picuture_file_extension = "png"
################################################ INPUT ###################################################


################################################ CODE ###################################################
#### take meta files from #Sample*Sept2021 and copy files to output_folder   #######


my_meta <- data.frame()

for (s in sample_names) {
  #s <- sample_names[1]
  sample_folder <- paste(s, "july2023", sep = "")
  meta <-
    read.csv(
      file = file.path(input_folder, sample_folder , encode_meta_file_name),
      header = TRUE,
      sep = "\t",
    )
  meta$folder <- sample_folder
  
  my_meta <- rbind(my_meta, meta)
}

# dete uninformative columsn
uniform_columns <- c()
for (cn in colnames(my_meta)) {
  #print(cn)
  different <- length(unique(my_meta[, cn]))
  #print(different)
  if (different <= 1) {
    uniform_columns <- c(uniform_columns, cn)
  }
}
my_meta <- my_meta[, setdiff(colnames(my_meta), uniform_columns)]

analyis_year <-
  substr(gsub(
    pattern = "s3://encode-public/",
    replacement = "",
    x = my_meta$s3_uri[my_meta$File.analysis.status == "released"]
  ),
  1,
  4)

my_meta <-
  my_meta[my_meta$File.analysis.status == "released", c(
    "Experiment.target",
    "Biosample.term.name",
    "File.accession",
    "Experiment.accession",
    "Experiment.date.released",
    "Biological.replicate.s.",
    "Technical.replicate.s."
  )]

# delete old analysis
my_meta$analyis.year <- analyis_year
# 
# for(Experiment.accession in unique(my_meta$Experiment.accession)){
#   # my_meta$Experiment.accession[1] <- Experiment.accession
#   dup_rows <- my_meta[Experiment.accession == my_meta$Experiment.accession,]
# }
#
df <- my_meta[,c("Experiment.target","Biosample.term.name","Experiment.accession", "Experiment.date.released")]
# repeated_analysis <- duplicated(df) | duplicated(df, fromLast = TRUE)
########################## last anayisy is news one checked by code above
my_meta <- my_meta[!duplicated(df, fromLast = TRUE),]

# get sample coverage
df_overview <- data.frame(target = unique(my_meta$Experiment.target))
for(s in sample_names){
  df_overview[,s] <- df_overview$target %in% my_meta$Experiment.target[my_meta$Biosample.term.name == s]
}
df_overview$sum <- apply(X = df_overview[,-1],MARGIN = 1,FUN = sum)

# only include targets with coverage in multiple samples 
my_meta <- my_meta[my_meta$Experiment.target %in% df_overview$target[df_overview$sum >= 2],]

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


for (r in 1:nrow(my_meta)) {
  f <- my_meta$File.accession[r]
  file_name  <- paste(f, "bed", sep = ".")
  target <-
    gsub(
      pattern = "-human",
      replacement = "",
      x = my_meta$Experiment.target[r],
      ignore.case = TRUE
    )
  out_file_name <- paste(my_meta$Biosample.term.name[r],
                         my_meta$Experiment.target[r],
                         file_name,
                         sep = "_")
  sample_folder <-
    paste(my_meta$Biosample.term.name[r], "july2023", sep = "")
  dir.create(file.path(output_folder, target), showWarnings = FALSE)
  file.copy(
    from = file.path(input_folder, sample_folder, file_name),
    to = file.path(output_folder, target, out_file_name)
  )
}

saveRDS(object = my_meta,
        file = file.path(output_folder,
                         paste("meta_summery", "RDS", sep = ".")))



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
