# Documentation -----------------------------------------------------------
##
##
##  input: 
##
##
##  output: 
##
##
##  v_01 21.07.2024
##  Author: Daniel Batyrev (HUJI 777634015)
##
# Set up Work Environment --------------------------------------------------
# Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/HumanEvo/HumanEvo/"
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
#################################################################################

df <- read.csv(file = "all_motifs_groups_all_with_CG.csv")

df.data <- df[,!names(df) %in% c("single.Motif.source","FIMO.gff","SpaMo","FIMO.Motif.source.file")]

# Create a new column with the extracted single Motif ID
df.data$database_Motif_ID <- ifelse(grepl("_HUMAN", df.data$single.Motif.ID, ignore.case = TRUE),
                                      sub("_HUMAN.*", "", df.data$single.Motif.ID),
                                      "")

# View the updated data frame
print(df.data[df.data$database_Motif_ID == df.data$Protein,])

df.data.exactMatch <- df.data[df.data$database_Motif_ID == df.data$Protein,]

# Load necessary library
library(ggplot2)

# Assuming df.data.exactMatch is already loaded
# Create a bar plot for FIMO.Group.Name
ggplot(df.data.exactMatch, aes(x = FIMO.Group.Name)) +
  geom_bar() +
  labs(title = "Histogram of FIMO Group Names",
       x = "FIMO Group Name",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for readability

