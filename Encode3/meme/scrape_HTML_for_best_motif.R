# Documentation -----------------------------------------------------------
##
##
##  input: 
##
##
##  output: 
##
##
##  v_01 07.07.2024
##  Author: Daniel Batyrev (HUJI 777634015)
##
# Set up Work Environment --------------------------------------------------
#Clear R working environment
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

library(rvest)

input_folder <- file.path("meme-chip_results","pooled_AGO1","centrimo_out")

# Load the HTML file
html_content <- read_html(file.path(input_folder,"centrimo.html"))

# Extract the table
table <- html_content %>% html_node("#motifs") %>% html_table()

# Extract concentration values
concentration_values <- as.numeric(table$concentration)

# Find the best concentration value
best_concentration <- max(concentration_values, na.rm = TRUE)

print(paste("The best concentration value is:", best_concentration))
