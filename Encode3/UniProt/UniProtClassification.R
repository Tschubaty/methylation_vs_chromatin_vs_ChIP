#################################################################
##  Uniprot classification
##  
##  input: 
##  ##
##  output:      
##  v01 - 21.06.2024
##  Author: Daniel Batyrev 777634015
#################################################################

#Clear R working environment 
#rm(list=ls())

this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"

setwd(this.dir)

file_path <- "idmapping_reviewed_true_AND_model_organ_2024_06_16.xlsx" 
UniProt <- openxlsx::read.xlsx(file_path)

# Combine all GO term entries into one big character vector
all_go_terms <- paste(UniProt$`Gene.Ontology.(GO)`, collapse = "; ")

# Split this big string into individual GO terms
go_terms_list <- unlist(strsplit(all_go_terms, "; "))

# Remove duplicates
unique_go_terms <- unique(go_terms_list)
# go terms have "DNA" inside
dna_binding <- unique_go_terms[grepl("DNA[- ]?binding", unique_go_terms, ignore.case = TRUE)]

# Function to check for any DNA-binding GO term in each protein's GO annotations
go_dna_binding <- function(go_string) {
  # Split the individual protein's GO string into a list of terms
  protein_go_terms <- unlist(strsplit(go_string, "; "))
  #print(protein_go_terms)
  # Check if any of these terms match the dna_binding terms found earlier
  paste(intersect(dna_binding, protein_go_terms),collapse = ";")
}

number_go_dna_binding <- function(go_string) {
  # Split the individual protein's GO string into a list of terms
  protein_go_terms <- unlist(strsplit(go_string, "; "))
  #print(protein_go_terms)
  # Check if any of these terms match the dna_binding terms found earlier
  length(intersect(dna_binding, protein_go_terms))
}

# Apply this function to each row in the UniProt DataFrame
UniProt$DNA_Binding_GO <- sapply(UniProt$`Gene.Ontology.(GO)`, go_dna_binding)
UniProt$number_DNA_Binding_GO <- sapply(UniProt$`Gene.Ontology.(GO)`, number_go_dna_binding)

UniProt[grepl(pattern = "Ago1",x = UniProt$From,ignore.case = T),]

UniProt$DNA_Binding_GO[grepl(pattern = "Ago1",x = UniProt$From,ignore.case = T)]
UniProt$DNA_Binding_GO[grepl(pattern = "PRPF4",x = UniProt$From,ignore.case = T)]

