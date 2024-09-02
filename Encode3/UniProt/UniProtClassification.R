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

# Combine all GO term entries into one big character vector
all_go_terms <- paste(UniProt$`Gene.Ontology.(molecular.function)`, collapse = "; ")


# Split this big string into individual GO terms
go_terms_list <- unlist(strsplit(all_go_terms, "; "))

# Remove duplicates
unique_go_terms <- unique(go_terms_list)
# go terms have "DNA" inside
dna_binding <- unique_go_terms[grepl("DNA[- ]?binding", unique_go_terms, ignore.case = TRUE)]

dna_binding <- c("cis-regulatory region sequence-specific DNA binding [GO:0000987]",                                     
"DNA-binding transcription factor activity [GO:0003700]" ,                                              
"DNA-binding transcription factor activity, RNA polymerase II-specific [GO:0000981]",                   
 "DNA-binding transcription repressor activity [GO:0001217]"  ,                                          
"DNA-binding transcription repressor activity, RNA polymerase II-specific [GO:0001227]",                
 "RNA polymerase II cis-regulatory region sequence-specific DNA binding [GO:0000978]",                   
 "sequence-specific double-stranded DNA binding [GO:1990837]",                                           
 "RNA polymerase II-specific DNA-binding transcription factor binding [GO:0061629]"   ,                  
 "DNA binding [GO:0003677]"   ,                                                                          
 "DNA-binding transcription activator activity, RNA polymerase II-specific [GO:0001228]" ,               
 "sequence-specific DNA binding [GO:0043565]" ,                                                          
 "DNA-binding transcription factor binding [GO:0140297]" ,                                               
 "RNA polymerase II transcription regulatory region sequence-specific DNA binding [GO:0000977]"  ,       
 "damaged DNA binding [GO:0003684]"  ,                                                                   
 "double-stranded DNA binding [GO:0003690]" ,                                                            
 "RNA polymerase II core promoter sequence-specific DNA binding [GO:0000979]" ,                          
 "core promoter sequence-specific DNA binding [GO:0001046]"                 ,                            
 "RNA polymerase I cis-regulatory region sequence-specific DNA binding [GO:0001165]"  ,                  
 "RNA polymerase I core promoter sequence-specific DNA binding [GO:0001164]" ,                           
 "DNA-binding transcription activator activity [GO:0001216]"  ,                                          
 "double-stranded methylated DNA binding [GO:0010385]" ,                                                 
 "hemi-methylated DNA-binding [GO:0044729]"  ,                                                           
 "DNA binding domain binding [GO:0050692]"    ,                             
 "intronic transcription regulatory region sequence-specific DNA binding [GO:0001161]"  ,                
 "single-stranded DNA binding [GO:0003697]"             ,                                                
 "sequence-specific single stranded DNA binding [GO:0098847]"    ,                                       
 "RNA polymerase II intronic transcription regulatory region sequence-specific DNA binding [GO:0001162]",
"nucleosomal DNA binding [GO:0031492]"               ,                                                  
 "RNA polymerase III type 3 promoter sequence-specific DNA binding [GO:0001006]")  


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

# UniProt[grepl(pattern = "Ago1",x = UniProt$From,ignore.case = T),]
# 
# UniProt$DNA_Binding_GO[grepl(pattern = "Ago1",x = UniProt$From,ignore.case = T)]
# UniProt$DNA_Binding_GO[grepl(pattern = "PRPF4",x = UniProt$From,ignore.case = T)]

# Load necessary libraries
library(dplyr)
library(stringr)

# File paths
uniprot_file <- "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/UniProt/idmapping_reviewed_true_AND_model_organ_2024_06_16.tsv"
motif_file <- "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/H12CORE_motifs.tsv"
output_file <- "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/UniProt/idmapping_with_H12CORE_motifs.tsv"

# Load the TSV files into data frames
uniprot_df <- read.delim(uniprot_file, stringsAsFactors = FALSE)
motif_df <- read.delim(motif_file, stringsAsFactors = FALSE)

# Preprocess the Motif column to extract the first part of the motif name (before the first dot)
motif_df$Motif_Short <- sapply(strsplit(motif_df$Motif, "\\."), `[`, 1)

# Function to find all matching motifs
find_motif_matches <- function(row, motifs, full_motifs) {
  from_matches <- full_motifs[motifs %in% row["From"]]
  entry_name_matches <- full_motifs[motifs %in% str_replace(row["Entry.Name"], "_HUMAN", "")]
  gene_names <- strsplit(row["Gene.Names"], " ")[[1]]
  gene_name_matches <- full_motifs[motifs %in% gene_names]
  
  # Combine all matches and remove duplicates
  all_matches <- unique(c(from_matches, entry_name_matches, gene_name_matches))
  
  # Return as a space-separated string
  if (length(all_matches) > 0) {
    return(paste(all_matches, collapse = " "))
  } else {
    return(NA)
  }
}

# Apply the function to each row in the UniProt DataFrame and add the results to a new column
uniprot_df$H12CORE_motif <- apply(uniprot_df, 1, function(row) find_motif_matches(row, motif_df$Motif_Short, motif_df$Motif))

# Save the updated DataFrame to a new TSV file
write.table(uniprot_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

print(sum(!is.na(uniprot_df$H12CORE_motif)))

cat("Output saved to", output_file, "\n")