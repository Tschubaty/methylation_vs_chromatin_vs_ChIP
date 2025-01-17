#################################################################
##  Chip-seq
##
##  input:  
##
##   E123 = K562
##   E118 = HepG2
##
##  v 
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

################################################ INPUT ###################################################

# Suggested colors for the biosamples
group.colors <- c(
  HepG2    = "#F8766D",
  # A warm reddish color
  K562     = "#00BFC4",
  # A cool cyan color
  GM12878  = "#A3A500",
  # A yellow-green color, good contrast with the others
  A549     = "#E76BF3"   # A vibrant purple color
)


# set up parallel processing
n_cores <- detectCores() - 1

CHR_NAMES <- paste0("chr", c(1:22))

# Define chromatin state names
chromatin_state_names <- c(
  "1_TssA",       # Active TSS
  "2_TssFlnk",    # Flanking Active TSS
  "3_TssFlnkU",   # Flanking TSS Upstream
  "4_TssFlnkD",   # Flanking TSS Downstream
  "5_Tx",         # Strong Transcription
  "6_TxWk",       # Weak Transcription
  "7_EnhG1",      # Genic Enhancer 1
  "8_EnhG2",      # Genic Enhancer 2
  "9_EnhA1",      # Active Enhancer 1
  "10_EnhA2",     # Active Enhancer 2
  "11_EnhWk",     # Weak Enhancer
  "12_ZNF_Rpts",  # ZNF Genes & Repeats
  "13_Het",       # Heterochromatin
  "14_TssBiv",    # Bivalent TSS
  "15_EnhBiv",    # Bivalent Enhancer
  "16_ReprPC",    # Repressed Polycomb
  "17_ReprPCWk",  # Weak Repressed Polycomb
  "18_Quies"      # Quiescent/Low Activity
)

# Define corresponding hex color codes for each chromatin state
chromatin_state_colors <- c(
  "#FF0000",  # 1_TssA       -> Red
  "#FF4500",  # 2_TssFlnk    -> Orange-Red
  "#FF9900",  # 3_TssFlnkU   -> Orange
  "#FFCC00",  # 4_TssFlnkD   -> Yellow-Orange
  "#00CC00",  # 5_Tx         -> Green
  "#006400",  # 6_TxWk       -> Dark Green
  "#FFD700",  # 7_EnhG1      -> Gold
  "#FFD700",  # 8_EnhG2      -> Gold (same as EnhG1)
  "#FFFF00",  # 9_EnhA1      -> Yellow
  "#FFDD00",  # 10_EnhA2     -> Yellow-Orange
  "#FFEA73",  # 11_EnhWk     -> Light Yellow
  "#9370DB",  # 12_ZNF_Rpts  -> Purple
  "#C0C0C0",  # 13_Het       -> Light Gray
  "#FF4500",  # 14_TssBiv    -> Orange-Red (similar to TssFlnk)
  "#FFDD00",  # 15_EnhBiv    -> Yellow-Orange (similar to EnhA2)
  "#808080",  # 16_ReprPC    -> Dark Gray
  "#A9A9A9",  # 17_ReprPCWk  -> Light Gray
  "#000000"   # 18_Quies     -> Black
)

# Create a named vector to map chromatin states to their colors
chromatin_state_map <- setNames(chromatin_state_colors, chromatin_state_names)

library(tidyr)
library(dplyr)

df <- data.frame()
for(s in names(group.colors )){
  print(s)
  df_s <- readRDS(file = file.path(this.dir,paste0(s,"_histone_summery.rds")))
  print(head(df_s))
  df_long <- df_s %>%
    pivot_longer(cols = -c(biosample, Chromatin_State, n_CpG),
                 names_to = "Variable",
                 values_to = "Value")
  df_long <- df_long %>%
    separate(Variable, into = c("Histone_Modification", "Experiment_ID"), sep = "_")
  df <- rbind(df, df_long)
}

library(ggplot2)
library(dplyr)

# Ensure Chromatin_State is a factor for proper ordering
df$Chromatin_State <- as.factor(df$Chromatin_State)

library(ggplot2)

# Create the heatmap with a 6-row, 3-column facet layout
p <- ggplot(df, aes(x = Histone_Modification, y = biosample, fill = Value)) +
  geom_tile(color = "white") +  # White grid lines
  facet_wrap(~ Chromatin_State, ncol = 3) +  # 3 columns, automatically adjusts rows
  scale_fill_gradient(low = "blue", high = "red") +  # Color scale
  theme_minimal() +  # Use a clean theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(title = "Histone Modification Heatmap",
       x = "Histone Modification",
       y = "Biosample",
       fill = "Value")

# Define output directory
output_dir <- file.path(this.dir, "plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)  # Create dir if not exists

# Save as a large SVG
ggsave(file.path(output_dir, "heatmap.svg"), plot = p, width = 15, height = 10)
