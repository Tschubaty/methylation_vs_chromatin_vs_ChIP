#################################################################
##  beta
##
##  input: everything in HG38
##
##          Encode3\
##
##
##  output:   Encode3\
##
##  v_1 29.09.2024
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
#################################### Libs ########################################
library(tidyr)
library(dplyr)
library(ggplot2)



################################## constants #####################################
start_script <- Sys.time()
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

# Define chromatin state colors
chromatin_state_colors_short <- c(
  "1" = "#FF0000",  "2" = "#FF4500",  "3" = "#FF9900",  "4" = "#FFCC00",
  "5" = "#00CC00",  "6" = "#006400",  "7" = "#FFD700",  "8" = "#FFD700",
  "9" = "#FFFF00",  "10" = "#FFDD00", "11" = "#FFEA73", "12" = "#9370DB",
  "13" = "#C0C0C0", "14" = "#FF4500", "15" = "#FFDD00", "16" = "#808080",
  "17" = "#A9A9A9", "18" = "#000000"
)


# Create a named vector to map chromatin states to their colors
chromatin_state_map <- setNames(chromatin_state_colors, chromatin_state_names)

# View the mapping (optional)
chromatin_state_map


########################### Figure 2E creation ###########################################



# Function to expand semicolon-separated ChIP targets into multiple rows
expand_chip_targets <- function(df, biosample) {
  df %>%
    separate_rows(!!sym(paste0("protein_hits_", biosample)), sep = ";") %>%
    rename(protein_hit = !!sym(paste0("protein_hits_", biosample))) %>%
    filter(!is.na(protein_hit))  # Remove NA values
}

# Process each biosample separately and combine results
biosamples <- c("A549", "GM12878", "HepG2", "K562")
summary_list <- list()


if (file.exists("final_summary.rds")) {
  final_summary <- readRDS("final_summary.rds")
} else {
  # INPUT #
  combined_methylation_data_compact_withpeaks <- readRDS(file.path(this.dir,"combined_methylation_data_compact_withpeaks.rds"))
  
  for (biosample in biosamples) {
    df_bio <- expand_chip_targets(combined_methylation_data_compact_withpeaks, biosample)
    
    summary_df <- df_bio %>%
      group_by(protein_hit) %>%
      mutate(total_cpgs = n()) %>%  # Total CpGs for each ChIP target in the biosample
      group_by(protein_hit, !!sym(paste0("Chromatin_State_", biosample))) %>%
      summarise(
        total_cpgs = first(total_cpgs),  # Keep total CpGs constant for each target
        count_cpgs = n(),  # Count CpGs in this chromatin state
        percentage_cpgs = count_cpgs / total_cpgs * 100
      ) %>%
      rename(chromatin_state = !!sym(paste0("Chromatin_State_", biosample))) %>%
      ungroup()
    
    summary_list[[biosample]] <- summary_df
  }
  
  # Combine results for all biosamples
  final_summary <- bind_rows(summary_list, .id = "biosample")
  
  # View results
  print(final_summary)
  
  saveRDS(object = final_summary,file = "final_summary.rds")
  
  saveRDS(final_summary, "final_summary.rds")
}



# Select top 5 protein hits for each chromatin state
top_hits <- final_summary %>%
  group_by(chromatin_state) %>%
  top_n(5, percentage_cpgs) %>%
  ungroup()


# Create boxplot
ggplot(final_summary, aes(x = factor(chromatin_state), y = percentage_cpgs, fill = biosample)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  # Boxplots without outlier points
  geom_jitter(data = top_hits, aes(color = biosample), width = 0.2, alpha = 0.6, size = 1.5) +  # Jittered points
  geom_text(data = top_hits, aes(label = protein_hit), vjust = -0.5, size = 3) + 
  scale_fill_manual(values = group.colors) +  # Apply custom fill colors
  scale_color_manual(values = group.colors) +  # Apply custom point colors
  labs(
    x = "Chromatin State",
    y = "Percentage of CpGs",
    title = "CpG Distribution Across Chromatin States for ChIP-seq Targets",
    fill = "Biosample",
    color = "Biosample"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(stringr)

# Extract short protein names before aggregation
final_summary2 <- final_summary %>%
  mutate(protein_hit_short = str_extract(protein_hit, "^[^_]+"))  # Extract first part before "_"

# Select top 5 protein hits for each chromatin state
top_hits <- final_summary2 %>%
  group_by(chromatin_state) %>%
  slice_max(order_by = percentage_cpgs, n = 10, with_ties = FALSE) %>%  # Selects top 5 per state
  ungroup()

# Create boxplot
p <- ggplot(final_summary2, aes(x = factor(chromatin_state), y = percentage_cpgs, fill = biosample)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +  # Boxplots without outlier points
  geom_jitter(data = top_hits, aes(color = biosample), width = 0.2, alpha = 0.6, size = 1.5) +  # Jittered points
  geom_text(data = top_hits, aes(label = protein_hit_short), vjust = -0.5, size = 3) +  # Use short label
  scale_fill_manual(values = group.colors) +  # Apply custom fill colors
  scale_color_manual(values = group.colors) +  # Apply custom point colors
  labs(
    x = "Chromatin State",
    y = "Percentage of CpGs",
    title = "CpG Distribution Across Chromatin States for ChIP-seq Targets",
    fill = "Biosample",
    color = "Biosample"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save figure
ggsave("CpG_Chromatin_State_Distribution.png", p, width = 10, height = 6, dpi = 300, limitsize = FALSE)

################# satte 1 
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)  # New package for non-overlapping labels

# Define custom colors
group.colors <- c(
  HepG2    = "#F8766D",  # Warm reddish
  K562     = "#00BFC4",  # Cool cyan
  GM12878  = "#A3A500",  # Yellow-green
  A549     = "#E76BF3"   # Vibrant purple
)

# Filter data to only include chromatin state 1
final_summary_filtered <- final_summary %>%
  filter(chromatin_state == 1) %>%
  mutate(protein_hit_short = str_extract(protein_hit, "^[^_]+"))  # Extract first part before "_"

# Select top 5 protein hits for chromatin state 1 (per biosample)
top_hits <- final_summary_filtered %>%
  group_by(biosample) %>%
  slice_max(order_by = percentage_cpgs, n = 5, with_ties = FALSE) %>%
  ungroup()

# Create boxplot with jittered points and properly aligned non-overlapping labels
p <- ggplot(final_summary_filtered, aes(x = factor(biosample), y = percentage_cpgs, fill = biosample)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, position = position_dodge(width = 0.75)) +  # Align boxplots
  geom_jitter(aes(color = biosample), alpha = 0.6, size = 1.5, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +  # Jitter within biosample columns
  geom_text_repel(data = top_hits, aes(label = protein_hit_short, color = biosample), 
                  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
                  size = 3, max.overlaps = Inf, box.padding = 0.3, force = 5) +  # Non-overlapping labels
  scale_fill_manual(values = group.colors) +  # Apply custom fill colors
  scale_color_manual(values = group.colors) +  # Apply custom point colors
  labs(
    x = "Biosample",
    y = "Percentage of CpGs",
    title = "CpG Distribution in Chromatin State 1 for ChIP-seq Targets",
    fill = "Biosample",
    color = "Biosample"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save figure
ggsave("CpG_Chromatin_State_1_Distribution.png", p, width = 10, height = 6, dpi = 300, limitsize = FALSE)


saveRDS(object = final_summary_filtered ,file = "final_summary_filtered.rds")



###### adding one plöot per proetien:
library(dplyr)
library(ggplot2)
library(ggrepel)
# plotting toggel
output_dir_plots <- file.path(this.dir, "protein_chromatin_state_plots")
# single output folder
dir.create(output_dir_plots, recursive = TRUE, showWarnings = FALSE)
if (FALSE) {

  
  prot_list <- sort(unique(final_summary2$protein_hit_short))
  
  for (prot in prot_list) {
    print(prot)
    lab_df <- final_summary2 %>% filter(protein_hit_short == prot)
    
    p <- ggplot(final_summary2,
                aes(
                  x = factor(chromatin_state),
                  y = percentage_cpgs,
                  fill = biosample
                )) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
      geom_jitter(
        data = lab_df,
        mapping = aes(color = biosample),
        width = 0.2,
        alpha = 0.6,
        size = 1.5
      ) +
      ggtitle(label = prot) +
      # geom_text_repel(
      #   data = lab_df,
      #   aes(label = protein_hit, color = biosample),
      #   size = 1.8,                 # smaller font (~6 pt)
      #   box.padding = 0.15,         # tighter box around labels
      #   point.padding = 0.1,        # tighter point padding
      #   max.overlaps = Inf,         # don’t drop labels
      #   segment.size = 0.2,         # thinner connecting lines
      #   segment.alpha = 0.5         # lighter connecting lines
      # )+
      scale_fill_manual(values = group.colors) +
      scale_color_manual(values = group.colors) +
      labs(
        x = "Chromatin State",
        y = "Percentage of CpGs",
        title = paste0("CpG Distribution Across Chromatin States – ", prot),
        fill = "Biosample",
        color = "Biosample"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      filename = file.path(
        output_dir_plots,
        paste0("CpG_Chromatin_State_Distribution_", prot, ".png")
      ),
      plot = p,
      width = 20,
      height = 12,
      dpi = 300,
      limitsize = FALSE
    )
  }
}
# library(dplyr)
# library(ggplot2)
# library(ggrepel)
# 
# # fixed order + offsets per biosample (4 lanes per state)
# lane_offset <- c(A549 = -0.30, GM12878 = -0.10, HepG2 = 0.10, K562 = 0.30)
# 
# final_summary2 <- final_summary2 %>%
#   mutate(
#     biosample = factor(biosample, levels = biosamples),
#     chrom_state_f = factor(chromatin_state, levels = sort(unique(chromatin_state))),
#     x_pos = as.numeric(chrom_state_f) + lane_offset[as.character(biosample)]
#   )
# 
# output_dir_plots <- file.path(this.dir,"protein_chromatin_state_plots")
# dir.create(output_dir_plots, recursive = TRUE, showWarnings = FALSE)
# 
# prot_list <- sort(unique(final_summary2$protein_hit_short))
# 
# for (prot in prot_list) {
#   print(prot)
#   lab_df <- final_summary2 %>%
#     filter(protein_hit_short == prot) %>%
#     mutate(x_pos = as.numeric(chrom_state_f) + lane_offset[as.character(biosample)])
#   
#   p <- ggplot(final_summary2, aes(x = x_pos, y = percentage_cpgs)) +
#     geom_boxplot(
#       aes(fill = biosample, group = interaction(chrom_state_f, biosample)),
#       width = 0.18, outlier.shape = NA, alpha = 0.5, na.rm = TRUE
#     ) +
#     geom_jitter(
#       data = lab_df, aes(color = biosample),
#       width = 0.03, height = 0, size = 1.5, alpha = 0.6, na.rm = TRUE
#     ) +
#     # geom_text_repel(
#     #   data = lab_df, aes(label = protein_hit, color = biosample),
#     #   size = 1.6, box.padding = 0.12, point.padding = 0.08,
#     #   segment.size = 0.2, segment.alpha = 0.5,
#     #   max.overlaps = Inf, na.rm = TRUE
#     # ) +
#     scale_fill_manual(values = group.colors, limits = biosamples, drop = FALSE) +
#     scale_color_manual(values = group.colors, limits = biosamples, drop = FALSE) +
#     scale_x_continuous(
#       breaks = seq_along(levels(final_summary2$chrom_state_f)),
#       labels = levels(final_summary2$chrom_state_f)
#     ) +
#     labs(
#       x = "Chromatin State", y = "Percentage of CpGs",
#       title = paste0("CpG Distribution Across Chromatin States – ", prot),
#       fill = "Biosample", color = "Biosample"
#     ) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
#   ggsave(
#     file.path(output_dir_plots, paste0("CpG_Chromatin_State_Distribution_", prot, ".png")),
#     p, width = 20, height = 12, dpi = 300, limitsize = FALSE
#   )
# }

#######################################################################################

# v01 unknown only one categroy 
# 
# summarize_chromatin_states <- function(file_path, biosample) {
#   df_ChIP <- read_bed_file(file_path)
#   
#   # Extract the column for the biosample's chromatin state
#   chromatin_col <- paste0("Chromatin_State_", biosample)
#   
#   # Replace missing or "NA" states with "Unknown"
#   df_ChIP[[chromatin_col]] <- ifelse(
#     is.na(df_ChIP[[chromatin_col]]) | df_ChIP[[chromatin_col]] == ".",
#     "Unknown",
#     df_ChIP[[chromatin_col]]
#   )
#   
#   # Summarize the distribution of chromatin states
#   summary <- df_ChIP %>%
#     mutate(!!sym(chromatin_col) := ifelse(
#       is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == ".",
#       "Unknown",  # Assign a label for missing or NA states
#       !!sym(chromatin_col)
#     )) %>%
#     group_by(!!sym(chromatin_col)) %>%
#     summarize(Count = n(), .groups = "drop") %>%
#     arrange(desc(Count))
#   
#   # Ensure all chromatin states are present, including "Unknown"
#   chromatin_states <- tibble(
#     Chromatin_State = c(names(chromatin_state_colors_short), "Unknown")
#   )
#   
#   summary <- chromatin_states %>%
#     left_join(summary %>% rename(Chromatin_State = !!sym(chromatin_col)), 
#               by = "Chromatin_State") %>%
#     mutate(Count = replace_na(Count, 0)) %>%
#     arrange(match(Chromatin_State, c(names(chromatin_state_colors_short), "Unknown")))
#   
#   # Add percentages column to the summary
#   summary <- summary %>%
#     mutate(
#       Percentage = (Count / sum(Count)) * 100
#     )
#   
#   
#   return(summary)
# }
# 
# 
# # Initialize a list to store summaries for all experiments
# all_experiment_summaries <- list()
# 
# # Iterate through all proteins and their respective experiments
# proteins <- list.files(path = input_ChIP_dir)
# for (protein in proteins) {
#   print(protein)
#   # skil log folder 
#   if (protein == "log") next
#   # Get all files for the protein
#   protein_dir <- file.path(input_ChIP_dir, protein)
#   experiment_files <- list.files(path = protein_dir, full.names = TRUE)
#   
#   for (file_path in experiment_files) {
#     # Extract experiment details from file name
#     file_name <- basename(file_path)
#     experiment_id <- strsplit(file_name, "_")[[1]][3]
#     biosample <- strsplit(file_name, "_")[[1]][1]
#     
#     # Summarize chromatin states for this experiment
#     summary <- summarize_chromatin_states(file_path, biosample)
#     
#     # Add experiment details
#     summary <- summary %>%
#       mutate(Protein = protein, ExperimentID = experiment_id, Biosample = biosample)
#     
#     # Append to the list
#     all_experiment_summaries[[paste0(protein, "_", experiment_id)]] <- summary
#   }
# }
# 
# # Combine all summaries into a single data frame
# final_summary <- bind_rows(all_experiment_summaries)
# 
# # Ensure Chromatin_State is ordered
# final_summary <- final_summary %>%
#   mutate(
#     Chromatin_State = factor(
#       Chromatin_State,
#       levels = c(names(chromatin_state_colors_short), "Unknown")
#     )
#   )
# 
# # Save the summary to a CSV file
# output_csv <- file.path("chromatin_state_summary.csv")
# write_csv(final_summary, output_csv)

################################################################################
# details unknown categtory 
summarize_chromatin_states <- function(file_path, biosample) {
  df_ChIP <- read_bed_file(file_path)
  
  # Extract the column for the biosample's chromatin state
  chromatin_col <- paste0("Chromatin_State_", biosample)
  
  # Categorize rows based on start_motif, start_cg, and chromatin state
  df_ChIP <- df_ChIP %>%
    mutate(
      !!sym(chromatin_col) := case_when(
        (is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == ".") & start_motif == -1 ~ "No motif in peak",
        (is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == ".") & start_motif != -1 & start_cg == 0 ~ "No CG in motif",
        (is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == ".") & start_motif != -1 & start_cg != 0 ~ "No State assignment for CG",
        TRUE ~ !!sym(chromatin_col)  # Keep original value otherwise
      )
    )
  
  # Summarize the distribution of chromatin states
  summary <- df_ChIP %>%
    group_by(!!sym(chromatin_col)) %>%
    summarize(Count = n(), .groups = "drop") %>%
    arrange(desc(Count))
  
  # Ensure all chromatin states are present
  chromatin_states <- tibble(
    Chromatin_State = c(names(chromatin_state_colors_short), 
                        "No motif in peak", "No CG in motif", "No State assignment for CG")
  )
  
  summary <- chromatin_states %>%
    left_join(summary %>% rename(Chromatin_State = !!sym(chromatin_col)), 
              by = "Chromatin_State") %>%
    mutate(Count = replace_na(Count, 0)) %>%
    arrange(match(Chromatin_State, c(names(chromatin_state_colors_short), 
                                     "No motif in peak", "No CG in motif", "No State assignment for CG")))
  
  # Add percentages column to the summary
  summary <- summary %>%
    mutate(
      Percentage = (Count / sum(Count)) * 100
    )
  
  return(summary)
}

## bugfinxin this seciton 

input_ChIP_dir <- file.path(dirname(this.dir), "meme", "fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3")

# Function to extract and clean column names from the first row of a text file
extract_colnames <- function(file_path) {
  # Read the first row
  first_row <- readLines(file_path, n = 1)
  
  # Remove the leading "#" and replace double tabs with single tab
  cleaned_row <- gsub("#", "", first_row)          # Remove leading #
  
  # Split the cleaned string into individual column names
  colnames <- unlist(strsplit(cleaned_row, "\t"))
  
  return(colnames)
}

# Function to read .bed files with dynamically extracted column names and correct column types
read_bed_file <- function(file_path) {
  # Extract the correct column names from the first row
  colnames <- extract_colnames(file_path)
  
  # Read the BED file, skipping the first line (header row)
  df <- readr::read_delim(
    file_path,
    delim = "\t",
    skip = 1,
    col_names = colnames,
    show_col_types = FALSE
  )
  
  return(df)
}

# Initialize a list to store summaries for all experiments
all_experiment_summaries <- list()

# Iterate through all proteins and their respective experiments
proteins <- list.files(path = input_ChIP_dir)
for (protein in proteins) {
  print(protein)
  # Skip log folder
  if (protein == "log") next
  # Get all files for the protein
  protein_dir <- file.path(input_ChIP_dir, protein)
  experiment_files <- list.files(path = protein_dir, full.names = TRUE)
  
  for (file_path in experiment_files) {
    # Extract experiment details from file name
    file_name <- basename(file_path)
    experiment_id <- strsplit(file_name, "_")[[1]][3]
    biosample <- strsplit(file_name, "_")[[1]][1]
    
    # Summarize chromatin states for this experiment
    summary <- summarize_chromatin_states(file_path, biosample)
    
    # Add experiment details
    summary <- summary %>%
      mutate(Protein = protein, ExperimentID = experiment_id, Biosample = biosample)
    
    # Append to the list
    all_experiment_summaries[[paste0(protein, "_", experiment_id)]] <- summary
  }
}

# Combine all summaries into a single data frame
final_summary <- bind_rows(all_experiment_summaries)

# Ensure Chromatin_State is ordered
final_summary <- final_summary %>%
  mutate(
    Chromatin_State = factor(
      Chromatin_State,
      levels = c(names(chromatin_state_colors_short), 
                 "No motif in peak", "No CG in motif", "No State assignment for CG")
    )
  )

# Save the summary to a CSV file
output_csv <- file.path("chromatin_state_summary.csv")
readr::write_csv(final_summary, output_csv)

###############################################################################
# detailed unknown cide 

library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

# Path to the saved file
output_csv <- file.path(output_dir_plots, "chromatin_state_summary.csv")

# Load the CSV file into a dataframe
final_summary <- read_csv(output_csv)


library(dplyr)

usable_summary <- final_summary %>%
  filter(as.character(Chromatin_State) %in% as.character(1:18)) %>%
  group_by(Protein, ExperimentID, Biosample) %>%
  summarize(
    Chromatin_State = "Usable_1_18",
    Count = sum(Count),
    Percentage = sum(Percentage),
    .groups = "drop"
  ) %>%
  select(Chromatin_State, Count, Percentage, Protein, ExperimentID, Biosample)

final_summary <- bind_rows(
  final_summary %>%
    mutate(Chromatin_State = as.character(Chromatin_State)),
  usable_summary
) %>%
  mutate(
    Chromatin_State = factor(
      Chromatin_State,
      levels = c(as.character(1:18),
                 "No motif in peak",
                 "No CG in motif",
                 "No State assignment for CG",
                 "Usable_1_18")
    )
  )


# Define chromatin state names and colors
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
  "#000000",  # 18_Quies     -> Black
  "#C2185B",  # No motif in peak -> Pink
  "#E91E63",  # No CG in motif -> Dodger Blue
  "#F48FB1",   # No State assignment for CG -> Lime Green,
  "#1E90FF"  # "Usable_1_18"
  
)

chromatin_state_labels <- c(
  "1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD",
  "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2",
  "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF_Rpts",
  "13_Het", "14_TssBiv", "15_EnhBiv", "16_ReprPC",
  "17_ReprPCWk", "18_Quies", 
  "No motif in peak", "No CG in motif", "No State assignment for CG","Usable_1_18"
)

# Ensure Chromatin_State is ordered
final_summary <- final_summary %>%
  mutate(
    Chromatin_State = factor(
      Chromatin_State,
      levels = c(as.character(1:18), "No motif in peak", "No CG in motif", "No State assignment for CG","Usable_1_18")
    )
  )

# Calculate z-scores and identify outliers
outlier_summary <- final_summary %>%
  group_by(Chromatin_State) %>%
  mutate(
    Mean_Percentage = mean(Percentage),
    SD_Percentage = sd(Percentage),
    Z_Score = (Percentage - Mean_Percentage) / SD_Percentage,
    Is_Outlier = abs(Z_Score) > 3,  # Outlier if z-score > threshold
    Label = ifelse(Is_Outlier, paste(Protein, Biosample, sep = " "), NA)
  ) %>%
  ungroup()

# Create final plot with outlier labels using ggrepel
ggplot(outlier_summary, aes(x = Chromatin_State, y = Percentage, color = Chromatin_State)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Jitter points
  # geom_text_repel(
  #   aes(label = Label), 
  #   size = 3, color = "black", na.rm = TRUE, 
  #   box.padding = 0.3, point.padding = 0.2, max.overlaps = Inf
  # ) +  # Repelling labels for outliers
  labs(
    title = "Chromatin State Percentage Distribution with Outliers",
    x = "Chromatin State",
    y = "Percentage"
  ) +
  scale_x_discrete(labels = chromatin_state_labels) +  # Set discrete x-axis labels with full names
  scale_color_manual(values = chromatin_state_colors) +  # Apply colors to Chromatin States
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none"  # Drop legend
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))  # Format y-axis as percentage

################################################################################
# library(ggplot2)
# library(ggrepel)
# library(dplyr)
# library(readr)
# 
# # Path to the saved file
# output_csv <- file.path(out_dir, "chromatin_state_summary.csv")
# 
# # Load the CSV file into a dataframe
# final_summary <- read_csv(output_csv)
# 
# # Define chromatin state names and colors
# chromatin_state_colors <- c(
#   "#FF0000",  # 1_TssA       -> Red
#   "#FF4500",  # 2_TssFlnk    -> Orange-Red
#   "#FF9900",  # 3_TssFlnkU   -> Orange
#   "#FFCC00",  # 4_TssFlnkD   -> Yellow-Orange
#   "#00CC00",  # 5_Tx         -> Green
#   "#006400",  # 6_TxWk       -> Dark Green
#   "#FFD700",  # 7_EnhG1      -> Gold
#   "#FFD700",  # 8_EnhG2      -> Gold (same as EnhG1)
#   "#FFFF00",  # 9_EnhA1      -> Yellow
#   "#FFDD00",  # 10_EnhA2     -> Yellow-Orange
#   "#FFEA73",  # 11_EnhWk     -> Light Yellow
#   "#9370DB",  # 12_ZNF_Rpts  -> Purple
#   "#C0C0C0",  # 13_Het       -> Light Gray
#   "#FF4500",  # 14_TssBiv    -> Orange-Red (similar to TssFlnk)
#   "#FFDD00",  # 15_EnhBiv    -> Yellow-Orange (similar to EnhA2)
#   "#808080",  # 16_ReprPC    -> Dark Gray
#   "#A9A9A9",  # 17_ReprPCWk  -> Light Gray
#   "#000000",  # 18_Quies     -> Black
#   "#808080"   # Unknown      -> Light Gray (default for unknown)
# )
# 
# chromatin_state_labels <- c(
#   "1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD",
#   "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2",
#   "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF_Rpts",
#   "13_Het", "14_TssBiv", "15_EnhBiv", "16_ReprPC",
#   "17_ReprPCWk", "18_Quies", "Unknown"
# )
# 
# # Ensure Chromatin_State is ordered
# final_summary <- final_summary %>%
#   mutate(
#     Chromatin_State = factor(
#       Chromatin_State,
#       levels = c(as.character(1:18), "Unknown")
#     )
#   )
# 
# # Calculate z-scores and identify outliers
# outlier_summary <- final_summary %>%
#   group_by(Chromatin_State) %>%
#   mutate(
#     Mean_Percentage = mean(Percentage),
#     SD_Percentage = sd(Percentage),
#     Z_Score = (Percentage - Mean_Percentage) / SD_Percentage,
#     Is_Outlier = abs(Z_Score) > 3,  # Outlier if z-score > threshold
#     Label = ifelse(Is_Outlier, paste(Protein, Biosample, sep = " "), NA)
#   ) %>%
#   ungroup()
# 
# # Create final plot with outlier labels using ggrepel
# ggplot(outlier_summary, aes(x = Chromatin_State, y = Percentage, color = Chromatin_State)) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Jitter points
#   geom_text_repel(
#     aes(label = Label), 
#     size = 3, color = "black", na.rm = TRUE, 
#     box.padding = 0.3, point.padding = 0.2, max.overlaps = Inf
#   ) +  # Repelling labels for outliers
#   labs(
#     title = "Chromatin State Percentage Distribution with Outliers",
#     x = "Chromatin State",
#     y = "Percentage"
#   ) +
#   scale_x_discrete(labels = chromatin_state_labels) +  # Set discrete x-axis labels with full names
#   scale_color_manual(values = chromatin_state_colors) +  # Apply colors to Chromatin States
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
#     legend.position = "none"  # Drop legend
#   ) +
#   scale_y_continuous(labels = scales::percent_format(scale = 1))  # Format y-axis as percentage
###################################### sibnglke protein ###########################
# Define output directory for plots

protein_plot_dir <- file.path(output_dir_plots, "protein_plots")
dir.create(protein_plot_dir, recursive = TRUE, showWarnings = FALSE)

proteins <- unique(final_summary$Protein)


if (FALSE) {
  # plot toggel
  
  
  for (protein in proteins) {
    protein_data <- final_summary %>%
      mutate(Label = ifelse(Protein == protein, paste(Protein, Biosample, sep = " "), NA))
    
    final_plot <- ggplot(protein_data,
                         aes(x = Chromatin_State, y = Percentage, color = Chromatin_State)) +
      geom_jitter(width = 0.2,
                  size = 2,
                  alpha = 0.7) +
      geom_text_repel(
        aes(label = Label),
        size = 3,
        color = "black",
        na.rm = TRUE,
        box.padding = 0.3,
        point.padding = 0.2,
        max.overlaps = Inf
      ) +
      labs(
        title = paste("Chromatin State Distribution for Protein:", protein),
        x = "Chromatin State",
        y = "Percentage"
      ) +
      scale_color_manual(values = chromatin_state_colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      scale_y_continuous(labels = scales::percent_format(scale = 1))
    
    plot_file <- file.path(protein_plot_dir,
                           paste0(protein, "_chromatin_state_plot.png"))
    ggsave(
      plot_file,
      final_plot,
      width = 24,
      height = 12,
      limitsize = FALSE
    )
    
    cat("Plot saved for protein:", protein, "\n")
  }
}
# get destributon fopr QC 

special_categories <- c(
  "No motif in peak",
  "No CG in motif",
  "No State assignment for CG",
  "Usable_1_18"
)

# colors exactly as you want
special_colors <- c(
  "No motif in peak" = "#C2185B",
  "No CG in motif" = "#E91E63",
  "No State assignment for CG" = "#F48FB1",
  "Usable_1_18" = "#1E90FF"
)

# adjustable number of bins
n_bins <- 100   # <- change this to control resolution
special_data <- final_summary %>%
  filter(as.character(Chromatin_State) %in% special_categories) %>%
  mutate(Chromatin_State = factor(as.character(Chromatin_State), levels = special_categories))

p <- ggplot(special_data, aes(x = Percentage, fill = Chromatin_State)) +
  geom_histogram(
    bins = n_bins,
    color = "black",
    alpha = 0.8
  ) +
  facet_wrap(~Chromatin_State, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = special_colors) +
  labs(
    title = paste("Histogram distributions (bins =", n_bins, ")"),
    x = "Percentage",
    y = "Count"
  ) +
  theme_minimal() +
  scale_x_continuous(labels = scales::percent_format(scale = 1))

print(p)

ggsave(
  file.path(output_dir_plots, paste0("special_histograms_bins_", n_bins, ".png")),
  p,
  width = 14,
  height = 10,
  limitsize = FALSE
)


#################################
# create subfolder once before loop
protein_plot_dir <- file.path(output_dir_plots, "protein_barplots")
dir.create(protein_plot_dir, recursive = TRUE, showWarnings = FALSE)

#protein_to_plot <- "ARNT"

for (protein_to_plot in proteins) {
  names(chromatin_state_colors) <- c(
    as.character(1:18),
    "No motif in peak",
    "No CG in motif",
    "No State assignment for CG",
    "Usable_1_18"
  )
  
  plot_order <- c("Chromatin states 1-18", "Special categories")
  
  protein_bar_data <- final_summary %>%
    mutate(
      Chromatin_State = as.character(Chromatin_State),
      Group = paste(Biosample, ExperimentID, sep = "\n"),
      PlotCategory = case_when(
        Chromatin_State %in% as.character(1:18) ~ "Chromatin states 1-18",
        Chromatin_State %in% c("No motif in peak", "No CG in motif", "No State assignment for CG") ~ "Special categories",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(Protein == protein_to_plot, !is.na(PlotCategory)) %>%
    mutate(
      PlotCategory = factor(PlotCategory, levels = plot_order),
      Chromatin_State = factor(
        Chromatin_State,
        levels = c(
          as.character(1:18),
          "No motif in peak",
          "No CG in motif",
          "No State assignment for CG"
        )
      ),
      XGroup = paste(PlotCategory, Group, sep = "___")
    )
  
  x_levels <- protein_bar_data %>%
    distinct(PlotCategory, Group, XGroup) %>%
    arrange(PlotCategory, Group) %>%
    pull(XGroup)
  
  protein_bar_data <- protein_bar_data %>%
    mutate(XGroup = factor(XGroup, levels = x_levels))
  
  x_labels_df <- protein_bar_data %>%
    distinct(XGroup, PlotCategory, Group) %>%
    arrange(PlotCategory, Group) %>%
    mutate(Label = paste0(as.character(PlotCategory), "\n", Group))
  
  x_labels <- x_labels_df$Label
  names(x_labels) <- x_labels_df$XGroup
  
  p <- ggplot(protein_bar_data,
              aes(x = XGroup, y = Count, fill = Chromatin_State)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_manual(values = chromatin_state_colors[c(as.character(1:18),
                                                        "No motif in peak",
                                                        "No CG in motif",
                                                        "No State assignment for CG")], drop = FALSE) +
    scale_x_discrete(labels = x_labels) +
    labs(
      title = paste("Absolute counts for protein:", protein_to_plot),
      x = NULL,
      y = "Count"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    legend.position = "right")
  
  print(p)
  
  # updated ggsave
  ggsave(
    filename = file.path(
      protein_plot_dir,
      paste0("barplot_", protein_to_plot, "_absolute_counts.png")
    ),
    plot = p,
    width = 20,
    height = 10,
    limitsize = FALSE
  )
  cat("Saved plot for protein:", protein_to_plot, "\n")
}  
  ##### filtering for QC
  
  motif_threshold <- 50
  cg_threshold <- 80
  
  special_categories <- c("No motif in peak",
                          "No CG in motif",
                          "No State assignment for CG",
                          "Usable_1_18")
  
  qc_data <- final_summary %>%
    filter(as.character(Chromatin_State) %in% special_categories) %>%
    mutate(Chromatin_State = as.character(Chromatin_State)) %>%
    
    group_by(Protein, ExperimentID, Biosample) %>%
    mutate(
      NoMotif = Percentage[Chromatin_State == "No motif in peak"][1],
      NoCG    = Percentage[Chromatin_State == "No CG in motif"][1],
      
      QC_Status = case_when(
        NoMotif >= motif_threshold ~ "FAIL_MOTIF",
        # priority 1
        NoCG >= cg_threshold ~ "FAIL_CG",
        # priority 2
        TRUE ~ "PASS"
      )
    ) %>%
    ungroup() %>%
    mutate(
      QC_Status = factor(QC_Status, levels = c("PASS", "FAIL_CG", "FAIL_MOTIF")),
      Chromatin_State = factor(Chromatin_State, levels = special_categories)
    )
  
  p <- ggplot(qc_data, aes(x = Percentage, fill = QC_Status)) +
    geom_histogram(bins = n_bins,
                   color = "black",
                   alpha = 0.9) +
    facet_wrap( ~ Chromatin_State, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(
      "PASS" = "#4CAF50",
      # green
      "FAIL_CG" = "#FF9800",
      # orange
      "FAIL_MOTIF" = "#D32F2F"  # red
    )) +
    labs(
      title = paste0(
        "QC reasons (NO Motif < ",
        motif_threshold,
        "%, NO CG < ",
        cg_threshold,
        "%)"
      ),
      x = "Percentage",
      y = "Count"
    ) +
    theme_minimal() +
    scale_x_continuous(labels = scales::percent_format(scale = 1))
  
  print(p)
  
  ggsave(
    file.path(
      output_dir_plots,
      paste0(
        "QC_reason_simple_motif_",
        motif_threshold,
        "_cg_",
        cg_threshold,
        ".png"
      )
    ),
    p,
    width = 14,
    height = 10,
    dpi = 300,
    limitsize = FALSE
  )
  


#########
low_cg_threshold <- 5     # very low No CG
high_motif_threshold <- 70  # very high No motif

interesting_cases <- final_summary %>%
  mutate(Chromatin_State = as.character(Chromatin_State)) %>%
  
  group_by(Protein, ExperimentID, Biosample) %>%
  summarize(
    NoMotif = Percentage[Chromatin_State == "No motif in peak"][1],
    NoCG    = Percentage[Chromatin_State == "No CG in motif"][1],
    .groups = "drop"
  ) %>%
  
  filter(
    NoCG < low_cg_threshold,
    NoMotif > high_motif_threshold
  ) %>%
  
  arrange(desc(NoMotif))

print(interesting_cases)


full_rows <- final_summary %>%
  semi_join(
    interesting_cases,
    by = c("Protein", "ExperimentID", "Biosample")
  )

print(full_rows, n = 200)



# Highest usable among FAIL_MOTIF
qc_data %>%
  filter(Chromatin_State == "Usable_1_18", QC_Status == "FAIL_MOTIF") %>%
  arrange(desc(Percentage)) %>%
  print(n = 50)


final_summary %>%
    filter(
        Protein == "ZBTB33")


qc_data %>%
      filter(Protein == "ZBTB33") 
