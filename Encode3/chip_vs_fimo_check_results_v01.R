#################################################################
##  Chip-seq vs methylation vs Chromatin
##
##  input: everything in HG38
##
##
##
##  output:   Encode3\simulation
##
##  v_1 08.11.2024
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
# Load necessary libraries
library(stringr)
library(readr)
library(foreach)
library(doParallel)
library(dplyr)
library(rlang)
library(ggplot2)
library(tidyr)
library(VennDiagram)
library(gridExtra)
##################################### INPUT ########################################



################################## constants #####################################
start_script <- Sys.time()
# Set input and output directories
input_folder <- file.path(
  this.dir,
  "meme",
  "fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3"
)
output_folder <- file.path(this.dir, "simulation")
RNA_seq_folder <- "RNA-seq"

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

# readLines(bed_file, n = 1)
#file_path <- bed_file
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

# Get the list of all files matching the pattern *_methylation.chromatinstate.fimo.bed
protein_folders <- list.dirs(path = input_folder,
                             recursive = FALSE ,
                             full.names = FALSE)
#List all protein folders, excluding "sbatch_scripts" and "log"
protein_folders <- protein_folders[!grepl("sbatch_scripts|log", protein_folders)]

summary_df <- readRDS(file = file.path(output_folder, "summary_df.rds"))

# Filter summary_df to keep only the rows with the lowest no_motif_no_CpG_ratio for each (protein, biosample, motif) combination
filtered_summary_df <- summary_df %>%
  group_by(protein, biosample, motif) %>%
  slice_min(no_motif_no_CpG_ratio, with_ties = FALSE) %>%
  ungroup()

# Function to calculate S_statistic for a given data frame
calculate_S_statistic <- function(df, biosample_comb) {
  mean_diff_per_sample <- df %>%
    mutate(diff_fRead = !!sym(paste0("fRead_", biosample_comb[1])) -!!sym(paste0("fRead_", biosample_comb[2]))) %>%  # Calculate difference
    group_by(sample) %>%  # Group by sample
    summarize(mean_diff = mean(diff_fRead, na.rm = TRUE))  # Calculate mean difference per sample
  
  # Calculate the S statistic
  S_statistic <- mean_diff_per_sample$mean_diff[biosample_comb[1] == mean_diff_per_sample$sample] -
    mean_diff_per_sample$mean_diff[biosample_comb[2] == mean_diff_per_sample$sample]
  
  return(S_statistic)
}

summary_test <- readRDS(file = file.path(output_folder, "summary_test.rds"))

# Loop over protein folders
for (protein in protein_folders) {
  # [1:1]
  # protein <- "CTCF"
  # protein <- "CEBPB"
  # protein <- "MAX"
  print(protein)
  
  protein_data <- filtered_summary_df[filtered_summary_df$protein == protein, ]
  # Filter the data for the specific protein
  
  # Loop over unique motifs
  unique_motifs <- unique(protein_data$motif)
  
  for (motif in unique_motifs) {
    # Filter the data for the current motif
    # debug # motif <- unique_motifs[1]
    print(motif)
    motif_data <- protein_data[protein_data$motif == motif, ]
    
    # Generate all combinations of two biosamples
    biosample_combinations <- combn(unique(motif_data$biosample), 2, simplify = FALSE)
    
    # Loop through each combination of two biosamples
    for (biosample_comb in biosample_combinations) {
      
      # biosample_comb <- biosample_combinations[[1]]
      
      print(biosample_comb)
      
      # Dynamically generate the chromatin state column names for both biosamples
      chromatin_state_col1 <- paste0("Chromatin_State_", biosample_comb[1])
      chromatin_state_col2 <- paste0("Chromatin_State_", biosample_comb[2])
      
      
      merged_df <- readRDS(file = file.path(
        output_folder,
        protein,
        paste(
          unlist(biosample_comb)[1],
          unlist(biosample_comb)[2],
          protein,
          motif,
          "best_comparison_all_data_points.rds",
          sep = "_"
        )
      ))
      
      biosample_both <- paste(unlist(biosample_comb)[1], unlist(biosample_comb)[2], sep = "_")
      
      df_complete <- merged_df[complete.cases(merged_df), ]
      
      # Filter rows where chromatin states are the same for both samples
      df_same_state <- df_complete %>%
        filter(!!sym(chromatin_state_col1) == !!sym(chromatin_state_col2))
      
      # Filter rows where chromatin states are "1" in both samples
      df_state1 <- df_complete %>%
        filter(!!sym(chromatin_state_col1) == "1" & !!sym(chromatin_state_col2) == "1")
      
      
      sample_counts <- table(df_complete$sample)
      sample_counts_same_state <- table(df_same_state$sample)
      sample_counts_state1 <- table(df_state1$sample)
     
      # Dynamically generate the column names for fRead values
      fRead_col1 <- paste0("fRead_", biosample_comb[1])
      fRead_col2 <- paste0("fRead_", biosample_comb[2])
      
      ### ploting 
      
      # Function to create a heatmap for a specific sample
      create_heatmap <- function(sample_name) {
        # Filter for the specified sample
        df_sample <- df_state1 %>% filter(sample == sample_name)
        
        # Define colors for the gradient, defaulting to grey if sample is not in group.colors
        light_color <- "#f7fbff"
        sample_color <- ifelse(sample_name %in% names(group.colors), group.colors[sample_name], "grey")
        intermediate_color <- colorRampPalette(c(light_color, sample_color))(3)[2]
        
        # Define custom colors and values for gradient
        custom_colors <- c(light_color, intermediate_color, sample_color)
        custom_values <- c(0, 0.01, 1)
        
        # Create the 2D heatmap plot with custom color gradient and overlay all points
        ggplot(df_sample, aes_string(x = fRead_col1, y = fRead_col2)) +
          stat_bin2d(binwidth = c(0.025, 0.025), aes(fill = ..count..)) +  # Set bin width to ensure consistent binning
          geom_point(color = sample_color, size = 1, alpha = 0.5) +
          labs(
            title = paste("ChIP in", sample_name),
            x = paste("fRead_", biosample_comb[1]),
            y = paste("fRead_", biosample_comb[2])
          ) +
          scale_fill_gradientn(colors = custom_colors, values = custom_values) +
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +  # Set axis limits without cutting off data
          theme_minimal()
      }
      
      # Create individual plots for each of the specified samples
      plot1 <- create_heatmap(biosample_comb[1])
      plot2 <- create_heatmap(paste(biosample_comb[1], biosample_comb[2], sep = "_"))
      plot3 <- create_heatmap(biosample_comb[2])
      
      # Combine the plots into a single grid with a super title
      combined_plot <- arrangeGrob(
        plot1, plot2, plot3,
        ncol = 3,
        top = textGrob(motif, gp = gpar(fontsize = 16, fontface = "bold"))
      )
      
      # Save the combined plot as an image
      output_file <- file.path(
        output_folder,
        protein,
        motif,
        paste(motif, biosample_comb[1], biosample_comb[2], "combined_heatmap.png", sep = "_")
      )
      
      ggsave(
        filename = output_file,
        plot = combined_plot,  # Use combined plot with title
        width = 15, height = 5  # Adjust width and height as needed
      )

      
      # Ensure that samples not in `group.colors` default to "grey"
      df_state1 <- df_state1 %>%
        mutate(color = ifelse(sample %in% names(group.colors), group.colors[sample], "grey"))
      
      # Define the output file path
      output_file <- file.path(
        output_folder,
        protein,
        motif,
        paste(motif, biosample_comb[1], biosample_comb[2], "scatter_plot_all.png", sep = "_")
      )
      
      # Create the plot
      scatter_plot <- ggplot(df_state1, aes_string(x = fRead_col1, y = fRead_col2)) +
        geom_point(aes(color = sample), size = 1, alpha = 0.6) +  # Map `sample` to the color legend
        scale_color_manual(values = c(group.colors, grey = "grey")) +  # Use custom colors with grey as default
        labs(
          title = motif,
          x = fRead_col1,
          y = fRead_col2,
          color = "Sample"  # Legend title
        ) +
        theme_minimal()
      
      # Save the plot as an image
      ggsave(
        filename = output_file,
        plot = scatter_plot,  # Use scatter plot with custom color legend
        width = 10, height = 8  # Adjust width and height as needed
      )
      
      # # Create the scatter plot
      # ggplot(df_state1, aes_string(x = fRead_col1, y = fRead_col2, color = "sample")) +
      #   geom_point() +
      #   labs(
      #     title = paste(protein,motif,"Scatter Plot of fRead_", biosample_comb[1], " vs fRead_", biosample_comb[2]),
      #     x = paste("fRead_", biosample_comb[1]),
      #     y = paste("fRead_", biosample_comb[2])
      #   ) +
      #   theme_minimal()
      # 
      # 
      # # Create the 2D histogram plot and build it to extract bin data
      # df_bins_plot <- ggplot(df_state1, aes_string(x = fRead_col1, y = fRead_col2)) +
      #   stat_bin2d(bins = 30)
      # 
      # # Use ggplot_build to extract binned data
      # df_bins <- ggplot_build(df_bins_plot)$data[[1]]
      # 
      # # Convert the binned data to a data frame and rename columns
      # df_3d <- df_bins %>%
      #   select(x, y, count) %>%
      #   rename(z = count)
      # 
      # # Apply a log10 transformation to the z values
      # df_3d <- df_3d %>%
      #   mutate(z = log10(z + 1))  # Adding 1 to avoid log(0) issues
      # 
      # # Define a custom color scale with transition points
      # custom_colorscale <- list(
      #   c(0, "blue"),     # 0% of scale (minimum)
      #   c(0.25, "green"), # 33% of scale
      #   c(0.50, "yellow"),# 66% of scale
      #   c(1, "red")       # 100% of scale (maximum)
      # )
      # 
      # # Create the 3D plot
      # plot_ly(data = df_3d, x = ~x, y = ~y, z = ~z, type = "mesh3d",
      #         intensity = ~z, colorscale = custom_colorscale,
      #         showscale = TRUE) %>%
      #   layout(
      #     title = paste(protein, motif, "3D Histogram of fRead_", biosample_comb[1], " vs fRead_", biosample_comb[2]),
      #     scene = list(
      #       xaxis = list(title = paste("fRead_", biosample_comb[1])),
      #       yaxis = list(title = paste("fRead_", biosample_comb[2])),
      #       zaxis = list(title = "Log(Frequency)", type = "log")
      #     )
      #   )
      # 
      # # Ensure 'sample' is a factor or character in df_state1 before binning
      # df_state1 <- df_state1 %>%
      #   mutate(sample = as.factor(sample))
      # 
      # # Bin the data into a 2D histogram and count occurrences
      # df_bins_plot <- ggplot(df_state1, aes_string(x = fRead_col1, y = fRead_col2, fill = "sample")) +
      #   stat_bin2d(bins = 30)
      # 
      # # Plot using ggplot2 with custom colors
      # df_bins_plot <- ggplot(df_state1, aes_string(x = fRead_col1, y = fRead_col2, fill = "sample")) +
      #   stat_bin2d(bins = 30) +
      #   scale_fill_manual(values = group.colors, na.value = "grey50") +  # Use grey for undefined colors
      #   labs(
      #     title = "2D Density Plot with Custom Colors",
      #     x = paste("fRead_", biosample_comb[1]),
      #     y = paste("fRead_", biosample_comb[2])
      #   ) +
      #   theme_minimal()
      # 
      # # Extract binned data
      # df_bins <- ggplot_build(df_bins_plot)$data[[1]]
      # df_3d <- df_bins %>%
      #   select(x, y, count, fill) %>%
      #   rename(z = count, sample = fill) %>%
      #   mutate(sample = as.factor(sample))  # Ensure sample is a factor
      # 
      # # Log transformation for z values
      # df_3d <- df_3d %>%
      #   mutate(z = log10(z + 1))
      # 
      # # Create the 3D plot with color by sample
      # plot_ly(data = df_3d, x = ~x, y = ~y, z = ~z, type = "mesh3d",
      #         color = sample,
      #         showscale = TRUE) %>%
      #   layout(
      #     title = "3D Histogram of fRead for Multiple Samples",
      #     scene = list(
      #       xaxis = list(title = paste("fRead_", biosample_comb[1])),
      #       yaxis = list(title = paste("fRead_", biosample_comb[2])),
      #       zaxis = list(title = "Log(Frequency)", type = "log")
      #     )
      #   )
      
      
      
      
    }
  }
}


# Define the unique samples
unique_samples <- unique(df_state1$sample)

# Bin data and calculate counts for each sample
binned_data <- unique_samples %>%
  lapply(function(sample) {
    # Filter for each sample
    sample_data <- df_state1 %>% filter(sample == !!sample)
    
    # Create 2D bins and count occurrences
    binned <- ggplot(sample_data, aes(x = !!sym(fRead_col1), y = !!sym(fRead_col2))) +
      stat_bin2d(bins = 30) %>%
      ggplot_build()
    
    # Extract binned data as a data frame and add sample and z level
    as.data.frame(binned$data[[1]]) %>%
      mutate(sample = sample) %>%
      rename(x = x, y = y, z = count)
  }) %>%
  bind_rows()

library(plotly)

# Define color and transparency for each sample
sample_colors <- c(
  "HepG2" = "rgba(248,118,109,0.5)",  # semi-transparent for HepG2
  "K562" = "rgba(0,191,196,0.5)",     # semi-transparent for K562
  "GM12878" = "rgba(163,165,0,0.5)",  # semi-transparent for GM12878
  "A549" = "rgba(231,107,243,0.5)"    # semi-transparent for A549
)

# Initialize plot
plot <- plot_ly()

# Add planes for each sample
for (sample in unique_samples) {
  # sample <- unique_samples[[1]]
  
  sample_data <- binned_data %>% filter(sample == !!sample)
  
  # Create a mesh3d plane with semi-transparent color for each sample
  plot <- plot %>%
    add_trace(
      data = sample_data,
      x = ~x,
      y = ~y,
      z = ~z,
      type = "mesh3d",
      opacity = 0.5,  # Set opacity for transparency
      color = I(sample_colors[sample]),
      showscale = FALSE
    )
}

# Layout adjustments
plot <- plot %>%
  layout(
    title = "3D Planes for fRead Counts by Sample",
    scene = list(
      xaxis = list(title = paste("fRead_", biosample_comb[1])),
      yaxis = list(title = paste("fRead_", biosample_comb[2])),
      zaxis = list(title = "Counts")
    )
  )

plot

######################################################




library(ggplot2)

# Specify the sample you want to plot
sample <- "A549"  # Replace with the desired sample name

# Filter df_state1 for the specified sample
df_sample <- df_state1 %>% filter(sample == !!sample)

# Define a custom color gradient for the fill scale
#custom_colors <- c("#f7fbff", "grey", group.colors[sample]) 


# Define colors for the gradient, with an intermediate color between light and the sample color
light_color <- "#f7fbff"
# Define the color for the sample, defaulting to grey if sample is not in group.colors
sample_color <- ifelse(sample %in% names(group.colors), group.colors[sample], "grey")
# Compute middle color
intermediate_color <- colorRampPalette(c(light_color, sample_color))(3)[2]  

# Define custom colors and values for gradient
custom_colors <- c(light_color, intermediate_color, sample_color)


custom_values <- c(0, 0.01, 1)  # Breakpoints for each color

# Create the 2D heatmap plot with custom color gradient
ggplot(df_sample, aes_string(x = fRead_col1, y = fRead_col2)) +
  stat_bin2d(bins = 30, aes(fill = ..count..)) +
  geom_point(color = sample_color, size = 1, alpha = 0.5)+
  labs(
    title = paste("2D Heatmap of Counts for Sample:", sample),
    x = paste("fRead_", biosample_comb[1]),
    y = paste("fRead_", biosample_comb[2])
  ) +
  scale_fill_gradientn(colors = custom_colors,
                       values = custom_values) +  # Apply custom color gradient
  theme_minimal()

####################################################################################


# Function to create a heatmap for a specific sample
create_heatmap <- function(sample_name) {
  # Filter for the specified sample
  df_sample <- df_state1 %>% filter(sample == sample_name)
  
  # Define colors for the gradient, defaulting to grey if sample is not in group.colors
  light_color <- "#f7fbff"
  sample_color <- ifelse(sample_name %in% names(group.colors), group.colors[sample_name], "grey")
  intermediate_color <- colorRampPalette(c(light_color, sample_color))(3)[2]
  
  # Define custom colors and values for gradient
  custom_colors <- c(light_color, intermediate_color, sample_color)
  custom_values <- c(0, 0.01, 1)
  
  # Create the 2D heatmap plot with custom color gradient and overlay all points
  ggplot(df_sample, aes_string(x = fRead_col1, y = fRead_col2)) +
    stat_bin2d(binwidth = c(0.025, 0.025), aes(fill = ..count..)) +  # Set bin width to 0.05 for both x and y axes
    geom_point(color = sample_color, size = 1, alpha = 0.5) +
    labs(
      title = paste("ChIP in ", sample_name),
      x = paste("fRead_", biosample_comb[1]),
      y = paste("fRead_", biosample_comb[2])
    ) +
    scale_fill_gradientn(colors = custom_colors, values = custom_values) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +  # Set axis limits without cutting off data
    theme_minimal()
}

# Create individual plots for each of the specified samples
plot1 <- create_heatmap(biosample_comb[1])
plot2 <- create_heatmap(paste(biosample_comb[1], biosample_comb[2], sep = "_"))
plot3 <- create_heatmap(biosample_comb[2])

# Arrange the plots in a single grid
grid.arrange(plot1, plot2, plot3, ncol = 3)

# Load necessary library for arranging
library(gridExtra)

# Arrange the plots in a single grid and save as an image
ggsave(
  filename = file.path(
    output_folder,
    protein,
    motif,
    paste(
      motif,
      unlist(biosample_comb)[1],
      unlist(biosample_comb)[2],
      "combined_heatmap.png",
      sep = "_"
    )
  ),  # File name and format
  plot = arrangeGrob(plot1, plot2, plot3, ncol = 3),  # Combine plots in one grid
  width = 15, height = 5  # Adjust width and height as needed
)

##########################################################



library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# Function to create a heatmap for a specific sample
create_heatmap <- function(sample_name) {
  # Filter for the specified sample
  df_sample <- df_state1 %>% filter(sample == sample_name)
  
  # Define colors for the gradient, defaulting to grey if sample is not in group.colors
  light_color <- "#f7fbff"
  sample_color <- ifelse(sample_name %in% names(group.colors), group.colors[sample_name], "grey")
  intermediate_color <- colorRampPalette(c(light_color, sample_color))(3)[2]
  
  # Define custom colors and values for gradient
  custom_colors <- c(light_color, intermediate_color, sample_color)
  custom_values <- c(0, 0.01, 1)
  
  # Create the 2D heatmap plot with custom color gradient and overlay all points
  ggplot(df_sample, aes_string(x = fRead_col1, y = fRead_col2)) +
    stat_bin2d(binwidth = c(0.025, 0.025), aes(fill = ..count..)) +  # Set bin width to ensure consistent binning
    geom_point(color = sample_color, size = 1, alpha = 0.5) +
    labs(
      title = paste("ChIP in", sample_name),
      x = paste("fRead_", biosample_comb[1]),
      y = paste("fRead_", biosample_comb[2])
    ) +
    scale_fill_gradientn(colors = custom_colors, values = custom_values) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +  # Set axis limits without cutting off data
    theme_minimal()
}

# Create individual plots for each of the specified samples
plot1 <- create_heatmap(biosample_comb[1])
plot2 <- create_heatmap(paste(biosample_comb[1], biosample_comb[2], sep = "_"))
plot3 <- create_heatmap(biosample_comb[2])

# Combine the plots into a single grid with a super title
combined_plot <- arrangeGrob(
  plot1, plot2, plot3,
  ncol = 3,
  top = textGrob(motif, gp = gpar(fontsize = 16, fontface = "bold"))
)

# Save the combined plot as an image
output_file <- file.path(
  output_folder,
  protein,
  motif,
  paste(motif, biosample_comb[1], biosample_comb[2], "combined_heatmap.png", sep = "_")
)

ggsave(
  filename = output_file,
  plot = combined_plot,  # Use combined plot with title
  width = 15, height = 5  # Adjust width and height as needed
)

################################################


