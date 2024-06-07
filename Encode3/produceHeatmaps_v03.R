# Documentation -----------------------------------------------------------
##
##
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
###########################################################################

bedColnames <- c(
"chr",
"start",
"end",
"name",
"score",
"strand",
"StartThick",
"EndThick",
"color",
"n_read",
"f_read",
"Reference_genotype",
"Sample_genotype",
"Quality_score")

# Get unique chromosomes
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

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(openxlsx)
library(grid)
library(cowplot)

file_path <- file.path("Chromatin","18State_model.xlsx")

# Read the Excel file
state_model <- read.xlsx(file_path, sheet = 1, colNames = TRUE)
state_model$STATE.NO. <- factor(x = ,state_model$STATE.NO.,levels = sort(unique(state_model$STATE.NO.)))

# Function to read and combine all chromosome files for a specific biosample
combine_chromosome_files <- function(biosample, input_dir) {
  combined_data <- do.call(rbind, lapply(CHR_NAMES, function(chr) {
    # chr <- CHR_NAMES[22]
    file <-
      file.path(input_dir,
                paste0(
                  chr,
                  "_cpgs_island_",
                  biosample,
                  "_chromatin_histones.bed"
                ))
    if (file.exists(file)) {
      print(file)
      # Read the header separately
      header <- readLines(file, n = 1)
      col_names <- strsplit(header, " ")[[1]]
      col_names[[1]] <- "chr"
      
      data <-
        read.table(
          file,
          header = FALSE,
          skip = 0,
          sep = "\t",
          stringsAsFactors = FALSE
        )
      colnames(data) <- col_names
      return(data)
    } else {
      stop("File not found: ", file)
      return(NULL)
    }
  }))
  
  if (is.null(combined_data)) {
    stop("No files found for the biosample: ", biosample)
  }
  
  return(combined_data)
}

biosamples <- c("K562", "GM12878", "HepG2", "A549")
for (biosample in biosamples) {
  # Main script
  #biosample <- "K562"  # Change this to the desired biosample
  input_dir <- "WGBS/byChr/histone_annotated"
  output_file <- paste0(biosample, "_chromatin_histone_heatmap.png")
  
 
  combined_fime_name <-
    file.path("Histones",
              paste0("chrALL_", biosample, "_chromatin_histones.rds"))
  if (file.exists(combined_fime_name)) {
    combined_data <- readRDS(combined_fime_name)
  } else{
    # Combine the data from all chromosome files
    combined_data <- combine_chromosome_files(biosample, input_dir)
    saveRDS(file = combined_fime_name, object = combined_data)
    print("combination complete")
  }
  
  colnames_state_summery <-
    c("biosample", names(combined_data)[6:length(names(combined_data))], "n_CpG")
  state_summery <-
    data.frame(matrix(
      ncol = length(colnames_state_summery),
      nrow = max(combined_data$Chromatin_State)
    ))
  colnames(state_summery) <- colnames_state_summery
  state_summery$biosample <- biosample
  
  
  for (state in 1:18) {
    print(state)
    state_data <-
      combined_data[combined_data$Chromatin_State == state,]
    state_summery$Chromatin_State[state] <- state
    state_summery$n_CpG[state] <- nrow(state_data)
    
    for (hist in names(combined_data)[7:length(names(combined_data))]) {
      print(hist)
      state_summery[state, hist] <- mean(state_data[, hist])
      
    }
  }
  state_summery$Chromatin_State <- factor(x = state_summery$Chromatin_State,levels = sort(unique(state_summery$Chromatin_State)))
  saveRDS(file = file.path("Histones", paste0(biosample, "_histone_summery.rds")), object = state_summery)
  print("saved  state_summery")
  

  
  state_summery_long <-
    melt(
      state_summery,
      id.vars = c("biosample", "Chromatin_State", "n_CpG"),
      variable.name = "Histone_Mark",
      value.name = "Occurrence_Probability"
    )
  
  state_summery_long$Histone_Mark <-
    factor(state_summery_long$Histone_Mark, levels = sort(as.character(
      unique(state_summery_long$Histone_Mark)
    )))
  
  heatmap_plot <-
    ggplot(
      state_summery_long,
      aes(x = Histone_Mark, y = Chromatin_State, fill = Occurrence_Probability)
    ) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(
      title = paste(
        "Histone Mark Presence Probability in Chromatin States for",
        biosample
      ),
      fill = "Probability"
    )
  
  ggsave(filename = file.path("Histones","plots",paste0(biosample, "_heatmap.png")),plot =  heatmap_plot)
  
}

# Merge state_summery_long with state_model to get DESCRIPTION and COLOR.CODE
merged_data <- state_summery_long %>%
  left_join(state_model, by = c("Chromatin_State" = "STATE.NO."))

# Convert COLOR.CODE to a usable color format
merged_data$COLOR.CODE <- sapply(merged_data$COLOR.CODE, function(x) {
  rgb_values <- as.numeric(unlist(strsplit(x, ",")))
  rgb(rgb_values[1], rgb_values[2], rgb_values[3], maxColorValue = 255)
})

# Create a custom theme for the y-axis labels
y_axis_labels <- merged_data$DESCRIPTION[1:18]
y_axis_colors <- merged_data$COLOR.CODE[1:18]
y_axis_df <- data.frame(
  DESCRIPTION = factor(y_axis_labels, levels = y_axis_labels),
  COLOR.CODE = y_axis_colors
)

heatmap_plot <- ggplot(merged_data, aes(x = Histone_Mark, y = DESCRIPTION, fill = Occurrence_Probability)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = paste("Histone Mark Presence Probability in Chromatin States for", biosample),
    fill = "Probability"
  )

#Function to create colored background for y-axis labels
background_grob <- function(y_axis_labels, y_axis_colors) {
  grobs <- list()
  for (i in seq_along(y_axis_labels)) {
    grobs[[i]] <- annotation_custom(
      grob = rectGrob(gp = gpar(fill = y_axis_colors[i], col = NA)),
      ymin = i - 0.5, ymax = i + 0.5, xmin = -2, xmax = 0.1
    )
  }
  grobs
}

# Add background rectangles for y-axis labels
for (bg in background_grob(y_axis_labels, y_axis_colors)) {
  heatmap_plot <- heatmap_plot + bg
}

# Print the plot
print(heatmap_plot)


# Function to create custom y-axis labels with background colors
custom_y_axis_labels <- function(labels, colors) {
  grobs <- list()
  for (i in seq_along(labels)) {
    grobs[[i]] <- gTree(children = gList(
      rectGrob(gp = gpar(fill = colors[i], col = NA)),
      textGrob(labels[i], x = unit(0.5, "npc"), gp = gpar(col = "black", fontsize = 10, fontface = "bold"))
    ), vp = viewport(y = unit(i, "npc") - unit(0.5, "npc"), height = unit(1, "npc") / length(labels)))
  }
  grobTree(children = do.call(gList, grobs))
}

# Add the custom y-axis labels with backgrounds to the plot
custom_labels_grob <- custom_y_axis_labels(y_axis_labels, y_axis_colors)
heatmap_grob <- ggplotGrob(heatmap_plot)

# Create a new plot with the custom y-axis labels
final_plot <- ggdraw() +
  draw_grob(custom_labels_grob, x = 0.05, y = 0.5, width = 0.1, height = 1) +
  draw_grob(heatmap_grob, x = 0.5, y = 0.5, width = 0.9, height = 1)

# Print the final plot
print(final_plot)