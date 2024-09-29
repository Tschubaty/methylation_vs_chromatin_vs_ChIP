#################################################################
##  Chip-seq vs methylation vs Chromatin
##
##  input: everything in HG38
##
##          Encode3\
##
##
##  output:   Encode3\
##
##  v_1 16.09.2024
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
library(foreach)
library(doParallel)
library(readr)
library(dplyr)
library(rlang) 
##################################### INPUT ########################################



################################## constants #####################################
start_script <- Sys.time()
# Set input and output directories
input_folder <- file.path(this.dir,"fimo_methylation")
output_folder <- file.path(this.dir,"simmulation")

RNA_seq_folder <- "RNA-seq"

sample_names <- c("K562", "HepG2","GM12878","A549")


# group.colors <-
#   c(HepG2 = "#F8766D",
#     K562 = "#00BFC4",
#     both = "#e4d00a")
#chip.states <- names(group.colors)

# set up parallel processing
n_cores <- detectCores() - 1
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
colnames_df <- c("Chromosome","Start","End","CpG_Island_Status","Strand",
                "Chromatin_State_A549","mRead_A549","nRead_A549","fRead_A549",
                "Chromatin_State_GM12878","mRead_GM12878","nRead_GM12878","fRead_GM12878",
                "Chromatin_State_HepG2","mRead_HepG2","nRead_HepG2","fRead_HepG2",
                "Chromatin_State_K562","mRead_K562","nRead_K562","fRead_K562",
                "peak_A549","peak_GM12878","peak_HepG2","peak_K562"
                )
##################################### function definitions ########################################


#################################### read INPUT ##################################################
# Function to extract and clean column names from the first row of a text file
extract_colnames <- function(file_path) {
  # Read the first row
  first_row <- readLines(file_path, n = 1)
  
  # Remove the leading "#" and replace double tabs with single tab
  cleaned_row <- gsub("#", "", first_row)          # Remove leading #
  cleaned_row <- gsub("\t\t", "\t", cleaned_row)   # Replace double tabs with single tab
  
  # Split the cleaned string into individual column names
  colnames <- unlist(strsplit(cleaned_row, "\t"))
  
  return(colnames)
}

# # Function to read .bed files with dynamic column names
# read_bed_file <- function(file_path) {
#   # Extract the correct column names from the first row
#   colnames <- extract_colnames(file_path)
#   
#   # Read the BED file, skipping the first line (header row)
#   df <- readr::read_delim(file_path, delim = "\t", skip = 1, col_names = colnames)
#   
#   return(df)
# }

# Function to read .bed files with dynamically extracted column names and correct column types
read_bed_file <- function(file_path) {
  # Extract the correct column names from the first row
  colnames <- extract_colnames(file_path)
  
  # Define column types: integer for all except the specified character columns
  col_types <- cols(
    Chromosome = col_character(),  # Chromosome as character
    Sample_Info = col_character(), # Sample info as character
    Strand = col_character(),      # Strand as character (e.g., "+" or "-")
    CpG_Island_Status = col_character(), # CpG Island status as character
    
    # Assume all other columns are integers, dynamically handle columns from the colnames
    .default = col_integer()       # All other columns are integers
  )
  
  # Read the BED file, skipping the first line (header row)
  df <- readr::read_delim(file_path, delim = "\t", skip = 1, col_names = colnames, col_types = col_types)
  
  return(df)
}

# Function to extract sample type from the filename
extract_sample_type <- function(file_name) {
  # Assuming file name format: PROTEIN_SAMPLE_methylation.chromatinstate.fimo.bed
  sample_type <- strsplit(basename(file_name), "_")[[1]][2]
  return(sample_type)
}

## Function to perform a full outer join and rename conflicting columns
full_outer_join_samples <- function(df_1, df_2, sample_1, sample_2) {
  
  # Rename conflicting columns in both dataframes before joining
  df_1 <- df_1 %>%
    rename(
      !!paste0("Chromatin_State_", sample_1) := Chromatin_State,
      !!paste0("Sample_Info_", sample_1) := Sample_Info,
      !!paste0("Score_", sample_1) := Score
    )
  
  df_2 <- df_2 %>%
    rename(
      !!paste0("Chromatin_State_", sample_2) := Chromatin_State,
      !!paste0("Sample_Info_", sample_2) := Sample_Info,
      !!paste0("Score_", sample_2) := Score
    )
  
  # Perform a full outer join on the specified columns, including CpG_Island_Status as an identifier
  df_joined <- full_join(
    df_1, df_2, 
    by = c("Chromosome", "Start", "End", "Strand", "CpG_Island_Status")
  )
  
  return(df_joined)
}

# List protein folders in the input directory
protein_folders <- list.dirs(input_folder, full.names = TRUE, recursive = FALSE)

# Loop through each protein folder
for (protein_folder in protein_folders[1:10]) {
  protein_name <- basename(protein_folder)  # Extract protein name
  bed_files <- list.files(protein_folder, pattern = "*_methylation.chromatinstate.fimo.bed", full.names = TRUE)
  print(protein_name)
  # Extract sample types for each file
  sample_types <- sapply(bed_files, extract_sample_type)
  
  # Generate all pairwise combinations of the sample types
  sample_pairs <- combn(sample_types, 2, simplify = FALSE)
  
  # # Process each .bed file
  # for (bed_file in bed_files) {
  #   print(basename(bed_file))
  #   #df <- read_bed_file(bed_file)
  #   
  #   # Perform further analysis on the dataframe `df`
  #   #print(head(df))
  # }
  
  for (pair in sample_pairs) {
    sample_1 <- pair[[1]]
    sample_2 <- pair[[2]]
    
    # Find the corresponding files for the two samples
    file_1 <- bed_files[grep(sample_1, bed_files)]
    file_2 <- bed_files[grep(sample_2, bed_files)]
    
    print(sprintf("Comparing samples: %s vs %s", sample_1, sample_2))
    
    # Read the data for both samples
    df_1 <- read_bed_file(file_1)
    df_2 <- read_bed_file(file_2)
    
    df_merged <- full_outer_join_samples(df_1, df_2, sample_1, sample_2)
    
    # Print the first few rows of the merged dataframe to verify
    print(head(df_merged))
    
  }
}



# df_rna <-
#   read.csv(file = file.path(RNA_seq_folder, "rna_seq_of_chip_sep_tartgets_summery.csv"))
# files <-
#   list.files(path = input_folder, pattern = "*.ChiP_vs_WGBS.*")
# 
# files <- file.path(input_folder , files)


# #debug files <- list.files(path = input_folder, pattern = "*random*")
# files[length(files) + 1] <-
#   "monte_carlo_results/data/script4.random.peaks_200bp_apart.RDS"
# 
# #debug
# files <- files[c(
#   which(grepl(x = files, pattern = "AGO1"))[1],
#   which(grepl(x = files, pattern = "CEBPB"))[1],
#   which(grepl(x = files, pattern = "IKZF1"))[1],
#   which(grepl(x = files, pattern = "MAFF"))[1],
#   which(grepl(x = files, pattern = "MAFK"))[1],
#   which(grepl(x = files, pattern = "NBN"))[1],
#   which(grepl(x = files, pattern = "PCBP2"))[1],
#   which(grepl(x = files, pattern = "PRPF4"))[1]
# )]

# files <- files[c(
#   which(grepl(x = files, pattern = "ZBTB33"))[1],
#   which(grepl(x = files, pattern = "AGO1"))[1],
#   which(grepl(x = files, pattern = "MAFF"))[1],
#   which(grepl(x = files, pattern = "PCBP2"))[1],
#   which(grepl(x = files, pattern = "CEBPB"))[1],
#   which(grepl(x = files, pattern = "MAFK"))[1],
#   which(grepl(x = files, pattern = "PRPF4"))[1],
#   which(grepl(x = files, pattern = "JUN"))[1],
#   which(grepl(x = files, pattern = "ATF1"))[1],
#   which(grepl(x = files, pattern = "YY1"))[1],
#   which(grepl(x = files, pattern = "CTCF"))[1]
# )]



processed_results_summery <- data.frame()


#n_cores <- 2
# path_file <- file.path(output_folder,sprintf("output_log start %s.txt",start_script))
# cl <- makeCluster(n_cores,outfile=path_file)
registerDoParallel(n_cores)

processed_results_summery <-
  foreach (f = 1:length(files), .combine = rbind) %dopar% {
    #for (f in 1:nrow(files)) {
    # f <- which(grepl(x = files,pattern = "ZBTB33"))[1]
    #f <- which(grepl(x = files,pattern = "GPBP1L1"))
    # f <- 1
    # f <- which(grepl(x = files,pattern = "ZZZ3"))[1]
    
    ## START execution
    start_time <- Sys.time()
    file_name <- files[f]
    protein_name <-
      gsub(pattern = "4. chromatin_annotated_v03/|.ChiP_vs_WGBS_vs_chromatin.RDS|monte_carlo_results/data/script4.|.peaks_200bp_apart.RDS|",
           replacement = "",
           x = file_name)
    
    # create folder
    target_folder <- file.path(output_folder, "targets")
    dir.create(file.path(output_folder, "targets"), showWarnings = FALSE)
    # create protein folder
    output_folder_protein <-
      file.path(target_folder, protein_name)
    dir.create(output_folder_protein, showWarnings = FALSE)
    # create plot folder
    output_folder_protein_plots <-
      file.path(output_folder_protein, "plots")
    dir.create(output_folder_protein_plots, showWarnings = FALSE)
    # create data folder
    output_folder_protein_data <-
      file.path(output_folder_protein, "data")
    dir.create(output_folder_protein_data, showWarnings = FALSE)
    
    summery_file_name <- file.path(output_folder_protein_data,
                                   paste("summery",
                                         protein_name,
                                         "RDS",
                                         sep = "."))
    
    if (file.exists(summery_file_name) & cluster) {
      df_temp <- readRDS(summery_file_name)
    } else{
      # read input
      df <- readRDS(file = file_name)
      # allocate output
      df_temp <- data.frame()
      
      ############# general plots #################
      
      
      # look only at data with same chromatin state
      df <-
        df[df$chromatin_state_HepG2 == df$chromatin_state_K562,]
      df$chromatin_state <-
        df$chromatin_state_HepG2
      df <-
        df[, !(names(df) %in% c("chromatin_state_HepG2" , "chromatin_state_K562"))]
      
      # look only at data with methylation values
      df <-
        df[!is.na(df$fRead_K562) & !is.na(df$fRead_HepG2) , ]
      #  restrict to unique peaks
      merge_colnames <-
        c(
          "chr",
          "merge_peak",
          "ChiP_state",
          "mRead_K562",
          "mRead_HepG2",
          "fRead_K562",
          "fRead_HepG2",
          "mean_n_K562",
          "mean_n_HepG2" ,
          "n_CPG" ,
          "chromatin_state",
          "hits",
          "name",
          "nRead_K562",
          "nRead_HepG2"
        )
      
      df <- df[!duplicated(df[, merge_colnames]), merge_colnames]
      
      
      # start plotting
      
      # ven diagram
      
      
      if (!cluster) {
        # if no data or random data
        if (protein_name == "random") {
          plot_ven <- ggplot2::ggplot()
        } else{
          ven_values <-  c(
            "HepG2" = sum(df$ChiP_state == "HepG2"),
            "K562" = sum(df$ChiP_state == "K562"),
            "K562&HepG2" = sum(df$ChiP_state == "both")
          )
          plot_ven <- plot(
            eulerr::euler(ven_values),
            fills = list(fill = group.colors),
            legend = list(side = "right"),
            quantities = list(cex = 2),
            labels = c()
          )
          
          ggplot2::ggsave(
            filename = paste(
              protein_name,
              "Venn_Diagramm_peaks",
              picuture_file_extension,
              sep = "."
            ),
            plot = plot_ven,
            device = picuture_file_extension,
            path =  output_folder_protein_plots,
            width = 1920,
            height = 1080,
            units = "px"
          )
        }
        
        # plot bio replicates
        plot_bio_rep <-
          ggplot2::ggplot(data = df,
                          mapping = ggplot2::aes(x = hits, fill = ChiP_state)) +
          ggplot2::geom_bar() +
          ggplot2::ggtitle(paste("peak classification according to sample: ", protein_name)) +
          ggplot2::scale_fill_manual(name = "ChiP_state",
                                     values = group.colors) +
          ggplot2::xlab("# replicated") +
          ggplot2::theme(
            legend.key = ggplot2::element_rect(fill = NA),
            text = ggplot2::element_text(size = 15),
            axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
          )
        
        ggplot2::ggsave(
          filename = paste(
            protein_name,
            "peak classification according to sample",
            picuture_file_extension,
            sep = "."
          ),
          plot = plot_bio_rep,
          device = picuture_file_extension,
          path =  output_folder_protein_plots,
          width = 1920,
          height = 1080,
          units = "px"
        )
        
        # extract all values
        rna_p <- df_rna[df_rna$gene_name == protein_name, ]
        #check if rna data present
        is_rna_data_present  <- nrow(rna_p) > 0
        if (is_rna_data_present) {
          ex <-
            unique(gsub(
              x =  grep(
                x = colnames(df_rna),
                pattern = "*._TPM",
                value = TRUE
              ),
              pattern = "_TPM.*",
              replacement = ""
            ))
          
          TPM <- as.vector(t(rna_p[paste(ex, "TPM", sep = "_")]))
          names(TPM) <- NULL
          TPM_ci_lower_bound <-
            as.vector(t(rna_p[paste(ex, "TPM_ci_lower_bound", sep = "_")]))
          names(TPM_ci_lower_bound) <- NULL
          TPM_ci_upper_bound <-
            as.vector(t(rna_p[paste(ex, "TPM_ci_upper_bound", sep = "_")]))
          names(TPM_ci_upper_bound) <- NULL
          FPKM = as.vector(t(rna_p[paste(ex, "FPKM", sep = "_")]))
          names(FPKM) <- NULL
          FPKM_ci_lower_bound <-
            as.vector(t(rna_p[paste(ex, "FPKM_ci_lower_bound", sep = "_")]))
          names(FPKM_ci_lower_bound) <- NULL
          FPKM_ci_upper_bound <-
            as.vector(t(rna_p[paste(ex, "FPKM_ci_upper_bound", sep = "_")]))
          names(FPKM_ci_upper_bound) <- NULL
          
          # build data_frame
          df_rna_p <- data.frame(
            experiment = ex,
            TPM = TPM,
            TPM_ci_lower_bound = TPM_ci_lower_bound,
            TPM_ci_upper_bound = TPM_ci_upper_bound,
            FPKM = FPKM,
            FPKM_ci_lower_bound = FPKM_ci_lower_bound,
            FPKM_ci_upper_bound = FPKM_ci_upper_bound
          )
          # plot_TPM <-
          #   ggplot2::ggplot(
          #     data = df_rna_p,
          #     mapping = ggplot2::aes(x = experiment, y = TPM, fill = experiment)
          #   ) +
          #   ggplot2::geom_bar(stat = "identity") +
          #   ggplot2::geom_errorbar(
          #     mapping = ggplot2::aes(ymin = TPM_ci_lower_bound, ymax = TPM_ci_upper_bound),
          #     width = 0.2,
          #     size = 1,
          #     color = "blue"
          #   ) +
          #   ggplot2::ggtitle(paste(protein_name, "expression in TPM"))
          
          
          plot_FPKM <-
            ggplot2::ggplot(
              data = df_rna_p,
              mapping = ggplot2::aes(
                x = experiment,
                y = FPKM,
                fill = experiment
              )
            ) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::geom_errorbar(
              mapping = ggplot2::aes(ymin = FPKM_ci_lower_bound, ymax = FPKM_ci_upper_bound),
              width = 0.2,
              size = 1,
              color = "black"
            ) +
            ggplot2::ggtitle(paste(protein_name, "expression in FPKM")) +
            ggplot2::theme(
              legend.key = ggplot2::element_rect(fill = NA),
              text = ggplot2::element_text(size = 15),
              axis.line = ggplot2::element_line(colour = "black"),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.background = ggplot2::element_blank(),
            )
        } else{
          plot_TPM <- ggplot2::ggplot()
          plot_FPKM <- ggplot2::ggplot()
          
        }
        
        ggplot2::ggsave(
          filename = paste(
            protein_name,
            "exression in FPKM",
            picuture_file_extension,
            sep = "."
          ),
          plot = plot_FPKM,
          device = picuture_file_extension,
          path =  output_folder_protein_plots,
          width = 1920,
          height = 1080,
          units = "px"
        )
      }
      ##############################m START chromatin speciifc ########################
      
      plot_list_chromatin_state_specific <- list()
      
      for (state_counter in 1:length(chromatin_state_names)) {
        #dubug  state_counter <- 1
        
        c_state <- chromatin_state_names[state_counter]
        #print(c_state)
        # only look at one chromatine state
        df_state <-
          df[!is.na(df$chromatin_state == c_state) &
               df$chromatin_state == c_state, ]
        
        # if no data or random data
        if (nrow(df_state) < 1 | protein_name == "random") {
          empty_plot <- ggplot2::ggplot()
          
          plot_list_chromatin_state_specific <-
            append(plot_list_chromatin_state_specific,
                   list(empty_plot,
                        empty_plot))
          
          df_temp <- rbind(
            df_temp,
            data.frame(
              protein_name = protein_name,
              chromatin_state = c_state,
              n_state = nrow(df_state),
              p_love = NA,
              p_hate = NA,
              statistic  = NaN,
              n_Chip_both = 0,
              n_Chip_HepG2 = 0,
              n_Chip_K562 = 0
            )
          )
          next
        }
        
        sim_number <- 100000
        
        meth_statistic_vec <-
          function(ChiP_state,
                   fRead_HepG2,
                   fRead_K562) {
            mean_K562_statistic_value <-
              mean(fRead_HepG2[which(ChiP_state == "K562")] -  fRead_K562[which(ChiP_state == "K562")], na.rm = TRUE)
            mean_HepG2_statistic_value <-
              mean(fRead_HepG2[which(ChiP_state == "HepG2")] - fRead_K562[which(ChiP_state == "HepG2")], na.rm = TRUE)
            statistic <-
              mean_HepG2_statistic_value  - mean_K562_statistic_value
            return(statistic)
          }
        
        #real_stat <-  meth_statistic(df_state)
        real_stat <-
          meth_statistic_vec(df_state$ChiP_state,
                             df_state$fRead_HepG2,
                             df_state$fRead_K562)
        
        
        file_name_mote_carlo <-
          file.path(
            output_folder_protein_data,
            paste(protein_name,
                  c_state,
                  "monte_carlo_values",
                  "RDS", sep = ".")
          )
        
        if (file.exists(file_name_mote_carlo)) {
          monte_carlo_values <- readRDS(file = file_name_mote_carlo)
        } else{
          fRead_HepG2 <- df_state$fRead_HepG2
          fRead_K562 <- df_state$fRead_K562
          
          monte_carlo_values <- sapply(
            1:sim_number,
            FUN =  function(x) {
              return(
                meth_statistic_vec(
                  ChiP_state = sample(df_state$ChiP_state),
                  fRead_HepG2 = fRead_HepG2,
                  fRead_K562 = fRead_K562
                )
              )
            }
          )
          
          saveRDS(object = monte_carlo_values,
                  file = file.path(
                    output_folder_protein_data,
                    paste(
                      protein_name,
                      c_state,
                      "monte_carlo_values",
                      "RDS",
                      sep = "."
                    )
                  ))
        }
        
        p_love = sum(real_stat < monte_carlo_values) / length(monte_carlo_values)
        p_hate = sum(real_stat > monte_carlo_values) / length(monte_carlo_values)
        
        if (!cluster) {
          plot_meth_statistic <-
            ggplot2::ggplot(data = data.frame(statistic  = monte_carlo_values)) +
            ggplot2::geom_histogram(mapping = ggplot2::aes(x = statistic),
                                    bins = 1000) +
            ggplot2::geom_vline(xintercept = real_stat, color = "red") +
            # ggplot2::ggtitle(c_state)+
            #   paste(
            #     protein_name,
            #     ": " ,
            #     c_state,
            #     sim_number,
            #     "perm. p_affin:",
            #     format.pval(p_love, eps = 1 / sim_number),
            #     "  p_avers:",
            #     format.pval(p_hate, eps = 1 / sim_number)
            #   )
          # ) +
          ggplot2::xlab("methylation statistic") +
            ggplot2::theme(
              legend.key = ggplot2::element_rect(fill = NA),
              text = ggplot2::element_text(size = 10),
              axis.line = ggplot2::element_line(colour = "black"),
              axis.text = ggplot2::element_text(size = 15),
              axis.title  = ggplot2::element_text(size = 20),
              panel.background = ggplot2::element_blank(),
            )
          
          ggplot2::ggsave(
            filename = paste(
              protein_name,
              c_state,
              "permutation test",
              picuture_file_extension,
              sep = "."
            ),
            plot = plot_meth_statistic,
            device = picuture_file_extension,
            path =  output_folder_protein_plots,
            width = 1920,
            height = 1080,
            units = "px"
          )
          
          
          plot_fReads <- ggplot2::ggplot(data = df_state) +
            ggplot2::geom_point(mapping =
                                  ggplot2::aes(
                                    x = fRead_HepG2,
                                    y = fRead_K562,
                                    color = ChiP_state
                                  )) +
            ggplot2::scale_color_manual(name = "ChiP_state",
                                        values = group.colors) +
            ggplot2::ggtitle(c_state)+
            # ggplot2::ggtitle(label = paste(
            #   protein_name,
            #   c_state ,
            #   "HepG2:",
            #   sum(df_state$ChiP_state == "HepG2"),
            #   "K562:",
            #   sum(df_state$ChiP_state == "K562"),
            #   "both:",
            #   sum(df_state$ChiP_state == "HepG2")
            # )) +
            ggplot2::scale_x_continuous(limits = c(0, 1),
                                        name = expression('methylation_value'["Hepg2"]),) +
            ggplot2::scale_y_continuous(limits = c(0, 1),
                                        name = expression('methylation_value'["K562"]),) +
            ggplot2::theme(
              legend.key = ggplot2::element_rect(fill = NA),
              text = ggplot2::element_text(size = 12),
              axis.text.x = ggplot2::element_text(face = "bold"),
              #color = group.colors["Hepg2"]), #, size = 22),
              axis.title.y = ggplot2::element_text(face = "bold"),
              #color = group.colors["K562"]), # size = 32,
              axis.line = ggplot2::element_line(colour = "black"),
              panel.background = ggplot2::element_blank(),
              legend.position = "none"
            )
          
          ggplot2::ggsave(
            filename = paste(
              protein_name,
              c_state,
              "methylation values",
              picuture_file_extension,
              sep = "."
            ),
            plot = plot_fReads,
            device = picuture_file_extension,
            path =  output_folder_protein_plots,
            width = 1920,
            height = 1080,
            units = "px"
          )
          
          
          # plot_density_methylation_values <-
          #   ggplot2::ggplot(data = df_state,
          #                   mapping = ggplot2::aes(x = fRead_HepG2,
          #                                          y = fRead_K562)) +
          #   stat_density_2d(geom =  "polygon",#  "raster", # , "point"
          #                   aes(alpha = (..level..) ^ 2, fill = ChiP_state))+
          #   ggplot2::geom_point(mapping =
          #                         ggplot2::aes(x = fRead_HepG2,
          #                                      y = fRead_K562,
          #                                      color = ChiP_state))
          #
          #
          #   ggplot2::ggplot(data = df_state,
          #                   mapping = ggplot2::aes(x = fRead_HepG2,
          #                                          y = fRead_K562)) +
          
          #                   sum(df_state$ChiP_state == "K562"),
          #                   "both:",
          #                   sum(df_state$ChiP_state == "HepG2")
          #     ))+
          #   ggplot2::geom_point( mapping = ggplot2::aes(color = ChiP_state))+
          #   ggplot2::stat_density_2d(geom = "polygon",
          #   bins = 10,
          #   mapping = ggplot2::aes(alpha = ..level..,
          #                          fill = ChiP_state,
          #                          group = ChiP_state))
          #
          #
          #   # ggplot2::scale_color_manual(name = "ChiP_state",
          #   #                             values = group.colors)+
          #
          # set.seed(123)
          # plot_data <-
          #   data.frame(
          #     X = c(rnorm(300, 3, 2.5), rnorm(150, 7, 2)),
          #     Y = c(rnorm(300, 6, 2.5), rnorm(150, 2, 2)),
          #     Label = c(rep('A', 300), rep('B', 150))
          #   )
          # ggplot(data = df_state, aes(x = fRead_HepG2,
          #                             y = fRead_K562, group = ChiP_state)) +
          #   stat_density_2d(geom = "polygon",
          #                   aes(alpha = (..level..) ^ 2, fill = ChiP_state),
          #                   bins = 10)+
          #   geom_point( aes(color = ChiP_state))+
          #   xlim(0, 0.02)+
          #   ylim(0, 0.02)
          #
          #
          #
          #
          #
          #     ggplot2::stat_density_2d(geom = "polygon",
          #                   aes(alpha = ..level.., fill = ChiP_state),
          #                   bins = 4,group = df_state$ChiP_state)
          #
          # +
          # ggplot2::annotate(
          #   "text",
          #   x = 0.1,
          #   y = 1.05,
          #   label = paste(sum(df_state$ChiP_state == "K562"), "K562 datapoints"),
          #   size = 4,
          #   color = group.colors["K562"]
          # ) +
          # ggplot2::annotate(
          #   "text",
          #   x = 0.5,
          #   y = 1.05,
          #   label = paste(sum(df_state$ChiP_state == "both"), "both datapoints"),
          #   size = 4,
          #   color = group.colors["both"]
          # ) +
          # ggplot2::annotate(
          #   "text",
          #   x = 0.9,
          #   y = 1.05,
          #   label = paste(sum(df_state$ChiP_state == "HepG2"), "HepG2 datapoints"),
          #   size = 4,
          #   color = group.colors["HepG2"]
          #)
          
          # plot_islands_chip <-
          #   ggplot2::ggplot(data = df_state,
          #                   mapping = ggplot2::aes(x = name, fill = ChiP_state)) +
          #   ggplot2::geom_bar() +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "CpG Island and Chip state"))
          #
          # plot_islands_hist <-
          #   ggplot2::ggplot(data = df_state,
          #                   mapping = ggplot2::aes(x = n_CPG , fill = name)) +
          #   ggplot2::geom_histogram(breaks = seq(from = 0, to = 100, by = 2)) +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "CpG Island and nCpG"))
          #
          # plot_CpG_quality_K562 <- ggplot2::ggplot(data = df_state) +
          #   ggplot2::geom_point(mapping = ggplot2::aes(x = n_CPG, y = mean_n_K562, color = ChiP_state)) +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "CpG_quality K562")) +
          #   ggplot2::ylim(0, 300)
          #
          #
          # plot_CpG_quality_HepG2 <- ggplot2::ggplot(data = df_state) +
          #   ggplot2::geom_point(mapping = ggplot2::aes(x = n_CPG, y = mean_n_HepG2, color = ChiP_state)) +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "CpG_quality HepG2")) +
          #   ggplot2::ylim(0, 300)
          #
          # plot_methylation_read_depth <-  ggplot2::ggplot(data = df_state) +
          #   ggplot2::geom_point(mapping = ggplot2::aes(x = fRead_HepG2, y = fRead_K562, color = n_CPG)) +
          #   #ggplot2::scale_colour_gradient(low = "yellow", high = "red", na.value = NA)+
          #   ggplot2::scale_colour_gradient2(
          #     #breaks = c(0,5,10,20),
          #     low = "red",
          #     mid = "white",
          #     high = "blue",
          #     midpoint = 25,
          #     space = "Lab",
          #     na.value = NA,
          #     guide = "colourbar",
          #     aesthetics = "colour"
          #   ) +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "methylation and read_depth"))
          #
          
          
          plot_list_chromatin_state_specific <-
            append(
              plot_list_chromatin_state_specific,
              list(plot_fReads, plot_meth_statistic)
            )
          
        } # end plotting
        
        
        df_temp <- rbind(
          df_temp,
          data.frame(
            protein_name = protein_name,
            chromatin_state = c_state,
            n_state = nrow(df_state),
            p_love = p_love,
            p_hate = p_hate,
            statistic  = real_stat,
            n_Chip_both = sum(df_state$ChiP_state == "both"),
            n_Chip_HepG2 = sum(df_state$ChiP_state == "HepG2"),
            n_Chip_K562 = sum(df_state$ChiP_state == "K562")
          )
        )
        
        saveRDS(object = df_state,
                file = file.path(
                  output_folder_protein_data,
                  paste(protein_name,
                        c_state,
                        "RDS",
                        sep = ".")
                ))
        
      } # end c_state loop
      
      saveRDS(object = df_temp,
              file = summery_file_name)
      
      # # lay <- rbind(
      # #   c(1, 2, 3, 4),
      # #   c(5, 6, 7, 8),
      # #   c(9, 9, 10, 11),
      # #   c(9, 9, 10, 11),
      # #   c(12, 12, 13, 14),
      # #   c(15, 15, 16, 17)
      # # )
      #
      # plot_list_general <- list(plot_ven, plot_bio_rep, plot_FPKM)
      #
      # # plot_bio_rep,
      # # plot_states,
      # # plot_islands_3corner,
      # # plot_islands_chip,
      # # plot_islands_hist,
      # #plot_CpG_quality_K562,
      # #plot_CpG_quality_HepG2,
      # #plot_methylation_read_depth,
      # #plot_methylation_state_read_depth
      # # plot_ven,
      # # plot_FPKM,
      #
      plot_list <- plot_list_chromatin_state_specific
      
      #append(plot_list_general,
      #       plot_list_chromatin_state_specific)
      
      final_plot <-
        gridExtra::grid.arrange(
          grobs = plot_list,
          shared_legend,
          top = protein_name, #paste(protein_name, Sys.Date()),
          ncol = 6
          # layout_matrix = lay
        )
      
      
      pic_file_name <-
        paste(protein_name,
              "overview",
              picuture_file_extension,
              sep = ".")
      
      ggplot2::ggsave(
        filename = pic_file_name,
        plot = final_plot,
        device = picuture_file_extension,
        path =  output_folder_protein_plots,
        width = 7680,
        height = 4320,
        units = "px"
      )
      
      print(paste("saved", pic_file_name))
    }
    print(sprintf(
      "%s Finished in %s on %s",
      protein_name,
      format(Sys.time() - start_time),
      sprintf("on %s.RDS",
              format(Sys.time(),
                     "%d-%b-%Y %H.%M"))
    ))
    return(df_temp)
    #processed_results_summery <- rbind(processed_results_summery,df_temp)
  }
#parallel::stopCluster(cl)
stopImplicitCluster()

saveRDS(object = processed_results_summery,
        file = file.path(
          this.dir,
          output_folder ,
          "meta",
          sprintf(
            "protein_sensetivity_summery %s.RDS",
            format(Sys.time(), "%d-%b-%Y %H.%M")
          )
        ))

print(sprintf(
  "%s Finished : %s %s",
  "Script",
  format(Sys.time() - start_script),
  sprintf(
    "protein_sensetivity_summery %s.RDS",
    format(Sys.time(),
           "%d-%b-%Y %H.%M")
  )
))