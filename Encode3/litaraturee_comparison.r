# Clear R working environment
rm(list = ls())
cluster <- FALSE

# Set working directory
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"
  picuture_file_extension <- "pdf"
} else {
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picuture_file_extension <- "svg"
}
setwd(this.dir)




library(readxl)

kaluscha_file <- file.path(
  this.dir,
  "literature",
  "41588_2022_1241_MOESM3_ESM.xlsx"
)

download.file(
  url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01241-6/MediaObjects/41588_2022_1241_MOESM3_ESM.xlsx",
  destfile = kaluscha_file,
  mode = "wb"
)

excel_sheets(kaluscha_file)

# Then inspect the antibody sheets/tables
lapply(excel_sheets(kaluscha_file), function(sheet) {
  message("SHEET: ", sheet)
  print(read_excel(kaluscha_file, sheet = sheet, n_max = 10))
})

# -------------------------------------------------------------------------
# Kaluscha et al. 2022 / Nature Genetics
# DOI: 10.1038/s41588-022-01241-6
#
# Extended Data Fig. 5b visual extraction.
# These are candidate enriched TF motifs in DNMT-TKO-specific ATAC-seq peaks.
# They are NOT all individually validated methylation-sensitive proteins.
# -------------------------------------------------------------------------

kaluscha2022_fig5b_motif_labels <- c(
  "EGR3",
  "NRF1",
  "ELK1",
  "ETS1",
  "ERG",
  "PAX1",
  "GMEB1",
  "ATF4",
  "JUN",
  "CREB1",
  "JUNB",
  "JUND",
  "FOSL2",
  "FOSB",
  "FOSL1",
  "CREB5",
  "ATF7",
  "ATF1",
  "CREM",
  "CREB3",
  "MLX",
  "CTCF",
  "E2F2",
  "YY1",
  "YY2",
  "NEUROG2",
  "TGIF2",
  "MYF6",
  "SRY",
  "FOXF2",
  "FOXP1",
  "FOXP2",
  "FOXK1",
  "FOXD1",
  "FOXO3",
  "FOXP3",
  "FOXG1",
  "FOXJ2",
  "ONECUT1",
  "ONECUT2",
  "CUX1",
  "CUX2",
  "DUX",
  "NFYB",
  "NFYA"
)

length(kaluscha2022_fig5b_motif_labels)
# -------------------------------------------------------------------------
# Short protein-level call table for appending to the big comparison table
# -------------------------------------------------------------------------

kaluscha2022_fig5b_results_raw <- tibble::tibble(
  protein = kaluscha2022_fig5b_motif_labels
) %>%
  dplyr::mutate(
    protein = stringr::str_to_upper(stringr::str_trim(protein)),
    
    Kaluscha2022_Fig5b_call =
      "candidate_enriched_motif_in_DNMT_TKO_specific_ATAC_peaks") %>%
  dplyr::distinct(protein, .keep_all = TRUE) %>%
  dplyr::arrange(protein)

kaluscha2022_fig5b_results_raw



library(readxl)

yin2017_file <- file.path(
  this.dir,
  "literature",
  "aaj2239_yin_sm_tables_s1-s6.xlsx"
)

yin2017_df <- read_excel(
  path = yin2017_file,
  sheet = 3,
  skip = 20,   # row 21 becomes the column-name row; data starts at row 22
  na = c("", "NA")
)

yin2017_calls <- yin2017_df %>%
  dplyr::transmute(
    protein = stringr::str_to_upper(stringr::str_trim(`TF name`)),
    Yin2017_call = stringr::str_trim(Call)
  ) %>%
  dplyr::filter(
    !is.na(protein),
    !is.na(Yin2017_call),
    Yin2017_call != ""
  ) %>%
  dplyr::distinct()

yin2017_calls

jams2022_file <- file.path(
  this.dir,
  "literature",
  "13059_2022_2713_MOESM4_ESM.xlsx"
)

jams2022_df <- readxl::read_excel(
  path = jams2022_file,
  sheet = 1,
  skip = 2,   # row 3 contains column names
  na = c("", "NA")
)

jams2022_calls <- jams2022_df %>%
  dplyr::transmute(
    protein = TF %>%
      stringr::str_trim() %>%
      stringr::str_remove("_HUMAN$"),
    call = stringr::str_trim(Call),
    cell_line = dplyr::coalesce(
      stringr::str_trim(`Cell line`),
      "unknown_cell_line"
    )
  ) %>%
  dplyr::filter(
    !is.na(protein),
    !is.na(call),
    call != ""
  ) %>%
  dplyr::distinct() %>%
  dplyr::group_by(protein) %>%
  dplyr::summarise(
    JAMS2022_call = {
      unique_calls <- sort(unique(call))
      
      if (length(unique_calls) == 1) {
        unique_calls
      } else {
        purrr::map_chr(unique_calls, function(current_call) {
          call_cell_lines <- sort(unique(
            cell_line[call == current_call]
          ))
          
          paste0(
            current_call,
            "__",
            paste(call_cell_lines, collapse = "_")
          )
        }) %>%
          paste(collapse = "__")
      }
    },
    .groups = "drop"
  )

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

miodownik_file <- file.path(
  this.dir,
  "literature",
  "pnas.2520814122.sd01.xlsx"
)

miodownik2025_sites <- read_excel(
  path = miodownik_file,
  sheet = "S1B",
  skip = 1
) %>%
  fill(protein, .direction = "down") %>%
  mutate(
    protein = str_to_upper(protein),
    direction = str_trim(`meth. Effect`),
    position_change = str_trim(`position/s`),
    reference_motif = `reference DNA`
  ) %>%
  filter(!is.na(direction)) %>%
  select(
    protein,
    direction,
    position_change,
    reference_motif,
    comments
  )

miodownik2025_simple <- miodownik2025_sites %>%
  dplyr::transmute(
    protein = stringr::str_to_upper(stringr::str_trim(protein)),
    direction = stringr::str_to_lower(stringr::str_trim(direction))
  ) %>%
  dplyr::filter(
    !is.na(protein),
    direction %in% c("positive", "negative")
  ) %>%
  dplyr::group_by(protein) %>%
  dplyr::summarise(
    Miodownik2025_methylation_direction = paste0(
      "positive_", sum(direction == "positive"),
      "__negative_", sum(direction == "negative")
    ),
    .groups = "drop"
  )

miodownik2025_simple


library(readr)
library(dplyr)
library(stringr)

luo2021_df <- read_csv(
  file.path(
    this.dir,
    "literature",
    "luo2021_figure1B_color_labels.csv"
  ),
  show_col_types = FALSE
)

luo2021_simple <- luo2021_df %>%
  dplyr::transmute(
    protein = stringr::str_to_upper(stringr::str_trim(TF)),
    
    Luo2021_peak_methylation_label = dplyr::recode(
      Luo2021_color_label,
      "high_methylation" = "high_methylation_peaks",
      "low_methylation" = "low_methylation_peaks",
      "intermediate" = "intermediate_methylation_peaks",
      "conflicting" = "conflicting_across_cell_lines"
    )
  )

str(luo2021_simple)

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# ---------------- Hu et al. 2013: Supplementary file 1B ----------------
# All proteins/cofactors present on the protein microarray

hu2013_s1b_file <- file.path(
  this.dir,
  "literature",
  "Hu2013_S1B.xlsx"
)

hu2013_tested_proteins <- read_excel(
  hu2013_s1b_file,
  col_names = FALSE
) %>%
  unlist(use.names = FALSE) %>%
  as.character() %>%
  str_trim() %>%
  .[!is.na(.) & . != ""] %>%
  str_to_upper() %>%
  unique() %>%
  tibble(protein = .)

# ---------------- Hu et al. 2013: Supplementary file 1C ----------------
# Proteins with at least one methylated-CpG motif hit

hu2013_s1c_file <- file.path(
  this.dir,
  "literature",
  "Hu2013_S1C.xlsx"
)

hu2013_s1c <- read_excel(
  hu2013_s1c_file,
  skip = 1,
  na = c("", "NA")
) %>%
  transmute(
    protein = str_to_upper(str_trim(`TF hit`)),
    Hu2013_motif_ids_raw = `Binding motif(s)`,
    Hu2013_subfamily = Subfamily
  ) %>%
  filter(!is.na(protein))

# One row per protein–methylated-motif hit
hu2013_methylated_motif_hits <- hu2013_s1c %>%
  mutate(
    motif_id = str_extract_all(Hu2013_motif_ids_raw, "M\\d+")
  ) %>%
  select(protein, Hu2013_subfamily, motif_id) %>%
  unnest_longer(motif_id) %>%
  distinct(protein, motif_id, .keep_all = TRUE)

# One row per assayed protein/cofactor
hu2013_protein_summary <- hu2013_tested_proteins %>%
  left_join(
    hu2013_methylated_motif_hits %>%
      group_by(protein) %>%
      summarise(
        Hu2013_n_mCpG_motif_hits = n_distinct(motif_id),
        Hu2013_mCpG_motif_ids = paste(sort(unique(motif_id)), collapse = ";"),
        Hu2013_subfamily = first(Hu2013_subfamily),
        .groups = "drop"
      ),
    by = "protein"
  ) %>%
  mutate(
    Hu2013_n_mCpG_motif_hits = coalesce(Hu2013_n_mCpG_motif_hits, 0L),
    
    Hu2013_no_mCpG_motif_hit =
      Hu2013_n_mCpG_motif_hits == 0,
    
    Hu2013_screen_status = case_when(
      Hu2013_n_mCpG_motif_hits > 0 ~ "mCpG_preferential_hit",
      Hu2013_n_mCpG_motif_hits == 0 ~ "assayed_no_mCpG_hit"
    )
  ) %>%
  arrange(desc(Hu2013_n_mCpG_motif_hits), protein)

# Counts of proteins with versus without methylated-motif hits
hu2013_overall_counts <- hu2013_protein_summary %>%
  count(Hu2013_screen_status, name = "n_proteins")

hu2013_overall_counts

# Inspect strongest methylated-motif hit proteins
hu2013_protein_summary %>%
  filter(Hu2013_n_mCpG_motif_hits > 0) %>%
  select(
    protein,
    Hu2013_n_mCpG_motif_hits,
    Hu2013_mCpG_motif_ids,
    Hu2013_subfamily
  ) %>%
  print(n = 50)


# Read every protein name from Hu2013 Supplementary file 1B,
# ignoring its original row/column layout
hu2013_s1b_all_proteins <- readxl::read_excel(
  hu2013_s1b_file,
  col_names = FALSE
) %>%
  unlist(use.names = FALSE) %>%
  as.character() %>%
  stringr::str_trim()

hu2013_s1b_all_proteins <- hu2013_s1b_all_proteins[
  !is.na(hu2013_s1b_all_proteins) &
    hu2013_s1b_all_proteins != ""
]

hu2013_s1b_all_proteins <- tibble::tibble(
  protein = stringr::str_to_upper(hu2013_s1b_all_proteins)
) %>%
  dplyr::distinct(protein)

# Only proteins not already present as positive 1C hits
hu2013_zero_hit_proteins <- hu2013_s1b_all_proteins %>%
  dplyr::anti_join(
    hu2013_protein_summary %>% dplyr::distinct(protein),
    by = "protein"
  ) %>%
  dplyr::mutate(
    Hu2013_n_mCpG_motif_hits = 0L,
    Hu2013_mCpG_motif_ids = NA_character_,
    Hu2013_subfamily = NA_character_,
    Hu2013_no_mCpG_motif_hit = TRUE,
    Hu2013_screen_status =
      "assayed_no_detectable_mCpG_preferential_hit"
  )

# Append the assay-negative S1B proteins to the existing positive-hit summary
hu2013_protein_summary <- dplyr::bind_rows(
  hu2013_protein_summary,
  hu2013_zero_hit_proteins
) %>%
  dplyr::arrange(
    dplyr::desc(Hu2013_n_mCpG_motif_hits),
    protein
  )

# Check the final counts
hu2013_protein_summary %>%
  dplyr::count(Hu2013_screen_status)

hu2013_simple <- hu2013_protein_summary %>%
  dplyr::transmute(
    protein = stringr::str_to_upper(stringr::str_trim(protein)),
    
    Hu2013_mCpG_screen = dplyr::case_when(
      Hu2013_n_mCpG_motif_hits > 0 ~ paste0(
        "mCpG_preferential_hit__",
        Hu2013_n_mCpG_motif_hits,
        "_motifs"
      ),
      
      Hu2013_n_mCpG_motif_hits == 0 ~
        "assayed_no_detectable_mCpG_preferential_hit__0_motifs",
      
      TRUE ~ NA_character_
    )
  )

str(hu2013_simple)


grau2023_s3 <- readxl::read_excel(
  file.path(
    this.dir,
    "literature",
    "Grau2023_Supplementary_Table_S3.xlsx"
  ),
  sheet = "Table_S3_clean",
  na = ""
)

str(grau2023_s3)


grau2023_s3_simple <- grau2023_s3 %>%
  dplyr::transmute(
    protein = stringr::str_to_upper(stringr::str_trim(protein)),
    
    Grau2023_methylation_model = dplyr::recode(
      Grau2023_methylation_model_code,
      "y"  = "methylation_improves_prediction",
      "n"  = "no_methylation_prediction_improvement",
      "i"  = "inconsistent_methylation_prediction_result",
      "NA" = "not_testable_single_cell_type"
    )
  )


song2021_calls <- readr::read_csv(
  file.path(
    this.dir,
    "literature",
    "song2021_symmetric_mC_calls.csv"
  ),
  show_col_types = FALSE
)

song2021_calls_short <- song2021_calls %>%
  dplyr::transmute(
    TF = Song2021_TF_raw,
    Call = Song2021_mC_call
  ) %>%
  dplyr::distinct()


stratified_test_qc_n50_best_experiments_FDR <- readxl::read_excel(
  "C:/Users/Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/QC_permutation_test/stratified_test_qc_n50_best_experiments_FDR.xlsx"
)


batyrev_results <- stratified_test_qc_n50_best_experiments_FDR %>%
  dplyr::filter(
    tested,
    !is.na(sig_direction)
  ) %>%
  dplyr::group_by(protein) %>%
  dplyr::summarise(
    n_significant_positive = sum(sig_direction == "significant_positive"),
    n_significant_negative = sum(sig_direction == "significant_negative"),
    n_not_significant      = sum(sig_direction == "not_significant"),
    .groups = "drop"
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    max_n = max(
      n_significant_positive,
      n_significant_negative,
      n_not_significant
    ),

    Batyrev_methylation_call = dplyr::case_when(

      # Tie involving both opposite significant directions:
      # positive vs negative, optionally also tied with non-significant
      n_significant_positive == max_n &&
        n_significant_negative == max_n &&
        max_n > 0 ~
        "mixed_or_conflicting_directions",

      # Positive wins outright OR ties only with non-significant
      n_significant_positive == max_n &&
        n_significant_positive > n_significant_negative ~
        "significant_higher_methylation",

      # Negative wins outright OR ties only with non-significant
      n_significant_negative == max_n &&
        n_significant_negative > n_significant_positive ~
        "significant_lower_methylation",

      # Non-significant has the sole majority
      n_not_significant == max_n ~
        "not_significant",

      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    protein,
    Batyrev_methylation_call,
    n_significant_positive,
    n_significant_negative,
    n_not_significant
  ) %>%
  dplyr::arrange(protein)

batyrev_results


# # -------------------------------------------------------------------------
# # Batyrev results:
# # One primary chromatin state per protein, selected only by CpG coverage.
# # No methylation result, rank statistic, or p-value is used for selection.
# # -------------------------------------------------------------------------
# 
# batyrev_results <- stratified_test_qc_n50_best_experiments_FDR %>%
#   dplyr::filter(
#     tested,
#     !is.na(sig_direction)
#   ) %>%
#   dplyr::mutate(
#     Batyrev_n_total_CpGs = n_sample1 + n_sample2 + n_both,
#     
#     # Used only as a coverage/balance tie-breaker
#     Batyrev_n_smallest_group = pmin(
#       n_sample1,
#       n_sample2,
#       n_both
#     )
#   ) %>%
#   dplyr::group_by(protein) %>%
#   dplyr::arrange(
#     dplyr::desc(Batyrev_n_total_CpGs),
#     dplyr::desc(Batyrev_n_smallest_group),
#     dplyr::desc(n_both),
#     chromatin_state,
#     .by_group = TRUE
#   ) %>%
#   dplyr::slice_head(n = 1) %>%
#   dplyr::ungroup() %>%
#   dplyr::transmute(
#     protein,
#     
#     Batyrev_methylation_call = dplyr::recode(
#       sig_direction,
#       "significant_positive" = "significant_higher_methylation",
#       "significant_negative" = "significant_lower_methylation",
#       "not_significant" = "not_significant"
#     ),
#     
#     Batyrev_selected_chromatin_state = chromatin_state,
#     Batyrev_selected_motif = motif,
#     
#     Batyrev_selected_biosample1 = biosample1,
#     Batyrev_selected_experiment_id1 = experiment_id1,
#     Batyrev_selected_biosample2 = biosample2,
#     Batyrev_selected_experiment_id2 = experiment_id2,
#     
#     Batyrev_selected_n_sample1 = n_sample1,
#     Batyrev_selected_n_sample2 = n_sample2,
#     Batyrev_selected_n_both = n_both,
#     Batyrev_selected_n_total_CpGs = Batyrev_n_total_CpGs,
#     Batyrev_selected_n_smallest_group = Batyrev_n_smallest_group,
#     
#     Batyrev_selected_S_statistic = observed_S_statistic,
#     Batyrev_selected_p_value_two_sided = p_value_two_sided,
#     Batyrev_selected_sig_direction = sig_direction
#   ) %>%
#   dplyr::arrange(protein)
# 
# batyrev_results

library(dplyr)
library(purrr)

# -------------------------------------------------------------------------
# Prepare one protein-keyed table per source
# No call labels are recoded or changed.
# -------------------------------------------------------------------------

our_results_raw <- batyrev_results %>%
  dplyr::transmute(
    protein,
    Batyrev_methylation_call,
    n_significant_positive,
    n_significant_negative,
    n_not_significant
  )

# our_results_raw <- batyrev_results %>%
#   dplyr::select(
#     protein,
#     Batyrev_methylation_call,
#     Batyrev_selected_chromatin_state,
#     Batyrev_selected_motif,
#     Batyrev_selected_n_total_CpGs,
#     Batyrev_selected_n_smallest_group
#   )

yin2017_results_raw <- yin2017_calls %>%
  dplyr::select(
    protein,
    Yin2017_call
  )

jams2022_results_raw <- jams2022_calls %>%
  dplyr::select(
    protein,
    JAMS2022_call
  )

miodownik2025_results_raw <- miodownik2025_simple %>%
  dplyr::select(
    protein,
    Miodownik2025_methylation_direction
  )

luo2021_results_raw <- luo2021_simple %>%
  dplyr::select(
    protein,
    Luo2021_peak_methylation_label
  )

hu2013_results_raw <- hu2013_simple %>%
  dplyr::select(
    protein,
    Hu2013_mCpG_screen
  )

grau2023_results_raw <- grau2023_s3_simple %>%
  dplyr::select(
    protein,
    Grau2023_methylation_model
  )

song2021_results_raw <- song2021_calls %>%
  dplyr::select(
    protein,
    Song2021_TF_raw,
    Song2021_mC_call
  )

# -------------------------------------------------------------------------
# Merge all sources: one row per protein
# -------------------------------------------------------------------------

all_methylation_results_raw <- list(
  our_results_raw,
  yin2017_results_raw,
  jams2022_results_raw,
  miodownik2025_results_raw,
  luo2021_results_raw,
  hu2013_results_raw,
  grau2023_results_raw,
  song2021_results_raw,
  kaluscha2022_fig5b_results_raw
) %>%
  purrr::reduce(
    dplyr::full_join,
    by = "protein"
  ) %>%
  dplyr::arrange(protein)

# -------------------------------------------------------------------------
# For novelty:
# Hu 2013 counts only when it found a positive mCpG-preferential hit.
# Hu "assayed_no_detectable..." does not count as previous coverage.
# -------------------------------------------------------------------------

hu2013_positive_results_raw <- hu2013_results_raw %>%
  dplyr::filter(
    stringr::str_detect(
      Hu2013_mCpG_screen,
      "^mCpG_preferential_hit"
    )
  ) %>%
  dplyr::select(protein)

external_proteins_positive <- dplyr::bind_rows(
  yin2017_results_raw %>% dplyr::select(protein),
  jams2022_results_raw %>% dplyr::select(protein),
  miodownik2025_results_raw %>% dplyr::select(protein),
  luo2021_results_raw %>% dplyr::select(protein),
  hu2013_positive_results_raw,
  grau2023_results_raw %>% dplyr::select(protein),
  song2021_results_raw %>% dplyr::select(protein)
) %>%
  dplyr::filter(
    !is.na(protein),
    protein != ""
  ) %>%
  dplyr::distinct(protein)

our_unique_proteins <- our_results_raw %>%
  dplyr::anti_join(
    external_proteins_positive,
    by = "protein"
  ) %>%
  dplyr::arrange(protein)

nrow(our_unique_proteins)

our_unique_proteins %>%
  print(n = Inf)

# writexl::write_xlsx(
#   stratified_test_qc_n50_best_experiments_FDR,
#   path = file.path(
#     this.dir,
#     "stratified_test_qc_n50_best_experiments_FDR.xlsx"
#   )
# )


# -------------------------------------------------------------------------
# Helper: compare each external source with Batyrev results
# -------------------------------------------------------------------------

compare_source_to_batyrev <- function(source_df, source_name) {
  
  source_df %>%
    dplyr::select(protein) %>%
    dplyr::filter(!is.na(protein), protein != "") %>%
    dplyr::distinct(protein) %>%
    dplyr::left_join(
      our_results_raw %>%
        dplyr::select(
          protein,
          Batyrev_methylation_call,
          n_significant_positive,
          n_significant_negative,
          n_not_significant
        ),
      by = "protein"
    ) %>%
    dplyr::mutate(
      source = source_name,
      present_in_batyrev = !is.na(Batyrev_methylation_call)
    )
}

source_overlap_all <- dplyr::bind_rows(
  compare_source_to_batyrev(yin2017_results_raw, "Yin2017"),
  compare_source_to_batyrev(jams2022_results_raw, "JAMS2022"),
  compare_source_to_batyrev(miodownik2025_results_raw, "Miodownik2025"),
  compare_source_to_batyrev(luo2021_results_raw, "Luo2021"),
  compare_source_to_batyrev(hu2013_positive_results_raw, "Hu2013_positive_hits"),
  compare_source_to_batyrev(grau2023_results_raw, "Grau2023"),
  compare_source_to_batyrev(song2021_results_raw, "Song2021"),
  compare_source_to_batyrev(kaluscha2022_fig5b_results_raw, "Kaluscha2022_Fig5b_candidates")
)

source_overlap_summary <- source_overlap_all %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(
    n_external_proteins = dplyr::n_distinct(protein),
    n_present_in_batyrev = sum(present_in_batyrev),
    n_not_present_in_batyrev = sum(!present_in_batyrev),
    .groups = "drop"
  ) %>%
  dplyr::arrange(source)

source_overlap_summary

source_overlap_call_counts <- source_overlap_all %>%
  dplyr::filter(present_in_batyrev) %>%
  dplyr::count(
    source,
    Batyrev_methylation_call,
    name = "n_proteins"
  ) %>%
  dplyr::arrange(source, dplyr::desc(n_proteins))

source_overlap_call_counts %>%
  print(n = Inf)

source_overlap_detail <- source_overlap_all %>%
  dplyr::filter(present_in_batyrev) %>%
  dplyr::arrange(source, protein)

source_overlap_detail %>%
  dplyr::select(
    source,
    protein,
    Batyrev_methylation_call,
    n_significant_positive,
    n_significant_negative,
    n_not_significant
  ) %>%
  print(n = Inf)


our_unique_summary <- our_unique_proteins %>%
  dplyr::count(
    Batyrev_methylation_call,
    name = "n_proteins"
  ) %>%
  dplyr::arrange(dplyr::desc(n_proteins))

our_unique_summary



source_overlap_detail %>%
  dplyr::filter(
    source == "Grau2023",
    Batyrev_methylation_call == "significant_higher_methylation"
  ) %>%
  dplyr::select(
    source,
    protein,
    Batyrev_methylation_call,
    n_significant_positive,
    n_significant_negative,
    n_not_significant
  ) %>%
  dplyr::arrange(protein) %>%
  print(n = Inf)

motif_significance_table <- motif_significance_table %>%
  mutate(
    net_significance_score =
      significant_positive - significant_negative,
    .after = total_comparisons
  )

print(motif_significance_table[motif_significance_table$significant_positive > 0,], n = Inf)



is_sig <- stratified_test_qc_n50_best_experiments_FDR$q_value_two_sided < 0.05
is_positive <- stratified_test_qc_n50_best_experiments_FDR$observed_S_statistic > 0  
  
is_interesting <- stratified_test_qc_n50_best_experiments_FDR$protein %in% stratified_test_qc_n50_best_experiments_FDR$protein[is_sig & is_positive]


print(x = stratified_test_qc_n50_best_experiments_FDR[is_interesting,],n = 100)
