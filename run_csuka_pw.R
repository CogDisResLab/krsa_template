# Perform confidence-scoring UKA on differential peptide data.
#
# Directional convention:
#   totalMeanLFC must represent the comparison encoded in the filename.
#   Positive = higher in the first group used to calculate the fold change.
#   Negative = higher in the second group.
#
# For a file containing "p-w":
#   positive totalMeanLFC = p - w
#   negative totalMeanLFC = lower in p than w
#
# Uses mapping database 0.6 and PamGene's score-based confidence table.
# FC-type UKA uses existing fold changes and does not calculate a replicate-based
# significance score. No additional log transformation is applied.

suppressPackageStartupMessages({
  library(purrr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(pgUpstream) # nolint: unused_import_linter.
  library(pgFCS)      # nolint: unused_import_linter.
  library(pgscales)   # nolint: unused_import_linter.
})

# ==============================================================================
# 1. IDENTIFY COMPARISON METADATA
# ==============================================================================

identify_dpp_comparison <- function(signal_path) {
  filename <- basename(signal_path)
  chip <- stringr::str_match(filename, "-(STK|PTK)\\.csv$")[, 2L]
  if (is.na(chip)) stop("Could not determine chip type from: ", signal_path)
  
  # Identify pairwise files by "p-w" anywhere in the filename.
  if (!str_detect(filename, fixed("p-w"))) {
    stop("Expected a pairwise DPP file containing 'p-w': ", signal_path)
  }
  
  comparison_name <- filename |>
    str_remove("^.*-dpp_?") |>
    str_remove("\\.csv$")
  
  tibble(
    signal_path = signal_path,
    run_prefix = str_remove(filename, "-dpp.*$"),
    comparison_name = comparison_name,
    chip = chip
  )
}

# ==============================================================================
# 2. PREPARE SIGNED DIFFERENTIAL PEPTIDE DATA
# ==============================================================================

prepare_signal_data <- function(
    signal_path,
    comparison_name) {
  
  raw_data <- readr::read_csv(
    signal_path,
    show_col_types = FALSE
  )
  required_columns <- c(
    "Peptide",
    "totalMeanLFC"
  )
  missing_columns <- setdiff(
    required_columns,
    colnames(raw_data)
  )
  if (length(missing_columns) > 0L) {
    stop(
      "Differential peptide file is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      "\nFile: ",
      signal_path
    )
  }
  
  signal_data <- raw_data |>
    dplyr::transmute(
      ID = gsub("[/ ]", "", as.character(Peptide)),
      
      # Preserve the original signed fold change exactly as written.
      # No absolute-value transformation and no sign reversal are applied.
      value = as.numeric(totalMeanLFC)
    ) |>
    dplyr::filter(
      !is.na(ID),
      ID != ""
    )
  if (any(!is.finite(signal_data$value))) {
    warning(
      "Missing or infinite totalMeanLFC values were found and will be removed from:\n",
      signal_path
    )
    
    signal_data <- signal_data |>
      dplyr::filter(
        is.finite(value)
      )
  }
  
  duplicate_peptides <- signal_data |>
    dplyr::count(
      ID,
      name = "n"
    ) |>
    dplyr::filter(
      n > 1L
    )
  
  if (nrow(duplicate_peptides) > 0L) {
    duplicate_conflicts <- signal_data |>
      dplyr::group_by(ID) |>
      dplyr::summarise(
        n_values = dplyr::n_distinct(value),
        .groups = "drop"
      ) |>
      dplyr::filter(
        n_values > 1L
      )
    
    if (nrow(duplicate_conflicts) > 0L) {
      stop(
        "Some peptides have multiple conflicting totalMeanLFC values in:\n",
        signal_path,
        "\nExamples: ",
        paste(
          utils::head(
            duplicate_conflicts$ID,
            10L
          ),
          collapse = ", "
        )
      )
    }
    
    signal_data <- signal_data |>
      dplyr::distinct(
        ID,
        value
      )
  }
  
  if (nrow(signal_data) == 0L) {
    stop(
      "No usable differential peptide values remained in:\n",
      signal_path
    )
  }
  
  message("Comparison: ", comparison_name)
  message(
    "Positive totalMeanLFC values are interpreted as higher in ",
    numerator_group,
    "; negative values are interpreted as higher in ",
    denominator_group,
    "."
  )
  
  direction_check <- signal_data |>
    dplyr::summarise(
      n_peptides = dplyr::n(),
      n_positive = sum(value > 0),
      n_negative = sum(value < 0),
      n_zero = sum(value == 0),
      mean_lfc = mean(value),
      median_lfc = median(value),
      minimum_lfc = min(value),
      maximum_lfc = max(value)
    )
  
  print(direction_check)
  
  signal_data |>
    dplyr::mutate(
      ID = factor(ID),
      colSeq = factor(comparison_name)
    ) |>
    dplyr::select(
      ID,
      value,
      colSeq
    )
}


# ==============================================================================
# 3. UKA FUNCTION
# ==============================================================================

filter_stk_antibody_recognition <- function(db_in, antibody_recognition, rescue) {
  if (!antibody_recognition || !"STK_ab_site" %in% names(db_in)) {
    return(db_in)
  }
  
  recognized <- db_in %>% filter(!is.na(STK_ab_site))
  
  if (!rescue) {
    return(recognized)
  }
  
  unrecognized_ids <- setdiff(unique(db_in$ID), unique(recognized$ID))
  if (length(unrecognized_ids) == 0) {
    return(recognized)
  }
  
  rescue_sites <- db_in %>%
    filter(ID %in% unrecognized_ids) %>%
    mutate(.dist_to_center = abs(Peptide_PhosSite - (nchar(as.character(Sequence)) + 1) / 2)) %>%
    arrange(ID, .dist_to_center, desc(PepProtein_SeqSimilarity)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    select(ID, Peptide_PhosSite)
  
  rescued_rows <- db_in %>%
    inner_join(rescue_sites, by = c("ID", "Peptide_PhosSite"))
  
  bind_rows(recognized, rescued_rows)
}

# Filter the mapping database and run PamGene UKA using fold changes.
perform_uka <- function(
    dpp_data,
    upstream_db,
    target_kinase_family,
    kinase_enrichment,
    minimum_sequence_homology = 0.9,
    minimum_phosphonet_score = 300,
    minimum_kl_score = 1.0,
    nperms = 500L,
    max_rank = 12L,
    minimum_set_size = 3L,
    antibody_recognition = TRUE,
    rescue_unrecognized_peptides = TRUE,
    remove_dual_specificity_kinases = TRUE) {
  
  subset_db <- upstream_db |>
    ungroup() |>
    mutate(ID = factor(gsub("[/ ]", "", ID))) |>
    filter(
      Database == "KL" |
        (Database %in% c("iviv", "PhosphoNET") &
           PepProtein_SeqSimilarity >= minimum_sequence_homology)
    ) |>
    filter(
      Database == "iviv" |
        (Database == "PhosphoNET" & Kinase_Score >= minimum_phosphonet_score) |
        (Database == "KL" & Kinase_Score >= minimum_kl_score)
    ) |>
    filter(Kinase_Rank <= max_rank, family == target_kinase_family)
  
  if (target_kinase_family == "STK") {
    subset_db <- filter_stk_antibody_recognition(
      subset_db, antibody_recognition, rescue_unrecognized_peptides
    )
  }
  if (nrow(subset_db) == 0L) stop("No database records remain after filtering.")
  if (length(intersect(as.character(dpp_data$ID), as.character(subset_db$ID))) == 0L) {
    stop("No peptide IDs match the filtered database.")
  }
  
  uka_result <- pgUpstream::pgCscoreAnalysis(
    dpp_data,
    dbFrame = subset_db,
    scores_df = scores_df,
    nPermutations = nperms,
    norm_func = "log10",
    ukaType = "FC",
    score_var = "Kinase_Score"
  )
  
  # FC analysis has no replicate-based significance score.
  uka_result$pPhenoScore <- rep(NA_real_, nrow(uka_result))
  
  uka_result <- uka_result |>
    filter(nFeatures >= minimum_set_size) |>
    left_join(kinase_enrichment, by = c("ClassName" = "Kinase_Name"))
  
  if (remove_dual_specificity_kinases) {
    if (target_kinase_family == "PTK") {
      uka_result <- uka_result |> filter(Kinase_group == "TYR")
    } else {
      uka_result <- uka_result |> filter(Kinase_group != "TYR")
    }
  }
  
  uka_result |>
    arrange(desc(combinedScore)) |>
    transmute(
      `Kinase Name` = ClassName,
      `Kinase Uniprot ID` = Kinase_UniprotID,
      `Kinase Group` = Kinase_group,
      `Kinase Family` = Kinase_family,
      `Significance Score` = pPhenoScore,
      `Specificity Score` = pFeatureScore,
      `Final score` = combinedScore,
      `Kinase Statistic` = NormalizedSetStat,
      `Kinase Change` = delta,
      `Peptide set size` = nFeatures,
      `Final score rank` = row_number(),
      contrast = as.character(dpp_data$colSeq[[1L]]),
      UKAType = "FC",
      UkaDbVersion = "0.6",
      ConfidenceMapping = "Score-based",
      InputSignalScale = "signed totalMeanLFC (unchanged)",
      ConfidenceNormalization = "log10",
      NumberOfPerms = nperms,
      Seed = 42L,
      MinimumSequenceHomology = minimum_sequence_homology,
      MinimumPNScore = minimum_phosphonet_score,
      MinimumKLScore = minimum_kl_score,
      MaxRank = max_rank,
      MinimumSetSize = minimum_set_size,
      AntibodyRecognition = antibody_recognition,
      RescueUnrecognizedPeptides = rescue_unrecognized_peptides,
      rmDualSpecKinases = remove_dual_specificity_kinases
    )
}

# ==============================================================================
# 4. LOCATE DIFFERENTIAL PEPTIDE FILES
# ==============================================================================

signal_files <- list.files(
  path = "results",
  pattern = "-dpp.*-(STK|PTK)\\.csv$",
  full.names = TRUE
)
# Keep only p-minus-w comparisons.
signal_files <- signal_files[
  str_detect(basename(signal_files), fixed("p-w"))
]

# ==============================================================================
# 5. RUN UKA ONLY WHEN PAIRWISE DPP FILES EXIST
# ==============================================================================

if (length(signal_files) == 0L) {
  message("No pairwise DPP files containing 'p-w' were found under results/. Skipping UKA.")
} else {
  
  # 5A. Load reference data.
  if (!"pgCscoreAnalysis" %in% getNamespaceExports("pgUpstream")) {
    stop("Your pgUpstream installation must export pgCscoreAnalysis.")
  }
  
  uka_db <- read_tsv(
    "reference_data/UKAdb_0.6.0_Tier2_260813_seqsim.txt",
    show_col_types = FALSE
  )
  scores_df <- read_csv(
    "reference_data/scoredf_calibrated.csv",
    show_col_types = FALSE
  )
  stopifnot(all(c("Database", "Kinase_Score", "cscore") %in% names(scores_df)))
  if ("Kinase_Rank" %in% names(scores_df)) stop("Use the Score-based confidence table.")
  
  # Get kinase details from the mapping database.
  kinase_enrichment <- uka_db |>
    select(Kinase_Name, Kinase_UniprotID, Kinase_group, Kinase_family) |>
    distinct()
  if (anyDuplicated(kinase_enrichment$Kinase_Name)) {
    stop("Conflicting kinase annotations; use the company's version-0.6 enrichment file.")
  }
  set.seed(42)
  
  
  # 5B. Identify comparisons from their filenames.
  comparison_metadata <- map_dfr(signal_files, identify_dpp_comparison)
  print(comparison_metadata)
  
  # 5C. Run UKA separately for each DPP file.
  uka_results <- purrr::pmap(
    comparison_metadata,
    function(signal_path, run_prefix, comparison_name, chip) {
      message("Running FC-type csUKA: ", comparison_name, " [", chip, "]")
      message("Input file: ", signal_path)
      
      prepared_data <- prepare_signal_data(signal_path, comparison_name)
      result <- perform_uka(
        dpp_data = prepared_data,
        upstream_db = uka_db,
        target_kinase_family = chip,
        kinase_enrichment = kinase_enrichment
      )
      
      output_path <- file.path(
        "results",
        paste0(run_prefix, "-csuka_fc_table_full_", comparison_name, ".csv")
      )
      write_csv(result, output_path)
      message("Wrote: ", output_path)
      result
    }
  )
  names(uka_results) <- basename(comparison_metadata$signal_path)
  
  # 5D. Summarize the directions of kinase changes.
  direction_summary <- map2_dfr(
    uka_results, comparison_metadata$signal_path,
    function(result, signal_path) {
      result |> summarise(
        input_file = basename(signal_path),
        n_kinases = n(),
        n_positive_change = sum(`Kinase Change` > 0, na.rm = TRUE),
        n_negative_change = sum(`Kinase Change` < 0, na.rm = TRUE),
        n_zero_change = sum(`Kinase Change` == 0, na.rm = TRUE),
        median_kinase_change = median(`Kinase Change`, na.rm = TRUE)
      )
    }
  )
  print(direction_summary)
}
