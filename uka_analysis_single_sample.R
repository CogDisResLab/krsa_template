# Perform UKA analysis on differential peptide data
#
# Directional convention:
#   totalMeanLFC must represent the comparison encoded in the filename.
#
# For a file containing "p-w":
#   positive totalMeanLFC = p - w
#   negative totalMeanLFC = lower in p than w
#
# This script never takes the absolute value, reverses the sign, or derives
# direction from factor ordering.

suppressPackageStartupMessages({
  library(purrr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(pgUpstream) # nolint: unused_import_linter.
  library(pgFCS)      # nolint: unused_import_linter.
  library(pgscales)   # nolint: unused_import_linter.
})

# ==============================================================================
# 1. IDENTIFY COMPARISON METADATA
# ==============================================================================

identify_dpp_comparison <- function(signal_path) {
  filename <- basename(signal_path)

  chip <- stringr::str_match(
    filename,
    "-(STK|PTK)\\.csv$"
  )[, 2L]

  if (is.na(chip)) {
    stop(
      "Could not determine chip type from filename:\n",
      signal_path,
      "\nExpected the filename to end in '-STK.csv' or '-PTK.csv'."
    )
  }

  comparison_name <- filename |>
    stringr::str_remove("^.*-dpp_?") |>
    stringr::str_remove("\\.csv$")

  if (!stringr::str_detect(
    comparison_name,
    stringr::fixed("p-w")
  )) {
    stop(
      "The selected differential peptide file does not contain the ",
      "expected directional comparison 'p-w':\n",
      signal_path,
      "\nComparison parsed as: ",
      comparison_name
    )
  }

  run_prefix <- filename |>
    stringr::str_remove("-dpp.*$")

  tibble::tibble(
    signal_path = signal_path,
    run_prefix = run_prefix,
    comparison_name = comparison_name,
    numerator_group = "p",
    denominator_group = "w",
    direction_label = "p - w",
    chip = chip
  )
}


# ==============================================================================
# 2. PREPARE SIGNED DIFFERENTIAL PEPTIDE DATA
# ==============================================================================

prepare_signal_data <- function(
    signal_path,
    numerator_group,
    denominator_group) {

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
      ID = as.character(Peptide),

      # Preserve the original signed fold change exactly as written.
      # No absolute-value transformation and no sign reversal are applied.
      value = as.numeric(totalMeanLFC)
    ) |>
    dplyr::filter(
      !is.na(ID),
      ID != ""
    )

  if (anyNA(signal_data$value)) {
    warning(
      "Missing totalMeanLFC values were found and will be removed from:\n",
      signal_path
    )

    signal_data <- signal_data |>
      dplyr::filter(
        !is.na(value)
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

  message(
    "Directional input: ",
    numerator_group,
    " - ",
    denominator_group
  )

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
      ID = factor(ID)
    ) |>
    dplyr::select(
      ID,
      value
    )
}


# ==============================================================================
# 3. UKA FUNCTION
# ==============================================================================

perform_uka <- function(
    dpp_data,
    upstream_db,
    target_kinase_family,
    kinase_enrichment,
    minimum_sequence_homology = 0.9,
    minimum_phosphonet_score = 300L,
    nperms = 500L,
    min_rank = 4L,
    max_rank = 12L,
    weight_iviv = 1L,
    weight_pnet = 1L,
    minimum_set_size = 3L) {

  if (!target_kinase_family %in% c("STK", "PTK")) {
    stop(
      "Unsupported kinase family: ",
      target_kinase_family,
      ". Expected 'STK' or 'PTK'."
    )
  }

  subset_db <- upstream_db |>
    dplyr::ungroup() |>
    dplyr::filter(
      PepProtein_SeqSimilarity >= minimum_sequence_homology,

      # Explicitly compare the database column to the function argument.
      .data$family == .env$target_kinase_family,

      Kinase_Rank <= max_rank
    ) |>
    dplyr::filter(
      Kinase_PKinase_PredictorVersion2Score >=
        minimum_phosphonet_score |
        Database == "iviv"
    )

  if (nrow(subset_db) == 0L) {
    stop(
      "No UKA database records remained for kinase family ",
      target_kinase_family,
      "."
    )
  }

  retained_families <- subset_db |>
    dplyr::distinct(family) |>
    dplyr::pull(family)

  if (
    length(retained_families) != 1L ||
      retained_families[[1L]] != target_kinase_family
  ) {
    stop(
      "Unexpected kinase families remained after filtering: ",
      paste(
        retained_families,
        collapse = ", "
      )
    )
  }

  uka_result <- pgScanAnalysis0_np(
    dpp_data,
    subset_db,
    nPermutations = nperms,
    scanRank = min_rank:max_rank,
    dbWeights = c(
      iviv = weight_iviv,
      PhosphoNET = weight_pnet
    )
  ) |>
    purrr::map(function(x) {
      x[["aResult"]][[1L]] |>
        dplyr::mutate(
          mxRank = x[["mxRank"]]
        )
    }) |>
    dplyr::bind_rows() |>
    dplyr::filter(
      nFeatures >= minimum_set_size
    ) |>
    makeSummary() |>
    dplyr::select(
      Kinase = ClassName,
      dplyr::everything()
    ) |>
    dplyr::left_join(
      kinase_enrichment,
      by = c(
        "Kinase" = "Kinase_Name"
      )
    ) |>
    dplyr::arrange(
      dplyr::desc(medianScore)
    ) |>
    dplyr::select(
      `Kinase Name` = Kinase,
      `Kinase Uniprot ID` = Kinase_UniprotID,
      `Kinase Group` = Kinase_group,
      `Kinase Family` = Kinase_family,
      `Mean Significance Score` = meanPhenoScore,
      `Mean Specificity Score` = meanFeatScore,
      `Median Final score` = medianScore,
      `Max Final score` = maxScore,
      `Median Kinase Statistic` = medianStat,
      `Mean Kinase Statistic` = meanStat,
      `SD Kinase Statistic` = sdStat,
      `Median Kinase Change` = medianDelta,
      `Mean peptide set size` = meanSetSize
    )

  uka_result
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
  stringr::str_detect(
    basename(signal_files),
    stringr::fixed("p-w")
  )
]


# ==============================================================================
# 5. RUN UKA ONLY WHEN P-W FILES EXIST
# ==============================================================================

if (length(signal_files) == 0L) {

  message(
    "No p-w differential peptide files ending in ",
    "'-STK.csv' or '-PTK.csv' were found under results/. ",
    "Skipping UKA analysis."
  )

} else {

  # ============================================================================
  # 5A. LOAD REFERENCE DATA
  # ============================================================================

  uka_db <- readRDS(
    file.path(
      "reference_data",
      "UKA_231031-86502-87102_UpstreamDb.rds"
    )
  )

  kinase_enrichment <- readr::read_csv(
    file.path(
      "reference_data",
      "UKA_Kinase_enrichment.csv"
    ),
    show_col_types = FALSE
  )


  # ============================================================================
  # 5B. IDENTIFY COMPARISONS
  # ============================================================================

  comparison_metadata <- purrr::map_dfr(
    signal_files,
    identify_dpp_comparison
  )

  print(
    comparison_metadata |>
      dplyr::select(
        comparison_name,
        numerator_group,
        denominator_group,
        direction_label,
        chip
      )
  )


  # ============================================================================
  # 5C. RUN UKA
  # ============================================================================

  uka_results <- purrr::pmap(
    comparison_metadata,
    function(
        signal_path,
        run_prefix,
        comparison_name,
        numerator_group,
        denominator_group,
        direction_label,
        chip) {

      message("")
      message(
        "Running UKA: ",
        direction_label,
        " [",
        chip,
        "]"
      )

      message(
        "Input file: ",
        signal_path
      )

      prepared_data <- prepare_signal_data(
        signal_path = signal_path,
        numerator_group = numerator_group,
        denominator_group = denominator_group
      )

      uka_result <- perform_uka(
        dpp_data = prepared_data,
        upstream_db = uka_db,
        target_kinase_family = chip,
        kinase_enrichment = kinase_enrichment
      )

      output_path <- file.path(
        "results",
        stringr::str_glue(
          "{run_prefix}-uka_table_full_{comparison_name}.csv"
        )
      )

      readr::write_csv(
        uka_result,
        output_path
      )

      message(
        "Wrote: ",
        output_path
      )

      uka_result
    }
  )

  names(uka_results) <- comparison_metadata$comparison_name


  # ============================================================================
  # 5D. DIRECTIONAL DIAGNOSTIC SUMMARY
  # ============================================================================

  direction_summary <- purrr::map2_dfr(
    uka_results,
    comparison_metadata$comparison_name,
    function(result, comparison_name) {
      result |>
        dplyr::summarise(
          comparison = comparison_name,
          n_kinases = dplyr::n(),
          n_positive_change = sum(
            `Median Kinase Change` > 0,
            na.rm = TRUE
          ),
          n_negative_change = sum(
            `Median Kinase Change` < 0,
            na.rm = TRUE
          ),
          n_zero_change = sum(
            `Median Kinase Change` == 0,
            na.rm = TRUE
          ),
          median_kinase_change = median(
            `Median Kinase Change`,
            na.rm = TRUE
          )
        )
    }
  )

  print(direction_summary)
}
