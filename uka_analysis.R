# Perform UKA analysis on KRSA signal files

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

# -------------------------------------------------------------------------
# 1. Reference data
# -------------------------------------------------------------------------

uka_db <- readRDS(
  file.path(
    "reference_data",
    "UKA_231031-86502-87102_UpstreamDb.rds"
  )
)

kinase_enrichment <- read_csv(
  file.path(
    "reference_data",
    "UKA_Kinase_enrichment.csv"
  ),
  show_col_types = FALSE
)

# -------------------------------------------------------------------------
# 2. Identify comparison metadata from a signal file
# -------------------------------------------------------------------------

identify_signal_comparison <- function(signal_path) {
  signal_data <- read_csv(
    signal_path,
    show_col_types = FALSE
  )

  required_columns <- c(
    "Group",
    "SampleName",
    "slope",
    "Peptide"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(signal_data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Signal file is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      "\nFile: ",
      signal_path
    )
  }

  groups <- signal_data |>
    pull(Group) |>
    as.character() |>
    unique()

  if (length(groups) != 2L) {
    stop(
      "UKA requires exactly two groups, but ",
      length(groups),
      " were found in:\n",
      signal_path,
      "\nGroups: ",
      paste(groups, collapse = ", ")
    )
  }

  # The reporting template writes comparisons as:
  # case versus control, with control conventionally named CTL_*.
  control_candidates <- groups[
    str_detect(
      groups,
      regex("(^CTL($|_)|control)", ignore_case = TRUE)
    )
  ]

  if (length(control_candidates) != 1L) {
    stop(
      "Could not identify exactly one control group in:\n",
      signal_path,
      "\nGroups found: ",
      paste(groups, collapse = ", "),
      "\nExpected one group beginning with 'CTL' or containing 'control'."
    )
  }

  control_group <- control_candidates[[1L]]
  case_group <- setdiff(groups, control_group)

  if (length(case_group) != 1L) {
    stop(
      "Could not identify exactly one case group in:\n",
      signal_path
    )
  }

  chip_match <- basename(signal_path) |>
  str_match("-(STK|PTK)\\.csv$")

  chip <- chip_match[1L, 2L]

  if (is.na(chip)) {
    stop(
      "Could not determine chip type from filename:\n",
      signal_path,
      "\nExpected the filename to end in '-STK.csv' or '-PTK.csv'."
    )
  }

  comparison_name <- basename(signal_path) |>
    str_remove("^.*-signal_") |>
    str_remove("\\.csv$")

  tibble(
    signal_path = signal_path,
    comparison_name = comparison_name,
    case_group = case_group,
    control_group = control_group,
    chip = chip
  )
}

# -------------------------------------------------------------------------
# 3. Prepare signal data for UKA
# -------------------------------------------------------------------------

prepare_signal_data <- function(
    signal_path,
    case_group,
    control_group) {
  signal_data <- read_csv(
    signal_path,
    show_col_types = FALSE
  ) |>
    filter(
      Group %in% c(case_group, control_group)
    ) |>
    transmute(
      # Explicitly order the factor so that UKA returns the directional
      # comparison as case minus control.
      #
      # Observed package behavior indicates that the UKA change is calculated
      # as the second factor level relative to the first factor level.
      grp = factor(
        Group,
        levels = c(control_group, case_group)
      ),
      colSeq = factor(SampleName),
      ID = factor(Peptide),
      value = as.numeric(slope)
    )

  if (anyNA(signal_data$grp)) {
    stop(
      "Missing phenotype assignments were created while preparing:\n",
      signal_path
    )
  }

  if (anyNA(signal_data$value)) {
    warning(
      "Missing slope values were detected in:\n",
      signal_path
    )
  }

  signal_data
}

# -------------------------------------------------------------------------
# 4. UKA function
# -------------------------------------------------------------------------

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

  subset_db <- upstream_db |>
    ungroup() |>
    filter(
      PepProtein_SeqSimilarity >= minimum_sequence_homology,
      .data$family == .env$target_kinase_family,
      Kinase_Rank <= max_rank
    ) |>
    filter(
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
    distinct(family) |>
    pull(family)

  if (
    length(retained_families) != 1L ||
      retained_families[[1L]] != target_kinase_family
  ) {
    stop(
      "Unexpected kinase families remained after filtering: ",
      paste(retained_families, collapse = ", ")
    )
  }

  uka_result <- pgScanAnalysis2g_np(
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
        mutate(mxRank = x[["mxRank"]])
    }) |>
    bind_rows() |>
    filter(nFeatures >= minimum_set_size) |>
    makeSummary() |>
    select(
      Kinase = ClassName,
      everything()
    ) |>
    left_join(
      kinase_enrichment,
      by = c("Kinase" = "Kinase_Name")
    ) |>
    arrange(desc(medianScore)) |>
    select(
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

# -------------------------------------------------------------------------
# 5. Locate signal files
# -------------------------------------------------------------------------

signal_files <- list.files(
  path = "results",
  pattern = "-signal_.*-(STK|PTK)\\.csv$",
  full.names = TRUE
)

# Exclude the p-w comparisons, as in the original script.
signal_files <- signal_files[
  !str_detect(
    basename(signal_files),
    fixed("p-w")
  )
]

if (length(signal_files) == 0L) {
  stop(
    "No UKA-compatible signal files were found under results/."
  )
}

comparison_metadata <- map_dfr(
  signal_files,
  identify_signal_comparison
)

print(
  comparison_metadata |>
    select(
      comparison_name,
      case_group,
      control_group,
      chip
    )
)

# -------------------------------------------------------------------------
# 6. Run UKA
# -------------------------------------------------------------------------

uka_results <- pmap(
  comparison_metadata,
  function(
      signal_path,
      comparison_name,
      case_group,
      control_group,
      chip) {
    message("")
    message(
      "Running UKA: ",
      case_group,
      " - ",
      control_group,
      " [",
      chip,
      "]"
    )

    prepared_data <- prepare_signal_data(
      signal_path = signal_path,
      case_group = case_group,
      control_group = control_group
    )

    message(
      "UKA phenotype levels: ",
      paste(
        levels(prepared_data$grp),
        collapse = " -> "
      )
    )

    signal_summary <- prepared_data |>
      group_by(grp) |>
      summarise(
        n_samples = n_distinct(colSeq),
        n_peptides = n_distinct(ID),
        mean_signal = mean(value, na.rm = TRUE),
        median_signal = median(value, na.rm = TRUE),
        .groups = "drop"
      )

    print(signal_summary)

    control_mean <- signal_summary |>
      filter(
        grp == control_group
      ) |>
      pull(mean_signal)

    case_mean <- signal_summary |>
      filter(
        grp == case_group
      ) |>
      pull(mean_signal)

    message(
      "Observed global mean difference, ",
      case_group,
      " - ",
      control_group,
      ": ",
      round(case_mean - control_mean, 4L)
    )

    uka_result <- perform_uka(
      dpp_data = prepared_data,
      upstream_db = uka_db,
      target_kinase_family = chip,
      kinase_enrichment = kinase_enrichment
    )

    output_path <- file.path(
      "results",
      str_glue(
        "{str_remove(basename(signal_path), '-signal_.*$')}",
        "-uka_table_full_",
        "{comparison_name}.csv"
      )
    )

    write_csv(
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

# -------------------------------------------------------------------------
# 7. Directional diagnostic summary
# -------------------------------------------------------------------------

direction_summary <- map2_dfr(
  uka_results,
  comparison_metadata$comparison_name,
  function(result, comparison_name) {
    result |>
      summarise(
        comparison = comparison_name,
        n_kinases = n(),
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

