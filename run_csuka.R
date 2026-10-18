# Perform confidence-scoring UKA on KRSA signal files
# Adapted from the supplied company helper.R and main.R.
# Requires csUKA-compatible pgUpstream and its dependencies.
# Retains all scores for kinases meeting the minimum peptide-set size.
# This is the company two-group Comparison model; it does not model donor pairing.

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


# -------------------------------------------------------------------------
# 1. Reference data
# -------------------------------------------------------------------------


input_signal_scale <- "log2" # Signals are already log2-transformed by KRSA.


uka_db <- read_tsv(
  "reference_data/UKAdb_0.6.0_Tier2_260813_seqsim.txt",
  show_col_types = FALSE
)

# PamGene score-to-confidence mapping table.
scores_df <- read_csv(
  "reference_data/scoredf_calibrated.csv",
  show_col_types = FALSE
)
stopifnot(all(c("Database", "Kinase_Score", "cscore") %in% names(scores_df)))
if ("Kinase_Rank" %in% names(scores_df)) {
  stop("Use the Score-based confidence table, without a Kinase_Rank column.")
}

# Use the annotations supplied with this database version.
kinase_enrichment <- uka_db |>
  select(Kinase_Name, Kinase_UniprotID, Kinase_group, Kinase_family) |>
  distinct()


# Match the company's default Random Number Ggenerator seed. 
set.seed(42)

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
      # Observed package behavior indicates that the UKA change is calculated as the second factor level relative to the first factor level.
      grp = factor(
        Group,
        levels = c(control_group, case_group)
      ),
      colSeq = factor(SampleName),
      ID = factor(gsub("[/ ]", "", Peptide)),
      value = as.numeric(slope)
    )
  if (any(!is.finite(signal_data$value))) {
    stop("Signal values must be finite. Check missing/non-numeric values in: ", signal_path)
  }
  if (input_signal_scale == "raw_slope") {
    if (any(signal_data$value <= 0)) {
      stop("Cannot log2-transform nonpositive raw slopes in: ", signal_path)
    }
    signal_data$value <- log2(signal_data$value)
  }
  if (anyNA(signal_data$grp)) {
    stop(
      "Missing phenotype assignments were created while preparing:\n",
      signal_path
    )
  }
  signal_data
}

# -------------------------------------------------------------------------
# 4. UKA function
# -------------------------------------------------------------------------

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

# Adaptation of the company's import_filter_db and Comparison csUKA call.
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

  # Explicit company convention: test = 1, control = 0.
  # Original preparation orders levels as control, case.
  expected_levels <- levels(dpp_data$grp)
  if (length(expected_levels) != 2L) stop("Exactly two groups are required.")
  dpp_data <- dpp_data |>
    mutate(grp = factor(
      ifelse(as.character(grp) == expected_levels[[2L]], "1", "0"),
      levels = c("0", "1")
    ))

  uka_result <- pgUpstream::pgCscoreAnalysis(
    dpp_data,
    dbFrame = subset_db,
    scores_df = scores_df,
    nPermutations = nperms,
    norm_func = "log10",
    ukaType = "Comparison",
    score_var = "Kinase_Score"
  ) |>
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
      contrast = paste(expected_levels[[2L]], "vs", expected_levels[[1L]]),
      UkaDbVersion = "0.6",
      ConfidenceMapping = "Score-based",
      InputSignalScale = input_signal_scale,
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

# -------------------------------------------------------------------------
# 5. Define comparisons and locate their signal files
# -------------------------------------------------------------------------

# First group = case; second group = reference/control.
# Edit the order to match your intended KRSA comparisons.
comparisons <- list(
  COMP1 = c("F_NHP", "F_NHP_HI"),
  COMP2 = c("M_NHP", "M_NHP_HI")
)

chip_type <- "STK"

signal_files <- list.files(
  path = "results",
  pattern = "-signal_.*-(STK|PTK)\\.csv$",
  full.names = TRUE
)

signal_files <- signal_files[
  !str_detect(basename(signal_files), fixed("p-w"))
] # Exclude the p-w comparisons, as in the original script.

comparison_metadata <- map_dfr(
  comparisons,
  function(comparison) {
    case_group <- comparison[[1L]]
    control_group <- comparison[[2L]]
    
    expected_suffix <- paste0(
      "-signal_",
      case_group, "_", control_group,
      "-", chip_type, ".csv"
    )
    
    matching_files <- signal_files[
      endsWith(basename(signal_files), expected_suffix)
    ]
    
    if (length(matching_files) != 1L) {
      stop(
        "Expected exactly one signal file ending in: ",
        expected_suffix,
        "\nFound: ", length(matching_files)
      )
    }
    
    signal_data <- read_csv(
      matching_files[[1L]],
      show_col_types = FALSE
    )
    
    if (!setequal(
      as.character(signal_data$Group),
      c(case_group, control_group)
    )) {
      stop("Signal-file groups do not match: ", expected_suffix)
    }
    
    tibble(
      signal_path = matching_files[[1L]],
      comparison_name = paste0(
        case_group, "_", control_group, "-", chip_type
      ),
      case_group = case_group,
      control_group = control_group,
      chip = chip_type
    )
  }
)

print(comparison_metadata)

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
      "Global mean difference on the analysis scale, ",
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
        "-csuka_table_full_",
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
          `Kinase Change` > 0,
          na.rm = TRUE
        ),
        n_negative_change = sum(
          `Kinase Change` < 0,
          na.rm = TRUE
        ),
        n_zero_change = sum(
          `Kinase Change` == 0,
          na.rm = TRUE
        ),
        median_kinase_change_across_kinases = median(
          `Kinase Change`,
          na.rm = TRUE
        )
      )
  }
)

print(direction_summary)

