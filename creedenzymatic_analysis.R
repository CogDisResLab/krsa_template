# Creedenzymatic Analysis

library(tidyverse)
library(creedenzymatic)

process_creedenzymatic <- function(krsa_path, uka_path, peptide_path) {
  krsa_data <- read_csv(krsa_path, show_col_types = FALSE) |>
    select(Kinase, Score = AvgZ) |>
    read_krsa(trns = "abs", sort = "desc")

  uka_data <- read_csv(uka_path, show_col_types = FALSE) |>
    select(Kinase = `Kinase Name`, Score = `Median Final score`) |>
    read_uka(trns = "abs", sort = "desc")

  peptide_data <- read_csv(peptide_path, show_col_types = FALSE) |>
    select(Peptide, Score = totalMeanLFC)

  kea3_data <- read_kea(
    peptide_data,
    sort = "asc",
    trns = "abs",
    method = "MeanRank",
    lib = "kinase-substrate"
  )

  ptmsea_data <- read_ptmsea(peptide_data)

  combined <- combine_tools(
    KRSA_df = krsa_data,
    UKA_df = uka_data,
    KEA3_df = kea3_data,
    PTM_SEA_df = ptmsea_data
  )

  combined
}

krsa_files <- list.files("results", "krsa", full.names = TRUE)

uka_files <- list.files("results", "uka", full.names = TRUE)

peptide_files <- list.files("results", "dpp", full.names = TRUE)

comparison_names <- krsa_files |>
  basename() |>
  str_replace(fixed("-krsa_table_"), "_") |>
  str_remove(fixed(".csv"))

result <- list(
  krsa_path = krsa_files,
  uka_path = uka_files,
  peptide_path = peptide_files
) |>
  pmap(process_creedenzymatic) |>
  set_names(comparison_names) |>
  imap_dfr(
    ~ write_csv(
      .x,
      file.path("results", str_glue("{.y}_creedenzymatic.csv"))
    ),
    .id = "Comparison"
  )
