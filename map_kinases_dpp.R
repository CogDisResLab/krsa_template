# Generate a list of kinases mapped to the LFC fo differentially phosphorylated peptides

library(tidyverse)
library(KRSA)

# Prepare the kinase data

ptk_coverage <- KRSA_coverage_PTK_PamChip_86402_v1 |>
  rename(Kinase = Kin, Peptide = Substrates)
stk_coverage <- KRSA_coverage_STK_PamChip_87102_v2 |>
  rename(Kinase = Kin, Peptide = Substrates)

coverage <- bind_rows(ptk_coverage, stk_coverage)

add_kinases <- function(dpp_file, coverage) {
  dpp_data <- read_csv(dpp_file, show_col_types = FALSE) |>
    select(Peptide, LFC = totalMeanLFC)

  dpp_data |>
    inner_join(coverage, by = "Peptide", relationship = "many-to-many") |>
    select(Kinase, Peptide, LFC)
}

dpp_files <- list.files("results", "dpp", full.names = TRUE) |>
  set_names(
    ~ basename(.x) |>
      str_remove(fixed(".csv")) |>
      str_remove(fixed("-dpp"))
  )

dpp_data <- dpp_files |>
  map(~ add_kinases(.x, coverage)) |>
  imap(~ write_csv(
    .x,
    file = file.path("results", str_glue("{.y}-mapped_kinases.csv"))
  ))
