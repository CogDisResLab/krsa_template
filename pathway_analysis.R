# Pathway Analysis STK - Su Data

library(tidyverse)
library(enrichR)
library(writexl)

# Load the reference data
stk_id <- readRDS("reference_data/stk_id_map.Rds")
stk_hgnc <- readRDS("reference_data/stk_hgnc_map.Rds")

stk_map <- stk_id |>
  inner_join(stk_hgnc)

# Copy this section for each comparison
#################

# Change the name of the file you are loading
sh_comparison <- read_csv("results/dpp_HSF1-SH_CTL-STK.csv") |>
  group_by(Peptide) |>
  filter(abs(LFC) == max(abs(LFC))) |>
  ungroup() |>
  select(Peptide, LFC) |>
  inner_join(stk_map) |>
  # Change the name of the file you are writing.
  write_csv("results/annotated_dpp_HSF1-SH_CTL-STK.csv")

sh_genes <- sh_comparison |>
  select(Gene, LFC) |>
  filter(LFC >= quantile(sh_comparison[["LFC"]], 0.90) |
    LFC <= quantile(sh_comparison[["LFC"]], 0.10)) |>
  pull(Gene)

dbs <- c(
  "GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
  "GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022"
)

sh_enriched <- enrichr(sh_genes, dbs) |>
  imap(~ write_csv(.x, str_glue("results/HSF1-SH_CTL-STK-{.y}-Pathways.csv")))

# Change the name of the file you are writing
write_xlsx(sh_enriched, "results/SH-STK-Pathways.xlsx")

#################
