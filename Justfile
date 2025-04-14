# Justfile for R Analysis Project with renv Integration

# Set options
set unstable

# Default target when no specific target is specified
default: all

# Render all Rmd files
render:
    #!/usr/bin/env Rscript
    rmd_files <- list.files(path = ".", pattern = "^((?!_).)*\\.Rmd$", full.names = TRUE)
    for (rmd_file in rmd_files) {
        cat("Rendering", rmd_file, "\n")
        system2("quarto", c("render", rmd_file))
    }

# Run specific R script analyses
uka:
    Rscript uka_analysis.R
    Rscript uka_analysis_single_sample.R

creeden:
    Rscript creedenzymatic_analysis.R
    Rscript generate_quartile_plots.R

# Workflow for statistical analysis
analysis: uka creeden

# Utility to generate new analysis from template
new-analysis NAME chiptype='STK' prefix='kinome':
    #!/usr/bin/env Rscript
    template_file <- "_template.Rmd"
    new_file <- paste0("{{NAME}}.Rmd")
    file.copy(template_file, new_file)

    rmd_lines <- readLines(new_file)
    rmd_lines <- sub("chip_type:.*$", paste0("chip_type: {{chiptype}}"), rmd_lines)
    rmd_lines <- sub("prefix:.*$", paste0("prefix: {{prefix}}"), rmd_lines)
    writeLines(rmd_lines, new_file)

# Clean up generated files and artifacts
clean:
    #!/usr/bin/env Rscript
    unlink(list.files(pattern = "\\.(html|docx|pdf)$"), recursive = TRUE)
    unlink(list.files("results/", pattern = "\\.csv$", full.names = TRUE), recursive = TRUE)
    unlink(list.files("figures/", pattern = "\\.(png|svg)$", full.names = TRUE), recursive = TRUE)
    unlink(list.files("datastore/", pattern = "\\.RData$", full.names = TRUE), recursive = TRUE)
    unlink(c("_cache", "_freeze"), recursive = TRUE)

# List available targets
list:
    @just --list

# Run all tasks in sequence
all: render analysis

# Aliases for common tasks
alias a := all
alias c := clean
alias n := new-analysis
