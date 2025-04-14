# KRSA Report and Kinome Analysis Template

## Overview

This repository provides a standardized, reproducible workflow for Kinome Analysis using tools from the [Cognitive Disorders Research Lab](https://cdrl-ut.org).

## Features

- Generate comprehensive KRSA reports for kinome array experiments
- Perform Universal Kinase Analysis (UKA)
- Create Creedenzymatic results
- Generate detailed quartile plots
- Standardize data and results management
- Ensure reproducible research workflows

## Requirements

### Software Prerequisites

1. **R** (Statistical Computing)
   - Download from [R Project Website](https://www.r-project.org/)
   - Recommended: Latest stable version

2. **Git** (Version Control)
   - Download from [Git Website](https://git-scm.com/downloads)
   - Ensures proper version tracking

3. **Just** (Task Runner)
   - Install via: The installation documentation [Just Documentation](https://just.systems/)
   - Manages project workflows

4. **Quarto** (Rendering Engine)
   - Install via: The installers on the [Quarto Website](https://quarto.org/docs/download/)
   - Handles the conversion of R Markdown files to Markdown with the R code
   evaluated.

5. **Text Editor/IDE**
   Recommended options:
   - [Visual Studio Code](https://code.visualstudio.com/)
   - [RStudio](https://posit.co/download/rstudio-desktop/)
   - [Neovim](https://neovim.io/)

6. **Additional Requirements**
   - Windows Only: [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
   - Optional: [Radian](https://github.com/randy3k/radian) (Enhanced R Console)

## Quickstart Guide

### 1. Repository Setup

1. Click "Use Template" to create a new repository
2. Clone the repository to your local machine

   ```bash
   git clone <your-repository-url>
   cd <repository-name>
   ```

### 2. Environment Preparation

1. Open R in the repository directory
2. Restore project dependencies

   ```r
   renv::restore()
   ```

### 3. Data Preparation

1. Place `SigmBg` and `SignalSaturation` files in the `kinome_data/` directory
2. Create a new analysis file from the template

   ```bash
   just new-analysis my_experiment
   ```

### 4. Configuration

1. Edit the newly created `.Rmd` file
2. Update file paths in the YAML frontmatter
3. Set chip type and prefix as needed

### 5. Run Analysis

Execute the complete workflow:

```bash
just render
just all
```

## Available Just Commands

| Command | Description |
|---------|-------------|
| `just list` | Show all available commands |
| `just new-analysis NAME CHIP_TYPE PREFIX` | Create a new analysis from template |
| `just render` | Render all Markdown reports |
| `just uka` | Run Universal Kinase Analysis |
| `just creeden` | Run Creedenzymatic analysis and generate quartile plots |
| `just clean` | Remove generated files and artifacts |

## Troubleshooting

- Ensure all dependencies are installed
- Check R and Just versions
- Verify data file locations

## License

Everything in this project is licensed under the MIT license except for the `UKA_*`
files in the `reference_data` folder.

## Contact

In case of any problems, please file an issue
