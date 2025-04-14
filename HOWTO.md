# KRSA Analysis: Detailed Workflow Guide

## Comprehensive Step-by-Step Procedure

### Preparation Phase

#### 1. Software Installation

1. Install R from [R Project Website](https://www.r-project.org/)
   - Verify installation: `R --version`
2. Install Git from [Git Website](https://git-scm.com/downloads)
   - Verify installation: `git --version`
3. Install [Just task runner](https://just.systems)
   Recommended Methods:
   - macOS: using `homebrew`

   ```bash
   brew install just
   ```

   - Linux: using `snap` or your distribution's package manager
   - Windows: Use one of the recommended methods from the
     [website](https://just.systems/man/en/packages.html#windows)
4. Install Quarto CLI from the [Quarto Website](https://quarto.org/docs/download/)
5. Install recommended IDE (VSCode, RStudio)
   Recommended options:
   - [Visual Studio Code](https://code.visualstudio.com/)
   - [RStudio](https://posit.co/download/rstudio-desktop/)
   - [Neovim](https://neovim.io/)

#### 2. Repository Setup

1. Navigate to the [template repository](https://www.github.com/CogDisResLab/krsa_template)
2. Click "Use Template" to create your own repository
3. Clone your new repository

   ```bash
   git clone https://github.com/your-username/your-repository.git
   cd your-repository
   ```

### Environment Configuration

#### 3. R Environment Preparation

> [!NOTE]
> This is the step that usually takes the most time. If in doubt,
try running `renv::hydrate` and `renv::update`

1. Open R in the repository directory
2. Initialize and restore project dependencies

   ```r
   # Restore project-specific packages
   renv::restore()
   
   # Verify installation
   renv::status()
   ```

### Data Preparation

#### 4. Data File Management

1. Locate your experimental data files:
   - Each experiment has two data files. Each file has a specific naming scheme
   that can be used to identify its contents:
      - **Signal Intensity File**: This file contains the signal data after
      subtracting the background intensity. These can be identified by the
      presence of the string `SigmBg` in the filename
      - **Signal Saturation File**: This file contains the data about
      he times when the detected signal intensity was higher than the maximum
      detectable signal. These can be identified by the presence
      of the string `SignalSaturation` in the filename.
2. Place files in `kinome_data/` directory
   - The structure in the `kinome_data/` directory is entirely up to personal
   preference. The convention is to separate by chips (STK and PTK) in different
   directories and then further use some kind of identifying prefix for the
   files that matches the analysis report prefix (see below).

### Analysis Setup

#### 5. Creating a New Analysis

1. Use Just to create a new analysis file

   ```bash
   # Basic usage
   just new-analysis my_experiment
   
   # Specify chip type (PTK or STK)
   just new-analysis my_experiment STK
   
   # Customize prefix
   just new-analysis my_experiment STK custom_prefix
   ```

#### 6. Configuration

1. Open the newly created `.Rmd` file
2. Edit YAML frontmatter:

   ```yaml
   params:
      title: "Add Title Here" # Add the title of the report
      subtitle: "Customer Name Here" # Add the customer/collaborator name
   ```

3. Update file paths for input data

   ```yaml
   params:
      signal_file: "kinome_data/your/file/name/here"
      saturation_file: "kinome_data/your/file/name/here"
   ```

4. Review analysis parameters

   ```yaml
   params:
      threshold: 2 # This is cosmetic and is used to mark significant Z scores
      # This is used to identify and manage multiple experiments in the same repository
      prefix: "kinome" 
      pairwise: FALSE # Set to TRUE if you are comparing one well to another well without replicates
   ```

5. Run the analysis until the group comparison chunk.

   In many cases, it is either hard or impossible to identify the group
   comparisons in advance. Similarly, unless you generated the data yourself,
   the naming scheme of the groups and sample names might not be what you expect.

   Thus, it becomes important to find out this information and set up the correct
   comparisons in this report.

   For this reason, we recommend running the report chunk by chunk until
   the *Group Comparison* chunk.

   This should give you the variable `groups` which has the extracted group
   names. This should allow you to identify the groups and comparisons.

6. Fill out the `comparisons` variable with the required comparisons.

   With the comparisons in hand, you can go ahead and edit the comparisons
   chunk to include the comparisons needed. It is important to note that
   the second group mentioned will be treated as a control group, against
   which the other group will be compared.

   ```{r}
   # Define Groups to be compared
   comparisons <- list(
     COMP1 = c(groups[[1L]], groups[[4L]]),
     COMP2 = c(groups[[1L]], groups[[3L]]),
     COMP3 = c(groups[[2L]], groups[[4L]]),
     COMP4 = c(groups[[4L]], groups[[3L]])
   )
   ```

   These groups are defined in terms of the `groups` variable. This ensures
   that the group names are sourced correctly.

7. Copy the `childA` chunk as many times as there are comparisons.
  
   The `childA` chunk is the main copy of what performs the group comparison
   analysis. You need to copy that as many times as needed. Once copied, each
   new copy would need two places changed. The `label` statement in the chunk
   options and the `random` variable in the chunk itself. These should be named
   sequential based on the alphabet as a suffix of `child`. Thus, the second
   chunk could be titled `childB` with the corresponding `random` variable being
   set to `B` as well. The following chunk would be `childC` and so on.

### Execution

#### 7. Run Complete Analysis

```bash
# Run entire workflow
just all

# Or run specific components
just render   # Generate reports
just uka      # Run Universal Kinase Analysis
just creeden  # Creedenzymatic analysis
```

It is important to run the steps in this order. The `render` step generates
the files consumed by both `uka` and `creeden` steps, while the `creeden`
step relies on the files generated by the `uka` step as well.

Failing to follow the order may result in cryptic errors.

### Output Management

#### 8. Interpreting Results

- Examine rendered PDF reports
- Check `results/` for CSV files
- Review `figures/` for generated plots

### Troubleshooting

#### 9. Common Issues

- **Dependencies**: `renv` can be unpredictable at times. In case of errors
   loading libraries, ensure that you've run `renv::hydrate()`
- **Incorrect Chip Type**: If you see an error implying that there are no
   peptides matching, you might have set the incorrect `chip_type` in the
   frontmatter. Ensure that the `chip_type` parameter is set correctly.

#### 10. Verification Steps

1. Confirm R and Just versions
2. Check `renv.lock` file
3. Verify data file integrity

## Best Practices

- Always use `renv` for dependency management
- Commit your `renv.lock` file to version control
- Keep input data files unchanged

## Getting Help

- Reach out to lab support
