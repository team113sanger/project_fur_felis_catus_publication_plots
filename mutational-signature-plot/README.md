# fur_mutational_signatures

This repository contains scripts for analyzing mutational signatures in feline cancer samples. These tools extract and identify mutational signatures from DNA sequencing data, helping to characterize the genomic landscape of various feline cancer types.

## Table of Contents
- [Repository Structure](#repository-structure)
- [Scripts Overview](#scripts-overview)
- [Requirements](#requirements)
- [Usage](#usage)
- [Notes](#notes)

## Repository Structure

```
├── README.md
├── analysis/                        # Output results directory
│   ├── 6555_2711/                   # Results by study ID
│   ├── 6711_2820/
│   └── ...
├── resources/                       # Reference data
│   ├── baitset/                     # Sequencing target regions
│   └── signatures/                  # Reference signature matrices
├── renv.lock                        # R environment lock file
├── scripts/                         # Analysis scripts
├── source_me.sh                     # Environment configuration
└── study_manifest.tsv               # Sample cohort definitions
```


## Scripts Overview

- **run_sigfit_extraction_and_fitting.R**: The main analysis script that:
  - Extracts de novo signatures from mutation catalogs
  - Calculates goodness-of-fit for different signature counts
  - Compares extracted signatures to reference signatures (COSMIC and custom)
  - Fits reference signatures to individual samples
  - Generates comprehensive visualizations

- **get_mutational_opportunities.R**: Generates an opportunities matrix for
  SigFit from a BED file. This normalizes mutation counts by trinucleotide
  context frequencies in the target regions.
- **make_SNV_Sigfit_table.strand192.R**: Processes MAF files to create tables
  with mutation context (trinucleotide context) for signature analysis,
  supporting both 96-context and 192-context (strand-specific) analyses.
- **submit_sigfit_jobs.sh**: Submits signature analysis jobs to a computing
  cluster (using LSF).
- **run_all_cohorts.sh**: Orchestrates signature analysis across all tumor
  cohorts defined in the study manifest.



## Requirements



These are the main dependencies:

- R (version 4.2.2)
- R packages including:
  - sigfit
  - BSgenome
  - dplyr
  - Biostrings
  - Cairo
  - GetoptLong

See the `renv.lock` file for the specific versions of R package dependencies used in the various scripts.

## Usage

1. Set up the environment for development on the Sanger HPC 
   ```bash
   source source_me.sh
   ```
   **NB** Otherwise you will need to install the dependencies manually, see the [`Requirements` section](#requirements) for details.

2. Run the analysis for all cohorts:
   ```bash
   bash scripts/run_all_cohorts.sh
   ```

Or for individual cohorts:
   ```bash
   bash scripts/submit_sigfit_jobs.sh -p PROJECT_DIR -g felcat9 -m MAF_FILE -c COSMIC_SIGNATURES -s SIGNAL_SIGNATURES
   ```

## Notes

- Two samples (CATD0653a and CATD294a) in Study 6712 (oral SCC) were moved to
  Study 6711 (cutaneous SCC) due to pathologist correction of tumor site. This
  move was done post-signature analysis, so the sample results remain in Study
  6712.

- The analysis uses both COSMIC signatures and custom "Signal" signatures for
  reference comparisons.