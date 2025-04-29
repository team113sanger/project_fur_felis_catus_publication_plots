# Oncoprints

This repository contains code for generating oncoprint plots to help visualise
both somatic mutations and copy number alterations across a given cancer cohort.
Somatic mutation data is inferred from MAF files, while copy number alteration
data is inferred from CNVKit genemetrics summary CSVs.

## Repository Structure

```
├── LICENSE
├── README.md
├── data
│   ├── FUR_cat.samples_to_exclude.nucleotide_variants.txt         # List of samples to exclude from analysis
│   ├── Felis_catus_9.0.gene_chromosome_mapping.txt                # Gene chromosome mapping file
│   ├── breakdown
│   │   ├── 6555_Lung_carcinoma_oncoprint_breakdown.csv
│   │   ├── <cohort + cancer type>_oncoprint_breakdown.csv
│   │   ├── ...
│   │   ├── B_cell_lymphoma_breakdown.csv
│   │   └── T_cell_lymphoma_breakdown.csv
│   ├── genemetrics_study_summary                                  # CNVKit genemetrics summary CSVs
│   │   ├── 6555_Lung_carcinoma.genemetrics_study_summary.csv
│   │   ├── <cohort + cancer type>_genemetrics_study_summary.csv
│   │   └── ...
│   └── maf_files                                                  # MAF files for somatic mutations
│       ├── keepPA_vaf_size_filt_matched_6555_2711.finalised.maf
│       ├── *_<cohort>_*.hotspot_mutations.maf
│       └── ...
├── notebooks
│   └── oncoprints.ipynb                                           # Notebook that creates result plots
├── results                                                        # Oncoprint result plots
│   ├── 6555_Lung_carcinoma_oncoprint.pdf
│   ├── <cohort + cancer typer>_oncoprint.pdf
│   ├── ...
│   ├── B_cell_lymphoma_oncoprint.pdf
│   └── T_cell_lymphoma_oncoprint.pdf
└── scripts
    └── oncoprints.sh
```

## Installation

1. Clone the repository:
   ```bash
   git clone git@gitlab.internal.sanger.ac.uk:DERMATLAS/fur/fur_penetrance_plots.git
   cd fur_penetrance_plots
   ```
2. Create and activate a virtual environment:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```


## Usage
To run the notebook, use the following command:

```bash
jupyter notebook notebooks/oncoprints.ipynb
```

## Requirements

- Python 3.8+
- pandas
- numpy
- matplotlib
- scipy
- jupyterlab

All required packages are listed in `requirements.txt`.
