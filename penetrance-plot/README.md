# Penetrance Plots

This repository contains code for generating and visualizing genome-wide DNA copy number alteration penetrance across cohorts of samples, based on CNVkit segmentation files.

## Repository Structure

```
├── data/
│   ├── FUR_cat.samples_to_exclude.nucleotide_variants.txt
│   └── cnvkit_segmentation_files/
│       └── <tumor_type>/<sample>.cns   # CNVkit segmentation files (.cns)
├── penetrance_plots/                    # Output PDF plots
├── penetrance_plots.ipynb               # Jupyter notebook with analysis functions
├── requirements.txt                     # Python dependencies
└── README.md                            # This file
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
jupyter notebook penetrance_plots.ipynb
```

## Requirements

- Python 3.8+
- pandas
- matplotlib

All required packages are listed in `requirements.txt`.

