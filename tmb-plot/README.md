# FUR Cat Tumour Mutational Burden Analysis

This repository contains a Jupyter notebook for computing and visualising tumor mutational burden (TMB) across the multiple FUR feline cancer cohorts. This notebook utlisises plotting code (in particular the `plotTMB` and `prepend` functions) from [TMB_plotter](https://github.com/AlexandrovLab/TMB_plotter).

This notebook also contains TMB calculations for samples with/without the SBS7 mutational signature.

## Requirements
All Python dependencies for this repository can be found in `requirements.txt` and can be installed using `pip install -r requirements.txt`

## Usage
To run the notebook, use the following command:
```bash
jupyter notebook tmb_analysis.ipynb
```

## Outputs
* `final_dataframe.csv` - filtered dataframe with the columns: Types,Mut_burden,log10BURDENpMB
* `TMB_plot.png` - Output TMB plot of per-cohort TMB distributions