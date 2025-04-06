# fur-cna-plots

A toolset for visualizing copy number alterations from CNVKit segmentation data via cohort-level, sample-level, and clustered heatmap plots.

## Summary
This repository provides scripts & notebooks that transform CNVKit segmentation files into figures for exploratory data analysis and publication.

 - Generates heatmaps that display the net fraction of gains and losses across genomic bins to identify recurrent copy number alterations across different tumor types.
 - Detailed segmentation plots for individual samples, highlighting specific genomic regions with significant copy number changes.
 - Clustered heatmaps focusing on focal (<1Mbp) copy number alterations with optional gene-chromosome annotations incorporated into the figure.
 - Bubble plots that visualize the amplifications and deletions withing genes against different tumor types.

The code for heatmap and segmentations plots, see `scripts/broad_level_cna_heatmap.py` and `scripts/focal_level_cna_heatmap.py`. 

The code for bubble plots, see `notebooks/broad_level_cna_heatmap.ipynb` and `notebooks/focal_level_cna_heatmap.ipynb`

## Installation

The code was developed and tested in Python 3.10, with dependencies listed in `requirements.txt`. To reproduce the environment:

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Usage of `scripts/broad_level_cna_heatmap.py`

```
usage: broad_level_cna_heatmap.py [-h] --mode {cohort,sample,both} --reference REFERENCE [--ignore_prefixes [IGNORE_PREFIXES ...]] [--bin_size BIN_SIZE] [--gain_threshold GAIN_THRESHOLD] [--loss_threshold LOSS_THRESHOLD] [--tumor TUMOR] [--tumor_dirs TUMOR_DIRS]
                                  [--exclude_file EXCLUDE_FILE] [--output_pdf OUTPUT_PDF]

Plot CNV cohort level (formerly heatmap), sample level (formerly segmentation), or both from CNVKit .cns files.

options:
  -h, --help            show this help message and exit
  --mode {cohort,sample,both}
                        Plot mode: 'cohort' for a multi-tumor cohort level plot, 'sample' for per-tumor sample level plot, or 'both' for both. In sample mode, if --tumor is not provided, all tumors are processed.
  --reference REFERENCE
                        Path to the reference FASTA file.
  --ignore_prefixes [IGNORE_PREFIXES ...]
                        Prefixes to ignore (default: AANG KZ MT).
  --bin_size BIN_SIZE   Bin size in base pairs (default: 1000000).
  --gain_threshold GAIN_THRESHOLD
                        Gain threshold (default: 0.585).
  --loss_threshold LOSS_THRESHOLD
                        Loss threshold (default: -0.4).
  --tumor TUMOR         Tumor type to process in sample mode. If not provided, all tumors are processed.
  --tumor_dirs TUMOR_DIRS
                        Optional JSON file mapping tumor types to glob patterns. If not provided, default paths are used.
  --exclude_file EXCLUDE_FILE
                        Path to a file with sample IDs (one per line) to exclude from plots.
  --output_pdf OUTPUT_PDF
                        Output PDF file (for cohort mode) or directory (for sample or both modes with multiple tumors). If not provided, plots are displayed interactively.
```

## Usage of `scripts/focal_level_cna_heatmap.py`

```
usage: focal_level_cna_heatmap.py [-h] [--cohort_json COHORT_JSON] --outdir OUTDIR [--focal_threshold FOCAL_THRESHOLD] [--gain_threshold GAIN_THRESHOLD] [--loss_threshold LOSS_THRESHOLD] [--min_samples MIN_SAMPLES] [--gene_chrom_map GENE_CHROM_MAP] [input_files ...]

Generate clustered heatmaps (vector-based via pcolormesh) of focal (<1Mbp) copy number alterations for one or multiple cohorts. Heatmaps are filtered by log₂FC thresholds and reported only for genes with alterations in a minimum number of samples. Optionally, annotate genes with
chromosome info (and add a colored annotation column) using a gene–chromosome mapping file. If a JSON file mapping cohorts to segmentation file patterns is provided via --cohort_json, a heatmap is generated for each cohort using the corresponding segmentation files. Output files
will be written to the directory specified by --outdir.

positional arguments:
  input_files           List of CNVKit segmentation files (tab-delimited with header) for single-cohort mode.

options:
  -h, --help            show this help message and exit
  --cohort_json COHORT_JSON
                        Optional JSON file mapping cohort IDs to a dict with key 'seg_files' (a list of file patterns).
  --outdir OUTDIR       Directory where output PDF files will be saved.
  --focal_threshold FOCAL_THRESHOLD
                        Maximum segment length (in bp) to be considered focal (default: 1000000).
  --gain_threshold GAIN_THRESHOLD
                        Minimum log₂FC to be considered a gain (default: 0.585).
  --loss_threshold LOSS_THRESHOLD
                        Maximum log₂FC to be considered a loss (default: -0.4).
  --min_samples MIN_SAMPLES
                        Minimum number of samples a gene must have a focal CNA (meeting gain/loss threshold) to be included (default: 1).
  --gene_chrom_map GENE_CHROM_MAP
                        Optional tab-delimited file mapping gene symbols to chromosome, start, and end. If provided, only genes present in this mapping are included and each gene's row is annotated.

```
