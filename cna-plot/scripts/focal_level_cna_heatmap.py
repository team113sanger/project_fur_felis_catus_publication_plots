#!/usr/bin/env python3
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib
from scipy.cluster.hierarchy import linkage, dendrogram
from pathlib import Path

# Unwanted gene symbols to filter out
UNWANTED_GENE_SYMBOLS = {"none", "unk", "incmpl", "cmpl"}

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate clustered heatmaps (vector-based via pcolormesh) of focal (<1Mbp) copy number "
                    "alterations for one or multiple cohorts. Heatmaps are filtered by log₂FC thresholds and reported "
                    "only for genes with alterations in a minimum number of samples. Optionally, annotate genes "
                    "with chromosome info (and add a colored annotation column) using a gene–chromosome mapping file. "
                    "If a JSON file mapping cohorts to segmentation file patterns is provided via --cohort_json, "
                    "a heatmap is generated for each cohort using the corresponding segmentation files. "
                    "Output files will be written to the directory specified by --outdir."
    )
    parser.add_argument(
        "input_files", nargs="*",
        help="List of CNVKit segmentation files (tab-delimited with header) for single-cohort mode."
    )
    parser.add_argument(
        "--cohort_json",
        help="Optional JSON file mapping cohort IDs to a dict with key 'seg_files' (a list of file patterns)."
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Directory where output PDF files will be saved."
    )
    parser.add_argument(
        "--focal_threshold", type=int, default=1000000,
        help="Maximum segment length (in bp) to be considered focal (default: 1000000)."
    )
    parser.add_argument(
        "--gain_threshold", type=float, default=0.585,
        help="Minimum log₂FC to be considered a gain (default: 0.585)."
    )
    parser.add_argument(
        "--loss_threshold", type=float, default=-0.4,
        help="Maximum log₂FC to be considered a loss (default: -0.4)."
    )
    parser.add_argument(
        "--min_samples", type=int, default=1,
        help="Minimum number of samples a gene must have a focal CNA (meeting gain/loss threshold) "
             "to be included (default: 1)."
    )
    parser.add_argument(
        "--gene_chrom_map",
        help="Optional tab-delimited file mapping gene symbols to chromosome, start, and end. "
             "If provided, only genes present in this mapping are included and each gene's row is annotated."
    )
    return parser.parse_args()

def parse_chrom_map(chrom_map_file):
    """
    Parses a tab-delimited file mapping gene symbols to chromosome.
    Uses only lines with 5 columns to build a mapping: { gene_symbol: Chrom }.
    """
    chrom_map = {}
    with open(chrom_map_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) == 5:
                gene_symbol, gene_id, chrom, start, end = fields
                if not gene_symbol.startswith("ENSF"):
                    chrom_map[gene_symbol] = chrom
            elif len(fields) == 4:
                # For this example, ignore lines with only gene IDs.
                pass
    return chrom_map

def process_file(file_path, focal_threshold, gain_threshold, loss_threshold):
    """
    Reads a segmentation file and returns a dictionary mapping each gene
    (from focal segments only, and meeting log₂FC thresholds) to the average log₂ copy number change.
    Only segments with log₂FC >= gain_threshold or log₂FC <= loss_threshold are considered.
    """
    df = pd.read_csv(file_path, sep="\t")
    df['segment_length'] = df['end'] - df['start']
    focal_df = df[df['segment_length'] < focal_threshold]
    
    gene_values = {}
    for _, row in focal_df.iterrows():
        log2_val = row['log2']
        if log2_val >= gain_threshold or log2_val <= loss_threshold:
            genes = [g.strip() for g in row['gene'].split(",") if g.strip()]
            for gene in genes:
                gene_values.setdefault(gene, []).append(log2_val)
    for gene in gene_values:
        gene_values[gene] = np.mean(gene_values[gene])
    return gene_values

def build_heatmap_df(input_files, focal_threshold, gain_threshold, loss_threshold, min_samples, chrom_map=None):
    """
    Builds a DataFrame with genes as rows and samples as columns.
    Each cell contains the average log₂ copy number change for that gene in that sample.
    Only segments meeting the log₂FC thresholds are considered.
    If chrom_map is provided, only genes present in the mapping are included.
    Only genes that appear in at least min_samples samples are retained.
    """
    sample_gene_data = {}
    for file in input_files:
        # Use the file stem (strip extension and any extra parts)
        sample_name = Path(file).stem.split(".")[0]
        gene_data = process_file(file, focal_threshold, gain_threshold, loss_threshold)
        sample_gene_data[sample_name] = gene_data

    if chrom_map is not None:
        all_genes = sorted({
            gene for gene_data in sample_gene_data.values() 
            for gene in gene_data 
            if gene.lower() not in UNWANTED_GENE_SYMBOLS and gene in chrom_map
        })
    else:
        all_genes = sorted({
            gene for gene_data in sample_gene_data.values() 
            for gene in gene_data 
            if gene.lower() not in UNWANTED_GENE_SYMBOLS
        })
    
    df_heat = pd.DataFrame(index=all_genes, columns=sample_gene_data.keys(), dtype=float)
    for sample, gene_dict in sample_gene_data.items():
        for gene, val in gene_dict.items():
            if gene.lower() not in UNWANTED_GENE_SYMBOLS:
                if chrom_map is None or gene in chrom_map:
                    df_heat.at[gene, sample] = val

    df_heat = df_heat[df_heat.count(axis=1) >= min_samples]
    return df_heat

def plot_clustered_heatmap(df_heat, output_file, chrom_map=None):
    """
    Clusters the data and plots a vector-based heatmap using pcolormesh.
    If a chromosome map is provided, a narrow annotation column is added as the first column
    of the heatmap. A vertical separator is drawn between the annotation and the data.
    """
    # Replace missing values with 0
    data = df_heat.fillna(0).values
    num_rows, num_cols = data.shape

    # Hierarchical clustering on rows and columns
    Z_rows = linkage(data, method='average', metric='euclidean')
    Z_cols = linkage(data.T, method='average', metric='euclidean')
    dendro_rows = dendrogram(Z_rows, orientation='left', no_plot=True)
    dendro_cols = dendrogram(Z_cols, no_plot=True)
    row_order = dendro_rows['leaves']
    col_order = dendro_cols['leaves']
    
    df_ordered = df_heat.iloc[row_order, :].iloc[:, col_order]
    data_ordered = df_ordered.values  # shape: (num_rows, num_cols)

    # Adjust figure size dynamically
    fig_width = max(12, 0.3 * num_cols + 4)
    fig_height = max(10, 0.3 * num_rows + 2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    
    # Prepare annotation column if chrom_map is provided
    if chrom_map is not None:
        unique_chroms = sorted(set(chrom_map.get(gene) for gene in df_ordered.index))
        chrom_to_idx = {chrom: i for i, chrom in enumerate(unique_chroms)}
        annotation_data = np.array([chrom_to_idx[chrom_map.get(gene)] for gene in df_ordered.index]).reshape(num_rows, 1)
        x_ann_edges = np.array([0, 1])
        x_data_edges = np.arange(1, num_cols + 2)
    else:
        x_data_edges = np.arange(0, num_cols + 1)
    
    y_edges = np.arange(num_rows + 1)
    
    # Plot the data heatmap (copy-number values)
    mesh_data = ax.pcolormesh(x_data_edges, y_edges, data_ordered, cmap=cm.RdBu_r,
                              norm=mcolors.Normalize(vmin=-np.nanmax(np.abs(data_ordered)),
                                                     vmax=np.nanmax(np.abs(data_ordered))),
                              shading='auto')
    
    # If annotation is provided, plot the annotation column and a vertical separator
    if chrom_map is not None:
        cmap_chrom_base = matplotlib.colormaps.get_cmap("tab20")
        cmap_chrom = cmap_chrom_base.resampled(len(unique_chroms))
        mesh_ann = ax.pcolormesh(x_ann_edges, y_edges, annotation_data, cmap=cmap_chrom,
                                 shading='auto')
        ax.axvline(x=1, color='black', linewidth=2)
    
    ax.invert_yaxis()
    y_tick_positions = np.arange(num_rows) + 0.5
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels(df_ordered.index, fontsize=8)
    
    if chrom_map is not None:
        x_tick_positions = np.arange(1, num_cols + 1) + 0.5
        ax.set_xticks(x_tick_positions)
        ax.set_xticklabels(df_ordered.columns, rotation=90, fontsize=8)
        ax.set_xlim(0, num_cols + 1)
    else:
        x_tick_positions = np.arange(num_cols) + 0.5
        ax.set_xticks(x_tick_positions)
        ax.set_xticklabels(df_ordered.columns, rotation=90, fontsize=8)
        ax.set_xlim(0, num_cols)
    
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    cbar = fig.colorbar(mesh_data, ax=ax)
    cbar.set_label("log₂ copy number change", fontsize=10)
    ax.set_title("Clustered Heatmap of Focal (<1Mbp) Copy Number Alterations", fontsize=12)
    
    if chrom_map is not None:
        from matplotlib.patches import Patch
        legend_handles = [Patch(facecolor=cmap_chrom(i), label=chrom) for i, chrom in enumerate(unique_chroms)]
        ax.legend(handles=legend_handles, title="Chromosome", loc="upper left", bbox_to_anchor=(1.01, 1))
    
    plt.savefig(output_file, format='pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"Clustered heatmap saved to {output_file}")

def main():
    args = parse_arguments()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Parse gene-chromosome map if provided
    chrom_map = None
    if args.gene_chrom_map:
        chrom_map = parse_chrom_map(args.gene_chrom_map)
    
    # Multiple-cohort mode if a JSON file is provided
    if args.cohort_json:
        with open(args.cohort_json, "r") as f:
            cohort_dict = json.load(f)
        for cohort_id, info in cohort_dict.items():
            seg_files = info.get("seg_files", [])
            files = []
            for pattern in seg_files:
                files.extend([str(p) for p in Path().glob(pattern)])
            if not files:
                print(f"No segmentation files found for cohort {cohort_id}. Skipping.")
                continue
            df_heat = build_heatmap_df(
                files,
                args.focal_threshold,
                args.gain_threshold,
                args.loss_threshold,
                args.min_samples,
                chrom_map=chrom_map
            )
            if df_heat.empty:
                print(f"No data available to plot for cohort {cohort_id} after filtering.")
                continue
            output_file = outdir / f"{cohort_id}_heatmap.pdf"
            plot_clustered_heatmap(df_heat, str(output_file), chrom_map=chrom_map)
    else:
        # Single-cohort mode, process input_files directly.
        if not args.input_files:
            print("No input files provided.")
            return
        files = []
        for pattern in args.input_files:
            files.extend([str(p) for p in Path().glob(pattern)])
        if not files:
            print("No segmentation files found from provided patterns.")
            return
        df_heat = build_heatmap_df(
            files,
            args.focal_threshold,
            args.gain_threshold,
            args.loss_threshold,
            args.min_samples,
            chrom_map=chrom_map
        )
        if df_heat.empty:
            print("No data available to plot after filtering.")
            return
        output_file = outdir / "heatmap.pdf"
        plot_clustered_heatmap(df_heat, str(output_file), chrom_map=chrom_map)

if __name__ == '__main__':
    main()
