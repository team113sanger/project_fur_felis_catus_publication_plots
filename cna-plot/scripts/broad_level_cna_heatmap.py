#!/usr/bin/env python

import argparse
import glob
import json
import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, dendrogram

def get_chromosome_lengths(reference_fasta, ignore_prefixes=None):
    """
    Reads a reference FASTA file and returns a dictionary of chromosome lengths.
    Sequences whose IDs (after removing an optional "chr" prefix) that start with any of
    the ignore_prefixes are skipped.
    """
    if ignore_prefixes is None:
        ignore_prefixes = []
    chrom_lengths = {}
    for record in SeqIO.parse(reference_fasta, "fasta"):
        chrom_id = record.id
        if chrom_id.lower().startswith("chr"):
            chrom_id = chrom_id[3:]
        if any(chrom_id.startswith(prefix) for prefix in ignore_prefixes):
            continue
        chrom_lengths[chrom_id] = len(record.seq)
    return chrom_lengths

def plot_heatmap(reference_fasta, tumor_dirs, ignore_prefixes, bin_size, gain_threshold, loss_threshold,
                 exclude_samples=None, output_pdf=None):
    # This function produces the cohort level plot.
    chrom_lengths = get_chromosome_lengths(reference_fasta, ignore_prefixes)
    bins = []
    bin_labels = []
    for chrom in sorted(chrom_lengths.keys()):
        length = chrom_lengths[chrom]
        for start in range(0, length, bin_size):
            end = min(start + bin_size, length)
            bins.append((chrom, start, end))
            bin_labels.append(f"{chrom}:{start}-{end}")
    
    tumor_heatmap = {}
    for tumor, pattern in tumor_dirs.items():
        files = glob.glob(pattern)
        if exclude_samples:
            files = [f for f in files if f.split('/')[-1].split('.')[0] not in exclude_samples]
        print(f"Processing {tumor}: {len(files)} samples found.")
        if len(files) == 0:
            print(f"No files found for tumor {tumor} after filtering. Skipping.")
            continue
        
        sample_matrix = np.zeros((len(files), len(bins)))
        for i, file in enumerate(files):
            df = pd.read_csv(file, sep='\t')
            df['chromosome'] = df['chromosome'].astype(str).str.replace('chr', '', regex=False)
            for j, (chrom, bin_start, bin_end) in enumerate(bins):
                segs = df[(df['chromosome'] == chrom) &
                          (df['end'] >= bin_start) &
                          (df['start'] <= bin_end)]
                if not segs.empty:
                    mean_log2 = segs['log2'].mean()
                    if mean_log2 > gain_threshold:
                        sample_matrix[i, j] = 1
                    elif mean_log2 < loss_threshold:
                        sample_matrix[i, j] = -1
                    else:
                        sample_matrix[i, j] = 0
                else:
                    sample_matrix[i, j] = 0
        prop_gain = np.sum(sample_matrix == 1, axis=0) / sample_matrix.shape[0]
        prop_loss = np.sum(sample_matrix == -1, axis=0) / sample_matrix.shape[0]
        net_score = prop_gain - prop_loss
        tumor_heatmap[tumor] = net_score

    heatmap_df = pd.DataFrame(tumor_heatmap).T
    heatmap_df.columns = bin_labels
    heatmap_df = heatmap_df.reindex(columns=bin_labels)
    
    chrom_order = sorted(chrom_lengths.keys())
    cum_offsets = {}
    current_offset = 0
    for chrom in chrom_order:
        cum_offsets[chrom] = current_offset
        current_offset += chrom_lengths[chrom]
    total_genome_length = current_offset

    genome_bins = []
    for chrom in sorted(chrom_lengths.keys()):
        length = chrom_lengths[chrom]
        for start in range(0, length, bin_size):
            end = min(start + bin_size, length)
            genome_bins.append((chrom, start, end))
    bin_positions = []
    for (chrom, start, end) in genome_bins:
        global_start = cum_offsets[chrom] + start
        global_end = cum_offsets[chrom] + end
        bin_positions.append((global_start, global_end))
    assert len(heatmap_df.columns) == len(bin_positions), "Mismatch in bin columns vs bin_positions."

    tumor_names = heatmap_df.index.tolist()
    Z = linkage(heatmap_df.values, method='average', metric='euclidean')
    dendro = dendrogram(Z, orientation='left', labels=tumor_names, no_plot=True)
    row_order = dendro['leaves']
    heatmap_df = heatmap_df.iloc[row_order, :]

    fig = plt.figure(figsize=(18, max(4, len(tumor_names)*0.6)))
    fig.subplots_adjust(left=0.3, wspace=0.05)
    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[1, 4])
    ax_dendro = fig.add_subplot(gs[0])
    ax_cohort = fig.add_subplot(gs[1])
    
    dendrogram(Z, orientation='left', labels=None, ax=ax_dendro, color_threshold=0)
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])

    norm = mcolors.Normalize(vmin=-0.4, vmax=0.4, clip=True)
    cmap = cm.RdBu_r

    for i, tumor_name in enumerate(heatmap_df.index):
        for j, col in enumerate(heatmap_df.columns):
            net_score = heatmap_df.iloc[i, j]
            global_start, global_end = bin_positions[j]
            width = global_end - global_start
            rect_color = cmap(norm(net_score))
            rect = plt.Rectangle((global_start, i), width, 0.8, color=rect_color)
            ax_cohort.add_patch(rect)
        ax_cohort.text(-total_genome_length*0.005, i + 0.4, tumor_name,
                        va='center', ha='right', fontsize=9)
    
    ax_cohort.set_yticks([])
    ax_cohort.set_ylim(0, len(heatmap_df))
    ax_cohort.set_xlim(0, total_genome_length)
    ax_cohort.set_xlabel("Genomic Coordinate")
    ax_cohort.set_title("Cohort Level Plot: Net Fraction of Gains/Losses per Tumor Type")
    
    for chrom in chrom_order:
        offset = cum_offsets[chrom]
        ax_cohort.axvline(offset, color='black', linewidth=1)
    ax_cohort.axvline(total_genome_length, color='black', linewidth=1)
    
    xticks = []
    xtick_labels = []
    for chrom in chrom_order:
        start = cum_offsets[chrom]
        end = start + chrom_lengths[chrom]
        xticks.append((start + end) / 2)
        xtick_labels.append(chrom)
    ax_cohort.set_xticks(xticks)
    ax_cohort.set_xticklabels(xtick_labels, rotation=90, fontsize=8)
    
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    fig.colorbar(sm, ax=ax_cohort, orientation='vertical', pad=0.02).set_label("Net (gain - loss) fraction")
    
    if output_pdf:
        plt.tight_layout()
        plt.savefig(output_pdf, format='pdf')
        print(f"Saved cohort level plot to {output_pdf}")
        plt.close()
    else:
        plt.show()

def plot_segmentation(reference_fasta, tumor_dirs, ignore_prefixes, bin_size, gain_threshold, loss_threshold,
                      tumor, exclude_samples=None, output_pdf=None):
    # This function produces the sample level plot.
    chrom_lengths = get_chromosome_lengths(reference_fasta, ignore_prefixes)
    chrom_order = sorted(chrom_lengths.keys())
    cum_offsets = {}
    current_offset = 0
    for chrom in chrom_order:
        cum_offsets[chrom] = current_offset
        current_offset += chrom_lengths[chrom]
    total_genome_length = current_offset

    if tumor not in tumor_dirs:
        print(f"Tumor type '{tumor}' not found in tumor directories mapping.")
        return
    pattern = tumor_dirs[tumor]
    files = glob.glob(pattern)
    if exclude_samples:
        files = [f for f in files if f.split('/')[-1].split('.')[0] not in exclude_samples]
    print(f"Processing {tumor}: {len(files)} samples found.")
    if len(files) == 0:
        print("No files found. Exiting.")
        return

    bins = np.arange(0, total_genome_length + bin_size, bin_size)
    n_bins = len(bins) - 1

    def create_profile_vector(df):
        profile = np.zeros(n_bins)
        for _, row in df.iterrows():
            chrom = row['chromosome']
            if chrom not in cum_offsets:
                continue
            seg_start = row['start']
            seg_end = row['end']
            global_start = cum_offsets[chrom] + seg_start
            global_end = cum_offsets[chrom] + seg_end
            start_bin = np.searchsorted(bins, global_start, side='right') - 1
            end_bin = np.searchsorted(bins, global_end, side='right') - 1
            profile[start_bin:end_bin+1] = row['log2']
        return profile

    profiles = []
    sample_names = []
    for file in files:
        df = pd.read_csv(file, sep='\t')
        df['chromosome'] = df['chromosome'].astype(str).str.replace('chr', '', regex=False)
        profiles.append(create_profile_vector(df))
        sample_names.append(file.split('/')[-1].split('.')[0])
    profiles = np.array(profiles)

    Z = linkage(profiles, method='average', metric='euclidean')
    dendro = dendrogram(Z, orientation='left', labels=sample_names, no_plot=True)
    order = dendro['leaves']

    all_log2_values = []
    for file in files:
        df = pd.read_csv(file, sep='\t')
        for val in df['log2'].values:
            if val <= loss_threshold or val >= gain_threshold:
                all_log2_values.append(val)
    tumor_min = min(all_log2_values) if all_log2_values else loss_threshold
    tumor_max = max(all_log2_values) if all_log2_values else gain_threshold
    norm = mcolors.Normalize(vmin=tumor_min, vmax=tumor_max)
    cmap = cm.RdBu_r

    fig = plt.figure(figsize=(18, max(4, len(files)*0.5)))
    fig.subplots_adjust(left=0.2, wspace=0.05)
    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[1, 4])
    ax_dendro = fig.add_subplot(gs[0])
    ax_sample = fig.add_subplot(gs[1])
    
    dendrogram(Z, orientation='left', labels=None, ax=ax_dendro, color_threshold=0)
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])

    for i, idx in enumerate(order):
        file = files[idx]
        df = pd.read_csv(file, sep='\t')
        df['chromosome'] = df['chromosome'].astype(str).str.replace('chr', '', regex=False)
        for _, row in df.iterrows():
            chrom = row['chromosome']
            if chrom not in cum_offsets:
                continue
            seg_start = row['start']
            seg_end = row['end']
            global_start = cum_offsets[chrom] + seg_start
            global_end = cum_offsets[chrom] + seg_end
            seg_width = global_end - global_start
            log2_val = row['log2']
            if loss_threshold < log2_val < gain_threshold:
                continue
            rect_color = cmap(norm(log2_val))
            rect = patches.Rectangle((global_start, i), seg_width, 0.8, color=rect_color)
            ax_sample.add_patch(rect)
        ax_sample.text(-total_genome_length*0.005, i+0.4, sample_names[idx],
                   va='center', ha='right', fontsize=8)

    ax_sample.set_yticks([])
    ax_sample.set_ylim(0, len(files))
    ax_sample.set_xlim(0, total_genome_length)
    ax_sample.set_xlabel("Genomic Coordinate")
    ax_sample.set_title(f"Sample Level Plot: Copy Number Segments: {tumor}")
    
    for chrom in chrom_order:
        offset = cum_offsets[chrom]
        ax_sample.axvline(offset, color='black', linewidth=1)
    ax_sample.axvline(total_genome_length, color='black', linewidth=1)
    
    xticks = []
    xtick_labels = []
    for chrom in chrom_order:
        start = cum_offsets[chrom]
        end = start + chrom_lengths[chrom]
        xticks.append((start + end) / 2)
        xtick_labels.append(chrom)
    ax_sample.set_xticks(xticks)
    ax_sample.set_xticklabels(xtick_labels, rotation=90, fontsize=8)
    
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax_sample, orientation='vertical', pad=0.02).set_label('log₂ copy number change')
    
    if output_pdf:
        plt.tight_layout()
        plt.savefig(output_pdf, format='pdf')
        print(f"Saved sample level plot for {tumor} to {output_pdf}")
        plt.close()
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(
        description="Plot CNV cohort level (formerly heatmap), sample level (formerly segmentation), or both from CNVKit .cns files."
    )
    parser.add_argument("--mode", choices=["cohort", "sample", "both"], required=True,
                        help="Plot mode: 'cohort' for a multi-tumor cohort level plot, 'sample' for per-tumor sample level plot, or 'both' for both. "
                             "In sample mode, if --tumor is not provided, all tumors are processed.")
    parser.add_argument("--reference", required=True,
                        help="Path to the reference FASTA file.")
    parser.add_argument("--ignore_prefixes", nargs='*', default=["AANG", "KZ", "MT"],
                        help="Prefixes to ignore (default: AANG KZ MT).")
    parser.add_argument("--bin_size", type=int, default=1000000,
                        help="Bin size in base pairs (default: 1000000).")
    parser.add_argument("--gain_threshold", type=float, default=0.585,
                        help="Gain threshold (default: 0.585).")
    parser.add_argument("--loss_threshold", type=float, default=-0.4,
                        help="Loss threshold (default: -0.4).")
    parser.add_argument("--tumor", type=str,
                        help="Tumor type to process in sample mode. If not provided, all tumors are processed.")
    parser.add_argument("--tumor_dirs", type=str, default="",
                        help="Optional JSON file mapping tumor types to glob patterns. If not provided, default paths are used.")
    parser.add_argument("--exclude_file", type=str, default="",
                        help="Path to a file with sample IDs (one per line) to exclude from plots.")
    parser.add_argument("--output_pdf", type=str, default="",
                        help="Output PDF file (for cohort mode) or directory (for sample or both modes with multiple tumors). "
                             "If not provided, plots are displayed interactively.")

    args = parser.parse_args()

    exclude_samples = set()
    if args.exclude_file:
        with open(args.exclude_file, 'r') as f:
            exclude_samples = {line.strip() for line in f if line.strip()}

    default_tumor_dirs = cancer_files = {
        '6555_Lung_carcinoma': 'data/CNVKit-segmentation/6555_Lung_carcinoma/*.call.median_centred.cns',
        '6712_Oral_SCC': 'data/CNVKit-segmentation/6712_Oral_SCC/*.call.median_centred.cns',
        '6841_Meningioma': 'data/CNVKit-segmentation/6841_Meningioma/*.call.median_centred.cns',
        '6945_Cholangiocarcinoma': 'data/CNVKit-segmentation/6945_Cholangiocarcinoma/*.call.median_centred.cns',
        '6982_Lymphoma': 'data/CNVKit-segmentation/6982_Lymphoma/*.call.median_centred.cns',
        '7040_BCC': 'data/CNVKit-segmentation/7040_BCC/*.call.median_centred.cns',
        '7098_Glioma': 'data/CNVKit-segmentation/7098_Glioma/*.call.median_centred.cns',
        '6711_Cutaneous_SCC': 'data/CNVKit-segmentation/6711_Cutaneous_SCC/*.call.median_centred.cns',
        '6713_Cutaneous_MCT': 'data/CNVKit-segmentation/6713_Cutaneous_MCT/*.call.median_centred.cns',
        '6864_Pancreatic_carcinoma': 'data/CNVKit-segmentation/6864_Pancreatic_carcinoma/*.call.median_centred.cns',
        '6973_OSA': 'data/CNVKit-segmentation/6973_OSA/*.call.median_centred.cns',
        '6990_Mammary_carcinoma': 'data/CNVKit-segmentation/6990_Mammary_carcinoma/*.call.median_centred.cns',
        '7097_CRC': 'data/CNVKit-segmentation/7097_CRC/*.call.median_centred.cns'
    }
    
    if args.tumor_dirs:
        with open(args.tumor_dirs, 'r') as f:
            tumor_dirs = json.load(f)
    else:
        tumor_dirs = default_tumor_dirs
    
    # Check the directories in the tumor_dirs mapping
    for _, pattern in tumor_dirs.items():
        # Each value has a trailing glob pattern, so we need to remove it
        directory = pattern.rsplit('/', 1)[0]
        if not os.path.isdir(directory):
            print(f"Directory {directory} does not exist. Please check the tumor_dirs mapping.")
            sys.exit(1)

    output_pdf = args.output_pdf if args.output_pdf else None

    if args.mode == "cohort":
        plot_heatmap(args.reference, tumor_dirs, args.ignore_prefixes,
                     args.bin_size, args.gain_threshold, args.loss_threshold,
                     exclude_samples=exclude_samples, output_pdf=output_pdf)
    elif args.mode == "sample":
        if args.tumor:
            plot_segmentation(args.reference, tumor_dirs, args.ignore_prefixes,
                              args.bin_size, args.gain_threshold, args.loss_threshold,
                              args.tumor, exclude_samples=exclude_samples, output_pdf=output_pdf)
        else:
            if output_pdf:
                if not os.path.isdir(output_pdf):
                    os.makedirs(output_pdf, exist_ok=True)
            for tumor in tumor_dirs.keys():
                tumor_output = os.path.join(output_pdf, f"segmentation_{tumor}.pdf") if output_pdf else None
                plot_segmentation(args.reference, tumor_dirs, args.ignore_prefixes,
                                  args.bin_size, args.gain_threshold, args.loss_threshold,
                                  tumor, exclude_samples=exclude_samples, output_pdf=tumor_output)
    elif args.mode == "both":
        # In both mode, output_pdf is treated as a directory.
        if output_pdf:
            if not os.path.isdir(output_pdf):
                os.makedirs(output_pdf, exist_ok=True)
            heatmap_pdf = os.path.join(output_pdf, "cohort_level_plot.pdf")
            plot_heatmap(args.reference, tumor_dirs, args.ignore_prefixes,
                         args.bin_size, args.gain_threshold, args.loss_threshold,
                         exclude_samples=exclude_samples, output_pdf=heatmap_pdf)
        else:
            plot_heatmap(args.reference, tumor_dirs, args.ignore_prefixes,
                         args.bin_size, args.gain_threshold, args.loss_threshold,
                         exclude_samples=exclude_samples, output_pdf=None)
        if args.tumor:
            seg_pdf = os.path.join(output_pdf, f"sample_level_plot_{args.tumor}.pdf") if output_pdf else None
            plot_segmentation(args.reference, tumor_dirs, args.ignore_prefixes,
                              args.bin_size, args.gain_threshold, args.loss_threshold,
                              args.tumor, exclude_samples=exclude_samples, output_pdf=seg_pdf)
        else:
            if output_pdf:
                if not os.path.isdir(output_pdf):
                    os.makedirs(output_pdf, exist_ok=True)
            for tumor in tumor_dirs.keys():
                tumor_output = os.path.join(output_pdf, f"sample_level_plot_{tumor}.pdf") if output_pdf else None
                plot_segmentation(args.reference, tumor_dirs, args.ignore_prefixes,
                                  args.bin_size, args.gain_threshold, args.loss_threshold,
                                  tumor, exclude_samples=exclude_samples, output_pdf=tumor_output)

if __name__ == "__main__":
    main()