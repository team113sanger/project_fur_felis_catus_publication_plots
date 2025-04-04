# Creation of cohort summaries with LOF variants and LOF + missense_variants in all genes and NHS genes

The objective is to generate a set of tables with the collated summary of all variants with a predicted Loss of Function (LOF) or LOF and Missense found in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label. Finally it proceeds to generate plots regarding the VAF distribution of the variants, per cohort and per mutation type.

## Get the MAF files with LOF and LOF + missense variants on genes by the NHS germline cancer risk testing

```bash 
PROJECT_DIR=/lustre/scratch127/casm/team113da/projects/fur_germline

# Run the script
cd ${PROJECT_DIR:?unset}

bash ${PROJECT_DIR:?unset}/scripts/run_risk_analysis.nhs.wmiss.sh

```

### Transform MAF files to XLSX

```bash
PROJECT_DIR=/lustre/scratch127/casm/team113da/projects/fur_germline

# Run the script
cd ${PROJECT_DIR:?unset}

for STUDY in `cut -f 2 ${PROJECT_DIR}/metadata/studies_reformat.tsv`; do
  cd ${PROJECT_DIR}/${STUDY}/risk_analysis/
  MAF=${PROJECT_DIR}/${STUDY}/risk_analysis/keepPA_vaf_size_filt_germline_${STUDY}.nhs.lof.wmiss.maf
  Rscript ${PROJECT_DIR:?unset}/scripts/MAF/tsv2xlsx.R ${MAF}
  cd ${PROJECT_DIR}
done

```
## Get the MAF files with LOF and LOF + missense variants on any genes

```bash 
PROJECT_DIR=/lustre/scratch127/casm/team113da/projects/fur_germline

# Run the script
cd ${PROJECT_DIR:?unset}

bash ${PROJECT_DIR:?unset}/scripts/run_risk_analysis.wmiss.sh

```

### Transform MAF files to XLSX

```bash
PROJECT_DIR=/lustre/scratch127/casm/team113da/projects/fur_germline

# Run the script
cd ${PROJECT_DIR:?unset}

for STUDY in `cut -f 2 ${PROJECT_DIR}/metadata/studies_reformat.tsv`; do
  cd ${PROJECT_DIR}/${STUDY}/risk_analysis/
  MAF=${PROJECT_DIR}/${STUDY}/risk_analysis/keepPA_vaf_size_filt_germline_${STUDY}.nhs.lof.wmiss.maf
  Rscript ${PROJECT_DIR:?unset}/scripts/MAF/tsv2xlsx.R ${MAF:?unset}
  cd ${PROJECT_DIR}
done

```

## Generate the cohort summaries with Tumour type label and variant counting 

```bash
PROJECT_DIR=/lustre/scratch127/casm/team113da/projects/fur_germline

# Run the script
cd ${PROJECT_DIR:?unset}

source ${PROJECT_DIR:?unset}/scripts/paper_tables/source_me.sh
Rscript ${PROJECT_DIR:?unset}/scripts/paper_tables/cohort_summary_tables.R

```

### Generate the XLSX files with the cohort summaries

```bash
PROJECT_DIR=/lustre/scratch127/casm/team113da/projects/fur_germline

# Run the script
cd ${PROJECT_DIR:?unset}

for MAF in ` ls -1 ${PROJECT_DIR}/cohort_analysis/cohort_tables/*.maf`; do
  cd ${PROJECT_DIR}/cohort_analysis/cohort_tables/
  Rscript ${PROJECT_DIR:?unset}/scripts/MAF/tsv2xlsx.R ${MAF:?unset}
  cd ${PROJECT_DIR}
done

```

## Output files

The list of files generated are following:

- **List of tables and `MAF` files generated**
  - **cohort_summary_numtab.tsv**: Table with the number of mutations found in each tumour type per consequence classification.
  - **FUR_keepPA_vaf_size_filt_germline_all_genes.risk.maf**: Table with the list of all variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types' cohorts and the addition of the tumour type label.
  - **FUR_keepPA_vaf_size_filt_germline_all_genes.risk.vafn0.25filt.maf**: Table with the list of all variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types' cohorts and the addition of the tumour type label. The variants were filtered by VAF > 0.25.
  - **FUR_keepPA_vaf_size_filt_germline_all_genes.risk.wmiss.maf**: Table with the list of all variants with a predicted Loss of Function (LOF) and Missense found in any of the samples across all of multiple Tumour types' cohorts and the addition of the tumour type label.
  - **FUR_keepPA_vaf_size_filt_germline_all_genes.risk.wmiss.vafn0.25filt.maf**: Table with the list of all variants with a predicted Loss of Function (LOF) and Missense found in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label. The variants were filtered by VAF > 0.25.
  - **FUR_keepPA_vaf_size_filt_germline.nhs.lof.maf**: Table with the list of all variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within all the genes whose human orthologue is present within the set of genes used by NHS germline cancer predisposition testing.
  - **FUR_keepPA_vaf_size_filt_germline.nhs.lof.vafn0.25filt.maf**: Table with the list of all variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within all the genes whose human orthologue is present within the set of genes used by NHS germline cancer predisposition testing.The variants were filtered by VAF > 0.25.
  - **FUR_keepPA_vaf_size_filt_germline.nhs.lof.wmiss.maf**: Table with the list of all variants with a predicted Loss of Function (LOF) and Missense found in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within NHS germline cancer predisposition orthologs.
  - **FUR_keepPA_vaf_size_filt_germline.nhs.lof.wmiss.vafn0.25filt.maf**: Table with the list of all variants with a predicted Loss of Function (LOF) and Missense found in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within NHS germline cancer predisposition orthologs. The variants are filtered by VAF > 0.25.

- **FUR_keepPA_vaf_size_filt_germline.nhs.lof.vafn0.25filt.lndepth.maf**: Table with the list of all variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within all the genes whose human orthologue is present within the set of genes used by NHS germline cancer predisposition testing.The variants were filtered by VAF > 0.25 and n_depth >6 reads.
  

- **Folders with the list of plots generated for each aggregated maf file**
  - **cohort_maf_lof**: Folder with the plots of the VAF distribution of the variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label.
  - **cohort_maf_lofwmis**: Folder with the plots of the VAF distribution of the variants with a predicted Loss of Function (LOF) and Missense found in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label.
  - **cohort_maf_nhs_lof**: Folder with the plots of the VAF distribution of the variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within all the genes whose human orthologue is present within the set of genes used by NHS germline cancer predisposition testing.
  - **cohort_maf_nhs_lof**: Folder with the plots of the VAF distribution of the variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within all the genes whose human orthologue is present within the set of genes used by NHS germline cancer predisposition testing.
  - **cohort_maf_nhs_lof_vaf_0.25filt**: Folder with the plots of the VAF distribution of the variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within all the genes whose human orthologue is present within the set of genes used by NHS germline cancer predisposition testing. The variants were filtered by VAF > 0.25 and minimum depth.
  - **cohort_maf_nhs_lof_vaf_0.25filt_lndepthf**: Folder with the plots of the VAF distribution of the variants with a predicted Loss of Function (LOF) in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within all the genes whose human orthologue is present within the set of genes used by NHS germline cancer predisposition testing. The variants were filtered by VAF > 0.25 and minimum coverage depth filter.
  
  - **cohort_maf_nhs_lofwmis**: Folder with the plots of the VAF distribution of the variants with a predicted Loss of Function (LOF) and Missense found in any of the samples across all of multiple Tumour types cohorts and the addition of the tumour type label within NHS germline cancer predisposition orthologs.

```bash
cd  ${PROJECT_DIR}/cohort_analysis/cohort_tables/
tree -L 3
.
├── cohort_summary_numtab.tsv
├── FUR_keepPA_vaf_size_filt_germline_all_genes.risk.maf
├── FUR_keepPA_vaf_size_filt_germline_all_genes.risk.maf.xlsx
├── FUR_keepPA_vaf_size_filt_germline_all_genes.risk.vafn0.25filt.maf
├── FUR_keepPA_vaf_size_filt_germline_all_genes.risk.vafn0.25filt.maf.xlsx
├── FUR_keepPA_vaf_size_filt_germline_all_genes.risk.wmiss.maf
├── FUR_keepPA_vaf_size_filt_germline_all_genes.risk.wmiss.maf.xlsx
├── FUR_keepPA_vaf_size_filt_germline_all_genes.risk.wmiss.vafn0.25filt.maf
├── FUR_keepPA_vaf_size_filt_germline_all_genes.risk.wmiss.vafn0.25filt.maf.xlsx
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.maf
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.maf.xlsx
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.vafn0.25filt.lndepthf.maf
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.vafn0.25filt.lndepthf.maf.xlsx
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.vafn0.25filt.maf
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.vafn0.25filt.maf.xlsx
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.wmiss.maf
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.wmiss.maf.xlsx
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.wmiss.vafn0.25filt.maf
├── FUR_keepPA_vaf_size_filt_germline.nhs.lof.wmiss.vafn0.25filt.maf.xlsx
└── plots
    ├── cohort_maf_lof
    │   ├── cohort_maf_lof_VAF_norm_histogram_per_mut_type.pdf
    │   └── cohort_maf_lof_VAF_norm_histogram_per_Tumour_type.pdf
    ├── cohort_maf_lof_vaf_0.25filt
    │   ├── cohort_maf_lof_vaf_0.25filt_VAF_norm_histogram_per_mut_type.pdf
    │   └── cohort_maf_lof_vaf_0.25filt_VAF_norm_histogram_per_Tumour_type.pdf
    ├── cohort_maf_lofwmis
    │   ├── cohort_maf_lofwmis_VAF_norm_histogram_per_mut_type.pdf
    │   └── cohort_maf_lofwmis_VAF_norm_histogram_per_Tumour_type.pdf
    ├── cohort_maf_lofwmis_vaf_0.25filt
    │   ├── cohort_maf_lofwmis_vaf_0.25filt_VAF_norm_histogram_per_mut_type.pdf
    │   └── cohort_maf_lofwmis_vaf_0.25filt_VAF_norm_histogram_per_Tumour_type.pdf
    ├── cohort_maf_nhs_lof
    │   ├── cohort_maf_nhs_lof_VAF_norm_histogram_per_Hugo_Symbol.pdf
    │   ├── cohort_maf_nhs_lof_VAF_norm_histogram_per_Hugo_Symbol_pvartype.pdf
    │   ├── cohort_maf_nhs_lof_VAF_norm_histogram_per_mut_type.pdf
    │   ├── cohort_maf_nhs_lof_VAF_norm_histogram_per_Tumour_type.pdf
    │   └── cohort_maf_nhs_lof_VAF_norm_ridgeline_per_Hugo_Symbol.pdf
    ├── cohort_maf_nhs_lof_vaf_0.25filt
    │   ├── cohort_maf_nhs_lof_vaf_0.25filt_VAF_norm_histogram_per_Hugo_Symbol.pdf
    │   ├── cohort_maf_nhs_lof_vaf_0.25filt_VAF_norm_histogram_per_Hugo_Symbol_pvartype.pdf
    │   ├── cohort_maf_nhs_lof_vaf_0.25filt_VAF_norm_histogram_per_mut_type.pdf
    │   ├── cohort_maf_nhs_lof_vaf_0.25filt_VAF_norm_histogram_per_Tumour_type.pdf
    │   └── cohort_maf_nhs_lof_vaf_0.25filt_VAF_norm_ridgeline_per_Hugo_Symbol.pdf
    ├── cohort_maf_nhs_lof_vaf_0.25filt_lndepthf
    │   ├── cohort_maf_nhs_lof_vaf_0.25filt_lndepthf_VAF_norm_histogram_per_Hugo_Symbol.pdf
    │   ├── cohort_maf_nhs_lof_vaf_0.25filt_lndepthf_VAF_norm_histogram_per_Hugo_Symbol_pvartype.pdf
    │   ├── cohort_maf_nhs_lof_vaf_0.25filt_lndepthf_VAF_norm_histogram_per_mut_type.pdf
    │   ├── cohort_maf_nhs_lof_vaf_0.25filt_lndepthf_VAF_norm_histogram_per_Tumour_type.pdf
    │   └── cohort_maf_nhs_lof_vaf_0.25filt_lndepthf_VAF_norm_ridgeline_per_Hugo_Symbol.pdf
    ├── cohort_maf_nhs_lofwmis
    │   ├── cohort_maf_nhs_lofwmis_VAF_norm_histogram_per_Hugo_Symbol.pdf
    │   ├── cohort_maf_nhs_lofwmis_VAF_norm_histogram_per_Hugo_Symbol_pvartype.pdf
    │   ├── cohort_maf_nhs_lofwmis_VAF_norm_histogram_per_mut_type.pdf
    │   ├── cohort_maf_nhs_lofwmis_VAF_norm_histogram_per_Tumour_type.pdf
    │   └── cohort_maf_nhs_lofwmis_VAF_norm_ridgeline_per_Hugo_Symbol.pdf
    └── cohort_maf_nhs_lofwmis_vaf_0.25filt
        ├── cohort_maf_nhs_lofwmis_vaf_0.25filt_VAF_norm_histogram_per_Hugo_Symbol.pdf
        ├── cohort_maf_nhs_lofwmis_vaf_0.25filt_VAF_norm_histogram_per_Hugo_Symbol_pvartype.pdf
        ├── cohort_maf_nhs_lofwmis_vaf_0.25filt_VAF_norm_histogram_per_mut_type.pdf
        ├── cohort_maf_nhs_lofwmis_vaf_0.25filt_VAF_norm_histogram_per_Tumour_type.pdf
        └── cohort_maf_nhs_lofwmis_vaf_0.25filt_VAF_norm_ridgeline_per_Hugo_Symbol.pdf

10 directories, 52 files

```

### Consequence considered as LOF

List of consequences considered as LOF:

```bash
transcript_ablation	1
splice_acceptor_variant	1
splice_donor_variant	1
stop_gained	1
frameshift_variant	1
stop_lost	1
start_lost	1
transcript_amplification	1
feature_elongation	1
feature_truncation	1

```

List of consequences considered as LOF + missense:

```bash
transcript_ablation	1
splice_acceptor_variant	1
splice_donor_variant	1
stop_gained	1
frameshift_variant	1
stop_lost	1
start_lost	1
transcript_amplification	1
feature_elongation	1
feature_truncation	1
missense_variant	2
```


