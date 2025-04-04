# fur_germline
A set of analysis scripts, metadata and results for running germline variant calling on FUR cats.

## Generation of Germline Variant calls

-  This project was setup with `scripts/populate_project.R`. Sample lists for each cohort were prepared from:
`/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_cat_cohort_files`
- Each cohort was processed using `https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/dermatlas_germlinepost_nf` version 0.2.1. Cohort analyses were parameterised using a
`<COHORT_DIR>/cohort_params.json` file. Each pipeline run was launched with a `run_cohort.sh` wrapper like so

```
cd 7098_3140
bsub < run_cohort.sh
```

## Post-processing of germline calls

Create MAFs from VCFs

```
PROJECTDIR=/lustre/scratch127/casm/team113da/projects/fur_germline
source ${PROJECTDIR}/scripts/QC/source_me.sh

for ID in `cut -f 2 ${PROJECTDIR}/metadata/studies_reformat.tsv`; do
  samplelist=/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_cat_cohort_files/${ID}/sample_lists/${ID}.samples_to_keep.nucleotide_variants.txt
  mkdir -p ${ID}/mafs
  cd ${ID}/mafs
  echo $PWD
  dir -1 ${PROJECTDIR}/${ID}/results/Final_joint_call/*.target.pass.vep.vcf.gz > vcfs_${ID}.list
  VCFLIST=$PWD/vcfs_${ID}.list

cat <<END > run_reformat.sh
#!/bin/bash
#BSUB -q normal
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000
#BSUB -oo reformat_out.o
#BSUB -eo reformat_out.e

${PROJECTDIR}/scripts/QC/somatic_variants_qc.sh -l ${VCFLIST} -m germline_${ID}.maf -s ${PROJECTDIR}/scripts -b felCat9 -a 99_Lives_AF -f filter2 -g -n &>reformat.log
END
  bsub < run_reformat.sh
  cd ${PROJECTDIR}
done
```

## Search for germline risk variants

### Loss-of-function variants
```
cd /lustre/scratch127/casm/team113da/projects/fur_germline

# This wrapper pulls out LOF variants from MAFs, generates QC plots
# and gene tile plots of variants present in at least 3 samples

bash ./scripts/run_risk_analysis.sh
```

### Variants in NHS cancer risk genes (DERMATLAS list)

BioMart (Ens v104) was used to extract cat<->human orthologs.
Genes in the NHS gene panel are matched to orthologs in cat.
Pull out genes in MAF using ENS IDs (not all cat genes have Hugo Symbols)

```
bash ./scripts/run_risk_analysis.nhs.sh
```

### Create XLSX files to share with the scientist
```
source scripts/MAF/source_me.sh

for f in `cut -f 2 metadata/studies_reformat.tsv`; do
  cd $f/risk_analysis
  echo $PWD
  for f in *maf */*/*tsv; do
    Rscript ../../scripts/MAF/tsv2xlsx.R $f
  done
  cd ../../
done
```

## Compare somatic and germline MAFs to get 2-hit genes

```
SOMATIC_DIR=/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_cat_maf/analysis
PROJECTDIR=/lustre/scratch127/casm/team113da/projects/fur_germline
source ${PROJECTDIR}/scripts/QC/source_me.sh

for f in `cut -f 2 metadata/studies_reformat.tsv`; do
  cd $f/risk_analysis;
  echo $PWD;
  maf="${PROJECTDIR}/${f}/mafs/keepPA_vaf_size_filt_germline_${f}.maf";
  out="${PROJECTDIR}/${f}/mafs/keepPA_vaf_size_filt_germline_${f}.2hit.maf"
  Rscript ${PROJECTDIR}/scripts/compare_mafs.R ${maf} ${SOMATIC_DIR}/${f}/keepPA_vaf_size_filt_matched_${f}.maf GERMLINE SOMATIC ${out};
  cd ../../;
done
```

### Sort MAFs and create xlsx files to post on gitlab for Louise

```
for f in `cat studies.list`; do
  cd $f
  head -n1 mafs/keepPA_vaf_size_filt_germline_${f}.2hit.maf > mafs/keepPA_vaf_size_filt_germline_${f}.2hit.sorted.maf
  tail -n +2 mafs/keepPA_vaf_size_filt_germline_${f}.2hit.maf | sort >> mafs/keepPA_vaf_size_filt_germline_${f}.2hit.sorted.maf
  Rscript ../scripts/MAF/tsv2xlsx.R mafs/keepPA_vaf_size_filt_germline_${f}.2hit.sorted.maf
  cd ../
done
```

### Create a summary MAF

```
head -n1 ./7040_3064/mafs/keepPA_vaf_size_filt_germline_7040_3064.2hit.sorted.maf > cohort_analysis/2-hit_genes/feline_cohort_keepPA_vaf_size_filt_germline.maf
grep -hv Hugo_Symbol */mafs/*sorted.maf >> cohort_analysis/2-hit_genes/feline_cohort_keepPA_vaf_size_filt_germline.maf

cd /lustre/scratch127/casm/team113da/projects/fur_germline/cohort_analysis/2-hit_genes
Rscript ../../scripts/MAF/tsv2xlsx.R feline_cohort_keepPA_vaf_size_filt_germline.maf 
```


### Generate cohort summary tables with LOF and LOF+missense variants

The information about the files generated can be found here. [README](./scripts/paper_tables/README.md)

### Generate the IGV plots from LOF variants on orthologues from Human Cancer predisposition genes
The information about the files generated can be found here. [README](./scripts/IGV_plots/README.md)

