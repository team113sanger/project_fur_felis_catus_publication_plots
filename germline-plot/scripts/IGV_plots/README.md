# IGV screen shots of the alignments that support the germline variants identified within cancer predisposition genes

# Overview
Documentation on the process followed to  generates IGV screen shots of the mapped reads on the genomic locations of the [`FUR_keepPA_vaf_size_filt_germline.nhs.lof.maf`](../../cohort_analysis/cohort_tables/FUR_keepPA_vaf_size_filt_germline.nhs.lof.maf) germline variants 

# Requirements

## Required software
The following software is required to be installed and present in the running path to run the scripts:
- IGV [v2.16.2](https://igv.org/doc/desktop/#DownloadPage/)
- R v4.4.0
- samtools [v1.19.2](https://github.com/samtools/samtools/releases/tag/1.19.2)
- bgzip v1.19.1 form [htslib v1.19.1 ](https://github.com/samtools/htslib/releases/tag/1.19.1)
- tabix v1.19.1 from [samtools v1.19.1](https://github.com/samtools/bcftools/releases/tag/1.19.1)

## Required files

To be able to run the scripts, you need to have the Felis catus FelCat9.0 genome and **ENSEMBLv104** annotation files. The files are available from ENSEMBL. The files are:

### _Felis catus_ genome and annotation files

To obtain the required files, please follow the instructions below, setting the `PROJECTDIR` variable to the directory where clone this repo into. The files will be saved in the following directory: `<PROJECTDIR>/resources/IGV/` 

```bash
PROJECTDIR=/lustre/scratch127/casm/team113da/projects/fur_germline

IGVRESOURCES=${PROJECTDIR}/resources/IGV
mkdir -p ${IGVRESOURCES}
cd ${IGVRESOURCES}  

# Donwlad the Felis catus genome from Ensembl:
wget ftp://ftp.ensembl.org/pub/release-104/fasta/felis_catus/dna/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.gz
gunzip -v -d Felis_catus.Felis_catus_9.0.dna.toplevel.fa.gz
# Index the genome file:
samtools faidx Felis_catus.Felis_catus_9.0.dna.toplevel.fa

# Download the Felis catus GTF file from Ensembl:
wget ftp://ftp.ensembl.org/pub/release-104/gtf/felis_catus/Felis_catus.Felis_catus_9.0.104.gtf.gz
gunzip -v -d Felis_catus.Felis_catus_9.0.104.gtf.gz
sort -k1,1 -k4,4n Felis_catus.Felis_catus_9.0.104.gtf > Felis_catus.Felis_catus_9.0.104.sorted.gtf
bgzip Felis_catus.Felis_catus_9.0.104.sorted.gtf
# Index the GTF file:
tabix -p gff Felis_catus.Felis_catus_9.0.104.sorted.gtf.gz

```

# Generate the batch files to plot the IGV screen shots

The batch files are generated using the `igv_plots.R` script. The script takes as input the following files:

- `--projdir` : The project directory where the repository got cloned into. The output files will be saved in the following directory: `<PROJECTDIR>/cohort_analysis/cohort_tables/igv_plots/batch_files/`
- `--bam_dir` :The directory where the BAM files are located with the following format: `<BAM_DIR>/Tumor_Sample_Barcodes/<PROJECTID>/<SAMPLEID>/<SAMPLEID>.sample.dupmarked.bam`
  - The BAM files should be indexed and the index files should be in the same directory as the BAM files
- `--study_list` : TSV file without header containing two columns "Tumour_type" and "StudyID_ProjectID"
- `--input_maf` : The MAF file with the following format:
  - The MAF file should be in the format of the MAF file. The MAF file should contain the following columns:
    - Chromosome
    - Start_Position
    - End_Position
    - Reference_Allele
    - Tumor_Seq_Allele2
    - Tumor_Sample_Barcode
- `--genome_file` : The genome file in FASTA format. The genome file should be indexed and the index file should be in the same directory as the genome file
- `--gtf_file` : The GTF file in GFF format. The GTF file should be indexed and the index file should be in the same directory as the GTF file


```bash 
PROJECTDIR=/lustre/scratch127/casm/team113da/projects/fur_germline

cd ${PROJECTDIR}
# Load the required modules:
source ${PROJECTDIR}/scripts/IGV_plots/source_me.sh

Rscript ${PROJECTDIR}/scripts/IGV_plots/igv_plots.R --projdir ${PROJECTDIR:?unset} \
  --bam_dir /lustre/scratch125/casm/team113da/users/bf14/samples \
  --study_list ${PROJECTDIR:?unset}/metadata/studies_reformat.tsv \
  --input_maf ${PROJECTDIR:?unset}/cohort_analysis/cohort_tables/FUR_keepPA_vaf_size_filt_germline.nhs.lof.vafn0.25filt.lndepthf.maf \
  --genome_file ${PROJECTDIR:?unset}/resources/IGV/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
  --gtf_file ${PROJECTDIR:?unset}/resources/IGV/Felis_catus.Felis_catus_9.0.104.sorted.gtf.gz


```
Output batch files will have the following format : <Hugo_Symbol>_<Chromosome>_<Start_Position>_<End_Position>_<Reference_Allele>_<Tumor_Seq_Allele2>.txt


# Plot the IGV screen shots

Once the batch files are generated, you can run IGV with the batch files. The batch files are located in the following directory: `<PROJECTDIR>/cohort_analysis/cohort_tables/igv_plots/batch_files/`

```bash
PROJECTDIR=/lustre/scratch127/casm/team113da/projects/fur_germline
# Load the required modules:
source ${PROJECTDIR}/scripts/IGV_plots/source_me.sh

cd ${PROJECTDIR}

for i in `ls -1 ${PROJECTDIR:?unset}/cohort_analysis/cohort_tables/igv_plots/batch_files/*.txt`; do
  echo "Processing ${i}" ;
  # Run IGV with the batch file:
  igv.sh -b ${i:?unset} ;
done

```

Output png files will have the following format : <Hugo_Symbol>_<Chromosome>_<Start_Position>_<End_Position>_<Reference_Allele>_<Tumor_Seq_Allele2>.png
