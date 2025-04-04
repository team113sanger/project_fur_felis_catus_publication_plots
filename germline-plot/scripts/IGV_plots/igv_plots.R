#!/usr/bin/env Rscript
 
library(here)
library(dplyr)
library(readr)
library(optparse)

#Get the options for tthe script
# Set the options list of the program
progname<- "igv_plots.R"
# Set the options list of the program
option_list<- list(
  make_option(c("--projdir"), action="store_true", default=NA, type="character", help="This full path to the project directory "),
  make_option(c("--bam_dir" ), action="store_true",type="character", default=NA, help="Full path to the BAM file directory. The files should be stored in the format <BAM_DIR>/Tumor_Sample_Barcodes/<PROJECTID>/<SAMPLEID>/<SAMPLEID>.sample.dupmarked.bam "),
  make_option(c("--study_list" ), action="store_true",type="character", default=NA, help="TSV file with two columns: Tumour_type   Study_project"),
  make_option(c("--input_maf" ), action="store_true",type="character", default=NA, help="Full path to the input MAF file"),
  make_option(c("--genome_file" ), action="store_true",type="character", default=NA, help="Full path to the reference genome uncompressed FASTA file "),
  make_option(c("--gtf_file" ), action="store_true",type="character", default=NA, help="Full path to the sorted, compressed and tabix indexed GTF file with the annotation "), 
  make_option(c("--use_here" ), action="store_true",type="character", default=FALSE, help="Optional, If set it assumes the project directory is taken by R here() function (running directory) ")    
)
#### Indicate that the file follows to the options used
parser<- OptionParser(usage = paste(progname, " [options] --projdir <PROJECTDIR> --bam_dir <BAMDIR> --study_list  <PROJECTDIR>/metadata/study_list.tsv --input_maf <input_maf> --genome_file /PATH/TO/genome.fa --gtf_file /PATH/TO/genome_annotation.gtf.gz ",sep=""), option_list=option_list)
arguments<- parse_args(parser, positional_arguments = 0) # no files after the options shoudl be provided
opt<-arguments$options

if(opt$use_here){
  # Configuration
  projdir<- here()
} else {
  # Configuration
  projdir<-opt$projdir
}
  output_dir <- file.path(projdir, "cohort_analysis", "cohort_tables", "igv_plots")
  screenshot_dir <- file.path(projdir, "cohort_analysis", "cohort_tables", "igv_plots", "screenshots")
  batch_dir <- file.path(projdir, "cohort_analysis", "cohort_tables", "igv_plots","batch_files")
  dir.create(screenshot_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(batch_dir, showWarnings = FALSE, recursive = TRUE)
########## INPUT files
#location of the BAM files 
#Format: <BAM_DIR>/Tumor_Sample_Barcodes/<PROJECTID>/<SAMPLEID>/<SAMPLEID>.sample.dupmarked.bam 
# Eg. <BAMDIR>/Tumor_Sample_Barcodes/2711/CATD0078b/CATD0078b.sample.dupmarked.bam 
bam_dir<- opt$bam_dir
input_file<- opt$input_maf
study_list<- opt$study_list


# Set the path to the GTF file
#gtf_file <- file.path(projdir, "resources", "IGV", "Felis_catus.Felis_catus_9.0.104.sorted.gtf.gz")
#genome_name <- file.path(projdir, "resources", "IGV","Felis_catus.Felis_catus_9.0.dna.toplevel.fa")
gtf_file<- opt$gtf_file
genome_name<- opt$genome_file

# Check the input files and resource files exist
if(!file.exists(input_file)){
  stop(paste0("Input MAF file doesn't exist"))
}
if(!file.exists(genome_name)){
  stop(paste0("Genome file doesn't exist"))
}
if(!file.exists(gtf_file)){
  stop(paste0("GTF file doesn't exist"))
}
if(!dir.exists(bam_dir)){
  stop(paste0("BAM directory doesn't exist"))
}


# Load input
maf <- read_tsv(input_file, col_types = cols())
studyl <- read_tsv(study_list, col_types = cols(), col_names = c("Tumour_type", "Study_project"))
studyl$ProjectID<- strsplit(studyl$Study_project, "_") %>% sapply(., function(x) x[2])
# Get the PRoject ID from the study list to the maf file
maf$ProjectID<- studyl$ProjectID[match(maf$Tumour_type , studyl$Tumour_type)]

# Group by region
grouped <- maf %>%
  group_by(Chromosome, Start_Position, End_Position, Reference_Allele,Tumor_Seq_Allele2,Hugo_Symbol) %>%
  summarise(Tumor_Sample_Barcodes = list(Tumor_Sample_Barcode), ProjectIDs=list(ProjectID), .groups = "drop")



# Generate batch commands
for (i in seq_len(nrow(grouped))) {
  
  Chromosome <- grouped$Chromosome[i]
  Start_Position <- grouped$Start_Position[i]
  End_Position <- grouped$End_Position[i]
  Reference_Allele <- grouped$Reference_Allele[i]
  Tumor_Seq_Allele2 <- grouped$Tumor_Seq_Allele2[i]
  Hugo_Symbol <- grouped$Hugo_Symbol[i]
  sample_list <- grouped$Tumor_Sample_Barcodes[[i]]
  project_list <- grouped$ProjectIDs[[i]]
  # Create a batch file for each region
  batch_file <- file.path(batch_dir,paste0(paste(Hugo_Symbol, Chromosome,  Start_Position, End_Position, Reference_Allele,Tumor_Seq_Allele2, sep="_"), "_igv_batch.txt"))
  # Create a snapshot name  
  snapshot_name <- paste0(paste(Hugo_Symbol, Chromosome,  Start_Position,
   End_Position, Reference_Allele,Tumor_Seq_Allele2, sep="_"), ".png")
  if(length(sample_list)!= length(project_list)){
    stop(paste0("Sample list and project list are not the same length"))
  }
  
  # Start fresh batch script
  cat("", file = batch_file)
  # Add genome and gtf file to batch script
  cat(
    "new\n",
    paste0("genome ", genome_name, "\n"),
    paste0("load ", gtf_file, "\n"),
    file = batch_file,
    append = TRUE
  )
  for (j in seq_len(length(sample_list))) {
    sample<-NULL
    project<-NULL
    sample <- sample_list[j]
    project <- project_list[j] 
    #Format: <BAM_DIR>/<PROJECTID>/<SAMPLEID>/<SAMPLEID>.sample.dupmarked.bam 
    # Eg. <BAMDIR>/2711/CATD0078b/CATD0078b.sample.dupmarked.bam 
    bam_path <- file.path(bam_dir,project,sample, paste0(sample, ".sample.dupmarked.bam"))
    # Check if the sample is in the bam_dir
    if(!file.exists(bam_path)){
      stop(paste0("BAM file doesn't exist for sample ", sample))
    }
    cat(
      paste0("load ", bam_path, "\n"),
      file = batch_file,
      append = TRUE
    )
  }
  # Add 23 bp UP and downstream of location the region to the batch script
  cat(
    paste0("snapshotDirectory ", screenshot_dir, "\n"),
    paste0("goto ", Chromosome, ":", (Start_Position-23), "-", (End_Position+23), "\n"),
    "sort base\n",
    "collapse\n",
    "colorBy FIRST_OF_PAIR_STRAND\n",
    paste0("snapshot ", snapshot_name, "\nexit\n"),
    file = batch_file,
    append = TRUE
  )
}
