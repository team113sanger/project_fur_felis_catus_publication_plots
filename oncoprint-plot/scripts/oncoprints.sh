#!/bin/bash
set -euo pipefail

# Color definitions
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print info messages in green
print_info() {
  echo -e "${GREEN}INFO: $1${NC}" >&2
}

# Function to print debug messages in blue
print_debug() {
  echo -e "${BLUE}DEBUG: $1${NC}" >&2
}

# List of cohorts
cohorts=(
  "6555_Lung_carcinoma"
  "6712_Oral_SCC"
  "6841_Meningioma"
  "6945_Cholangiocarcinoma"
  "6982_Lymphoma"
  "7040_BCC"
  "7098_Glioma"
  "6711_Cutaneous_SCC"
  "6713_Cutaneous_MCT"
  "6864_Pancreatic_carcinoma"
  "6973_OSA"
  "6990_Mammary_carcinoma"
  "7097_CRC"
)

print_info "Starting script execution"

# This script location
script_dir=$(dirname "$(realpath "$0")")
print_debug "Script directory: $script_dir"

# MAF directory path
maf_dir="${script_dir}/../data/maf_files"
if [[ ! -d "$maf_dir" ]]; then
  print_info "MAF directory not found: $maf_dir"
  exit 1
fi
print_debug "MAF directory: $maf_dir"

# Mapping file path
mapping_file="${script_dir}/../data/Felis_catus_9.0.gene_chromosome_mapping.txt"
if [[ ! -f "$mapping_file" ]]; then
  print_info "Mapping file not found: $mapping_file"
  exit 1
fi
print_debug "Mapping file: $mapping_file"

# Study summary directory path
study_summary_dir="${script_dir}/../data/genemetrics_study_summary"
if [[ ! -d "$study_summary_dir" ]]; then
  print_info "Study summary directory not found: $study_summary_dir"
  exit 1
fi
print_debug "Study summary directory: $study_summary_dir"

# Exclude file
exclude_file="${script_dir}/../data/FUR_cat.samples_to_exclude.nucleotide_variants.txt"
if [[ ! -f "$exclude_file" ]]; then
  print_info "Exclude file not found: $exclude_file"
  exit 1
fi
print_debug "Exclude file: $exclude_file"

# Breakdown csv file by cohort directory
breakdown_dir="${script_dir}/../data/breakdown"
if [[ ! -d "$breakdown_dir" ]]; then
  print_info "Breakdown directory not found: $breakdown_dir"
  exit 1
fi
print_debug "Breakdown directory: $breakdown_dir"

total_cohorts=${#cohorts[@]}
print_info "Processing $total_cohorts cohorts"

# Loop over each cohort
for i in "${!cohorts[@]}"; do
  cohort=${cohorts[$i]}
  current=$((i+1))
  
  print_info "Processing cohort $current/$total_cohorts: $cohort"
  
  # Extract the numeric ID (the part before the first underscore)
  id=$(echo "$cohort" | cut -d'_' -f1)
  # Extract the cancer type (everything after the first underscore)
  cancer=$(echo "$cohort" | cut -d'_' -f2-)
  
  print_debug "Cohort ID: $id, Cancer type: $cancer"
  
  # Search for the keepPA MAF file that contains the cohort ID
  maf_file=$(find "$maf_dir" -type f -name "*keepPA*${id}*.maf" | head -n 1)
  if [[ -z "$maf_file" ]]; then
    print_info "MAF file for cohort ID ${id} not found. Skipping ${cohort}..."
    continue
  fi
  print_debug "MAF file: $maf_file"

  # Search for the breakdown file that contains the cohort ID
  breakdown_file=$(find "$breakdown_dir" -type f -name "*${id}*.csv" | head -n 1)
  if [[ -z "$breakdown_file" ]]; then
    print_info "Breakdown file for cohort ID ${id} not found. Skipping ${cohort}..."
    continue
  fi
  print_debug "Breakdown file: $breakdown_file"

  # Find study summary file path by searching for the cohort ID
  study_summary=$(find "$study_summary_dir" -type f -name "*${id}*.csv" | head -n 1)
  if [[ -z "$study_summary" ]]; then
    print_info "Study summary file for cohort ID ${id} not found. Skipping ${cohort}..."
    continue
  fi
  print_debug "Study summary file: $study_summary"
  
  # Define the output file name
  output_file="foo/${cohort}.oncoprint.pdf"
  print_debug "Output file: $output_file"
  
  # Build the command as a multi-line string with proper line continuation
  cmd="python3 -m fur_cnvkit oncoprint \\
  ${maf_file} \\
  ${study_summary} \\
  ${cancer} \\
  --output_file ${output_file} \\
  --num_recurrent_samples 10 \\
  --gene_chrom_file ${mapping_file} \\
  --gene_sort position \\
  --cluster_samples \\
  --breakdown_file ${breakdown_file} \\
  --exclude_file ${exclude_file}"

  print_info "Running command for ${cohort}:"
  echo -e "${BLUE}$cmd${NC}" # Print the command in blue
  
  # Execute the command
  eval "$cmd"
  
  print_info "Completed processing cohort $current/$total_cohorts: $cohort"
done

print_info "Script execution completed"