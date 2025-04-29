#!/bin/bash
set -e

source source_me.sh

MANIFEST="study_manifest.tsv"
OUTPUT_DIR="/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_mutational_signatures/analysis"
MAF_DIR="/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_hotspot_mutations/analysis/finalised_mafs"
SIGFIT_SCRIPT="scripts/submit_sigfit_jobs.sh"
SIGNATURES_DIR="/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_mutational_signatures/resources/signatures"
COSMIC_SIG="${SIGNATURES_DIR}/COSMIC_v3.4_SBS_GRCh38.sigfit.txt"
REFSIG="${SIGNATURES_DIR}/RefSig_SBS_v2.03.reformatted.tsv"
GENOME="felcat9"

mkdir -p "$OUTPUT_DIR"

# Read manifest and run sigfit for each study
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r INDICATION SEQSCAPE_ID CANAPPS_ID; do
    STUDY_DIR="${SEQSCAPE_ID}_${CANAPPS_ID}"
    FULL_PATH="${OUTPUT_DIR}/${STUDY_DIR}"

    # Create study dir
    mkdir -p "$FULL_PATH"

    # Copy over maf
    MAF_FILENAME="keep_vaf_size_filt_matched_${SEQSCAPE_ID}_${CANAPPS_ID}.finalised.maf"
    cp "${MAF_DIR}/$MAF_FILENAME" "$FULL_PATH/"

    # Run submit_sigfit_jobs.sh 
    bash "$SIGFIT_SCRIPT" -p "$FULL_PATH" -g "$GENOME" -m "$MAF_FILENAME" -c "$COSMIC_SIG" -s "$REFSIG"

    echo "Submitted study: ${SEQSCAPE_ID}_${CANAPPS_ID}"
done

echo "done"