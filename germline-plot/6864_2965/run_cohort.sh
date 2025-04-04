#!/bin/bash
#BSUB -q oversubscribed
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000
#BSUB -oo nf_out.o
#BSUB -eo nf_out.e

PARAMS_FILE="/lustre/scratch127/casm/team113da/projects/fur_germline/6864_2965/cohort_params.json"

module load nextflow-23.10.0
module load /software/modules/ISG/singularity/3.11.4

# nextflow pull "https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/dermatlas_germlinepost_nf"

nextflow run "https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/dermatlas_germlinepost_nf" \
-r 0.2.1 \
-resume \
-params-file $PARAMS_FILE \
-profile farm22 \
-c /lustre/scratch127/casm/team113da/projects/fur_germline/scripts/nextflow.config
