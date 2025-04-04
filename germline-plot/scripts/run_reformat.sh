#!/bin/bash
#BSUB -q normal
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000
#BSUB -oo reformat_out.o
#BSUB -eo reformat_out.e

/lustre/scratch127/casm/team113da/projects/fur_germline/scripts/QC/germline_gatk_qc.sh -l /lustre/scratch127/casm/team113da/projects/fur_germline/vcfs_6945_3142.list -m germline_6945_3142.maf -s /lustre/scratch127/casm/team113da/projects/fur_germline/scripts -b felCat9 -r -a 99_Lives_AF -f filter2 -g &>reformat.log
