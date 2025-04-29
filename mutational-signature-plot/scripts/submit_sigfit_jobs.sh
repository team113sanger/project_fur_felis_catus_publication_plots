#!/bin/bash

usage () {

	echo -e "\nUsage: $0 [options]\n
	-p \$PROJECTDIR 
	-g genome 
	-m MAF file 
	-c Cosmic signaure matrix
	-s Signal signature matrix
	where\n
	genome is one of: grch38, canfam3, mm10, felcat9, bostau9
	\$PROJECTDIR is the path to your directory with scripts and resources directories
	Optional opportunities file is required for non-human analysis or custom human target panels\n"

}

while getopts ":p:g:m:o:c:s:h" flag; do
        case "${flag}" in
                p) PROJECTDIR=$OPTARG
                   ;;
                g) genome=$OPTARG
                   ;;
                m) maf=$OPTARG
                   ;;
                c) cosmic=$OPTARG
                   ;;
                s) signal=$OPTARG
                   ;;
                : )
                   echo "Invalid option: $OPTARG requires an argument" 1>&2
                   ;;
                h | *)
                   usage
                   exit 0
                   ;;
    esac
done
shift $((OPTIND -1))


if [[ -z $PROJECTDIR || -z $genome || -z $maf ]]; then
	usage
	exit 1
fi

sigfit_script=scripts/run_sigfit_extraction_and_fitting.R
table_script=scripts/make_SNV_Sigfit_table.strand192.R

if [[ ! -f $sigfit_script || ! -f $cosmic || ! -f $signal || ! -f $table_script ]]; then
	usage
	echo "ERROR: one or more files does not exist:"
	echo "$sigfit_script"
	echo "$cosmic"
	echo "$signal"
	exit 1
fi

maf=$PROJECTDIR/$maf
cohort_size=$(awk -F'\t' 'NR>1 {print $11}' "${maf}" | sort | uniq | wc -l)

echo "Cohort size is $cohort_size"

Rscript scripts/get_mutational_opportunities.R $cohort_size felcat9 resources/baitset/S3250994_Feline_HSA_Jan2020_146.bed

mv felcat9_opportunities.RData $PROJECTDIR

echo "MAF is $maf"

Rscript $table_script $maf $genome

mv sigfit_var_table.strand192.tsv $PROJECTDIR

mkdir -p $PROJECTDIR/logs

cmd="Rscript $sigfit_script --file $PROJECTDIR/sigfit_var_table.strand192.tsv --output_dir $PROJECTDIR --cosmic $cosmic --signal $signal --opp_file $PROJECTDIR/felcat9_opportunities.RData"

bsub -e $PROJECTDIR/logs/extract_%J.e -o $PROJECTDIR/logs/extract_%J.o -M40000 -R 'select[mem>40000] rusage[mem=40000]' -q normal -n 6 "$cmd"