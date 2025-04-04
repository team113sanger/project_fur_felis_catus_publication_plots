#!/bin/bash

PROJECTDIR=$PWD
SCRIPTDIR=${PROJECTDIR}/scripts
LISTDIR=/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_cat_cohort_files


for STUDY in `cut -f 2 ${PROJECTDIR}/metadata/studies_reformat.tsv`; do
    MAF_IN=${PROJECTDIR}/${STUDY}/mafs/keepPA_vaf_size_filt_germline_${STUDY}.maf
    MAF_OUT=${PROJECTDIR}/${STUDY}/risk_analysis/keepPA_vaf_size_filt_germline_${STUDY}.risk.maf
    MAF_OUTWMISS=${PROJECTDIR}/${STUDY}/risk_analysis/keepPA_vaf_size_filt_germline_${STUDY}.risk.wmiss.maf
    SAMPLELIST=${LISTDIR}/${STUDY}/sample_lists/${STUDY}.samples_to_keep.nucleotide_variants.txt

    # Calculate the height of the plot based on sample list
    plot_height=`wc -l ${SAMPLELIST} | awk '{print $1}'`
    check=$(echo $plot_height / 5 | bc -l | perl -ne 's/\S+\.(\S+)/$1/;print')
    if [[ "$plot_height" -le 5 ]]; then
        plot_height=1
    elif [[ "$check" -gt 0 ]]; then
        let plot_height=1+$(echo $plot_height / 5 | bc)
    else
        let plot_height=$plot_height/5
    fi
    plot_height=$(echo $plot_height*1.5 | bc -l)

    # Create output directory
    echo "### ${MAF_IN} ###"
    mkdir -p  ${PROJECTDIR}/${STUDY}/risk_analysis/

    # Create MAFs with putative risk variants
    head -n 1 ${MAF_IN} > ${MAF_OUT}
    for f in `grep -vw 2 ${PROJECTDIR}/resources/germline/vep_csq_priority.txt | cut -f 1`; do
        grep -w $f ${MAF_IN} >> ${MAF_OUT}
    done

    # Create MAFs with putative risk variants and missense variants
    head -n 1 ${MAF_IN} > ${MAF_OUTWMISS}
    for f in `grep -vw 2 ${PROJECTDIR}/resources/germline/vep_csq_priority.txt | cut -f 1`; do
        grep -w $f ${MAF_IN} >> ${MAF_OUTWMISS}
    done
    grep -w missense_variant ${MAF_IN} >> ${MAF_OUTWMISS}

    # Create QC plots
    mkdir -p ${PROJECTDIR}/${STUDY}/risk_analysis/QC
    cd ${PROJECTDIR}/${STUDY}/risk_analysis/QC
    echo $PWD
    echo "Creating QC plots"
    Rscript ${SCRIPTDIR}/MAF/plot_vaf_vs_depth_from_maf.R --file $MAF_OUT --samplefile ${SAMPLELIST} --width 10 --height $plot_height -ncol 5 --germline
    #bash ${PROJECTDIR}/scripts/QC/somatic_variants_qc.sh -m $MAF_IN -b GRCh38 -p -g -n -s ${SCRIPTDIR} -f filter2
    
    # Create oncoplots
    echo $PWD
    echo "Creating tile plot"
    mkdir -p ${PROJECTDIR}/${STUDY}/risk_analysis/TILE_PLOT
    cd ${PROJECTDIR}/${STUDY}/risk_analysis/TILE_PLOT
    Rscript ${SCRIPTDIR}/MAF/maketileplot_from_maf.R -a ${MAF_OUT} -s ${SAMPLELIST} -c 3  --sortbyfrequency -w 8 -t 5

    cd ${PROJECTIDR}
done



