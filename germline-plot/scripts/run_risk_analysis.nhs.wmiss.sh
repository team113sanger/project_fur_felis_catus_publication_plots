#!/bin/bash
# Wrapper script to take the maf files filterd them and generate 

PROJECTDIR=$PWD
#SCRIPTDIR=/lustre/scratch127/casm/team113da/projects/fur_germline/scripts
SCRIPTDIR=${PROJECTDIR}/scripts
LISTDIR=/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_cat_cohort_files
NHSLIST=${PROJECTDIR}/resources/ensembl_v104/nhs_genes_ensv104_cat_human_ortholog.list

for STUDY in `cut -f 2 ${PROJECTDIR}/metadata/studies_reformat.tsv`; do
    MAF_IN=${PROJECTDIR}/${STUDY}/mafs/keepPA_vaf_size_filt_germline_${STUDY}.maf
    MAF_OUT1=${PROJECTDIR}/${STUDY}/risk_analysis/keepPA_vaf_size_filt_germline_${STUDY}.nhs.lof.maf
    MAF_OUT2=${PROJECTDIR}/${STUDY}/risk_analysis/keepPA_vaf_size_filt_germline_${STUDY}.nhs.all.maf
    MAF_OUT3=${PROJECTDIR}/${STUDY}/risk_analysis/keepPA_vaf_size_filt_germline_${STUDY}.nhs.lof.wmiss.maf
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

    # Create MAFs with NHS gene orthologs (including missense mutations) and 'voi=yes' variants
    head -n 1 ${MAF_IN} > ${MAF_OUT2}
    for f in `cut -f 2 -d "/" ${NHSLIST}`; do
        grep -w $f ${MAF_IN} | grep -w yes >> ${MAF_OUT2}
    done

    # Create MAFs with putative risk variants in NHS gene orthologs
    head -n 1 ${MAF_IN} > ${MAF_OUT1}
    for f in `grep -vw 2 ${PROJECTDIR}/resources/germline/vep_csq_priority.txt | cut -f 1`; do
        grep -w $f ${MAF_OUT2} >> ${MAF_OUT1}
    done

    # Create MAFs with putative risk variants in NHS gene orthologs
    head -n 1 ${MAF_IN} > ${MAF_OUT3}
    for f in `grep -vw 2 ${PROJECTDIR}/resources/germline/vep_csq_priority.txt | cut -f 1`; do
        grep -w $f ${MAF_OUT2} >> ${MAF_OUT3}
    done
    grep -w missense_variant ${MAF_OUT2} >> ${MAF_OUT3}


    # Create QC plots
    mkdir -p ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/ALL/QC
    cd ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/ALL/QC
    echo $PWD
    echo "Creating QC plots"
    Rscript ${SCRIPTDIR}/MAF/plot_vaf_vs_depth_from_maf.R --file ${MAF_OUT2} --samplefile ${SAMPLELIST} --width 10 --height $plot_height -ncol 5 --germline
    #bash ${PROJECTDIR}/scripts/QC/somatic_variants_qc.sh -m $MAF_IN -b GRCh38 -p -g -n -s ${SCRIPTDIR} -f filter2
    
    mkdir -p ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/LOF/QC
    cd ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/LOF/QC
    echo $PWD
    echo "Creating QC plots"
    Rscript ${SCRIPTDIR}/MAF/plot_vaf_vs_depth_from_maf.R --file ${MAF_OUT1} --samplefile ${SAMPLELIST} --width 10 --height $plot_height -ncol 5 --germline
    #bash ${PROJECTDIR}/scripts/QC/somatic_variants_qc.sh -m $MAF_IN -b GRCh38 -p -g -n -s ${SCRIPTDIR} -f filter2

    mkdir -p ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/LOF_wmiss/QC
    cd ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/LOF_wmiss/QC
    echo $PWD
    echo "Creating QC plots"
    Rscript ${SCRIPTDIR}/MAF/plot_vaf_vs_depth_from_maf.R --file ${MAF_OUT3} --samplefile ${SAMPLELIST} --width 10 --height $plot_height -ncol 5 --germline
    
    # Create oncoplots NHS ALL
    mkdir -p ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/ALL/TILE_PLOT
    cd ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/ALL/TILE_PLOT
    echo $PWD
    echo "Creating tile plot"
    Rscript ${SCRIPTDIR}/MAF/maketileplot_from_maf.R -a ${MAF_OUT2} -s ${SAMPLELIST} -c 3 --sortbyfrequency -w 8 -t 5

    # Create oncoplots NHS LOF
    mkdir -p ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/LOF/TILE_PLOT
    cd ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/LOF/TILE_PLOT
    echo $PWD
    echo "Creating tile plot"
    Rscript ${SCRIPTDIR}/MAF/maketileplot_from_maf.R -a ${MAF_OUT1} -s ${SAMPLELIST} -c 3 --sortbyfrequency -w 8 -t 5

    # Create oncoplots NHS LOF with missense
    mkdir -p ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/LOF_wmiss/TILE_PLOT
    cd ${PROJECTDIR}/${STUDY}/risk_analysis/NHS/LOF_wmiss/TILE_PLOT
    echo $PWD
    echo "Creating tile plot"
    Rscript ${SCRIPTDIR}/MAF/maketileplot_from_maf.R -a ${MAF_OUT3} -s ${SAMPLELIST} -c 3 --sortbyfrequency -w 8 -t 5

    cd ${PROJECTIDR}
done



