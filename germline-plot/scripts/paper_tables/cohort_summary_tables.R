#!/usr/bin/env Rscript
# This script is used to generate the cohort summary tables 

library(here)
library(data.table)
library(ggplot2)
library(ggridges)
library(paletteer)
library(RColorBrewer)

## Functions to check if a file exists
filecheck<-function(file){
        if(file.exists(file)){
                message(paste0("File exists: ", file))
        }else{
                stop(paste0("File does not exist: ", file))
        }
        return(file)
}
## Generate a VAf plot per gin
plot_vafdists<-function(maf, prefix, outdir, plot_per_gene=TRUE){
        #Check if outdir exists
        if(!dir.exists(outdir)){
                dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        }
        ntums<-length(unique(maf$Tumour_type))
        ngmutgenes<-length(unique(maf$Hugo_Symbol))
        message(paste0("Number of Tumour types: ", ntums))
        message(paste0("Number of genes: ", ngmutgenes))
        #PLOT VAF dist per Variant type 
        vaf_histogram_pmtype <- ggplot(maf, aes(x = VAF_norm, fill=Variant_Type)) +
                geom_histogram(binwidth = 0.01, color = "black", alpha = 0.5) +
                labs(title = paste0(prefix, "\n", "Histogram of VAF_norm Distribution per Tumour_type"),
                        x = "VAF_norm",
                        y = "Frequency") +
                scale_fill_paletteer_d("ggthemes::Color_Blind")+
                theme_bw() +
                theme(strip.text.y = element_text(size = 7, face= "bold"),
                        legend.position = "none")+
                facet_wrap(~Variant_Type, ncol = 1, scales = "free_y")

        # Save the histogram plot
        ggsave(file.path(outdir,paste0(prefix,  "_VAF_norm_histogram_per_mut_type.pdf")),
        vaf_histogram_pmtype, width = 10, height = 9)
        #PLOT VAF dist per tumour type
        vaf_histogram_pttype <- ggplot(maf, aes(x = VAF_norm, fill=Tumour_type)) +
                geom_histogram(binwidth = 0.01, fill = "#006BA4FF", color = "black", alpha = 0.5) +
                labs(title = paste0(prefix, "\n", "Histogram of VAF_norm Distribution per Tumour_type"),
                        x = "VAF_norm",
                        y = "Frequency") +
                theme_bw() +
                theme(strip.text.y = element_text(size = 7, face= "bold"),
                        legend.position = "none")+
                facet_wrap(~Tumour_type, ncol = 1, scales = "free_y")

        # Save the histogram plot
        ggsave(file.path(outdir,paste0(prefix,  "_VAF_norm_histogram_per_Tumour_type.pdf")),
        vaf_histogram_pttype, width = 10, height = (15+ (ntums*0.1)))
        
        if(plot_per_gene){
                # Function to plot the VAF distribution per gene
                vaf_histogram <- ggplot(maf, aes(x = VAF_norm)) +
                                geom_histogram(binwidth = 0.01, fill = "#006BA4FF", color = "black", alpha = 0.5) +
                                labs(title = paste0(prefix, "\n", "Histogram of VAF_norm Distribution per Hugo_Symbol"),
                                        x = "VAF_norm",
                                        y = "Frequency") +
                                theme_bw() +
                                theme(strip.text.y = element_text(size = 7, face= "bold"),
                                        legend.position = "none") +
                                facet_grid(Hugo_Symbol ~ ., scales = "free_y")

                # Save the histogram plot
                hheight<- NULL
                ifelse((15+ (ngmutgenes*0.1))>40, hheight<-40, hheight<-(15+ (ngmutgenes*0.1)))
                ggsave(file.path(outdir,paste0(prefix,  "_VAF_norm_histogram_per_Hugo_Symbol.pdf")),
                        vaf_histogram, width = 10, height=hheight)
                
                # Function to plot the VAF distribution per gene per Variant type
                vaf_histogram_pgpvt <- ggplot(maf, aes(x = VAF_norm, fill=Variant_Type)) +
                                geom_histogram(binwidth = 0.01, color = "black", alpha = 0.5) +
                                labs(title = paste0(prefix, "\n", "Histogram of VAF_norm Distribution per Hugo_Symbol"),
                                        x = "VAF_norm",
                                        y = "Frequency") +
                                theme_bw() +
                                scale_fill_paletteer_d("ggthemes::Color_Blind")+
                                theme(strip.text.y = element_text(size = 7, face= "bold"),
                                        legend.position = "right") +
                                facet_grid(Hugo_Symbol ~ ., scales = "free_y")

                # Save the histogram plot
                hheight<- NULL
                ifelse((15+ (ngmutgenes*0.1))>40, hheight<-40, hheight<-(15+ (ngmutgenes*0.1)))
                ggsave(file.path(outdir,paste0(prefix,  "_VAF_norm_histogram_per_Hugo_Symbol_pvartype.pdf")),
                        vaf_histogram_pgpvt, width = 10, height=hheight)
                        
                # Ridgeline plot
                vaf_ridgeline <- ggplot(maf, aes(x = VAF_norm, y = Hugo_Symbol, fill = ..x..)) +
                        geom_density_ridges_gradient(scale=0.5,rel_min_height = 0.01, 
                                jittered_points = TRUE, point_shape = "|",point_size = 5,
                                position = position_points_jitter(width = 0.04, height = 0, seed=42)) +
                        scale_fill_gradientn(colors = brewer.pal(11, "RdYlBu")) +
                        labs(title = paste0(prefix, "\n", "Ridgeline Plot of VAF_norm Distribution per Hugo_Symbol"),
                        x = "VAF_norm",
                        y = "Hugo Symbol",
                        fill = "VAF_norm") +
                        xlim(0,1) +
                        theme_ridges() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))

                # Save the ridgeline plot
                ggsave(file.path(outdir, paste0(prefix, "_VAF_norm_ridgeline_per_Hugo_Symbol.pdf")),
                vaf_ridgeline, width = 10, height = (20+ (ngmutgenes*0.1)))
        }

}



########### Variables
# Define the project directory and variables
projdir<-here()
message(paste0("Project directory: ", projdir))
resdir<-file.path(projdir, "cohort_analysis", "cohort_tables")
plotsdir<-file.path(projdir, "cohort_analysis", "cohort_tables","plots")
dir.create(plotsdir, showWarnings = FALSE, recursive = TRUE)


# Read the list of studies
studies<-fread(file.path(projdir, "metadata", "studies_reformat.tsv"), 
        header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(studies)<-c("Tumour_type","studyIDs" )
summary_numtab<-NULL

#Final cohort summary tables
cohort_maf_lof<-NULL
cohort_maf_lofwmis<-NULL
cohort_maf_nhs_lof<-NULL
cohort_maf_nhs_lofwmis<-NULL
# Loop through the studies
for(i in 1:nrow(studies)){
        study<-NULL
        tumour_type<-NULL
        maf_lof<-NULL
        maf_lofwmis<-NULL
        maf_nhs_lof<-NULL
        maf_nhs_lofwmis<-NULL
        #Assign the variables
        study<-studies$studyIDs[i]
        tumour_type<-studies$Tumour_type[i]
        message(paste0("Processing study: ", study))
        maf_lof<-fread(filecheck(file.path(projdir,study, "risk_analysis", 
                paste0("keepPA_vaf_size_filt_germline_", study,".risk.maf"))), header=TRUE, sep="\t")
        maf_lofwmis<-fread(filecheck(file.path(projdir,study, "risk_analysis", 
                paste0("keepPA_vaf_size_filt_germline_", study,".risk.wmiss.maf"))), header=TRUE, sep="\t")
        maf_nhs_lof<-fread(filecheck(file.path(projdir,study, "risk_analysis", 
                paste0("keepPA_vaf_size_filt_germline_", study,".nhs.lof.maf"))), header=TRUE, sep="\t")
        maf_nhs_lofwmis<-fread(filecheck(file.path(projdir,study, "risk_analysis", 
                paste0("keepPA_vaf_size_filt_germline_", study,".nhs.lof.wmiss.maf"))), header=TRUE, sep="\t")
        #Add tumour type to the maf tables
        maf_lof$Tumour_type<-tumour_type
        maf_lofwmis$Tumour_type<-tumour_type
        maf_nhs_lof$Tumour_type<-tumour_type
        maf_nhs_lofwmis$Tumour_type<-tumour_type

        #Count number of mutations
        all_genes_lof_nmut<-length(unique(paste(maf_lof$Hugo_Symbol, maf_lof$Chromosome, maf_lof$Reference_Allele, maf_lof$Consequence, maf_lof$Tumor_Sample_Barcode, sep="_")))
        all_genes_lofwmiss_nmut<-length(unique(paste(maf_lofwmis$Hugo_Symbol, maf_lofwmis$Chromosome, maf_lofwmis$Reference_Allele, maf_lofwmis$Consequence, maf_lofwmis$Tumor_Sample_Barcode, sep="_")))
        nhs_genes_lof_nmut<-length(unique(paste(maf_nhs_lof$Hugo_Symbol, maf_nhs_lof$Chromosome, maf_nhs_lof$Reference_Allele, maf_nhs_lof$Consequence, maf_nhs_lof$Tumor_Sample_Barcode, sep="_")))
        nhs_genes_lofwmiss_nmut<-length(unique(paste(maf_nhs_lofwmis$Hugo_Symbol, maf_nhs_lofwmis$Chromosome, maf_nhs_lofwmis$Reference_Allele, maf_nhs_lofwmis$Consequence, maf_nhs_lofwmis$Tumor_Sample_Barcode, sep="_")))

        #Add the tables to the final tables
        cohort_maf_lof<-rbind(cohort_maf_lof, maf_lof)
        cohort_maf_lofwmis<-rbind(cohort_maf_lofwmis, maf_lofwmis)
        cohort_maf_nhs_lof<-rbind(cohort_maf_nhs_lof, maf_nhs_lof)
        cohort_maf_nhs_lofwmis<-rbind(cohort_maf_nhs_lofwmis, maf_nhs_lofwmis)

        #Add the summary to the summary table
        summary_numtab<-rbind(summary_numtab, data.frame(Tumour_type=tumour_type, Study=study, 
                All_genes_lof_nmut=all_genes_lof_nmut, All_genes_lofwmiss_nmut=all_genes_lofwmiss_nmut, 
                NHS_genes_lof_nmut=nhs_genes_lof_nmut, NHS_genes_lofwmiss_nmut=nhs_genes_lofwmiss_nmut))
}
    

#Write the cohort summary tables
message("Writing the cohort summary tables")
fwrite(cohort_maf_lof, 
        file=file.path(resdir, paste0("FUR_keepPA_vaf_size_filt_germline_all_genes.risk.maf")),
        sep="\t", quote=FALSE, row.names=FALSE)
fwrite(cohort_maf_lofwmis, 
        file=file.path(resdir, paste0("FUR_keepPA_vaf_size_filt_germline_all_genes.risk.wmiss.maf")),
        sep="\t", quote=FALSE, row.names=FALSE)
fwrite(cohort_maf_nhs_lof, 
        file=file.path(resdir, paste0("FUR_keepPA_vaf_size_filt_germline.nhs.lof.maf")),
        sep="\t", quote=FALSE, row.names=FALSE)
fwrite(cohort_maf_nhs_lofwmis, 
        file=file.path(resdir, paste0("FUR_keepPA_vaf_size_filt_germline.nhs.lof.wmiss.maf")),
        sep="\t", quote=FALSE, row.names=FALSE)
# Write filterded cohort summary tables by VAF>=0.25
message("Writing the cohort summary tables filtered by VAF>=0.25")
fwrite(cohort_maf_lof[cohort_maf_lof$VAF_norm>=0.25,], 
        file=file.path(resdir, paste0("FUR_keepPA_vaf_size_filt_germline_all_genes.risk.vafn0.25filt.maf")),
        sep="\t", quote=FALSE, row.names=FALSE)
fwrite(cohort_maf_lofwmis[cohort_maf_lofwmis$VAF_norm>=0.25,],
        file=file.path(resdir, paste0("FUR_keepPA_vaf_size_filt_germline_all_genes.risk.wmiss.vafn0.25filt.maf")),
        sep="\t", quote=FALSE, row.names=FALSE)
fwrite(cohort_maf_nhs_lof[cohort_maf_nhs_lof$VAF_norm>=0.25,],  
        file=file.path(resdir, paste0("FUR_keepPA_vaf_size_filt_germline.nhs.lof.vafn0.25filt.maf")),
        sep="\t", quote=FALSE, row.names=FALSE)
fwrite(cohort_maf_nhs_lofwmis[cohort_maf_nhs_lofwmis$VAF_norm>=0.25,],
        file=file.path(resdir, paste0("FUR_keepPA_vaf_size_filt_germline.nhs.lof.wmiss.vafn0.25filt.maf")),
        sep="\t", quote=FALSE, row.names=FALSE)

###### Filter all variants for final table
#Variants with VAF_norm>=0.25
cohort_maf_nhs_loff<-cohort_maf_nhs_lof[cohort_maf_nhs_lof$VAF_norm>=0.25,]
# Variants with low depth to be removed, using 10% of the eveness of coverage used as threshold 6 reads as 60X
cohort_maf_nhs_loff_ldepf<-cohort_maf_nhs_loff[cohort_maf_nhs_loff$n_depth>6,]
fwrite(cohort_maf_nhs_loff_ldepf,  
        file=file.path(resdir, paste0("FUR_keepPA_vaf_size_filt_germline.nhs.lof.vafn0.25filt.lndepthf.maf")),
        sep="\t", quote=FALSE, row.names=FALSE)


# Write the summary table
fwrite(summary_numtab, 
        file=file.path(resdir, "cohort_summary_numtab.tsv"), sep="\t", quote=FALSE, row.names=FALSE)



# Plot the VAF distributions
message("Plotting the VAF distributions")
plot_vafdists(maf=cohort_maf_lof, prefix="cohort_maf_lof", outdir=file.path(plotsdir, "cohort_maf_lof"), plot_per_gene=FALSE)
plot_vafdists(maf=cohort_maf_lofwmis, prefix="cohort_maf_lofwmis", outdir=file.path(plotsdir, "cohort_maf_lofwmis"), plot_per_gene=FALSE)
plot_vafdists(maf=cohort_maf_nhs_lof, prefix="cohort_maf_nhs_lof", outdir=file.path(plotsdir, "cohort_maf_nhs_lof"))
plot_vafdists(maf=cohort_maf_nhs_lofwmis, prefix="cohort_maf_nhs_lofwmis", outdir=file.path(plotsdir, "cohort_maf_nhs_lofwmis"))

message("Plotting the VAF distributions of VAF>=0.25 filtered data")
plot_vafdists(maf=cohort_maf_lof[cohort_maf_lof$VAF_norm>=0.25,], 
        prefix="cohort_maf_lof_vaf_0.25filt", outdir=file.path(plotsdir, "cohort_maf_lof_vaf_0.25filt"), plot_per_gene=FALSE)
plot_vafdists(maf=cohort_maf_lofwmis[cohort_maf_lofwmis$VAF_norm>=0.25,], prefix="cohort_maf_lofwmis_vaf_0.25filt", outdir=file.path(plotsdir, "cohort_maf_lofwmis_vaf_0.25filt"), plot_per_gene=FALSE)
plot_vafdists(maf=cohort_maf_nhs_lof[cohort_maf_nhs_lof$VAF_norm>=0.25,], prefix="cohort_maf_nhs_lof_vaf_0.25filt", outdir=file.path(plotsdir, "cohort_maf_nhs_lof_vaf_0.25filt"))
plot_vafdists(maf=cohort_maf_nhs_lofwmis[cohort_maf_nhs_lofwmis$VAF_norm>=0.25,], prefix="cohort_maf_nhs_lofwmis_vaf_0.25filt", outdir=file.path(plotsdir, "cohort_maf_nhs_lofwmis_vaf_0.25filt"))
plot_vafdists(maf=cohort_maf_nhs_loff_ldepf, prefix="cohort_maf_nhs_lof_vaf_0.25filt_lndepthf", outdir=file.path(plotsdir, "cohort_maf_nhs_lof_vaf_0.25filt_lndepthf") )
message("Done")