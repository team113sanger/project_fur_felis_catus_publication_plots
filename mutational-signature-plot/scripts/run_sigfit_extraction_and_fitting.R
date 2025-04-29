.libPaths(c("/software/team113/dermatlas/R/R-4.2.2/lib/R/library", 
			"/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_mutational_signatures/renv/library/R-4.2/x86_64-pc-linux-gnu"))

library(sigfit)
library(Cairo)
library(dplyr)
library(GetoptLong)

# This script take a catalogue of variants and runs
# Sigfit to extract signatures.
# The catalogue can be created using the script 'make_SNV_Sigfit_table.R'
# Rscript make_SNV_Sigfit_table.R your_file.maf GRCh38


############### Functions ###############

### Rerun extraction using many more iterations to get better results

rerun_extraction <- function(mut_catalogue, best) {

	# Run extraction 
	samples_extr_best = extract_signatures(counts = mut_catalogue,
											nsignatures = best,
											seed = 9242,
											iter = 20000)
#	print("Best estimate")
#	head(samples_extr_best)

	# Retrieve signatures
	signatures_rerun <- retrieve_pars(samples_extr_best, "signatures")

	# "Re-fit these signatures back to the original catalogues,
	# as this sometimes results in exposure estimates that are
	# more accurate than the ones obtained by signature extraction."

	print("Fitting")
	samples_extr_best_refit <- fit_signatures(counts = mut_catalogue,
										 signatures = signatures_rerun,
										 seed = 1844,
										 iter = 20000)

	# Write outputs to tables
	print("Getting signature")
	genome_signatures <- samples_extr_best_refit$data$signatures

	print("Getting exposure and recon")
	exposures <- retrieve_pars(samples_extr_best_refit, "exposures")
	reconstructions <- retrieve_pars(samples_extr_best_refit, "reconstructions")

	write_to_tsv(genome_signatures, paste0("extraction_genome_signatures_", best, "sigs.tsv"))
	write_to_tsv(exposures$mean, paste0("extraction_signature_exposures_", best, "sigs.tsv"))
	write_to_tsv(reconstructions$mean, paste0("extraction_reconstructed_signatures_", best, "sigs.tsv"))

	# plot everything relative to exome
	plot_all(samples_extr_best_refit, exp_cex_names = 0.9, out_path = ".")
	save(list = ls(all.names = TRUE), file = paste0("sigfit_extraction", best, "sigs.RData"))

}


### Write to tsv file

write_to_tsv <- function(df, filename) {
	write.table(df, filename, quote = F, sep = "\t", row.names = T, col.names = NA)
}

 
### Write to tsv file, with no rownames

write_to_tsv_norownames <- function(df, filename) {
	write.table(df, filename, quote = F, sep = "\t", row.names = F)
}


### Calculate cosine similarity

cosine <- function (x, y) {
	# x and y are vectors
	a <- sum( x * y)
	b <- sqrt(sum(x * x))
	c <- sqrt(sum(y * y))

	cos.sim <- a / ( b * c)
	return(cos.sim)
}


### Run cosine similarity check

get_cosine_sim <- function(file1, file2, ref_sigs) {
	query <- read.table(file1, header = T, sep = "\t", row.names = 1, check.names = F)
	ref <- read.table(file2, header = T, sep = "\t", row.names = 1, check.names = F)

	# Transpose the data frames
	query <- t(query)
	ref <- t(ref)

	# Sort to have motifs in the same order in both files
	order <- rownames(query)

	query <- data.frame(query[order, ])
	ref <- ref[order, ]

	# Print signatures to file
	write_to_tsv(query, "query_signatures.tsv")
	write_to_tsv(ref, paste0(ref_sigs, "_signatures.txt"))

	# Run each pairwise comparison. Compare each signature in query
	# to each signature in the reference

	sim_df <- data.frame(Query = character(), Reference = character(), Cosine_sim = numeric())
	highest_sim <- data.frame(Query = character(), Reference = character(), Cosine_sim = numeric())
	index <- 0

	for (i in 1:ncol(query)) {
		x <- query[ ,i]
		name <- colnames(query)[i]
		#print(x)
		vals <- c()
		for (j in 1:ncol(ref)) {
			y <- ref[ ,j]
			sim <- rbind(cosine(x, y))
			cat(paste("Sim is ", sim, " ", colnames(ref)[j]), "\n")
			index <- index + 1
			sim_df[index, ] <- c(name, colnames(ref)[j], sim)
			vals <- c(vals, sim)	
		}

		# Get the highest cosine similarity
		# If the highest and second highest are within 0.03, print out both

		tmp <- data.frame(vals,colnames(ref)[1:ncol(ref)])
		tmp <- tmp[order(tmp$vals, decreasing = TRUE),]
		if (tmp[1, 1] - 0.03 < tmp[2, 1]) {
			cat(paste("Max is ",tmp[1, 1]," and ", tmp[2, 1]),"\n")
			cat(paste(name, tmp[1, 2], tmp[2, 2], tmp[1, 1], tmp[2, 1], sep = "\t"),"\n")
			highest_sim[i, ] <- c(name, tmp[1, 2], tmp[1,1], tmp[2, 2], tmp[2, 1])

		} else {
			cat(paste("Max is ", tmp[1 ,1]), "\n")
			cat(paste(name, tmp[1, 2], tmp[1, 1], sep = "\t"), "\n")
			highest_sim[i, ] <- c(name, tmp[1, 2], tmp[1,1], "N/A", "N/A")
		}
	}

	save(list = ls(all.names = TRUE), file = paste0("cosine_sim_", ncol(query), "sigs.RData"))
	# Write results to file

	write_to_tsv_norownames(sim_df, paste0("cosine_sim_all_", ref_sigs, ".tsv"))
	write_to_tsv_norownames(highest_sim, paste0("cosine_sim_highest_", ref_sigs, ".tsv"))
}


### Signature fitting

run_sig_fitting <- function(muts, signatures, sample, prefix) {

	signatures <- signatures
	sample <- sample
	prefix <- prefix
	# Run fitting against specific signatures
	samples_fit <- fit_signatures(counts = muts,
										 signatures = signatures,
										 seed = 48573,
										 iter = 20000)

	# Retreive exposures

	exposures <- retrieve_pars(samples_fit, "exposures")
	#reconstructions <- retrieve_pars(samples_fit, "reconstructions")
	# Plot spectrum, exposures and reconstructions
	plot_all(mcmc_samples = samples_fit, 
		out_path = ".",
		prefix = paste0(prefix, "_fitting_", sample))

	# Write out exposures to a file
	write_to_tsv(exposures$mean, paste0(prefix, "_fitting_exposures_", sample, ".tsv"))
	#write_to_tsv(reconstructions$mean, paste0(prefix, "_fitting_recontructions_", sample, ".tsv"))

	# Now run refitting base on signatures with mean exposure > 0.1
	# Get index of signatures present with exposure > 0.1
	exp <- t(exposures$mean)
	exp_index <- which(exp > 0.1)
	if (length(exp_index) < 2) {
		sig <- rownames(exp)[exp_index]
	   print(paste("Not refitting for", sample, ", 1 signature present:", sig))
	} else {
		samples_fit_2 <- fit_signatures(counts = muts,
											 signatures = signatures[exp_index, ],
											 seed = 48573,
											 iter = 20000)

		exposures_2 <- retrieve_pars(samples_fit_2, "exposures")
		#reconstructions_2 <- retrieve_pars(samples_fit_2, "reconstructions")

		plot_all(mcmc_samples = samples_fit_2, 
			out_path = ".",
			prefix = paste0(prefix, "_Refitting_", sample))

		# Write out exposures to a file
		write_to_tsv(exposures_2$mean, paste0(prefix, "_Refitting_exposures_", sample, ".tsv"))
		#write_to_tsv(reconstructions_2$mean, paste0(prefix, "_Refitting_recontructions_", sample, ".tsv"))
		# Save Rdata
		save(list = ls(all.names = TRUE), file = paste0(prefix, "_sigfit_fitting_", sample, ".RData"))
	}
}



############### Main ###############

# Get command line paramaters

opp_file <- NA
opps <- NA

spec = "

Run SigFit with a custom refdb and optional covariates file.

  Usage: Rscript run_sigfit_extraction_and_fitting.R [options]

  Options:
    <file=s> SigFit variants file with sequence context [required]
    <output_dir=s> Output directory [required]
    <cosmic_file=s> Cosmic signature file [required]
    <signal_file=s> Signal signature file [required]
    <opp_file=s> Opportunities matrix file [optional]

"

GetoptLong(spec, template_control = list(opt_width = 21))


# Check that file exists and output directory exists

for (f in c(file, cosmic_file, signal_file)) {
	if (! file.exists(f)) {
		stop(paste("File", f, "does not exist"))
	}
}

if (! dir.exists(output_dir)) {
	stop(paste("Directory", output_dir, "does not exist"))
} else if (! is.na(opp_file) && ! file.exists(opp_file)) {
	stop(paste("Opportunity file", opp_file, "does not exist"))
}

setwd(output_dir)


# Read the input file

mutations <- read.table(file, header = T, sep = "\t", check.names = F)
mutations_nostrand <- mutations %>% select(-Strand)

# Load the opportunity file if provided

if (! is.na(opp_file)) {
	load(opp_file)
	print("Loaded opportunity file")
	head(opps)
}

# Build mutation catalogues, 96 trinucleotide types and 192 types (including strand info)

catalogue_exome <- build_catalogues(mutations)
catalogue_nostrand_exome <- build_catalogues(mutations_nostrand)

print("Catalogue file read")

# Get samples with minimum 100 mutations
catalogue_exome_100 <- as.data.frame(catalogue_nostrand_exome) %>% filter(rowSums(.) >= 100) %>% as.matrix()
if (nrow(catalogue_exome_100) == 0) {
	print("NO SAMPLES WITH AT LEAST 100 MUTATIONS, SKIPPING SIGNATURE FITTING.")
}


# Write out catalog files from exome

write_to_tsv(catalogue_exome, "exome_catalogue_192.tsv")
write_to_tsv(catalogue_nostrand_exome, "exome_catalogue_96.tsv")
write_to_tsv(catalogue_exome_100, "exome_catalogue_96_min100.tsv")


# Get the number of samples to adjust mfrow and pdf height

num_samples <- length(unique(mutations$Sample))
print(paste("Number of samples is", num_samples))


# Plot the overall spectrum
pdf("snv_spectrum_192_exome.pdf", width = 20, height = num_samples * 5)
par(mfrow = c(num_samples, 1), oma = c(7,7,7,2), mar = c(5,5,8,5))
plot_spectrum(catalogue_exome)
dev.off()

pdf("snv_spectrum_96_exome.pdf", width = 20, height = num_samples * 5)
par(mfrow = c(num_samples, 1), oma = c(7,7,7,2), mar = c(5,5,8,5))
plot_spectrum(catalogue_nostrand_exome)
dev.off()


# Convert catalogue to human genome now so we can compare to COSMIC genomes and
# get exposures later
print("Converting 192")
if (is.na(opp_file)) {
	catalogue <- convert_signatures(catalogue_exome,
					opportunities_from="human-exome",
					opportunities_to="human-genome")
} else {
	print("Using custom opportunities file")
	catalogue <- convert_signatures(catalogue_exome,
					opportunities_from=opps,
					opportunities_to="human-genome")
}

catalogue <- catalogue  * rowSums(catalogue_exome)

print("Converting 96")
if (is.na(opp_file)) {
	catalogue_nostrand <- convert_signatures(catalogue_nostrand_exome,
					opportunities_from="human-exome",
					opportunities_to="human-genome")
} else {
	print("Using custom opportunities file")
	catalogue_nostrand <- convert_signatures(catalogue_nostrand_exome,
					opportunities_from=opps,
					opportunities_to="human-genome")
}

catalogue_nostrand <- catalogue_nostrand  * rowSums(catalogue_nostrand_exome)

if (nrow(catalogue_exome_100) > 1) {
	print("Converting 96 min100")
	if (is.na(opp_file)) {
		catalogue_100 <- convert_signatures(catalogue_exome_100,
						opportunities_from="human-exome",
						opportunities_to="human-genome")
	} else {
		print("Using custom opportunities file")
		catalogue_100 <- convert_signatures(catalogue_exome_100,
						opportunities_from=opps,
						opportunities_to="human-genome")
	}

	catalogue_100 <- catalogue_100 * rowSums(catalogue_exome_100)
	write_to_tsv(catalogue_100, "genome_converted_catalogue_96_min100.tsv")
}

write_to_tsv(catalogue, "genome_converted_catalogue_96.tsv")
write_to_tsv(catalogue_nostrand, "genome_converted_catalogue_192.tsv")

samples_extr_nostrand = extract_signatures(counts = catalogue_nostrand,
									nsignatures = 1:8,
									seed = 3485,
									iter = 5000)


cosine_gof_96 <- data.frame(calculate_gof(samples_extr_nostrand, stat = "cosine"))
l2_gof_96 <- data.frame(calculate_gof(samples_extr_nostrand, stat = "L2"))

write_to_tsv_norownames(cosine_gof_96, "gof_cosine_96.tsv")
write_to_tsv_norownames(l2_gof_96, "gof_L2_96.tsv")


### Draw a goodness-of-fit plot

pdf("gof_plot_cosine_96.pdf", width = 8, height = 8)
plot_gof(samples_extr_nostrand, stat = "cosine")
dev.off()

pdf("gof_plot_L2_96.pdf", width = 8, height = 8)
plot_gof(samples_extr_nostrand, stat = "L2")
dev.off()



# Retrieve the best estimate, and also run best-1 and best+1 signatures

best <- cosine_gof_96$best[1]
print(c("BEST cosine:", best))

save.image(file = "sigfit_extraction.RData")

which = 1

for (numsigs in c(best, best + 1, best + 2)) {
	dirname_suffix = "sigs"
	if (which == 1) {
		dirname_suffix = paste0(dirname_suffix, "_best_estimate")
	}
	sig_output_dir <- file.path(output_dir, paste0(numsigs, dirname_suffix))
	dir.create(file.path(sig_output_dir))
	setwd(file.path(sig_output_dir))

	rerun_extraction(catalogue_nostrand, numsigs)

	# Calculate cosine similarity with known signatures

	get_cosine_sim(paste0("extraction_genome_signatures_", numsigs, "sigs.tsv"), cosmic_file, "cosmic")
	get_cosine_sim(paste0("extraction_genome_signatures_", numsigs, "sigs.tsv"), signal_file, "signal")
	which <- which + 1
}

save.image(file = "sigfit_extraction.RData")


# Run signature fitting on samples with >= 100 mutations and input cosmic file

if (exists("catalogue_100")) {
	for (name in row.names(catalogue_100)) {
		catalogue_subset <- as.data.frame(catalogue_100) %>% filter(rownames(catalogue_100) == name) %>% as.matrix()
		sample <- row.names(catalogue_subset)[1]

		fit_output_dir <- file.path(output_dir, paste0("sig_fitting_", sample))
		dir.create(file.path(fit_output_dir))
		setwd(file.path(fit_output_dir))

		ref_sigs_cosmic <- read.table(cosmic_file, header = T, sep = "\t", row.names = 1, check.names = F)
		run_sig_fitting(catalogue_subset, ref_sigs_cosmic, sample, "cosmic")

		ref_sigs_signal <- read.table(signal_file, header = T, sep = "\t", row.names = 1, check.names = F)
		run_sig_fitting(catalogue_subset, ref_sigs_signal, sample, "signal")
	}
}

save.image(file = "sigfit_extraction.RData")

print("Done!")
