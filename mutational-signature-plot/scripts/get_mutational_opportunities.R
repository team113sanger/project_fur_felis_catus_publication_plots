suppressMessages(library(BSgenome))
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))

# Create an opportunity matrix to be used as input to Sigfit

usage <- "\n\tThis script generates an oppportunities matrix for SigFit for a given BED file.

	Usage:

	 Rscript get_mutational_opportunities.R nsamples genome bedfile

	where 

	nsamples = the number of samples in your cohort
	genome is one of: canfam3, felcat9, mm10, grch38, bostau9
	bedfile is the path to a BED formatted file of target regions\n\n"


args <- commandArgs(trailingOnly = T)

nsamples <- args[1]
genome <- args[2]
bedfile <- args[3]

#available.genomes()

bases = c("A", "T", "C", "G", "a", "t", "c", "g")

first <- c("ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT")
second <- c("ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT")

# 96-element vector, ordered
context_order <- c(first, first, first, second, second, second)

# ---------- Functions ---------- #

get_trinuc_freq <- function(seqs) {
	trinuc_freq <- trinucleotideFrequency(seqs)
	sums <- colSums(trinuc_freq)
	df_sums <- data.frame(cbind(Context = names(sums), Sums = unname(sums)))
	# Get rev. comp where necessary and merge counts
	for (row in 1:nrow(df_sums)) {
		if (grepl("\\SG\\S", df_sums[row, 'Context']) || grepl("\\SA\\S", df_sums[row, 'Context'])) {
			df_sums[row, 'Context'] <- reverseComplement(DNAStringSet(df_sums[row, 'Context']))
		}
	}
	df_sums <- df_sums %>% group_by(Context) %>% summarise(Sums = sum(as.numeric(Sums))) %>% data.frame()
	return(df_sums)
}

make_opp_matrix <- function(df_sums, nsamples, output) {
	# Create a 96-element data frame, ordered
	df_sums_first <- df_sums %>% filter(Context %in% first) %>% arrange(match(Context, first))
	df_sums_second <- df_sums %>% filter(Context %in% second) %>% arrange(match(Context, first))

	context_96 <- rbind(df_sums_first, df_sums_first, df_sums_first, df_sums_second, df_sums_second, df_sums_second)
	print(context_96)
	opps <- matrix(rep(context_96$Sums, nsamples), nrow = as.numeric(nsamples), byrow = TRUE)
	colnames(opps) <- context_96$Context
	#print(opps)
	save(opps, file = paste0(output, "_opportunities.RData"))
}

# ---------- Main ---------- #

# Check the input parameters

if (is.na(nsamples) || is.na(genome) || is.na(bedfile)) {
	cat(usage)
	stop("Inputs missing")
	
} else if (! file.exists(bedfile)) {
	stop(paste("BED file", bedfile, "does not exist"))
}

# Load the requested genome (these must be pre-installed)

if (genome == "canfam3") {
	library(BSgenome.CanFam3.NCBI.canfam3.1)
} else if (genome == "felcat9") {
	library(BSgenome.Fcatus.NCBI.felCat9.0)
} else if (genome == "bostau9") {
	library(BSgenome.BosTau9.NCBI.bostau9)
} else if (genome == "mm10") {
	library(BSgenome.Mmusculus.UCSC.mm10)
} else if (genome == "grch38") {
	library(BSgenome.Hsapiens.NCBI.GRCh38)
} else {
	stop(paste("Genome", genome, "is not recognised"))
}

# Get tricleotide counts from BED file regions

bedfile <- read.table(bedfile, header = F, sep = "\t")
colnames(bedfile) <- c("Chromosome", "Start", "End")

# BED file is 0-based; add 1 to first coord
bedfile$Start <- bedfile$Start + 1

seqs <- DNAStringSet()

if (genome == "canfam3") {
	trinuc_sums <- data.frame()
	seqs <- getSeq(
			Cfamiliaris3.1,
			as.character(bedfile$Chromosome),
			start = bedfile$Start, 
			end = bedfile$End
		)
} else if (genome == "felcat9") {
	seqs <- getSeq(
			Fcatus9.0,
			as.character(bedfile$Chromosome),
			start = bedfile$Start, 
			end = bedfile$End
		)
} else if (genome == "bostau9") {
	print("Checking bostau9")
	seqs <- getSeq(
			Btaurus9,
			as.character(bedfile$Chromosome),
			start = bedfile$Start, 
			end = bedfile$End
		)
	str(seqs)
} else if (genome == "mm10") {
	seqs <- getSeq(
			Mmusculus,
			as.character(bedfile$Chromosome),
			start = bedfile$Start, 
			end = bedfile$End
		)
} else if (genome == "grch38") {
	seqs <- getSeq(
			Hsapiens,
			as.character(bedfile$Chromosome),
			start = bedfile$Start, 
			end = bedfile$End
		)
}

# Get the trinucleotide sums and make the opportunity matrix

trinuc_sums <- get_trinuc_freq(seqs)
make_opp_matrix(trinuc_sums, nsamples, genome)
print(trinuc_sums)



