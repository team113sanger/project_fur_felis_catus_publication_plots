suppressMessages(library(BSgenome))
suppressMessages(library(dplyr))
suppressMessages(library(readr))

# Input format for multi-sample comparison in Sigfit (matrix or data frame)
#
#      Sample     Ref Alt Trinuc
# [1,] "Sample 1" "C" "T" "TCT" 
# [2,] "Sample 1" "T" "C" "ATT" 
# [3,] "Sample 1" "C" "T" "ACT" 
# [4,] "Sample 1" "C" "T" "CCG" 

args <- commandArgs(trailingOnly = T)

file <- args[1]
genome <- args[2]

bases = c("A", "T", "C", "G", "a", "t", "c", "g")

#maf_table <- read.table(file, sep = "\t", header = T, check.names = F)
maf_table <- read_tsv(file)

maf_table <- maf_table %>% select(Tumor_Sample_Barcode, Reference_Allele, Tumor_Seq_Allele2, Chromosome, Start_Position, TRANSCRIPT_STRAND) %>%
				mutate(Pos = Start_Position - 1, .before = "TRANSCRIPT_STRAND") %>% filter(Reference_Allele %in% bases & Tumor_Seq_Allele2 %in% bases)

head(maf_table)


if (genome == "canfam3") {
	library(BSgenome.CanFam3.NCBI.canfam3.1)
	ref = "Cfamiliaris3.1"
	print(Cfamiliaris3.1)
	seqs <- getSeq(
			Cfamiliaris3.1, 
			as.character(maf_table$Chromosome),
			start = maf_table$Pos, 
			width = 3
		)
} else if (genome == "felcat9") {
	library(BSgenome.Fcatus.NCBI.felCat9.0)
	seqs <- getSeq(
		Fcatus9.0,
		as.character(maf_table$Chromosome),
		start = maf_table$Pos, 
		width = 3
	)
} else if (genome == "bostau9") {
	print("Getting bostau9 seqs")
	library(BSgenome.BosTau9.NCBI.bostau9)
	seqs <- getSeq(
			Btaurus9,
			as.character(maf_table$Chromosome),
			start = maf_table$Pos, 
			width = 3
		)
} else if (genome == "mm10") {
	library(BSgenome.Mmusculus.UCSC.mm10)
	seqs <- getSeq(
			Mmusculus,
			as.character(maf_table$Chromosome),
			start = maf_table$Pos, 
			width = 3
		)
} else if (genome == "grch38") {
	library(BSgenome.Hsapiens.NCBI.GRCh38)
	maf_table$Chromosome <- sub("chr", "", maf_table$Chromosome)
	print(head(maf_table))
	seqs <- getSeq(
			Hsapiens,
			as.character(maf_table$Chromosome),
			start = maf_table$Pos, 
			width = 3
		)
}


# Get the mutation context (1bp on either side)


# Print out table for Sigfit

var_table <- cbind(maf_table, seqs) %>% select(Tumor_Sample_Barcode, Reference_Allele, Tumor_Seq_Allele2, seqs, TRANSCRIPT_STRAND)
colnames(var_table) <- c("Sample", "Ref", "Alt", "Trinuc", "Strand")

write.table(var_table, file = "sigfit_var_table.strand192.tsv", quote = F, sep = "\t", row.names = F)