library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

maf1_file <- args[1]
maf2_file <- args[2]
cat1 <- args[3]
cat2 <- args[4]
outfile <- args[5]

# Read in files

maf1 <- read.table(maf1_file, header = T, sep = "\t", comment.char = "")
maf2 <- read.table(maf2_file, header = T, sep = "\t", comment.char = "")

# Check MAF headers
headers1 <- colnames(maf1)
headers2 <- colnames(maf2)

if (!identical(headers1, headers2)) {
#    stop("MAF hedaers do not match")
    print("Column headers differ; using common columns only")
    cols_to_use <- intersect(headers1, headers2)
    maf1 <- maf1 %>% select(matches(cols_to_use))
    maf2 <- maf2 %>% select(matches(cols_to_use))
#    print(colnames(maf1))
#    print(colnames(maf2))
#    print(cols_to_use)
}

# Create a 'Patient' column

maf1 <- maf1 %>%
    mutate("Patient" = str_sub(Tumor_Sample_Barcode, end = -2))

maf2 <- maf2 %>%
    mutate("Patient" = str_sub(Tumor_Sample_Barcode, end = -2))

# Get common genes

common_genes_1 <- maf1 %>%
    semi_join(maf2, by = c("Hugo_Symbol", "Gene", "Patient")) %>%
    mutate("Category" = cat1)

common_genes_2 <- maf2 %>%
    semi_join(maf1, by = c("Hugo_Symbol", "Gene", "Patient")) %>%
    mutate("Category" = cat2)

# Print out appended maf

merged_maf <- rbind(common_genes_1, common_genes_2)

write.table(merged_maf, file = outfile, row.names = F, quote = F, sep = "\t")

