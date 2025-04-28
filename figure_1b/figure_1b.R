###############################################
## 0) Load Necessary Libraries
###############################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(grid)

###############################################
## 1) Define Parameters and MAF File Paths
###############################################
top_n <- 5  # Number of top genes to display per tumor type

maf_files <- list.files(
  path       = "./data", 
  pattern    = "*.maf", 
  full.names = TRUE
)

tumor_types <- c(
  "Lung adenocarcinoma", "Cutaneous squamous cell carcinoma",
  "Oral squamous cell carcinoma", "Mast cell tumour",
  "Meningioma", "Pancreatic carcinoma",
  "Cholangiocarcinoma", "Osteosarcoma", "Lymphoma",
  "Mammary gland carcinoma", "Basal cell carcinoma",
  "Colorectal carcinoma", "Glioma"
)

###############################################
## 2) Function: Calculate Top N Genes per Tumor Type
###############################################
calculate_top_n_genes <- function(maf_file, tumor_type, top_n) {
  maf_data <- read.delim(maf_file, stringsAsFactors = FALSE)
  
  mutation_proportion <- maf_data %>%
    filter(Hugo_Symbol != "-") %>%
    group_by(Hugo_Symbol) %>%
    summarise(Mutated_Samples = n()) %>%
    mutate(Tumor_Type = tumor_type) %>%
    arrange(desc(Mutated_Samples)) %>%
    slice_head(n = top_n)
  
  return(mutation_proportion)
}

###############################################
## 3) Get Top Genes Across All Tumor Types
###############################################
top_genes_data <- map2_df(maf_files, tumor_types, 
                          ~calculate_top_n_genes(.x, .y, top_n))

# Unique list of top genes across all tumor types
unique_top_genes <- unique(top_genes_data$Hugo_Symbol)

###############################################
## 4) Calculate Variant Classifications (mutation_types_data)
###############################################
# Read all MAF files, filter for our unique_top_genes, 
# and count unique samples by gene & mutation type.

calculate_mutation_types <- function(maf_file) {
  maf_data <- read.delim(maf_file, stringsAsFactors = FALSE)
  
  mutation_types <- maf_data %>%
    filter(Hugo_Symbol %in% unique_top_genes) %>%
    group_by(Hugo_Symbol, Variant_Classification) %>%
    summarise(
      Mutated_Samples = n_distinct(Tumor_Sample_Barcode),
      .groups         = "drop"
    )
  
  return(mutation_types)
}

# Combine all MAF data
maf_data_combined <- map_df(maf_files, read.delim, stringsAsFactors = FALSE)

# Total unique samples in the cohort
total_samples_cohort <- maf_data_combined %>%
  pull(Tumor_Sample_Barcode) %>%
  n_distinct()

# Summarize mutation types across all files
mutation_types_data <- map_df(maf_files, calculate_mutation_types) %>%
  group_by(Hugo_Symbol, Variant_Classification) %>%
  summarise(
    Total_Mutations   = sum(Mutated_Samples),
    .groups           = "drop"
  ) %>%
  mutate(Sample_Proportion = Total_Mutations / total_samples_cohort)

###############################################
## 5) Determine Gene Order 
###############################################
gene_order_by_total_muts <- mutation_types_data %>%
  group_by(Hugo_Symbol) %>%
  summarise(Total_Mutations = sum(Total_Mutations, na.rm = TRUE)) %>%
  arrange(desc(Total_Mutations)) %>%
  pull(Hugo_Symbol)

###############################################
## 6) Calculate Heatmap Data (gene x tumor type)
###############################################
calculate_mutation_proportion_for_gene_list <- function(maf_file, tumor_type, gene_list) {
  maf_data <- read.delim(maf_file, stringsAsFactors = FALSE)
  
  mutation_proportion <- maf_data %>%
    filter(Hugo_Symbol %in% gene_list) %>%
    group_by(Hugo_Symbol) %>%
    summarise(Mutated_Samples = n_distinct(Tumor_Sample_Barcode)) %>%
    mutate(
      Mutation_Proportion = Mutated_Samples / n_distinct(maf_data$Tumor_Sample_Barcode),
      Tumor_Type          = tumor_type
    ) %>%
    complete(Hugo_Symbol = gene_list, fill = list(Mutation_Proportion = 0))
  
  return(mutation_proportion)
}

# Mutation proportions for all top genes across all tumor types
mutation_data <- map2_df(maf_files, tumor_types, 
                         ~calculate_mutation_proportion_for_gene_list(.x, .y, unique_top_genes))

# Make a consistent data frame for the heatmap
heatmap_data <- mutation_data %>%
  rename(Gene = Hugo_Symbol) %>%
  complete(Tumor_Type, Gene, fill = list(Mutation_Proportion = NA)) %>%
  filter(!is.na(Tumor_Type))

###############################################
## 7) Pivot Heatmap Data to a Wide Matrix
###############################################
# rows = Gene, columns = Tumor_Type, values = Mutation_Proportion
heatmap_matrix <- heatmap_data %>%
  select(Gene, Tumor_Type, Mutation_Proportion) %>%
  pivot_wider(names_from = Tumor_Type, values_from = Mutation_Proportion) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

###############################################
## 8) Prepare Stacked Bars Data
###############################################
# 8a) Pivot mutation_types_data to wide form: 
#     rows = gene, columns = variant classification, value = Sample_Proportion
mutation_types_data$Hugo_Symbol <- factor(mutation_types_data$Hugo_Symbol, 
                                          levels = unique_top_genes) 
stacked_df <- mutation_types_data %>%
  select(Hugo_Symbol, Variant_Classification, Sample_Proportion) %>%
  pivot_wider(
    names_from  = Variant_Classification,
    values_from = Sample_Proportion,
    values_fill = 0
  )

# 8b) Turn into a matrix
stacked_mat <- as.matrix(stacked_df[, -1])
rownames(stacked_mat) <- stacked_df$Hugo_Symbol

# 8c) Colors for variant classes
variant_colors <- c(
  "Frame_Shift_Del"        = "orange", 
  "Frame_Shift_Ins"        = "brown", 
  "In_Frame_Del"           = "blue",
  "In_Frame_Ins"           = "red",
  "Missense_Mutation"      = "springgreen3",
  "Nonsense_Mutation"      = "black",
  "Nonstop_Mutation"       = "pink",
  "Splice_Site"            = "purple",
  "Targeted_Region"        = "cyan",
  "Translation_Start_Site" = "yellow"
)

# Keep only columns actually present in stacked_mat
variant_order <- intersect(colnames(stacked_mat), names(variant_colors))
variant_colors <- variant_colors[variant_order]

###############################################
## 9) Reorder Genes by Total Stacked Height 
##    (Largest sum -> top row)
###############################################
# 9a) Compute row sums
row_totals <- rowSums(stacked_mat, na.rm = TRUE)

# 9b) Sort genes by decreasing total
gene_order_by_stack <- names(sort(row_totals, decreasing = TRUE))

# 9c) Reorder stacked_mat and heatmap_matrix
#     so that row 1 is the highest stacked bar
stacked_mat     <- stacked_mat[gene_order_by_stack, , drop = FALSE]
heatmap_matrix  <- heatmap_matrix[gene_order_by_stack, , drop = FALSE]

###############################################
## 10) Custom Annotation for Stacked Bars
###############################################
my_stacked_bar_fun <- function(index) {
  n <- length(index)
  pushViewport(viewport(xscale = c(0, 0.4), yscale = c(0.5, n + 0.5)))
  
  for (i in seq_along(index)) {
    row_i <- index[i]
    y <- n - i + 1  
    
    gene_name  <- rownames(stacked_mat)[row_i]
    row_values <- stacked_mat[gene_name, ]
    
    x_left <- 0
    for (class in variant_order) {
      portion  <- row_values[class]
      x_right  <- x_left + portion
      
      grid.rect(
        x      = unit((x_left + x_right) / 2, "native"),
        y      = unit(y, "native"),
        width  = unit(x_right - x_left, "native"),
        height = unit(0.8, "native"),
        gp     = gpar(fill = variant_colors[class], col = NA)
      )
      x_left <- x_right
    }
  }
  
  # Add x-axis
  grid.xaxis(
    at = seq(0, 0.4, by = 0.1), 
    label = seq(0, 0.4, by = 0.1),
    gp = gpar(fontsize = 8),
    main = TRUE
  )
  
  popViewport()
}


# Legend for stacked bars
lgd_varclass <- Legend(
  labels    = variant_order,
  legend_gp = gpar(fill = variant_colors),
  title     = "Alteration Type"
)

###############################################
## 11) Build the Row Annotation 
###############################################
ha_stacked <- rowAnnotation(
  "." = AnnotationFunction(
    fun   = my_stacked_bar_fun,
    width = unit(4, "cm"),
    which = 'row'
  )
)

###############################################
## 12) Build the Main Heatmap
###############################################
# Define a color ramp for the heatmap
col_fun <- colorRamp2(
  breaks = c(0.001, 0.5, 1), 
  colors = c("lightblue1", "cyan3", "blue4")
)

ht_main <- Heatmap(
  matrix          = heatmap_matrix,
  name            = "Mutation\nProportion",
  col             = col_fun,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  show_row_dend   = FALSE,
  show_column_dend= FALSE,
  row_names_side  = "left",
  column_title    = "Tumor Type",
  row_title       = "Gene",
  row_names_gp    = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  height          = unit(15, "cm"),
  width           = unit(8, "cm"),
  rect_gp         = gpar(col = "black", lwd = 0.5),
  heatmap_legend_param = list(
    title = "Mutation Proportion",
    at = c(0, 0.5, 1),
    labels = c(">0", "0.5", "1")
  )
)

###############################################
## 13) Combine the Heatmap + Stacked Bars, Draw
###############################################
ht_with_stacked <- ht_main + ha_stacked

draw(
  ht_with_stacked,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE,
  annotation_legend_list = list(lgd_varclass)
)

write.table(heatmap_matrix, file = "heatmap_matrix.csv", sep = ",", row.names = TRUE, col.names = TRUE)
