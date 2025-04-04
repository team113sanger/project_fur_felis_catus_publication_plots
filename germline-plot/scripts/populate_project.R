library(tidyverse)
library(fs)
library(here)
library(glue)
library(jsonlite)


samples_to_keep <- fs::dir_ls("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_cat_cohort_files",   
                    recurse = TRUE, glob = "*samples_to_keep.nucleotide_variants.txt")

passed_samples <- read_tsv(samples_to_keep, col_names = c("sample"), id = "source") |> 
mutate(cohort = str_remove(basename(source),pattern = ".samples_to_keep.nucleotide_variants.txt")) |> 
select(-source) 




 tmp <- passed_samples |> 
 filter(!grepl(sample, pattern= "a")) |> 
 group_by(cohort) |> 
 group_split()
tmp[[1]] |> print(n = 100)


study_paths <- fs::dir_ls("/lustre/scratch125/casm/team113da/users/bf14/samples", 
type = "directory") |>
  as_tibble() |>
  mutate(study_id = basename(value)) |>
  filter(grepl(study_id, pattern = "[0-9]+")) |>
  mutate(CanApps_ID = as.numeric(study_id))


project_dir <- here()
studies <- read_tsv("metadata/study_manifest.tsv", )
studies

study_info <- studies |>
  left_join(study_paths, by = c("CanApps_ID")) |>
  rowwise() |>
  mutate(cohort_id = glue("{SeqScape_ID}_{CanApps_ID}")) |>
  mutate(path = glue("{project_dir}/{cohort_id}"))

study_info[["path"]]

output_dirs <- study_info |> pull(path)
walk(output_dirs, dir.create)
walk(paste0(output_dirs, "/results"), dir.create)

file_ids <- paste0(output_dirs, "/params.json")



items <- map(study_info[["value"]], ~dir_ls(.x, recurse = TRUE, glob = "*.bam"))


sample_manifests <- tibble(object = purrr::flatten(items)) |> 
unnest(cols = c(object)) |>
mutate(study = str_match(object, pattern = ".*/samples/([0-9]+)/CATD*")[,2], .before = "object") |> 
mutate(sample = str_match(object, pattern = ".*/samples/.*/(CATD.*)/CAT")[,2], .before = "object") |> 
mutate(object_index = paste0(object, ".bai")) 

sample_manifests

grouped_manifests <- sample_manifests |> 
left_join(study_info, by = c("study" = "study_id")) |> 
group_by(path) |>
filter(!grepl(sample, pattern = "a")) |>
select(sample, object, object_index)  

manifests <- group_split(grouped_manifests, .keep = FALSE)
names(manifests) <- group_keys(grouped_manifests)[["path"]]


filtered_manifests <- map2(manifests, tmp[1:13], ~filter(.x, sample %in% .y$sample))



# map2(filtered_manifests, basename(names(filtered_manifests)),
#      ~write_tsv(.x, file = paste0(.y, "/sample_manifest.tsv")))



map_files <- fs::dir_ls("/lustre/scratch127/casm/team113da/projects/fur_germline",   
                    recurse = TRUE, glob = "*results/sample_map.txt")

original_set <- map(map_files, read_tsv, col_names = c("sample", "object"))


filter_missing <- function(x,y){
differences <- setdiff(x[["sample"]], y[["sample"]])
x |> filter(!sample %in% differences)}



order_1 <- basename(dirname(dirname(names(original_set))))
order_2 <- map_chr(tmp[1:13], ~select(.x, cohort) |> unique() |> pull(cohort))
sum(order_1 == order_2)

filtered_sample_maps <- map2(original_set,tmp[1:13], filter_missing) 
imap(filtered_sample_maps, ~write_tsv(.x, file = .y,col_names = FALSE))


parameters <- map2(
  output_dirs, as.list(study_info[["cohort_id"]]),
  ~ toJSON(list(
      study_id = .y,
      chrom_list = "/lustre/scratch127/casm/team113da/projects/fur_germline/metadata/cat_contigs.txt",
      post_process_only = TRUE,
      summarise_results = FALSE,
      tsv_file  = glue("{project_dir}/{.y}/sample_manifest.tsv"),
      samples_to_process  = -1,
      run_mode  = "sort_inputs",
      run_markDuplicates  = TRUE,
      run_coord_sort_cram  = TRUE,
      run_deepvariant  = FALSE,
      run_haplotypecaller  = TRUE,
      reference_genome  = "/lustre/scratch124/casm/team78pipelines/reference/Cat/Felis_catus_9.0/genome.fa",
      baitset = "/lustre/scratch125/casm/team113da/projects/FUR/FUR_bases/FUR_base_cat/metadata/references/baitset/DNA/S3250994_Feline_HSA_Jan2020_146.bed",
      vep_cache = "/lustre/scratch125/casm/team113da/users/bf14/git/ensembl-vep/cache",
      custom_files = "/lustre/scratch125/casm/team113da/users/bf14/variant_caller_benchmarking/VEP/SRA_99_Lives.snp_indel.filt.54cats.short.alt.no_header.vcf.gz{,.tbi}",
      custom_args = "SRA_99_Lives.snp_indel.filt.54cats.short.alt.no_header.vcf.gz,99_Lives,vcf,exact,0,AF",
      species = "felis_catus",
      db_version = "104",
      assembly = "Felis_catus_9.0",
      outdir  = glue("{project_dir}/{.y}/results")
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  )
)
file_ids <- paste0(output_dirs, "/cohort_params.json")
map2(parameters, file_ids, ~ write_file(x = .x, file =   .y))


