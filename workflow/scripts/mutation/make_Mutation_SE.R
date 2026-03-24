snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(fs)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(dplyr)
})

message("Loading mutation calls: ", snk$input$raw)
mut_df <- read_sheet_tsv(snk$input$raw) |>
  normalize_character_columns() |>
  filter(
    !is.na(sample),
    nzchar(sample),
    !is.na(gene),
    nzchar(gene)
  )

mutation_categories <- c(
  "MutKnownFunctional",
  "MutLikelyFunctional",
  "MutNovel"
)
mutation_calls <- mut_df |>
  filter(category %in% mutation_categories) |>
  distinct(sample, gene)

sample_ids <- sort(unique(mut_df$sample))
gene_ids <- sort(unique(mutation_calls$gene))
mutation_mat <- matrix(
  0,
  nrow = length(gene_ids),
  ncol = length(sample_ids),
  dimnames = list(gene_ids, sample_ids)
)

if (nrow(mutation_calls) > 0) {
  mutation_mat[cbind(
    match(mutation_calls$gene, gene_ids),
    match(mutation_calls$sample, sample_ids)
  )] <- 1
}
storage.mode(mutation_mat) <- "double"

row_df <- data.frame(
  gene = gene_ids,
  stringsAsFactors = FALSE,
  row.names = gene_ids
)
col_df <- sample_annotation_from_ids(sample_ids)

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = mutation_mat),
  rowData = S4Vectors::DataFrame(row_df),
  colData = S4Vectors::DataFrame(col_df)
)
S4Vectors::metadata(se) <- list(
  datatype = "mutation",
  annotation = "mutation"
)

dir_create(path_dir(snk$output$se))
saveRDS(se, snk$output$se)
write_matrix_tsv(
  matrix = mutation_mat,
  row_ids = row_df$gene,
  row_id_name = "gene",
  output_path = snk$output$matrix
)
