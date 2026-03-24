snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(fs)
  library(SummarizedExperiment)
})

message("Loading RNASeq matrix: ", snk$input$raw)
rna <- wide_assay_tsv_to_se(
  path = snk$input$raw,
  datatype = "RNASeq",
  feature_col = "gene",
  row_data_name = "gene"
)
se <- rna$se

dir_create(path_dir(snk$output$se))
saveRDS(se, snk$output$se)
write_matrix_tsv(
  matrix = rna$matrix,
  row_ids = rna$row_data$gene,
  row_id_name = "gene",
  output_path = snk$output$matrix
)
