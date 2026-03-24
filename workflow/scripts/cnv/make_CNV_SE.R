snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(fs)
  library(SummarizedExperiment)
})

message("Loading CNV matrix: ", snk$input$raw)
cnv <- wide_assay_tsv_to_se(
  path = snk$input$raw,
  datatype = "cnv",
  feature_col = "feature",
  row_data_name = "feature"
)
se <- cnv$se

dir_create(path_dir(snk$output$se))
saveRDS(se, snk$output$se)
write_matrix_tsv(
  matrix = cnv$matrix,
  row_ids = cnv$row_data$feature,
  row_id_name = "feature",
  output_path = snk$output$matrix
)
