suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
})
snakemake@source("helpers.R")
snk <- parse_snakemake()

sample_df <- data.table::fread(snk$input$sampleMetadata) |>
  as.data.frame(stringsAsFactors = FALSE)
rownames(sample_df) <- sample_df$sampleid

assay_paths <- unname(unlist(snk$input$assays))
assays <- lapply(assay_paths, readRDS)
assay_names <- vapply(
  seq_along(assays),
  function(idx) {
    meta <- S4Vectors::metadata(assays[[idx]])
    if (!is.null(meta$datatype) && nzchar(meta$datatype)) {
      return(meta$datatype)
    }
    path_ext_remove(path_file(assay_paths[[idx]]))
  },
  character(1)
)
names(assays) <- make.unique(assay_names)

sample_maps <- lapply(names(assays), function(assay_name) {
  se <- assays[[assay_name]]
  col_df <- as.data.frame(
    SummarizedExperiment::colData(se),
    stringsAsFactors = FALSE
  )
  names(col_df) <- clean_names_simple(names(col_df))

  primary_col <- first_present(
    c("sampleid", "patient_id", "patient.id", "sample_id", "biobase_id"),
    names(col_df)
  )

  if (is.na(primary_col)) {
    stop(
      "No primary mapping column found for assay ",
      assay_name,
      call. = FALSE
    )
  }

  data.frame(
    assay = assay_name,
    primary = as.character(col_df[[primary_col]]),
    colname = colnames(se),
    stringsAsFactors = FALSE
  )
})
sample_map <- do.call(rbind, sample_maps)

missing_primary <- setdiff(unique(sample_map$primary), rownames(sample_df))
if (length(missing_primary)) {
  missing_df <- data.frame(
    sampleid = missing_primary,
    patient.id = missing_primary,
    tissue = NA_character_,
    tissue.name = NA_character_,
    stringsAsFactors = FALSE,
    row.names = missing_primary
  )
  for (col_name in setdiff(names(sample_df), names(missing_df))) {
    missing_df[[col_name]] <- NA
  }
  missing_df <- missing_df[, names(sample_df), drop = FALSE]
  sample_df <- rbind(sample_df, missing_df)
}

mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = assays,
  colData = S4Vectors::DataFrame(sample_df),
  sampleMap = sample_map,
  metadata = list(dataset = "PDXE")
)

dir_create(path_dir(snk$output$mae))
saveRDS(mae, snk$output$mae)
