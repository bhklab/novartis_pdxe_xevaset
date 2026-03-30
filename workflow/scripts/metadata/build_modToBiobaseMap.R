snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fs)
  library(SummarizedExperiment)
})

model_df <- data.table::fread(snk$input$modelInfo) |>
  as.data.frame(stringsAsFactors = FALSE)

extract_mapping <- function(se_path, datatype) {
  se <- readRDS(se_path)
  col_df <- as.data.frame(
    SummarizedExperiment::colData(se),
    stringsAsFactors = FALSE
  )
  col_df$colname <- colnames(se)
  names(col_df) <- clean_names_simple(names(col_df))

  patient_col <- first_present(
    c("patient_id", "patient.id", "sampleid", "sample_id", "biobase_id"),
    names(col_df)
  )

  if (is.na(patient_col)) {
    stop(
      "Could not determine patient/sample column for datatype ",
      datatype,
      call. = FALSE
    )
  }

  biobase_col <- first_present(
    c("biobase_id", "sampleid", "patient_id", "patient.id", "sample_id"),
    names(col_df)
  )
  biobase_ids <- if (is.na(biobase_col)) {
    col_df$colname
  } else {
    col_df[[biobase_col]]
  }

  tibble::tibble(
    patient.id = as.character(col_df[[patient_col]]),
    biobase.id = as.character(biobase_ids),
    mDataType = datatype
  )
}

mapping_df <- bind_rows(
  extract_mapping(snk$input$rnaseq, "RNASeq"),
  extract_mapping(snk$input$cnv, "cnv"),
  extract_mapping(snk$input$mutation, "mutation")
) |>
  inner_join(model_df |> distinct(model.id, patient.id), by = "patient.id") |>
  select(model.id, biobase.id, mDataType) |>
  distinct() |>
  arrange(mDataType, model.id, biobase.id)

dir_create(path_dir(snk$output$modToBiobaseMap))
data.table::fwrite(mapping_df, snk$output$modToBiobaseMap, sep = "\t")
