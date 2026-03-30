suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(readxl)
})
snakemake@source("../helpers.R")
snk <- parse_snakemake()

sheet_cfg <- snk$config[["metadata"]][["supplementary"]][["sheets"]]

read_sheet <- function(sheet_name) {
  readxl::read_xlsx(snk$input$workbook, sheet = sheet_name) |>
    as.data.frame(stringsAsFactors = FALSE, check.names = FALSE) |>
    normalize_character_columns()
}

normalize_wide_sheet <- function(df, first_col_name) {
  names(df)[[1]] <- first_col_name
  df[[first_col_name]] <- trimws(as.character(df[[first_col_name]]))
  df
}

normalize_treatment_sheet <- function(df) {
  names(df) <- clean_names_simple(names(df))
  df <- normalize_character_columns(df)
  if ("treatment" %in% names(df)) {
    df$treatment <- gsub('"', "", df$treatment, fixed = TRUE)
  }
  df
}

rnaseq_df <- read_sheet(sheet_cfg[["rnaseq"]][["name"]]) |>
  normalize_wide_sheet("gene")
cnv_df <- read_sheet(sheet_cfg[["cnv"]][["name"]]) |>
  normalize_wide_sheet("feature")
mutation_df <- read_sheet(sheet_cfg[["mutation"]][["name"]])
names(mutation_df) <- clean_names_simple(names(mutation_df))
mutation_df <- normalize_character_columns(mutation_df)
pct_raw_df <- read_sheet(sheet_cfg[["pct_raw"]][["name"]]) |>
  normalize_treatment_sheet()
pct_curve_metrics_df <- read_sheet(
  sheet_cfg[["pct_curve_metrics"]][["name"]]
) |>
  normalize_treatment_sheet()

dir_create(unique(path_dir(unlist(snk$output, use.names = FALSE))))
data.table::fwrite(rnaseq_df, snk$output$rnaseq, sep = "\t")
data.table::fwrite(cnv_df, snk$output$cnv, sep = "\t")
data.table::fwrite(mutation_df, snk$output$mutation, sep = "\t")
data.table::fwrite(pct_raw_df, snk$output$pctRaw, sep = "\t")
data.table::fwrite(
  pct_curve_metrics_df,
  snk$output$pctCurveMetrics,
  sep = "\t"
)
