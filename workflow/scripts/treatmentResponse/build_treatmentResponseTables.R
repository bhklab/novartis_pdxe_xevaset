snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fs)
})

experiment_df <- data.table::fread(snk$input$experiment) |>
  as.data.frame(stringsAsFactors = FALSE)
curve_metrics <- data.table::fread(snk$input$curveMetrics) |>
  as.data.frame(stringsAsFactors = FALSE)
model_df <- data.table::fread(snk$input$modelInfo) |>
  as.data.frame(stringsAsFactors = FALSE)
drug_df <- data.table::fread(snk$input$treatmentMetadata) |>
  as.data.frame(stringsAsFactors = FALSE)
control_drug <- snk$config[["treatmentResponse"]][["control_drug"]]

experiment_df <- experiment_df |>
  filter(
    model.id %in% model_df$model.id,
    drug %in% c(drug_df$drug.id, control_drug)
  ) |>
  arrange(model.id, time)

observed_treatments <- sort(intersect(
  unique(model_df$drug[model_df$drug != control_drug]),
  unique(experiment_df$drug[experiment_df$drug != control_drug])
))

dropped_treatments <- sort(setdiff(
  unique(drug_df$drug.id),
  c(observed_treatments, control_drug)
))
if (length(dropped_treatments) > 0) {
  message(
    "Dropping ",
    length(dropped_treatments),
    " treatment metadata row(s) not represented in both model and experiment tables: ",
    paste(dropped_treatments, collapse = ", ")
  )
}

experiment_df <- experiment_df |>
  filter(drug %in% c(observed_treatments, control_drug))

curve_metrics <- curve_metrics |>
  filter(treatment %in% observed_treatments) |>
  arrange(model, treatment)

dir_create(unique(path_dir(c(
  snk$output$experiment,
  snk$output$curveMetrics
))))
data.table::fwrite(experiment_df, snk$output$experiment, sep = "\t")
data.table::fwrite(curve_metrics, snk$output$curveMetrics, sep = "\t")
