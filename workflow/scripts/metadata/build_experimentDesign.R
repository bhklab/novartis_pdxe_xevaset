snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(data.table)
  library(fs)
})

model_df <- data.table::fread(snk$input$modelInfo) |>
  as.data.frame(stringsAsFactors = FALSE)

control_drug <- snk$config[["treatmentResponse"]][["control_drug"]]
drugs <- setdiff(unique(model_df$drug), control_drug)

exp_design <- list()
dropped_batches <- list()
for (patient_id in unique(model_df$patient.id)) {
  for (drug_id in drugs) {
    batch <- list(
      batch.name = paste(patient_id, drug_id, sep = "."),
      control = model_df$model.id[
        model_df$patient.id == patient_id & model_df$drug == control_drug
      ],
      treatment = model_df$model.id[
        model_df$patient.id == patient_id & model_df$drug == drug_id
      ]
    )

    if (length(batch$control) > 0 && length(batch$treatment) > 0) {
      exp_design[[batch$batch.name]] <- batch
    } else {
      dropped_batches[[length(dropped_batches) + 1L]] <- data.frame(
        batch.name = batch$batch.name,
        patient.id = patient_id,
        drug.id = drug_id,
        has_control = length(batch$control) > 0,
        has_treatment = length(batch$treatment) > 0,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(dropped_batches) > 0) {
  dropped_df <- do.call(rbind, dropped_batches)
  message(
    "Dropping ",
    nrow(dropped_df),
    " incomplete experiment-design batch(es)."
  )
  message(
    "  missing control: ",
    sum(!dropped_df$has_control)
  )
  message(
    "  missing treatment: ",
    sum(!dropped_df$has_treatment)
  )
  message(
    "  first dropped batches: ",
    paste(utils::head(dropped_df$batch.name, 10), collapse = ", ")
  )
}

dir_create(path_dir(snk$output$expDesign))
saveRDS(exp_design, snk$output$expDesign)
