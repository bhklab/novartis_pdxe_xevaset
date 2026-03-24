if (exists("snakemake")) {
  snakemake@source("helpers.R")
  snk <- parse_snakemake()
} else {
  source("workflow/scripts/helpers.R")

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 8) {
    stop(
      paste(
        "Usage:",
        "Rscript workflow/scripts/build_XevaSet.R",
        "<multi_assay_experiment>",
        "<model_info>",
        "<treatment_metadata>",
        "<experiment>",
        "<curve_metrics>",
        "<exp_design>",
        "<mod_to_biobase_map>",
        "<output_xeva>"
      ),
      call. = FALSE
    )
  }

  snk <- list(
    input = list(
      multiAssayExperiment = args[[1]],
      modelInfo = args[[2]],
      treatmentMetadata = args[[3]],
      experiment = args[[4]],
      curveMetrics = args[[5]],
      expDesign = args[[6]],
      modToBiobaseMap = args[[7]]
    ),
    output = list(xeva = args[[8]])
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(Biobase)
})

if (!requireNamespace("Xeva", quietly = TRUE)) {
  stop(
    "The Xeva package is not available in the current R library paths.",
    call. = FALSE
  )
}

model_df <- data.table::fread(snk$input$modelInfo) |>
  as.data.frame(stringsAsFactors = FALSE)
drug_df <- data.table::fread(snk$input$treatmentMetadata) |>
  as.data.frame(stringsAsFactors = FALSE)
experiment_df <- data.table::fread(snk$input$experiment) |>
  as.data.frame(stringsAsFactors = FALSE)
mod_map <- data.table::fread(snk$input$modToBiobaseMap) |>
  as.data.frame(stringsAsFactors = FALSE)
exp_design <- readRDS(snk$input$expDesign)
mae <- readRDS(snk$input$multiAssayExperiment)

control_drug <- "untreated"
if (
  !is.null(snk$config) &&
    !is.null(snk$config[["treatmentResponse"]]) &&
    !is.null(snk$config[["treatmentResponse"]][["control_drug"]])
) {
  control_drug <- as.character(snk$config[["treatmentResponse"]][["control_drug"]])
}

retain_drugs_with_data <- function(model_df, drug_df, experiment_df, control_drug) {
  retained_drugs <- sort(Reduce(intersect, list(
    unique(as.character(drug_df$drug.id)),
    unique(as.character(model_df$drug)),
    unique(as.character(experiment_df$drug))
  )))
  dropped_drugs <- sort(setdiff(
    unique(as.character(drug_df$drug.id)),
    c(retained_drugs, control_drug)
  ))

  if (length(dropped_drugs) > 0) {
    message(
      "Dropping ",
      length(dropped_drugs),
      " treatment metadata row(s) not represented in both model and experiment tables: ",
      paste(dropped_drugs, collapse = ", ")
    )
  }

  list(
    model = model_df[model_df$drug %in% retained_drugs, , drop = FALSE],
    drug = drug_df[drug_df$drug.id %in% retained_drugs, , drop = FALSE],
    experiment = experiment_df[experiment_df$drug %in% retained_drugs, , drop = FALSE]
  )
}


collect_warnings <- function(expr) {
  warnings_seen <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings_seen))
}


safe_batch_response <- function(
    object,
    batch_id,
    metric,
    min_time = 10,
    treatment_only = FALSE,
    max_time = NULL,
    concurrent_time = TRUE,
    vol_normal = FALSE,
    log_volume = FALSE
) {
  attempts <- list(
    list(label = "default", concurrent.time = concurrent_time, impute.value = TRUE),
    list(label = "no_impute", concurrent.time = concurrent_time, impute.value = FALSE)
  )

  last_error <- NA_character_
  last_warnings <- character()

  for (attempt in attempts) {
    response_attempt <- collect_warnings(
      tryCatch(
        Xeva::response(
          object,
          batch = batch_id,
          res.measure = metric,
          treatment.only = treatment_only,
          max.time = max_time,
          impute.value = attempt$impute.value,
          min.time = min_time,
          concurrent.time = attempt$concurrent.time,
          vol.normal = vol_normal,
          log.volume = log_volume,
          verbose = FALSE
        ),
        error = function(err) err
      )
    )

    last_warnings <- response_attempt$warnings
    if (!inherits(response_attempt$value, "error")) {
      return(list(
        ok = TRUE,
        result = response_attempt$value,
        method = attempt$label,
        warnings = response_attempt$warnings
      ))
    }

    last_error <- conditionMessage(response_attempt$value)
  }

  list(
    ok = FALSE,
    result = NULL,
    method = NA_character_,
    warnings = last_warnings,
    error = last_error
  )
}


set_batch_metric <- function(
    object,
    metric,
    min_time = 10,
    treatment_only = FALSE,
    max_time = NULL,
    concurrent_time = TRUE,
    vol_normal = FALSE,
    log_volume = FALSE
) {
  sensitivity <- methods::slot(object, "sensitivity")
  batch_ids <- Xeva::batchInfo(object)

  metric_cols <- switch(
    metric,
    angle = c("slope.control", "slope.treatment", "angle"),
    abc = c("auc.control", "auc.treatment", "abc"),
    TGI = c("TGI"),
    stop("Unsupported batch metric: ", metric, call. = FALSE)
  )

  sensitivity$batch[, metric_cols] <- NA_real_

  fallback_batches <- character()
  failed_batches <- character()

  for (batch_id in batch_ids) {
    batch_result <- safe_batch_response(
      object = object,
      batch_id = batch_id,
      metric = metric,
      min_time = min_time,
      treatment_only = treatment_only,
      max_time = max_time,
      concurrent_time = concurrent_time,
      vol_normal = vol_normal,
      log_volume = log_volume
    )

    if (!batch_result$ok) {
      failed_batches <- c(
        failed_batches,
        sprintf("%s (%s)", batch_id, batch_result$error)
      )
      next
    }

    if (identical(batch_result$method, "no_impute")) {
      fallback_batches <- c(fallback_batches, batch_id)
    }

    if (metric == "angle") {
      metric_values <- c(
        batch_result$result$control$value,
        batch_result$result$treatment$value,
        batch_result$result$value
      )
      if (length(metric_values) != length(metric_cols)) {
        failed_batches <- c(
          failed_batches,
          sprintf("%s (empty %s value)", batch_id, metric)
        )
        next
      }
      sensitivity$batch[batch_id, metric_cols] <- as.numeric(metric_values)
    } else if (metric == "abc") {
      metric_values <- c(
        batch_result$result$control$value,
        batch_result$result$treatment$value,
        batch_result$result$value
      )
      if (length(metric_values) != length(metric_cols)) {
        failed_batches <- c(
          failed_batches,
          sprintf("%s (empty %s value)", batch_id, metric)
        )
        next
      }
      sensitivity$batch[batch_id, metric_cols] <- as.numeric(metric_values)
    } else if (metric == "TGI") {
      metric_value <- batch_result$result$value
      if (length(metric_value) == 0) {
        failed_batches <- c(
          failed_batches,
          sprintf("%s (empty %s value)", batch_id, metric)
        )
        next
      }
      sensitivity$batch[batch_id, metric_cols] <- as.numeric(metric_value[[1]])
    }
  }

  if (length(fallback_batches) > 0) {
    message(
      "Recovered ",
      metric,
      " for ",
      length(fallback_batches),
      " batch(es) by disabling interpolation after Xeva::response() failed with imputation enabled."
    )
  }

  if (length(failed_batches) > 0) {
    warning(
      paste0(
        "Unable to compute ",
        metric,
        " for ",
        length(failed_batches),
        " batch(es). First failures: ",
        paste(utils::head(failed_batches, 10), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  methods::slot(object, "sensitivity") <- sensitivity
  object
}


create_xeva <- function(molecular_profiles) {
  Xeva::createXevaSet(
    name = "PDXE xevaSet",
    model = model_df,
    drug = drug_df,
    experiment = experiment_df,
    expDesign = exp_design,
    molecularProfiles = molecular_profiles,
    modToBiobaseMap = mod_map
  )
}

filtered_inputs <- retain_drugs_with_data(
  model_df,
  drug_df,
  experiment_df,
  control_drug = control_drug
)
model_df <- filtered_inputs$model
drug_df <- filtered_inputs$drug
experiment_df <- filtered_inputs$experiment

tissue_missing <- is.na(model_df$tissue) | !nzchar(trimws(as.character(model_df$tissue)))
tissue_name_missing <- is.na(model_df$tissue.name) | !nzchar(trimws(as.character(model_df$tissue.name)))
oncotree_missing <- rep(FALSE, nrow(model_df))
if ("oncotree.mainType" %in% colnames(model_df)) {
  oncotree_missing <- is.na(model_df$oncotree.mainType) |
    !nzchar(trimws(as.character(model_df$oncotree.mainType)))
}

missing_annotation_rows <- model_df[
  tissue_missing | tissue_name_missing | oncotree_missing,
  ,
  drop = FALSE
]

if (nrow(missing_annotation_rows) > 0) {
  stop(
    paste0(
      "Model metadata contains missing tissue annotations for model.id: ",
      paste(utils::head(missing_annotation_rows$model.id, 10), collapse = ", ")
    ),
    call. = FALSE
  )
}

pdxe <- tryCatch(
  create_xeva(mae),
  error = function(err) {
    message(
      "Falling back to ExpressionSet list for Xeva::createXevaSet(): ",
      conditionMessage(err)
    )
    exp_list <- lapply(
      MultiAssayExperiment::experiments(mae),
      se_to_expression_set
    )
    create_xeva(exp_list)
  }
)

for (metric in c("mRECIST", "slope", "AUC", "angle", "abc", "TGI")) {
  message("Computing response metric: ", metric)
  if (metric %in% c("angle", "abc", "TGI")) {
    pdxe <- set_batch_metric(pdxe, metric = metric)
  } else {
    pdxe <- tryCatch(
      Xeva::setResponse(pdxe, res.measure = metric, verbose = TRUE),
      error = function(err) {
        warning(
          paste0(
            "Skipping response metric ",
            metric,
            ": ",
            conditionMessage(err)
          ),
          call. = FALSE
        )
        pdxe
      }
    )
  }
}

required_batch_cols <- c(
  "batch.name",
  "slope.control",
  "slope.treatment",
  "angle",
  "auc.control",
  "auc.treatment",
  "abc",
  "TGI"
)
missing_batch_cols <- setdiff(
  required_batch_cols,
  colnames(methods::slot(pdxe, "sensitivity")$batch)
)
if (length(missing_batch_cols) > 0) {
  stop(
    "The rebuilt XevaSet is missing expected batch sensitivity columns: ",
    paste(missing_batch_cols, collapse = ", "),
    call. = FALSE
  )
}

dir_create(path_dir(snk$output$xeva))
saveRDS(pdxe, snk$output$xeva)
