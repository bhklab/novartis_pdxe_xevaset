if (exists("snakemake")) {
  snakemake@source("helpers.R")
  snk <- parse_snakemake()
} else {
  source("workflow/scripts/helpers.R")

  args <- commandArgs(trailingOnly = TRUE)
  if (!(length(args) %in% c(8, 9))) {
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
        "<output_xeva>",
        "[output_tsv_dir]"
      ),
      call. = FALSE
    )
  }

  default_tsv_dir <- fs::path(
    path_dir(args[[8]]),
    paste0(path_ext_remove(path_file(args[[8]])), "_tsv")
  )

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
    output = list(
      xeva = args[[8]],
      tsvDir = if (length(args) >= 9) args[[9]] else default_tsv_dir
    )
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

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop(
    "The jsonlite package is not available in the current R library paths.",
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
  control_drug <- as.character(snk$config[["treatmentResponse"]][[
    "control_drug"
  ]])
}

retain_drugs_with_data <- function(
  model_df,
  drug_df,
  experiment_df,
  control_drug
) {
  retained_drugs <- sort(Reduce(
    intersect,
    list(
      unique(as.character(drug_df$drug.id)),
      unique(as.character(model_df$drug)),
      unique(as.character(experiment_df$drug))
    )
  ))
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
    experiment = experiment_df[
      experiment_df$drug %in% retained_drugs,
      ,
      drop = FALSE
    ]
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
    list(
      label = "default",
      concurrent.time = concurrent_time,
      impute.value = TRUE
    ),
    list(
      label = "no_impute",
      concurrent.time = concurrent_time,
      impute.value = FALSE
    )
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


serialize_tsv_value <- function(x) {
  if (is.null(x) || !length(x)) {
    return(NA_character_)
  }

  if (inherits(x, "sessionInfo")) {
    x <- paste(capture.output(print(x)), collapse = " | ")
  }

  if (inherits(x, "POSIXt")) {
    return(format(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  }

  if (inherits(x, "Date")) {
    return(as.character(x))
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.data.frame(x)) {
    x <- lapply(x, function(col) {
      if (is.factor(col)) {
        return(as.character(col))
      }
      col
    })
  }

  value <- if (is.list(x) || length(x) > 1) {
    jsonlite::toJSON(
      x,
      auto_unbox = TRUE,
      null = "null",
      dataframe = "rows",
      na = "null"
    )
  } else {
    as.character(x[[1]])
  }

  gsub("[\r\n\t]+", " ", value)
}


normalize_for_tsv <- function(df) {
  df[] <- lapply(
    df,
    function(col) {
      if (is.factor(col)) {
        return(as.character(col))
      }

      if (inherits(col, "POSIXt")) {
        return(format(col, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
      }

      if (inherits(col, "Date")) {
        return(as.character(col))
      }

      if (is.list(col)) {
        return(vapply(col, serialize_tsv_value, character(1)))
      }

      col
    }
  )
  df
}


write_tsv_export <- function(df, output_path, rownames_col = NULL) {
  out_df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)

  if (!is.null(rownames_col)) {
    out_df <- data.frame(
      stats::setNames(
        data.frame(rownames = rownames(out_df), stringsAsFactors = FALSE),
        rownames_col
      ),
      out_df,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out_df <- normalize_for_tsv(out_df)
  dir_create(path_dir(output_path))
  data.table::fwrite(out_df, output_path, sep = "\t")
  out_df
}


batch_definitions_to_df <- function(exp_design) {
  data.table::rbindlist(
    lapply(
      exp_design,
      function(batch) {
        as.data.frame(as.list(batch), stringsAsFactors = FALSE)
      }
    ),
    fill = TRUE
  ) |>
    as.data.frame(stringsAsFactors = FALSE)
}


build_experiment_exports <- function(experiment_list) {
  measurement_rows <- lapply(
    experiment_list,
    function(model_obj) {
      measurement_df <- as.data.frame(model_obj@data, stringsAsFactors = FALSE)
      measurement_df$model.id <- model_obj@model.id
      measurement_df$measurement_index <- seq_len(nrow(measurement_df))
      measurement_df[,
        c(
          "model.id",
          "measurement_index",
          setdiff(colnames(measurement_df), c("model.id", "measurement_index"))
        ),
        drop = FALSE
      ]
    }
  )

  metadata_rows <- lapply(
    experiment_list,
    function(model_obj) {
      row <- list()

      for (slot_name in slotNames(model_obj)) {
        if (identical(slot_name, "data")) {
          next
        }

        slot_value <- methods::slot(model_obj, slot_name)
        if (
          is.list(slot_value) &&
            !is.data.frame(slot_value) &&
            length(slot_value) > 0 &&
            !is.null(names(slot_value)) &&
            all(nzchar(names(slot_value)))
        ) {
          for (value_name in names(slot_value)) {
            row[[paste(
              slot_name,
              value_name,
              sep = "."
            )]] <- serialize_tsv_value(
              slot_value[[value_name]]
            )
          }
        } else {
          row[[slot_name]] <- serialize_tsv_value(slot_value)
        }
      }

      as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
    }
  )

  list(
    metadata = data.table::rbindlist(metadata_rows, fill = TRUE) |>
      as.data.frame(stringsAsFactors = FALSE),
    measurements = data.table::rbindlist(measurement_rows, fill = TRUE) |>
      as.data.frame(stringsAsFactors = FALSE)
  )
}


annotation_to_df <- function(annotation) {
  do.call(
    rbind,
    lapply(
      names(annotation),
      function(field_name) {
        data.frame(
          field = field_name,
          value = serialize_tsv_value(annotation[[field_name]]),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}


export_xeva_tsv <- function(pdxe, output_dir) {
  if (dir_exists(output_dir)) {
    dir_delete(output_dir)
  }
  dir_create(output_dir)

  manifest_entries <- list()
  add_manifest_entry <- function(
    component,
    export_path,
    description,
    exported_df
  ) {
    manifest_entries[[length(manifest_entries) + 1]] <<- data.frame(
      component = component,
      relative_path = fs::path_rel(export_path, start = output_dir),
      rows = nrow(exported_df),
      columns = ncol(exported_df),
      description = description,
      stringsAsFactors = FALSE
    )
  }

  export_component <- function(
    component,
    relative_path,
    description,
    df,
    rownames_col = NULL
  ) {
    export_path <- fs::path(output_dir, relative_path)
    exported_df <- write_tsv_export(
      df = df,
      output_path = export_path,
      rownames_col = rownames_col
    )
    add_manifest_entry(component, export_path, description, exported_df)
  }

  export_component(
    component = "annotation",
    relative_path = "annotation.tsv",
    description = "XevaSet annotation metadata, including creation info.",
    df = annotation_to_df(methods::slot(pdxe, "annotation"))
  )
  export_component(
    component = "model",
    relative_path = "model.tsv",
    description = "Model metadata from slot(pdxe, 'model').",
    df = methods::slot(pdxe, "model")
  )
  export_component(
    component = "drug",
    relative_path = "drug.tsv",
    description = "Treatment metadata from slot(pdxe, 'drug').",
    df = methods::slot(pdxe, "drug")
  )
  export_component(
    component = "sensitivity_model",
    relative_path = "sensitivity_model.tsv",
    description = "Model-level response summaries from slot(pdxe, 'sensitivity')$model.",
    df = methods::slot(pdxe, "sensitivity")$model
  )
  export_component(
    component = "sensitivity_batch",
    relative_path = "sensitivity_batch.tsv",
    description = "Batch-level response summaries from slot(pdxe, 'sensitivity')$batch.",
    df = methods::slot(pdxe, "sensitivity")$batch
  )
  export_component(
    component = "batch_definitions",
    relative_path = "batch_definitions.tsv",
    description = "Treatment-control batch definitions from slot(pdxe, 'expDesign').",
    df = batch_definitions_to_df(methods::slot(pdxe, "expDesign"))
  )
  export_component(
    component = "mod_to_biobase_map",
    relative_path = "mod_to_biobase_map.tsv",
    description = "Mapping between Xeva molecular profile identifiers and Biobase annotations.",
    df = methods::slot(pdxe, "modToBiobaseMap")
  )

  experiment_exports <- build_experiment_exports(methods::slot(
    pdxe,
    "experiment"
  ))
  export_component(
    component = "experiment_model_metadata",
    relative_path = "experiment_model_metadata.tsv",
    description = "Flattened per-model metadata from slot(pdxe, 'experiment').",
    df = experiment_exports$metadata
  )
  export_component(
    component = "experiment_measurements",
    relative_path = "experiment_measurements.tsv",
    description = "Long-form tumor-volume time-course data from slot(pdxe, 'experiment').",
    df = experiment_exports$measurements
  )

  molecular_profiles <- methods::slot(pdxe, "molecularProfiles")
  for (profile_name in names(molecular_profiles)) {
    profile <- molecular_profiles[[profile_name]]

    export_component(
      component = paste0("molecular_profile_assay:", profile_name),
      relative_path = fs::path(
        "molecular_profiles",
        paste0(profile_name, ".tsv")
      ),
      description = paste0(
        "Expression matrix for molecular profile '",
        profile_name,
        "'."
      ),
      df = data.frame(
        as.data.frame(Biobase::exprs(profile), check.names = FALSE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      rownames_col = "feature_rowname"
    )
  }

  write_tsv_export(
    df = do.call(rbind, manifest_entries),
    output_path = fs::path(output_dir, "file_manifest.tsv")
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

tissue_missing <- is.na(model_df$tissue) |
  !nzchar(trimws(as.character(model_df$tissue)))
tissue_name_missing <- is.na(model_df$tissue.name) |
  !nzchar(trimws(as.character(model_df$tissue.name)))
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

if (!is.null(snk$output$tsvDir) && nzchar(snk$output$tsvDir)) {
  export_xeva_tsv(pdxe, snk$output$tsvDir)
}
