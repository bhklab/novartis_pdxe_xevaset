suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fs)
})
snakemake@source("../helpers.R")
snk <- parse_snakemake()

tissue_name_map <- unlist(
  snk$config[["metadata"]][["tissue_name_map"]],
  use.names = TRUE
)

message("Loading normalized PCT tables")

pct_raw <- read_sheet_tsv(snk$input$pctRaw) |>
  normalize_character_columns()

curve_metrics <- read_sheet_tsv(snk$input$pctCurveMetrics) |>
  normalize_character_columns()

required_raw <- c(
  "model",
  "tumor_type",
  "treatment",
  "volume_mm3",
  "body_weight_g",
  "days_post_t0"
)
missing_raw <- setdiff(required_raw, colnames(pct_raw))
if (length(missing_raw)) {
  stop(
    "Missing required PCT raw columns: ",
    paste(missing_raw, collapse = ", "),
    call. = FALSE
  )
}

required_metrics <- c(
  "model",
  "treatment",
  "treatment_target",
  "treatment_type",
  "response_category"
)
missing_metrics <- setdiff(required_metrics, colnames(curve_metrics))
if (length(missing_metrics)) {
  stop(
    "Missing required PCT curve metric columns: ",
    paste(missing_metrics, collapse = ", "),
    call. = FALSE
  )
}

best_response_col <- first_present(
  c("best_response", "bestresponse"),
  colnames(curve_metrics)
)
day_best_response_col <- first_present(
  c("day_best_response", "day_bestresponse"),
  colnames(curve_metrics)
)
best_avg_response_col <- first_present(
  c("best_avg_response", "bestavgresponse"),
  colnames(curve_metrics)
)
day_best_avg_response_col <- first_present(
  c("day_best_avg_response", "day_bestavgresponse"),
  colnames(curve_metrics)
)
time_to_double_col <- first_present(
  c("time_to_double", "timetodouble"),
  colnames(curve_metrics)
)
day_last_col <- first_present(c("day_last"), colnames(curve_metrics))
response_category_col <- first_present(
  c("response_category", "responsecategory"),
  colnames(curve_metrics)
)

pct_raw <- pct_raw |>
  transmute(
    model = as.character(model),
    tumor_type = as.character(tumor_type),
    treatment = normalize_treatment_name(
      gsub('"', "", as.character(treatment), fixed = TRUE)
    ),
    volume_mm3 = as.numeric(volume_mm3),
    body_weight_g = as.numeric(body_weight_g),
    days_post_t0 = as.numeric(days_post_t0),
    tvol_difference = if ("tvol_difference" %in% colnames(pct_raw)) {
      as.numeric(tvol_difference)
    } else {
      NA_real_
    },
    bw_difference = if ("bw_difference" %in% colnames(pct_raw)) {
      as.numeric(bw_difference)
    } else {
      NA_real_
    }
  ) |>
  filter(!is.na(model), !is.na(treatment), !is.na(days_post_t0))

duplicate_timepoints <- pct_raw |>
  count(model, treatment, days_post_t0, name = "replicate_count") |>
  filter(replicate_count > 1)

if (nrow(duplicate_timepoints) > 0) {
  duplicate_groups <- duplicate_timepoints |>
    distinct(model, treatment)
  message(
    "Aggregating ",
    nrow(duplicate_timepoints),
    " duplicated patient/treatment/day measurement set(s) across ",
    nrow(duplicate_groups),
    " patient-treatment group(s) before timecourse splitting."
  )
}

experiment_split <- pct_raw |>
  group_split(model, treatment) |>
  lapply(split_timecourse) |>
  data.table::rbindlist(fill = TRUE) |>
  as.data.frame(stringsAsFactors = FALSE)

experiment_df <- experiment_split |>
  transmute(
    model.id = model.id,
    drug = treatment,
    time = days_post_t0,
    volume = volume_mm3,
    body.weight = body_weight_g
  )

model_df <- experiment_split |>
  distinct(model.id, tumor_type, model, treatment) |>
  transmute(
    model.id = model.id,
    tissue = tumor_type,
    tissue.name = dplyr::coalesce(
      unname(tissue_name_map[tumor_type]),
      tumor_type
    ),
    patient.id = model,
    drug = treatment
  ) |>
  filter(
    !is.na(model.id),
    !is.na(tissue),
    !is.na(tissue.name),
    !is.na(patient.id),
    !is.na(drug)
  )

curve_metrics <- curve_metrics |>
  transmute(
    model = as.character(model),
    treatment = normalize_treatment_name(as.character(treatment)),
    treatment.target = as.character(treatment_target),
    treatment.type = as.character(treatment_type),
    best.response = as.numeric(.data[[best_response_col]]),
    day.best.response = as.numeric(.data[[day_best_response_col]]),
    best.avg.response = as.numeric(.data[[best_avg_response_col]]),
    day.best.avg.response = as.numeric(.data[[day_best_avg_response_col]]),
    time.to.double = as.numeric(.data[[time_to_double_col]]),
    day.last = as.numeric(.data[[day_last_col]]),
    response.category = as.character(.data[[response_category_col]])
  )

drug_df <- curve_metrics |>
  distinct(treatment, treatment.target, treatment.type) |>
  transmute(
    drug.id = treatment,
    standard.name = treatment,
    treatment.target = treatment.target,
    treatment.type = treatment.type
  ) |>
  filter(!is.na(drug.id))

dir_create(unique(path_dir(c(
  snk$output$model,
  snk$output$drug,
  snk$output$experiment,
  snk$output$curveMetrics
))))
data.table::fwrite(model_df, snk$output$model, sep = "\t")
data.table::fwrite(drug_df, snk$output$drug, sep = "\t")
data.table::fwrite(experiment_df, snk$output$experiment, sep = "\t")
data.table::fwrite(curve_metrics, snk$output$curveMetrics, sep = "\t")
