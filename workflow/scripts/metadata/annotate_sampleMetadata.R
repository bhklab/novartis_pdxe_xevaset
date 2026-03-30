snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(AnnotationGx)
  library(data.table)
  library(dplyr)
  library(fs)
})

default_oncotree_code_map <- c(
  "BRCA" = "BREAST",
  "CRC" = "BOWEL",
  "CM" = "MEL",
  "GC" = "STOMACH",
  "NSCLC" = "LUNG",
  "PDAC" = "PANCREAS"
)

resolve_oncotree_code_map <- function(config_map = NULL) {
  if (is.null(config_map)) {
    return(default_oncotree_code_map)
  }

  if (
    is.null(names(config_map)) ||
      length(names(config_map)) != length(config_map) ||
      any(is.na(names(config_map))) ||
      !all(nzchar(names(config_map)))
  ) {
    return(default_oncotree_code_map)
  }

  config_names <- names(config_map)
  config_map <- normalize_missing_character(unname(config_map))
  names(config_map) <- config_names
  config_map <- config_map[!is.na(config_map)]

  if (!length(config_map)) {
    return(default_oncotree_code_map)
  }

  config_map
}

extract_first_nested_value <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }

  if (is.list(x)) {
    x <- unlist(x, recursive = TRUE, use.names = FALSE)
  }

  x <- normalize_missing_character(x)
  x <- x[!is.na(x)]

  if (!length(x)) {
    return(NA_character_)
  }

  x[[1]]
}

build_oncotree_annotation_tbl <- function(project_to_oncotree_code_map) {
  required_codes <- sort(unique(unname(project_to_oncotree_code_map)))

  oncotree_raw <- tryCatch(
    as.data.frame(
      AnnotationGx::getOncotreeTumorTypes(),
      stringsAsFactors = FALSE
    ),
    error = function(err) {
      stop(
        paste0(
          "AnnotationGx::getOncotreeTumorTypes failed: ",
          conditionMessage(err)
        ),
        call. = FALSE
      )
    }
  )

  missing_codes <- setdiff(required_codes, oncotree_raw$code)
  if (length(missing_codes) > 0) {
    stop(
      paste0(
        "OncoTree annotations were not found for configured tissue code(s): ",
        paste(missing_codes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  oncotree_lookup <- oncotree_raw |>
    filter(code %in% required_codes) |>
    transmute(
      oncotree.code = normalize_missing_character(code),
      oncotree.name = normalize_missing_character(name),
      oncotree.mainType = normalize_missing_character(mainType),
      oncotree.umlsCode = vapply(
        externalReferences,
        extract_first_nested_value,
        character(1)
      ),
      oncotree.nciTissueCode = vapply(
        tissue,
        extract_first_nested_value,
        character(1)
      )
    ) |>
    distinct(oncotree.code, .keep_all = TRUE)

  annotation_tbl <- data.frame(
    tissue = names(project_to_oncotree_code_map),
    oncotree.code = unname(project_to_oncotree_code_map),
    stringsAsFactors = FALSE
  ) |>
    left_join(oncotree_lookup, by = "oncotree.code")

  missing_annotation_rows <- annotation_tbl |>
    filter(
      is.na(oncotree.code) |
        is.na(oncotree.name) |
        is.na(oncotree.mainType) |
        is.na(oncotree.umlsCode) |
        is.na(oncotree.nciTissueCode)
    )

  if (nrow(missing_annotation_rows) > 0) {
    stop(
      paste0(
        "Configured tissue code(s) are missing required OncoTree annotations: ",
        paste(missing_annotation_rows$tissue, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  annotation_tbl
}

first_non_missing <- function(x) {
  x <- normalize_missing_character(x)
  x <- x[!is.na(x)]
  if (length(x)) {
    return(x[[1]])
  }
  NA_character_
}

coalesce_missing_character <- function(x, fallback) {
  x <- normalize_missing_character(x)
  fallback <- normalize_missing_character(fallback)
  dplyr::coalesce(x, fallback)
}

config_oncotree_code_map <- NULL
if (
  !is.null(snk$config) &&
    !is.null(snk$config[["metadata"]]) &&
    !is.null(snk$config[["metadata"]][["oncotree_tissue_code_map"]])
) {
  config_oncotree_code_map <- unlist(
    snk$config[["metadata"]][["oncotree_tissue_code_map"]],
    use.names = TRUE
  )
}

project_to_oncotree_code_map <- resolve_oncotree_code_map(
  config_map = config_oncotree_code_map
)
message("Retrieving OncoTree branch annotations from AnnotationGx")
oncotree_annotation_tbl <- build_oncotree_annotation_tbl(
  project_to_oncotree_code_map = project_to_oncotree_code_map
)

model_df <- data.table::fread(snk$input$model) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  distinct(model.id, .keep_all = TRUE) |>
  mutate(
    tissue = normalize_missing_character(tissue),
    tissue.name = normalize_missing_character(tissue.name)
  )

patient_backfill <- model_df |>
  group_by(patient.id) |>
  summarise(
    tissue.patient = first_non_missing(tissue),
    tissue.name.patient = first_non_missing(tissue.name),
    .groups = "drop"
  )

model_df <- model_df |>
  left_join(patient_backfill, by = "patient.id") |>
  mutate(
    tissue = coalesce_missing_character(tissue, tissue.patient),
    tissue.name = coalesce_missing_character(tissue.name, tissue.name.patient)
  ) |>
  left_join(oncotree_annotation_tbl, by = "tissue") |>
  select(-tissue.patient, -tissue.name.patient) |>
  arrange(model.id)

invalid_rows <- model_df |>
  filter(
    is.na(tissue) |
      is.na(tissue.name) |
      is.na(oncotree.code) |
      is.na(oncotree.name) |
      is.na(oncotree.mainType) |
      is.na(oncotree.umlsCode) |
      is.na(oncotree.nciTissueCode)
  )

if (nrow(invalid_rows) > 0) {
  stop(
    paste0(
      "Model metadata still has missing tissue annotations after patient-level backfill for model.id: ",
      paste(utils::head(invalid_rows$model.id, 10), collapse = ", ")
    ),
    call. = FALSE
  )
}

sample_df <- model_df |>
  transmute(
    sampleid = patient.id,
    patient.id = patient.id,
    tissue = tissue,
    tissue.name = tissue.name,
    oncotree.code = oncotree.code,
    oncotree.name = oncotree.name,
    oncotree.mainType = oncotree.mainType,
    oncotree.umlsCode = oncotree.umlsCode,
    oncotree.nciTissueCode = oncotree.nciTissueCode
  ) |>
  group_by(sampleid, patient.id) |>
  summarise(
    tissue = first_non_missing(tissue),
    tissue.name = first_non_missing(tissue.name),
    oncotree.code = first_non_missing(oncotree.code),
    oncotree.name = first_non_missing(oncotree.name),
    oncotree.mainType = first_non_missing(oncotree.mainType),
    oncotree.umlsCode = first_non_missing(oncotree.umlsCode),
    oncotree.nciTissueCode = first_non_missing(oncotree.nciTissueCode),
    .groups = "drop"
  ) |>
  arrange(sampleid)

dir_create(unique(path_dir(c(
  snk$output$modelInfo,
  snk$output$sampleMetadata
))))
data.table::fwrite(model_df, snk$output$modelInfo, sep = "\t")
data.table::fwrite(sample_df, snk$output$sampleMetadata, sep = "\t")
