suppressPackageStartupMessages({
  library(fs)
})

parse_snakemake <- function() {
  if (!exists("snakemake")) {
    stop("This script is intended to be executed by Snakemake.", call. = FALSE)
  }

  if (length(snakemake@log) > 0) {
    log_file <- snakemake@log[[1]]
    dir_create(path_dir(log_file))
    sink(log_file, append = FALSE, split = TRUE)
    message_connection <- file(log_file, open = "at")
    sink(message_connection, append = TRUE, type = "message")
  }

  list(
    input = snakemake@input,
    output = snakemake@output,
    params = snakemake@params,
    threads = snakemake@threads,
    config = snakemake@config
  )
}


clean_names_simple <- function(x) {
  x <- gsub("%", "", x, fixed = TRUE)
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("(^_+|_+$)", "", x)
  tolower(x)
}


normalize_missing_character <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_
  x
}


normalize_treatment_name <- function(x) {
  x <- normalize_missing_character(x)
  x[x == "5FFU"] <- "5FU"
  x
}


first_non_missing_value <- function(x) {
  if (is.character(x)) {
    x <- normalize_missing_character(x)
  }
  idx <- which(!is.na(x))
  if (!length(idx)) {
    return(NA)
  }
  x[[idx[[1]]]]
}


mean_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}


aggregate_timecourse_by_day <- function(df, day_col = "days_post_t0") {
  if (!nrow(df)) {
    return(df)
  }

  day_values <- sort(unique(df[[day_col]]))
  rows <- lapply(
    day_values,
    function(day_value) {
      idx <- which(df[[day_col]] == day_value)
      row <- df[idx[[1]], , drop = FALSE]

      for (col_name in colnames(df)) {
        values <- df[[col_name]][idx]
        if (col_name == day_col) {
          row[[col_name]] <- day_value
        } else if (is.numeric(df[[col_name]])) {
          row[[col_name]] <- mean_or_na(values)
        } else {
          row[[col_name]] <- first_non_missing_value(values)
        }
      }

      row
    }
  )

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


abbreviate_treatment <- function(drug_name) {
  pieces <- strsplit(drug_name, "\\+", fixed = FALSE)[[1]]
  pieces <- trimws(pieces)
  codes <- vapply(
    pieces,
    function(piece) {
      letters <- strsplit(piece, "", fixed = TRUE)[[1]]
      if (length(letters) <= 4) {
        return(paste0(letters, collapse = ""))
      }
      paste0(
        letters[1],
        letters[2],
        letters[length(letters) - 1],
        letters[length(letters)]
      )
    },
    character(1)
  )
  paste(codes, collapse = ".")
}


split_timecourse <- function(df) {
  df <- df[order(df$days_post_t0), , drop = FALSE]
  if (anyDuplicated(df$days_post_t0)) {
    df <- aggregate_timecourse_by_day(df)
  }
  base_model_id <- paste0(
    gsub("-", ".", df$model[[1]], fixed = TRUE),
    ".",
    abbreviate_treatment(df$treatment[[1]])
  )

  if (nrow(df) <= 1) {
    df$model.id <- base_model_id
    return(df)
  }

  segment <- cumsum(c(0, diff(df$days_post_t0) < 0))
  if (max(segment) == 0) {
    df$model.id <- base_model_id
  } else {
    df$model.id <- paste0(base_model_id, ".", segment)
  }
  df
}


expression_set_to_se <- function(eset, datatype) {
  assay_mat <- Biobase::exprs(eset)

  row_df <- as.data.frame(Biobase::fData(eset), stringsAsFactors = FALSE)
  if (!nrow(row_df)) {
    row_df <- data.frame(
      feature_id = rownames(assay_mat),
      stringsAsFactors = FALSE
    )
  }
  if (is.null(rownames(row_df)) || any(rownames(row_df) == "")) {
    rownames(row_df) <- rownames(assay_mat)
  }

  col_df <- as.data.frame(Biobase::pData(eset), stringsAsFactors = FALSE)
  if (!nrow(col_df)) {
    col_df <- data.frame(
      sampleid = colnames(assay_mat),
      stringsAsFactors = FALSE
    )
  }
  if (is.null(rownames(col_df)) || any(rownames(col_df) == "")) {
    rownames(col_df) <- colnames(assay_mat)
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = as.matrix(assay_mat)),
    rowData = S4Vectors::DataFrame(row_df),
    colData = S4Vectors::DataFrame(col_df)
  )
  S4Vectors::metadata(se) <- list(
    datatype = datatype,
    annotation = datatype
  )
  se
}


se_to_expression_set <- function(se) {
  assay_mat <- SummarizedExperiment::assay(se)
  row_df <- as.data.frame(
    SummarizedExperiment::rowData(se),
    stringsAsFactors = FALSE
  )
  col_df <- as.data.frame(
    SummarizedExperiment::colData(se),
    stringsAsFactors = FALSE
  )

  if (!nrow(row_df)) {
    row_df <- data.frame(
      feature_id = rownames(assay_mat),
      stringsAsFactors = FALSE
    )
  }
  if (!nrow(col_df)) {
    col_df <- data.frame(
      sampleid = colnames(assay_mat),
      stringsAsFactors = FALSE
    )
  }

  rownames(row_df) <- rownames(assay_mat)
  rownames(col_df) <- colnames(assay_mat)

  Biobase::ExpressionSet(
    assayData = as.matrix(assay_mat),
    phenoData = Biobase::AnnotatedDataFrame(col_df),
    featureData = Biobase::AnnotatedDataFrame(row_df)
  )
}


first_present <- function(candidates, choices) {
  match <- candidates[candidates %in% choices]
  if (!length(match)) {
    return(NA_character_)
  }
  match[[1]]
}


read_sheet_tsv <- function(path) {
  data.table::fread(
    path,
    sep = "\t",
    data.table = FALSE,
    check.names = FALSE
  )
}


normalize_character_columns <- function(df) {
  df[] <- lapply(
    df,
    function(col) {
      if (is.character(col)) {
        return(trimws(col))
      }
      col
    }
  )
  df
}


sample_annotation_from_ids <- function(sample_ids) {
  sample_ids <- as.character(sample_ids)
  data.frame(
    sampleid = sample_ids,
    patient.id = sample_ids,
    biobase.id = sample_ids,
    stringsAsFactors = FALSE,
    row.names = sample_ids
  )
}


wide_assay_tsv_to_se <- function(
  path,
  datatype,
  feature_col,
  row_data_name = feature_col,
  assay_name = "exprs"
) {
  df <- read_sheet_tsv(path)

  if (!(feature_col %in% names(df))) {
    stop(
      "Missing feature column '",
      feature_col,
      "' in ",
      path,
      call. = FALSE
    )
  }

  feature_ids <- trimws(as.character(df[[feature_col]]))
  keep_rows <- !is.na(feature_ids) & nzchar(feature_ids)
  df <- df[keep_rows, , drop = FALSE]
  feature_ids <- feature_ids[keep_rows]

  sample_ids <- setdiff(names(df), feature_col)
  if (!length(sample_ids)) {
    stop("No sample columns found in ", path, call. = FALSE)
  }

  assay_df <- df[, sample_ids, drop = FALSE]
  assay_df[] <- lapply(
    assay_df,
    function(col) as.numeric(as.character(col))
  )
  assay_mat <- as.matrix(assay_df)
  storage.mode(assay_mat) <- "double"
  rownames(assay_mat) <- make.unique(feature_ids, sep = "__")
  colnames(assay_mat) <- sample_ids

  row_df <- data.frame(
    stringsAsFactors = FALSE,
    row.names = rownames(assay_mat)
  )
  row_df[[row_data_name]] <- feature_ids
  col_df <- sample_annotation_from_ids(sample_ids)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = setNames(list(assay_mat), assay_name),
    rowData = S4Vectors::DataFrame(row_df),
    colData = S4Vectors::DataFrame(col_df)
  )
  S4Vectors::metadata(se) <- list(
    datatype = datatype,
    annotation = datatype
  )

  list(se = se, matrix = assay_mat, row_data = row_df)
}


write_matrix_tsv <- function(matrix, row_ids, row_id_name, output_path) {
  out_df <- data.frame(
    row_ids = row_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out_df)[[1]] <- row_id_name
  out_df <- cbind(out_df, as.data.frame(matrix, check.names = FALSE))
  data.table::fwrite(out_df, output_path, sep = "\t")
}


brainarray_package_info <- function(
  cel_file,
  version,
  organism,
  annotation_source
) {
  platform <- affyio::read.celfile.header(cel_file, info = "full")$cdfName
  platform <- gsub("cdf$", "", platform)
  platform <- gsub("stv1$", "st", platform)
  platform <- gsub("stv2$", "st", platform)
  package_name <- paste0(platform, organism, annotation_source, "cdf")
  file_name <- paste0(package_name, "_", version, ".tar.gz")
  list(
    package_name = package_name,
    file_name = file_name,
    url = sprintf(
      "https://mbni.org/customcdf/%s/%s.download/%s",
      version,
      annotation_source,
      file_name
    )
  )
}
