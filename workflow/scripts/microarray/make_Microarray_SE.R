snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(data.table)
  library(affyio)
  library(affy)
  library(fs)
  library(R.utils)
  library(biomaRt)
  library(SummarizedExperiment)
  library(S4Vectors)
})

micro_cfg <- snk$config[["molecularProfiles"]][["microarray"]]
if (grepl("\\.rds$", snk$input$raw_tar, ignore.case = TRUE)) {
  obj <- readRDS(snk$input$raw_tar)
  se <- if (methods::is(obj, "SummarizedExperiment")) {
    obj
  } else {
    expression_set_to_se(obj, "microarray")
  }
  expr_matrix <- SummarizedExperiment::assay(se)
} else {
  tmp_root <- path(
    snk$config[["directories"]][["procdata"]],
    "_tmp",
    "microarray",
    micro_cfg[["geo_accession"]]
  )
  cel_dir <- path(tmp_root, "cel")
  dir_create(cel_dir)

  message("Extracting microarray archive: ", snk$input$raw_tar)
  utils::untar(snk$input$raw_tar, exdir = cel_dir)

  gz_files <- dir_ls(
    cel_dir,
    regexp = "(?i)\\.cel\\.gz$",
    type = "file"
  )
  if (length(gz_files)) {
    for (gz_file in gz_files) {
      target <- path_ext_remove(gz_file)
      if (!file_exists(target)) {
        R.utils::gunzip(
          gz_file,
          destname = target,
          overwrite = TRUE,
          remove = FALSE
        )
      }
    }
  }

  cel_files <- dir_ls(
    cel_dir,
    regexp = "(?i)\\.cel$",
    type = "file"
  )
  if (!length(cel_files)) {
    stop(
      "No CEL files found after extracting ",
      snk$input$raw_tar,
      call. = FALSE
    )
  }

  pkg_info <- brainarray_package_info(
    cel_file = cel_files[[1]],
    version = micro_cfg[["brainarray"]][["version"]],
    organism = micro_cfg[["brainarray"]][["organism"]],
    annotation_source = micro_cfg[["brainarray"]][["annotation_source"]]
  )

  if (!requireNamespace(pkg_info$package_name, quietly = TRUE)) {
    install.packages(pkg_info$url, repos = NULL, type = "source")
  }
  library(pkg_info$package_name, character.only = TRUE)

  message("Normalizing microarray CEL files with ", pkg_info$package_name)
  raw_affy <- affy::ReadAffy(
    filenames = cel_files,
    cdfname = pkg_info$package_name
  )
  norm_affy <- affy::expresso(
    raw_affy,
    bg.correct = FALSE,
    normalize = TRUE,
    normalize.method = "quantiles",
    pmcorrect.method = "pmonly",
    summary.method = "medianpolish"
  )

  expr_matrix <- Biobase::exprs(norm_affy)
  rownames(expr_matrix) <- gsub("_at$", "", rownames(expr_matrix))
  colnames(expr_matrix) <- path_ext_remove(path_file(colnames(expr_matrix)))
  colnames(expr_matrix) <- sub("_.*$", "", colnames(expr_matrix))

  sdrf <- data.table::fread(snk$input$sdrf, sep = "\t", data.table = FALSE)
  colnames(sdrf) <- clean_names_simple(colnames(sdrf))

  required_sdrf <- c(
    "array_data_file",
    "characteristics_model_number",
    "characteristics_passage",
    "characteristics_primary_site"
  )
  missing_sdrf <- setdiff(required_sdrf, colnames(sdrf))
  if (length(missing_sdrf)) {
    stop(
      "Missing SDRF columns: ",
      paste(missing_sdrf, collapse = ", "),
      call. = FALSE
    )
  }

  sample_df <- sdrf |>
    transmute(
      biobase.id = sub(
        "_.*$",
        "",
        path_ext_remove(path_file(array_data_file))
      ),
      patient.id = paste0("X-", as.character(characteristics_model_number)),
      passage = as.character(characteristics_passage),
      tumor.type = as.character(characteristics_primary_site)
    ) |>
    filter(biobase.id %in% colnames(expr_matrix)) |>
    distinct(biobase.id, .keep_all = TRUE) |>
    arrange(patient.id, passage, biobase.id)

  expr_matrix <- expr_matrix[, sample_df$biobase.id, drop = FALSE]

  ensembl <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = micro_cfg[["ensembl"]][["dataset"]],
    host = micro_cfg[["ensembl"]][["host"]]
  )

  gene_map <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = unique(rownames(expr_matrix)),
    mart = ensembl
  ) |>
    as.data.frame(stringsAsFactors = FALSE) |>
    filter(hgnc_symbol != "") |>
    distinct(ensembl_gene_id, .keep_all = TRUE) |>
    distinct(hgnc_symbol, .keep_all = TRUE)

  common_ids <- intersect(rownames(expr_matrix), gene_map$ensembl_gene_id)
  gene_map <- gene_map[
    match(common_ids, gene_map$ensembl_gene_id),
    ,
    drop = FALSE
  ]
  expr_matrix <- expr_matrix[common_ids, , drop = FALSE]
  rownames(expr_matrix) <- gene_map$hgnc_symbol

  row_df <- data.frame(
    geneName = gene_map$hgnc_symbol,
    ensembl.id = gene_map$ensembl_gene_id,
    stringsAsFactors = FALSE,
    row.names = gene_map$hgnc_symbol
  )
  row.names(sample_df) <- sample_df$biobase.id

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = expr_matrix),
    rowData = S4Vectors::DataFrame(row_df),
    colData = S4Vectors::DataFrame(sample_df)
  )
}
S4Vectors::metadata(se) <- list(
  datatype = "microarray",
  annotation = "microarray"
)

dir_create(unique(path_dir(c(snk$output$se, snk$output$matrix))))
saveRDS(se, snk$output$se)
data.table::fwrite(
  data.table::as.data.table(expr_matrix, keep.rownames = "geneName"),
  snk$output$matrix,
  sep = "\t"
)
