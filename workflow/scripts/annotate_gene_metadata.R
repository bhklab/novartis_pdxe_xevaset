snakemake@source("helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(httr)
  library(jsonlite)
  library(S4Vectors)
  library(SummarizedExperiment)
})

gencode_cfg <- snk$config[["annotationgx"]][["gencode"]]
species <- gencode_cfg[["species"]]
if (is.null(species) || !nzchar(species)) {
  species <- "homo_sapiens"
}
batch_size <- gencode_cfg[["batch_size"]]
if (is.null(batch_size)) {
  batch_size <- 1000
}
batch_size <- as.integer(batch_size)
request_timeout_seconds <- gencode_cfg[["request_timeout_seconds"]]
if (is.null(request_timeout_seconds)) {
  request_timeout_seconds <- 120
}
request_timeout_seconds <- as.integer(request_timeout_seconds)

profile_inputs <- list(
  RNASeq = list(
    path = snk$input$rnaseq,
    output = snk$output$rnaseq,
    feature_col = "gene"
  ),
  cnv = list(
    path = snk$input$cnv,
    output = snk$output$cnv,
    feature_col = "feature"
  ),
  mutation = list(
    path = snk$input$mutation,
    output = snk$output$mutation,
    feature_col = "gene"
  )
)

empty_gene_table <- function() {
  data.table(
    query_type = character(),
    query_value = character(),
    gene_id = character(),
    gene_version = numeric(),
    symbol = character(),
    description = character(),
    biotype = character(),
    seq_region_name = character(),
    start = integer(),
    end = integer(),
    strand = integer(),
    strand_char = character(),
    assembly_name = character(),
    canonical_transcript = character(),
    object_type = character(),
    species = character(),
    release = integer()
  )
}

strand_to_char <- function(x) {
  dplyr::case_when(
    as.integer(x) == 1L ~ "+",
    as.integer(x) == -1L ~ "-",
    TRUE ~ "*"
  )
}

field <- function(record, name, default = NA) {
  value <- record[[name]]
  if (is.null(value) || !length(value)) {
    return(default)
  }
  value[[1]]
}

get_current_release_info <- function(species) {
  response <- httr::GET(
    "https://rest.ensembl.org/info/data",
    query = list("content-type" = "application/json"),
    httr::accept_json(),
    httr::timeout(request_timeout_seconds)
  )
  httr::stop_for_status(response)
  payload <- jsonlite::fromJSON(
    httr::content(response, as = "text", encoding = "UTF-8")
  )

  data.table(
    species = species,
    ensembl_release = as.integer(payload$releases[[1]]),
    queried_at = as.character(Sys.time())
  )
}

release_info <- get_current_release_info(species)
ensembl_release <- release_info$ensembl_release[[1]]

record_to_row <- function(record, query_value, query_type) {
  if (is.null(record) || !length(record)) {
    return(NULL)
  }

  strand <- field(record, "strand", NA_integer_)
  data.table(
    query_type = query_type,
    query_value = query_value,
    gene_id = field(record, "id", NA_character_),
    gene_version = field(record, "version", NA_real_),
    symbol = field(record, "display_name", NA_character_),
    description = field(record, "description", NA_character_),
    biotype = field(record, "biotype", NA_character_),
    seq_region_name = field(record, "seq_region_name", NA_character_),
    start = field(record, "start", NA_integer_),
    end = field(record, "end", NA_integer_),
    strand = strand,
    strand_char = strand_to_char(strand),
    assembly_name = field(record, "assembly_name", NA_character_),
    canonical_transcript = field(record, "canonical_transcript", NA_character_),
    object_type = field(record, "object_type", NA_character_),
    species = field(record, "species", NA_character_),
    release = ensembl_release
  )
}

lookup_ids <- function(ids) {
  ids <- sort(unique(na.omit(ids)))
  if (!length(ids)) {
    return(empty_gene_table())
  }

  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  rows <- lapply(seq_along(batches), function(batch_idx) {
    batch_ids <- batches[[batch_idx]]
    message(
      "GENCODE ID batch ",
      batch_idx,
      "/",
      length(batches),
      " (",
      length(batch_ids),
      " IDs)"
    )
    response <- httr::POST(
      "https://rest.ensembl.org/lookup/id",
      query = list("content-type" = "application/json"),
      body = list(ids = batch_ids, expand = 0),
      encode = "json",
      httr::accept_json(),
      httr::timeout(request_timeout_seconds)
    )
    httr::stop_for_status(response)
    payload <- jsonlite::fromJSON(
      httr::content(response, as = "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    )

    rbindlist(
      lapply(names(payload), function(id) {
        record_to_row(payload[[id]], id, "id")
      }),
      fill = TRUE
    )
  })

  rbindlist(rows, fill = TRUE)
}

lookup_symbols <- function(symbols) {
  symbols <- sort(unique(na.omit(symbols)))
  if (!length(symbols)) {
    return(empty_gene_table())
  }

  batches <- split(symbols, ceiling(seq_along(symbols) / batch_size))
  rows <- lapply(seq_along(batches), function(batch_idx) {
    batch_symbols <- batches[[batch_idx]]
    message(
      "GENCODE symbol batch ",
      batch_idx,
      "/",
      length(batches),
      " (",
      length(batch_symbols),
      " symbols)"
    )
    response <- httr::POST(
      paste0("https://rest.ensembl.org/lookup/symbol/", species),
      query = list("content-type" = "application/json"),
      body = list(symbols = batch_symbols, expand = 0),
      encode = "json",
      httr::accept_json(),
      httr::timeout(request_timeout_seconds)
    )
    httr::stop_for_status(response)
    payload <- jsonlite::fromJSON(
      httr::content(response, as = "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    )

    rbindlist(
      lapply(names(payload), function(symbol) {
        record_to_row(payload[[symbol]], symbol, "symbol")
      }),
      fill = TRUE
    )
  })

  rbindlist(rows, fill = TRUE)
}

is_ensembl_gene_id <- function(x) {
  grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", x)
}

feature_query_table <- function(se, profile_name, feature_col) {
  row_df <- as.data.frame(rowData(se), stringsAsFactors = FALSE)
  feature_values <- normalize_missing_character(row_df[[feature_col]])
  query_type <- ifelse(is_ensembl_gene_id(feature_values), "id", "symbol")
  query_value <- ifelse(
    query_type == "id",
    sub("\\.[0-9]+$", "", feature_values),
    feature_values
  )

  data.table(
    profile = profile_name,
    row_index = seq_along(feature_values),
    feature_rowname = rownames(row_df),
    feature_column = feature_col,
    feature_value = feature_values,
    query_type = query_type,
    query_value = query_value
  )
}

gencode_columns <- c(
  "gencode.query.type" = "query_type",
  "gencode.query.value" = "query_value",
  "gencode.mapped" = "gencode.mapped",
  "gencode.gene.id" = "gene_id",
  "gencode.gene.version" = "gene_version",
  "gencode.symbol" = "symbol",
  "gencode.description" = "description",
  "gencode.biotype" = "biotype",
  "gencode.seq.region" = "seq_region_name",
  "gencode.start" = "start",
  "gencode.end" = "end",
  "gencode.strand" = "strand_char",
  "gencode.assembly" = "assembly_name",
  "gencode.canonical.transcript" = "canonical_transcript",
  "gencode.object.type" = "object_type",
  "gencode.species" = "species",
  "gencode.release" = "release"
)

se_by_profile <- lapply(profile_inputs, function(profile) readRDS(profile$path))
feature_dt <- rbindlist(
  lapply(names(profile_inputs), function(profile_name) {
    profile <- profile_inputs[[profile_name]]
    feature_query_table(
      se_by_profile[[profile_name]],
      profile_name,
      profile$feature_col
    )
  }),
  fill = TRUE
)

id_queries <- sort(unique(na.omit(feature_dt[query_type == "id"]$query_value)))
symbol_queries <- sort(unique(na.omit(
  feature_dt[query_type == "symbol"]$query_value
)))

message(
  "Annotating ",
  length(id_queries),
  " unique Ensembl gene IDs with current GENCODE"
)
id_hits <- lookup_ids(id_queries)

message(
  "Annotating ",
  length(symbol_queries),
  " unique gene symbols with current GENCODE"
)
symbol_hits <- lookup_symbols(symbol_queries)

gene_hits <- rbindlist(list(id_hits, symbol_hits), fill = TRUE)
gene_hits <- gene_hits[!duplicated(paste(query_type, query_value))]

annotated_dt <- merge(
  feature_dt,
  gene_hits,
  by = c("query_type", "query_value"),
  all.x = TRUE,
  sort = FALSE
)
setorder(annotated_dt, profile, row_index)
annotated_dt[, gencode.mapped := !is.na(gene_id) & nzchar(gene_id)]

unmapped_rows <- list()

for (profile_name in names(profile_inputs)) {
  profile <- profile_inputs[[profile_name]]
  se <- se_by_profile[[profile_name]]
  row_df <- as.data.frame(rowData(se), stringsAsFactors = FALSE)
  source_row_df <- row_df
  profile_annotation <- annotated_dt[profile == profile_name]
  setorder(profile_annotation, row_index)

  for (new_col in names(gencode_columns)) {
    row_df[[new_col]] <- profile_annotation[[gencode_columns[[new_col]]]]
  }

  rowData(se) <- S4Vectors::DataFrame(row_df)
  se_metadata <- metadata(se)
  se_metadata$gencode <- as.list(release_info[1, ])
  metadata(se) <- se_metadata

  dir_create(path_dir(profile$output))
  saveRDS(se, profile$output)

  unmapped <- profile_annotation[gencode.mapped == FALSE]
  if (nrow(unmapped)) {
    source_cols <- as.data.table(
      source_row_df[unmapped$row_index, , drop = FALSE]
    )
    setnames(
      source_cols,
      names(source_cols),
      paste0("source.", names(source_cols))
    )
    unmapped_rows[[profile_name]] <- cbind(
      unmapped[,
        .(
          profile,
          feature_rowname,
          feature_column,
          feature_value,
          query_type,
          query_value
        )
      ],
      source_cols
    )
  }
}

unmapped_genes <- if (length(unmapped_rows)) {
  rbindlist(unmapped_rows, fill = TRUE)
} else {
  data.table(
    profile = character(),
    feature_rowname = character(),
    feature_column = character(),
    feature_value = character(),
    query_type = character(),
    query_value = character(),
    source.gene = character(),
    source.feature = character()
  )
}
dir_create(path_dir(snk$output$unmapped))
data.table::fwrite(unmapped_genes, snk$output$unmapped, sep = "\t", na = "")
