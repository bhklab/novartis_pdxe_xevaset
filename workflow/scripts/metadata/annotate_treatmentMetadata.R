snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(AnnotationGx)
  library(data.table)
  library(fs)
})

thread_count <- 1L
if (!is.null(snk$threads) && length(snk$threads) > 0) {
  thread_count <- as.integer(snk$threads[[1]])
}

options(
  mc.cores = thread_count,
  annotationgx.mc.cores = thread_count,
  log_level = "WARN"
)

unichem_source_spec <- data.table::data.table(
  source_name = c(
    "chembl",
    "drugbank",
    "chebi",
    "pharmgkb",
    "lincs",
    "clinicaltrials",
    "nih_ncc",
    "fdasrs",
    "rxnorm"
  ),
  output_col = c(
    "unichem.ChEMBL",
    "unichem.DrugBank",
    "unichem.ChEBI",
    "unichem.PharmGKB",
    "unichem.LINCS",
    "unichem.ClinicalTrials",
    "unichem.NIH_NCC",
    "unichem.FDASRS",
    "unichem.RxNorm"
  )
)

component_value_cols <- c(
  "pubchem.CID",
  "pubchem.Title",
  "pubchem.URL",
  "pubchem.MolecularFormula",
  "pubchem.InChIKey",
  "pubchem.CanonicalSMILES",
  unichem_source_spec$output_col,
  "chembl.pref_name",
  "chembl.parent_molecule_chembl_id",
  "chembl.target_chembl_id",
  "chembl.target_pref_name",
  "chembl.target_type",
  "chembl.target_organism",
  "chembl.record_id",
  "chembl.mechanism_of_action",
  "chembl.mechanism_comment",
  "chembl.action_type"
)

normalize_value <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_
  x
}

normalize_display_name <- function(x) {
  x <- normalize_value(x)
  uppercase_idx <- !is.na(x) & grepl("[A-Z]", x) & !grepl("[a-z]", x)
  x[uppercase_idx] <- tools::toTitleCase(tolower(x[uppercase_idx]))
  x
}

collapse_unique_values <- function(x, sep = ",") {
  x <- normalize_value(x)
  x <- unique(x[!is.na(x)])
  if (!length(x)) {
    return(NA_character_)
  }
  paste(x, collapse = sep)
}

split_components <- function(drug_id) {
  if (is.na(drug_id) || !nzchar(trimws(drug_id))) {
    return(character())
  }

  components <- trimws(strsplit(as.character(drug_id), " + ", fixed = TRUE)[[
    1
  ]])
  components[nzchar(components)]
}

split_csv_values <- function(x) {
  x <- normalize_value(x)
  x <- x[!is.na(x)]

  if (!length(x)) {
    return(character())
  }

  values <- trimws(strsplit(x[[1]], ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

collapse_component_values <- function(values) {
  values <- normalize_value(values)

  if (!length(values)) {
    return(NA_character_)
  }

  if (length(values) == 1L) {
    return(values[[1]])
  }

  values[is.na(values)] <- "NA"
  paste(values, collapse = "; ")
}

collapse_display_names <- function(primary_values, fallback_values) {
  fallback_values <- normalize_value(fallback_values)
  if (!length(fallback_values)) {
    return(NA_character_)
  }

  primary_values <- normalize_value(primary_values)
  if (length(primary_values) < length(fallback_values)) {
    primary_values <- c(
      primary_values,
      rep(NA_character_, length(fallback_values) - length(primary_values))
    )
  }

  display_values <- primary_values
  missing_idx <- is.na(display_values)
  display_values[missing_idx] <- fallback_values[missing_idx]
  display_values <- display_values[!is.na(display_values)]

  if (!length(display_values)) {
    return(NA_character_)
  }

  paste(display_values, collapse = " + ")
}

build_pubchem_url <- function(cid) {
  cid <- normalize_value(cid)
  if (is.na(cid)) {
    return(NA_character_)
  }

  sprintf("https://pubchem.ncbi.nlm.nih.gov/compound/%s", cid)
}

chunk_values <- function(x, chunk_size = 50L) {
  if (!length(x)) {
    return(list())
  }

  split(x, ceiling(seq_along(x) / chunk_size))
}

is_trueish_value <- function(x) {
  if (!length(x) || is.null(x) || all(is.na(x))) {
    return(FALSE)
  }

  if (is.logical(x)) {
    return(isTRUE(x[[1]]))
  }

  x_chr <- trimws(as.character(x[[1]]))
  if (!nzchar(x_chr)) {
    return(FALSE)
  }

  tolower(x_chr) %in% c("1", "true", "t", "yes")
}

chembl_name_cache <- new.env(parent = emptyenv())

resolve_chembl_ids_from_names <- function(candidate_names) {
  candidate_names <- normalize_value(candidate_names)
  candidate_names <- unique(candidate_names[!is.na(candidate_names)])

  if (!length(candidate_names)) {
    return(character())
  }

  for (candidate_name in candidate_names) {
    cache_key <- tolower(candidate_name)
    if (exists(cache_key, envir = chembl_name_cache, inherits = FALSE)) {
      cached_ids <- get(cache_key, envir = chembl_name_cache, inherits = FALSE)
      if (length(cached_ids)) {
        return(cached_ids)
      }
      next
    }

    resolved_ids <- tryCatch(
      {
        response <- AnnotationGx::queryChemblAPI(
          resource = "molecule",
          field = "pref_name",
          filter_type = "iexact",
          value = candidate_name
        )

        molecules <- response[["molecules"]]
        if (is.null(molecules) || !length(molecules)) {
          character()
        } else {
          molecule_dt <- data.table::as.data.table(molecules)
          normalize_value(molecule_dt$molecule_chembl_id)
        }
      },
      error = function(err) {
        warn_annotation_issue(
          paste0(
            "AnnotationGx::queryChemblAPI(molecule, ",
            candidate_name,
            ")"
          ),
          err
        )
        character()
      }
    )

    resolved_ids <- unique(resolved_ids[!is.na(resolved_ids)])
    assign(cache_key, resolved_ids, envir = chembl_name_cache)

    if (length(resolved_ids)) {
      return(resolved_ids)
    }
  }

  character()
}

warn_annotation_issue <- function(step, err) {
  warning(
    paste0(step, " failed: ", conditionMessage(err)),
    call. = FALSE
  )
}

fetch_chembl_records <- function(
  resource,
  response_field,
  field,
  ids,
  step,
  chunk_size = 10L
) {
  ids <- normalize_value(ids)
  ids <- unique(ids[!is.na(ids)])

  if (!length(ids)) {
    return(data.table::data.table())
  }

  tryCatch(
    data.table::rbindlist(
      lapply(
        chunk_values(ids, chunk_size = chunk_size),
        function(id_chunk) {
          response <- AnnotationGx::queryChemblAPI(
            resource = resource,
            field = field,
            filter_type = "in",
            value = paste(id_chunk, collapse = ",")
          )

          records <- response[[response_field]]
          if (is.null(records) || !length(records)) {
            return(data.table::data.table())
          }

          data.table::as.data.table(records)
        }
      ),
      fill = TRUE,
      use.names = TRUE
    ),
    error = function(err) {
      warn_annotation_issue(step, err)
      data.table::data.table()
    }
  )
}

build_preferred_chembl_names <- function(chembl_ids) {
  chembl_ids <- normalize_value(chembl_ids)
  chembl_ids <- unique(chembl_ids[!is.na(chembl_ids)])

  if (!length(chembl_ids)) {
    return(data.table::data.table(
      chembl_id = character(),
      `chembl.pref_name` = character()
    ))
  }

  chembl_molecules <- fetch_chembl_records(
    resource = "molecule",
    response_field = "molecules",
    field = "molecule_chembl_id",
    ids = chembl_ids,
    step = "AnnotationGx::queryChemblAPI(molecule)"
  )

  if (!nrow(chembl_molecules)) {
    return(data.table::data.table(
      chembl_id = chembl_ids,
      `chembl.pref_name` = NA_character_
    ))
  }

  for (field in c(
    "molecule_chembl_id",
    "pref_name",
    "therapeutic_flag",
    "availability_type"
  )) {
    if (!(field %in% names(chembl_molecules))) {
      chembl_molecules[, (field) := NA_character_]
    }
  }

  chembl_molecules[,
    molecule_chembl_id := normalize_value(molecule_chembl_id)
  ]
  chembl_molecules[, pref_name := normalize_display_name(pref_name)]

  preferred_name_map <- chembl_molecules[,
    .(parent_pref_name = collapse_unique_values(pref_name)),
    by = .(chembl_id = molecule_chembl_id)
  ]

  chembl_form_links <- fetch_chembl_records(
    resource = "molecule_form",
    response_field = "molecule_forms",
    field = "parent_chembl_id",
    ids = chembl_ids,
    step = "AnnotationGx::queryChemblAPI(molecule_form)"
  )

  preferred_form_name_map <- data.table::data.table(
    chembl_id = character(),
    preferred_form_name = character()
  )

  if (nrow(chembl_form_links) > 0) {
    for (field in c(
      "molecule_chembl_id",
      "parent_chembl_id",
      "is_parent"
    )) {
      if (!(field %in% names(chembl_form_links))) {
        chembl_form_links[, (field) := NA_character_]
      }
    }

    chembl_form_links[,
      molecule_chembl_id := normalize_value(molecule_chembl_id)
    ]
    chembl_form_links[,
      parent_chembl_id := normalize_value(parent_chembl_id)
    ]
    chembl_form_links <- chembl_form_links[
      !vapply(is_parent, is_trueish_value, logical(1))
    ]

    form_chembl_ids <- unique(stats::na.omit(
      chembl_form_links$molecule_chembl_id
    ))
    chembl_form_molecules <- fetch_chembl_records(
      resource = "molecule",
      response_field = "molecules",
      field = "molecule_chembl_id",
      ids = form_chembl_ids,
      step = "AnnotationGx::queryChemblAPI(molecule)"
    )

    if (nrow(chembl_form_molecules) > 0) {
      for (field in c(
        "molecule_chembl_id",
        "pref_name",
        "therapeutic_flag",
        "availability_type"
      )) {
        if (!(field %in% names(chembl_form_molecules))) {
          chembl_form_molecules[, (field) := NA_character_]
        }
      }

      chembl_form_molecules[,
        molecule_chembl_id := normalize_value(molecule_chembl_id)
      ]
      chembl_form_molecules[, pref_name := normalize_display_name(pref_name)]

      chembl_form_candidates <- merge(
        chembl_form_links[, .(parent_chembl_id, molecule_chembl_id)],
        chembl_form_molecules,
        by = "molecule_chembl_id",
        all.x = TRUE,
        sort = FALSE
      )
      chembl_form_candidates <- chembl_form_candidates[!is.na(pref_name)]

      if (nrow(chembl_form_candidates) > 0) {
        chembl_form_candidates[,
          therapeutic_rank := as.integer(vapply(
            therapeutic_flag,
            is_trueish_value,
            logical(1)
          ))
        ]
        chembl_form_candidates[,
          availability_rank := as.integer(
            normalize_value(availability_type) == "1"
          )
        ]
        chembl_form_candidates[, source_order := seq_len(.N)]
        data.table::setorder(
          chembl_form_candidates,
          parent_chembl_id,
          -therapeutic_rank,
          -availability_rank,
          source_order
        )

        preferred_form_name_map <- chembl_form_candidates[,
          .(preferred_form_name = pref_name[[1]]),
          by = .(chembl_id = parent_chembl_id)
        ]
      }
    }
  }

  preferred_name_map <- merge(
    preferred_name_map,
    preferred_form_name_map,
    by = "chembl_id",
    all.x = TRUE,
    sort = FALSE
  )
  preferred_name_map[,
    `chembl.pref_name` := data.table::fcoalesce(
      preferred_form_name,
      parent_pref_name
    )
  ]

  preferred_name_map[, .(chembl_id, `chembl.pref_name`)]
}

message("Loading treatment metadata: ", snk$input$drug)
drug_df <- data.table::fread(snk$input$drug)
treatment_df <- unique(drug_df, by = "drug.id")
treatment_df[,
  `:=`(
    standard.name = data.table::fcoalesce(standard.name, drug.id),
    treatment.target = data.table::fcoalesce(treatment.target, "NA"),
    treatment.type = data.table::fcoalesce(treatment.type, "single")
  )
]
data.table::setorder(treatment_df, drug.id)
treatment_df[, component_list := lapply(drug.id, split_components)]

components <- unique(unlist(treatment_df$component_list, use.names = FALSE))
components <- sort(components)
message("Annotating ", length(components), " unique treatment component(s)")

component_annotations <- data.table::data.table(component = components)

cid_map <- tryCatch(
  {
    mapped <- data.table::as.data.table(
      AnnotationGx::mapCompound2CID(components, first = TRUE)
    )
    mapped[,
      .(
        component = normalize_value(name),
        pubchem.CID = normalize_value(cids)
      )
    ]
  },
  error = function(err) {
    warn_annotation_issue("AnnotationGx::mapCompound2CID", err)
    data.table::data.table(
      component = components,
      pubchem.CID = NA_character_
    )
  }
)
component_annotations <- merge(
  component_annotations,
  cid_map,
  by = "component",
  all.x = TRUE,
  sort = FALSE
)

valid_cids <- unique(stats::na.omit(component_annotations$pubchem.CID))
valid_cid_ints <- suppressWarnings(as.integer(valid_cids))
valid_cid_ints <- unique(stats::na.omit(valid_cid_ints))
message("Mapped ", length(valid_cids), " component(s) to PubChem CID")

pubchem_annotations <- data.table::data.table(`pubchem.CID` = character())
if (length(valid_cids) > 0) {
  pubchem_annotations <- tryCatch(
    {
      props <- data.table::as.data.table(
        AnnotationGx::getPubchemCompound(
          ids = valid_cid_ints,
          from = "cid",
          to = "property",
          properties = c(
            "Title",
            "MolecularFormula",
            "InChIKey",
            "CanonicalSMILES"
          )
        )
      )

      if (
        "ConnectivitySMILES" %in%
          names(props) &&
          !"CanonicalSMILES" %in% names(props)
      ) {
        data.table::setnames(
          props,
          "ConnectivitySMILES",
          "CanonicalSMILES"
        )
      }

      for (col_name in c(
        "CID",
        "Title",
        "MolecularFormula",
        "InChIKey",
        "CanonicalSMILES"
      )) {
        if (!(col_name %in% names(props))) {
          props[, (col_name) := NA_character_]
        }
      }

      props[,
        .(
          `pubchem.CID` = normalize_value(CID),
          `pubchem.Title` = normalize_value(Title),
          `pubchem.MolecularFormula` = normalize_value(MolecularFormula),
          `pubchem.InChIKey` = normalize_value(InChIKey),
          `pubchem.CanonicalSMILES` = normalize_value(CanonicalSMILES)
        )
      ]
    },
    error = function(err) {
      warn_annotation_issue("AnnotationGx::getPubchemCompound", err)
      data.table::data.table(`pubchem.CID` = character())
    }
  )
}

component_annotations <- merge(
  component_annotations,
  unique(pubchem_annotations, by = "pubchem.CID"),
  by = "pubchem.CID",
  all.x = TRUE,
  sort = FALSE
)
component_annotations[,
  `pubchem.URL` := vapply(
    `pubchem.CID`,
    build_pubchem_url,
    character(1)
  )
]

unichem_annotations <- data.table::data.table(`pubchem.CID` = character())
unichem_sources <- tryCatch(
  AnnotationGx::getUnichemSources(TRUE),
  error = function(err) {
    warn_annotation_issue("AnnotationGx::getUnichemSources", err)
    NULL
  }
)

if (!is.null(unichem_sources) && length(valid_cids) > 0) {
  pubchem_source_id <- unichem_sources[Name == "pubchem", SourceID][1]

  if (!is.na(pubchem_source_id)) {
    unichem_results <- tryCatch(
      AnnotationGx::queryUnichemCompound(
        type = "sourceID",
        compound = valid_cids,
        sourceID = pubchem_source_id
      ),
      error = function(err) {
        warn_annotation_issue("AnnotationGx::queryUnichemCompound", err)
        NULL
      }
    )

    if (!is.null(unichem_results)) {
      unichem_annotations <- data.table::rbindlist(
        lapply(
          valid_cids,
          function(cid) {
            row <- as.list(rep(
              NA_character_,
              length(unichem_source_spec$output_col)
            ))
            names(row) <- unichem_source_spec$output_col
            row[["pubchem.CID"]] <- cid

            result <- unichem_results[[cid]]
            if (is.null(result) || inherits(result, "unichem_error")) {
              return(data.table::as.data.table(row))
            }

            mappings <- tryCatch(
              data.table::as.data.table(result$External_Mappings),
              error = function(...) data.table::data.table()
            )

            if (!nrow(mappings)) {
              return(data.table::as.data.table(row))
            }

            mappings <- mappings[Name %in% unichem_source_spec$source_name]
            if (!nrow(mappings)) {
              return(data.table::as.data.table(row))
            }

            mappings <- mappings[,
              .(value = collapse_unique_values(compoundID)),
              by = Name
            ]
            mappings <- merge(
              mappings,
              unichem_source_spec,
              by.x = "Name",
              by.y = "source_name",
              all.x = TRUE,
              sort = FALSE
            )

            for (idx in seq_len(nrow(mappings))) {
              row[[mappings$output_col[[idx]]]] <- mappings$value[[idx]]
            }

            data.table::as.data.table(row)
          }
        ),
        fill = TRUE,
        use.names = TRUE
      )
    }
  } else {
    warning(
      "UniChem source ID for PubChem was not found; skipping UniChem annotation.",
      call. = FALSE
    )
  }
}

component_annotations <- merge(
  component_annotations,
  unique(unichem_annotations, by = "pubchem.CID"),
  by = "pubchem.CID",
  all.x = TRUE,
  sort = FALSE
)
for (col_name in unichem_source_spec$output_col) {
  if (!(col_name %in% names(component_annotations))) {
    component_annotations[, (col_name) := NA_character_]
  }
}
component_annotations[,
  `unichem.ChEMBL` := vapply(
    seq_len(.N),
    function(idx) {
      current_ids <- split_csv_values(`unichem.ChEMBL`[[idx]])
      if (length(current_ids)) {
        return(collapse_unique_values(current_ids))
      }

      fallback_ids <- resolve_chembl_ids_from_names(c(
        `pubchem.Title`[[idx]],
        component[[idx]]
      ))
      collapse_unique_values(fallback_ids)
    },
    character(1)
  )
]

chembl_annotations <- data.table::data.table(component = character())
chembl_component_map <- data.table::rbindlist(
  lapply(
    seq_len(nrow(component_annotations)),
    function(idx) {
      chembl_ids <- split_csv_values(
        component_annotations$`unichem.ChEMBL`[[idx]]
      )
      if (!length(chembl_ids)) {
        return(NULL)
      }

      data.table::data.table(
        component = component_annotations$component[[idx]],
        chembl_id = chembl_ids
      )
    }
  ),
  fill = TRUE,
  use.names = TRUE
)

if (nrow(chembl_component_map) > 0) {
  chembl_pref_name_map <- build_preferred_chembl_names(
    unique(chembl_component_map$chembl_id)
  )
  chembl_annotations <- merge(
    chembl_component_map,
    chembl_pref_name_map,
    by = "chembl_id",
    all.x = TRUE,
    sort = FALSE
  )[,
    .(
      `chembl.pref_name` = collapse_unique_values(`chembl.pref_name`)
    ),
    by = component
  ]

  chembl_fields <- c(
    "parent_molecule_chembl_id",
    "target_chembl_id",
    "record_id",
    "mechanism_of_action",
    "mechanism_comment",
    "action_type"
  )

  chembl_raw <- tryCatch(
    data.table::as.data.table(
      AnnotationGx::getChemblMechanism(unique(chembl_component_map$chembl_id))
    ),
    error = function(err) {
      warn_annotation_issue("AnnotationGx::getChemblMechanism", err)
      data.table::data.table()
    }
  )

  if (nrow(chembl_raw) > 0) {
    chembl_raw[,
      molecule_chembl_id := normalize_value(molecule_chembl_id)
    ]
    chembl_raw[,
      target_chembl_id := normalize_value(target_chembl_id)
    ]
    for (field in chembl_fields) {
      if (!(field %in% names(chembl_raw))) {
        chembl_raw[, (field) := NA_character_]
      }
    }

    chembl_by_id <- chembl_raw[,
      lapply(.SD, collapse_unique_values),
      by = molecule_chembl_id,
      .SDcols = chembl_fields
    ]

    chembl_mechanism_annotations <- merge(
      chembl_component_map,
      chembl_by_id,
      by.x = "chembl_id",
      by.y = "molecule_chembl_id",
      all.x = TRUE,
      sort = FALSE
    )[,
      lapply(.SD, collapse_unique_values),
      by = component,
      .SDcols = chembl_fields
    ]

    data.table::setnames(
      chembl_mechanism_annotations,
      chembl_fields,
      paste0("chembl.", chembl_fields)
    )

    chembl_target_annotations <- data.table::data.table(component = character())
    chembl_target_link <- merge(
      chembl_component_map,
      chembl_raw[, .(
        chembl_id = molecule_chembl_id,
        target_chembl_id = target_chembl_id
      )],
      by = "chembl_id",
      all.x = TRUE,
      sort = FALSE
    )
    chembl_target_ids <- unique(stats::na.omit(
      chembl_target_link$target_chembl_id
    ))

    if (length(chembl_target_ids) > 0) {
      chembl_target_raw <- tryCatch(
        data.table::rbindlist(
          lapply(
            chunk_values(chembl_target_ids, chunk_size = 10L),
            function(target_id_chunk) {
              response <- AnnotationGx::queryChemblAPI(
                resource = "target",
                field = "target_chembl_id",
                filter_type = "in",
                value = paste(target_id_chunk, collapse = ",")
              )

              if (is.null(response$targets) || !length(response$targets)) {
                return(data.table::data.table())
              }

              data.table::as.data.table(response$targets)
            }
          ),
          fill = TRUE,
          use.names = TRUE
        ),
        error = function(err) {
          warn_annotation_issue("AnnotationGx::queryChemblAPI(target)", err)
          data.table::data.table()
        }
      )

      if (nrow(chembl_target_raw) > 0) {
        chembl_target_fields <- c(
          "target_chembl_id",
          "pref_name",
          "target_type",
          "organism"
        )

        for (field in chembl_target_fields) {
          if (!(field %in% names(chembl_target_raw))) {
            chembl_target_raw[, (field) := NA_character_]
          }
        }

        chembl_target_raw[,
          target_chembl_id := normalize_value(target_chembl_id)
        ]
        chembl_target_by_id <- chembl_target_raw[,
          .(
            `chembl.target_pref_name` = collapse_unique_values(pref_name),
            `chembl.target_type` = collapse_unique_values(target_type),
            `chembl.target_organism` = collapse_unique_values(organism)
          ),
          by = target_chembl_id
        ]

        chembl_target_annotations <- merge(
          chembl_target_link,
          chembl_target_by_id,
          by = "target_chembl_id",
          all.x = TRUE,
          sort = FALSE
        )[,
          lapply(.SD, collapse_unique_values),
          by = component,
          .SDcols = c(
            "chembl.target_pref_name",
            "chembl.target_type",
            "chembl.target_organism"
          )
        ]
      }
    }

    chembl_mechanism_annotations <- merge(
      chembl_mechanism_annotations,
      chembl_target_annotations,
      by = "component",
      all.x = TRUE,
      sort = FALSE
    )

    chembl_annotations <- merge(
      chembl_annotations,
      chembl_mechanism_annotations,
      by = "component",
      all.x = TRUE,
      sort = FALSE
    )
  }
}

component_annotations <- merge(
  component_annotations,
  chembl_annotations,
  by = "component",
  all.x = TRUE,
  sort = FALSE
)

for (col_name in component_value_cols) {
  if (!(col_name %in% names(component_annotations))) {
    component_annotations[, (col_name) := NA_character_]
  }
}

treatment_annotations <- data.table::rbindlist(
  lapply(
    seq_len(nrow(treatment_df)),
    function(idx) {
      components_for_treatment <- treatment_df$component_list[[idx]]
      component_rows <- component_annotations[
        match(components_for_treatment, component)
      ]

      row <- list(
        drug.id = treatment_df$drug.id[[idx]],
        `annotationgx.components` = if (length(components_for_treatment)) {
          paste(components_for_treatment, collapse = "; ")
        } else {
          NA_character_
        },
        `annotationgx.componentCount` = as.integer(length(
          components_for_treatment
        )),
        `annotationgx.standardizedName` = collapse_display_names(
          primary_values = component_rows$`chembl.pref_name`,
          fallback_values = data.table::fcoalesce(
            component_rows$`pubchem.Title`,
            components_for_treatment
          )
        )
      )

      for (col_name in component_value_cols) {
        row[[col_name]] <- collapse_component_values(
          component_rows[[col_name]]
        )
      }

      data.table::as.data.table(row)
    }
  ),
  fill = TRUE,
  use.names = TRUE
)

treatment_df[, component_list := NULL]
final_treatment_metadata <- cbind(
  treatment_df,
  treatment_annotations[, !"drug.id"]
)
final_treatment_metadata[,
  standard.name := data.table::fcoalesce(
    `annotationgx.standardizedName`,
    standard.name
  )
]

dir_create(path_dir(snk$output$treatmentMetadata))
data.table::fwrite(
  final_treatment_metadata,
  snk$output$treatmentMetadata,
  sep = "\t"
)
