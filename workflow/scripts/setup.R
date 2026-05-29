BiocManager::install(
  c("Xeva", "Biostrings"),
  ask = FALSE,
  update = FALSE
)

annotationgx_ref <- Sys.getenv(
  "ANNOTATIONGX_REF",
  "v0.0.0.9097"
)

remotes::install_github(
  paste0("bhklab/AnnotationGx@", annotationgx_ref),
  upgrade = "never"
)
