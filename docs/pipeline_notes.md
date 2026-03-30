# Pipeline Notes

This documents notable decisions made during curation.

- `5FFU` is normalized to `5FU`. This appears to be a typo, not a separate treatment.
- Treatment metadata is filtered to only drugs represented in both model and experiment data. This drops orphan entries such as `abraxane + gemcitabine`.
- Treatment metadata is enriched with `AnnotationGx` where possible, including PubChem identifiers and titles, UniChem cross-references, ChEMBL mechanisms, and ChEMBL target names.
- `untreated` is retained as the control treatment in the XevaSet drug table.
- `expDesign` only keeps complete batches with both control and treatment arms. Dropped incomplete batches are logged during the build.
- Duplicate same-day timecourse measurements are aggregated before splitting. Timecourses are only split on true day resets, not repeated day values.
- Missing `tissue` and `tissue.name` values are backfilled from other rows for the same patient before `AnnotationGx` is used to assign OncoTree branch codes, branch names, main types, and ontology identifiers. Remaining blanks are treated as build failures.
- Batch-level `angle`, `abc`, and `TGI` metrics can fail when Xeva interpolation breaks on sparse trajectories. The pipeline retries those batches with `impute.value = FALSE`.
- A small number of sparse batches still return empty values for `angle`, `abc`, and `TGI`, so those remain `NA` in the final XevaSet.
