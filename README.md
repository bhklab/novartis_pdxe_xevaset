# Novartis PDXE Snakemake Pipeline

This project builds the Novartis PDXE `XevaSet` using Snakemake and `pixi` for environment management.
The paper supplemental workbook is the only raw data source used by the workflow.

## Pipeline Design Notes

- The workflow is driven by a single source of truth: the Novartis supplemental workbook configured in `config/pipeline.yaml`. The download step accepts either the published URL or a local file path, which keeps development and reruns straightforward.
- Workbook extraction is separated from downstream curation. Relevant sheets are first materialized as tabular intermediates under `data/procdata/supplementary/`, then normalized into annotated model, sample, treatment, experiment, and curve-metric tables.
- RNA-seq, CNV, and mutation data are processed independently into modality-specific `SummarizedExperiment` objects. Their feature metadata is then enriched with current GENCODE annotations before integration.
- The pipeline builds treatment-response tables and experiment-design metadata explicitly before final assembly. This makes the `XevaSet` construction step easier to inspect and debug because key inputs already exist as named intermediate artifacts.
- Final outputs are layered: `results/PDXE_MultiAssayExperiment.rds` contains the harmonized molecular assays, `results/Xeva_PDXE.rds` adds the metadata and treatment-response structures required for downstream `Xeva` analysis, `results/Xeva_PDXE_tsv/` mirrors the final `XevaSet` contents as TSV exports for non-R workflows, and `results/unmapped_genes.tsv` lists molecular features that did not map to GENCODE.

## Setup

```bash
pixi install
pixi run setup
```

The setup task installs AnnotationGx from the remote `bhklab/AnnotationGx` repository. Set `ANNOTATIONGX_REF` only if a different remote ref should be used.

## Run Pipeline

Modify the pipeline config at `config/pipeline.yaml` as desired, and then run the snakemake pipeline with:

```bash
pixi run snakemake --cores <N>
```

## Outputs

- `results/PDXE_MultiAssayExperiment.rds`
- `results/Xeva_PDXE.rds`
- `results/Xeva_PDXE_tsv/`
- `results/unmapped_genes.tsv`

## Run QC notebook

There is a pixi task setup to quickly knit the QC notebook to HTML in `qc/`.

```bash
pixi run qc
```
