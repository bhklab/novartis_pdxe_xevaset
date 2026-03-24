# Novartis PDXE Snakemake Pipeline

This project builds the Novartis PDXE `XevaSet` using Snakemake and `pixi` for environment management.
The paper supplemental workbook is the only raw data source used by the workflow.

## Pipeline Design Notes

- The workflow is driven by a single source of truth: the Novartis supplemental workbook configured in `config/pipeline.yaml`. The download step accepts either the published URL or a local file path, which keeps development and reruns straightforward.
- Workbook extraction is separated from downstream curation. Relevant sheets are first materialized as tabular intermediates under `data/procdata/supplementary/`, then normalized into annotated model, sample, treatment, experiment, and curve-metric tables.
- RNA-seq, CNV, and mutation data are processed independently into modality-specific `SummarizedExperiment` objects. That keeps assay-specific logic isolated while giving the later integration steps a consistent interface.
- The pipeline builds treatment-response tables and experiment-design metadata explicitly before final assembly. This makes the `XevaSet` construction step easier to inspect and debug because key inputs already exist as named intermediate artifacts.
- Final outputs are layered: `results/PDXE_MultiAssayExperiment.rds` contains the harmonized molecular assays, and `results/Xeva_PDXE.rds` adds the metadata and treatment-response structures required for downstream `Xeva` analysis.

## Setup

```bash
pixi install
pixi run setup
```

## Run Pipeline

Modify the pipeline config at `config/pipeline.yaml` as desired, and then run the snakemake pipeline with:

```bash
pixi run snakemake --cores <N>
```

## Outputs

- `results/PDXE_MultiAssayExperiment.rds`
- `results/Xeva_PDXE.rds`
