from pathlib import Path

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

metadata_cfg = config["metadata"]
supplementary_name = metadata_cfg["supplementary"]["file"]
sheet_cfg = metadata_cfg["supplementary"]["sheets"]


rule downloadSupplementaryWorkbook:
    output:
        workbook=rawdata / supplementary_name,
    log:
        logs / "metadata" / "downloadSupplementaryWorkbook.log",
    params:
        url=metadata_cfg["supplementary"]["url"],
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.workbook})
        mkdir -p $(dirname {log})
        if [[ "{params.url}" =~ ^https?://|^ftp:// ]]; then
            curl -L "{params.url}" -o "{output.workbook}" > "{log}" 2>&1
        else
            cp "{params.url}" "{output.workbook}" > "{log}" 2>&1
        fi
        """


rule extractSupplementarySheets:
    input:
        workbook=rules.downloadSupplementaryWorkbook.output.workbook,
    output:
        rnaseq=procdata / "supplementary" / sheet_cfg["rnaseq"]["file"],
        cnv=procdata / "supplementary" / sheet_cfg["cnv"]["file"],
        mutation=procdata / "supplementary" / sheet_cfg["mutation"]["file"],
        pctRaw=procdata / "supplementary" / sheet_cfg["pct_raw"]["file"],
        pctCurveMetrics=procdata
        / "supplementary"
        / sheet_cfg["pct_curve_metrics"]["file"],
    log:
        logs / "metadata" / "extractSupplementarySheets.log",
    script:
        scripts / "metadata" / "extractSupplementarySheets.R"


rule preprocessMetadata:
    input:
        pctRaw=rules.extractSupplementarySheets.output.pctRaw,
        pctCurveMetrics=rules.extractSupplementarySheets.output.pctCurveMetrics,
    output:
        model=procdata / "metadata" / "preprocessed_model.tsv",
        drug=procdata / "metadata" / "preprocessed_drug.tsv",
        experiment=procdata / "metadata" / "preprocessed_experiment.tsv",
        curveMetrics=procdata / "metadata" / "preprocessed_curve_metrics.tsv",
    log:
        logs / "metadata" / "preprocessMetadata.log",
    script:
        scripts / "metadata" / "preprocessMetadata.R"


rule annotate_treatmentMetadata:
    input:
        drug=rules.preprocessMetadata.output.drug,
        curveMetrics=rules.preprocessMetadata.output.curveMetrics,
    output:
        treatmentMetadata=procdata
        / "metadata"
        / "annotations"
        / "PDXE_treatmentMetadata_annotated.tsv",
    log:
        logs / "metadata" / "annotate_treatmentMetadata.log",
    script:
        scripts / "metadata" / "annotate_treatmentMetadata.R"


rule annotate_sampleMetadata:
    input:
        model=rules.preprocessMetadata.output.model,
    output:
        modelInfo=procdata
        / "metadata"
        / "annotations"
        / "PDXE_modelMetadata_annotated.tsv",
        sampleMetadata=procdata
        / "metadata"
        / "annotations"
        / "PDXE_sampleMetadata_annotated.tsv",
    log:
        logs / "metadata" / "annotate_sampleMetadata.log",
    script:
        scripts / "metadata" / "annotate_sampleMetadata.R"


rule build_experimentDesign:
    input:
        modelInfo=rules.annotate_sampleMetadata.output.modelInfo,
    output:
        expDesign=procdata / "metadata" / "annotations" / "PDXE_expDesign.rds",
    log:
        logs / "metadata" / "build_experimentDesign.log",
    script:
        scripts / "metadata" / "build_experimentDesign.R"


rule build_modToBiobaseMap:
    input:
        modelInfo=rules.annotate_sampleMetadata.output.modelInfo,
        rnaseq=procdata / "rnaseq" / "RNASeq_SE.rds",
        cnv=procdata / "cnv" / "CNV_SE.rds",
        mutation=procdata / "mutation" / "Mutation_SE.rds",
    output:
        modToBiobaseMap=procdata
        / "metadata"
        / "annotations"
        / "PDXE_modToBiobaseMap.tsv",
    log:
        logs / "metadata" / "build_modToBiobaseMap.log",
    script:
        scripts / "metadata" / "build_modToBiobaseMap.R"
