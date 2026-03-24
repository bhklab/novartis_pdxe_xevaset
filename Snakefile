from pathlib import Path


configfile: "config/pipeline.yaml"


rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])


include: "workflow/rules/metadata.smk"
include: "workflow/rules/rnaseq.smk"
include: "workflow/rules/cnv.smk"
include: "workflow/rules/mutation.smk"
include: "workflow/rules/treatmentResponse.smk"


rule all:
    input:
        xeva=results / "Xeva_PDXE.rds",
    localrule: True


rule build_MultiAssayExperiment:
    input:
        sampleMetadata=rules.annotate_sampleMetadata.output.sampleMetadata,
        assays=[
            rules.make_RNASeq_SE.output.se,
            rules.make_CNV_SE.output.se,
            rules.make_Mutation_SE.output.se,
        ],
    output:
        mae=results / "PDXE_MultiAssayExperiment.rds",
    log:
        logs / "build_MultiAssayExperiment.log",
    script:
        "workflow/scripts/build_MultiAssayExperiment.R"


rule build_XevaSet:
    input:
        multiAssayExperiment=rules.build_MultiAssayExperiment.output.mae,
        modelInfo=rules.annotate_sampleMetadata.output.modelInfo,
        sampleMetadata=rules.annotate_sampleMetadata.output.sampleMetadata,
        treatmentMetadata=rules.annotate_treatmentMetadata.output.treatmentMetadata,
        experiment=rules.build_treatmentResponseTables.output.experiment,
        curveMetrics=rules.build_treatmentResponseTables.output.curveMetrics,
        expDesign=rules.build_experimentDesign.output.expDesign,
        modToBiobaseMap=rules.build_modToBiobaseMap.output.modToBiobaseMap,
    output:
        xeva=results / "Xeva_PDXE.rds",
    log:
        logs / "build_XevaSet.log",
    script:
        "workflow/scripts/build_XevaSet.R"
