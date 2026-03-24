from pathlib import Path

procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")


rule build_treatmentResponseTables:
    input:
        experiment=rules.preprocessMetadata.output.experiment,
        curveMetrics=rules.preprocessMetadata.output.curveMetrics,
        modelInfo=rules.annotate_sampleMetadata.output.modelInfo,
        treatmentMetadata=rules.annotate_treatmentMetadata.output.treatmentMetadata,
    output:
        experiment=procdata / "treatmentResponse" / "PDXE_experiment.tsv",
        curveMetrics=procdata / "treatmentResponse" / "PDXE_curve_metrics.tsv",
    log:
        logs / "treatmentResponse" / "build_treatmentResponseTables.log",
    script:
        scripts / "treatmentResponse" / "build_treatmentResponseTables.R"
