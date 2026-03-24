from pathlib import Path

procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")


rule make_RNASeq_SE:
    input:
        raw=rules.extractSupplementarySheets.output.rnaseq,
    output:
        se=procdata / "rnaseq" / "RNASeq_SE.rds",
        matrix=procdata / "rnaseq" / "RNASeq_expression.tsv",
    log:
        logs / "rnaseq" / "make_RNASeq_SE.log",
    script:
        scripts / "rnaseq" / "make_RNASeq_SE.R"
