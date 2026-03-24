from pathlib import Path

procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")


rule make_CNV_SE:
    input:
        raw=rules.extractSupplementarySheets.output.cnv,
    output:
        se=procdata / "cnv" / "CNV_SE.rds",
        matrix=procdata / "cnv" / "CNV_expression.tsv",
    log:
        logs / "cnv" / "make_CNV_SE.log",
    script:
        scripts / "cnv" / "make_CNV_SE.R"
