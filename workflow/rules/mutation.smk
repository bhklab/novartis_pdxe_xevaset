from pathlib import Path

procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")


rule make_Mutation_SE:
    input:
        raw=rules.extractSupplementarySheets.output.mutation,
    output:
        se=procdata / "mutation" / "Mutation_SE.rds",
        matrix=procdata / "mutation" / "Mutation_expression.tsv",
    log:
        logs / "mutation" / "make_Mutation_SE.log",
    script:
        scripts / "mutation" / "make_Mutation_SE.R"
