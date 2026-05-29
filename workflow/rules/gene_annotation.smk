from pathlib import Path

procdata = Path(config["directories"]["procdata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")


rule annotate_gene_metadata:
    input:
        rnaseq=rules.make_RNASeq_SE.output.se,
        cnv=rules.make_CNV_SE.output.se,
        mutation=rules.make_Mutation_SE.output.se,
    output:
        rnaseq=procdata / "gene_annotation" / "RNASeq_SE.rds",
        cnv=procdata / "gene_annotation" / "CNV_SE.rds",
        mutation=procdata / "gene_annotation" / "Mutation_SE.rds",
        unmapped=results / "unmapped_genes.tsv",
    log:
        logs / "gene_annotation" / "annotate_gene_metadata.log",
    script:
        scripts / "annotate_gene_metadata.R"
