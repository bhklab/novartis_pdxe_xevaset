from pathlib import Path

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

micro_cfg = config["molecularProfiles"]["microarray"]
raw_tar_name = micro_cfg["raw_tar_file"]
sdrf_name = micro_cfg["sdrf_file"]


rule download_MicroarrayRaw:
    output:
        raw_tar=rawdata / "microarray" / raw_tar_name,
    log:
        logs / "microarray" / "download_MicroarrayRaw.log",
    params:
        url=micro_cfg["raw_tar_url"],
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.raw_tar})
        mkdir -p $(dirname {log})
        if [[ "{params.url}" =~ ^https?://|^ftp:// ]]; then
            curl -L "{params.url}" -o "{output.raw_tar}" > "{log}" 2>&1
        else
            cp "{params.url}" "{output.raw_tar}" > "{log}" 2>&1
        fi
        """


rule download_MicroarraySDRF:
    output:
        sdrf=rawdata / "microarray" / sdrf_name,
    log:
        logs / "microarray" / "download_MicroarraySDRF.log",
    params:
        url=micro_cfg["sdrf_url"],
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.sdrf})
        mkdir -p $(dirname {log})
        if [[ "{params.url}" =~ ^https?://|^ftp:// ]]; then
            curl -L "{params.url}" -o "{output.sdrf}" > "{log}" 2>&1
        else
            cp "{params.url}" "{output.sdrf}" > "{log}" 2>&1
        fi
        """


rule make_Microarray_SE:
    input:
        raw_tar=rules.download_MicroarrayRaw.output.raw_tar,
        sdrf=rules.download_MicroarraySDRF.output.sdrf,
    output:
        se=procdata / "microarray" / "Microarray_SE.rds",
        matrix=procdata / "microarray" / "Microarray_expression.tsv",
    log:
        logs / "microarray" / "make_Microarray_SE.log",
    script:
        scripts / "microarray" / "make_Microarray_SE.R"
