from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "workflow/config/pipeline.yaml"
conda_env: "../envs/rnaseq.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

version = config["GDSC_version"]
release = config["GDSC_release"]

################################################################################
## RNA-SEQ

rnaseq = config["molecularProfiles"]['rnaseq']

rule download_RNASEQ:
    input:
        processed = HTTP.remote(rnaseq['processed']['url'])
    output:
        processed = rawdata / "rnaseq/rnaseq_all_20220624.zip",
    shell:
        """
        mv {input.processed} {output.processed} 
        """

rule make_RNASEQ_SE:
    input:
        processed = rules.download_RNASEQ.output.processed,
        sampleMetadata = procdata / metadata / f"{version}_{release}_preprocessed_sampleMetadata.tsv",
        geneAnnotation = procdata / metadata / "preprocessed_geneAnnotation.tsv"
    output:
        tpm = procdata / "rnaseq" / "rnaseq_tpm.tsv",
        fpkm = procdata / "rnaseq" / "rnaseq_fpkm.tsv",
        read_count = procdata / "rnaseq" / "rnaseq_read_count.tsv",
        metadata = procdata / "rnaseq" / "rnaseq_metadata.json",
        rse_list = results / "rnaseq" / "rnaseq_rse_list.RDS"
    log:
        logs / "rnaseq/preprocess_RNASEQ.log"
    conda:
        conda_env
    threads:
        3
    script:
        scripts / "rnaseq/make_RNASEQ_SE.R"
