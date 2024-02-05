from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "workflow/config/pipeline.yaml"

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

version = config["GDSC_version"]
release = config["GDSC_release"]

fusion_conda_env = "../envs/fusion.yaml"
fusions = config['molecularProfiles']['fusion']

rule download_FUSION:
    input:
        gene_fusions = HTTP.remote(fusions['gene_fusions']['url']),
    output:
        gene_fusions = rawdata / "fusion/Fusions_20230725.zip"
    log:
        "logs/fusion/download_FUSION.log"
    shell:
        """
        mv {input.gene_fusions} {output.gene_fusions} > {log} 2>&1
        """

rule make_FUSION_SE:
    input:
        gene_fusions = rules.download_FUSION.output.gene_fusions,
        sampleMetadata = procdata / metadata / f"{version}_{release}_preprocessed_sampleMetadata.tsv",
        geneAnnotation = procdata / metadata / "preprocessed_geneAnnotation.tsv"
    output:
        rse_list = results / "data/fusion/fusion_rse_list.RDS",
        fusion = procdata / "fusion/fusion.tsv",
        metadata = procdata / "fusion/fusion_metadata.json"
    log:
        logs / "fusion/make_FUSION_SE.log"
    conda:
        fusion_conda_env
    script:
        scripts / "fusion/make_FUSION_SE.R"
