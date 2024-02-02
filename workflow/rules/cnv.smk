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

cnv = config['molecularProfiles']['cnv']
cnv_conda_env = "../envs/cnv.yaml"

rule download_CNV_WESData:
    input:
        WES = HTTP.remote(cnv['WES_CNV']['url'])
    output:
        WES_zipped = "rawdata/cnv/WES_pureCN_CNV_genes.zip"
    log:
        "logs/cnv/downloadCNV_WESData.log"
    shell:
        """
        mv {input.WES} {output.WES_zipped} > {log} 2>&1
        """

rule make_CNV_SE:
    input:
        WES_zipped = rules.download_CNV_WESData.output.WES_zipped,
        sampleMetadata = procdata / metadata / f"{version}_{release}_preprocessed_sampleMetadata.tsv",
        geneAnnotation = procdata / metadata / "preprocessed_geneAnnotation.tsv"
    output:
        CNV_se = "results/data/cnv/CNV_SE.qs",
    log:
        "logs/cnv/preprocessCNV.log",
    conda:
        cnv_conda_env,
    threads:
        6
    script:
        "../scripts/cnv/preprocess_CNV.R"



# rule downloadCNV_WGSData:
#     input:
#         WGS = HTTP.remote(molecularProfiles['cnv']['WGS_CNV']['url']),
#     output:
#         WGS_genes = "rawdata/cnv/WGS_purple_CNV_genes_20230303.csv",
#         WGS_category = "rawdata/cnv/WGS_purple_genes_cn_category_20230303.csv",
#         WGS_total_cnv = "rawdata/cnv/WGS_purple_genes_total_copy_number_20230303.csv",
#     shell:
#         "unzip -d $(dirname {output.WGS_genes}) {input.WGS}; rm {input.WGS}"