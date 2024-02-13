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


mutation_conda_env = "../envs/mutation.yaml"
mutation = config['molecularProfiles']['mutation']

rule download_MUTATION_processed:
    input:
        all_mutations = HTTP.remote(mutation['all_mutations']['url']),
        mutation_genes = HTTP.remote(mutation['Genes_Metadata']['url']),
    output:
        all_mutations = rawdata / "mutation" / "mutations_all_20230202.zip",
        mutation_genes = rawdata / "mutation" / "driver_mutations_20221208.csv",
    log:
        "logs/mutation/download_MUTATION_processed.log",
    shell:
        """
        mv {input.all_mutations} {output.all_mutations} && \
        mv {input.mutation_genes} {output.mutation_genes} > {log} 2>&1
        """

rule make_MUTATION_SE:
    input:
        all_mutations = rules.download_MUTATION_processed.output.all_mutations, 
        mutation_genes = rules.download_MUTATION_processed.output.mutation_genes,
        sampleMetadata = procdata / metadata / f"GDSC_{release}_preprocessed_sampleMetadata.tsv",
        geneAnnotation = procdata / metadata / "preprocessed_geneAnnotation.tsv"
    output:
        rse_list = results / "data/mutation/mutation_rse_list.RDS",
        protein = procdata / "mutation" / "mutation_protein.tsv",
        rna = procdata / "mutation" / "mutation_rna.tsv",
        cdna = procdata / "mutation" / "mutation_cdna.tsv",
        vaf = procdata / "mutation" / "mutation_vaf.tsv",
        effect = procdata / "mutation" / "mutation_effect.tsv",
        metadata = procdata / "mutation" / "mutation_metadata.json",
    log:
        logs / "mutation/preprocess_MUTATION.log",
    conda:
        mutation_conda_env,
    threads:
        6
    script:
        scripts / "mutation/make_MUTATION_SE.R"

# Until figure out how to process VCF files for mutations, use the processed data
# rule download_MUTATION_VCF:
#     input:
#         WGS_VCF = HTTP.remote(mutation['WGS_VCF']['url']),
#         WES_VCF = HTTP.remote(mutation['WES_VCF']['url']),
#     output:
#         WGS_VCF = "rawdata/mutation/mutations_wgs_vcf_20221123.zip",
#         WES_VCF = "rawdata/mutation/mutations_wes_vcf_20221010.zip",
#     shell:
#         """
#         mv {input.WGS_VCF} {output.WGS_VCF} && \
#         mv {input.WES_VCF} {output.WES_VCF}
#         """

