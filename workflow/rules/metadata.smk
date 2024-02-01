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


rule preprocess_metadata:
    input:
        treatmentMetadata = expand(
            procdata / metadata / "{version}_{release}_preprocessed_treatmentMetadata.tsv", 
            version = version, release = release),
        sampleMetadata = expand(
            procdata / metadata / "{version}_{release}_preprocessed_sampleMetadata.tsv", 
            version = version, release = release),
        geneAnnotation = procdata / metadata / "preprocessed_geneAnnotation.tsv"

rule downloadSampleMetadata:
    input:
        sampleMetadata = lambda wc:
            HTTP.remote(config["metadata"][wc.release]['sampleMetadata']['url']),
    output:
        sampleMetadata = metadata / "{version}_{release}_sampleMetadata.xlsx",
    shell:
        """
        mv {input.sampleMetadata} {output.sampleMetadata} 
        """

rule downloadTreatmentMetadata:
    input:
        treatmentMetadata = lambda wc:
            HTTP.remote(config["metadata"][wc.release]['treatmentMetadata']['url']),
    output:
        treatmentMetadata = metadata / "{version}_{release}_treatmentMetadata.csv",
    shell:
        """
        mv {input.treatmentMetadata} {output.treatmentMetadata}
        """

rule downloadMappingFiles:
    input:
        sampleMappingFile = HTTP.remote(config["metadata"]['sampleMappingFile']['url']),
        treatmentMappingFile = HTTP.remote(config["metadata"]['treatmentMappingFile']['url']),
        GDSC_to_CCLE_treatmentMappingFile = HTTP.remote(config["metadata"]['GDSC_to_CCLE']['treatmentMappingFile']['url']),
        GDSC_to_CCLE_sampleMappingFile = HTTP.remote(config["metadata"]['GDSC_to_CCLE']['sampleMappingFile']['url']),
    output:
        sampleMappingFile = metadata / "sampleMappingFile.xlsx",
        treatmentMappingFile = metadata / "treatmentMappingFile.csv",
        GDSC_to_CCLE_treatmentMappingFile = metadata / "GDSC_to_CCLE_treatmentMappingFile.xlsx",
        GDSC_to_CCLE_sampleMappingFile = metadata / "GDSC_to_CCLE_sampleMappingFile.xlsx",
    shell:
        """
        mv {input.sampleMappingFile} {output.sampleMappingFile};
        mv {input.treatmentMappingFile} {output.treatmentMappingFile};
        mv {input.GDSC_to_CCLE_treatmentMappingFile} {output.GDSC_to_CCLE_treatmentMappingFile};
        mv {input.GDSC_to_CCLE_sampleMappingFile} {output.GDSC_to_CCLE_sampleMappingFile}
        """

rule downloadCellModelPassportsMetadata:
    input:
        sampleAnnotation = HTTP.remote(config["molecularProfiles"]['cellModelPassports']['sampleAnnotation']['url']),
        geneAnnotation = HTTP.remote(config["molecularProfiles"]['cellModelPassports']['geneAnnotation']['url']),
    output:
        CMP_sampleAnnotation = metadata / "cellModelPassports_sampleAnnotation.csv",
        CMP_geneAnnotation = metadata / "cellModelPassports_geneAnnotation.csv",
    shell:
        """
        # input files are .csv.gz, so we need to unzip them
        gunzip -c {input.sampleAnnotation} > {output.CMP_sampleAnnotation};
        gunzip -c {input.geneAnnotation} > {output.CMP_geneAnnotation}
        """

rule preprocess_sampleMetadata:
    input:
        sampleMetadata = metadata / "{version}_{release}_sampleMetadata.xlsx",
        GDSC_to_CCLE_sampleMappingFile = metadata / "GDSC_to_CCLE_sampleMappingFile.xlsx",
        CMP_sampleAnnotation = metadata / "cellModelPassports_sampleAnnotation.csv"
    output:
        sampleMetadata = procdata / metadata / "{version}_{release}_preprocessed_sampleMetadata.tsv",
    log: logs / metadata / "{version}_{release}_preprocess_sampleMetadata.log"
    script:
        scripts / metadata / "preprocess_sampleMetadata.R"

rule preprocess_treatmentMetadata:
    input:
        treatmentMetadata = metadata / "{version}_{release}_treatmentMetadata.csv",
        GDSC_to_CCLE_treatmentMappingFile = metadata / "GDSC_to_CCLE_treatmentMappingFile.xlsx"
    output:
        treatmentMetadata = procdata / metadata / "{version}_{release}_preprocessed_treatmentMetadata.tsv",
    log: logs / metadata / "{version}_{release}_preprocess_treatmentMetadata.log"
    script:
        scripts / metadata / "preprocess_treatmentMetadata.R"


# This rule is a wrapper that retrieves the Ensembl annotation for a given species and build.
rule download_EnsemblAnnotation:
    output:
        ensemblAnnotation = metadata / "Ensembl/{}_{}annotation.gtf".format(
            config["metadata"]["referenceGenome"]["build"],
            config["metadata"]["referenceGenome"]["release"],
        ),
    params:
        species = config["metadata"]["referenceGenome"]["species"],
        release = config["metadata"]["referenceGenome"]["release"],
        build = config["metadata"]["referenceGenome"]["build"],
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    log:
        logs / metadata/ "download_EnsemblAnnotation.log",
    wrapper:
        "v3.3.6/bio/reference/ensembl-annotation"

rule preprocess_geneAnnotation:
    input:
        geneAnnotation = metadata / "cellModelPassports_geneAnnotation.csv",
        ensemblAnnotation = rules.download_EnsemblAnnotation.output.ensemblAnnotation
    output:
        geneAnnotation = procdata / metadata / "preprocessed_geneAnnotation.tsv",
    log: logs / metadata / "preprocess_geneAnnotation.log"
    script:
        scripts / metadata / "preprocess_geneAnnotation.R"

