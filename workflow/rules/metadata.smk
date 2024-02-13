from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "workflow/config/pipeline.yaml"
conda_env = "../envs/metadata.yaml"

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

version = config["GDSC_version"]
release = config["GDSC_release"]

annotationGx_docker = "docker://jjjermiah/annotationgx-r:0.0.0.9016"
cellosaurus_docker = "docker://quay.io/biocontainers/r-cellosaurus:0.8.1--r43hdfd78af_0"


rule preprocess_metadata:
    input:
        treatmentMetadata_annot = expand(
            procdata / metadata / "GDSC_{release}_treatmentMetadata_annotated.tsv",
            version = version, release = release),
        sampleMetadata_annot = expand(
            procdata / metadata / "GDSC_{release}_sampleMetadata_mappedCellosaurus.tsv", 
            version = version, release = release),
        geneAnnotation = procdata / metadata / "preprocessed_geneAnnotation.tsv"

rule downloadSampleMetadata:
    input:
        sampleMetadata = lambda wc:
            HTTP.remote(config["metadata"][wc.release]['sampleMetadata']['url']),
    output:
        sampleMetadata = metadata / "GDSC_{release}_sampleMetadata.xlsx",
    shell:
        """
        mv {input.sampleMetadata} {output.sampleMetadata} 
        """

rule downloadTreatmentMetadata:
    input:
        treatmentMetadata = lambda wc:
            HTTP.remote(config["metadata"][wc.release]['treatmentMetadata']['url']),
    output:
        treatmentMetadata = metadata / "GDSC_{release}_treatmentMetadata.csv",
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
        sampleMetadata = metadata / "GDSC_{release}_sampleMetadata.xlsx",
        GDSC_to_CCLE_sampleMappingFile = metadata / "GDSC_to_CCLE_sampleMappingFile.xlsx",
        CMP_sampleAnnotation = metadata / "cellModelPassports_sampleAnnotation.csv"
    output:
        sampleMetadata = procdata / metadata / "GDSC_{release}_preprocessed_sampleMetadata.tsv",
    log: logs / metadata / "GDSC_{release}_preprocess_sampleMetadata.log"
    conda: conda_env
    script:
        scripts / metadata / "preprocess_sampleMetadata.R"

rule preprocess_treatmentMetadata:
    input:
        treatmentMetadata = metadata / "GDSC_{release}_treatmentMetadata.csv",
        GDSC_to_CCLE_treatmentMappingFile = metadata / "GDSC_to_CCLE_treatmentMappingFile.xlsx"
    output:
        treatmentMetadata = procdata / metadata / "GDSC_{release}_preprocessed_treatmentMetadata.tsv",
    log: logs / metadata / "GDSC_{release}_preprocess_treatmentMetadata.log"
    conda: conda_env
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
    conda: conda_env
    script:
        scripts / metadata / "preprocess_geneAnnotation.R"

##### ANNOTATIONS

rule map_treatments_to_PubChemCID:
    input:
        treatmentMetadata = procdata / metadata / "GDSC_{release}_preprocessed_treatmentMetadata.tsv",
    output:
        treatment_CIDS = procdata / metadata / "annotation" / "GDSC_{release}_treatmentMetadata_mapped_PubChem.tsv",
    log: logs / metadata / "GDSC_{release}_map_treatments_to_PubChemCID.log"
    threads:
        8
    container: 
        annotationGx_docker
    script:
        scripts / metadata / "map_treatments_to_PubChemCID.R"

# rule annotate_PubChemCIDS:
#     input:
#         treatment_CIDS = procdata / metadata / "annotation" / "GDSC_{release}_treatmentMetadata_MappedCIDS.tsv",
#     output:
#         annotated_CIDs = procdata / metadata / "annotation" / "GDSC_{release}_CIDS_{annotationType}.tsv",
#     log: logs / metadata / "GDSC_{release}_CIDS_{annotationType}.log"
#     threads:
#         8
#     container: 
#         annotationGx_docker
#     script:
#         scripts / metadata / "annotate_PubChemCIDS.R"

rule annotate_ChEMBL:
    input:
        annotated_CIDS = procdata / metadata / "annotation" / "GDSC_{release}_treatmentMetadata_mapped_PubChem.tsv",
    output:
        annotated_ChEMBL = procdata / metadata / "annotation" / "GDSC_{release}_ChEMBL_annotated.tsv",
    log: logs / metadata / "GDSC_{release}_ChEMBL_annotated.log"
    threads:
        8
    container: 
        annotationGx_docker
    script:
        scripts / metadata / "annotate_ChEMBL.R"


# PubChemAnnotations = ['ChEMBL ID', 'NSC Number', 'Drug Induced Liver Injury', 'CAS', 'ATC Code']
rule annotate_TreatmentData:
    input:
        annotated_ChEMBL = procdata / metadata / "annotation" / "GDSC_{release}_ChEMBL_annotated.tsv",
        treatmentMetadata = procdata / metadata / "GDSC_{release}_treatmentMetadata_mapped_PubChem.tsv",
    output:
        annotated_treatmentMetadata = results / "data" / "metadata" / "GDSC_{release}_treatmentMetadata_annotated.tsv",
    log: logs / metadata / "GDSC_{release}_treatmentMetadata_annotated.log"
    container: 
        annotationGx_docker
    script:
        scripts / metadata / "combine_annotated_treatmentData.R"


rule getCellosaurusObject:
    output:
        cellosaurus_object = "metadata/cellosaurus.RDS",
    container: 
        cellosaurus_docker
    script:
        scripts / metadata / "getCellosaurus/getCellosaurusObject.R"

rule annotate_SampleMetadata:
    input:
        sampleMetadata = procdata / metadata / "GDSC_{release}_preprocessed_sampleMetadata.tsv",
        cellosaurus_object = "metadata/cellosaurus.RDS", 
    output:
        sample_Cellosaurus_file = results / "data" / "metadata" / "GDSC_{release}_sampleMetadata_mappedCellosaurus.tsv",
    container: 
        cellosaurus_docker
    script:
        scripts / metadata / "getCellosaurus/mapCellosaurus.R"