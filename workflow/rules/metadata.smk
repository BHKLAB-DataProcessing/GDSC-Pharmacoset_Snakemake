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

annotationGx_docker = "docker://bhklab/annotationgx-r:0.0.0.9095"

################################################################################################
# MOTHER RULES
################################################################################################
rule run_annotate_treatmentMetadata:
    input:
        treatment_CIDS = expand(
            procdata / metadata / "annotation" / "GDSC_{release}_treatmentMetadata_annotated.tsv",
            release = release), 

rule preprocess_metadata:
    input:
        treatmentMetadata_annot = expand(
            results / "data" / "metadata" / "GDSC_{release}_treatmentMetadata_annotated.tsv",
            release = release),
        sampleMetadata_annot = expand(
            results / "data" / "metadata" / "GDSC_{release}_sampleMetadata_mappedCellosaurus.tsv",
            release = release),
        geneAnnotation = procdata / metadata / "preprocessed_geneAnnotation.tsv"

rule download_ONLY:
    input:
        sampleMetadata = metadata / "GDSC_{release}_sampleMetadata.xlsx",
        CMP_sampleAnnotation = metadata / "cellModelPassports_sampleAnnotation.csv",
        CMP_geneAnnotation = metadata / "cellModelPassports_geneAnnotation.csv",
        ensemblAnnotation = metadata / "Ensembl/{}_{}annotation.gtf".format(
            config["metadata"]["referenceGenome"]["build"],
            config["metadata"]["referenceGenome"]["release"],
        ),


################################################################################################
# DOWNLOAD RULES
################################################################################################
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

################################################################################################
# SAMPLE METADATA RULES
################################################################################################

rule preprocess_sampleMetadata:
    input:
        sampleMetadata = metadata / "GDSC_{release}_sampleMetadata.xlsx",
        CMP_sampleAnnotation = metadata / "cellModelPassports_sampleAnnotation.csv"
    output:
        sampleMetadata = procdata / metadata / "GDSC_{release}_preprocessed_sampleMetadata.tsv",
    log: logs / metadata / "GDSC_{release}_preprocess_sampleMetadata.log"
    conda: conda_env
    script:
        scripts / metadata / "preprocess_sampleMetadata.R"


rule annotate_SampleMetadata:
    input:
        sampleMetadata = procdata / metadata / "GDSC_{release}_preprocessed_sampleMetadata.tsv"
    output:
        sample_Cellosaurus_file = results / "data" / "metadata" / "GDSC_{release}_sampleMetadata_mappedCellosaurus.tsv",
    container: 
        annotationGx_docker
    threads:
        4
    script:
        scripts / metadata / "getCellosaurus/mapCellosaurus.R"

################################################################################################
# GENE ANNOTATION RULES
################################################################################################
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

################################################################################################
# TREATMENT METADATA RULES
################################################################################################

rule preprocess_treatmentMetadata:
    input:
        treatmentMetadata = metadata / "GDSC_{release}_treatmentMetadata.csv",
    output:
        treatmentMetadata = procdata / metadata / "GDSC_{release}_preprocessed_treatmentMetadata.tsv",
    log: logs / metadata / "GDSC_{release}_preprocess_treatmentMetadata.log"
    conda: conda_env
    script:
        scripts / metadata / "preprocess_treatmentMetadata.R"

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

rule annotate_treatmentMetadata:
    input:
        annotated_CIDS = procdata / metadata / "annotation" / "GDSC_{release}_treatmentMetadata_mapped_PubChem.tsv",
    output:
        annotated_treatmentMetadata = procdata / metadata / "annotation" / "GDSC_{release}_treatmentMetadata_annotated.tsv",
    log: 
        logs / metadata / "GDSC_{release}_annotate_treatmentMetadata.log"
    threads:
        8
    container: 
        annotationGx_docker
    script:
        scripts / metadata / "annotate_treatmentMetadata.R"
