configfile: "config/pipeline.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
version = config["GDSC_version"]
release = config["GDSC_release"]

include: "workflow/rules/metadata.smk"
include: "workflow/rules/treatmentResponse.smk"
include: "workflow/rules/rnaseq.smk"
# include: "workflow/rules/microarray.smk"
include: "workflow/rules/cnv.smk"
include: "workflow/rules/mutation.smk"
include: "workflow/rules/fusion.smk"

rule get_PharmacoSet:
    input:
        expand(
            results / "data" / "{version}_{release}/{version}_{release}_PharmacoSet.RDS",
                # version = ["GDSC1", "GDSC2"],#version,
                # release = [8.4, 8.5]#release       
                version = version,
                release = release       
            )

rule build_PharmacoSet:
    input: 
        treatmentResponseExperiment =
            results / "data/treatmentResponse" / "{version}_{release}_treatmentResponseExperiment_fitted.RDS",
        multiAssayExperiment = 
            results / "data" / "{version}_{release}/{version}_{release}_MultiAssayExperiment.RDS",
        annotated_treatmentMetadata =
            procdata / "metadata" / "annotation" / "GDSC_{release}_treatmentMetadata_annotated.tsv",
        annotated_sampleMetadata =
            procdata / "metadata" / "annotation" / "GDSC_{release}_sampleMetadata_mappedCellosaurus.tsv",
    output:
        results / "data" / "{version}_{release}/{version}_{release}_PharmacoSet.RDS",
    log:
        logs / "{version}_{release}_build_PharmacoSet.log",
    conda: "workflow/envs/PharmacoSet.yaml"
    script:
        "scripts/build_PharmacoSet.R"

rule build_MultiAssayExperiment:
    input:
        summarizedExperiment_lists = [
            # rules.make_MICROARRAY_SE.output.microarray_SE,
            rules.make_RNASEQ_SE.output.rse_list,
            rules.make_CNV_SE.output.rse_list,
            rules.make_MUTATION_SE.output.rse_list,
            rules.make_FUSION_SE.output.rse_list,
        ],
        sampleMetadata = rules.annotate_SampleMetadata.output.sampleMetadata,
    output:
        mae = results / "data" / "{version}_{release}/{version}_{release}_MultiAssayExperiment.RDS",
    log:
        logs / "{version}_{release}_build_MultiAssayExperiment.log",
    conda: 
        "workflow/envs/PharmacoSet.yaml"
    script:
        "scripts/build_MultiAssayExperiment.R"

# rule get_Pharmacoset:
#     input:
#         expand(
#             results / "data" / "{version}_{release}/{version}_{release}_Pharmacoset.RDS",
#                 # version = ["GDSC1", "GDSC2"],#version,
#                 # release = [8.4, 8.5]#release       
#                 version = version,
#                 release = release       
#             )

# rule build_Pharmacoset:
#     input:
#         treatmentResponseExperiment = 
#             results / "data/treatmentResponse" / "{version}_{release}_treatmentResponseExperiment_fitted.RDS",
#         annotated_treatmentMetadata = 
#             procdata / "metadata" / "annotation" / "GDSC_{release}_treatmentMetadata_annotated.tsv",
#         annotated_sampleMetadata = 
#             procdata /  "metadata" / "annotation" / "GDSC_{release}_sampleMetadata_mappedCellosaurus.tsv",
#         geneAnnotation = 
#             procdata / metadata / "preprocessed_geneAnnotation.tsv",
#         summarizedExperimentLists = [
#             rules.make_MICROARRAY_SE.output.microarray_SE,
#             rules.make_RNASEQ_SE.output.rse_list,
#             rules.make_CNV_SE.output.rse_list,
#             rules.make_MUTATION_SE.output.rse_list,
#             rules.make_FUSION_SE.output.rse_list,
#     ],
#     output:
#         results / "data" / "{version}_{release}/{version}_{release}_Pharmacoset.RDS",
#     log:
#         logs / "{version}_{release}_build_Pharmacoset.log",
#     conda: "envs/snakemake.yaml"
#     script:
#         "scripts/build_PharmacoSet.R"
