from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

version = config["GDSC_version"]
release = config["GDSC_release"]
treatmentResponse = config["treatmentResponse"]

conda: "workflow/envs/snakemake.yaml"
rule process_allTreatmentResponse:
    input:
        preprocessed = expand(
            results / "data/treatmentResponse" / "{version}_{release}_treatmentResponseExperiment_fitted.RDS",
            version = ["GDSC1", "GDSC2"],
            release = ["8.4", "8.5",]
        )

rule download_treatmentResponse:
    input:
        rawdata = lambda wc: 
            HTTP.remote(treatmentResponse[wc.release][wc.version]["rawdata"]["url"]),
        processed = lambda wc: 
            HTTP.remote(treatmentResponse[wc.release][wc.version]["processed"]["url"])
    output:
        rawdata = rawdata / "treatmentResponse/release_{release}/{version}_public_raw_data.csv",
        processed = rawdata / "treatmentResponse/release_{release}/{version}_fitted_dose_response.xlsx",
    log:
        logs / "treatmentResponse/{version}_{release}/download_treatmentResponse.log"
    shell:
        """
        mv {input.rawdata} {output.rawdata} && \
        mv {input.processed} {output.processed} > {log} 2>&1
        """

rule build_treatmentResponseExperiment:
    input:
        rawdata = rawdata / "treatmentResponse/release_{release}/{version}_public_raw_data.csv",
        processed = rawdata / "treatmentResponse/release_{release}/{version}_fitted_dose_response.xlsx",
        treatmentMetadata = procdata / metadata / "GDSC_{release}_preprocessed_treatmentMetadata.tsv",
        sampleMetadata = procdata / metadata / "GDSC_{release}_preprocessed_sampleMetadata.tsv", 
    output:
        tre = results / "data/treatmentResponse" / "{version}_{release}_treatmentResponseExperiment.RDS",
        raw = procdata / "treatmentResponse" / "{version}_{release}_treatmentResponse_raw.tsv",
        published_profiles = procdata / "treatmentResponse" / "{version}_{release}_treatmentResponse_published_profiles.tsv",
    log:
        logs / "treatmentResponse/{version}_{release}/preprocess_treatmentResponse.log"
    conda:
        "../envs/treatmentResponse.yaml"
    threads:
        8
    resources:
        mem_mb = 16000
    script:
        scripts / "treatmentResponse/build_treatmentResponseExperiment.R"

rule fit_reatmentResponseExperiment:
    input:
        tre = results / "data/treatmentResponse" / "{version}_{release}_treatmentResponseExperiment.RDS",
    output:
        tre_fit = results / "data/treatmentResponse" / "{version}_{release}_treatmentResponseExperiment_fitted.RDS",
    log:
        logs / "treatmentResponse/{version}_{release}/fit_treatmentResponse.log"
    conda:
        "../envs/treatmentResponse.yaml"
    threads:
        30 
    resources:
        mem_mb = 96000
    script:
        scripts / "treatmentResponse/fit_treatmentResponseExperiment.R"