from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "workflow/config/pipeline.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])

version = config["GDSC_version"]
release = config["GDSC_release"]
treatmentResponse = config["treatmentResponse"]

rule process_allTreatmentResponse:
    input:
        preprocessed = f"procdata/treatmentResponse/{version}_{release}_treatmentResponse_preprocessed.qs"

rule download_treatmentResponse:
    input:
        rawdata = lambda wc: 
            HTTP.remote(treatmentResponse[wc.release][wc.version]["rawdata"]["url"]),
        processed = lambda wc: 
            HTTP.remote(treatmentResponse[wc.release][wc.version]["processed"]["url"])
    output:
        rawdata = "rawdata/treatmentResponse/release_{release}/{version}_public_raw_data.csv",
        processed = "rawdata/treatmentResponse/release_{release}/{version}_fitted_dose_response.xlsx",
    log:
        "logs/treatmentResponse/{version}_{release}/download_treatmentResponse.log"
    shell:
        """
        mv {input.rawdata} {output.rawdata} && \
        mv {input.processed} {output.processed} > {log} 2>&1
        """

rule preprocess_treatmentResponse:
    input:
        rawdata = "rawdata/treatmentResponse/release_{release}/{version}_public_raw_data.csv",
        processed = "rawdata/treatmentResponse/release_{release}/{version}_fitted_dose_response.xlsx",
        treatmentMetadata = procdata / metadata / "{version}_{release}_preprocessed_treatmentMetadata.tsv",
    output:
        preprocessed = "procdata/treatmentResponse/{version}_{release}_treatmentResponse_preprocessed.qs"
    log:
        "logs/treatmentResponse/{version}_{release}/preprocess_treatmentResponse.log"
    threads:
        4
    script:
        "../scripts/treatmentResponse/preprocess_treatmentResponse.R"