from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "workflow/config/pipeline.yaml"
conda: "workflow/envs/snakemake.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

version = config["GDSC_version"]
release = config["GDSC_release"]
treatmentResponse = config["treatmentResponse"]

rule process_allTreatmentResponse:
    input:
        preprocessed = expand(
            "procdata/treatmentResponse/{version}_{release}_treatmentResponse_preprocessed.qs",
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

rule preprocess_treatmentResponse:
    input:
        rawdata = rawdata / "treatmentResponse/release_{release}/{version}_public_raw_data.csv",
        processed = rawdata / "treatmentResponse/release_{release}/{version}_fitted_dose_response.xlsx",
        treatmentMetadata = procdata / metadata / "{version}_{release}_preprocessed_treatmentMetadata.tsv",
        sampleMetadata = procdata / metadata / "{version}_{release}_preprocessed_sampleMetadata.tsv", 
    output:
        tre = "procdata/treatmentResponse/{version}_{release}_treatmentResponse_preprocessed.qs",
        raw = procdata / "treatmentResponse/{version}_{release}_treatmentResponse_raw.tsv",
        published_profiles = procdata / "treatmentResponse/{version}_{release}_treatmentResponse_published_profiles.tsv",
    log:
        logs / "treatmentResponse/{version}_{release}/preprocess_treatmentResponse.log"
    conda:
        "../envs/treatmentResponse.yaml"
    threads:
        10
    script:
        scripts / "treatmentResponse/preprocess_treatmentResponse.R"