from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import json

HTTP = HTTPRemoteProvider()

configfile: "workflow/config/pipeline.yaml"


rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")


HTTP = HTTPRemoteProvider()

# NEED TO USE CONTAINER INSTEAD
# conda_env = "../envs/microarray.yaml"
microarray_container = "docker://jjjermiah/gdsc_microarray:0.2"

################################################################################
## MICROARRAY 
# todo: this is a pretty exhaustive method to download the data. 
# think of a better solution to download the data. 
rule download_MicroArrayMetadata:
    input:
        srdf = HTTP.remote(config["molecularProfiles"]['microarray']['metadata_srdf']['url']),
        fileList = HTTP.remote(config["molecularProfiles"]['microarray']['metadata_json']['url']),
    output:
        srdf = rawdata / "microarray/E-MTAB-3610.sdrf.txt",
        filelist = rawdata / "microarray/metadata/E-MTAB-3610_filelist.json",
    shell:
        """
        mv {input.srdf} {output.srdf} && mv {input.fileList} {output.filelist}
        """

# This checkpoint rule is used to parse the metadata file and save it as a json file
# The json file is used to download the expression data in downloadExpressionData
checkpoint load_MicroArrayMetadata:
    input:
        filelist = rawdata / "microarray/metadata/E-MTAB-3610_filelist.json",
    output:
        microarrayFiles = rawdata / "microarray/metadata/E-MTAB-3610_expressionFiles.json",
    run:
        import json
        with open(input.filelist) as f:
            filelist = json.load(f)
        
        microarrayFiles = dict()

        for i, file in enumerate(filelist):
            samplename = file['attributes'][0]['value']
            microarrayFiles[samplename] = dict()

            microarrayFiles[samplename]['filename'] = file['path']
            microarrayFiles[samplename]['size'] = file['size']
            microarrayFiles[samplename]['description'] = file['attributes'][1]['value']
       
        # save to json file
        with open(output.microarrayFiles, 'w') as f:
            json.dump(microarrayFiles, f, indent=4)


# input function to get the files from the FTP server and save them locally
def getMicroArrayFiles(wildcards):
    basepath = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/610/E-MTAB-3610/Files/"
    
    with checkpoints.load_MicroArrayMetadata.get().output[0].open() as f:
        arrayFiles = json.load(f)
    
    # get the first 10 keys of the dictionary
    files = dict(list(arrayFiles.items()))

    # append filename to basepath for each sample and return 
    ftpFilePaths = [basepath + arrayFiles[sample]['filename'] for sample in files]

    # return HTTP.remote(ftpFilePaths)
    return [rawdata / "microarray" / arrayFiles[sample]['filename'] for sample in files]

rule make_MICROARRAY_SE:
    input: 
        CELfiles = getMicroArrayFiles,
        CEL_metadata = rawdata / "microarray/E-MTAB-3610.sdrf.txt",
        CEL_FileList = rawdata / "microarray/metadata/E-MTAB-3610_expressionFiles.json",
        sampleMetadata = procdata / metadata / f"GDSC_{release}_preprocessed_sampleMetadata.tsv",
    output:
        microarray_SE = results / "data/microarray/microarray_SE.RDS",
        microarray_expr = procdata / "microarray/microarray_expr.tsv",
        metadata = procdata / "microarray/microarray_metadata.json",
    log: 
        logs / "microarray" / "microarray_SE.log"
    container: 
        microarray_container
    threads:
        4
    script:
        scripts / "microarray" / "make_MICROARRAY_SE.R"

# This rule is primarily to parallelize the download of the ~1018 .CEL files.
# It would be called by the preprocess_MicroArray rule when the input Function getMicroArrayFiles returns a list
# of files that it needs at the rawdata / "microarray" directory.
# the sample wildcard is then used to construct the FTP path to the file and download it using wget
rule download_MicroArrayCEL: 
    output:
        rawdata / "microarray" / "{sample}.cel"
    log:
        logs / "microarray" / "download_MicroArrayCEL" / "Download_{sample}.log"
    retries: 5 # Sometimes it randomly fails to download a file
    shell:
        """
        ftpFilePath="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/610/E-MTAB-3610/Files"
        wget -O {output} $ftpFilePath/{wildcards.sample}.cel > {log} 2>&1
        """