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


HTTP = HTTPRemoteProvider()


################################################################################
## MICROARRAY 
# todo: this is a pretty exhaustive method to download the data. 
# think of a better solution to download the data. 
rule download_MicroArrayMetadata:
    input:
        srdf = HTTP.remote(molecularProfiles['microarray']['metadata_srdf']['url']),
        fileList = HTTP.remote(molecularProfiles['microarray']['metadata_json']['url']),
        tsv = HTTP.remote(molecularProfiles['microarray']['metadata_tsv']['url'])
    output:
        srdf = rawdata / "microarray/E-MTAB-3610.sdrf.txt",
        filelist = rawdata / "microarray/E-MTAB-3610_filelist.json",
        tsv = rawdata / "microarray/E-MTAB-3610.tsv",
    shell:
        """
        mv {input.srdf} {output.srdf} && mv {input.fileList} {output.filelist} && mv {input.tsv} {output.tsv}
        """

# This checkpoint rule is used to parse the metadata file and save it as a json file
# The json file is used to download the expression data in downloadExpressionData
checkpoint load_MicroArrayMetadata:
    input:
        filelist = rawdata / "microarray/E-MTAB-3610_filelist.json",
    output:
        microarrayFiles = rawdata / "microarray/E-MTAB-3610_expressionFiles.json",
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
    return [rawdata / "microarray/" + arrayFiles[sample]['filename'] for sample in files]

rule preprocess_MicroArray:
    input: 
        CELfiles = getMicroArrayFiles,
        CEL_metadata = rawdata / "microarray/E-MTAB-3610.sdrf.txt",
        CEL_FileList = rawdata / "microarray/E-MTAB-3610_expressionFiles.json",
    output:
        preprocessedMicroArray = procdata / "microarray/preprocessedMicroArray.qs",
    threads:
        4
    script:
        scripts / "microarray/preprocess_MICROARRAY.R"

rule download_MicroArrayFILE: 
    output:
        rawdata / "microarray/{sample}.cel"
    log:
        logs / "microarray/Download_{sample}.log"
    retries: 5
    shell:
        """
        ftpFilePath="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/610/E-MTAB-3610/Files"
        wget -O {output} $ftpFilePath/{wildcards.sample}.cel > {log} 2>&1
        """