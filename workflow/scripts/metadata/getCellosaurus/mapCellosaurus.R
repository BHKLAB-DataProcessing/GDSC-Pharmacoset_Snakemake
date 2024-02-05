## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image()
    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)

}

sampleMetadata <- read.table(INPUT[['sampleMetadata']], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cellosaurus_object <- INPUT[['cellosaurus_object']]
object <- readRDS(cellosaurus_object)

gdsc_mapped <- Cellosaurus::mapCells(object, sampleMetadata$GDSC.sampleid)

str(gdsc_mapped)