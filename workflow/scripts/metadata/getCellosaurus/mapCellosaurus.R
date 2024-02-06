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

# sampleMetadata <- data.table::fread("procdata/metadata/GDSC2_8.4_preprocessed_sampleMetadata.tsv", sep="\t", header=T)
sampleMetadata <- data.table::fread(INPUT[['sampleMetadata']], sep="\t", header=T)
cellosaurus_object <- INPUT[['cellosaurus_object']]
object <- readRDS(cellosaurus_object)
# /home/bioinf/bhklab/jermiah/repos/GDSC-Pharmacoset_Snakemake/metadata/cellosaurus.RDS
gdsc_mapped <- data.table::as.data.table(Cellosaurus::mapCells(object, unique(sampleMetadata[,GDSC.sampleid])),keep.rownames = T)
gdsc_mapped


cellosaurus_DT <- data.table::as.data.table(object)
columns <- c(
    "cellLineName", "synonyms", "accession", "misspellings",
    # "diseases", 
    "category", "samplingSite", "isCancer", 
    "oncotreeName", "oncotreeTissue", "oncotreeLevel",
    "depmapId", "atccId", "ncbiTaxonomyId", "sangerModelId",
    "sexOfCell", "ageAtSampling"
)

cellosaurus_DT <- cellosaurus_DT[, ..columns]
names(cellosaurus_DT) <- paste0("cellosaurus.", names(cellosaurus_DT))
annotated_sampleMetadata <- merge(
    sampleMetadata[], 
    cellosaurus_DT[!is.na(cellosaurus.sangerModelId), ], 
    by.x = "CMP.model_id", by.y = "cellosaurus.sangerModelId", all.x = TRUE, allow.cartesian = TRUE)

data.table::fwrite(annotated_sampleMetadata, OUTPUT[['sample_Cellosaurus_file']], sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
