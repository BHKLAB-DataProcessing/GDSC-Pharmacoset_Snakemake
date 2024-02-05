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

annotationType <- WILDCARDS$annotationType
compound_nameToCIDS <- data.table::fread(INPUT$treatment_CIDS, header = TRUE, sep = "\t")

### EXTERNAL ANNOTATIONS
stopifnot(annotationType %in% c('ChEMBL ID', 'NSC Number', 'Drug Induced Liver Injury', 'CAS', 'ATC Code'))

# Create a BPParam object
bpParam <- BiocParallel::MulticoreParam(workers = THREADS, progressbar = TRUE, stop.on.error = FALSE)

# Create a partial function to pass to bplapply
partialFunction <- function(CID, annotationType) {
    data.table::as.data.table(
        AnnotationGx::getPubChemAnnotation(CID, annotationType = annotationType))
}

# create a partial function that runs bptry and bplapply
bp_func <- function(CID, annotationType, BPPARAM) {
    result <- BiocParallel::bptry(
        BiocParallel::bplapply(
            CID, 
            partialFunction,
            annotationType = annotationType,
            BPPARAM = BPPARAM
        )
    )
    # Get the successful runs
    successful <- result[BiocParallel::bpok(result)]
    failed <- result[!BiocParallel::bpok(result)] # unused

    data.table::rbindlist(successful)[!duplicated(cid),]
}




# Get the CID to annotation mapping
annotated_CIDs <- bp_func(
    compound_nameToCIDS[, pubchem.CID], 
    annotationType, 
    bpParam
)

# Merge back with compound_nameToCIDS
annotated_CIDs <- merge(
    compound_nameToCIDS, 
    annotated_CIDs, 
    by.x = "pubchem.CID", 
    by.y = "cid",
    all.x = TRUE
)


data.table::fwrite(
    annotated_CIDs,
    OUTPUT$annotated_CIDs,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
)