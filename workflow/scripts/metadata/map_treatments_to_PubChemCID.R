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

treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata)


compound_nameToCIDS <- AnnotationGx::getPubChemCompound(
    treatmentMetadata[1:20, GDSC.treatmentid],
    from='name',
    to='cids',
    batch = FALSE,
    verbose = FALSE,
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS, progressbar = TRUE, stop.on.error = FALSE)
)


compound_nameToCIDS <- compound_nameToCIDS[!is.na(compound_nameToCIDS$name) & !duplicated(compound_nameToCIDS$name),]
names(compound_nameToCIDS) <- c("GDSC.treatmentid", "pubchem.CID")


data.table::fwrite(compound_nameToCIDS, OUTPUT$treatment_CIDS, quote = FALSE, sep = "\t", row.names = FALSE)


# propertiesFromCID <- 
#     AnnotationGx::getPubChemCompound(
#         compound_nameToCIDS[, pubchem.CID], 
#         from='cid', 
#         to='property', 
#         properties=c('Title', 'MolecularFormula', 'InChIKey', 'CanonicalSMILES'),
#         BPPARAM = BiocParallel::MulticoreParam(workers = THREADS, progressbar = TRUE, stop.on.error = FALSE))
# names(propertiesFromCID) <- paste0("pubchem.", names(propertiesFromCID))


# CIDtoSynonyms <- 
#     AnnotationGx::getPubChemCompound(
#         compound_nameToCIDS[, PubChem.CID], 
#         from='cid', 
#         to='synonyms',
#         BPPARAM = BiocParallel::MulticoreParam(workers = THREADS, progressbar = TRUE, stop.on.error = FALSE))
# names(CIDtoSynonyms) <- c("pubchem.CID", "pubchem.Synonyms")

# treatments_annotated_properties <- Reduce(
#     function(x, y) merge(x, y, by = "pubchem.CID", all = TRUE), list(compound_nameToCIDS, propertiesFromCID, CIDtoSynonyms))
