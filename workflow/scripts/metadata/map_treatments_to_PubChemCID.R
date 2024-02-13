## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image("map_treatments_to_PubChem.RData")
    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)

}else{
    load("map_treatments_to_PubChem.RData")
}
library(data.table)
message("Starting map_treatments_to_PubChem.R")
message("INPUT: ", INPUT)
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata)
message("treatmentMetadata: ", paste0(capture.output(str(treatmentMetadata)), collapse = "\n"))


message("Running AnnotationGx::getPubchemCompound on ",nrow(treatmentMetadata), " treatments" )
compound_nameToCIDS <- AnnotationGx::getPubchemCompound(
    treatmentMetadata[1:10, GDSC.treatmentid],
    from='name',
    to='cids',
)

compound_nameToCIDS <- compound_nameToCIDS[!is.na(compound_nameToCIDS$name) & !duplicated(compound_nameToCIDS$name),]
names(compound_nameToCIDS) <- c("GDSC.treatmentid", "pubchem.CID")

# annotations <- c('ChEMBL ID', 'NSC Number', 'Drug Induced Liver Injury', 'CAS', 'ATC Code')
# pubchem.ChEMBL.ID	pubchem.CID	pubchem.NSC.Number	pubchem.DILI.Status	pubchem.CAS.Number	pubchem.ATC.Code
pseudo_annotate <- function(cid, heading){
    result <- AnnotationGx::annotatePubchemCompound(cid = cid, heading = heading)
    if(length(result) == 0 || is.na(result) || is.null(result)){
        return(NA)
    }else{
        return(result)
    }
}

compound_nameToCIDS[, 'pubchem.ChEMBL.ID' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'ChEMBL ID')), by = pubchem.CID]
compound_nameToCIDS[, 'pubchem.NSC.Number' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'NSC Number')), by = pubchem.CID]
compound_nameToCIDS[, 'pubchem.DILI.Status' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'Drug Induced Liver Injury')), by = pubchem.CID]
compound_nameToCIDS[, 'pubchem.CAS.Number' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'CAS')), by = pubchem.CID]
compound_nameToCIDS[, 'pubchem.ATC.Code' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'ATC Code')), by = pubchem.CID]


annotated_treatments <- merge(treatmentMetadata, compound_nameToCIDS, by.x = "GDSC.treatmentid", by.y = "GDSC.treatmentid", all.x = TRUE)

# compound_nameToCIDS[, 'pubchem.ChEMBL.ID' := lapply(pubchem.CID, function(x) AnnotationGx::annotatePubchemCompound(cid = x, heading = 'ChEMBL ID')), by = pubchem.CID]
# compound_nameToCIDS[, 'pubchem.NSC.Number' := lapply(pubchem.CID, function(x) AnnotationGx::annotatePubchemCompound(cid = x, heading = 'NSC Number')), by = pubchem.CID]
# compound_nameToCIDS[, 'pubchem.DILI.Status' := lapply(pubchem.CID, function(x) AnnotationGx::annotatePubchemCompound(cid = x, heading = 'Drug Induced Liver Injury')), by = pubchem.CID]
# compound_nameToCIDS[, 'pubchem.CAS.Number' := lapply(pubchem.CID, function(x) AnnotationGx::annotatePubchemCompound(cid = x, heading = 'CAS')), by = pubchem.CID]
# compound_nameToCIDS[, 'pubchem.ATC.Code' := lapply(pubchem.CID, function(x) AnnotationGx::annotatePubchemCompound(cid = x, heading = 'ATC Code')), by = pubchem.CID]
# compound_nameToCIDS


# compound_nameToCIDS[1:10, pubchem.NSC.Number := AnnotationGx::annotatePubchemCompound(cid = pubchem.CID, heading = 'NSC Number')]
# compound_nameToCIDS[1:10, pubchem.DILI.Status := AnnotationGx::annotatePubchemCompound(cid = pubchem.CID, heading = 'Drug Induced Liver Injury')]
# compound_nameToCIDS[1:10, pubchem.CAS.Number := AnnotationGx::annotatePubchemCompound(cid = pubchem.CID, heading = 'CAS')]
# compound_nameToCIDS[1:10, pubchem.ATC.Code := AnnotationGx::annotatePubchemCompound(cid = pubchem.CID, heading = 'ATC Code')]


data.table::fwrite(annotated_treatments, OUTPUT$treatment_CIDS, quote = FALSE, sep = "\t", row.names = FALSE)

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
