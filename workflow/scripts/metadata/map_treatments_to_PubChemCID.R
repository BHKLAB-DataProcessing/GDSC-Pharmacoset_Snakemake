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
    treatmentMetadata[1:100, GDSC.treatmentid],
    from='name',
    to='cids',
)
compound_nameToCIDS[cids == "PUGREST.NotFound", cids := NA_character_]
compound_nameToCIDS <- compound_nameToCIDS[!is.na(compound_nameToCIDS$name) & !duplicated(compound_nameToCIDS$name) & !is.na(cids),]
names(compound_nameToCIDS) <- c("GDSC.treatmentid", "pubchem.CID")

### Get properties from CID
properties=c('Title', 'MolecularFormula', 'InChIKey', 'CanonicalSMILES')
message("Getting the following properties from PubChem: ", paste(properties, collapse= " "), " for ", nrow(compound_nameToCIDS), " compounds")
pubchemProperties <- compound_nameToCIDS[, AnnotationGx::getPubchemCompound(ids = pubchem.CID, from = 'cid', to = 'property', properties= properties)]
names(pubchemProperties) <- paste0("pubchem.", names(pubchemProperties))
pubchemProperties[, pubchem.CID := as.character(pubchem.CID)]

# set each pubchem.CID column to character


pubchem_annotated <- merge(compound_nameToCIDS, pubchemProperties, by= "pubchem.CID", all.x = TRUE)

pseudo_annotate <- function(cid, heading){
    result <- AnnotationGx::annotatePubchemCompound(cid = cid, heading = heading)
    if(length(result) == 0 || is.na(result) || is.null(result)){
        return(NA)
    }else{
        return(result)
    }
}


# annotations <- c('ChEMBL ID', 'NSC Number', 'Drug Induced Liver Injury', 'CAS', 'ATC Code')
message("Annotating with ChEMBL ID")
pubchem_annotated[, 'pubchem.ChEMBL.ID' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'ChEMBL ID')), by = pubchem.CID]

message("Annotating with NSC Number")
pubchem_annotated[, 'pubchem.NSC.Number' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'NSC Number')), by = pubchem.CID]

message("Annotating with Drug Induced Liver Injury")
pubchem_annotated[, 'pubchem.DILI.Status' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'Drug Induced Liver Injury')), by = pubchem.CID]

message("Annotating with CAS Number")
pubchem_annotated[, 'pubchem.CAS.Number' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'CAS')), by = pubchem.CID]

message("Annotating with ATC Code")
pubchem_annotated[, 'pubchem.ATC.Code' := lapply(pubchem.CID, function(x) pseudo_annotate(x, 'ATC Code')), by = pubchem.CID]

message("Annotating with Synonyms")
annotated_treatments <- merge(treatmentMetadata, pubchem_annotated, by.x = "GDSC.treatmentid", by.y = "GDSC.treatmentid", all.x = TRUE)

message("Writing to ", OUTPUT$treatment_CIDS)
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
