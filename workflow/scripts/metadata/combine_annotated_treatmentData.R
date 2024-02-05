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



annotated_CIDS <- 
    lapply(INPUT$annotated_CIDS, function(file) data.table::fread(file, header = TRUE, sep = "\t"))
annotated_ChEMBL <- data.table::fread(INPUT$annotated_ChEMBL, header = TRUE, sep = "\t")
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata, header = TRUE, sep = "\t")

# merge all the annotated_dt on "cid" column 
annotated_DT <- Reduce(function(x, y) merge(x, y, by = c("pubchem.CID", "GDSC.treatmentid"), all = TRUE), annotated_CIDS)
# [1] "cid"                       "ChEMBL ID"                
# [3] "NSC Number"                "Drug Induced Liver Injury"
# [5] "CAS"                       "ATC Code"

oldNames <- c("ChEMBL ID", "NSC Number", "Drug Induced Liver Injury", "CAS", "ATC Code")
newNames <- c("pubchem.ChEMBL.ID", "pubchem.NSC.Number", "pubchem.DILI.Status", "pubchem.CAS.Number", "pubchem.ATC.Code")
data.table::setnames(annotated_DT, oldNames, newNames)


##### rename chembl
# [1] "molecule.chembl.id"        "action.type"              
# [3] "mechanism.of.action"       "molecular.mechanism"      
# [5] "mechanism.comment"         "parent.molecule.chembl.id"
# [7] "target.chembl.id" 
newNames <- c(
    "molecule.chembl.id", "chembl.ActionType", "chembl.MechanismOfAction", "chembl.MolecularMechanism", 
    "chembl.MechanismComment", "chembl.ParentMoleculeChEMBL.ID", "chembl.TargetChEMBL.ID")

oldNames <- c(
    "molecule.chembl.id", "action.type", "mechanism.of.action", 
    "molecular.mechanism", "mechanism.comment", "parent.molecule.chembl.id", "target.chembl.id")

data.table::setnames(annotated_ChEMBL, oldNames, newNames)

annotations_combined_DT <- merge(annotated_DT, annotated_ChEMBL, by.y = "molecule.chembl.id", by.x = "pubchem.ChEMBL.ID", all = TRUE)


# merge the annotations with the treatmentMetadata
final_annotated_treatmentMetadata <- merge(treatmentMetadata, annotations_combined_DT, by = "GDSC.treatmentid", all.x = TRUE)

dir.create(dirname(file.path(OUTPUT$annotated_treatmentMetadata)), recursive = TRUE)
data.table::fwrite(final_annotated_treatmentMetadata, OUTPUT$annotated_treatmentMetadata, quote = FALSE, sep = "\t", row.names = FALSE)
