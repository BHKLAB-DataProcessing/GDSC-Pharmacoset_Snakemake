## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image("combine_annotated_treatmentData.RData")
    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)

}
q()


annotated_CIDS <- 
    lapply(INPUT$annotated_CIDS, function(file) data.table::fread(file, header = TRUE, sep = "\t"))
annotated_ChEMBL <- data.table::fread(INPUT$annotated_ChEMBL, header = TRUE, sep = "\t")
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata, header = TRUE, sep = "\t")

# # merge all the annotated_dt on "cid" column 
# annotated_DT <- Reduce(function(x, y) merge(x, y, by = c("pubchem.CID", "GDSC.treatmentid"), all = TRUE), annotated_CIDS)
# # [1] "cid"                       "ChEMBL ID"                
# # [3] "NSC Number"                "Drug Induced Liver Injury"
# # [5] "CAS"                       "ATC Code"

# oldNames <- c("ChEMBL ID", "NSC Number", "Drug Induced Liver Injury", "CAS", "ATC Code")
# newNames <- c("pubchem.ChEMBL.ID", "pubchem.NSC.Number", "pubchem.DILI.Status", "pubchem.CAS.Number", "pubchem.ATC.Code")
# data.table::setnames(annotated_DT, oldNames, newNames)


##### rename chembl
# [1] "Molecule Chembl ID"        "Action Type"              
# [3] "Mechanism Of Action"       "Molecular Mechanism"      
# [5] "Mechanism Comment"         "Parent Molecule Chembl ID"
# [7] "Target Chembl ID"  
oldNames <- c("Molecule Chembl ID", "Action Type", "Mechanism Of Action", "Molecular Mechanism", "Mechanism Comment", "Parent Molecule Chembl ID", "Target Chembl ID")
newNames <- c(
    "molecule.chembl.id", "chembl.ActionType", "chembl.MechanismOfAction", "chembl.MolecularMechanism", 
    "chembl.MechanismComment", "chembl.ParentMoleculeChEMBL.ID", "chembl.TargetChEMBL.ID")


data.table::setnames(annotated_ChEMBL, oldNames, newNames)

# annotations_combined_DT <- merge(annotated_DT, annotated_ChEMBL, by.y = "molecule.chembl.id", by.x = "pubchem.ChEMBL.ID", all = TRUE)
# merge the annotations with the treatmentMetadata
final_annotated_treatmentMetadata <- merge(treatmentMetadata, annotated_ChEMBL,  by.x = "pubchem.ChEMBL.ID", by.y = "molecule.chembl.id", all.x = TRUE)

dir.create(dirname(file.path(OUTPUT$annotated_treatmentMetadata)), recursive = TRUE)
data.table::fwrite(final_annotated_treatmentMetadata, OUTPUT$annotated_treatmentMetadata, quote = FALSE, sep = "\t", row.names = FALSE)
