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
load("combine_annotated_treatmentData.RData")
library(data.table)
annotated_ChEMBL <- data.table::fread(INPUT$annotated_ChEMBL, header = TRUE, sep = "\t")
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata, header = TRUE, sep = "\t")

# names(annotated_ChEMBL)
#  [1] "action_type"               "binding_site_comment"     
#  [3] "direct_interaction"        "disease_efficacy"         
#  [5] "max_phase"                 "mec_id"                   
#  [7] "mechanism_comment"         "mechanism_of_action"      
#  [9] "mechanism_refs"            "molecular_mechanism"      
# [11] "molecule_chembl_id"        "parent_molecule_chembl_id"
# [13] "record_id"                 "selectivity_comment"      
# [15] "site_id"                   "target_chembl_id"         
# [17] "variant_sequence"  

name_map <- list(
    "molecule_chembl_id" = "molecule.chembl.id",
    "action_type" = "chembl.ActionType",
    "mechanism_of_action" = "chembl.MechanismOfAction",
    "molecular_mechanism" = "chembl.MolecularMechanism",
    "mechanism_comment" = "chembl.MechanismComment",
    "parent_molecule_chembl_id" = "chembl.ParentMoleculeChEMBL.ID",
    "target_chembl_id" = "chembl.TargetChEMBL.ID"
)
fields <- names(name_map)
annotated_ChEMBL <- annotated_ChEMBL[, ..fields]
data.table::setnames(annotated_ChEMBL, old = names(annotated_ChEMBL), new = unlist(name_map))


# annotations_combined_DT <- merge(annotated_DT, annotated_ChEMBL, by.y = "molecule.chembl.id", by.x = "pubchem.ChEMBL.ID", all = TRUE)
# merge the annotations with the treatmentMetadata
final_annotated_treatmentMetadata <- merge(treatmentMetadata, annotated_ChEMBL,  by.x = "pubchem.ChEMBL.ID", by.y = "molecule.chembl.id", all.x = TRUE)

dir.create(dirname(file.path(OUTPUT$annotated_treatmentMetadata)), recursive = TRUE)
data.table::fwrite(final_annotated_treatmentMetadata, OUTPUT$annotated_treatmentMetadata, quote = FALSE, sep = "\t", row.names = FALSE)
