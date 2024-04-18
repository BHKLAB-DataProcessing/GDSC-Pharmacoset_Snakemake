## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
# This snippet is run at the beginning of a snakemake run to setup the env
# Helps to load the workspace if the script is run independently or debugging
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads

    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(
            file = snakemake@log[[1]], 
            append = FALSE, 
            type = c("output", "message"), 
            split = TRUE
    )

    # Assuming that this script is named after the rule
    # Saves the workspace to "resources/"annotate_treatmentMetadata"
    # file.path("resources", paste0(snakemake@rule, ".RData")) |> 
        # save.image()
}else{
    # If the snakemake object does not exist, load the workspace
    # file.path("resources", "annotate_treatmentMetadata.RData") |>
        # load()
}


###############################################################################
# Load INPUT
###############################################################################
annotatedCIDs <- data.table::fread(INPUT$annotated_CIDS)
message("annotatedCIDs: ", paste0(capture.output(str(annotatedCIDs)), collapse = "\n"))



###############################################################################
# Main Script
###############################################################################

unichem_sources <- AnnotationGx::getUnichemSources(T)
data.table::setkey(unichem_sources, Name)

sources_of_interest <- c("chembl", "drugbank", "chebi", "phamgkb", "lincs", "clinicaltrials", "nih_ncc", "fdasrs", "pharmgkb", "rxnorm")

sourceID <- unichem_sources[Name == "pubchem", SourceID]

message("\n\nAnnotating with unichem...")
unichem_mappings <- lapply(annotatedCIDs$pubchem.CID, function(x){
  tryCatch({
    result <- AnnotationGx::queryUnichemCompound(type = "sourceID", compound = x, sourceID = sourceID)

    subset <- result$External_Mappings[Name %in% sources_of_interest, .(compoundID, Name)]
    # make Name the column names and the values the compoundID 
    subset$cid <- x
    dcast(subset, cid ~ Name, value.var = "compoundID", fun.aggregate = list)
  }, error = function(e) NULL)
  } 
  ) |> data.table::rbindlist(fill = T)
show(unichem_mappings)

# for each column, if its a list then make it a string with a comma separator
for(col in names(unichem_mappings)){
  if(is.list(unichem_mappings[[col]])){
    unichem_mappings[[col]] <- sapply(unichem_mappings[[col]], function(x) paste(x, collapse = ","))
  }
}

# Rename columns like drugbank to unichem.DrugBank etc, using the unichem_sources "NameLabel" column 
names(unichem_mappings) <- paste(
    "unichem", 
    unichem_sources[names(unichem_mappings), gsub(" ", "_", NameLabel)], 
    sep = "."
    )

all_annotated_treatmentMetadata <- merge(
    treatment_annotations, 
    unichem_mappings, 
    by.x = "pubchem.CID", 
    by.y = "unichem.NA",        # dw about name being NA itll get removed  by this merge
    all.x = T
)

## ------------------------------------------------------------------------- ##


annotated_treatmentMetadata <- copy(all_annotated_treatmentMetadata)
message("\n\nAnnotating with ChEMBL using Unichem-obtained ChEMBL IDs")
chembl_mechanisms_dt <- AnnotationGx::getChemblMechanism(
    annotated_treatmentMetadata$unichem.ChEMBL
)

chembl_cols_of_interest <- c(
        "molecule_chembl_id",  "parent_molecule_chembl_id", "target_chembl_id", "record_id", 
        "mechanism_of_action", "mechanism_comment", "action_type"
    )

annotated_treatmentMetadata <- merge(
    annotated_treatmentMetadata, 
    chembl_mechanisms_dt[, ..chembl_cols_of_interest], 
    by.x = "unichem.ChEMBL",
    by.y = "molecule_chembl_id", 
    all.x = TRUE
    )

data.table::setnames(
    annotated_treatmentMetadata, 
    chembl_cols_of_interest, 
    paste0("chembl.", chembl_cols_of_interest), 
    skip_absent = TRUE)

annotated_treatmentMetadata <- annotated_treatmentMetadata[!duplicated(pubchem.CID),]


final_treatmentMetadata <- merge(
    treatmentMetadata, 
    all_annotated_treatmentMetadata, 
    by.x = "GDSC.treatmentid", 
    by.y = "pubchem.name", 
    all.x = TRUE

)



###############################################################################
# Save OUTPUT 
###############################################################################

message("\n\nWriting out cleaned data to ", OUTPUT$annotated_treatmentMetadata)

data.table::fwrite(
    final_treatmentMetadata, 
    file = OUTPUT$annotated_treatmentMetadata,
    quote = FALSE, 
    sep = "\t"
)

