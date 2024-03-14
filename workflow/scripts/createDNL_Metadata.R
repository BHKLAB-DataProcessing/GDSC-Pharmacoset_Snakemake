# Script that creates a metadata file for the Data Nutrition Label (DNL) for the GDSC2 dataset

# Load libraries
library(data.table)

# Load GDSC2 data
pset<- readRDS("results/data/GDSC2_8.5_Pharmacoset.RDS") 

treatmentResponse <- slot(pset, "treatmentResponse")

treatmentResponse_metadata <- treatmentResponse@metadata[-which(names(treatmentResponse@metadata) == "sessionInfo")]

molecularProfiles <- slot(pset, "molecularProfiles")
molecularProfiles_metadata <- lapply(names(molecularProfiles), function(x){
    se <- molecularProfiles[[x]]
    metadata <- se@metadata
    metadata <- metadata[-which(names(metadata) == "sessionInfo")]
    return(metadata)
})


treatment <- data.table::as.data.table(slot(pset, "treatment"))
sample <- data.table::as.data.table(slot(pset, "sample"))

sampleMetadata <- {
    numSamples <- nrow(sample)
    sampleAnnotations <- names(sample)
    list(
        numSamples = numSamples,
        sampleAnnotations = sampleAnnotations
    )
}


treatmentMetadata <- {
    numTreatments <- nrow(treatment)
    treatmentAnnotations <- names(treatment)
    list(
        numTreatments = numTreatments,
        treatmentAnnotations = treatmentAnnotations
    )
}

metadata_list <- list(
    sampleMetadata = sampleMetadata,
    treatmentMetadata = treatmentMetadata,
    molecularProfiles = molecularProfiles_metadata,
    treatmentResponse = treatmentResponse_metadata
)


# metadata_list <- lapply(names(metadata_list), function(x){
#     # molecularProfiles[[x]]@metadata
#     # return metadata without sessionInfo
#     metadata <- molecularProfiles[[x]]@metadata
#     metadata <- metadata[-which(names(metadata) == "sessionInfo")]
#     return(metadata)

# })

# metadata_list |> jsonlite::toJSON(x=_, pretty = T)

# Save metadata
save_path <- paste0("results/metadata/", name(pset), "_metadata.json")
dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)

jsonlite::write_json(metadata_list, save_path, pretty = T)
