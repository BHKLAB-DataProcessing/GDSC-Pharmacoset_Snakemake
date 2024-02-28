# Script that creates a metadata file for the Data Nutrition Label (DNL) for the GDSC2 dataset

# Load libraries
library(data.table)

# Load GDSC2 data
gdsc2 <- readRDS("results/data/GDSC2_8.5_Pharmacoset.RDS") 

treatmentResponse <- slot(gdsc2, "treatmentResponse")

treatmentResponse_metadata <- treatmentResponse@metadata[-which(names(treatmentResponse@metadata) == "sessionInfo")]

molecularProfiles <- slot(gdsc2, "molecularProfiles")
molecularProfiles_metadata <- lapply(names(molecularProfiles), function(x){
    se <- molecularProfiles[[x]]
    metadata <- se@metadata
    metadata <- metadata[-which(names(metadata) == "sessionInfo")]
    return(metadata)
})


treatment <- data.table::as.data.table(slot(gdsc2, "treatment"))
sample <- data.table::as.data.table(slot(gdsc2, "sample"))

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
dir.create("results/metadata", recursive = TRUE, showWarnings = FALSE)

jsonlite::write_json(metadata_list, "results/metadata/GDSC2_metadata.json")
