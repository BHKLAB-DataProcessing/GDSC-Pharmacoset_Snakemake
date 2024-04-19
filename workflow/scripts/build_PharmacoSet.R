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
    # Saves the workspace to "resources/"build_PharmacoSet"
    file.path("resources", paste0(snakemake@rule, ".RData")) |> 
        save.image()
}else{
    # If the snakemake object does not exist, load the workspace
    file.path("resources", "build_PharmacoSet.RData") |>
        load()
}

library(PharmacoGx)
library(data.table)
###############################################################################
# Load INPUT
###############################################################################
message("Reading in SampleMetadata")
sampleMetadata <- data.table::fread(INPUT$annotated_sampleMetadata)

message("Reading in TreatmentMetadata")
treatmentMetadata <- data.table::fread(INPUT$annotated_treatmentMetadata)

message("Reading in TreatmentResponseExperiment")
tre <- readRDS(INPUT$treatmentResponseExperiment)

message("Reading in MultiAssayExperiment")
mae <- readRDS(INPUT$multiAssayExperiment)

###############################################################################
# Main Script
###############################################################################

## ------------------------------------------------------------------------- ##
# SANITY CHECKS

# Need to make sure that the sampleMetadata and the colnames of the summarized
# experiments are the same

allSamples <- sapply(names(mae), function(se) colnames(mae[[se]])) |>
    unlist() |>
    unique()

# Check if all the samples in the MultiAssayExperiment are in the sampleMetadata
stopifnot(all(allSamples %in% sampleMetadata$sampleid))

treSamples <- tre@colData$sampleid |> 
    unique() 
stopifnot(all(treSamples %in% sampleMetadata$sampleid))

allTreatments <- unique(tre@rowData$treatmentid)

# Check if all the treatments in the TreatmentResponseExperiment are in the treatmentMetadata
stopifnot(all(allTreatments %in% treatmentMetadata$treatmentid))

## ------------------------------------------------------------------------- ##

treatmentMetadata <- treatmentMetadata[!duplicated(treatmentid),]
treatment_df <-  data.frame(
    treatmentMetadata, 
    row.names = treatmentMetadata$treatmentid)

sampleMetadata <- sampleMetadata[!duplicated(sampleid),]
sample_df <- data.frame(
    sampleMetadata, 
    row.names = sampleMetadata$sampleid
)

name <- paste(WILDCARDS$version, WILDCARDS$release, sep = "v")
print(paste("Name:", name))

pset <- PharmacoGx::PharmacoSet2(
    name = name,
    treatment = treatment_df,
    sample = sample_df,
    molecularProfiles = mae,
    treatmentResponse = tre,
    curation = list(
        sample = sample_df,
        treatment = treatment_df,
        tissue = data.frame()
    ),
    datasetType = "sensitivity"
)



###############################################################################
# Save OUTPUT 
###############################################################################
# OUTPUT[[1]]

message("Saving pset to ", OUTPUT[[1]])
dir.create(dirname(OUTPUT[[1]]), showWarnings = FALSE, recursive = TRUE)

saveRDS(pset, OUTPUT[[1]])
