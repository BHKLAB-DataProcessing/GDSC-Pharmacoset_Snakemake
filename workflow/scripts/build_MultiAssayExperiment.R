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

    # Assuming that this script is named after _treatmentMetadatathe rule
    # Saves the workspace to "resources/"build_MultiAssayExperiment"
    file.path(paste0(snakemake@rule, ".RData")) |> 
        save.image()
}else{
    # If the snakemake object does not exist, load the workspace
    file.path("build_MultiAssayExperiment.RData") |>
        load()
}
library(MultiAssayExperiment)
library(data.table)
###############################################################################
# Load INPUT
###############################################################################
print(paste("Loading: ", INPUT$summarizedExperiment_lists, sep = "\n\t"))
summarizedExperimentLists <- unlist(lapply(INPUT$summarizedExperiment_lists, readRDS))
show(summarizedExperimentLists)

print(paste("Loading: ", INPUT$sampleMetadata, sep = "\n\t"))
sampleMetadata <- data.table::fread(INPUT$sampleMetadata)
sampleMetadata[, sampleid := GDSC.sampleid]
sampleMetadata
###############################################################################
# Main Script
###############################################################################
data.table::setcolorder(sampleMetadata, c("sampleid"))
str(sampleMetadata)

se_list <- sapply(summarizedExperimentLists, function(SE){
    if(!all(colnames(SE) %in% sampleMetadata$sampleid)){
        print("Not all colnames of the summarized experiment are in the sample metadata")
        subset_SE <- SE[, colnames(SE) %in% sampleMetadata$sampleid]

        # remove duplicated columns
        subset_SE <- subset_SE[, !duplicated(colnames(subset_SE))]
    
        return(subset_SE)
    }
    else{
        subset_SE <- SE[, !duplicated(colnames(SE))]
        return(subset_SE)
    }
})

stopifnot(all(sapply(se_list, function(x){
    # make sure colnames of each SE is in sampleMetadat$GDSC.sampleid
    all(colnames(x) %in% sampleMetadata$sampleid)
})))

sampleMetadata <- sampleMetadata[!duplicated(sampleMetadata$sampleid), ]
colData <- data.frame(
    sampleid = sampleMetadata$sampleid,
    batchid = rep(NA, nrow(sampleMetadata)),
    row.names = sampleMetadata$sampleid
)

se_list_wColData <- sapply(seq_along(se_list), function(i) {
    SE <- se_list[[i]]
    # if SE does not have a colData slot, create one
    if(ncol(slot(SE, "colData")) == 0){
        print(paste("SE :", names(se_list)[i], " does not have a colData slot. Adding one."))
        colData <- DataFrame(
            sampleid = colnames(SE),
            batchid = rep(NA, ncol(SE)),
            row.names = colnames(SE)
        )
        slot(SE, "colData") <- colData
    }
    return(SE)
})
names(se_list_wColData) <- names(se_list)

# Create a sample map for each experiment in the ExperimentList
sampleMapList <- lapply(se_list_wColData, function(se){
    data.frame(
        primary = colnames(se),
        colname = colnames(se),
        stringsAsFactors = FALSE
    )
})

ExpList <- MultiAssayExperiment::ExperimentList(se_list_wColData)

# Create a MultiAssayExperiment object with the ExperimentList, column data, and sample map
mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = ExpList,
    colData = colData,
    sampleMap = MultiAssayExperiment::listToMap(sampleMapList)
)
message("MultiAssayExperiment:\n", paste(capture.output(show(mae)), collapse = "\n\t"))

###############################################################################
# Save OUTPUT 
###############################################################################
# OUTPUT$mae

message("Saving mae to: ", OUTPUT$mae)
dir.create(dirname(OUTPUT$mae), showWarnings = FALSE, recursive = TRUE)
saveRDS(mae, OUTPUT$mae)
