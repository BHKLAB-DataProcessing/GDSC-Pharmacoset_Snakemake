#' RULE: build_PharmacoSet
#' AUTHOR: Jermiah Joseph
#' DATE: 01-15-2024
#' This script takes in the following files:
#' - INPUT$metadata
#' - INPUT$summarizedExperiments
#' and outputs the following files:
#' - OUTPUT$pset
#' 
#' Libraries Used:
#' - MultiAssayExperiment
#' - log4r
#' - BiocParallel
#' - qs
#' - SummarizedExperiment
#' - CoreGx
#' - PharmacoGx
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

suppressPackageStartupMessages(library(PharmacoGx))


# 0.2 Read in the metadata
# --------------------------

treatmentMetadata <- data.table::fread(INPUT$annotated_treatmentMetadata)
sampleMetadata <- data.table::fread(INPUT$annotated_sampleMetadata)
geneAnnotation  <- data.table::fread(INPUT$geneAnnotation)


# 0.3 Read in the summarized experiments
# --------------------------------------
# Read the summarized experiments
print(paste("Loading: ", INPUT$summarizedExperimentLists, sep = "\n\t"))
summarizedExperimentLists <- unlist(lapply(INPUT$summarizedExperimentLists, readRDS))

# 0.4 Read in treatmentResponseExperiment
# ----------------------------------------
print(paste("Loading: ", INPUT$treatmentResponseExperiment, sep = "\n\t"))
tre <- readRDS(INPUT$treatmentResponseExperiment)

# 1.0 Build MultiAssayExperiment
# ------------------------------
# Extract unique sample IDs from the summarized experiments

se_list <- sapply(summarizedExperimentLists, function(x){
    if(!all(colnames(x) %in% sampleMetadata$GDSC.sampleid)){
        print("Not all colnames of the summarized experiment are in the sample metadata")
        new_x <- x[, colnames(x) %in% sampleMetadata$GDSC.sampleid]

        # remove duplicated columns
        new_x <- new_x[, !duplicated(colnames(new_x))]
    
        return(new_x)
    }
    else{
        new_x <- x[, !duplicated(colnames(x))]
        return(new_x)
    }
})
stopifnot(all(sapply(se_list, function(x){
    # make sure colnames of each SE is in sampleMetadat$GDSC.sampleid
    all(colnames(x) %in% sampleMetadata$GDSC.sampleid)
})))
summarizedExperimentLists <- se_list
summarizedExperimentLists <- sapply(summarizedExperimentLists, function(x){
    x@colData <- DataFrame(
        sampleid = colnames(x),
        batchid = rep(NA, ncol(x)),
        row.names = colnames(x)
    )
    x
})
# Remove duplicate sample IDs
sample <- sampleMetadata[!duplicated(GDSC.sampleid), ]

# convert sample into a data frame with the rownames being the sample IDs
# and ordered by the sample IDs
sample <- data.frame(sample, row.names = sample$GDSC.sampleid)
sample <- sample[order(rownames(sample)), ]

print(sprintf("Total number of samples across all experiments: %d", nrow(sample)))


treatment <- treatmentMetadata[!duplicated(treatmentMetadata$GDSC.treatmentid), ]
treatment <- data.frame(treatment, row.names = treatment$GDSC.treatmentid)
treatment <- treatment[order(rownames(treatment)), ]
print(sprintf("Total number of treatments: %d", nrow(treatment)))


# Create a data frame for the column data, including sample IDs and batch IDs
colData <- data.frame(
    sampleid = sample$GDSC.sampleid,
    batchid = rep(NA, length(sample$GDSC.sampleid)),
    row.names = sample$GDSC.sampleid
)
print(sprintf("Column data has %d rows and %d columns", nrow(colData), ncol(colData)))

# Create an ExperimentList object from the filtered summarized experiments
ExpList <- MultiAssayExperiment::ExperimentList(summarizedExperimentLists)
print(paste("ExperimentList:", capture.output(show(ExpList)), sep = ""))

# Create a sample map for each experiment in the ExperimentList
sampleMapList <- lapply(summarizedExperimentLists, function(se){
    data.frame(
        primary = colnames(se),
        colname = colnames(se),
        stringsAsFactors = FALSE
    )
})
names(sampleMapList) <- names(ExpList)
print(paste("Sample map list:", capture.output(str(sampleMapList)), sep = ""))

# Convert the sample map list to a single sample map
mae_sampleMap <- MultiAssayExperiment::listToMap(sampleMapList)
print(paste("Sample map:\n", capture.output(str(mae_sampleMap)), sep = ""))

# Create a MultiAssayExperiment object with the ExperimentList, column data, and sample map
mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = ExpList,
    colData = colData,
    sampleMap = MultiAssayExperiment::listToMap(sampleMapList)
)
print(paste("MultiAssayExperiment:\n", capture.output(show(mae)), sep = ""))


name <- paste(WILDCARDS$version, WILDCARDS$release, sep = "v")
print(paste("Name:", name))

print("Creating PharmacoSet object")
pset <- PharmacoGx::PharmacoSet2(
    name = name,
    treatment = treatment,
    sample = sample,
    molecularProfiles = mae,
    treatmentResponse = tre,
    perturbation = list(),
    curation = list(sample = data.frame(), treatment = data.frame(), tissue = data.frame()),
    datasetType = "sensitivity"
)

print(paste("PharmacoSet object:", capture.output(show(pset)), sep = ""))

print(paste("Saving PharmacoSet object to", OUTPUT[[1]]))
dir.create(dirname(OUTPUT[[1]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(pset, file = OUTPUT[[1]])


