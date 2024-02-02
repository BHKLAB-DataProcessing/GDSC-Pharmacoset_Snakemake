#' RULE: make_CNV_SE
#' AUTHOR: Jermiah Joseph
#' DATE: 01-15-2024
#' This script takes in the following files:
#' - INPUT$preprocessedCNV
#' and outputs the following files:
#' - OUTPUT$CNV_se
#' 
#' Libraries Used:
#' - data.table
#' - SummarizedExperiment
#' - log4r
#' - BiocParallel
#' - qs
#' - GenomicRanges

## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}

library(data.table)
suppressPackageStartupMessages(library(GenomicRanges))

# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(
      log4r::file_appender(LOGFILE, append = TRUE),
      log4r::console_appender()
    )
)


# make a function to easily log messages to the logger
info <- function(...) log4r::info(logger, paste0(..., collapse = " "))

info("Starting make_CNV_SE.R\n")

# 0.2 Read in the input files
# ---------------------------
info("Loading ", INPUT$preprocessedCNV, " ")
preproc <- qs::qread(INPUT$preprocessedCNV, nthreads = THREADS)

assays_ <- preproc$assays
rowRanges_ <- preproc$GRanges
metadata_ <- preproc$metadata

# Append to metadata as much information as possible 
metadata_$annotation <- "cnv"
metadata_$date_created <- Sys.time()
metadata_$sessionInfo <- capture.output(sessionInfo())

info("Matrices:\n", paste(capture.output(str(assays_)), collapse = "\n"))


# This code block creates a list of SummarizedExperiment objects, 
# one for each assay in the 'assays_' object.
# Iterate over each assay in the 'assays_' object
rse_list <- lapply(names(assays_), function(assay_name){
    # Get the assay data for the current assay
    assay <- assays_[[assay_name]]
    
    # Get the sample IDs and gene IDs for the assay
    sampleids <- colnames(assay)
    geneids <- rownames(assay)
    # Create the colData data frame, which contains information about each sample
    colData <- data.frame(
        sampleid = sampleids,
        batchid = rep(NA, length(sampleids)),
        row.names = sampleids
    )

    rowRanges <- rowRanges_[rowRanges_$symbol %in% geneids,]
    # Create a SummarizedExperiment object for the current assay
    SummarizedExperiment::SummarizedExperiment(
        assays = list(exprs = assays_[[assay_name]]),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata_
    )
})

# Assign names to the rse_list based on the assay names
names(rse_list) <- names(assays_)

rse_list <- list(
    "cnv" = rse_list
)

info("SummarizedExperiment:\n", paste(capture.output(rse_list),collapse = "\n"))

dir.create(dirname(OUTPUT$CNV_se), recursive = TRUE, showWarnings = FALSE)
qs::qsave(rse_list, file = OUTPUT$CNV_se)
