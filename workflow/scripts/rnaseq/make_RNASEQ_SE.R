# RULE: preprocess_RNASEQ
# AUTHOR: Jermiah Joseph
# DATE: 01-15-2024

## This script processes RNA-seq data for the GDSC-Pharmacoset project.
## It performs the following steps:
## 1. Parse the Snakemake object to retrieve input, output, wildcards, and thread information.
## 2. Load required packages: data.table, GenomicRanges.
## 3. Read in metadata and gene annotation files.
## 4. Unzip the processed RNA-seq data.
## 5. Read in the RNA-seq data files and add a source column.
## 6. Create metadata from configuration information and RNA-seq data.
## 7. Subset the RNA-seq data to include only samples from GDSC metadata.
##         - data is taken from cell model passports, which has a lot more samples than GDSC
## 8. Replace cell model passport (CMP) model IDs with sample IDs.
## 9. Subset the data for only genes in GDSC gene annotation.
## 10. Convert the data to a matrix.
## 11. Create genomic ranges from the gene annotation.
## 12. Create a SummarizedExperiment object for each dataset type (fpkm, read_count, tpm).
## 13. Save the SummarizedExperiment objects to output files.
## 14. Save the metadata to a JSON file.

# Required packages: data.table, GenomicRanges
## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)
    
}

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))



# 0.2 read in metadata
# --------------------
sample <- fread(INPUT$sampleMetadata)
geneAnnot <- fread(INPUT$geneAnnotation)
# 0.3 read in rnaseq data
# -----------------------
allDir <- paste0(dirname(INPUT$processed), "/all")
dir.create(allDir, recursive = TRUE, showWarnings = FALSE)
print(paste0("Unzipping ", INPUT$processed, " into ", allDir))
unzip(INPUT$processed, exdir = allDir)

# list.files(allDir)
# [1] "rnaseq_all_data_20220624.csv"   "rnaseq_fpkm_20220624.csv"      
# [3] "rnaseq_read_count_20220624.csv" "rnaseq_tpm_20220624.csv" 

# file names without date or extension
parsedNames <- lapply(list.files(allDir), 
    function(x) gsub("^rnaseq_|_\\d{8}\\.csv", "", x)
)
# for each csv in allDir, read in the csv into a data.table and add source col
rnaseq_data <- lapply(file.path(allDir, list.files(allDir)), function(file){
    # if basename  of file has "all" then return NULL
    if(grepl("all", basename(file))){
        print(paste0("Skipping ", file))
        return(NULL)
    }

    print(paste0("Reading in ", file))
    data <- data.table::fread(file, header = TRUE, sep = ",", showProgress = F)
    data$file <- file
    data
})
names(rnaseq_data) <- parsedNames

# NOTE: we are not going to use the "all" data due to the resources required
dataset_types <- c("fpkm", "read_count", "tpm")

# Create metadata from config information and rnaseq_data
metadata <- lapply(dataset_types, function(x) {
    list(
        data_source = snakemake@config$molecularProfiles$rnaseq$processed,
        filename = basename(unique(rnaseq_data[[x]][["file"]])),
        annotation = "rnaseq",
        datatype = x,
        gene_annotation = lapply(snakemake@config$metadata$referenceGenome, as.character),
        date = Sys.Date(),
        sessionInfo = capture.output(sessionInfo())    
)})
names(metadata) <- dataset_types

# 1.0 Subset rnaseq data to only include samples from GDSC metadata
# -----------------------------------------------------------------

print("Subsetting rnaseq data to only include samples from GDSC metadata")


rse_list <- BiocParallel::bplapply(
    dataset_types, 
    function(x){ 
        print(sprintf("Subsetting %s", x))

        datatype_dt <- 
            rnaseq_data[[x]][
                -(1:4),
                colnames(rnaseq_data[[x]]) %in% c("model_id", sample$CMP.model_id), 
                with = FALSE]

        # get the rows of the sample Metadata that are in the dataset
        sample_subset <- 
            sample[CMP.model_id %in% colnames(datatype_dt), .(GDSC.sampleid, CMP.model_id)]
        
        # replace the model id (i.e SIDM00003) with sampleid (i.e M14)
        data.table::setnames(
            x = datatype_dt,
            old = c("model_id", sample_subset$CMP.model_id), 
            new = c("gene_id", sample_subset$GDSC.sampleid))

        print(sprintf(
            "Subsetting %s for only genes in GDSC gene annotation", x))

        matched_genes <- 
            geneAnnot[CMP.GENE_ID %in% datatype_dt$gene_id & !is.na(gene_name), .(CMP.GENE_ID, gene_name)]
        
        # replace the gene_id with gene_name using matched_genes
        dt <- merge(datatype_dt, matched_genes, by.x = "gene_id", by.y = "CMP.GENE_ID", all = FALSE)

        # set column order alphabetically and remove gene_id [-1]
        ordered_cols <- c("gene_name", sort(colnames(datatype_dt)[-1]))
        dt <- dt[, ..ordered_cols]

        print(sprintf("Converting %s to matrix", x))
        mtx <- as.matrix(
            dt[, !c("gene_name"), with = FALSE],
            rownames = dt[["gene_name"]]
        )

        # Convert all values in mtx from character to numeric
        mtx <- apply(mtx, 2, as.numeric)

        print(sprintf(
            "Matrix %s has %d rows and %d columns", x, nrow(mtx), ncol(mtx)))

        rowRanges <- GenomicRanges::makeGRangesFromDataFrame(
            geneAnnot[CMP.GENE_ID %in% matched_genes$CMP.GENE_ID], 
            keep.extra.columns = TRUE)

        metadata[[x]]$samples <- ncol(mtx)
        metadata[[x]]$genes <- nrow(mtx)

        rse <- SummarizedExperiment::SummarizedExperiment(
            assays = list(exprs = mtx),
            rowRanges = rowRanges,
            metadata = metadata[[x]])

        return(rse)
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS)
)
names(rse_list) <- paste0("rnaseq.", dataset_types)


print(paste("Saving Output Files to", OUTPUT$rse_list))
saveRDS(rse_list, file = OUTPUT$rse_list)

result <- lapply(dataset_types, function(x){
    print(paste0("Writing ", x, " to ", OUTPUT[[x]]))
    dir.create(dirname(OUTPUT[[x]]), recursive = TRUE, showWarnings = FALSE)

    write.table(
        x = SummarizedExperiment::assay(rse_list[[paste0("rnaseq.", x)]]), 
        file = OUTPUT[[x]], 
        quote = FALSE, 
        sep = "\t", 
        row.names = TRUE)
})

print(paste("Saving metadata to", OUTPUT$metadata))
jsonlite::write_json(metadata, OUTPUT$metadata)

