# AUTHOR: Jermiah Joseph
# CREATED: 01-08-2024
# This script takes in the following files:
# 1. sampleMetadata: A data.table containing the sample metadata
# 2. geneAnnotation: A data.table containing the gene annotation
# 3. WES_zipped: A zipped file containing the WES CNV data
# The script then filters the WES CNV data to only include samples in the sampleMetadata
# and genes in the geneAnnotation. The script then creates a SummarizedExperiment for each
# column in the WES CNV data and saves the SummarizedExperiment to a file.

# NOTES:
# 1. Categorisation of Total Copy Number values
# source(https://depmap.sanger.ac.uk/documentation/datasets/copy-number/)
    # The total copy number values have been categorised (CNA Call) 
    # using the following calculation:

    # Val = round( 2 * 2^log2(C/Ploidy) )

    # if Val == 0: Category = 'Deletion'
    # if Val == 1: Category = 'Loss'
    # if Val == 2: Category = 'Neutral'
    # if Val == 3: Category = 'Gain'
    # if Val >= 4: Category = 'Amplification'

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

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))

    
# 0.1 Read in the input files
# ---------------------------
sampleMetadata <- fread(INPUT$sampleMetadata)
geneAnnot <- fread(INPUT$geneAnnotation)

WES_CNV_dir <- paste0(dirname(INPUT$WES_zipped), "/WES_CNV")
dir.create(WES_CNV_dir, recursive = TRUE, showWarnings = FALSE)
print(paste0("Unzipping ", INPUT$WES_zipped, " into ", WES_CNV_dir))
unzip(INPUT$WES_zipped, exdir = WES_CNV_dir)

inputFilesNames <- file.path(WES_CNV_dir, list.files(WES_CNV_dir))

# 0.3 Read in the raw CNV data
# -----------------------
# Read in large file "WES_pureCN_CNV_genes_20221213.csv" :
# script was written to handle both WES and WGS, only to find out that WGS 
# does not contain data for any samples in GDSC sample metadata
# rawdata/cnv/WES_pureCN_CNV_genes_20221213.csv

gene_files <- inputFilesNames[grepl("genes", inputFilesNames)]

# Get the largest file in the list of gene_files
file <- gene_files[which.max(file.size(gene_files))]


print(paste("Loading", file," "))
df <- data.table::fread(
    file, header = TRUE, showProgress = FALSE,
    sep = ",", stringsAsFactors = FALSE, nThread = THREADS)

print(sprintf("Loaded %s with %d rows and %d columns", file, nrow(df), ncol(df)))

print(sprintf("Subsetting %s to only samples in GDSC sample metadata", file))
df <- df[model_id %in% sampleMetadata[, CMP.model_id]]

print(sprintf("Subsetting %s to only genes in GDSC gene annotation", file))
# drop the symbol column 
df_subset <- merge(
    df[gene_id %in% geneAnnot$CMP.GENE_ID, !c("symbol"), with = FALSE], 
    data.table::as.data.table(geneAnnot)[,.(gene_name, CMP.GENE_ID)], 
    by.x = "gene_id", by.y = "CMP.GENE_ID", all.x = TRUE)

print(sprintf("data.table now has %d rows and %d columns", nrow(df_subset), ncol(df_subset)))
print(sprintf("Total number of model_ids: %d", uniqueN(df_subset[, model_id])))
print(sprintf("Total number of gene_ids: %d", uniqueN(df_subset[, gene_id])))

# 1. Build data structures for each datatype:
# -------------------------------------------
wes_gene_dt <- merge(df_subset,
    sampleMetadata[, .(GDSC.sampleid, CMP.model_id)],
    by.x = "model_id", by.y = "CMP.model_id", all.x = TRUE)

# ------------------------------
# print(paste("\n",capture.output(wes_gene_dt[, .N, by = source])))
#    source        N
# 1: Sanger 16884273
# 2:  Broad  1843646

cols <- c(
    "GDSC.sampleid", "gene_name",
    "total_copy_number", "cn_category", 
    "seg_mean", "gene_mean", "num_snps", "gatk_mean_log2_copy_ratio", "source")
dt <- wes_gene_dt[, ..cols, with = FALSE]

# subset assay_dt to only include rows with source == "Sanger" 
source_ <- "Sanger"
assay_dt <- dt[source == source_, !c("source"), with = FALSE]

# Generic metadata list to be used (and updated for each assay)
metadata <- list(
    annotation = "cnv", 
    data_source = snakemake@config$molecularProfiles$cnv,
    gene_annotation = snakemake@config$metadata$referenceGenome,
    class = "RangedSummarizedExperiment",
    filename = file,
    date = Sys.Date()
)


# for each column that isnt GDSC.sampleid or symbol, create a matrix with genes rows
# and samples as columns and set the rownames to the gene_id
# the `dcast` function to reshape the data from 
# long to wide format, with genes as rows and samples as columns.
# The resulting data frame `assay_dt.t` contains the `col` values 
# (i.e total_copy_number, cn_category, etc) for each gene-sample combination.
# If the `col` values are of type character, the function `first` 
# is used to aggregate the values, otherwise the mean is calculated.
assayNames <- cols[!cols %in% c("GDSC.sampleid", "gene_name", "source")]
rse_list <- BiocParallel::bplapply(
    assayNames,
    function(col){
        
        message(paste("Casting ", col))
        assay_dt.t <- dcast(
            unique(assay_dt[!is.na(gene_name), c("GDSC.sampleid", "gene_name", col), with = FALSE]),
            gene_name ~ GDSC.sampleid,
            value.var = col,
            fun.aggregate = if(is.character(assay_dt[[col]])) dplyr::first else mean
        )
        assay_dt.t  <- assay_dt.t[gene_name %in% geneAnnot$gene_name,]

        print(paste("Converting ", col, " to matrix"))
        mtx <- as.matrix(
            assay_dt.t[, !c("gene_name"), with = FALSE],
            rownames = assay_dt.t[["gene_name"]])
            
        # now NA rownames
        mtx <- mtx[!is.na(rownames(mtx)),]
        rowRanges <- GenomicRanges::makeGRangesFromDataFrame(
            geneAnnot[gene_name %in% rownames(mtx) & !is.na(gene_name),], 
            keep.extra.columns = TRUE
        )

        print(sprintf(
            "Matrix %s has %d rows and %d columns", col, nrow(mtx), ncol(mtx)))

        tmpMetadata <- c(
            metadata,
            list(
                datatype = col,
                numSamples = ncol(mtx),
                numGenes = nrow(mtx)
            )
        )

        colData <- data.frame(
            sampleid = colnames(mtx),
            batchid = rep(NA, ncol(mtx))
        )

        print(paste("Creating SummarizedExperiment for ", col))
        rse <- SummarizedExperiment::SummarizedExperiment(
            assays = list(exprs = mtx),
            colData = colData,
            rowRanges = rowRanges,
            metadata = tmpMetadata
        )
        
        print(paste("Writing ", col, " to ", OUTPUT[[col]]))
        write.table(
            mtx, 
            file = OUTPUT[[col]],
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)
        
        return(rse)
    },
    BPPARAM = BiocParallel::MulticoreParam(
        workers = THREADS, timeout=3600)
)
names(rse_list) <- paste0("cnv.", assayNames)

# Get numSamples in each matrix
numSamples <- max(sapply(rse_list, function(x) ncol(SummarizedExperiment::assay(x))))


# For each rse, add the metadata
for (i in seq_along(rse_list)){
    rse_list[[i]]@metadata <- metadata
}

# 3. Save Output
# # --------------
print("Saving Output Files")

# make output directory if it doesnt exist
dir.create(dirname(OUTPUT$rse_list), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(metadata, OUTPUT$metadata)

saveRDS(rse_list, file = OUTPUT$rse_list)

