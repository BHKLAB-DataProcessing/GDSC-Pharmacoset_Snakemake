
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)

}

library(data.table)
suppressPackageStartupMessages(library(GenomicRanges))

# 0.2 read in metadata
# -------------------- 
sampleMetadata <- fread(INPUT$sampleMetadata)
geneAnnot <- fread(INPUT$geneAnnotation)

# 0.3 read in mutation data
# -------------------------
dir <- paste0(dirname(INPUT$all_mutations), "/all")
allDir <- paste0(dirname(INPUT$all_mutations), "/all")
dir.create(allDir, recursive = TRUE, showWarnings = FALSE)
unzip(INPUT$all_mutations, exdir = allDir)

# 0.4 read mutation gene metadata
# -------------------------------
print(paste("Loading", file.path(allDir, list.files(allDir))))
mut_data <- data.table::fread(
    file.path(allDir, list.files(allDir)), 
    header = TRUE, sep = ",", showProgress = F)

# 0.5 read mutation gene metadata
# -------------------------------
print(paste("Loading", INPUT$mutation_genes))
mutGenesAnnot <- data.table::fread(INPUT$mutation_genes, header = TRUE, sep = ",")

# 1.0 Subset mutation data to only include samples from GDSC metadata
# -------------------------------------------------------------------
print("Subsetting mutation data to only include samples from GDSC metadata")
data.table::setkey(mut_data, model_id)
mut_dt <- unique(merge(
    sampleMetadata[, .(GDSC.sampleid, CMP.model_id)], 
    mut_data, 
    by.x = "CMP.model_id", by.y = "model_id"))


# Setup metadata for SummarizedExperiment object
# --------------------------------------------------
metadata <- list(
    annotation = "mutation",
    class = "RangedSummarizedExperiment",
    data_source = snakemake@config$molecularProfiles$mutation$all_mutations,
    filename = basename(INPUT$all_mutations),
    date_created = Sys.Date()
)


# Create Assays object from mut_dt
# -------------------------------------------------
# isolate only columns to use for assay
assay_cols <- c("protein_mutation", "rna_mutation",
    "cdna_mutation", "vaf", "effect")
cols <- c("GDSC.sampleid", "gene_symbol", assay_cols)

# use only sanger data
assay_dt <- mut_dt[source == "Sanger", ..cols]

matrices <- BiocParallel::bplapply(assay_cols, function(x){
    print(sprintf("Subsetting df for %s", x))
    cols_ <- c("GDSC.sampleid", "gene_symbol", x)
    assay_mtx <- assay_dt[, ..cols_]

    # replace "-" in the x col with "wt"
    assay_mtx[, (x) := lapply(.SD, function(x) gsub("-", "wt", x)), .SDcols = x]

    print(sprintf("casting df for %s", x))
    # dcast assay_mtx to wide format
    assay_mtx <- data.table::dcast(
        unique(assay_mtx), gene_symbol ~ GDSC.sampleid, 
        value.var = x, fun.aggregate = dplyr::first, fill = "wt")

    print(sprintf("Converting %s into matrix", x))
    # create a matrix with gene_symbol as rownames and GDSC.sampleid as colnames
    mtx <- as.matrix(assay_mtx[, -c("gene_symbol"), with = FALSE], 
        rownames = assay_mtx[["gene_symbol"]])
    
    print(sprintf(
        "Matrix %s has %d rows and %d columns", x, nrow(mtx), ncol(mtx)))

    return(mtx)
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS))
names(matrices) <- assay_cols


print("Creating SummarizedExperiment objects")
rse_list <- lapply(names(matrices), function(matrix_name){

    x <- matrices[[matrix_name]]
    # subset geneAnnot
    assay <- x[rownames(x) %in% geneAnnot$gene_name,]
    geneAnnot_sub <- unique(geneAnnot[gene_name %in% rownames(assay)], by = "gene_name")

    dim(assay)
    dim(geneAnnot_sub)
    # create GRanges object
    gr <- GenomicRanges::makeGRangesFromDataFrame(
        df = geneAnnot_sub, keep.extra.columns=TRUE, na.rm=TRUE,
        start.field = "start", end.field = "end", seqnames.field = "seqnames",)
    
    colData <- data.frame(
        sampleid = colnames(assay),
        batchid = rep(NA, ncol(assay))
    )

    metadata$datatype = matrix_name
    metadata$numSamples = ncol(assay)
    metadata$numGenes = nrow(assay)
    metadata$gene_annotation = INPUT$geneAnnotation

    rse <- SummarizedExperiment::SummarizedExperiment(
        assays = list(exprs = assay),
        rowRanges = gr,
        colData = colData,
        metadata = metadata
    )
})

names(rse_list) <- paste0("mut.", gsub("_mutation", "", names(matrices)))

print("Done creating SummarizedExperiment objects")
print("Writing output to disk")

# 5.0 Save Output
# ---------------
message("Saving metadata to: ", OUTPUT$metadata)
jsonlite::write_json(metadata, OUTPUT$metadata)

message("Saving RSE list to: ", OUTPUT$rse_list)
dir.create(dirname(OUTPUT$rse_list), recursive = TRUE, showWarnings = FALSE)
saveRDS(rse_list, file = OUTPUT$rse_list)


# write assay matrices to disk
message("Writing assay matrices to disk")
outputs_written <- lapply(names(rse_list), function(x) {
    assay <- SummarizedExperiment::assay(rse_list[[x]], "exprs")
    filename <- gsub("mut.", "", x)
    message("Writing ", filename, " to ", OUTPUT[[filename]])
    write.table(
        assay,
        file = OUTPUT[[filename]],
        quote = FALSE,
        sep = "\t",
        row.names = TRUE
    )
})
