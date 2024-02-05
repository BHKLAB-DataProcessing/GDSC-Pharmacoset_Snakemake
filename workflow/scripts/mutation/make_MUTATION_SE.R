
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    # save.image()
    
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

# 2.0 Create GRanges object from mutGenesAnnot
# --------------------------------------------
# Unsure if this is needed, keeping for future reference. 
# GRanges_dt <- data.table::as.data.table(geneAnnot)
# GRanges_dt <- GRanges_dt[, 
#     .(CMP_gene_id, ensembl_gene_id, entrez_id, strand, seqnames, 
#     hgnc_id, refseq_id, uniprot_id, gene_type, source)]

# mutGenesAnnot <- merge(
#     mutGenesAnnot[, !c("strand")], 
#     GRanges_dt, 
#     by.x = "gene_id", by.y = "CMP_gene_id")

# mutGenesAnnot <- mutGenesAnnot[gene_id %in% mut_dt$gene_id]

# gr <- GenomicRanges::makeGRangesFromDataFrame(
#     df = mutGenesAnnot, keep.extra.columns=TRUE, na.rm=TRUE,
#     start.field = "chr_start", end.field = "chr_end", seqnames.field = "seqnames",)

# 3.0 Create Assays object from mut_dt
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
    mtx
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS))
names(matrices) <- assay_cols




# 4.0 Setup metadata for SummarizedExperiment object
# --------------------------------------------------
metadata <- list(
    data_source = snakemake@config$molecularProfiles$mutation$SUMMARY,
    filename.data = INPUT$all_mutations,
    filename.gene_metadata = INPUT$mutation_genes)

geneAnnot 

print("Creating SummarizedExperiment objects")

rse_list <- lapply(matrices, function(x) {
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

    # create SummarizedExperiment object
    rse <- SummarizedExperiment::SummarizedExperiment(
        assays = list(exprs = assay),
        rowRanges = gr,
        metadata = metadata,
        colData = colData)

})
names(rse_list) <- paste0("mut.", gsub("_mutation", "", names(matrices)))

print("Done creating SummarizedExperiment objects")
print("Writing output to disk")

# 5.0 Save Output
# ---------------
jsonlite::write_json(metadata, OUTPUT$metadata)

dir.create(dirname(OUTPUT$rse_list), recursive = TRUE, showWarnings = FALSE)
saveRDS(rse_list, file = OUTPUT$rse_list)

outputs_written <- lapply(names(rse_list), function(x) {
    assay <- SummarizedExperiment::assay(rse_list[[x]], "exprs")
    filename <- gsub("mut.", "", x)
    write.table(
        assay,
        file = OUTPUT[[filename]],
        quote = FALSE,
        sep = "\t",
        row.names = TRUE)
})