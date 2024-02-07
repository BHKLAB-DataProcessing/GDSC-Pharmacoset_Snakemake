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

library(data.table)
suppressPackageStartupMessages(library(GenomicRanges))

# 0.2 read in metadata
# -------------------- 
sampleMetadata <- fread(INPUT$sampleMetadata)
geneAnnot <- fread(INPUT$geneAnnotation)


# 0.2 read gene fusions
# ---------------------
zipDir <- dirname(INPUT$gene_fusions)
unzipDir <- file.path(dirname(zipDir), "gene_fusions")
unzip(INPUT$gene_fusions, exdir = unzipDir)

print(paste("Loading", list.files(unzipDir), " "))
dt <- data.table::fread(file.path(unzipDir, list.files(unzipDir)))


cols <- c(
    "model_id", "model_name", "tissue",
    "chr_3prime", "chr_5prime", 
    "gene_id_3prime", "gene_id_5prime", 
    "gene_symbol_3prime", "gene_symbol_5prime")

# 1.0 subset gene fusions to only genes in GDSC sample metadata
# ----------------------------------------------------------------
fusion_dt <- dt[model_id %in% sampleMetadata[, CMP.model_id], ..cols]

fusion_dt <- fusion_dt[gene_id_3prime %in% geneAnnot$CMP.GENE_ID]
fusion_dt <- fusion_dt[gene_id_5prime %in% geneAnnot$CMP.GENE_ID]

fusion_dt$status <- "fusion"

dt_t <- data.table::dcast(
    fusion_dt[, .(gene_symbol_5prime, model_name, gene_symbol_3prime, status)], 
    gene_symbol_3prime + gene_symbol_5prime ~ model_name, 
    value.var = "status", 
    fun.aggregate = dplyr::first,
    fill = "wt",
    sep = "_"
)

# make a matrix from dt_t with a combination of 
# gene_symbol_3prime and gene_symbol_5prime as rownames and
# model_id as colnames
mtx <- as.matrix(
    dt_t[, -c("gene_symbol_3prime", "gene_symbol_5prime")],
    rownames = paste0(dt_t$gene_symbol_3prime, "_", dt_t$gene_symbol_5prime))

# match the colnames of mtx to the GDSC.SampleName of 
# sampleMetadata[, .(GDSC.sampleid, GDSC.Sample_Name)]
# and then set the colnames of mtx to the corresponding GDSC.sampleid
matched_names <- match(colnames(mtx), sampleMetadata$GDSC.Sample_Name)
colnames(mtx) <- sampleMetadata$GDSC.sampleid[matched_names]

metadata <- list(
    data_source = snakemake@config$molecularProfiles$fusion,
    filename = basename(INPUT$gene_fusions),
    annotation = "fusion",
    samples = ncol(mtx),
    date = Sys.Date(),
    sessionInfo = capture.output(sessionInfo())
)


dt_t[, rownames := paste0(gene_symbol_3prime, "_", gene_symbol_5prime)]
rowData <- unique(merge(
    dt_t[, c("gene_symbol_3prime", "gene_symbol_5prime", "rownames")],
    fusion_dt[, 
        c("chr_3prime", "chr_5prime", 
        "gene_id_3prime", "gene_id_5prime", 
        "gene_symbol_3prime", "gene_symbol_5prime")],
    by = c("gene_symbol_3prime", "gene_symbol_5prime")
))

print("Creating SummarizedExperiment object")

rowData <- data.frame(
    rowData,
    stringsAsFactors = FALSE,
    row.names = rowData$rownames)

outputFiles <- list(
    "assays" = mtx,
    "rowData" = rowData,
    "metadata" = metadata)

rse_list <- list(
    "fusion" = SummarizedExperiment::SummarizedExperiment(
        assays = list(fusions = mtx),
        colData = data.frame(
            sampleid = colnames(mtx),
            batchid = rep(NA, ncol(mtx))
        ),
        rowData = rowData,
        metadata = metadata))


# 0.3 save the output
print(paste("Saving output to", OUTPUT$rse_list))
dir.create(dirname(OUTPUT$rse_list), showWarnings = FALSE, recursive = TRUE)
saveRDS(rse_list, file = OUTPUT$rse_list)

write.table(
    mtx,
    file = OUTPUT$fusion,
    quote = FALSE,
    sep = "\t",
    row.names = TRUE)

jsonlite::write_json(metadata, OUTPUT$metadata)
