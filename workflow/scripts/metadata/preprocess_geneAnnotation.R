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

CMP_geneAnnotation <- fread(INPUT$geneAnnotation, header = TRUE)
names(CMP_geneAnnotation) <- paste0("CMP.", toupper(names(CMP_geneAnnotation)))


ensemblAnnot <- rtracklayer::import(INPUT$ensemblAnnotation)
ensembl_dt <- as.data.table(ensemblAnnot)

# This is only used for gencode
# [,c("gene_id", "transcript_id", "havana_transcript", "exon_id") 
#   := lapply(.SD, function(x) gsub("\\.\\d+$", "", x)), 
#     .SDcols = c("gene_id", "transcript_id", "havana_transcript", "exon_id")]

# remove columns that are not needed
cols <- c(
    "seqnames", "start", "end", "strand", 
    "gene_id", "gene_name", "gene_source")
ensembl_dt <- ensembl_dt[type == "gene", ..cols]
names(ensembl_dt) <- paste0("ENSEMBL.", toupper(names(ensembl_dt)))

mergedGeneAnnotation <- merge(CMP_geneAnnotation, ensembl_dt, by.x = "CMP.ENSEMBL_GENE_ID", by.y = "ENSEMBL.GENE_ID", all.x = TRUE)

# for each column in x, if al the values are NA, remove the column
mergedGeneAnnotation <- x[, .SD, .SDcols = colnames(x)[colSums(is.na(x)) < nrow(x)]]

# replace any "" with NA
mergedGeneAnnotation[mergedGeneAnnotation == ""] <- NA

# ENSEMBL.SEQNAMES : seqnames
# ENSEMBL.START : start
# ENSEMBL.END : end
# ENSEMBL.STRAND : strand
# CMP.ENSEMBL_GENE_ID : gene_id
# ENSEMBL.GENE_NAME : gene_name
# convert the above columns from left to right, and set those as the first columns
oldCols <- c("ENSEMBL.SEQNAMES", "ENSEMBL.START", "ENSEMBL.END", "ENSEMBL.STRAND", "CMP.ENSEMBL_GENE_ID", "ENSEMBL.GENE_NAME") 
newCols <- c("seqnames", "start", "end", "strand", "gene_id", "gene_name")
setnames(mergedGeneAnnotation, oldCols, newCols)
mergedGeneAnnotation <- mergedGeneAnnotation[, c(newCols, setdiff(colnames(mergedGeneAnnotation), newCols)), with = FALSE][order(gene_id)]

# write to OUTPUT$geneAnnotation
fwrite(mergedGeneAnnotation, OUTPUT$geneAnnotation, quote = TRUE, sep = "\t", na = "NA", row.names = FALSE)

