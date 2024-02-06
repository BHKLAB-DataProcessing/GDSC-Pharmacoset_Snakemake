# 
# AUTHOR: Jermiah Joseph
# CREATED: 01-08-2024
# DESCRIPTION:
#   This script takes in a directory of .cel files and a sample annotation file
#   and performs RMA normalization on the .cel files. The normalized expression
#   matrix is then collapsed by gene symbol and the sample annotations are
#   merged with the expression matrix. The gene annotations are also merged with
#   the expression matrix. The final output is a list of 3 items:
#       1. expr_annot: a data.table containing the gene annotations
#       2. sample_annotations: a data.table containing the sample annotations
#       3. expr: a data.table containing the normalized expression matrix


# PACKAGE DEPENDENCIES:
#  - affy
#  - readr
#  - data.table
#  - BiocParallel
#  - WGCNA
#  - qs
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
# if (!requireNamespace("preprocessCore", quietly = TRUE)) {
#     BiocManager::install("preprocessCore", configure.args = "--disable-threading", ask = FALSE)
# }

# exit
celDirectory <- INPUT$CELfiles
sample_annotations <- INPUT$CEL_metadata

# ----------------------------- MICROARRAY DATA ----------------------------- ##
# Get a list of all the .cel files in the rawdata expression directory
# Read in the .cel files using the affy package and perform RMA normalization
# Extract the expression matrix from the normalized data
# Remove the ".cel" extension from the column names of the expression matrix

# fileList <- list.files(celDirectory, pattern = ".cel", full.names = TRUE)
message("Reading in the .cel files and performing RMA normalization")
affyFiles <- affy::read.affybatch(celDirectory[1:10L])

eSet <- affy::rma(affyFiles)
mtx <- affy::exprs(eSet)                            
colnames(mtx) <- sub(pattern = ".cel", "", colnames(mtx), fixed = TRUE) 

## --------------------------- SAMPLE ANNOTATIONS --------------------------- ##
# NOTE: These sample annotations are specific to the affymetrix data
# The rawdata files use Assay IDs specific to this technology and have provided annotations 
# to map these IDs back to the sample name 
sample_annotations <- data.table::fread(sample_annotations, header=TRUE, sep="\t")

# Get the common column names between the matrix and sample annotations
common_colnames <- intersect(colnames(mtx), sample_annotations$`Assay Name`)

# Subset the matrix and sample annotations to only include the common columns
mtx <- mtx[, colnames(mtx) %in% common_colnames]
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` %in% common_colnames, ]
sample_annotations <- sample_annotations[match(colnames(mtx), sample_annotations$`Assay Name`), ]

# Check if the column names of the matrix and sample annotations match
all.equal(colnames(mtx), sample_annotations$`Assay Name`)

## ----------------------------- GENE ANNOTATIONS ----------------------------- ##
library(hgu219.db)
k <- keys(hgu219.db, keytype="PROBEID")

columns <- c("SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID", "UNIPROT", "PMID") 

annotations <- BiocParallel::bplapply(columns, function(col){
        data <- data.table::as.data.table(
            select(hgu219.db, keys=k, columns=c(col), keytype="PROBEID"))
        data <- data[!duplicated(data$PROBEID), ]
        data
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS)
)

# merge all data.tables on "PROBEID"
gene_annotations <- Reduce(
    function(x, y) merge(x, y, by = "PROBEID", all = TRUE), annotations)

# # Get common gene names
common_genes <- intersect(rownames(mtx), gene_annotations$PROBEID)

# # Subset and match both matrices
mtx <- mtx[rownames(mtx) %in% common_genes, ]
gene_annotations <- gene_annotations[gene_annotations$PROBEID %in% common_genes, ]
gene_annotations <- gene_annotations[match(rownames(mtx), gene_annotations$PROBEID), ]
all.equal(rownames(mtx), gene_annotations$PROBEID)  # Check if matching worked


mtx_collapsed <- WGCNA::collapseRows(datET = mtx, rowGroup = gene_annotations$SYMBOL, rowID = rownames(mtx))$datETcollapsed
colnames(mtx_collapsed) <- sample_annotations$`Characteristics[cell line]`

expr <- data.table::as.data.table(mtx_collapsed, keep.rownames = "gene.symbol")

# subset gene_annotations on SYMBOL with the values in expr$gene.symbol
# only keep the rows in gene_annotations that have a match in expr

expr_annot <- gene_annotations[gene_annotations$SYMBOL %in% expr$gene.symbol, .(SYMBOL, GENENAME, ENSEMBL, ENTREZID, UNIPROT, PMID)]


# ----------------------------- OUTPUT ----------------------------- ##


# combined expr_annot and sample_annotations and mtx into a list of 3 items
outputList <- list(expr_annot, sample_annotations, expr)
names(outputList) <- c("expr_annot", "sample_annotations", "expr")


# make output directory if it doesnt exist
dir.create(dirname(OUTPUT[[1]]), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputList, OUTPUT[[1]], nthreads = THREADS)
