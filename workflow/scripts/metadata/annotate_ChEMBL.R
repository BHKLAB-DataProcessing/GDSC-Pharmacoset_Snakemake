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

chembl_dt <- data.table::fread(INPUT$annotated_CIDS, header = TRUE, sep = "\t")

annotated_ChEMBL_dt <- AnnotationGx::getChemblMechanism(as.character(chembl_dt[,pubchem.ChEMBL.ID]))

data.table::fwrite(
    annotated_ChEMBL_dt,
    OUTPUT$annotated_ChEMBL,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
)


