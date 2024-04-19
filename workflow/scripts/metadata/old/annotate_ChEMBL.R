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
message("Reading annotated ChEMBL file: ", INPUT$annotated_CIDS)
chembl_dt <- data.table::fread(INPUT$annotated_CIDS, header = TRUE, sep = "\t")


# split every value of pubchem.ChEMBL.ID column by "; " and take the first value
chembl_mechs <- BiocParallel::bplapply(chembl_dt[pubchem.ChEMBL.ID!="",]$pubchem.ChEMBL.ID, function(ID) {
    id <- strsplit(ID, "; ")[[1]][1]
    AnnotationGx::getChemblMechanism(id)
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = 8)
) |> data.table::rbindlist()



# if any column has a list in any of its entries, then collapse the list into a string with ;
list_columns <- sapply(chembl_mechs, is.list)
for (i in names(list_columns)[list_columns]) {
    chembl_mechs[[i]] <- sapply(chembl_mechs[[i]], function(x) x[1])
}
chembl_mechs <- chembl_mechs[!duplicated(molecule_chembl_id),] 

message("Writing annotated ChEMBL data to file: ", OUTPUT$annotated_ChEMBL)
data.table::fwrite(
    chembl_mechs,
    OUTPUT$annotated_ChEMBL,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
)


