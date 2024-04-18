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

# cleanCharacterStrings function from utils.R
snakemake@source("utils.R")

library(data.table)

# 1.0 Read Input Data
# -------------------

treatmentMetadata <- fread(INPUT$treatmentMetadata)

names(treatmentMetadata) <- paste0("GDSC.", names(treatmentMetadata))

treatmentMetadata[, GDSC.treatmentid := cleanCharacterStrings(GDSC.DRUG_NAME)]



# get all the column names, and reorder so that GDSC.treatment is the first column and CCLE.treatment is the second column
treatmentMetadata <- 
    treatmentMetadata[, 
    c("GDSC.treatmentid", setdiff(names(treatmentMetadata), c("GDSC.treatmentid"))),
    with = FALSE]

# write treatmentMetadata to OUTPUT$treatmentMetadata csv file
fwrite(treatmentMetadata, OUTPUT$treatmentMetadata, sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = TRUE)
