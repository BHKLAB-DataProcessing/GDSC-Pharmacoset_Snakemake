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
    save.image("preprocess_treatmentMetadata.RData")
}

# cleanCharacterStrings function from utils.R
snakemake@source("utils.R")

library(data.table)

# 1.0 Read Input Data
# -------------------

treatmentMetadata <- fread(INPUT$treatmentMetadata)

names(treatmentMetadata) <- paste0("GDSC.", names(treatmentMetadata))

treatmentMetadata[, GDSC.treatmentid := cleanCharacterStrings(GDSC.DRUG_NAME)]

# For some reason, the GDSC to CCLE treatment mapping file is organized so that the actual information starts at row 3
# and ends at row 25 and the columns are B to E.... need to do some hard coding here too lol
GDSC_to_CCLE_treatmentMappingFile <- readxl::read_excel(
    path = INPUT$GDSC_to_CCLE_treatmentMappingFile, 
    sheet = 1, range = "B10:E25", col_names = TRUE, na = "NA")
GDSC_to_CCLE_treatmentMappingFile <- 
    as.data.table(GDSC_to_CCLE_treatmentMappingFile)[, .(GDSC.DRUG_ID = `GDSC1000 drug ids`, CCLE.treatmentid = `CCLE name`)]



treatmentMetadata <- merge(treatmentMetadata, GDSC_to_CCLE_treatmentMappingFile, by = "GDSC.DRUG_ID", all.x = TRUE)

# get all the column names, and reorder so that GDSC.treatment is the first column and CCLE.treatment is the second column
treatmentMetadata <- 
    treatmentMetadata[, 
    c("GDSC.treatmentid", "CCLE.treatmentid", setdiff(names(treatmentMetadata), c("GDSC.treatmentid", "CCLE.treatmentid"))),
    with = FALSE]

# write treatmentMetadata to OUTPUT$treatmentMetadata csv file
fwrite(treatmentMetadata, OUTPUT$treatmentMetadata, sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = TRUE)
