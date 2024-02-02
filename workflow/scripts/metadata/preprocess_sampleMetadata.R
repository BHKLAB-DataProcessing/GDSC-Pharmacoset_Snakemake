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

# cleanCharacterStrings function from utils.R

snakemake@source("utils.R")

library(data.table)

# 1.0 Read Input Data
# -------------------
cat("Reading sample metadata file: ", INPUT$sampleMetadata, sep="\n\t")
raw_sampleMetadata <- readxl::read_excel(INPUT$sampleMetadata, sheet = 1, col_names = TRUE, na = "NA")
raw_sampleMetadata <- as.data.table(raw_sampleMetadata)

# read GDSC to CCLE sample mapping file
# note: the excel file is for some reason organized so that the actual information starts at row 9
# and ends at row 398 and the columns are C to E.... need to do some hard coding here lol 
cat("Reading GDSC to CCLE sample mapping file: ", INPUT$GDSC_to_CCLE_sampleMappingFile, sep="\n\t")
GDSC_to_CCLE_sampleMappingFile <- readxl::read_excel(
    path = INPUT$GDSC_to_CCLE_sampleMappingFile,  sheet = 1, range = "C9:E398",
    col_names = TRUE, na = "NA")

GDSC_to_CCLE_sampleMappingFile <- as.data.table(GDSC_to_CCLE_sampleMappingFile
    )[, .(GDSC.COSMIC_identifier = `GDSC1000 cosmic id`, CCLE.sampleid = `CCLE name`, GDSC.Sample_Name = `GDSC1000 name`)]
# ONLY 390 SAMPLES ARE MAPPED USING THIS FILE

# CMP SAMPLE ANNOTAITON
CMP_sampleAnnotation <- fread(INPUT$CMP_sampleAnnotation, header = TRUE)
# subset the CMP_sampleAnnotation to only the samples in the sampleMetadata
# and columns of interest
cols <- c(
    'model_id', 'sample_id', 'model_name', 'tissue', 'cancer_type', 
    'cancer_type_ncit_id', 'sample_site', 'COSMIC_ID', 'BROAD_ID', 'CCLE_ID',
    'gender', 'ethnicity', 'age_at_sampling', "RRID"
)
CMP_sampleAnnotation <- CMP_sampleAnnotation[
    model_name %in% raw_sampleMetadata$`Sample Name`,
    ..cols]
# rename GDSC.sampleid = sample_id, GDSC.Sample_Name = model_name, CCLE.sampleid = CCLE_ID
names(CMP_sampleAnnotation) <- ifelse(
    !grepl("^GDSC", names(CMP_sampleAnnotation)), 
    paste0("GDSC.", names(CMP_sampleAnnotation)), 
    gsub("^GDSC", "GDSC.", names(CMP_sampleAnnotation)))
setnames(CMP_sampleAnnotation, "GDSC.sample_id", "CMP.sampleid")
setnames(CMP_sampleAnnotation, "GDSC.model_id", "CMP.model_id")
setnames(CMP_sampleAnnotation, "GDSC.model_name", "GDSC.sampleid")
setnames(CMP_sampleAnnotation, "GDSC.CCLE_ID", "CCLE.sampleid")
CMP_sampleAnnotation[, GDSC.sampleid := cleanCharacterStrings(GDSC.sampleid)]

# 2.0 Preprocess Sample Metadata
# ------------------------------

cat("Cleaning sample metadata file", sep="\n\t")
sampleMetadata <- copy(raw_sampleMetadata)
# clean column names using data.table
# remove \r and \n from column names
setnames(sampleMetadata, gsub("[\r\n]", "", names(sampleMetadata)))
setnames(sampleMetadata, gsub(" ", "_", names(sampleMetadata)))
# drop the row where  GDSC.Sample_Name == "TOTAL:""
sampleMetadata <- sampleMetadata[Sample_Name != "TOTAL:", ]

# Renames the columns of the sampleMetadata dataframe by adding a prefix "GDSC." 
# to the column names that do not already start with "GDSC."
# The grepl() function is used to check if the column names start with "GDSC".
# The gsub() function is used to replace the "GDSC" prefix with "GDSC." if it already has.
names(sampleMetadata) <- ifelse(
    !grepl("^GDSC", names(sampleMetadata)), 
    paste0("GDSC.", names(sampleMetadata)), 
    gsub("^GDSC", "GDSC.", names(sampleMetadata)))

# clean the sample names in the sampleMetadata dataframe
sampleMetadata[, GDSC.sampleid := cleanCharacterStrings(GDSC.Sample_Name)]


# NOTE: NO LONGER MERGING missingSamples INTO THE SAMPLEMETADATA. 
# leaving here for reference
# This step is for mapping to CCLE:
# there are two tables, GDSC_to_CCLE_sampleMappingFile and CMP_sampleAnnotation
# We want to combine the two tables where the CMP_sampleAnnotation is the main table 
# and the GDSC_to_CCLE_sampleMappingFile is the secondary table
# We only want the rows from the secondary table where the GDSC.COSMIC_identifier is NOT in the main table's 'COSMIC_ID' column 
missingSamples <- GDSC_to_CCLE_sampleMappingFile[!GDSC.COSMIC_identifier %in% CMP_sampleAnnotation$GDSC.COSMIC_ID, ] 


# Now we want to merge the CMP_sampleAnnotation with the sampleMetadata dataframe
# need to convert the GDSC.COSMIC_identifier to double
sampleMetadata[, GDSC.COSMIC_ID := as.numeric(GDSC.COSMIC_identifier)][, GDSC.COSMIC_identifier := NULL]
CMP_sampleAnnotation[, GDSC.COSMIC_ID := as.numeric(GDSC.COSMIC_ID)]
sampleMetadata <- merge(sampleMetadata, CMP_sampleAnnotation[, !c("GDSC.sampleid")], by = "GDSC.COSMIC_ID", all.x = TRUE)

# get all the column names, and reorder so that GDSC.sampleid is the first column and CCLE.sampleid is the second column

firstcols <- c("GDSC.sampleid", "GDSC.Sample_Name", "CCLE.sampleid", "GDSC.BROAD_ID", "GDSC.RRID", "CMP.model_id", "CMP.sampleid")
columnNames <- c(firstcols, names(sampleMetadata)[!names(sampleMetadata) %in% firstcols])
sampleMetadata <- sampleMetadata[, columnNames, with = FALSE]

cat("Writing sample metadata file: ", OUTPUT$sampleMetadata, sep="\n\t")
fwrite(sampleMetadata, OUTPUT$sampleMetadata, sep = "\t", quote = FALSE, na = "NA", row.names = FALSE, col.names = TRUE)


