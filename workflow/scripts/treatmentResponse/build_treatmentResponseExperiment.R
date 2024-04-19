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

# 0.1 Startup 
# -----------
# load data.table, suppressPackageStartupMessages unless there is an error
suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
snakemake@source("../metadata/utils.R")

# 1.0 Read input data into data.table
# -----------------------------------
# input:
#     rawdata = "rawdata/treatmentResponse/release_{release}/{version}_public_raw_data.csv",
#     processed = "rawdata/treatmentResponse/release_{release}/{version}_fitted_dose_response.xlsx",
#     treatmentMetadata = procdata / metadata / "{version}_{release}_preprocessed_treatmentMetadata.tsv",
rawdata <- fread(INPUT$rawdata)
processed <- as.data.table(readxl::read_xlsx(INPUT$processed))
treatmentMetadata <- fread(INPUT$treatmentMetadata)
sampleMetadata <- fread(INPUT$sampleMetadata)

# 2.0 Subset data
# ---------------
# This line subsets for treatments that are in the treatmentMetadata and have a valid TAG

rawdata_s <- rawdata[, GDSC.sampleid := ..sampleMetadata[match(SANGER_MODEL_ID, CMP.model_id), GDSC.sampleid]]

rawdata_st <- rawdata_s[, GDSC.treatmentid := ..treatmentMetadata[match(DRUG_ID, GDSC.DRUG_ID), GDSC.treatmentid]]

subsetted_rawdata <- rawdata_st[(!is.na(DRUG_ID) & !is.na(GDSC.treatmentid)) | 
        grepl("^(L|R|A|N|B)", TAG)]

## ------------------- GDSC SPECIFIC STEP ------------------- ##
# In the essence of using the processed GDSC data as much as possible
# we will normalize the rawdata using the gdscIC50 package
# neg_control_TAGS <- subsetted_rawdata[grepl("^NC", TAG), unique(TAG)]
# if(length(neg_control_TAGS) > 1){
#     print(paste0("More than one negative control found: ", paste0(neg_control_TAGS, collapse = ", ")))
#     print(paste0("Setting all negative controls to: ", neg_control_TAGS[1]))
#     subsetted_rawdata[TAG %in% neg_control_TAGS, TAG := neg_control_TAGS[1]]
#     neg_control_tag <- neg_control_TAGS[1]
# }else if(length(neg_control_TAGS) == 0){
#     stop("No negative control found")
# }else{
#     cat(paste0("Negative control found: ", neg_control_TAGS))
#     neg_control_tag <- neg_control_TAGS[1]
# }

neg_control_tag <- ifelse(WILDCARDS$version == "GDSC1", "NC-0", "NC-1")

# OLD METHODOLOGY
# gdsc_sens[, Viability := {
#   control <- .SD[TAG == ..control.column, ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
#   background <- .SD[TAG.1 == "B", ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
#   Viability <- (INTENSITY - background) / (control - background) * 100
#   list(Viability = Viability)
# }, .(MASTER_CELL_ID, BARCODE)]


# rawdata_st[, Viability := {
#     control <- .SD[TAG == ..neg_control_tag, ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
#     background <- .SD[grepl("^B", TAG), ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
#     Viability <- (INTENSITY - background) / (control - background) * 100
#     list(Viability = Viability)
#     }, .(MASTER_CELL_ID, BARCODE)]


if (!require("gdscIC50", quietly = TRUE))
    pak::pkg_install("cancerrxgene/gdscIC50")


subsetted_rawdata <- gdscIC50::removeFailedDrugs(subsetted_rawdata)

subsetted_rawdata <- gdscIC50::removeMissingDrugs(subsetted_rawdata)

normData <- gdscIC50::normalizeData(subsetted_rawdata, trim=TRUE, neg_control=neg_control_tag, pos_control="B")

# the normalization removes the treatmentid and sampleid columns
normData <- merge(normData, unique(subsetted_rawdata[, .(MASTER_CELL_ID, GDSC.sampleid)]), by = "MASTER_CELL_ID", all.x = T)
normData <- merge(normData, unique(subsetted_rawdata[, .(DRUG_ID, GDSC.treatmentid)]), by.x = "DRUG_ID_lib", by.y = "DRUG_ID", all.x = T)


assay <- normData[
    !(is.na(GDSC.sampleid)), 
    .(GDSC.sampleid, GDSC.treatmentid, Dose = CONC, Viability = normalized_intensity, BARCODE)]


assay[!is.na(GDSC.treatmentid), .N, by=.(GDSC.sampleid, GDSC.treatmentid)][order(-N)]$GDSC.treatmentid |> unique() -> drug_exp_count

print(paste0("Merging treatmentMetadata with processed data"))
procData <- merge(processed, treatmentMetadata[, .(GDSC.DRUG_ID, GDSC.treatmentid)], by.x = "DRUG_ID",  by.y = "GDSC.DRUG_ID", all.x = T)
procData <- procData[, GDSC.sampleid := cleanCharacterStrings(CELL_LINE_NAME)][!is.na(GDSC.sampleid) & GDSC.sampleid %in%  sampleMetadata$GDSC.sampleid]


# 3.0 
# -------------------
print(paste0("Subsetting normData to not include missing values"))
subset_normData <- unique(assay[
  !is.na(Viability) & 
    !is.na(GDSC.treatmentid) & 
    !is.na(Dose) & 
    !is.na(GDSC.sampleid)])


# CoreGx Functions require that the main terminology is "treatmentid" and "sampleid"
# we will add these columns to the subset_normData

subset_normData[, 'treatmentid' := GDSC.treatmentid]
subset_normData[, 'sampleid' := GDSC.sampleid]


unique_treatments <- unique(subset_normData$GDSC.treatmentid)
unique_samples <- unique(subset_normData$GDSC.sampleid)

print(paste0("Number of unique treatments: ", length(unique_treatments)))
print(paste0("Number of unique samples: ", length(unique_samples)))

# subset to only first 30 treatments
subset_normData <- subset_normData[GDSC.treatmentid %in% unique_treatments[1:30]]

# subset to only first 100 samples
subset_normData <- subset_normData[GDSC.sampleid %in% unique_samples[1:100]]


print("Constructing tre")
TREDataMapper <- CoreGx::TREDataMapper(rawdata=subset_normData)
CoreGx::rowDataMap(TREDataMapper) <- list(
    id_columns = (c("treatmentid", "Dose")),
    mapped_columns = c())

CoreGx::colDataMap(TREDataMapper) <- list(
    id_columns = c("sampleid", "BARCODE"),
    mapped_columns = c())

CoreGx::assayMap(TREDataMapper) <- list(
    sensitivity = list(
        c("treatmentid", "BARCODE", "Dose", "sampleid"),
        c("Viability")))

gdsc_tre <- CoreGx::metaConstruct(TREDataMapper)

published_profiles <- procData[
    GDSC.sampleid %in% gdsc_tre$sensitivity$sampleid & 
    GDSC.treatmentid %in% gdsc_tre$sensitivity$treatmentid,
    .(sampleid = GDSC.sampleid, treatmentid = GDSC.treatmentid,  LN_IC50, AUC, RMSE, Z_SCORE)]

CoreGx::assay(gdsc_tre, "profiles_published") <- published_profiles

# ADD METADATA
# ------------ 
metadata <- list(
    data_source = snakemake@config$treatmentResponse[[WILDCARDS$release]][[WILDCARDS$version]],
    filename = list(
        description = "The filename used in the pipeline to create this object. Likely renamed from the original file",
        raw = basename(INPUT$rawdata),
        processed = basename(INPUT$processed)
    ),
    date = Sys.Date()
)
CoreGx::metadata(gdsc_tre) <- metadata
######
# 4.0 Save output
# ---------------
dir.create(dirname(OUTPUT$raw), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUTPUT$published_profiles), showWarnings = FALSE, recursive = TRUE)
fwrite(subset_normData, OUTPUT$raw, sep = "\t", quote = F)
fwrite(published_profiles, OUTPUT$published_profiles, sep = "\t", quote = F)
saveRDS(gdsc_tre, OUTPUT$tre)



