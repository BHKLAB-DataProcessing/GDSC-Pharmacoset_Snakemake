## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads

    # setup logger if log file is provided
    ifelse(length(snakemake@log)>0, 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE), NULL)

    save.image()

}


# 0.1 Startup 
# -----------
# load data.table, suppressPackageStartupMessages unless there is an error
suppressPackageStartupMessages(library(data.table, quietly = TRUE))


# 1.0 Read input data into data.table
# -----------------------------------
rawdata <- fread(INPUT$rawdata)
processed <- as.data.table(read_xlsx(INPUT$processed))




# 2.0 Subset data
# ---------------
rawdata_s <- rawdata[!is.na(DRUG_ID) | grepl("^(L|R|A|N|B)", TAG)] 

## ------------------- GDSC SPECIFIC STEP ------------------- ##
# In the essence of using the processed GDSC data as much as possible
# we will normalize the rawdata using the gdscIC50 package

if (!require("gdscIC50", quietly = TRUE))
        pak::pkg_install("cancerrxgene/gdscIC50")

# drop all rows where TAG starts with L or FAIL
dt <- rawdata_s[!(grepl("^L", TAG) ) | !TAG == 'FAIL']

neg_control_TAGS <- dt[grepl("^NC", TAG), unique(TAG)]
if(length(neg_control_TAGS) > 1){
    print(paste0("More than one negative control found: ", paste0(neg_control_TAGS, collapse = ", ")))
    print(paste0("Setting all negative controls to: ", neg_control_TAGS[1]))
    dt[TAG %in% neg_control_TAGS, TAG := neg_control_TAGS[1]]
    neg_control_tag <- neg_control_TAGS[1]
}else if(length(neg_control_TAGS) == 0){
    stop("No negative control found")
}else{
    cat(paste0("Negative control found: ", neg_control_TAGS))
    neg_control_tag <- neg_control_TAGS[1]
}


normData <- suppressWarnings(
    gdscIC50::normalizeData(dt, trim=T, neg_control=neg_control_tag, pos_control="B"))

# normalizeData drops the DRUG_NAME column
print(paste0("Merging DRUG_NAME back into normData"))
normData <- merge(normData, treatment[,.(DRUG_ID, DRUG_NAME)], by.x = "DRUG_ID_lib",  by.y = "DRUG_ID", all.x = T)


