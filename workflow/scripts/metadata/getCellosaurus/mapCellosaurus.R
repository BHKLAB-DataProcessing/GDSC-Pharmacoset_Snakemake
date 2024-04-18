## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image("mapCellosaurus.RData")
    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)

}
load("mapCellosaurus.RData")

# Set this for mapping 
options("mc.cores" = THREADS)

options("log_level" = "INFO")   # AnnotationGx logging level

sampleMetadata <- data.table::fread(INPUT[['sampleMetadata']], sep="\t", header=T)

# The main columns in sampleMetadata are:
# GDSC.sampleid GDSC.Sample_Name CCLE.sampleid GDSC.BROAD_ID GDSC.RRID CMP.model_id CMP.sampleid GDSC.COSMIC_ID

# First we will map the GDSC.COSMIC_ID to the Cellosaurus accession
sampleMetadata[, cellosaurus.acc := {
    mapped <- AnnotationGx::mapCell2Accession(
        as.character(GDSC.COSMIC_ID), 
        from = "dr",
        BPPARAM = BPPARAM
        )
    return(mapped$ac)
    }
]
field_map <- list(
    "id" = "id",
    "ac" = "accession",
    "sy" = "synonyms",
    "misspelling" = "misspellings",
    "di" = "diseases",
    "ca" = "category",
    "sx" = "sexOfCell",
    "ag" = "ageAtSampling",
    "derived-from-site" = "samplingSite"
)

fields <- names(field_map)

sampleMetadata[, paste0("cellosaurus.", fields) := {
    mapped <- AnnotationGx::mapCell2Accession(
        as.character(GDSC.COSMIC_ID), 
        from = "dr",
        to = fields
    )
    return(mapped[, 1:length(fields)])
    }
]
save.image("mapCellosaurus.RData")

# rename using map
for (i in names(field_map)) {
    data.table::setnames(sampleMetadata, paste0("cellosaurus.", i), paste0("cellosaurus.", field_map[[i]]))
}




data.table::fwrite(sampleMetadata, OUTPUT[['sample_Cellosaurus_file']], sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
