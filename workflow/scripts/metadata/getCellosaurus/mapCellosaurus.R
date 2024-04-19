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

message("Mapping GDSC.COSMIC_ID to Cellosaurus accession")
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

show(sampleMetadata)
message("Number of samples: ", nrow(sampleMetadata))

message("Annotating Cellosaurus accessions")
annotated_accessions <- AnnotationGx::annotateCellAccession(
    accessions = sampleMetadata$cellosaurus.acc,
)
message("Number of annotated accessions: ", nrow(annotated_accessions))

message("Number of unique categories: ")
annotated_accessions[, .N, by = "category"]

message("Number of unique sexOfCell: ")
annotated_accessions[, .N, by = "sexOfCell"]

annotated_accessions[, synonyms := sapply(synonyms, function(x) paste(x, collapse = "; "))]
annotated_accessions[, diseases := sapply(diseases, function(x) paste(x, collapse = "; "))]
annotated_accessions[, c("crossReferences", "hierarchy", "comments") := NULL]

names(annotated_accessions) <- paste0("cellosaurus.", names(annotated_accessions))
annotated_accessions <- unique(annotated_accessions)
show(annotated_accessions)

final_annotated <- merge(
    annotated_accessions,
    sampleMetadata, 
    by.x = "cellosaurus.accession",
    by.y = "cellosaurus.acc",
    all.x = TRUE
) |> unique()

final_annotated[, sampleid := "GDSC.sampleid"]

message("Writing output to: ", OUTPUT[['sampleMetadata']])

data.table::fwrite(
    final_annotated, 
    OUTPUT[['sampleMetadata']], 
    sep="\t", 
    quote = FALSE, 
    row.names = FALSE, 
    col.names = TRUE
)
