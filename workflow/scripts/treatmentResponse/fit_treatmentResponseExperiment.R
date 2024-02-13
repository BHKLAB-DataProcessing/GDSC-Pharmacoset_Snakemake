## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image("build_treatmentResponseExperiment.RData")
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)
}


suppressPackageStartupMessages(library(PharmacoGx))
suppressPackageStartupMessages(library(CoreGx))


## ------------------- Load Data ------------------- ##
tre <- readRDS(INPUT$tre)

tre_fit <- tre |> CoreGx::endoaggregate(
    {  # the entire code block is evaluated for each group in our group by
        # 1. fit a log logistic curve over the dose range
        fit <- PharmacoGx::logLogisticRegression(Dose, Viability,
            viability_as_pct=FALSE)
        # 2. compute curve summary metrics
        ic50 <- PharmacoGx::computeIC50(Dose, Hill_fit=fit)
        aac <- PharmacoGx::computeAUC(Dose, Hill_fit=fit)
        # 3. assemble the results into a list, each item will become a
        #   column in the target assay.
        list(
            HS=fit[["HS"]],
            E_inf = fit[["E_inf"]],
            EC50 = fit[["EC50"]],
            Rsq=as.numeric(unlist(attributes(fit))),
            aac_recomputed=aac,
            ic50_recomputed=ic50
        )
    },
    assay="sensitivity",
    target="profiles_recomputed",
    enlist=FALSE,  # this option enables the use of a code block for aggregation
    by=c("treatmentid", "sampleid"),
    nthread=THREADS  # parallelize over multiple cores to speed up the computation
)

saveRDS(tre_fit, OUTPUT$tre_fit)