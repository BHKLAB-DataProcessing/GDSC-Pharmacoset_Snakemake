
FROM  bioconductor/bioconductor_docker:devel 

# Install any needed packages specified in requirements.txt
RUN R -e "BiocManager::install(c('affy', 'affyio', 'BiocParallel', 'hgu219.db', 'hgu219cdf', 'WGCNA', 'MultiAssayExperiment'))"
RUN R -e "install.packages(c('data.table', 'qs'))"
