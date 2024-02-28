#!env bash
set -o pipefail

# install BiocManager
R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
R -e "BiocManager::install(c('affy', 'affyio', 'BiocParallel', 'hgu219.db', 'hgu219cdf', 'WGCNA', 'MultiAssayExperiment'))"

R -e 'BiocManager::install("preprocessCore", configure.args = "--disable-threading", force=TRUE)'