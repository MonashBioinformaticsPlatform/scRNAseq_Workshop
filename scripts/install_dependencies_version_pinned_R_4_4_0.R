# For R 4.4.0
#
install.packages("Seurat", version = "5.1.0")
install.packages("clustree", version = "2.1.6")
install.packages("harmony")
install.packages(c("dplyr", "tidyverse", "remotes", "R.utils", "RColorBrewer", "hdf5r"))
install.packages("bookdown")
install.packages("pander")

## Install required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.19")
BiocManager::install(c('SingleR', 'celldex',
                       'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment','edgeR'),
                     update=FALSE, ask=FALSE)

devtools::install_github("hhoeflin/hdf5r@dc4774c")
devtools::install_github('immunogenomics/presto@1.0.0')
