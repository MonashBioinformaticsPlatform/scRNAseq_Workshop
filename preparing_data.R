library(tidyverse)
library(Seurat)
library(scDblFinder)

set.seed(1)

# Download data from Gene Expression Omnibus (GEO) ====

dir.create("data", showWarnings=FALSE)

# Download from GEO can be slow.
options(timeout=3600)

if (!file.exists("data/GSE96583_RAW.tar") || 
    file.info("data/GSE96583_RAW.tar")$size != 76195840) {
  download.file(
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE96583&format=file",
    "data/GSE96583_RAW.tar")
}

if (!file.exists("data/GSE96583_batch2.total.tsne.df.tsv.gz") ||
    file.info("data/GSE96583_batch2.total.tsne.df.tsv.gz")$size != 756342) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.total.tsne.df.tsv.gz",
    "data/GSE96583_batch2.total.tsne.df.tsv.gz")
}

if (!file.exists("data/GSE96583_batch2.genes.tsv.gz") ||
    file.info("data/GSE96583_batch2.genes.tsv.gz")$size != 277054) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.genes.tsv.gz",
    "data/GSE96583_batch2.genes.tsv.gz")
}


untar("data/GSE96583_RAW.tar", exdir="data")


# Load the data ====

cells <- read.table("data/GSE96583_batch2.total.tsne.df.tsv.gz", sep="\t")

counts <- cbind(
  ReadMtx("data/GSM2560248_2.1.mtx.gz", "data/GSM2560248_barcodes.tsv.gz", "data/GSE96583_batch2.genes.tsv.gz"),
  ReadMtx("data/GSM2560249_2.2.mtx.gz", "data/GSM2560249_barcodes.tsv.gz", "data/GSE96583_batch2.genes.tsv.gz"))

colnames(counts) <- rownames(cells)

# This dataset doesn't have Mito values. We have added Mitochondria values to illustrate what is the behavior o cells under stress  -
  # The logic is that if from all features expressed in a cell the proportion of mito is more than 5% - 10% is because the cell is stress. However this hard threshold miss the biological reasons 
  # why mito can be expressed. If the cell has other features expressed and mito could be a biological reaction. if the cell has low expression of other features and also high mito, the cell may be stressed or broken. 
  #  To illustrate this cases in the same dataset, adding mito values, is not trivial.

 

mt_rows <- grep("^MT-", rownames(counts))
n_cells <- dim(counts)[2]

# Generate simulated values (e.g., Poisson distribution) for mitochondrial genes
set.seed(123)
simulated_counts <- matrix(rpois(length(mt_rows) * n_cells, lambda = 3), nrow = length(mt_rows))

# Set 30% of simulated mitochondrial counts to zero to add sparsity
simulated_counts[sample(length(simulated_counts), length(simulated_counts) * 0.3)] <- 0

# Update the existing matrix with these simulated mitochondrial counts
counts[mt_rows, ] <- simulated_counts

# get kang2018.rds using the next lines
#download.file(
#"https://bioinformatics.erc.monash.edu/home/lper0012/SingleCellWorkshopData/data.tar","data.tar")

#untar("data.tar")

kang<-readRDS("~/data/tasks/training/SingleCellWorkshopData/data/kang2018.rds")

#colnames(kang)

#use the same cell names:
  
  
 # counts[,colnames(kang)]
  
  
cnts<- CreateSeuratObject(counts[,colnames(kang)], meta.data=cells)

cnts <- PercentageFeatureSet(cnts , pattern = "^MT-", col.name = "percent.mt", assay = "RNA")


VlnPlot(cnts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3 )


library(scrattch.io)
#counts 
write_dgCMatrix_h5(counts[,colnames(kang)], cols_are = "cell_names", "filtered_feature_bc_matrix.h5",
                               ref_name = "human" )

write.table(kang@meta.data,'metadata.txt')
