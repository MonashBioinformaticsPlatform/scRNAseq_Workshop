

# Load data ================


## Setup the Seurat Object --------


## The data set --------
#
# The dataset used in this workshop is a *modified* version derived from this
# study ([see here](https://pubmed.ncbi.nlm.nih.gov/29227470/)). It has been
# adapted to introduce additional complexity for instructional purposes. *Please
# refrain from drawing any biological conclusions from this data as it does not
# represent real experimental results*.
#
# This dataset represents human peripheral blood mononuclear cells (PBMCs),
# pooled from eight individual donors. Genetic differences among donors enable
# the identification of some cell doublets, enhancing data complexity. It
# includes two single-cell sequencing batches, one of which was stimulated with
# IFN-beta. Additionally, mitochondrial expression levels have been introduced to
# demonstrate how to interpret and apply mitochondrial thresholds for filtering
# purposes.
#
#### Note: What does the data look like?
#
# What do the input files look like? It varies, but this
# is the output of the CellRanger pipleine, described
# [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/gex-outputs)

├── analysis
│   ├── clustering
│   ├── diffexp
│   ├── pca
│   ├── tsne
│   └── umap
├── cloupe.cloupe
├── filtered_feature_bc_matrix
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── filtered_feature_bc_matrix.h5.  --> matrix to read in h5 format
├── metrics_summary.csv
├── molecule_info.h5
├── possorted_genome_bam.bam
├── possorted_genome_bam.bam.bai
├── raw_feature_bc_matrix
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── raw_feature_bc_matrix.h5
└── web_summary.html

###
#
# We start by reading in the data. There are several options for
# loading the data. The Read10X() function reads in the output of the
# [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
# pipeline from 10X, returning a unique molecular identified (UMI) count matrix.
# The values in this matrix represent the number of molecules for each feature
# (i.e. gene; row) that are detected in each cell (column).
#
# We next use the count matrix to create a Seurat object. The object
# serves as a container that contains both data (like the count matrix) and
# analysis (like PCA, or clustering results) for a single-cell dataset. For a
# technical discussion of the Seurat object structure, check out the [GitHub
# Wiki](https://github.com/satijalab/seurat/wiki). For example, the count matrix
# is stored in pbmc@assays$RNA@counts.

library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)


## Different ways of loading the data --------
#
# Example 1. Load your data using the path to the folder:
# filtered_feature_bc_matrix that is in the output folder of your cellranger run
# using the Read10X function.

## Load the PBMC dataset
# pbmc.data <- Read10X(data.dir = "outs/filtered_feature_bc_matrix")
## Initialize the Seurat object with the raw (non-normalized data).
# seurat_object <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200)


# Example 2. Load your data directing the ReadMtx function to each of the
# relevant files in the filtered_feature_bc_matrix folder in the outputs from your
# cellranger run. MTX is a simple text format for sparse matrices.


# expression_matrix <- ReadMtx(
#   mtx = "outs/filtered_feature_bc_matrix/count_matrix.mtx.gz", features = "outs/filtered_feature_bc_matrix/features.tsv.gz",
#   cells = "outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
# )
#
# seurat_object <- CreateSeuratObject(counts = expression_matrix)


# Example 3. Load your data using the Read10X_h5 function to each of the relevant
# HDF5 files. HDF5 is an efficient binary format.

pbmc.data <- Read10X_h5("data/filtered_feature_bc_matrix.h5")
metadata <- read.table("data/metadata.txt")



seurat_object <- CreateSeuratObject(counts = pbmc.data ,
                                 assay = "RNA", project = 'pbmc')



# Use AddMetaData to add new meta data to object
seurat_object  <- AddMetaData(object = seurat_object, metadata = metadata)



#   **What does data in a count matrix look like?**

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D","TCL1A","MS4A1"), 1:30]

# The . values in the matrix represent 0s (no molecules detected). Since most
# values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation
# whenever possible. This results in significant memory and speed savings for
# Drop-seq/inDrop/10x data.

dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size / sparse.size


# (PART) Single Cell Analysis ================


# QC Filtering ================
#
# The steps below encompass the standard pre-processing workflow for scRNA-seq
# data in Seurat. These represent the selection and filtration of cells based on
# QC metrics, data normalization and scaling, and the detection of highly variable
# features.


## QC and selecting cells for further analysis --------
#
#### Why do we need to do this?
#
# Low quality cells can add noise to your results leading you to the wrong
# biological conclusions. Using only good quality cells helps you to avoid this.
# Reduce noise in the data by filtering out low quality cells such as dying or
# stressed cells (high mitochondrial expression) and cells with few features that
# can reflect empty droplets.
#
####
#
# Seurat allows you to easily explore QC metrics and filter cells
# based on any user-defined criteria. A few QC metrics [commonly
# used](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/) by the community
# include
#
# * The number of unique genes detected in each cell.
#     + Low-quality cells or empty droplets will often have very few genes
#     + Cell doublets or multiplets may exhibit an aberrantly high gene count
# * Similarly, the total number of molecules detected within a cell (correlates
#   strongly with unique genes)
# * The percentage of reads that map to the mitochondrial genome
#     + Low-quality / dying cells often exhibit extensive mitochondrial
#       contamination
#     + We calculate mitochondrial QC metrics with the PercentageFeatureSet()
#       function, which calculates the percentage of counts originating from a set
#       of features
#     + We use the set of all genes starting with MT- as a set of mitochondrial
#       genes

# The $ operator can add columns to object metadata.
# This is a great place to stash QC stats
seurat_object$percent.mt <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

#### Challenge: The meta.data slot in the Seurat object
#
# Where are QC metrics stored in Seurat?
#
# * The number of unique genes and total molecules are automatically calculated
#   during CreateSeuratObject()
#     + You can find them stored in the object meta data
#
# What do you notice has changed within the meta.data table now that we have
# calculated mitochondrial gene proportion?
#
####
#
# In the example below, we visualize QC metrics. We will use these to filter
# cells.

# Visualize QC metrics as a violin plot
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships,
# but can be used for anything calculated by the object,
# i.e. columns in object metadata, PC scores etc.
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#### Challenge: Red Blood Cells
#
# Genes "HBA1", "HBA2", and "HBB" are components of hemoglobin in red blood cells.
#
# * Use PercentageFeatureSet, passing these genes to the "features" argument, to
#   find cells that might be red blood cells.
# * How do cells high in these genes differ from other cells, in terms of number
#   of features or total count?
# * Should we remove these cells?
#
####
#
# Let's look again at the number of features (genes) to the percent mitochondrial
# genes plot.

VlnPlot(seurat_object, features = "nFeature_RNA")

# Zoom in to the max and min
VlnPlot(seurat_object, features = "nFeature_RNA") + scale_y_continuous(limits = c(1000,3000))
VlnPlot(seurat_object, features = "nFeature_RNA", y.max =1000)

FeatureScatter(seurat_object, feature1 = "nFeature_RNA", feature2 = "percent.mt")

# You can check different thresholds of mito percentage.

#Number of cells that would be left after filters
# Proportion of cells with less than 5% mito

mean(seurat_object$percent.mt < 5)

# Proportion of cells with less than 2% mito

mean(seurat_object$percent.mt < 2)

# Ok, let's go with these filters:
#
# * We filter cells to have >200 unique features
# * We filter cells that have >5% mitochondrial counts
#
# Let's apply this and subset our data. This will remove the cells we think are of
# poor quality.

seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt < 5)

# Let's replot the feature scatters and see what they look like.

FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# We also wondered if cells with high counts might be doublets. Should we also
# filter cells with very high counts? With this data, we know for certain some of
# the doublets!

FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="multiplets")


# Normalisation ================
#
#### Why do we need to do this?
#
# The sequencing depth can be different per cell. This can bias the counts of
# expression showing higher numbers for more sequenced cells leading to the wrong
# biological conclusions. To correct this the feature counts are normalized.
#
####
#
# After removing unwanted cells from the dataset, the next step is to normalize
# the data. By default, we employ a global-scaling normalization method
# "LogNormalize" that normalizes the feature expression measurements for each
# cell by the total expression, multiplies this by a scale factor (10,000 by
# default), and log-transforms the result. Normalized values are stored in
# seurat_object$RNA@data.

seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1e4)

# For clarity, in this previous line of code (and in future commands), we provide
# the default values for certain parameters in the function call. However, this
# isn't required and the same behavior can be achieved with:

seurat_object <- NormalizeData(seurat_object)

# There are other options for normalization such as
# [SCTtransform]( https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)
# which was popularized in 2019, however Log base normalization
# continued to be preferred as they perform better. [See
# here](https://www.nature.com/articles/s41592-023-01814-1) the for more details.


# PCAs and UMAPs ================


## Identification of highly variable features (feature selection) --------
#
#### Why do we need to do this?
#
# Identifying the most variable features allows retaining the real biological
# variability of the data and reduce noise in the data.
#
####
#
# We next calculate a subset of features that exhibit high cell-to-cell
# variation in the dataset (i.e, they are highly expressed in some
# cells, and lowly expressed in others). The Seurat authors and
# [others](https://www.nature.com/articles/nmeth.2645) have found that focusing
# on these genes in downstream analysis helps to highlight biological signal in
# single-cell datasets.
#
# The procedure in Seurat is described in detail
# [here](https://doi.org/10.1016/j.cell.2019.05.031), and improves on previous
# versions by directly modeling the mean-variance relationship inherent in
# single-cell data, and is implemented in the FindVariableFeatures() function.
# By default, Seurat returns 2,000 features per dataset. These will be used in
# downstream analysis, like PCA.

seurat_object <- FindVariableFeatures(seurat_object, selection.method = 'vst', nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


## Scaling the data --------
#
#### Why do we need to do this?
#
# Highly expresed genes can overpower the signal of other less expresed genes with
# equal importance. Within the same cell the assumption is that the underlying RNA
# content is constant. Aditionally, If variables are provided in vars.to.regress,
# they are individually regressed against each feature, and the resulting
# residuals are then scaled and centered. This step allows controling for cell
# cycle and other factors that may bias your clustering.
#
####
#
# Next, we apply a linear transformation ('scaling') that is a standard
# pre-processing step prior to dimensional reduction techniques like PCA. The
# ScaleData() function:
#
# * Shifts the expression of each gene, so that the mean expression across cells
#   is 0
# * Scales the expression of each gene, so that the variance across cells is 1
#     + This step gives equal weight in downstream analyses, so that
#       highly-expressed genes do not dominate
# * The results of this are stored in seurat_object$RNA@scale.data

all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)


# Dimensionality reduction ================
#
#### Why do we need to do this?
#
# Imagine each gene represents a dimension - or an axis on a plot. We could
# plot the expression of two genes with a simple scatterplot. But a genome has
# thousands of genes - how do you collate all the information from each of those
# genes in a way that allows you to visualise it in a 2 dimensional image. This
# is where dimensionality reduction comes in, we calculate meta-features that
# contains combinations of the variation of different genes. From thousands of
# genes, we end up with 10s of meta-features
#
####


## Perform linear dimensional reduction --------
#
# Next we perform PCA on the scaled data. By default, only the previously
# determined variable features are used as input, but can be defined using
# features argument if you wish to choose a different subset.

seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Seurat provides several useful ways of visualizing both cells and features that
# define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()

# Examine and visualize PCA results a few different ways
print(seurat_object$pca, dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat_object, dims = 1:2, reduction = 'pca')
DimPlot(seurat_object, reduction = 'pca')

# In particular DimHeatmap() allows for easy exploration of the primary sources
# of heterogeneity in a dataset, and can be useful when trying to decide which PCs
# to include for further downstream analyses. Both cells and features are ordered
# according to their PCA scores. Setting cells to a number plots the 'extreme'
# cells on both ends of the spectrum, which dramatically speeds plotting for large
# datasets. Though clearly a supervised analysis, we find this to be a valuable
# tool for exploring correlated feature sets.

DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE)


## Determine the 'dimensionality' of the dataset --------
#
# To overcome the extensive technical noise in any single feature for scRNA-seq
# data, Seurat clusters cells based on their PCA scores, with each PC essentially
# representing a 'metafeature' that combines information across a correlated
# feature set. The top principal components therefore represent a robust
# compression of the dataset. However, how many components should we choose to
# include? 10? 20? 100?
#
# -----
#
# *Note*: The Seurat developers suggest using a JackStraw
#  resampling test to determine this. See [Macosko *et
#  al*](http://www.cell.com/abstract/S0092-8674(15)00549-8),
#  and the original [seurat_object3
#  vignette](https://satijalab.org/seurat/articles/seurat_object3k_tutorial.html#determine-the-dimensionality-of-the-dataset-1).
#  We're going to use an Elbow Plot instead here, because its much quicker.
#
# -----
#
# An alternative heuristic method generates an 'Elbow plot': a ranking of
# principle components based on the percentage of variance explained by each
# one (ElbowPlot() function). In this example, we can observe an 'elbow' around
# PC9-10, suggesting that the majority of true signal is captured in the first 10
# PCs.

ElbowPlot(seurat_object)

# Identifying the true dimensionality of a dataset -- can be challenging/uncertain
# for the user. We therefore suggest these three approaches to consider. The
# first is more supervised, exploring PCs to determine relevant sources of
# heterogeneity, and could be used in conjunction with GSEA for example. The
# second implements a statistical test based on a random null model, but is
# time-consuming for large datasets, and may not return a clear PC cutoff. The
# third is a heuristic that is commonly used, and can be calculated instantly. In
# this example, all three approaches yielded similar results, but we might have
# been justified in choosing anything between PC 7-12 as a cutoff.
#
# We chose 10 here, but encourage users to consider the following:
#
# * In the original version of this vignette with the PBMC3k dataset, genes
#   strongly associated with PCs 12 and 13 defined rare immune subsets (i.e. MZB1
#   is a marker for plasmacytoid DCs). However, these groups are so rare, they
#   are difficult to distinguish from background noise for a dataset of this size
#   without prior knowledge.
# * We encourage users to repeat downstream analyses with a different number
#   of PCs (10, 15, or even 50!). As you will observe, the results often do not
#   differ dramatically.
# * We advise users to err on the higher side when choosing this parameter. For
#   example, performing downstream analyses with only 5 PCs does significantly and
#   adversely affect results.


## Run non-linear dimensional reduction (UMAP/tSNE) --------
#
# Seurat offers several non-linear dimensional reduction techniques, such as tSNE
# and UMAP, to visualize and explore these datasets. The goal of these algorithms
# is to learn the underlying manifold of the data in order to place similar
# cells together in low-dimensional space. Cells within the graph-based clusters
# determined above should co-localize on these dimension reduction plots. As input
# to the UMAP and tSNE, we suggest using the same PCs as input to the clustering
# analysis.

seurat_object <- RunUMAP(seurat_object, dims = 1:10)

DimPlot(seurat_object, reduction = 'umap')

#### Challenge: PC genes
#
# You can plot gene expression on the UMAP with the FeaturePlot() function.
#
# Try out some genes that were highly weighted in the principal component
# analysis. How do they look?
#
####


## Save --------
#
# You can save the object at this point so that it can easily be loaded back
# in with readRDS() without having to rerun the computationally intensive steps
# performed above, or easily shared with collaborators.

saveRDS(seurat_object, file = "seurat_object_tutorial_saved.rds")

# Tip: For faster saving and loading, try the "qs" package.


# Data set integration with Harmony ================
#
### Why do we need to do this?
#
# You can have data coming from different samples, batches or experiments and you
# will need to combine them.
#
###
#
# * ind identifies a cell as coming from one of 8 individuals.
# * stim identifies a cell as control or stimulated with IFN-beta.
# * cell contains the cell types identified by the creators of this data set.
# * multiplets classifies cells as singlet or doublet.

DimPlot(seurat_object, reduction="umap", group.by="ind")
DimPlot(seurat_object, reduction="umap", group.by="stim")

seurat_object<- FindNeighbors(seurat_object, reduction="pca", dims=1:10)
seurat_object <- FindClusters(seurat_object, resolution=0.25)
seurat_object$pca_clusters <- seurat_object$seurat_clusters

DimPlot(seurat_object, reduction="umap", group.by="pca_clusters")

# There is a big difference between unstimulated and stimulated cells. This
# has split cells of the same type into pairs of clusters. If the difference
# was simply uniform, we could regress it out (e.g. using ScaleData(...,
# vars.to.regress="stim")). However, as can be seen in the PCA plot, the
# difference is not uniform and we need to do something cleverer.
#
# We will use [Harmony](https://github.com/immunogenomics/harmony), which can
# remove non-uniform effects. We will try to remove both the small differences
# between individuals and the large difference between the unstimulated and
# stimulated cells.
#
# Harmony operates only on the PCA scores. The original gene expression levels
# remain unaltered.

library(harmony)

seurat_object <- RunHarmony(seurat_object, c("stim", "ind"), reduction="pca",reduction.save="harmony")

# This has added a new set of reduced dimensions to the Seurat object,
# seurat_object$harmony which is a modified version of the existing
# seurat_object$pca reduced dimensions. The PCA plot shows a large difference
# between 'ctrl' and 'stim', but this has been removed in the harmony reduction.

DimPlot(seurat_object, reduction="pca", group.by="stim")
DimPlot(seurat_object, reduction="harmony", group.by="stim")

# We can use harmony the same way we used the pca reduction to compute a UMAP
# layout or to find clusters.

seurat_object <- RunUMAP(seurat_object, reduction="harmony", dims=1:10, reduction.name="umap_harmony")

DimPlot(seurat_object, reduction="umap_harmony", group.by="stim")

seurat_object <- FindNeighbors(seurat_object, reduction="harmony", dims=1:10)
seurat_object <- FindClusters(seurat_object, resolution=0.25)
seurat_object$harmony_clusters <- seurat_object$seurat_clusters

DimPlot(seurat_object, reduction="umap_harmony", group.by="harmony_clusters")
DimPlot(seurat_object, reduction="umap", group.by="harmony_clusters")

# Having found a good set of clusters, we would usually perform differential
# expression analysis on the original data and include batches/runs/individuals
# as predictors in the linear model. In this example we could now compare
# un-stimulated and stimulated cells within each cluster. A particularly nice
# statistical approach that is possible here would be to convert the counts
# to pseudo-bulk data for the eight individuals, and then apply a bulk RNA-Seq
# differential expression analysis method. However there is still the problem that
# unstimulated and stimulated cells were processed in separate batches.


# Clustering ================
#
#### Why do we need to do this?
#
# Clustering the cells will allow you to visualise the variability of your data,
# can help to segregate cells into cell types.
#
####


## Cluster cells --------
#
# Seurat v3 applies a graph-based clustering approach,
# building upon initial strategies in ([Macosko *et
# al*](http://www.cell.com/abstract/S0092-8674(15)00549-8)). Importantly, the
# *distance metric* which drives the clustering analysis (based on previously
# identified PCs) remains the same. However, our approach to partitioning
# the cellular distance matrix into clusters has dramatically improved. Our
# approach was heavily inspired by recent manuscripts which applied graph-based
# clustering approaches to scRNA-seq data [[SNN-Cliq, Xu and Su, Bioinformatics,
# 2015]](http://bioinformatics.oxfordjournals.org/content/early/2015/02/10/bioinformatics.btv088.abstract)
# and CyTOF data [[PhenoGraph, Levine *et al*., Cell,
# 2015]](http://www.ncbi.nlm.nih.gov/pubmed/26095251). Briefly, these methods
# embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph,
# with edges drawn between cells with similar feature expression patterns, and
# then attempt to partition this graph into highly interconnected 'quasi-cliques'
# or 'communities'.
#
# As in PhenoGraph, we first construct a KNN graph based on the euclidean distance
# in PCA space, and refine the edge weights between any two cells based on the
# shared overlap in their local neighborhoods (Jaccard similarity). This step is
# performed using the FindNeighbors() function, and takes as input the previously
# defined dimensionality of the dataset (first 10 PCs).
#
# To cluster the cells, we next apply modularity optimization techniques such
# as the Louvain algorithm (default) or SLM [[SLM, Blondel *et al*., Journal of
# Statistical Mechanics]](http://dx.doi.org/10.1088/1742-5468/2008/10/P10008),
# to iteratively group cells together, with the goal of optimizing the standard
# modularity function. The FindClusters() function implements this procedure, and
# contains a resolution parameter that sets the 'granularity' of the downstream
# clustering, with increased values leading to a greater number of clusters. We
# find that setting this parameter between 0.4-1.2 typically returns good results
# for single-cell datasets of around 3K cells. Optimal resolution often increases
# for larger datasets. The clusters can be found using the Idents() function.

seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(seurat_object), 5)

# Check out the clusters.

DimPlot(seurat_object,reduction = "umap_harmony")
# Equivalent to
# DimPlot(seurat_object,reduction="umap", group.by="seurat_clusters")
# DimPlot(seurat_object,reduction="umap", group.by="RNA_snn_res.0.5")

#### Challenge: Try different cluster settings
#
# Run FindNeighbours and FindClusters again, with a different number of dimensions
# or with a different resolution. Examine the resulting clusters using DimPlot.
#
# To maintain the flow of this tutorial, please put the output of this exploration
# in a different variable, such as seurat_object2!
#
####


## Choosing a cluster resolution --------
#
# Its a good idea to try different resolutions when clustering to identify the
# variability of your data.

resolution = 2
seurat_object <- FindClusters(object = seurat_object, reduction = "umap_harmony", resolution = seq(0.1, resolution, 0.1),
    dims = 1:10)

# the different clustering created
names(seurat_object@meta.data)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_object), 5)

# Plot a clustree to decide how many clusters you have and what resolution capture
# them.

library(clustree)
clustree(seurat_object, prefix = "RNA_snn_res.",show_axis=TRUE) + theme(legend.key.size = unit(0.05, "cm"))

# Name cells with the corresponding cluster name at the resolution you pick. This
# case we are happy with 0.5.

# The name of the cluster is prefixed with 'RNA_snn_res' and the number of the resolution
Idents(seurat_object) <- seurat_object$RNA_snn_res.0.5

# Plot the UMAP with colored clusters with Dimplot

DimPlot(seurat_object, reduction = "umap_harmony", label = TRUE, repel = TRUE, label.box = TRUE) + NoLegend()


# Cluster Markers ================
#
#### Why do we need to do this?
#
# Single cell data helps to segragate cell types. Use markers to identify cell
# types. warning: In this example the cell types/markers are well known and making
# this step easy, but in reality this step needs the experts curation.
#
####


## Finding differentially expressed features (cluster biomarkers) --------
#
# Seurat can help you find markers that define clusters via differential
# expression. By default, it identifies positive and negative markers of a single
# cluster (specified in ident.1), compared to all other cells. FindAllMarkers()
# automates this process for all clusters, but you can also test groups of
# clusters vs. each other, or against all cells.
#
# The min.pct argument requires a feature to be detected at a minimum percentage
# in either of the two groups of cells, and the thresh.test argument requires
# a feature to be differentially expressed (on average) by some amount between
# the two groups. You can set both of these to 0, but with a dramatic increase
# in time - since this will test a large number of features that are unlikely
# to be highly discriminatory. As another option to speed up these computations,
# max.cells.per.ident can be set. This will downsample each identity class to have
# no more cells than whatever this is set to. While there is generally going to
# be a loss in power, the speed increases can be significant and the most highly
# differentially expressed features will likely still rise to the top.

# find all markers of cluster 2
cluster2.markers <- FindMarkers(seurat_object, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(seurat_object, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_object.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

# Seurat has several tests for differential expression which can be set with the
# test.use parameter (see our [DE vignette](de_vignette.html) for details). For
# example, the ROC test returns the 'classification power' abs(AUC-0.5)*2 for any
# individual marker, ranging from 0 = random to 1 = perfect.

cluster0.markers <- FindMarkers(seurat_object, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# We include several tools for visualizing marker expression. VlnPlot() (shows
# expression probability distributions across clusters), and FeaturePlot()
# (visualizes feature expression on a tSNE or PCA plot) are our most commonly
# used visualizations. We also suggest exploring RidgePlot(), CellScatter(), and
# DotPlot() as additional methods to view your dataset.

VlnPlot(seurat_object, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(seurat_object, features = c("NKG7", "PF4"), slot = 'counts', log = TRUE)

FeaturePlot(seurat_object, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A"))
FeaturePlot(seurat_object, features = c("FCGR3A", "LYZ", "PPBP", "CD8A"))

#   **Other useful plots**
# These are ridgeplots, cell scatter plots and dotplots. Replace FeaturePlot with
# the other functions.

RidgePlot(seurat_object, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A"))
RidgePlot(seurat_object, features = c("FCGR3A", "LYZ", "PPBP", "CD8A"))

# For CellScatter plots, will need the cell id of the cells you want to look at.
# You can get this from the cell metadata (seurat_object@meta.data).

head( seurat_object@meta.data )

CellScatter(seurat_object, cell1 = "AGGGCGCTATTTCC-1", cell2 = "GGAGACGATTCGTT-1")

CellScatter(seurat_object, cell1 = "GGAGACGATTCGTT-1", cell2 = "TACGAGACCTATTC-1")




# DotPlots

DotPlot(seurat_object, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features. In
# this case, we are plotting the top 10 markers (or all markers if less than 10)
# for each cluster.

maxcells=1500
top10 <- seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(subset(seurat_object, downsample = maxcells), features = top10$gene) + NoLegend()



## Use makers to label or find a cluster --------
#
# If you know markers for your cell types, use AddModuleScore to label them.

genes_markers <- list(Naive_CD4_T = c("IL7R", "CCR7"))

seurat_object <- AddModuleScore(object = seurat_object, features = genes_markers, ctrl = 5, name = "Naive_CD4_T",
    search = TRUE)


# notice the name of the cluster has a 1 at the end
names(seurat_object@meta.data)

# label that cell type
seurat_object$cell_label = NA
seurat_object$cell_label[seurat_object$Naive_CD4_T1 > 1] = "Naive_CD4_T"
Idents(seurat_object) = seurat_object$cell_label

# plot
# Using a custom colour scale
FeaturePlot(seurat_object, features = "Naive_CD4_T1", label = TRUE, repel = TRUE, reduction = "umap_harmony") + scale_colour_gradientn(colours = c("lightblue","beige","red"))


## Assigning cell type identity to clusters --------
#
# Fortunately in the case of this dataset, we can use canonical markers to easily
# match the unbiased clustering to known cell types:
#
# | Markers | Cell Type
# |---------------|----------
# | IL7R, CCR7 | Naive CD4+ T
# | CD14, LYZ | CD14+ Mono
# | IL7R, S100A4 | Memory CD4+
# | MS4A1 | B
# | CD8A | CD8+ T
# | FCGR3A, MS4A7 | FCGR3A+ Mono
# | GNLY, NKG7 | NK
# | FCER1A, CST3 | DC
# | PPBP | Platelet



DimPlot(seurat_object,group.by = "RNA_snn_res.0.2",reduction = "umap_harmony")|FeaturePlot(seurat_object, features = c( "MS4A1"),reduction = "umap_harmony")

#### Challenge: Match cluster numbers with cell labels
#
# Use the markers provided and the resolution 0.2 to identity the cell labels
#
#   **code ideas?**

Idents(seurat_object) <- seurat_object$RNA_snn_res.0.2
# this is not the proper order. Make sure the labels are in the same order that the numers they should replace
new.cluster.ids <- c("Naive CD4 T","B cells","CD14+ Monocytes","CD4 T cells","CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes","Megakaryocytes" )
names(new.cluster.ids) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
DimPlot(seurat_object, reduction = 'umap_harmony', label = TRUE, pt.size = 0.5) + NoLegend()

####
#
#### save the plots

#### save the seurat object

saveRDS(seurat_object, file = "seurat_object3k_final.rds")


# (PART) Futher Analysis ================


# SingleR ================

#install.packages("BiocManager")
#BiocManager::install(c("SingleCellExperiment","SingleR","celldex"),ask=F)
library(SingleCellExperiment)
library(SingleR)
library(celldex)

# In this workshop we have focused on the Seurat package. However, there
# is another whole ecosystem of R packages for single cell analysis within
# Bioconductor. We won't go into any detail on these packages in this
# workshop, but there is good material describing the object type online :
# [OSCA](https://robertamezquita.github.io/orchestratingSingleCellAnalysis/data-infrastructure.html).
#
# For now, we'll just convert our Seurat object into an object called
# SingleCellExperiment. Some popular packages from Bioconductor that work with
# this type are Slingshot, Scran, Scater.

sce <- as.SingleCellExperiment(seurat_object)
sce

# We will now use a package called SingleR to label each cell. SingleR uses a
# reference data set of cell types with expression data to infer the best label
# for each cell. A convenient collection of cell type reference is in the celldex
# package which currently contains the follow sets:

ls('package:celldex')

# In this example, we'll use the HumanPrimaryCellAtlasData set, which contains
# high-level, and fine-grained label types. Lets download the reference dataset

# This too is a sce object,
# colData is equivalent to seurat's metadata
ref.set <- celldex::HumanPrimaryCellAtlasData()

# The "main" labels.

unique(ref.set$label.main)

# An example of the types of "fine" labels.

head(unique(ref.set$label.fine))

# Now we'll label our cells using the SingleCellExperiment object, with the above
# reference set.

pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.main)

# Keep any types that have more than 10 cells to the label, and put those labels
# back on our Seurat object and plot our on our umap.

lbls.keep <- table(pred.cnts$labels)>10
seurat_object$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')
DimPlot(seurat_object, reduction='umap_harmony', group.by='SingleR.labels')

# Compare cell labels by different annotation methods:

DimPlot(seurat_object,group.by = "RNA_snn_res.0.2",reduction = "umap_harmony")

DimPlot(seurat_object,group.by = "cell",reduction = "umap_harmony")





# Differential Expression ================
#
# There are many different methods for calculating differential expression between
# groups in scRNAseq data. There are a number of review papers worth consulting on
# this topic.
#
# There is the [Seurat differential expression
# Vignette](https://satijalab.org/seurat/archive/v3.1/de_vignette.html) which
# walks through the variety implemented in Seurat.
#
# There is also a good discussion of useing [pseudobulk
# approaches](http://bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html#creating-pseudo-bulk-samples)
# which is worth checking out if youre planning differential expression analyses.

head(seurat_object@meta.data)

# How cells from each condition do we have?

table(seurat_object$stim)

# How many cells per individuals per group?

table(seurat_object$ind, seurat_object$stim)

# And for each sample, how many of each cell type has been classified?

table(paste(seurat_object$ind,seurat_object$stim), seurat_object$cell)


## Prefiltering --------
#
#### Why do we need to do this?
#
# If expression is below a certain level, it will be almost impossible to see any
# differential expression.
#
####
#
# When doing differential expression, you generally ignore genes with low
# expression.
# In single cell datasets, there are many genes like this. Filtering here to make
# our dataset smaller so it runs quicker, and there is less aggressive correction
# for multiple hypotheses.
#
# How many genes before filtering?

seurat_object

# How many copies of each gene are there?

total_per_gene <- rowSums(GetAssayData(seurat_object, assay='RNA', slot='counts'))
hist(log10(total_per_gene))

# Lets keep only those genes with at least 50 copies across the entire experiment.

seurat_object<- seurat_object[total_per_gene >= 50, ]

# How many genes after filtering?

seurat_object

# We might like to see the effect of IFN-beta stimulation on each cell type
# individually. For the purposes of this workshop, just going to test one cell
# type; CD14+ Monocytes
#
# An easy way is to subset the object.

# Set idents to 'cell' column.
Idents(seurat_object) <- seurat_object$cell
DimPlot(seurat_object, reduction = "umap_harmony")
seurat_object_celltype <- seurat_object[, seurat_object$cell == "CD14+ Monocytes" ]
DimPlot(seurat_object_celltype, reduction = "umap_harmony")


##  Default Wilcox test --------
#
# To run this test, we change the Idents to the factor(column) we want to test. In
# this case, that's 'stim'.

# Change Ident to Condition
Idents(seurat_object_celltype) <- seurat_object_celltype$stim

# default, wilcox test
de_result_wilcox <- FindMarkers(seurat_object_celltype,
            ident.1 = 'stim',
            ident.2 = 'ctrl',
            logfc.threshold = 0, # Give me ALL results
            min.pct = 0
            )

# Add average expression for plotting
de_result_wilcox$AveExpr<- rowMeans(seurat_object_celltype[rownames(de_result_wilcox),])

# Look at the top differentially expressed genes.

head(de_result_wilcox)

p1 <- ggplot(de_result_wilcox, aes(x=AveExpr, y=avg_log2FC, col=p_val_adj < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) +
  theme_bw() +
  ggtitle("Wilcox Test")


p2 <- ggplot(de_result_wilcox, aes(x=avg_log2FC, y=-log10(p_val), col=p_val_adj < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) +
  theme_bw() +
  ggtitle("Wilcox Test (Volcano)")

p1 + p2


## Seurat Negative binomial --------
#
# Negative binonial test is run almost the same way - just need to specify it
# under 'test.use'


# Change Ident to Condition
Idents(seurat_object_celltype) <- seurat_object_celltype$stim

# default, wilcox test
de_result_negbinom <- FindMarkers(seurat_object_celltype,
            test.use="negbinom", # Choose a different test.
            ident.1 = 'stim',
            ident.2 = 'ctrl',
            logfc.threshold = 0, # Give me ALL results
            min.pct = 0
)

# Add average expression for plotting
de_result_negbinom$AveExpr<- rowMeans(seurat_object_celltype[rownames(de_result_negbinom),])

# Look at the top differentially expressed genes.

head(de_result_negbinom)

p1 <- ggplot(de_result_negbinom, aes(x=AveExpr, y=avg_log2FC, col=p_val_adj < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) +
  theme_bw() +
  ggtitle("Negative Bionomial Test")


p2 <- ggplot(de_result_negbinom, aes(x=avg_log2FC, y=-log10(p_val), col=p_val_adj < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) +
  theme_bw() +
  ggtitle("Negative Bionomial Test (Volcano)")

p1 + p2


## Pseudobulk --------
#
# Pseudobulk analysis is an option where you have biological replicates. It is
# essentially pooling the individual cell counts and treating your expreiment like
# a bulk RNAseq.
#
# First, you need to build a pseudobulk matrix - the AggregateExpression()
# function can do this, once you set the 'Idents' of your seurat object to your
# grouping factor (here, thats a combination of individual+treatment called
# 'sample', instead of the 'stim' treatment column).

# Tools for bulk differential expression
library(limma)
library(edgeR)


# Change idents to ind for grouping.
seurat_object_celltype$sample <- factor(paste(seurat_object_celltype$stim, seurat_object_celltype$ind, sep="_"))
Idents(seurat_object_celltype) <- seurat_object_celltype$sample

# THen pool together counts in those groups
# AggregateExperssion returns a list of matricies - one for each assay requested (even just requesting one)
pseudobulk_matrix_list <- AggregateExpression( seurat_object_celltype,  slot = 'counts', assays='RNA' )
pseudobulk_matrix      <- pseudobulk_matrix_list[['RNA']]
colnames(pseudobulk_matrix) <- as.character(colnames(pseudobulk_matrix)) # Changes colnames to simple text
pseudobulk_matrix[1:5,1:4]

# Now it looks like a bulk RNAseq experiment, so treat it like one.
#
# We can use the popular limma package for differential expression. Here is one
# [tutorial](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html),
# and the hefty reference manual is hosted by
# [bioconductor](https://bioconductor.org/packages/release/bioc/html/limma.html).
#
# In brief, this code below constructs a linear model for this experiment that
# accounts for the variation in individuals and treatment. It then tests for
# differential expression between 'stim' and 'ctrl' groups.

dge <- DGEList(pseudobulk_matrix)
dge <- calcNormFactors(dge)

# Remove _ or - and everything after it - yeilds stim group
stim <- gsub("[-_].*","",colnames(pseudobulk_matrix))

# Removing everything before the _ or - for the individual, then converting those numerical ind explictiy to text. Else limma will treat them as numbers!
ind  <- as.character(gsub(".*[-_]","",colnames(pseudobulk_matrix)))

design <- model.matrix( ~0 + stim + ind)
vm  <- voom(dge, design = design, plot = FALSE)
fit <- lmFit(vm, design = design)

contrasts <- makeContrasts(stimstim - stimctrl, levels=coef(fit))
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)

de_result_pseudobulk <- topTable(fit, n = Inf, adjust.method = "BH")
de_result_pseudobulk <- arrange(de_result_pseudobulk , adj.P.Val)

# Look at the significantly differentially expressed genes:

head(de_result_pseudobulk)

p1 <- ggplot(de_result_pseudobulk, aes(x=AveExpr, y=logFC, col=adj.P.Val < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) +
  theme_bw() +
  ggtitle("Pseudobulk")


p2 <- ggplot(de_result_pseudobulk, aes(x=logFC, y=-log10(P.Value), col=adj.P.Val < 0.05)) +
  geom_point() +
  scale_colour_manual(values=c('TRUE'="red",'FALSE'="black")) +
  theme_bw() +
  ggtitle("Pseudobulk Test (Volcano)")

p1 + p2


#### Discussion
#
# These methods give different results. How would you decide which to use? How
# could you check an individual gene?
#
####


# Cell cycle Assignment ================
#
# In some datasets, the phase of cell cycle that a cell is in (G1/G2M/S) can
# account for
# alot of the observed transcriptomic variation. There may be clustering by phase,
# or
# separation in the UMAP by phase.
#
# Seurat provides a simple method for assigning cell cycle state to each cell.
# Other methods are available.
#
# More information about assigning cell cycle states to cells is in the [cell
# cycle vignette](https://satijalab.org/seurat/articles/cell_cycle_vignette.html)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Use those lists with the cell cycle scoring function in Seurat.
seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes)

# Which adds S.Score, G2M.Score and Phase calls to the metadata.

head(seurat_object@meta.data)

# We can then check the cell phase on the UMAP. In this dataset, phase isn't
# driving the clustering, and would not require any further handling.

DimPlot(seurat_object, reduction = 'umap_harmony', group.by = "Phase")

# Where a bias _is_ present, your course of action depends on the task at
# hand. It might involve 'regressing out' the cell cycle variation when
# scaling data ScaleData(kang, vars.to.regress="Phase"), omitting cell-cycle
# dominated clusters, or just accounting for it in your differential expression
# calculations.
#
# If you are working with non-human data, you will need to convert these gene
# lists, or find new cell cycle associated genes in your species.
