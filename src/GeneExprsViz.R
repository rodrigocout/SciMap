#http://www.bioconductor.org/packages/release/bioc/vignettes/Nebulosa/inst/doc/introduction.html#

#Tutorial from Nebulosa package

library("Nebulosa")
library("scater")
library("scran")
library("DropletUtils")
library("BiocFileCache")

bfc <- BiocFileCache(ask = FALSE)
data_file <- bfcrpath(bfc, file.path(
  "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell",
  "pbmc3k",
  "pbmc3k_filtered_gene_bc_matrices.tar.gz"
))

untar(data_file, exdir = tempdir())
pbmc <- read10xCounts(file.path(tempdir(),
  "filtered_gene_bc_matrices",
  "hg19"
))

rownames(pbmc) <- uniquifyFeatureNames(rowData(pbmc)[["ID"]],
                                       rowData(pbmc)[["Symbol"]])


#Quality control
# remove features that are not expressed in at least 3 cells
i <- rowSums(counts(pbmc) > 0)
is_expressed <- i > 3
pbmc <- pbmc[is_expressed, ]

#And cells not expressing at least one UMI in at least 200 genes.

i <- colSums(counts(pbmc) > 0)
is_expressed <- i > 200
pbmc <- pbmc[,is_expressed]

# remove outlier cells based on the number of genes being expressed in each cell, library size, and expression of mitochondrial genes using the #perCellQCMetrics and quickPerCellQC functions from the scater package.

is_mito <- grepl("^MT-", rownames(pbmc))
qcstats <- perCellQCMetrics(pbmc, subsets = list(Mito = is_mito))
qcfilter <- quickPerCellQC(qcstats, percent_subsets = c("subsets_Mito_percent"))

#Data normalization

logcounts(pbmc) <- log1p(counts(pbmc) / colSums(counts(pbmc)) * 1e4)


Dimensionality reduction

dec <- modelGeneVar(pbmc)
top_hvgs <- getTopHVGs(dec, n = 3000) #top 3000 most highly-variable genes

set.seed(66)
pbmc <- runPCA(pbmc, scale = TRUE, subset_row = top_hvgs)

#UMAP

pbmc <- runUMAP(pbmc, dimred = "PCA")

#Visualize data with Nebulosa

#Plot with one gene

plot_density(pbmc, "CD4")

plotUMAP(pbmc, colour_by = "CD4")

#Plot multiple genes
p3 <- plot_density(pbmc, c("CD8A", "CCR7"))
p3 + plot_layout(ncol = 1)


p4 <- plot_density(pbmc, c("CD8A", "CCR7"), joint = TRUE)
p4 + plot_layout(ncol = 1)

#--
plotUMAP(pbmc, colour_by = "label", text_by = "label")

#--
p_list <- plot_density(pbmc, c("CD8A", "CCR7"), joint = TRUE, combine = FALSE)
p_list[[length(p_list)]]
#--














