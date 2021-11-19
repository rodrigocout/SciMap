
# Load data ----
library(dplyr)
library(Seurat)
library(patchwork)

ind <- 5

data_list <- list.dirs(path = "C:/Data/phd/code/scimap/SciMap/data", full.names = TRUE, recursive = FALSE)
data_names <- basename(data_list)

rm_data <- c()

for(k in 1:length(data_names)) {
  nm <- data_names[k]
  nm_start <- substr(nm,1,1)
  if(nm_start != "D") {
    rm_data <- c(rm_data, k)
  }
}

data_list <- data_list[- rm_data]
data_names <- data_names[- rm_data]

# choose dataset
data_dir <- data_list[ind]

pbmc.data <- Read10X(data.dir = data_dir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k")

# get percentage of mitochondrial data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


# trim data
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# extract features ----

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 7)

# normalization, normalized data is in pbmc[["RNA"]]@data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# find features with high variance accross cells
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# scale data, pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# infl. of PCAs
pbmc <- JackStraw(pbmc, num.replicate = 100, dims = 30)
pbmc <- ScoreJackStraw(pbmc, dims = 1:30)
JackStrawPlot(pbmc, dims = 1:30)

# plot results ----

dm <- 20
# find clusters

pbmc <- FindNeighbors(pbmc, dims = 1:dm)

pbmc <- FindClusters(pbmc, resolution = 0.5)

head(Idents(pbmc), 5)

# UMAP analysis

pbmc <- RunUMAP(pbmc, dims = 1:dm)

# find markers

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
m_i_genes <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)

m_i_genes <- m_i_genes[7]
m_i_genes <- m_i_genes[[1]]

new.cluster.ids <- m_i_genes
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()