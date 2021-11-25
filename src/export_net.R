# Load data ----
library(dplyr)
library(Seurat)
library(patchwork)
library(R.matlab)
library(WGCNA)
library(data.table)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

ind <- 1
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
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "ipsc")

# get percentage of mitochondrial data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


# trim data
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# extract features ----

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 7)

# normalization, normalized data is in pbmc[["RNA"]]@data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# find gene names
all.genes <- rownames(pbmc)

genes <- list(all.genes)

# save data
writeMat("C:/Data/phd/code/D1_data.mat",x=pbmc[["RNA"]]@data)

write.csv(genes, "C:/Data/phd/code/genes_D1.csv")
