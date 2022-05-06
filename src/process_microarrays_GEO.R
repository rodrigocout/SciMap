library(GEOquery)
library(limma)
library(affy)

# load series and platform data from GEO
gset <- getGEO("GSE96036", GSEMatrix =TRUE)
#Check
show(gset)
show(pData(phenoData(gset[[1]]))[1:5,c(1,6,10)])
#---
#Make a phenotype table
Pheno <- pData(phenoData(gset[[1]]))[,c(1,6,10)]
#-----------
#Process
if (length(gset) > 1) idx <- grep("GPL21185", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#Make a gene expression matrix
gexp <- exprs(gset)
#Check 
head(gexp)
# log2 transform
Log2gexp <- log2(gexp)

#Normalize the log2 transform dataset
NormGexp <- normalizeBetweenArrays(Log2gexp) 

#------------
#Probe to gene
genes <- rownames(gexp)
library("biomaRt")
ensembl = useMart(biomart= "ensembl",dataset="hsapiens_gene_ensembl")

#Check the microarray platform
tables <- listAttributes(ensembl)
tables[grep('agilent', tables[,1]),]

#Make a data frame
agilent.df<-getBM(attributes = c("hgnc_symbol","external_gene_name", "agilent_sureprint_g3_ge_8x60k_v2"), 
                  filters=c("agilent_sureprint_g3_ge_8x60k_v2"),values=genes, mart=ensembl)
#check 
head(agilent.df)
#Make probes rownames
rownames(agilent.df) <- make.names(agilent.df$agilent_sureprint_g3_ge_8x60k_v2, unique=TRUE)

#Merge with the IDs
Gexp <- merge(agilent.df, gexp, by = "row.names")
#Check
head(Gexp)
#Process
colnames(Gexp)[1] <- c('agilent_probes')
Gexp$agilent_sureprint_g3_ge_8x60k_v2 <- NULL

#----------
Log2NormGexp <- merge(agilent.df, NormGexp, by = "row.names") #Here I changed to NormGexp
head(Log2Gexp)
#Process
colnames(Log2NormGexp)[1] <- c('agilent_probes')
Log2NormGexp$agilent_sureprint_g3_ge_8x60k_v2 <- NULL




write.table(Gexp, "Raw_expression_Dataset_GSE96036_Kawataetal2019.txt", 
            sep = '\t', quote = F, row.names = F)

write.table(Log2NormGexp, "Log2Normalized_expression_Dataset_GSE96036_Kawataetal2019.txt", 
            sep = '\t', quote = F, row.names = F)

write.table(Pheno,"Phenotype_Dataset_GSE96036_Kawataetal2019.txt", 
            sep = '\t', quote = F)
