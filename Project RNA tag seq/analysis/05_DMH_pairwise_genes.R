
# pairwise comparisons of normalized expression of important genes 
# selected by large logFC, MEcolor strength, and theoretical importance 

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(CorLevelPlot)
library(gridExtra)
library(clusterProfiler)
library(enrichplot)
library(biomaRt)
library(AnnotationDbi)
library(annotables)
grcm38 <- grcm38
library(ggpubr)
library(tidyverse)


#### making normalized counts - shouldn't need to be done again ####
data <- read.csv("ham_brain_data/DMH_counts.csv")

phenoData <- read.csv("ham_brain_data/DMH_id.csv")

data[1:10, 1:10]
head(phenoData)

# prepare data

data <- data %>% 
  gather(key = "samples", value = "counts", -X) %>% 
  rename(gene = X) %>% 
  inner_join(., phenoData, by = c("samples" = "X")) %>% 
  select(1, 3, 4) %>% 
  spread(key = "subject", value = "counts") %>% 
  column_to_rownames(var = "gene") %>%
  slice(-1:-5)


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK


table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detected as outliers
data <- data[gsg$goodGenes == TRUE,]


# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)
# no outliers removed???? -- wait -- come back after retesting DMH WGCNA with some "outliers" removed


# # pca - method 2 for finding outliers 
# 
# pca <- prcomp(t(data))
# pca.dat <- pca$x
# 
# pca.var <- pca$sdev^2
# pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
# 
# pca.dat <- as.data.frame(pca.dat)
# 
# ggplot(pca.dat, aes(PC1, PC2)) +
#   geom_point() +
#   geom_text(label = rownames(pca.dat)) +
#   labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
#        y = paste0('PC2: ', pca.var.percent[2], ' %'))
# #KM193 again seems like a potential outlier 


# exclude outlier samples
samples.to.be.excluded <- c('KM193')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
### ***  note in results that we exclude KM193 because of higher variance compared to all other samples 



# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

phenoData <- phenoData %>% 
  column_to_rownames(var = "subject")

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# fixing column names in colData
names(colData)

# selecting relevant info
colData <- colData %>% 
  select(3, 10, 11, 13, 14, 15, 20, 22)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))



# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not specifying model "because we need this DSeq data set to perform variance stabilizing transformation"


## remove all genes with counts < 10 in more than 3 samples -- here, same as DEseq2 guide

dds75 <- dds[rowSums(counts(dds) >= 10) >= 3,]
nrow(dds75) #


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

norm.counts2 <- norm.counts
norm.counts2 <- as.data.frame(norm.counts2)
norm.counts3 <- tibble::rownames_to_column(norm.counts2, "subject")
colData2 <- tibble::rownames_to_column(colData, "subject")
normcounts.coldata.DMH <- cbind(colData2, norm.counts3)

write.csv(normcounts.coldata.DMH, "results/normcounts_coldata_DMH.csv")


#### shouldn't need to do above again ####

normcounts.coldata.DMH <- read.csv("results/normcounts_coldata_DMH.csv")
normcounts.coldata.DMH$X <- NULL
DMH_limma_results1 <- readRDS("results/DMH_limma_results.RDS")
DMH_limma_results1 <- DMH_limma_results1 %>% 
  select(1,2,3,9,10)