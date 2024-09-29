
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


# pca - method 2 for finding outliers

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
#KM193, 194, 197  outliers - not removing for normalization (??)


# exclude outlier samples
# samples.to.be.excluded <- c('KM193','KM194', 'KM197')
samples.to.be.excluded <- c(NA)
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

dds75 <- dds[rowSums(counts(dds) >= 5) >= 3,]
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
DMH_MEdarkred_limma <- read.csv("results/DMH_MEdarkred_limma.csv") #LOWER expression in Stress
DMH_MEdarkred_limma$X <- NULL
DMH_MEskyblue_limma <- read.csv("results/DMH_MEskyblue_limma.csv") #HIGHER expression in Stress
DMH_MEskyblue_limma$X <- NULL
DMH_MEdarkgrey_limma <- read.csv("results/DMH_MEdarkgrey_limma.csv") #HIGHER expression in Stress
DMH_MEdarkgrey_limma$X <- NULL




### top 50 biggest logFC genes 
DMH_limma_results1 %>% 
  arrange(logFC) %>% 
  head(25)

DMH_limma_results1 %>% 
  arrange(-logFC) %>% 
  head(25)




# visualize some top logFCs

# Gabpb2 - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, Gabpb2, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="C. DMH - Gabpb2",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 


# Izumo1 - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, Izumo1, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="D. DMH - Izumo1",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



# Npw - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, Npw, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="E. DMH - Npw",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



# Lrp2 - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, Lrp2, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="F. DMH - Lrp2",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 





### High MM genes from MEs (c-H)

# Arhgdig - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, Arhgdig, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="C. DMH - Arhgdig",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



# Rpl36 - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, Rpl36, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="D. DMH - Rpl36",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 


# Cntnap1 - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, Cntnap1, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="E. DMH - Cntnap1",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 


# Fgfr1 - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, Fgfr1, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="F. DMH - Fgfr1",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 


# Fcer1g - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, Fcer1g, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="G. DMH - Fcer1g",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 


# C1qc - save as w400 h600
ggplot(normcounts.coldata.DMH, aes(group, C1qc, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="H. DMH - C1qc",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 





