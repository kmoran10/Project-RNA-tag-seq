
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
data <- read.csv("ham_brain_data/LH_counts.csv")

phenoData <- read.csv("ham_brain_data/LH_id.csv")

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
# plot(htree)
# KM193 potentially seems like an outlier sample


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
normcounts.coldata.LH <- cbind(colData2, norm.counts3)

write.csv(normcounts.coldata.LH, "results/normcounts_coldata_LH.csv")


#### shouldn't need to do above again ####

normcounts.coldata.LH <- read.csv("results/normcounts_coldata_LH.csv")
normcounts.coldata.LH$X <- NULL
LH_limma_results1 <- readRDS("results/LH_limma_results.RDS")
LH_limma_results1 <- LH_limma_results1 %>% 
  select(1,2,3,9,10)
LH_MEyellow_limma <- read.csv("results/LH_MEyellow_limma.csv") #LOWER expression in Stress
LH_MEyellow_limma$X <- NULL
LH_MEmagenta_limma <- read.csv("results/LH_MEmagenta_limma.csv") #HIGHER expression in Stress
LH_MEmagenta_limma$X <- NULL


### top 50 biggest logFC genes 
LH_limma_results1 %>% 
  arrange(logFC) %>% 
  head(25)

LH_limma_results1 %>% 
  arrange(-logFC) %>% 
  head(25)

#visualize some favorites 

ggplot(normcounts.coldata.LH, aes(x=group,y=Hcrt, fill=group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = .1)) +
  stat_compare_means(method = "t.test") + 
  ggtitle("Hcrt")


# hcrt - save as w400 h600
ggplot(normcounts.coldata.LH, aes(group, Hcrt, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="C. LH - Hcrt",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 





# ggplot(normcounts.coldata.LH, aes(x=Hcrt,y=wt.rate.pct, color=group)) +
#   geom_smooth(method="lm") +
#   geom_point()
# 
# ggplot(normcounts.coldata.LH, aes(x=Hcrt,y=food.eff, color=group)) +
#   geom_smooth(method="lm") +
#   geom_point()
# 
# ggplot(normcounts.coldata.LH, aes(x=Hcrt,y=food.total, color=group)) +
#   geom_smooth(method="lm") +
#   geom_point()
# 
# ggplot(normcounts.coldata.LH, aes(x=Hcrt,y=wt.day.15, color=group)) +
#   geom_smooth(method="lm") +
#   geom_point()
# 
# ggplot(normcounts.coldata.LH, aes(x=Hcrt,y=total.fat, color=group)) +
#   geom_smooth(method="lm") +
#   geom_point()
# 
# ggplot(normcounts.coldata.LH, aes(x=Hcrt,y=g.hoarded, color=group)) +
#   geom_smooth(method="lm") +
#   geom_point()
# 
# ggscatter(normcounts.coldata.LH, x="Hcrt", y="wt.rate.pct", add="reg.line")+
#   stat_cor(label.x = 8, size = 6) +
#   scale_fill_manual(values=c("white", "darkgray")) +
#   labs(title="Hcrt & rate of weight gain") +
#   theme_classic() +
#   theme(axis.line = element_line(colour = 'black', size = 1),
#         axis.ticks = element_line(colour = "black", size = 1),
#         legend.position="none",
#         axis.text=element_text(size=16),
#         axis.title=element_text(size=18,face="bold"),
#         plot.title = element_text(size=24, hjust = 0.4)) 




# Mc3r - save as w400 h600
ggplot(normcounts.coldata.LH, aes(group, Mc3r, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="D. LH - Mc3r",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 




# Npvf - save as w400 h600
ggplot(normcounts.coldata.LH, aes(group, Npvf, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="E. LH - Npvf",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



# Gpr50 - save as w400 h600
ggplot(normcounts.coldata.LH, aes(group, Gpr50, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="F. LH - Gpr50",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 




#### Pairwise of High ME MM genes - yellow


# visualize 

ggplot(normcounts.coldata.LH, aes(group, Mog, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="C. LH - Mog",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



ggplot(normcounts.coldata.LH, aes(group, Cldn11, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="D. LH - Cldn11",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 




#### Pairwise of High ME MM genes - magenta

# visualize 

ggplot(normcounts.coldata.LH, aes(group, Cartpt  , fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="E. LH - Cartpt",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



ggplot(normcounts.coldata.LH, aes(group, Basp1, fill=group)) +
  stat_compare_means(method = "t.test", size = 6, label.x = 0.8) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size = 4, position=position_jitter(0.15)) + 
  scale_fill_manual(values=c("blue3", "orange2")) +
  labs(title="F. LH - Basp1",x="Group", y = "Normalized gene expression") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 




