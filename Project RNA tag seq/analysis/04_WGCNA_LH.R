
# WGCNA for LH
# BiocManager::install("")

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
library(tidyverse)
source("functions/gettop10GO.R")


allowWGCNAThreads()          # allow multi-threading (optional)

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
plot(htree)
# KM193 potentially seems like an outlier sample


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
#KM193 again seems like a potential outlier 


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


## remove all genes with counts < 15 in more than 75% of samples (19*0.75=14.25)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 15,]
nrow(dds75) # 8409 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()




# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
power

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)
# selecting a power of 12 since it seems to maximize around here and has minimal mean connectivity
# 4 or 5 is where it crossed over R^2 of .8, but mean connectivity was still relatively high.
# 10 could potentially work if this R^2 is "excessively high" for some reason. 



# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 12 
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000, #suitable for 16gigs of RAM in PC
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25, #threshold of merging simmilar modules 
                          numericLabels = FALSE, #to set model eigengene labels as colors 
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor




# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)




# grey module = all genes that doesn't fall into other modules were assigned to the grey module





# 6A. Relate modules to traits --------------------------------------------------
# module trait associations



# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(stress_bin = ifelse(grepl('Stress', group), 1, 0)) %>% 
  select(9)


# binarize categorical variables 

# colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))
# 
# severity.out <- binarizeCategoricalColumns(colData$severity,
#                                            includePairwise = FALSE,
#                                            includeLevelVsAll = TRUE,
#                                            minCount = 1)
# 
# 
# traits <- cbind(traits, severity.out)
### NOT DOING ABOVE, NOT RELEVANT

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
module.trait.corr.pvals

# View(module.trait.corr.pvals)


# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')



### NEED TO adjust this step depending on dimensions of heatmap.data
names(heatmap.data)

heatmap.data2 <- heatmap.data
module.gene.mapping <- as.data.frame(bwnet$colors)
MEallcolors <- tibble::rownames_to_column(module.gene.mapping, "symbol")
colnames(MEallcolors)[2] <- "MEcolor"
color_counts <- table(MEallcolors$MEcolor)
column_names <- colnames(heatmap.data2)

new_column_names <- sapply(column_names, function(col) {
  # Extract the color part (e.g., "MEblack" -> "black")
  color <- sub("ME", "", col)
  
  # If the color exists in the color_counts, add the count to the column name
  if (color %in% names(color_counts)) {
    paste0(col, "(", color_counts[color], ")")
  } else {
    col # Leave the column name unchanged if the color is not found
  }
})

colnames(heatmap.data2) <- new_column_names

names(heatmap.data2)

## IN THE FOLLOWING: POSITIVE VALUES MEAN ME EXPRESSION IN HIGHER IN TRAIT CODED WITH 1 COMPARED TO TRAIT CODED WITH 0 
### SAVE 550x700

CorLevelPlot(heatmap.data2,
             x = names(heatmap.data2)[15:15], #trait data
             y = names(heatmap.data2)[1:14], #ME data
             col = c("blue3", "skyblue", "white", "#f7ab5e", "orange2"),
             main = "A. LH - WGCNA 
Module Eigengenes",
             fontCorval = 1,
             cexCorval = 1.6,
             cexLabY = 1.5)

## so for LH, modules magenta and yellow are significantly altered by stress

module.gene.mapping <- as.data.frame(bwnet$colors)

## FROM HERE, CAN PERFORM FURTHER ANALYSIS ON GENES IN IMPORTANT MODULES - the following lists the genes in relevant modules that are selected from above analysis to be modules that are differentially expressed between groups


## THIS ME has LOWER expression in Stress
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames()

lmy <- lm(MEyellow ~ stress_bin, heatmap.data)
summary(lmy)
  


## THIS ME has HIGHER expression in Stress
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'magenta') %>% 
  rownames()

lmm <- lm(MEmagenta ~ stress_bin, heatmap.data)
summary(lmm)


# 6B. Intramodular analysis: Identifying driver genes ---------------
#"highly connected intramodular hub genes"


# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure[1:10,1:10]
module.membership.measure.pvals[1:10,1:10] # just a check that this ran right


# Calculate the gene significance and associated p-values
#correlate expression data with trait of interest - FOR ME, JUST stress_bin

gene.signf.corr <- cor(norm.counts, traits$stress_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)
#top 25 genes in the LH sig associated with stress experience 
### NEED TO DO SOMETHING SIMILAR TO THIS *JUST* WITHIN SIG MODULES - GET "HIGHEST MM GENES"   ### basically take module.membership.measure.pvals, flip orientation, filter only relevant module, then arrange(V1) 



# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.




#### FOLLOWUP ANALYSIS OF VARIOUS THINGS FOUND ABOVE 


library(clusterProfiler)
library(enrichplot)
library(biomaRt)
library(AnnotationDbi)
library(annotables)
grcm38 <- grcm38
source("functions/gettop10GO.R")

#1st - pull LH_limma_results1 <- readRDS("results/LH_limma_results.RDS")
#2nd - filter genes that are only in modules of interest - call them MEcolor
#3rd - attach LH_limma_results1 to filteres MEcolor datasets
#4th - do GO analysis of merged MEcolor-limma datasets 
#5th - and highest MM ranking of MEcolor datasets


## steps 1-3
LH_limma_results1 <- readRDS("results/LH_limma_results.RDS")
MEallcolors <- tibble::rownames_to_column(module.gene.mapping, "symbol")
colnames(MEallcolors)[2] <- "MEcolor"

gene.signf.corr2 <- as.data.frame(gene.signf.corr)
gene.signf.corr2 <- tibble::rownames_to_column(gene.signf.corr2, "symbol")
colnames(gene.signf.corr2)[2] <- "gene.signif.corr"

gene.signf.corr.pvals2 <- as.data.frame(gene.signf.corr.pvals)
gene.signf.corr.pvals2 <- tibble::rownames_to_column(gene.signf.corr.pvals2, "symbol")
colnames(gene.signf.corr.pvals2)[2] <- "gene.signif.corr.pval"

MEallcolors2 <- left_join(MEallcolors, gene.signf.corr2, by = "symbol") %>% 
  left_join(., gene.signf.corr.pvals2, by = "symbol")

write.csv(MEallcolors2, "results/MEallcolors_LH.csv")


MEyellow <- MEallcolors2 %>% 
  filter(MEcolor == "yellow")

MEyellow.limma <- LH_limma_results1 %>% 
  left_join(MEyellow, by = "symbol") %>%
  filter(!is.na(MEcolor)) %>% 
  select(1,2,3,4,10,11,12,13)

write.csv(MEyellow.limma, "results/LH_MEyellow_limma.csv")


MEmagenta <- MEallcolors2 %>% 
  filter(MEcolor == "magenta")

MEmagenta.limma <- LH_limma_results1 %>% 
  left_join(MEmagenta, by = "symbol") %>%
  filter(!is.na(MEcolor)) %>% 
  select(1,2,3,4,10,11,12,13)

write.csv(MEmagenta.limma, "results/LH_MEmagenta_limma.csv")

## GO ANALYSIS OF YELLOW MODULE
gettop10GO(MEyellow.limma, my_showCategory) %>% 
  mutate(comparison = "Control - Stress") -> GOterms_LH_yellow

write.csv(GOterms_LH_yellow, "results/GOterms_LH_yellow.csv")

## HIGHEST MM OF YELLOW MODULE
MEyellow.limma %>% 
  arrange(gene.signif.corr.pval) %>% 
  head(5)




## GO ANALYSIS OF MAGENTA MODULE
gettop10GO(MEmagenta.limma, my_showCategory) %>% 
  mutate(comparison = "Control - Stress") -> GOterms_LH_magenta

write.csv(GOterms_LH_magenta, "results/GOterms_LH_magenta.csv")

## HIGHEST MM OF MAGENTA MODULE
MEmagenta.limma %>% 
  arrange(gene.signif.corr.pval) %>% 
  head(5)





#### Module boxplots 
source("functions/geom_boxjitter.R")



MEdf <- as.data.frame(module_eigengenes) %>% 
  tibble::rownames_to_column(., "subject")

cdgrp <- as.data.frame(colData) %>% 
  tibble::rownames_to_column(., "subject") %>% 
  select(subject, group)

MEdf2 <- left_join(MEdf, cdgrp, by="subject")

# yellow  - save 325x500
MEdf2 %>%
  ggplot(aes(group,MEyellow, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  scale_fill_manual(values=c("blue3", "darkorange")) +
  labs(title="Yellow:
146 genes",x="Group", y = "Module Eigengene") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
                legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



# magenta 325 x 500
MEdf2 %>%
  ggplot(aes(group,MEmagenta, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  scale_fill_manual(values=c("blue3", "darkorange")) +
  labs(title="Magenta:
37 genes",x="Group", y = "Module Eigengene") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 




