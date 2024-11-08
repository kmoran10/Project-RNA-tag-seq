
# WGCNA for ARC --------- TESTING RESULTS IF REMOVE "OUTLIER-Y" SAMPLES
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

data <- read.csv("ham_brain_data/ARC_counts.csv")

phenoData <- read.csv("ham_brain_data/ARC_id.csv")

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
# this is a weird one. two somewhat distinct trees. But maybe 196 and 207 are outliers?



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
# here, seems like 210, 193, and maybe 196?


#so 196 out for sure
# but we'll start with 196, 193, and 210 out

# exclude outlier samples
samples.to.be.excluded <- c("KM193","KM196","KM210")
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
### ***  note removed outliers 



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
nrow(dds75) # 9193 genes


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
##### THIS PART NEEDS TO ADJUST PER ROI ANALYSIS 
# 6-10 could work. even as low as 4?
  # start with 6 - 3 sig diff modules out of 11
  # check 10 - 2 sig dif MEs out of 12. 
  # hm, lets try 4
  




# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 6  ### SELECTED BY ABOVE 
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
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

#View(module.trait.corr.pvals)


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
# power 6 - 11 MEs? some def too big tho? -- BUT 3 are sig diff? 
# power 10 - 12 MEs - only 2 sig difs?
# power 4 only 13 MEs??? - 2 sig dif.
# stick with 6 i guess



### NEED TO adjust this step depending on dimensions of heatmap.data

## IN THE FOLLOWING: POSITIVE VALUES MEAN ME EXPRESSION IS HIGHER IN TRAIT CODED WITH 1 COMPARED TO TRAIT CODED WITH 0 - specifically positive = higher in stressed
### SAVE 550x800
CorLevelPlot(heatmap.data2,
             x = names(heatmap.data2)[12:12], #trait data
             y = names(heatmap.data2)[1:11], #ME data
             col = c("blue3", "skyblue", "white", "#f7ab5e", "orange2"),
             main = "E. ARC -
WGCNA Module Eigengenes",
             fontCorval = 1,
             cexCorval = 1.6,
             cexLabY = 1.5) 

## so for ARC, modules yellow, purple, and blue are significantly altered by stress

module.gene.mapping <- as.data.frame(bwnet$colors)




## FROM HERE, CAN PERFORM FURTHER ANALYSIS ON GENES IN IMPORTANT MODULES - the following lists the genes in relevant modules that are selected from above analysis to be modules that are differentially expressed between groups


## THIS ME has LOWER expression in Stress
module.gene.mapping %>%
  filter(`bwnet$colors` == 'yellow') %>%
  rownames()

module.gene.mapping %>%
  filter(`bwnet$colors` == 'purple') %>%
  rownames

## THIS ME has HIGHER expression in Stress
module.gene.mapping %>%
  filter(`bwnet$colors` == 'blue') %>%
  rownames()



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
#top 25 genes in the ARC sig associated with stress experience
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

#1st - pull ARC_limma_results1 <- readRDS("results/ARC_limma_results.RDS")
#2nd - filter genes that are only in modules of interest - call them MEcolor
#3rd - attach ARC_limma_results1 to filteres MEcolor datasets
#4th - do GO analysis of merged MEcolor-limma datasets
#5th - and highest MM ranking of MEcolor datasets


## steps 1-3
ARC_limma_results1 <- readRDS("results/ARC_limma_results.RDS")

gene.signf.corr2 <- as.data.frame(gene.signf.corr)
gene.signf.corr2 <- tibble::rownames_to_column(gene.signf.corr2, "symbol")
colnames(gene.signf.corr2)[2] <- "gene.signif.corr"

gene.signf.corr.pvals2 <- as.data.frame(gene.signf.corr.pvals)
gene.signf.corr.pvals2 <- tibble::rownames_to_column(gene.signf.corr.pvals2, "symbol")
colnames(gene.signf.corr.pvals2)[2] <- "gene.signif.corr.pval"

MEallcolors2 <- left_join(MEallcolors, gene.signf.corr2, by = "symbol") %>%
  left_join(., gene.signf.corr.pvals2, by = "symbol")

write.csv(MEallcolors2, "results/MEallcolors_ARC.csv")

module.membership.measure2 <- as.data.frame(module.membership.measure)
module.membership.measure2 <- tibble::rownames_to_column(module.membership.measure2, "module")



### down in stressed  yellow, purple
### up in stressed    blue



MEyellow <- MEallcolors2 %>%
  filter(MEcolor == "yellow")

MEyellow.limma <- ARC_limma_results1 %>%
  left_join(MEyellow, by = "symbol") %>%
  filter(!is.na(MEcolor)) %>%
  select(1,2,3,4,10,11,12,13)

write.csv(MEyellow.limma, "results/ARC_MEyellow_limma.csv")


MEpurple <- MEallcolors2 %>%
  filter(MEcolor == "purple")

MEpurple.limma <- ARC_limma_results1 %>%
  left_join(MEpurple, by = "symbol") %>%
  filter(!is.na(MEcolor)) %>%
  select(1,2,3,4,10,11,12,13)

write.csv(MEpurple.limma, "results/ARC_MEpurple_limma.csv")


MEblue <- MEallcolors2 %>%
  filter(MEcolor == "blue")

MEblue.limma <- ARC_limma_results1 %>%
  left_join(MEblue, by = "symbol") %>%
  filter(!is.na(MEcolor)) %>%
  select(1,2,3,4,10,11,12,13)

write.csv(MEblue.limma, "results/ARC_MEblue_limma.csv")




## GO ANALYSIS OF yellow MODULE
gettop10GO(MEyellow.limma, my_showCategory) %>%
  mutate(comparison = "Control - Stress") -> GOterms_ARC_yellow

View(GOterms_ARC_yellow)

write.csv(GOterms_ARC_yellow, "results/GOterms_ARC_yellow.csv")

## HIGHEST MM OF yellow MODULE
MMyellow <- module.membership.measure2 %>% 
  filter(module == "MEyellow") %>% 
  gather(., gene, MM)
MMyellow <- MMyellow[-1,]
MMyellow %>% 
  arrange(desc(MM)) %>% 
  head(10)




## GO ANALYSIS OF purple MODULE
gettop10GO(MEpurple.limma, my_showCategory) %>%
  mutate(comparison = "Control - Stress") -> GOterms_ARC_purple

View(GOterms_ARC_purple)

write.csv(GOterms_ARC_purple, "results/GOterms_ARC_purple.csv")

## HIGHEST MM OF purple MODULE
MMpurple <- module.membership.measure2 %>% 
  filter(module == "MEpurple") %>% 
  gather(., gene, MM)
MMpurple <- MMpurple[-1,]
MMpurple %>% 
  arrange(desc(MM)) %>% 
  head(10)




## GO ANALYSIS OF blue MODULE
gettop10GO(MEblue.limma, my_showCategory) %>%
  mutate(comparison = "Control - Stress") -> GOterms_ARC_blue

View(GOterms_ARC_blue)

write.csv(GOterms_ARC_blue, "results/GOterms_ARC_blue.csv")

## HIGHEST MM OF blue MODULE
MMblue <- module.membership.measure2 %>% 
  filter(module == "MEblue") %>% 
  gather(., gene, MM)
MMblue <- MMblue[-1,]
MMblue %>% 
  arrange(desc(MM)) %>% 
  head(10)
