
## following the DESeq2 guide for all data together


library(tidyverse)
library(DESeq2)


cts_with_unassigned <- read.csv("ham_brain_data/all_counts.csv", row.names = "X")
cts <- cts_with_unassigned[-c(32317:32345),]
# cts <- cts[,-1] # not removing 193 just yet


coldata <- read.csv("ham_brain_data/all_id.csv", row.names = "X")
# coldata <- coldata[-1,] # not removing 193 just yet

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))



#### ROI checking first - do all samples from all ROIs go with their ROIs

dds <- DESeqDataSetFromMatrix(countData = cts,
                          colData = coldata,
                          design= ~ ROI)


dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients


#removing 0s
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

res <- results(dds, name = "group_Stress_vs_Control")

res

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="group_Stress_vs_Control", type="apeglm")
resLFC

resOrdered <- res[order(res$pvalue),]
resOrdered

summary(res)

res05 <- results(dds, alpha=0.05)
summary(res05)

library(IHW)

resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult



plotMA(res, ylim=c(-2,2))

plotMA(resLFC, ylim=c(-2,2))

# USE BELOW idx fxn to "pluck" and identify genes from plot by clicking on them in plot

# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]


#Alternative shrinkage estimators

resNorm <- lfcShrink(dds, coef=2, type="normal")
library(ashr)
resAsh <- lfcShrink(dds, coef=2, type="ashr")


par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")


dev.off()

plotCounts(dds, gene=which.min(res$padj), intgroup="group")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="group", 
                returnData=TRUE)
ggplot(d, aes(x=group, y=count)) + 
  geom_boxplot() +
  geom_point(position=position_jitter(w=0.1,h=0))


mcols(res)$description


write.csv(as.data.frame(resOrdered), 
          file="results/ARC_bygroup_DESeq_results.csv")




#### Count data transformations


#variance stabilizing transformations
vsd <- vst(dds, blind=FALSE)
#regularized logarithm
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)


# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))


#### Data quality assessment by sample clustering and visualization

# Heatmap of the count matrix

coldata$group <- as.factor(coldata$group)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("group","mother")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# not sure what use this is.

## Heatmap of the sample-to-sample distances

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# not sure what use this is either..... 


### Principal component plot of the samples (PCA)

plotPCA(vsd, intgroup=c("group"))

mat <- assay(vsd)
mm <- model.matrix(~group, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup=c("group"))



pcaData <- plotPCA(vsd, intgroup=c("ROI"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=ROI, label = rownames(pcaData))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + geom_text() + labs(title="ROI")








