
## following the DESeq2 guide for my ARC data
## removing KM193 because outlier and was noted as problematic in punching. 


library(tidyverse)
library(DESeq2)


cts_with_unassigned <- read.csv("ham_brain_data/ARC_counts.csv", row.names = "X")
cts <- cts_with_unassigned[-c(32317:32345),]
cts <- cts[,-1]

coldata <- read.csv("ham_brain_data/ARC_id.csv", row.names = "X")
coldata <- coldata[-1,]

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))


### "Quick Start"

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ group)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name = "group_Stress_vs_Control")
# or to shrink log fold changes association with condition:
# res.s <- lfcShrink(dds, coef="group_Stress_vs_Control", type="apeglm")
# lots of warning/errors


### Deeper

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



pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group, label = rownames(pcaData))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + geom_text() + labs(title="ARC")




### more about shrinkage estimators

resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)


resAshT <- lfcShrink(dds, coef=2, type="ashr", lfcThreshold=1)
plotMA(resAshT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)



### Tests of log2 fold change above or below a threshold

par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()











