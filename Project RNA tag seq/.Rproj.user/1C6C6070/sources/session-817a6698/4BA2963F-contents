
setwd("Z:/DelvilleLab/Kevin Lab/Project RNA tag seq/limma example")


#### limma example from guide 

library(limma)
library(Glimma)
library(edgeR)
library(tidyverse)
library(scales)


library(tidyverse)
library(DESeq2)
library(pheatmap)
library(annotables)
library(Mus.musculus)


# url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
# utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
# utils::untar("GSE63310_RAW.tar", exdir = ".")
# files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
#            "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
#            "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
# for(i in paste(files, ".gz", sep=""))
#   R.utils::gunzip(i, overwrite=TRUE)

files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
           "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
           "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
           "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)

x <- readDGE(files, columns=c(1,3))
class(x)
dim(x)

samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames

colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples


geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)

genes <- genes[!duplicated(genes$ENTREZID),]

x$genes <- genes
x

# 5 Data pre-processing

# 5.1 Transformations from the raw-scale
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)


L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)


summary(lcpm)


# filtering low counts

table(rowSums(x$counts==0)==9)


keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)



lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")



# 5.3 Normalising gene expression distributions

x <- calcNormFactors(x, method = "TMM") #trimmed mean of M-values (TMM)
x$samples$norm.factors


x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors


lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")



# 5.4 Unsupervised clustering of samples

# one of the most important exploratory plots to examine for gene expression analyses is the multi-dimensional scaling (MDS) plot

#seems similar to PCA

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")


glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=F) #change to T to explore 


#### 6 Differential expression analysis

#### more or less, all of the above is stuff i explored with DESeq. Here is where we can start, using our two matrices of data


# 6.1 Creating a design matrix and contrasts

design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design


contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP, 
  BasalvsML = Basal - ML, 
  LPvsML = LP - ML, 
  levels = colnames(design))
contr.matrix


# 6.2 Removing heteroscedascity from count data

par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v


vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")


# 6.3 Fitting linear models for comparisons of interest --- It can be seen from this plot that the variance is no longer dependent on the mean expression level.

# 6.4 Examining the number of DE genes

summary(decideTests(efit))

tfit <- treat(vfit, lfc=1) # get log-FCs
dt <- decideTests(tfit)
summary(dt)


# Genes that are DE in multiple comparisons can be extracted using the results from decideTests, where 0s represent genes that are not DE, 1s represent genes that are up-regulated, and -1s represent genes that are down-regulated.
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)

vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

write.fit(tfit, dt, file="results.txt")


# 6.5 Examining individual DE genes from top to bottom

basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)

head(basal.vs.lp)
head(basal.vs.ml)


# 6.6 Useful graphical representations of differential expression results

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=F) #change to T to explore



library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column") # Heatmap of log-CPM values for top 100 genes DE in basal versus LP










