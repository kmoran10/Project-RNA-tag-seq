### PCA HAM!

library(limma)
# library(Mus.musculus)  --- prob useless
library(DESeq2)
library(edgeR)
library(tidyverse)


## Getting PCA from the DEG results not LIMMA

# Expression values
dlNorm <-  read.csv("ham_brain_data/LH_counts.csv", row.names = 1)

#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]



