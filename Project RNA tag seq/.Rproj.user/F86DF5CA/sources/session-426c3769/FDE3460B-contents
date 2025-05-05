
### PCA

library(limma)
library(Mus.musculus)
library(DESeq2)
library(edgeR)
library(tidyverse)


## Getting PCA from the DEG results not LIMMA

# Expression values
dlNorm <-  read.csv("ms_brain_data/LH_counts.csv", row.names = 1)

#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]


## only 20ish genes... VAST majority seem unmapped.  
## so, something bad is happening during mapping. back to the pod
