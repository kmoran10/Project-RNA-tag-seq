
library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(tidyverse)
grcm38 # mouse genes


source("functions/gettop10GO.R")


my_logFC_threshold = 0.2

y1a <- readRDS("results/LH_limma_results.RDS")


gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "Control - Stress") -> top10_GOterms_LH


write.csv(top10_GOterms_LH,"results/results_tables/top100_GOterms_LH.csv", row.names = F)




y2a <- readRDS("results/DMH_limma_results.RDS")

gettop10GO(y2a, my_showCategory) %>% 
  mutate(comparison = "Control - Stress") -> top10_GOterms_DMH


write.csv(top10_GOterms_DMH,"results/results_tables/top100_GOterms_DMH.csv", row.names = F)





y3a <- readRDS("results/ARC_limma_results.RDS")

gettop10GO(y3a, my_showCategory) %>% 
  mutate(comparison = "Control - Stress") -> top10_GOterms_ARC


write.csv(top10_GOterms_ARC,"results/results_tables/top100_GOterms_ARC.csv", row.names = F)





#### 
library(tidyverse)

LH <- readRDS("results/LH_limma_results.RDS")
DMH <- readRDS("results/DMH_limma_results.RDS")
ARC <- readRDS("results/ARC_limma_results.RDS")

LH  %>% filter(symbol == "Hcrt")
DMH %>% filter(symbol == "Npvf")
ARC %>% filter(symbol == "Npy")
