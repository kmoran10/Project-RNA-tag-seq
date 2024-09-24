

#### limma 

library(limma)
library(Glimma)
library(edgeR)

library(biomaRt)
library(AnnotationDbi)
library(annotables)
grcm38 <- grcm38
library(tidyverse)


cts_with_unassigned <- read.csv("ham_brain_data/ARC_counts.csv", row.names = "X")
dlNorm  <- cts_with_unassigned[-c(32317:32345),]
dlNorm  <- dlNorm [,-1]

coldata <- read.csv("ham_brain_data/ARC_id.csv", row.names = "X")
coldata <- coldata[-1,]

dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]

dlNorm  <- dlNorm [, rownames(coldata)]
all(rownames(coldata) == colnames(dlNorm))


####

# normalize all coldata values

columns_to_normalize <- 10:22

normalize_z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

for (col in columns_to_normalize) {
  col_name <- names(coldata)[col]
  new_col_name <- paste0(col_name, ".Z")
  
  coldata[new_col_name] <- normalize_z_score(coldata[, col])
}


####

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0 = DGEList(d, group = coldata$group)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)

dge.dl$samples$group

coldata %>% 
  dplyr::select(subject, group) -> var_info



coldata$group %>% 
  factor(.,levels = c("Control","Stress")) -> group.dl

design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames


v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)


contr.matrix <- makeContrasts(group.dlControl-group.dlStress,
                              levels = design.dl)


vfit.dl2 <- contrasts.fit(vfit.dl, contr.matrix)

efit.dl2 = eBayes(vfit.dl2)

p.dl.limma2 = efit.dl2[["p.value"]]
head(p.dl.limma2)

saveRDS(v.dl, "results/limma_vdl_ARC.RDS")



### random sampling

R = 5000
set.seed(312)

#to store pvalues in
p.dl.rand = vector('list',length = R)

# to store "t" values (coefficients)
p.dl.rand.t = vector('list',length = R)

for( g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  
  group.dl.rand = sample(group.dl)
  
  # Model
  design.dl.rand = model.matrix(~0 + group.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  
  vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contr.matrix)
  
  efit.dl.rand2 = eBayes(vfit.dl.rand2)
  
  p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
  p.dl.rand.t[[g]] = efit.dl.rand2[["t"]]
  
}


q.dl <- Reduce(`+`, lapply(p.dl.rand, \(x) {
  (x < p.dl.limma2)
}))

head(q.dl)

q.dl = q.dl / R

q.dl = as.data.frame(q.dl)


efit.dl2[["p.value"]] <- q.dl

sum(duplicated(row.names(efit.dl2$coefficients)))



saveRDS(q.dl,("results/limma_vdl_cutoff5_2000_tworand_ARC.RDS"))

####

q.dl <- readRDS("results/limma_vdl_cutoff5_2000_tworand_ARC.RDS")


##### Analysis pulling genes out for each contrast 
tmp1 <- contrasts.fit(efit.dl2, coef = 1) # 
# tmp2 <- contrasts.fit(efit.dl2, coef = 2) # no 2 since design had no intercept?




ARC_limma_results <- topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('symbol') %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)



ARC_limma_results1 <- ARC_limma_results %>% 
  left_join(grcm38, by = "symbol") %>% 
  filter(!is.na(symbol)) %>% 
  filter(!is.na(entrez)) %>%
  select(symbol,logFC,P.Value,adj.P.Val,chr,entrez,start,end,biotype,description)




saveRDS(ARC_limma_results1,"results/ARC_limma_results.RDS")


ARC_limma_results1 <- readRDS("results/ARC_limma_results.RDS")

ARC_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC>0.2),
            Down = sum(logFC<0.2)) %>% 
  mutate(.,Total = Up + Down) 

hist(ARC_limma_results1$logFC)












