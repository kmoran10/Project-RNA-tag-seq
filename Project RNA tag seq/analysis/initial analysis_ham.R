
### initial cleaning and checking raw counts of genes per subject and ROI. HAM GENOME!!

library(limma)
library(tidyverse)
library(scales)


library(tidyverse)
library(DESeq2)
library(pheatmap)
library(annotables)



rawcounts <- read.delim("raw_data/combined_using_Maur2_ref.gff.txt",header=F)

colnames(rawcounts) <- c("gene","counts", "id")

bc <- rawcounts


bcx <- bc %>%  pivot_wider(names_from = id, values_from = counts)
bcx1 <- bcx %>% column_to_rownames(.,var = "gene")

write.csv(bcx1, "raw_data/rnaseq_counts_ham.csv", row.names = F)

id <- read.csv("raw_data/id.csv")



bc1 <- bc %>%
  mutate(subject = str_replace(id, "([A-Za-z]+)(\\d+)", "KM\\2")) %>% 
  mutate(ROI = str_extract(id, "[A-Za-z]+"))

df1 <- bc1 %>%  pivot_wider(names_from = gene, values_from = counts)

a <- merge(id, df1)
a1 <- a %>% relocate(id, .before=subject) %>% relocate(ROI, .after=subject) %>% select(c(1:23))
idxx <- a1 %>% column_to_rownames(.,var = "id")


bcx1 <- bcx1[,rownames(idxx)]
all(rownames(idxx) == colnames(bcx1)) #check


#### count checks


colSums(bcx1[,]) %>% 
  as_tibble_row() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'id') %>% 
  dplyr::rename(genecounts = V1) -> genecounts

genecounts <- genecounts %>%
  full_join(a1) %>% na.omit() 

genecounts %>% 
  ggplot(aes(genecounts)) +
  geom_histogram(bins = 54,alpha =0.5,color = 'grey') +
  theme_classic() 

counts <- genecounts %>%
  ggplot(aes(id,genecounts))+
  geom_bar(stat = 'identity')+
  coord_flip()

ggsave("imgs/braincounts_byid_ham.png",counts, height = 5,width= 5, dpi = 150)

#

source("functions/geom_boxjitter.R")

counts2 <- genecounts %>%
  ggplot(aes(group,genecounts, color = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.6)) +
  theme_classic()+
  theme(legend.position = "none")


ggsave("imgs/braincounts_bygroup_ham.png",counts2, height = 5,width= 5, dpi = 150)


counts3 <- genecounts %>%
  ggplot(aes(group,genecounts, color = ROI))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  theme_classic()


ggsave("imgs/braincounts_bygroupandROI_ham.png",counts3, height = 5,width= 5, dpi = 150)



ROIsplit <- split(idxx, idxx$ROI)

ARC_id <- ROIsplit$ARC
DMH_id <- ROIsplit$DMH
LH_id <- ROIsplit$LH


write.csv(ARC_id, "ham_brain_data/ARC_id.csv", row.names = T)
write.csv(DMH_id, "ham_brain_data/DMH_id.csv", row.names = T)
write.csv(LH_id, "ham_brain_data/LH_id.csv", row.names = T)

write.csv(idxx, "ham_brain_data/all_id.csv", row.names = T)



ARC_counts <- bcx1 %>%
  select(contains("ARC"))

DMH_counts <- bcx1 %>%
  select(contains("DMH"))

LH_counts <- bcx1 %>%
  select(contains("LH"))


write.csv(ARC_counts, "ham_brain_data/ARC_counts.csv", row.names = T)
write.csv(DMH_counts, "ham_brain_data/DMH_counts.csv", row.names = T)
write.csv(LH_counts, "ham_brain_data/LH_counts.csv", row.names = T)

write.csv(bcx1, "ham_brain_data/all_counts.csv", row.names = T)

