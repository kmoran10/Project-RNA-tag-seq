

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
library(organism, character.only = TRUE)
library(DOSE)
library(EnhancedVolcano)
library(tidyverse)
grcm38 # mouse genes


my_logFC_threshold = 0.2

LH <- readRDS("results/LH_limma_results.RDS")

y1a <- readRDS("results/LH_limma_results.RDS")

dc <- y1a %>% mutate(contrast = "Control vs. Stressed") %>% mutate(log10 = -log10(P.Value))


dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.5)%>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.5) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 



vp_LH <- ggplot(data = dc, 
                  aes(x = logFC, 
                      y = log10, 
                      colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange2", "grey","blue3"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-2.5,2.5),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-.75,.75),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black",vjust =1, hjust =.45)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1, vjust = -.55)+
  labs(title = "A. Differential Gene Expression in LH",
       x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))+
  theme_bw() +
  annotate(geom="text", x=2.5, y=.5, label=paste0(" ", "\n", "Lower in Stressed"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.5, y=.5, label=paste0(" ", "\n", "Higher in Stressed"),
           color="black", size = 5)+
  scale_x_continuous(limits = c(-3.5,3.5),breaks = c(-3,-2,-1,0,1,2,3))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )

vp_LH #875 x 550


### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 254 
dc %>% filter(between(logFC, .2, .75))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 2 05
dc %>% filter(between(logFC, .75, 1.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 45 
dc %>% filter(between(logFC, 1.5, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #4 

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 271 
dc %>% filter(between(logFC, -0.75, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 203
dc %>% filter(between(logFC, -1.5, -0.75)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 56
dc %>% filter(between(logFC, -3, -1.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 12 

#top 10 genes down in stressed
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes up in stressed 
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
################



 y1a %>% filter(symbol == "Nmu")
 y1a %>% filter(entrez == "11604")
 
 
 
 #875 x 550
ggplot(data = dc, 
        aes(x = logFC, 
            y = log10, 
            colour=diffexpressed)) +
   geom_point(alpha=0.25, size=3.5) +
   scale_color_manual(values=c("orange2", "grey","blue3"))+
   xlim(c(-2.5, 2.5)) +
   ylim(0,4)+
   # geom_vline(xintercept=c(-2.5,2.5),lty=4,col="black",lwd=0.8) +
   geom_vline(xintercept=c(-.75,.75),lty=4,col="black",lwd=0.8) +
   geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
   geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
   geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black",vjust =1, hjust =.8)+
   geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1, vjust = -.8)+
   labs(title = "A. Differential Gene Expression in LH",
        x="log2 Fold Change",
        y=bquote(~-Log[10]~italic(eFDR)))+
   theme_bw() +
   annotate(geom="text", x=2.5, y=.5, label=paste0(" ", "\n", "Lower in Stressed"),
            color="black", size = 5)+
   annotate(geom="text", x=-2.5, y=.5, label=paste0(" ", "\n", "Higher in Stressed"),
            color="black", size = 5)+
   scale_x_continuous(limits = c(-3.5,3.5),breaks = c(-3,-2,-1,0,1,2,3))+
   theme(axis.text.x = element_text(vjust = 1,size = 20),
         # axis.ticks = element_blank(),
         axis.text.y = element_text(hjust = 0.5,size = 20),
         axis.text = element_text(color="#3C3C3C",size = 20),
         axis.title = element_text(size = 20),   
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = "none",
         strip.background = element_blank(),
         strip.text.x = element_text(size = 20),
         text = element_text(size = 20))








LH <- readRDS("results/LH_limma_results.RDS")

top25sigdiff_genes_LH <- LH %>% 
  filter(logFC >= 1.25 | logFC <= -1.25) %>%
  filter(P.Value < 0.05) %>% 
  arrange(desc(abs(logFC))) %>% 
  select(symbol, logFC, P.Value, chr, description) %>% 
  head(25)

write.csv(top25sigdiff_genes_LH,"results/results_tables/top25sigdiff_genes_LH.csv", row.names = F)

