######################################################################################################
# Comparing Ki-67 RNAseqV2 gene expression levels in TCGA
# Script author: David Chen
# Script maintainer: David Chen
# Date: 03/14/18
# Notes:
######################################################################################################

rm(list=ls())
source("~/repos/Repos_for_Manuscript_Code/BRCA1-like_analyses/helper_functions.R");

## Tumor annotation:
my_samples <- loadReceptorPositiveTumors(receptorPosOnly=TRUE); 
my_samples$group[my_samples$SVM_BRCA1=="BRCA1-like"] <- 1;
my_samples$group[my_samples$SVM_BRCA1=="non-BRCA1-like"] <- 0;
my_samples$group <- as.factor(my_samples$group);

## Load & subset RNA-seq data:
tcgaGeneExpr <- loadTCGAGeneLeveExpr(); 
tcgaGeneExpr <- tcgaGeneExpr[ , colnames(tcgaGeneExpr) %in% my_samples$patients];

## Subset to gene of interest:
ki67 <- subset(tcgaGeneExpr, rownames(tcgaGeneExpr) %in% "MKI67");
ki67 <- as.data.frame(t(ki67)); 
ki67$patients <- rownames(ki67);
ki67 <- merge(ki67, my_samples[ , c("patients","group","Age","Stage","ER","PR","HER2")], by="patients");

## Statistical test:
fitki67 <- lm(
  MKI67 ~ group + Age + Stage + ER + PR + HER2, 
  data = ki67
);
summary(fitki67)

## Data visualization:
table(ki67$group)
ki67$BRCAness[ki67$group==1] <- "BRCA1-like (n=159)";
ki67$BRCAness[ki67$group==0] <- "non-BRCA1-like (n=578)";

png("~/Downloads/BRCA1ness_figures/051618_Figure3C.png", res=300, units="in", height=8.27, width=6);
ggplot(ki67, aes(x=BRCAness, y=MKI67)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0) +
  geom_point(position=position_jitter(width=0.25), alpha=0.3) + 
  theme_classic() +
  theme(axis.text.x=element_text(size=21,color="black"), axis.text.y=element_text(size=21,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=30,color="black") ) +    
  labs(y="MKI67 expression") +
  
  geom_segment(aes(x=1, y=10000, xend=2, yend=10000), size=0.3, inherit.aes=FALSE) +
  geom_segment(aes(x=1, y=9000, xend=1, yend=10000), size=0.3, inherit.aes=FALSE) +
  geom_segment(aes(x=2, y=8000, xend=2, yend=10000), size=0.3, inherit.aes=FALSE) +
  annotate("text", x=1.5, y=10500, label="*** P = 1.21e-07", size=7)
dev.off();
