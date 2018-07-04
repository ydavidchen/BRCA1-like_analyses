######################################################################################################
# Comparing DNMT1/3A/3B RNAseq gene expression levels in TCGA
# Script author: David Chen
# Script maintainer: David Chen
# Notes:
######################################################################################################

rm(list=ls())
library(ggplot2)
library(reshape2)
source("~/repos/Repos_for_Manuscript_Code/BRCA1-like_analyses/helper_functions.R");
my_samples <- loadReceptorPositiveTumors(receptorPosOnly=TRUE); 
tcgaGeneExpr <- loadTCGAGeneLeveExpr(); 
tcgaGeneExpr <- tcgaGeneExpr[ , colnames(tcgaGeneExpr) %in% my_samples$patients];

## Subset to gene of interest for statistical test & data visualization:
plt.DNMTs <- subset(tcgaGeneExpr, rownames(tcgaGeneExpr) %in% c("DNMT1", "DNMT3A", "DNMT3B"));
plt.DNMTs <- as.data.frame(t(plt.DNMTs));
plt.DNMTs$patients <- rownames(plt.DNMTs);
plt.DNMTs <- merge(plt.DNMTs, my_samples[ , c("patients","TNBC","SVM_BRCA1","Age","Stage","ER", "PR", "HER2")],by="patients");
plt.DNMTs$BRCAness[plt.DNMTs$SVM_BRCA1=="BRCA1-like"] <- 1;
plt.DNMTs$BRCAness[plt.DNMTs$SVM_BRCA1=="non-BRCA1-like"] <- 0;
plt.DNMTs$BRCAness <- as.factor(plt.DNMTs$BRCAness);
table(plt.DNMTs$BRCAness)

## Statistical tests:
fitDNMT1 <- lm(
  DNMT1 ~ BRCAness + Age + Stage + ER + PR + HER2, 
  data = plt.DNMTs
);
summary(fitDNMT1)

fitDNMT3A <- lm(
  DNMT3A ~ BRCAness + Age + Stage + ER + PR + HER2, 
  data = plt.DNMTs
);
summary(fitDNMT3A)

fitDNMT3B <- lm(
  DNMT3B ~ BRCAness + Age + Stage + ER + PR + HER2, 
  data = plt.DNMTs
);
summary(fitDNMT3B)

## Data visualization:
plt.DNMTs$TNBC <- paste(plt.DNMTs$TNBC, "(159 BRCA1-like, 578 non-BRCA1-like)");
plt.DNMTs <- melt(plt.DNMTs[ , c("patients","SVM_BRCA1","TNBC","DNMT1","DNMT3A","DNMT3B")]);
ggplot(plt.DNMTs, aes(x=variable, y=value, color=SVM_BRCA1)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0, position=position_dodge(width=0.775)) + #must hide outliers if jitters are to be shown
  geom_point(aes(color=SVM_BRCA1), position=position_jitterdodge(jitter.width=0.25), alpha=0.25) + # add jitter
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  scale_y_continuous(limits = c(0,8000))+
  labs(y="mRNA transcript count") +
  theme_classic() +
  theme(axis.text.x=element_text(face="italic",size=16,color="black"), axis.text.y=element_text(size=16,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=16,color="black"),
        legend.position="top", legend.title=element_blank(),legend.text=element_text(size=14,color="black",face="bold"),
        strip.text.x=element_text(size=12,colour="black",face="bold")) +  
  
  facet_wrap( ~ Subtype, ncol=1) +
  geom_segment(aes(x=0.81, y=7000, xend=1.19, yend=7000), size=0.3, inherit.aes=F) + #DNMT1
  geom_segment(aes(x=1.81, y=5000, xend=2.19, yend=5000), size=0.3, inherit.aes=F) + #DNMT3A
  geom_segment(aes(x=2.81, y=3000, xend=3.19, yend=3000), size=0.3, inherit.aes=F) + #DNMT3B
  geom_segment(data=data.frame(x=0.81,y=6600,xend=0.81,yend=7000,Subtype="Non-TNBC (159 BRCA1-like, 578 non-BRCA1-like)"), aes(x=x,y=y,xend=xend,yend=yend),size=0.3,inherit.aes=F) +
  geom_segment(data=data.frame(x=1.19,y=4200,xend=1.19,yend=7000,Subtype="Non-TNBC (159 BRCA1-like, 578 non-BRCA1-like)"), aes(x=x,y=y,xend=xend,yend=yend),size=0.3,inherit.aes=F) +
  geom_segment(data=data.frame(x=1.81,y=3000,xend=1.81,yend=5000,Subtype="Non-TNBC (159 BRCA1-like, 578 non-BRCA1-like)"), aes(x=x,y=y,xend=xend,yend=yend),size=0.3,inherit.aes=F) +
  geom_segment(data=data.frame(x=2.19,y=2000,xend=2.19,yend=5000,Subtype="Non-TNBC (159 BRCA1-like, 578 non-BRCA1-like)"), aes(x=x,y=y,xend=xend,yend=yend),size=0.3,inherit.aes=F) +
  geom_segment(data=data.frame(x=2.81,y=2100,xend=2.81,yend=3000,Subtype="Non-TNBC (159 BRCA1-like, 578 non-BRCA1-like)"), aes(x=x,y=y,xend=xend,yend=yend),size=0.3,inherit.aes=F) +
  geom_segment(data=data.frame(x=3.19,y=1200,xend=3.19,yend=3000,Subtype="Non-TNBC (159 BRCA1-like, 578 non-BRCA1-like)"), aes(x=x,y=y,xend=xend,yend=yend),size=0.3,inherit.aes=F) +  geom_text(data=data.frame(x=1,y=7500,label="P = 4.06e-10 ***",Subtype="Non-TNBC (159 BRCA1-like, 578 non-BRCA1-like)"), aes(x=x,y=y,label=label),size=5,inherit.aes=F) +
  geom_text(data=data.frame(x=2,y=5500,label="P = 2.16e-05 ***",Subtype="Non-TNBC (159 BRCA1-like, 578 non-BRCA1-like)"), aes(x=x,y=y,label=label),size=5,inherit.aes=F) +
  geom_text(data=data.frame(x=3,y=3500,label="P = 3.84e-09 ***",Subtype="Non-TNBC (159 BRCA1-like, 578 non-BRCA1-like)"), aes(x=x,y=y,label=label),size=5,inherit.aes=F)

