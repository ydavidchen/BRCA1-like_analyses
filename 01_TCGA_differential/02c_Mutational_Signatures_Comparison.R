####################################################################################################
# Correlation of WTSI mutational signature with BRCA1ness in TCGA
# Script author: David Chen
# Date: 03/14/2018
# Notes:
####################################################################################################

rm(list=ls())

library(ggplot2)
library(matrixStats)
library(reshape2)
library(doParallel); registerDoParallel(detectCores() - 1)

## Clinical annotation for study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);
my_samples <- subset(my_samples, TNBC=="Non-TNBC");

## Load mutational signatures:
rosenthalTCGA <- read.table("~/Dropbox (Christensen Lab)/Pan-cancer-analyses/TCGA_pan_cancer_signature_Rosenthal.txt", header=T);
rosenthalTCGA$id <- gsub(".", "-", rosenthalTCGA$id, fixed=TRUE);
rosenthalTCGA <- subset(rosenthalTCGA, cancer == "BRCA");
colnames(rosenthalTCGA)[1] <- "patients"; #for merging

## Merge:
comp.mat <- merge(my_samples, rosenthalTCGA, by="patients");
comp.mat <- subset(comp.mat, method == "WTSI");

## For statistical tests:
comp.mat$group[comp.mat$SVM_BRCA1=="BRCA1-like"] <- 1;
comp.mat$group[comp.mat$SVM_BRCA1=="non-BRCA1-like"] <- 0;
comp.mat$group <- as.factor(comp.mat$group);

## For data visualization:
comp.mat$BRCAness[comp.mat$SVM_BRCA1=="BRCA1-like"] <- "BRCA1-like (n=151)";
comp.mat$BRCAness[comp.mat$SVM_BRCA1=="non-BRCA1-like"] <- "non-BRCA1-like (n=505)";

#----------------------------------------Statistical test: Signature 3----------------------------------------
fitSig3 <- lm(
  Signature.3 ~ group + Age + Stage + ER + PR + HER2, 
  data = comp.mat
);
summary(fitSig3)

png("~/Downloads/BRCA1ness_figures/Figure2B.png", res=300, units="in", height=8.27, width=5.84);
ggplot(comp.mat, aes(x=BRCAness, y=Signature.3, color=BRCAness)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0, position=position_dodge(width=0.775)) + #must hide outliers if jitters are to be shown
  geom_point(aes(color=BRCAness), position=position_jitterdodge(jitter.width=0.25), alpha=0.25) + # add jitter
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  labs(y="Somatic Substitution Signature 3") +
  theme_classic() +
  theme(axis.text.x=element_text(size=16,color="black"), axis.text.y=element_text(size=16,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=16,color="black"),
        legend.position="top", legend.title=element_blank(),legend.text=element_text(size=14,color="black",face="bold"),
        strip.text.x=element_text(size=12,colour="black",face="bold")) +
  
  geom_segment(aes(x=1, y=1, xend=2, yend=1), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1, y=0.9, xend=1, yend=1), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2, y=0.8, xend=2, yend=1), size=0.3, inherit.aes=F) +
  annotate("text", x=1.5, y=1.03, label="P = 2.46e-05 ***", size=7)
dev.off();

#----------------------------------------Statistical test: Signature 1----------------------------------------
fitSig1 <- lm(
  Signature.1 ~ group + Age + Stage + ER + PR + HER2, 
  data = comp.mat
);
summary(fitSig1)

png("~/Downloads/BRCA1ness_figures/FigureS6.png", res=300, units="in", height=8.27, width=5.84);
ggplot(comp.mat, aes(x=BRCAness, y=Signature.1, color=BRCAness)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0, position=position_dodge(width=0.775)) + #must hide outliers if jitters are to be shown
  geom_point(aes(color=BRCAness), position=position_jitterdodge(jitter.width=0.25), alpha=0.25) + # add jitter
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  labs(y="Somatic Mutational Signature 1") +
  theme_classic() +
  theme(axis.text.x=element_text(size=16,color="black"), axis.text.y=element_text(size=16,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=16,color="black"),
        legend.position="top", legend.title=element_blank(),legend.text=element_text(size=14,color="black",face="bold"),
        strip.text.x=element_text(size=12,colour="black",face="bold")) +
  
  geom_segment(aes(x=1, y=1.02, xend=2, yend=1.02), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1, y=0.95, xend=1, yend=1.02), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2, y=1.01, xend=2, yend=1.02), size=0.3, inherit.aes=F) +
  annotate("text", x=1.5, y=1.04, label="P = 0.0056 **", size=7)
dev.off();

