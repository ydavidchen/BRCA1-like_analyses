# Comparison of Somatic Mutational Signatures by SVM BRCA1-like status in TCGA
# Script author: David Chen
# Script maintainer: David Chen
# Notes:

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
source("plot_themes.R");

my_samples <- loadReceptorPositiveTumors(receptorPosOnly=TRUE);
rosenthalTCGA <- loadMutationalSignatures();
comp.mat <- merge(my_samples, rosenthalTCGA, by="patients");
comp.mat <- subset(comp.mat, method == "WTSI");
comp.mat$group[comp.mat$SVM_BRCA1=="BRCA1-like"] <- 1;
comp.mat$group[comp.mat$SVM_BRCA1=="non-BRCA1-like"] <- 0;
comp.mat$group <- as.factor(comp.mat$group);
table(comp.mat$group)
comp.mat$BRCAness[comp.mat$SVM_BRCA1=="BRCA1-like"] <- "BRCA1-like \n (n=151)";
comp.mat$BRCAness[comp.mat$SVM_BRCA1=="non-BRCA1-like"] <- "non-BRCA1-like \n (n=505)";

#----------------------------------------Statistical test: Signature 3----------------------------------------
fitSig3 <- lm(
  Signature.3 ~ group + Age + Stage + ER + PR + HER2, 
  data = comp.mat
);
summary(fitSig3)

png("~/Downloads/Final_revision/Figure3B.png", res=300, units="in", height=8.27, width=6.5);
ggplot(comp.mat, aes(x=BRCAness, y=Signature.3)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0) +
  geom_point(position=position_jitter(width=0.25), alpha=0.3) + 
  labs(y="Mutational Signature 3") +
  myBoxplotTheme +
  
  geom_segment(aes(x=1, y=1, xend=2, yend=1), size=0.3, inherit.aes=FALSE) +
  geom_segment(aes(x=1, y=0.9, xend=1, yend=1), size=0.3, inherit.aes=FALSE) +
  geom_segment(aes(x=2, y=0.8, xend=2, yend=1), size=0.3, inherit.aes=FALSE) +
  annotate("text", x=1.5, y=1.03, label="*** P = 2.46e-05", size=10)
dev.off();

#----------------------------------------Statistical test: Signature 1----------------------------------------
fitSig1 <- lm(
  Signature.1 ~ group + Age + Stage + ER + PR + HER2, 
  data = comp.mat
);
summary(fitSig1)

myBoxplotTheme$legend.position <- "none";
ggplot(comp.mat, aes(x=BRCAness, y=Signature.1, color=BRCAness)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0, position=position_dodge(width=0.775)) +
  geom_point(aes(color=BRCAness), position=position_jitterdodge(jitter.width=0.25), alpha=0.25) +
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  labs(y="Mutational Signature 1") +
  myBoxplotTheme +
  geom_segment(aes(x=1, y=1.02, xend=2, yend=1.02), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1, y=0.95, xend=1, yend=1.02), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2, y=1.01, xend=2, yend=1.02), size=0.3, inherit.aes=F) +
  annotate("text", x=1.5, y=1.04, label="** P = 0.0056", size=7)

