# Verifying expression of miRNA showed differential methylation in BRCA1-like tumors
# Script author: David Chen
# Script maintainer: David Chen
# Notes:
# 1. miR124-2 / miR-124a2: HYPERmethylated in BRCA1-like; modulates proliferation (PMID:25731732)

rm(list=ls());
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
source("plot_themes.R");

## Data loading & subsetting:
my_samples <- loadReceptorPositiveTumors(receptorPosOnly=TRUE);
miRExpr <- loadTCGAmiRNA();

my_miRs <- subset(miRExpr, rownames(miRExpr) %in% c("hsa-mir-124-2"));
my_miRs <- as.data.frame(t(my_miRs));
my_miRs$patients <- rownames(my_miRs);
my_miRs <- merge(
  my_miRs, 
  my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER","PR","HER2")], 
  by="patients"
);
table(my_miRs$SVM_BRCA1)
my_miRs$group[my_miRs$SVM_BRCA1 == "BRCA1-like"] <- 1;
my_miRs$group[my_miRs$SVM_BRCA1 == "non-BRCA1-like"] <- 0;
table(my_miRs$group)

## Statistical modeling:
fit_miR124a2 <- lm(
  `hsa-mir-124-2` ~ group + Age + Stage + ER + PR + HER2,
  data = my_miRs
);
summary(fit_miR124a2)

## Data visualization:
my_miRs$BRCAness[my_miRs$SVM_BRCA1 == "BRCA1-like"] <- "BRCA1-like \n (n=106)";
my_miRs$BRCAness[my_miRs$SVM_BRCA1 == "non-BRCA1-like"] <- "non-BRCA1-like \n (n=400)";
myBoxplotTheme$legend.position <- "none"; 
ggplot(my_miRs, aes(x=BRCAness, y=`hsa-mir-124-2`, color=BRCAness)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0, position=position_dodge(width=0.775)) +
  geom_point(aes(color=BRCAness), position=position_jitterdodge(jitter.width=0.25), alpha=0.5) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  labs(y="hsa-mir-124-2 miRseq Z-score") +
  geom_segment(aes(x=1, y=0, xend=2, yend=0), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1, y=-0.03, xend=1, yend=0), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2, y=-0.02, xend=2, yend=0), size=0.3, inherit.aes=F) +
  annotate("text", x=1.5, y=0.005, label="* P = 0.021", size=6) +
  myBoxplotTheme
