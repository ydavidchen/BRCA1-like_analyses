# Comparison of Ki-67 microarray gene expression by BRCA1-like status in METABRIC
# Script author: David Chen
# Script maintainer: David Chen
# Notes:

rm(list=ls()); 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
source("plot_themes.R");
sample_clinical <- loadMETABRICtumors(receptorPosOnly=TRUE); 
arrayExpr <- loadMETABRICarrayExpr(); 

## Load gene expression data & subset to Ki67:
plt.Ki67 <- arrayExpr[rownames(arrayExpr)=="MKI67", , drop=FALSE];
plt.Ki67 <- as.data.frame(t(plt.Ki67));
plt.Ki67$SAMPLE_ID <- rownames(plt.Ki67);
plt.Ki67 <- merge(
  plt.Ki67,
  sample_clinical[ , c("SAMPLE_ID","SVM_BRCA1","AGE_AT_DIAGNOSIS","Stage","ER_STATUS","PR_STATUS","HER2_STATUS")],
  by = "SAMPLE_ID"
);
plt.Ki67$group[plt.Ki67$SVM_BRCA1=="BRCA1-like"] <- 1;
plt.Ki67$group[plt.Ki67$SVM_BRCA1=="non-BRCA1-like"] <- 0;
plt.Ki67$group <- as.factor(plt.Ki67$group);
table(plt.Ki67$group)

## Statistical modeling:
fit_other <- lm(
  MKI67 ~ group + AGE_AT_DIAGNOSIS + Stage + ER_STATUS + PR_STATUS + HER2_STATUS,
  data = plt.Ki67
);
summary(fit_other)

## Data visualization:
plt.Ki67$BRCAness[plt.Ki67$SVM_BRCA1=="BRCA1-like"] <- "BRCA1-like \n (n=101)";
plt.Ki67$BRCAness[plt.Ki67$SVM_BRCA1=="non-BRCA1-like"] <- "non-BRCA1-like \n (n=1,046)";
myBoxplotTheme$legend.position <- "none";
myBoxplotTheme$plot.title$size <- 20;
ggplot(plt.Ki67, aes(x=BRCAness, y=MKI67, color=BRCAness)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0, position=position_dodge()) +
  geom_point(aes(color=BRCAness), position=position_jitterdodge(jitter.width=0.25), alpha=0.25) +
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  labs(y="MKI67 expression (microarray)", title="METABRIC") +
  myBoxplotTheme +
  geom_segment(aes(x=1, y=7.3, xend=2, yend=7.3), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1, y=7, xend=1, yend=7.3), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2, y=6.8, xend=2, yend=7.3), size=0.3, inherit.aes=F) +
  annotate("text", x=1.5, y=7.4, label="*** P = 2.49e-12", size=7)
