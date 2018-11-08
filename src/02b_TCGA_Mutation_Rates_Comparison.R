# Comparison of Mutational Rates per Mb with SVM BRCA1-like Status
# Script author: David Chen
# Script maintainer: David Chen
# Notes:

rm(list=ls()); 
library(lmtest);
library(MASS);
library(matrixStats);
library(sandwich);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
source("plot_themes.R");

## Load clinical & mutational data:
my_samples <- loadReceptorPositiveTumors(receptorPosOnly=TRUE);
mutMeasure <- loadMutationMeasures(); 
comp.mat <- merge(
  my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER", "PR", "HER2","TNBC")], 
  mutMeasure, 
  by = "patients"
);
table(comp.mat$SVM_BRCA1)

## Statistical modeling:
comp.mat$group[comp.mat$SVM_BRCA1 == "BRCA1-like"] <- 1;
comp.mat$group[comp.mat$SVM_BRCA1 == "non-BRCA1-like"] <- 0;
comp.mat$group <- as.factor(comp.mat$group);
fit.rates <- rlm(
  `Mutation.Rate...Mbp.` ~ group + Age + Stage + ER + PR + HER2,
  data = subset(comp.mat, TNBC=="Non-TNBC")
);
summary(fit.rates)
svar.rlm <- sandwich(fit.rates);
coeftest(fit.rates, vcov.=svar.rlm, type="HC0");

## Data visualization:
comp.mat$BRCAness[comp.mat$SVM_BRCA1=="BRCA1-like"] <- "BRCA1-like \n (n=138)";
comp.mat$BRCAness[comp.mat$SVM_BRCA1=="non-BRCA1-like"] <- "non-BRCA1-like \n (n=453)";

png("~/Downloads/BRCA1ness_figures/Figure3A.png", res=300, units="in", height=8.27, width=6);
ggplot(comp.mat, aes(x=BRCAness, y=Mutation.Rate...Mbp.)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0) +
  geom_point(position=position_jitter(width=0.25), alpha=0.3) + 
  labs(y="Mutation rate per Mb") +
  myBoxplotTheme +
  
  geom_segment(aes(x=1, y=11.5, xend=2, yend=11.5), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1, y=10, xend=1, yend=11.5), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2, y=8, xend=2, yend=11.5), size=0.3, inherit.aes=F) +
  annotate("text", x=1.5, y=11.95, label="** P = 0.0015", size=7)
dev.off();
