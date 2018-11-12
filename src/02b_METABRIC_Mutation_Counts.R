# Compare METABRIC Mutational Burden in Receptor-positive Tumors
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

cBIO_MUT_COUNT_PATH <- "~/Dropbox (Christensen Lab)/Breast_Cancer_Data_Sets/METABRIC_data_set_cBio/cBio_METABRIC_mutation_count_query.csv";

sample_clinical <- loadMETABRICtumors(receptorPosOnly=TRUE); 
compMat <- read.csv(cBIO_MUT_COUNT_PATH, stringsAsFactors=FALSE);
compMat <- merge(compMat, sample_clinical, by="SAMPLE_ID");

## Statistical test: Poisson regression
compMat$group[compMat$SVM_BRCA1=="BRCA1-like"] <- 1;
compMat$group[compMat$SVM_BRCA1=="non-BRCA1-like"] <- 0;
compMat$group <- as.factor(compMat$group);
table(compMat$group)

fit <- rlm(
  Mutation.Count ~ group + AGE_AT_DIAGNOSIS + Stage + ER_STATUS + PR_STATUS + HER2_STATUS,
  data = compMat
);
summary(fit)
svar.rlm <- sandwich(fit);
coeftest(fit, vcov.=svar.rlm, type="HC0");

## Calculate mean by group:
mean(compMat$Mutation.Count[compMat$group==1], na.rm=TRUE)
sd(compMat$Mutation.Count[compMat$group==1], na.rm=TRUE)

mean(compMat$Mutation.Count[compMat$group==0], na.rm=TRUE)
sd(compMat$Mutation.Count[compMat$group==0], na.rm=TRUE)

## Data visualization:
compMat$SVM_BRCA1[compMat$SVM_BRCA1=="BRCA1-like"] <- "BRCA1-like \n (n=104)";
compMat$SVM_BRCA1[compMat$SVM_BRCA1=="non-BRCA1-like"] <- "non-BRCA1-like \n (n=1085)";
png("~/Downloads/BRCA1ness_figures/METABRIC_mut_count.png", res=300, units="in", height=8.27, width=6);
ggplot(compMat, aes(x=SVM_BRCA1, y=Mutation.Count)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0) +
  geom_point(position=position_jitter(width=0.25), alpha=0.3) + 
  labs(y="Mutation count", title="METABRIC") +
  myBoxplotTheme +
  geom_segment(aes(x=1, y=25, xend=2, yend=25), size=0.3, inherit.aes=FALSE) +
  annotate("text", x=1.5, y=25.5, label="P = 0.30", size=7)
dev.off();
