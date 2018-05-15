##########################################################################################################
# Relation of Total Mutational Rates with SVM BRCA1-like Status
# Script author: David Chen
# Date: 02/23/18; 05/08/18 (revision)
# Notes:
##########################################################################################################

rm(list=ls())

library(gdata)
library(ggplot2)
library(lmtest)
library(MASS)
library(matrixStats)
library(reshape2)
library(sandwich)
library(doParallel); registerDoParallel(detectCores() - 1)

## Clinical annotation for the study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);

## Load mutational count & rates data (Kandoth et al. 2013):
mutMeasure <- read.xls("~/Dropbox (Christensen Lab)/Pan-cancer-analyses/Kandoth2013Nature/Supplementary_Table_3a.xls", stringsAsFactors=F);
mutMeasure <- subset(mutMeasure, Cancer.Type == "brca");
mutMeasure$patients <- substr(mutMeasure$TCGA.ID, 1, 12);
stopifnot( anyDuplicated(mutMeasure$patients) == 0 )

sum(my_samples$patients %in% mutMeasure$patients)

## Join with sample annotation:
comp.mat <- merge(
  my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER", "PR", "HER2","TNBC")], 
  mutMeasure, 
  by = "patients"
);
table(comp.mat$SVM_BRCA1[comp.mat$TNBC=="TNBC"])
table(comp.mat$SVM_BRCA1[comp.mat$TNBC=="Non-TNBC"])
table(comp.mat$SVM_BRCA1)

## Statistical test:
comp.mat$group[comp.mat$SVM_BRCA1=="BRCA1-like"] <- 1;
comp.mat$group[comp.mat$SVM_BRCA1=="non-BRCA1-like"] <- 0;
comp.mat$group <- as.factor(comp.mat$group);

fit.rates <- rlm(
  `Mutation.Rate...Mbp.` ~ group + Age + Stage + ER + PR + HER2,
  data = subset(comp.mat, TNBC=="Non-TNBC")
);
summary(fit.rates)
svar.rlm <- sandwich(fit.rates);
coeftest(fit.rates, vcov.=svar.rlm, type="HC0");
## z test of coefficients:
##   
##   Estimate Std. Error z value  Pr(>|z|)    
##  (Intercept)        0.5566147  0.2321946  2.3972  0.016521 *  
##  group1             0.3407573  0.1070674  3.1826  0.001459 ** 
##  Age                0.0112294  0.0024814  4.5255 6.025e-06 ***
##  StageStage III-IV  0.1874977  0.0900624  2.0819  0.037355 *  
##  ERPositive        -0.1422425  0.2219066 -0.6410  0.521522    
##  PRPositive        -0.1794163  0.1265139 -1.4182  0.156146    
##  HER2Positive       0.3206599  0.1020799  3.1413  0.001682 ** 
##  ---
##  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Data visualization:
comp.mat$BRCAness[comp.mat$SVM_BRCA1=="BRCA1-like"] <- "BRCA1-like (n=205)";
comp.mat$BRCAness[comp.mat$SVM_BRCA1=="non-BRCA1-like"] <- "non-BRCA1-like (n=457)";
png("~/Downloads/BRCA1ness_figures/Figure2B.png", res=300, units="in", height=8.27, width=5.84);
ggplot(comp.mat, aes(x=BRCAness, y=Mutation.Rate...Mbp., color=BRCAness)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0) +
  geom_point(aes(color=BRCAness), position=position_jitterdodge(dodge.width=0.5), alpha=0.3) + 
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=20,color="black"), axis.text.y=element_text(size=20,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=22,color="black"),
        legend.position="top", legend.title=element_blank(),legend.text=element_text(size=18,color="black",face="bold"),
        strip.text.x=element_text(size=12,colour="black",face="bold")) +    
  labs(y="Mutation rate per mega base-pair (Mb)") +
  geom_segment(aes(x=1, y=11.5, xend=2, yend=11.5), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1, y=10, xend=1, yend=11.5), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2, y=8, xend=2, yend=11.5), size=0.3, inherit.aes=F) +
  annotate("text", x=1.5, y=12, label="P = 0.017**", size=7)
dev.off();
