######################################################################################################
# Comparison of Horvath DNA methylation age by BRCA1-like status
# Script author: David Chen
# Date: 03/08/18; 04/24/18
# Notes:
# 1. Motivation: to follow-up with the differential methylation in BRCA1-like tumors.
######################################################################################################

rm(list=ls())
library(gdata)
library(ggplot2)
library(reshape2)
library(wateRmelon)

## Load data & results:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/betas.TCGAC450k.RData");
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/030418_DMRcate_Non-TNBCs.RData");

## Check if all Non-TNBCs:
table(my_samples$TNBC)

## Mutual subsetting
betas_TCGA450k <- betas_TCGA450k[ , colnames(betas_TCGA450k) %in% my_samples$patients];
dim(betas_TCGA450k)

## Execute methylation age calculation:
DNAmAge <- agep(betas_TCGA450k);
DNAmAge <- as.data.frame(DNAmAge);
colnames(DNAmAge) <- "Horvath";
DNAmAge$patients <- rownames(DNAmAge);

## Merge into sample annotation:
my_samples <- merge(my_samples, DNAmAge, by="patients"); 

#-------------------------------------------------Age acceleration & visualization-------------------------------------------------
## Calculate age acceleration (horvath vs. age regression residual; see Johnson et al)
hist(my_samples$Horvath); 

Clock_regress <- lm(Horvath ~ Age, data=my_samples);  #order of regression matters
summary(Clock_regress)$coefficients[2, 4] #Exact P-value
## 1.339065e-08
cor.test(my_samples$Age, my_samples$Horvath, method="spearman")
## S = 3786100, p-value = 4.444e-09
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##   rho 
## 0.3195721 
my_samples$`Methylation age acceleration (years)` <- resid(Clock_regress);

## Scatter plot with color:
my_samples$BRCAness[my_samples$SVM_BRCA1=="BRCA1-like"] <- "BRCA1-like (n=68)";
my_samples$BRCAness[my_samples$SVM_BRCA1=="non-BRCA1-like"] <- "non-BRCA1-like (n=254)";
table(my_samples$BRCAness) #check

## Show age acceleration:
ggplot(my_samples, aes(x=Age, y=Horvath, color=BRCAness)) +
  geom_point(aes(size=`Methylation age acceleration (years)`, alpha=BRCAness)) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_alpha_manual(values=c(0.65, 0.55)) +
  labs(x="Chronological Age (years)", y="Methylation Age (years)") +
  geom_abline(slope=1, intercept=0, color="darkgray") +
  
  scale_x_continuous(limit=c(10,130), breaks=seq(25,100,by=25)) +
  scale_y_continuous(limit=c(10,130), breaks=seq(25,125,by=25)) +
  scale_size_continuous(breaks=seq(-30, 30, by=15)) +
  
  theme_bw() +
  theme(axis.text=element_text(size=16,color="black"),
        axis.title=element_text(size=16,color="black"),
        legend.position="right", legend.title=element_text(size=14,color="black",face="bold"),legend.text=element_text(size=14,color="black",face="bold"),
        strip.text.x=element_text(size=12,colour="black",face="bold")) +
  annotate("text", 110, 25, label="Spearman's rho = 0.32", size=6) +
  annotate("text", 110, 15, label="P = 4.44E-09", size=6) 

## 2 scatter plots:
png("~/Downloads/BRCA1ness_figures/Figure3C.png", res=300, units="in", height=8.27, width=8.27);
ggplot(my_samples, aes(x=Age, y=Horvath, color=BRCAness)) +
  geom_point(aes(alpha=BRCAness), size=3) +
  geom_smooth(aes(fill=BRCAness), method="lm", alpha=0.3, se=FALSE) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_alpha_manual(values=c(0.9, 0.8)) +
  labs(x="Chronological Age (years)", y="Methylation Age (years)") +
  geom_abline(slope=1, intercept=0, color="lightgray", linetype="dashed") +
  scale_x_continuous(limit=c(25,130), breaks=seq(25,100,by=25)) +
  scale_y_continuous(limit=c(25,130), breaks=seq(25,125,by=25)) +
  scale_size_continuous(breaks=seq(-30, 30, by=15)) +
  theme_bw() +
  theme(axis.text=element_text(size=20,color="black"), axis.title=element_text(size=20,color="black"),
        legend.position="top", legend.title=element_blank(),legend.text=element_text(size=14,color="black",face="bold"),
        strip.text.x=element_text(size=12,colour="black",face="bold"))
dev.off();

#-------------------------------------------------Linear model controlling for covariates-------------------------------------------------
## Statistical tests:
my_samples$group[my_samples$SVM_BRCA1=="BRCA1-like"] <- 1;
my_samples$group[my_samples$SVM_BRCA1=="non-BRCA1-like"] <- 0;
my_samples$group <- as.factor(my_samples$group);
table(my_samples$group)
table(my_samples$SVM_BRCA1)

horvath.comp <- lm(Horvath ~ group + Stage + ER + PR + HER2, data=my_samples);
summary(horvath.comp)

chrono.comp <- lm(Age ~ group + Stage + ER + PR + HER2, data=my_samples);
summary(chrono.comp)

## Data visualization by group: 
plt.age <- melt(my_samples[ , c("patients", "BRCAness", "Age", "Horvath")]);
plt.age$variable <- gsub("Age", "Chronological Age", plt.age$variable);
plt.age$variable <- gsub("Horvath", "Methylation Age", plt.age$variable);

png("~/Downloads/BRCA1ness_figures/Figure3C_right.png", res=300, units="in", height=8.27, width=8.27);
ggplot(plt.age, aes(x=variable, y=value, color=BRCAness)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0, position=position_dodge(width=0.775)) + #must hide outliers if jitters are to be shown
  geom_point(position=position_jitterdodge(jitter.width=0.25), alpha=0.3) + # add jitter
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  scale_y_continuous(breaks=seq(25, 125, by=25), limits=c(20,154)) +
  theme_classic() +
  theme(axis.text=element_text(size=20,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=20,color="black"),
        legend.position="top", legend.title=element_blank(),legend.text=element_text(size=14,color="black",face="bold"),
        strip.text.x=element_text(size=12,colour="black",face="bold")) +
  
  labs(y="Age (years)") +
  
  geom_segment(aes(x=1.81, y=150, xend=2.19, yend=150), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1.81, y=120, xend=1.81, yend=150), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2.19, y=140, xend=2.19, yend=150), size=0.3, inherit.aes=F) +
  annotate("text", x=2, y=153, label="P = 0.0182", size=7) +
  
  geom_segment(aes(x=0.81, y=100, xend=1.19, yend=100), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=0.81, y=95, xend=0.81, yend=100), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1.19, y=95, xend=1.19, yend=100), size=0.3, inherit.aes=F) +
  annotate("text", x=1, y=103, label="n.s.", size=7)
dev.off();
