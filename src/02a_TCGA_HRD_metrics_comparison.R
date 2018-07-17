##########################################################################################################
# Relation of HR Deficiency Metrics with SVM BRCA1-like Score
# Script author: David Chen
# Script maintainer: David Chen
# Notes:
##########################################################################################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
source("plot_themes.R");

my_samples <- loadReceptorPositiveTumors(receptorPosOnly=FALSE);
MarquardHRD <- loadHRDmetrics();
comp.mat <- merge(
  my_samples[ , c("patients","BRCA1_prob","SVM_BRCA1","Age","Stage","ER", "PR", "HER2","TNBC")], 
  MarquardHRD, 
  by = "patients"
);
table(comp.mat$SVM_BRCA1[comp.mat$TNBC=="TNBC"])
table(comp.mat$SVM_BRCA1[comp.mat$TNBC=="Non-TNBC"])
comp.mat$Subtype[comp.mat$TNBC=="TNBC"] <- "TNBC (64 BRCA1-like, 7 non-BRCA1-like)";
comp.mat$Subtype[comp.mat$TNBC=="Non-TNBC"] <- "Non-TNBC (143 BRCA1-like, 503 non-BRCA1-like)";

#-------------------------------------------Correlative analyses-------------------------------------------
## HRD-LOH Score:
stopifnot( sum(is.na(comp.mat$HRD.LOH)) == 0 ) #ensure no missing value in outcome
fitHRD <- lm(HRD.LOH ~ BRCA1_prob, data=comp.mat);
summary(fitHRD)

## LST Score
stopifnot( sum(is.na(comp.mat$LST)) == 0 ) #ensure no missing value in outcome
fitLST <- lm(LST ~ BRCA1_prob, data = comp.mat);
summary(fitLST)

## Data visualization:
ggplot(comp.mat, aes(x=BRCA1_prob, y=HRD.LOH)) +
  geom_jitter() +
  geom_smooth(method="lm", se=FALSE) +
  labs(x="BRCA1-like probability", y="HRD-LOH Index") +
  myScatterTheme + 
  annotate(geom="text",x=0.25,y=40,label="P < 2.2e-16 ***",size=6) +
  annotate(geom="text",x=0.25,y=37,label="R-squared = 0.27",size=6) +
  annotate(geom="text",x=0.25,y=34,label="Coef. = 11.5",size=6)

ggplot(comp.mat, aes(x=BRCA1_prob, y=LST)) +
  geom_jitter() +
  geom_smooth(method="lm", se=FALSE) +
  labs(x="BRCA1-like probability", y="Large Scale Transition (LST) score") +
  theme_classic() +
  theme(axis.text.x=element_text(size=15,color="black"), axis.title.x=element_text(size=15,color="black"),
        axis.text.y=element_text(size=15,color="black"), axis.title.y=element_text(size=15,color="black"),
        strip.text.x=element_text(size=12,colour="black",face="bold"),
        legend.position="top", legend.title=element_blank(), legend.text=element_text(size=15,color="black")) +
  annotate(geom="text",x=0.25,y=40,label="P < 2.2e-16 ***",size=6) +
  annotate(geom="text",x=0.25,y=37,label="R-squared = 0.31",size=6) +
  annotate(geom="text",x=0.25,y=34,label="Coef. = 14.2",size=6)
