#######################################################################################################
# Exploratory analysis of patient survival by Kaplan-Meier Method in TCGA
# Script author: David Chen
# Script maintainer: David Chen
# Notes:
#######################################################################################################

rm(list=ls())
library(gridExtra)
source("~/repos/Repos_for_Manuscript_Code/BRCA1-like_analyses/helper_functions.R"); 
source("~/repos/Repos_for_Manuscript_Code/BRCA1-like_analyses/plot_themes.R"); 

## Load study population & merge w/ survival data:
my_samples <- loadReceptorPositiveTumors(receptorPosOnly=TRUE);
metaSurv <- loadTCGASurvivalMeta();
sum(metaSurv$patients %in% my_samples$patients)

clin.breast <- merge(
  metaSurv[ , c("patients", "OS.time","OS", "PFI.time","PFI", "DFI.time","DFI")], 
  my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER","PR","HER2","TNBC","PAM50","PAM50lite")], 
  by = "patients"
);

## For survival analysis, also remove HER2+ tumors:
table(clin.breast$HER2, useNA="always")
clin.breast <- subset(clin.breast, HER2=="Negative");

## Exclude tumors without stage
my_samples <- subset(my_samples, !is.na(Stage));

## Survival modeling:
clin.breast$group <- 1 * (clin.breast$SVM_BRCA1=="BRCA1-like");
sModels <- computeKaplanMeier(clin.breast);
sModels$OS
sModels$PFS
sModels$DFS

## Data visualization:
myPlots <- drawKaplanMeier(sModels, pTheme=mySurvTheme);
grid.arrange(
  myPlots$OS$plot, 
  myPlots$PFS$plot,
  myPlots$DFS$plot,
  ncol = 3
)
