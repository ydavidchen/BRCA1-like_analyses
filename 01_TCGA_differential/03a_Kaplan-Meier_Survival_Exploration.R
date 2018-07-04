#######################################################################################################
# Explore Patient Survival
# Script author: David Chen
# Date: 07/04/18
# Notes:
# 1. Survival data (OS, PFS, DFS) were from Liu et al. 2018 Cell
#######################################################################################################

rm(list=ls())
library(gridExtra)
library(matrixStats)
source("~/repos/Repos_for_Manuscript_Code/BRCA1-like_analyses/helper_functions.R"); 
source("~/repos/Repos_for_Manuscript_Code/BRCA1-like_analyses/plot_themes.R"); 

## Load study population & merge w/ survival data:
my_samples <- loadReceptorPositiveTumors(receptorPosOnly=TRUE);

metaSurv <- loadSurvivalMeta();
sum(metaSurv$patients %in% my_samples$patients)

clin.breast <- merge(
  metaSurv[ , c("patients", "OS.time","OS", "PFI.time","PFI", "DFI.time","DFI")], 
  my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER","PR","HER2","TNBC","PAM50","PAM50lite")], 
  by = "patients"
);

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
