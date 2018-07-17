#######################################################################################################
# Exploratory analysis of patient survival by Kaplan-Meier Method in TCGA
# Script author: David Chen
# Script maintainer: David Chen
# Notes:
#######################################################################################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
source("plot_themes.R");

## Load study population & merge w/ survival data:
my_samples <- loadReceptorPositiveTumors(receptorPosOnly=TRUE);
metaSurv <- loadTCGASurvivalMeta();
sum(metaSurv$patients %in% my_samples$patients)

clin.breast <- merge(
  metaSurv[ , c("patients", "OS.time","OS", "PFI.time","PFI", "DFI.time","DFI")], 
  my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER","PR","HER2","TNBC","PAM50","PAM50lite")], 
  by = "patients"
);

## For survival analysis, keep only HER2- tumors:
calcDescrpStatsByStrata(clin.breast)
clin.breast <- subset(clin.breast, HER2=="Negative");

## Exclude tumors without stage:
my_samples <- subset(my_samples, ! is.na(Stage));

## Survival modeling:
clin.breast$group <- 1 * (clin.breast$SVM_BRCA1 == "BRCA1-like"); #conver to 1/0
myPlots <- computeMultipleKaplanMeier(data=clin.breast, pTheme=mySurvTheme);
gridExtra::grid.arrange(
  myPlots$OS$plot, 
  myPlots$PFS$plot,
  myPlots$DFS$plot,
  ncol = 3
)
