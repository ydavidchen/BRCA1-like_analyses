# TCGA breast tumor characteristics by SVM-BRCA1-like status
# Script author: David Chen
# Script maintainer: David Chen
# Notes:

rm(list=ls());
library(tableone);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
source("plot_themes.R");
TSIZE <- 5; 
BIN_COLORS <- c("dimgray","lightblue"); 

## Clinical annotation for study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);

mySampNoExclu <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/031318_Complete_Sample_Sheet_No_Exclusion.txt", header=T, sep="\t", stringsAsFactors=F);
mySampNoExclu$includedInBRCA1Analysis <- mySampNoExclu$patients %in% my_samples$patients;
mySampNoExclu$anyRecepPos <- mySampNoExclu$ER=="Positive" |  mySampNoExclu$PR=="Positive" | mySampNoExclu$HER2=="Positive";

#---------------------------------------Distribution of BRCA1-like probability---------------------------------------
table(mySampNoExclu$anyRecepPos, useNA="ifany");
labelRP <- "Receptor Positive (n = 754)"; 
labelTN <- "TNBC (n = 100)";

mySampNoExclu$recepPos[mySampNoExclu$anyRecepPos] <- labelRP;
mySampNoExclu$recepPos[! mySampNoExclu$anyRecepPos] <- labelTN;
mySampNoExclu$recepPos[is.na(mySampNoExclu$anyRecepPos)] <- "Unknown";

contTab <- table(
  BRCAness = mySampNoExclu$SVM_BRCA1, 
  Subtype = mySampNoExclu$recepPos
);
contTab

png("~/Downloads/BRCA1ness_figures/Figure2_TCGA.png", res=300, units="in", height=6, width=11.69);
ggplot(subset(mySampNoExclu, recepPos != "Unknown"), aes(x=reorder(patients, BRCA1_prob), y=BRCA1_prob)) +
  geom_bar(aes(fill=recepPos), stat="identity", show.legend=FALSE) +
  geom_hline(yintercept=0.50, linetype="dashed") + 
  labs(x="\n TCGA Sample", y="SVM BRCA1-like probability \n", title="TCGA") +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(breaks=seq(0,1,0.2), expand=c(0,0), limits=c(0,1.05)) +
  scale_fill_manual(values=BIN_COLORS) +
  myWaterfallTheme  +
  facet_wrap(~ recepPos, scale="free") +
  
  geom_text(data=data.frame(x=150, y=0.6, label="159 (21.1%) \n BRCA1-like", recepPos=labelRP), aes(x,y,label=label), size=TSIZE, inherit.aes=FALSE) +
  geom_text(data=data.frame(x=150, y=0.4, label="595 (78.9%) \n non-BRCA1-like", recepPos=labelRP), aes(x,y,label=label), size=TSIZE, inherit.aes=FALSE) +
  geom_text(data=data.frame(x=20, y=0.6, label="88 (88.0%) \n BRCA1-like", recepPos=labelTN), aes(x,y,label=label), size=TSIZE, inherit.aes=FALSE) +
  geom_text(data=data.frame(x=20, y=0.4, label="12 (12.0%) \n non-BRCA1-like", recepPos=labelTN), aes(x,y,label=label), size=TSIZE, inherit.aes=FALSE)
dev.off();

#---------------------------------------------Table One: Overall---------------------------------------------
## Table One
myFactorVars <- c("Stage","ER","PR","HER2","TNBC","PAM50", "race", "anyRecepPos", "includedInBRCA1Analysis");
myVars <- c("Age", myFactorVars);
myTableOne <- CreateTableOne(
  data = mySampNoExclu, 
  vars = myVars, 
  strata = "SVM_BRCA1",
  factorVars = myFactorVars,
  smd = TRUE,
  test = TRUE,
  testExact = fisher.test,
  argsExact = list(workspace = 2*10^5), #default
  testNormal = oneway.test,
  argsNormal = list(var.equal=TRUE),
  includeNA = TRUE
);

## Visualize & export:
y <- print(myTableOne, showAllLevels=TRUE);
# write.csv(y, file="~/Downloads/051118_TCGA-BRCA_No-exclusion_TableOne.csv",quote=F);
# WriteXLS(mySampNoExclu, ExcelFileName="~/Downloads/031318_Sample_Sheet_TCGA_NoExclu.xls");

## Response to Reviewer #2 comment 6:
myTableOneERpos <- CreateTableOne(
  data = subset(mySampNoExclu, ER=="Positive"), 
  vars = myVars, 
  strata = "SVM_BRCA1",
  factorVars = myFactorVars,
  smd = TRUE,
  test = TRUE,
  testExact = fisher.test,
  argsExact = list(workspace = 2*10^5), #default
  testNormal = oneway.test,
  argsNormal = list(var.equal=TRUE),
  includeNA = TRUE
);
yERpos <- print(myTableOneERpos, showAllLevels=TRUE);
# write.csv(yERpos, file="~/Downloads/101818_TCGA-BRCA_ER+BRCA.csv",quote=FALSE);

#---------------------------------------------Table One: Subset---------------------------------------------
## Table One
myFactorVars <- c("Stage","ER","PR","HER2");
myVars <- c("Age", myFactorVars);
myTableOne <- CreateTableOne(
  data = subset(mySampNoExclu, TNBC=="Non-TNBC" & includedInBRCA1Analysis),
  vars = myVars, 
  strata = "SVM_BRCA1",
  factorVars = myFactorVars,
  smd = FALSE, 
  test = FALSE,
  includeNA = TRUE
);

## Visualize & export:
y <- print(myTableOne, showAllLevels=TRUE); 
# write.csv(y, file="~/Downloads/BRCA1ness_figures/041118_TCGA_Subset_Table_one.csv",quote=F);

