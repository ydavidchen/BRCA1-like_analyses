# Explore clinical & survival data of METABRIC for Table 1/S1
# Script author: David Chen
# Script maintainer: David Chen
# Notes:

rm(list=ls());
library(ggplot2);
library(matrixStats);
library(reshape2);
library(tableone);
library(WriteXLS);
library(doParallel); registerDoParallel(detectCores() - 1);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("plot_themes.R");
TSIZE <- 5; 
BIN_COLORS <- c("dimgray","lightblue"); 

## Load clinical data:
sample_clinical <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030718_METABRIC_study_population.txt",sep="\t",stringsAsFactors=F);

sampNoExclu <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/031318_METABRIC_NoExclusion.txt",sep="\t",stringsAsFactors=F);
sampNoExclu$isIncludedInBRCA1analysis <- sampNoExclu$SAMPLE_ID %in% sample_clinical$SAMPLE_ID;
sampNoExclu$anyRecepPos <- sampNoExclu$ER_STATUS=="+" | sampNoExclu$PR_STATUS=="+" | sampNoExclu$HER2_STATUS=="+";

#---------------------------------------Distribution of BRCA1-like probability---------------------------------------
## Calculate proportions for showing on plot:
table(sampNoExclu$anyRecepPos, useNA="ifany");
labelRP <- "Receptor Positive (n = 1,651)"; 
labelTN <- "TNBC (n = 317)";

sampNoExclu$recepPos[sampNoExclu$anyRecepPos] <- labelRP;
sampNoExclu$recepPos[! sampNoExclu$anyRecepPos] <- labelTN;

contTab <- table(
  BRCAness = sampNoExclu$SVM_BRCA1, 
  Subtype = sampNoExclu$recepPos
);
contTab

png("~/Downloads/BRCA1ness_figures/Figure2_METABRIC.png", res=300, units="in", height=6, width=11.69);
ggplot(sampNoExclu, aes(x=reorder(SAMPLE_ID, BRCA1_Prob), y=BRCA1_Prob)) + #fill=DNAm_avail
  geom_bar(aes(fill=recepPos), stat="identity", show.legend=FALSE) + #customized
  geom_hline(yintercept=0.50, linetype="dashed") + 
  labs(x="\n METABRIC Sample", y="SVM BRCA1-like probability \n", title="METABRIC") +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(breaks=seq(0,1,0.2), expand=c(0,0), limits=c(0,1.05)) +
  scale_fill_manual(values=BIN_COLORS) +
  geom_hline(yintercept=0.50, linetype="dashed") + 
  myWaterfallTheme +
  facet_wrap(~ recepPos, scale="free")  +
  
  geom_text(data=data.frame(x=400, y=0.6, label="143 (8.7%) \n BRCA1-like", recepPos=labelRP), aes(x,y,label=label), size=TSIZE, inherit.aes=FALSE) +
  geom_text(data=data.frame(x=400, y=0.4, label="1508 (91.3%) \n non-BRCA1-like", recepPos=labelRP), aes(x,y,label=label), size=TSIZE, inherit.aes=FALSE) +
  geom_text(data=data.frame(x=50, y=0.6, label="200 (63.1%) \n BRCA1-like", recepPos=labelTN), aes(x,y,label=label), size=TSIZE, inherit.aes=FALSE) +
  geom_text(data=data.frame(x=50, y=0.4, label="117 (36.9%) \n non-BRCA1-like", recepPos=labelTN), aes(x,y,label=label), size=TSIZE, inherit.aes=FALSE)
dev.off();

#---------------------------------------Table One on overall---------------------------------------
myFactorVars <- c("Stage","ER_STATUS","PR_STATUS","HER2_STATUS","CLAUDIN_SUBTYPE","anyRecepPos","isIncludedInBRCA1analysis");
myVars <- c("AGE_AT_DIAGNOSIS", myFactorVars);
myTableOne <- CreateTableOne(
  data = sampNoExclu, 
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
y <- print(myTableOne, showAllLevels=TRUE); 
# write.csv(y, file="~/Downloads/BRCA1ness_figures/041118_METABRIC__Table_one.csv",quote=F);

#---------------------------------------------Response to BCR Reviewer #2 request---------------------------------------------
myTableOneERpos <- CreateTableOne(
  data = subset(sampNoExclu, ER_STATUS=="+"),
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
# write.csv(yERpos, file="~/Downloads/101818_METABRIC_ERpos_table1.csv",quote=FALSE);

#---------------------------------------------Response to BCR Reviewer #2, comment 7---------------------------------------------
png("~/Downloads/BRCA1ness_figures/Age_distribution_METABRIC.png", res=300, units="in", height=6, width=11.69);
ggplot(sampNoExclu, aes(AGE_AT_DIAGNOSIS, color=SVM_BRCA1, fill=SVM_BRCA1)) +
  geom_density(alpha=0.25) +
  labs(title="METABRIC", x="Age") +
  scale_y_continuous(limits=c(0,0.04)) +
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  facet_wrap(~ SVM_BRCA1) +
  myScatterTheme
dev.off();
