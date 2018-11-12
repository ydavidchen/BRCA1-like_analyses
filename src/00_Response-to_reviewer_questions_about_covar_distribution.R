# Response to Reviewer Comments related to datasets and covariate distributions
# Script author: David Chen
# Script maintainer: David Chen
# Notes:

rm(list=ls());
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");
source("plot_themes.R");

TCGA_ALL_TUMOR_PATH  <- "~/repos/BRCA1ness_by_SVM/annotations_and_backups/031318_Complete_Sample_Sheet_No_Exclusion.txt"; 
METABRIC_ALL_TUMOR_PATH <- "~/repos/BRCA1ness_by_SVM/annotations_and_backups/031318_METABRIC_NoExclusion.txt"; 

combine_TCGA_and_metabric <- function(path.tcga, path.metabric) {
  ## TCGA:
  mySampNoExclu_tcga <- read.csv(path.tcga, sep="\t", header=TRUE, stringsAsFactors=FALSE);
  mySampNoExclu_tcga$anyBRCAalt <- mySampNoExclu_tcga$BRCA1.mut | mySampNoExclu_tcga$BRCA2.mut | mySampNoExclu_tcga$BRCA1meth; 
  mySampNoExclu_tcga <- mySampNoExclu_tcga[ , c("patients","BRCA1_prob","SVM_BRCA1","Age","ajcc_pathologic_tumor_stage","Stage","anyBRCAalt")]; 
  mySampNoExclu_tcga$misclassif <- mySampNoExclu_tcga$anyBRCAalt & mySampNoExclu_tcga$SVM_BRCA1=="non-BRCA1-like";
  mySampNoExclu_tcga$dataset <- "TCGA"; 
  
  ## METABRIC:
  sampNoExclu_metabric <- read.csv(path.metabric, sep="\t", stringsAsFactors=FALSE);
  sampNoExclu_metabric$anyBRCAalt <- sampNoExclu_metabric$BRCA1mut | sampNoExclu_metabric$BRCA2mut; 
  sampNoExclu_metabric <- sampNoExclu_metabric[ , c("SAMPLE_ID","BRCA1_Prob","SVM_BRCA1","AGE_AT_DIAGNOSIS","TUMOR_STAGE","Stage","anyBRCAalt")]; 
  sampNoExclu_metabric$misclassif <- sampNoExclu_metabric$anyBRCAalt & sampNoExclu_metabric$SVM_BRCA1=="non-BRCA1-like";
  sampNoExclu_metabric$dataset <- "METABRIC";
  
  ## Combine data sets to increase statistical power:
  colnames(sampNoExclu_metabric) <- colnames(mySampNoExclu_tcga);
  combinedTCGA_METABRIC <- rbind(mySampNoExclu_tcga, sampNoExclu_metabric);
  
  combinedTCGA_METABRIC$ternaryStage[combinedTCGA_METABRIC$ajcc_pathologic_tumor_stage %in% c("1","Stage I","Stage IA","Stage IB")] <- "Stage I";
  combinedTCGA_METABRIC$ternaryStage[combinedTCGA_METABRIC$ajcc_pathologic_tumor_stage %in% c("2","Stage II","Stage IIA","Stage IIB")] <- "Stage II";
  combinedTCGA_METABRIC$ternaryStage[combinedTCGA_METABRIC$ajcc_pathologic_tumor_stage %in% c("3","4","Stage III","Stage IIIA", "Stage IIIB","Stage IIIC","Stage IV")] <- "Stage III-IV";
  
  return(combinedTCGA_METABRIC)
}

eval_BRCAalt_vs_BRCAness <- function(sampAnnot) {
  contTab <- table(
    BRCAalt = sampAnnot$anyBRCAalt,
    BRCAness = sampAnnot$SVM_BRCA1
  );
  contTab <- contTab[c(2,1), ];
  print(contTab);
  print(fisher.test(contTab));
}

## Overall summaries:
combinedTCGA_METABRIC <- combine_TCGA_and_metabric(TCGA_ALL_TUMOR_PATH, METABRIC_ALL_TUMOR_PATH);
eval_BRCAalt_vs_BRCAness(combinedTCGA_METABRIC)

contTabStage <- table(
  Stage = combinedTCGA_METABRIC$Stage, 
  BRCA.alt = combinedTCGA_METABRIC$anyBRCAalt,
  useNA = "no"
);
contTabStage <- contTabStage[ , c(2,1)];
contTabStage
fisher.test(contTabStage)

#---------------------------------------------Response to Reviewer #1, comment g---------------------------------------------
contTabRace <- matrix(c(60,61, (201+20+1), (488+36+0)), byrow=TRUE, ncol=2);
contTabRace
sum(contTabRace)
fisher.test(contTabRace)

#---------------------------------------------Response to Reviewer #2, comment 7---------------------------------------------
## Stage, based on numbers reported in Supplemental Table 3:
contTabStageByDataset <- matrix(c((63+165),(26+102),(238+473),(220+1106)), byrow=TRUE, ncol=2);
contTabStageByDataset
sum(contTabStageByDataset)
fisher.test(contTabStageByDataset)

contTabStageTCGA <- matrix(c(63,165,238,473), byrow=TRUE, ncol=2);
contTabStageTCGA
fisher.test(contTabStageTCGA)

contTabStageMETABRIC <- matrix(c(26,102,220,1106), byrow=TRUE, ncol=2);
contTabStageMETABRIC
fisher.test(contTabStageMETABRIC)

## Age
fitAgeByDataset <- lm(Age ~ dataset, data=combinedTCGA_METABRIC);
summary(fitAgeByDataset)

ggplot(combinedTCGA_METABRIC, aes(x=dataset, Age)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0) +
  geom_point(position=position_jitter(width=0.25), alpha=0.3) + 
  labs(y="Age in years") +
  geom_segment(aes(x=1, y=98, xend=2, yend=98), size=0.3, inherit.aes=FALSE) +
  annotate("text", x=1.5, y=100.5, label="*** P = 2.52E-07", size=7) + 
  myBoxplotTheme

#------------------------------------------Response to BCR Reviewer #2 Comment 5------------------------------------------
## Restrict to non-BRCA1-like tumors
nonbrca_joined <- subset(combinedTCGA_METABRIC, SVM_BRCA1 == "non-BRCA1-like");
table(nonbrca_joined$anyBRCAalt, useNA="ifany")

## Compare continuous variables:
## 1) BRCA1-like probability scores:
fit.prob <- lm(BRCA1_prob ~ anyBRCAalt, data=nonbrca_joined);
summary(fit.prob);

## 2) Subject chronological age:
fit.age  <- lm(Age ~ anyBRCAalt, data=nonbrca_joined);
summary(fit.age); 

## Data visualization:
nonbrca_joined$BRCA.alteration <- ifelse(nonbrca_joined$anyBRCAalt, "BRCA1/2 altered \n (n=64)", "No evidence of \n BRCA1/2 alteration (n=2210)");
ggplot(nonbrca_joined, aes(x=BRCA.alteration, BRCA1_prob)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0) +
  geom_point(position=position_jitter(width=0.25), alpha=0.3) + 
  labs(y="BRCA1-like probability score", title="All non-BRCA1-like tumors (n=2,274)") +
  myBoxplotTheme

ggplot(nonbrca_joined, aes(x=BRCA.alteration, Age)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0) +
  geom_point(position=position_jitter(width=0.25), alpha=0.3) + 
  labs(y="Age (years)", title="All non-BRCA1-like tumors (n=2,274)") +
  myBoxplotTheme

## Compare tumor stage:
ct_tcga_nbl <- table(
  Stage = nonbrca_joined$ternaryStage,
  BRCA.alteration = nonbrca_joined$anyBRCAalt,
  useNA = "no"
);
ct_tcga_nbl <- ct_tcga_nbl[ , c(2,1)];
ct_tcga_nbl
fisher.test(ct_tcga_nbl)
