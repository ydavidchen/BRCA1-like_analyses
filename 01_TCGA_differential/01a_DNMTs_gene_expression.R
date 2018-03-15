######################################################################################################
# Comparing DNMT1/3A/3B RNAseqV2 gene expression levels
# Script author: David Chen
# Notes:
######################################################################################################

rm(list=ls())

library(data.table)
library(gdata)
library(matrixStats)
library(splitstackshape)

## Clinical annotation for study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);
table(my_samples$TNBC, useNA="ifany")

#---------------------------------------Load and clean RNA-seq data set---------------------------------------
rnaMedianExpr <- fread("~/Dropbox (Christensen Lab)/Breast_Cancer_Data_Sets/All_TCGA_breast_RNAseq/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt");
class(rnaMedianExpr) <- "data.frame";
rnaMedianExpr <- rnaMedianExpr[-1, ]; 
rownames(rnaMedianExpr) <- rnaMedianExpr$`Hybridization REF`; 
rnaMedianExpr <- rnaMedianExpr[ , grepl("01A", colnames(rnaMedianExpr)) | grepl("01B", colnames(rnaMedianExpr))];
rnaMedianExpr$`Hybridization REF` <- NULL;
colnames(rnaMedianExpr) <- substr(colnames(rnaMedianExpr), 1, 12); 

rnaMedianExpr <- rnaMedianExpr[ , colnames(rnaMedianExpr) %in% my_samples$patients]; #restrict

rnaMapping <- data.frame(geneName=rownames(rnaMedianExpr)); 
rnaMapping <- cSplit(rnaMapping, "geneName", sep="|", drop=F); 
class(rnaMapping) <- "data.frame"; 
rownames(rnaMapping) <- rnaMapping$geneName; 
rnaMapping$geneName <- NULL; 
colnames(rnaMapping) <- c("geneName", "entrezID");

rnaMedianExpr <- merge(rnaMedianExpr, rnaMapping[ , "geneName", drop=F], by="row.names");
rnaMedianExpr <- subset(rnaMedianExpr, geneName != "?");
if( sum(duplicated(rnaMedianExpr$geneName) != 0 ) ){
  print("Duplicate genes exist!")
  rnaMedianExpr <- subset(rnaMedianExpr, ! duplicated(geneName));
}
rownames(rnaMedianExpr) <- rnaMedianExpr$geneName; 
rnaMedianExpr$Row.names <- rnaMedianExpr$geneName <- NULL; 
rnaMedianExpr <- as.matrix(rnaMedianExpr); 
if(mode(rnaMedianExpr) != "numeric") { mode(rnaMedianExpr) <- "numeric" };
sum(rnaMedianExpr == 0)
hist(log2(rnaMedianExpr+1)) #check

#------------------------------------------------Subset to gene of interest------------------------------------------------
## Subset to gene of interest:
mRNA.chosen <- subset(rnaMedianExpr, rownames(rnaMedianExpr) %in% c("DNMT1", "DNMT3A", "DNMT3B"));
plt.DNMTs <- data.frame( t(mRNA.chosen) ); 
plt.DNMTs$patients <- rownames(plt.DNMTs);
plt.DNMTs <- merge(
  plt.DNMTs,
  my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER", "PR", "HER2","TNBC")],
  by = "patients"
);
plt.DNMTs$BRCAness[plt.DNMTs$SVM_BRCA1=="BRCA1-like"] <- 1;
plt.DNMTs$BRCAness[plt.DNMTs$SVM_BRCA1=="non-BRCA1-like"] <- 0;
plt.DNMTs$BRCAness <- as.factor(plt.DNMTs$BRCAness);
table(plt.DNMTs$BRCAness)

#------------------------------------------------Statistical tests------------------------------------------------
fitDNMT1_other <- lm(
  DNMT1 ~ BRCAness + Age + Stage + ER + PR + HER2, 
  data = subset(plt.DNMTs, TNBC=="Non-TNBC")
);
summary(fitDNMT1_other)

fitDNMT3A_other <- lm(
  DNMT3A ~ BRCAness + Age + Stage + ER + PR + HER2, 
  data = subset(plt.DNMTs, TNBC=="Non-TNBC")
);
summary(fitDNMT3A_other)

fitDNMT3B_other <- lm(
  DNMT3B ~ BRCAness + Age + Stage + ER + PR + HER2, 
  data = subset(plt.DNMTs, TNBC=="Non-TNBC")
);
summary(fitDNMT3B_other)
