######################################################################################################
# Comparing Ki-67 RNAseqV2 gene expression levels
# Script author: David Chen
# Notes: MKI67 gene
######################################################################################################

rm(list=ls())

library(data.table)
library(gdata)
library(matrixStats)
library(splitstackshape)

my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/031318_Complete_Sample_Sheet_No_Exclusion.txt", header=T, sep="\t", stringsAsFactors=F);
my_samples <- subset(my_samples, TNBC=="Non-TNBC");

my_samples$group[my_samples$SVM_BRCA1=="BRCA1-like"] <- 1;
my_samples$group[my_samples$SVM_BRCA1=="non-BRCA1-like"] <- 0;
my_samples$group <- as.factor(my_samples$group);

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

## Subset to gene of interest:
ki67 <- subset(rnaMedianExpr, rownames(rnaMedianExpr) %in% "MKI67");
ki67 <- data.frame( t(ki67) ); 
ki67$patients <- rownames(ki67);

my_samples <- merge(
  my_samples,
  ki67,
  by = "patients"
);

#------------------------------------------------Statistical test------------------------------------------------
fitki67_other <- lm(
  MKI67 ~ group + Age + Stage + ER + PR + HER2, 
  data = my_samples
);
summary(fitki67_other)

