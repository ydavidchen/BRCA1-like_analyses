######################################################################################################
# Comparing Ki-67 RNAseqV2 gene expression levels
# Script author: David Chen
# Date: 03/14/18
# Notes:
######################################################################################################

rm(list=ls())

library(data.table)
library(gdata)
library(ggplot2)
library(matrixStats)
library(RColorBrewer); gradient_cols <- brewer.pal(12, "Paired"); 
library(reshape2)
library(splitstackshape)
library(survival)
library(survminer)
library(doParallel); registerDoParallel(detectCores() - 1)

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

table(my_samples$SVM_BRCA1)
my_samples$BRCAness[my_samples$SVM_BRCA1=="BRCA1-like"] <- "BRCA1-like (n=159)";
my_samples$BRCAness[my_samples$SVM_BRCA1=="non-BRCA1-like"] <- "non-BRCA1-like (n=595)";

png("~/Downloads/BRCA1ness_figures/Figure2C.png", res=300, units="in", height=8.27, width=5.84);
ggplot(my_samples, aes(x=BRCAness, y=MKI67, color=BRCAness)) +
  geom_boxplot(outlier.size=0, outlier.shape=0, outlier.alpha=0, position=position_dodge(width=0.775)) + #must hide outliers if jitters are to be shown
  geom_point(aes(color=BRCAness), position=position_jitterdodge(jitter.width=0.25), alpha=0.25) + # add jitter
  scale_fill_manual(values=c("mediumorchid3","darkolivegreen3")) +
  scale_color_manual(values=c("mediumorchid3","darkolivegreen3")) + 
  labs(y="MKI67 expression") + #title="TCGA"
  theme_classic() +
  theme(axis.text.x=element_text(size=14,color="black"), axis.text.y=element_text(size=20,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=20,color="black"), 
        title=element_text(size=20,color="black",face="bold"), legend.position="none",
        strip.text.x=element_text(size=12,colour="black",face="bold")) +
  
  geom_segment(aes(x=1, y=10000, xend=2, yend=10000), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1, y=9000, xend=1, yend=10000), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2, y=8000, xend=2, yend=10000), size=0.3, inherit.aes=F) +
  annotate("text", x=1.5, y=10500, label="P = 2.04e-07 ***", size=7)
dev.off();
