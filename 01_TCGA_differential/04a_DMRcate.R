######################################################################################################
# DMRcate analysis of TCGA Non-TNBCs
# Script author: David Chen
# Notes:
######################################################################################################

rm(list=ls())

library(DMRcate)
library(matrixStats)
library(doParallel); registerDoParallel(detectCores() - 1)

## 450K annotation:
library(IlluminaHumanMethylation450kanno.ilmn12.hg19);
annot.450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19); 
sexProbes  <- annot.450k$Name[annot.450k$chr %in% c("chrX", "chrY")];
XProbes <- annot.450k$Name[annot.450k$chr=="chrX"]; #chrX probes

#-----------------------------------------------Load data for the study population-----------------------------------------------
## Load clinical data & define study population:
my_samples <- read.csv("~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", header=T, sep="\t", stringsAsFactors=F);
my_samples <- subset(my_samples, TNBC=="Non-TNBC");
sum(duplicated(my_samples$patients))

## Drop samples with missing covariates for which need to be adjusted:
my_samples <- subset(my_samples, ! (is.na(Age) | is.na(Stage) | is.na(ER) | is.na(PR)  | is.na(HER2) ) );

## Main variable of interest:
my_samples$group[my_samples$SVM_BRCA1=="BRCA1-like"] <- 1;
my_samples$group[my_samples$SVM_BRCA1=="non-BRCA1-like"] <- 0;
my_samples$group <- as.factor(my_samples$group);
table(my_samples$group)
table(my_samples$SVM_BRCA1)

## Load beta-values, apply constraints, and convert to M-values:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/betas.TCGAC450k.RData");
betas_TCGA450k[betas_TCGA450k <= 0] <- 0.000001;
betas_TCGA450k[betas_TCGA450k >= 1] <- 0.999999; 
breast_mvals <- logit2(betas_TCGA450k);

## Mutually subset clinical & 450K data:
breast_mvals <- breast_mvals[ , colnames(breast_mvals) %in% my_samples$patients];
my_samples <- subset(my_samples, patients %in% colnames(breast_mvals));

table(rowSums(is.na(breast_mvals))) #ensure no missing value

## Separate chrX & autosomal probes:
tcga_auto_mval <- breast_mvals[! rownames(breast_mvals) %in% sexProbes, ];
tcga_chrX_mval <- breast_mvals[rownames(breast_mvals) %in% XProbes, ];
rm(betas_TCGA450k, breast_mvals, targets) #reduce memory usage

#----------------------------------------------Execute DMRcate----------------------------------------------
## Create design matrix:
if(identical(colnames(tcga_auto_mval), colnames(tcga_chrX_mval)) &
   identical(colnames(tcga_auto_mval), my_samples$patients) ){
  print("Conditions met! Proceed to assemble a common design matrix to be used for autosomes and chrX...")
  my_design <- model.matrix( ~ group + Age + Stage + ER + PR + HER2, data=my_samples) #same as limma setup 
} else {
  stop("Sample IDs failed to match! Design matrix not created.")
}

which(colnames(my_design) == "group1")

tcga_annot_auto <- cpg.annotate(
  datatype = "array",
  object = tcga_auto_mval,
  what = "M",
  arraytype = "450K",
  analysis.type = "differential", 
  design = my_design, 
  coef = which(colnames(my_design) == "group1"),
  fdr = 0.05 #individual probe FDR
);
## Your contrast returned 5996 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

tcga_annot_chrX <- cpg.annotate(
  datatype = "array", 
  object = tcga_chrX_mval,
  what = "M",
  arraytype = "450K",
  analysis.type = "differential", 
  design = my_design, 
  coef = which(colnames(my_design) == "group1"),
  fdr = 0.05 #individual probe FDR
);
## Your contrast returned 917 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

dmrcate_auto <- dmrcate(tcga_annot_auto, min.cpgs=10, C=2);
dmrcate_chrX <- dmrcate(tcga_annot_chrX, min.cpgs=10, C=2);

#----------------------------------------------Extract results and export----------------------------------------------
## Extract DMRs:
autosomal_ranges <- extractRanges(dmrcate_auto, genome="hg19");
chrX_ranges <- extractRanges(dmrcate_chrX, genome="hg19");
res.dmrcate <- rbind(
  as.data.frame(autosomal_ranges), 
  as.data.frame(chrX_ranges)  
);
rownames(res.dmrcate) <- NULL;
res.dmrcate$isRegionSignif <- res.dmrcate$Stouffer < 0.05;

sum(res.dmrcate$seqnames!="chrX" & res.dmrcate$isRegionSignif) #152
sum(res.dmrcate$seqnames=="chrX" & res.dmrcate$isRegionSignif) #50

## Extract DMPs/CpGs:
## (Note: Do not subset based on statistical significance yet: CpG universe info needed)
dmps <- rbind(dmrcate_auto$input, dmrcate_chrX$input);
sum(dmps$indfdr < 0.05) #6913

## Export objects:
rm("annot.450k","sexProbes","XProbes","tcga_auto_mval","tcga_chrX_mval"); 
save(
  list = ls(),
  file = "~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/030418_DMRcate_Non-TNBCs.RData", 
  compress = TRUE
); 