################################################################################################
# Differentially methylated CpGs in TCGA Non-TNBCs
# Script author: David Chen
# Notes:
################################################################################################

rm(list=ls())

library(DMRcate)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(matrixStats)
library(pheatmap)
library(splitstackshape)
library(WriteXLS)

library(RColorBrewer)
gradient_cols <- brewer.pal(12, "Paired")
cols <- brewer.pal(8, "Set1"); #use in numeric order

## Load data & results:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/betas.TCGAC450k.RData");
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/030418_DMRcate_Non-TNBCs.RData");

## 450K annotation:
library(IlluminaHumanMethylation450kanno.ilmn12.hg19);
annot.450k <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19));
annot.450k$isEnhancer <- ifelse(annot.450k$Enhancer=="TRUE" | annot.450k$Phantom != "", "Yes", "No");
annot.450k$isPromoter <- ifelse(grepl("TSS", annot.450k$UCSC_RefGene_Group), "Yes", "No");

#----------------------------------------------Heat map of individually significant CpGs----------------------------------------------
beta_cut <- 0.125;
sum(dmps$betafc > beta_cut); sum(dmps$betafc < -beta_cut)

## Ranked mean beta FC plot:
ranked_DMPs <- sort(dmps$betafc[dmps$indfdr < 0.05]);
plot(
  ranked_DMPs, 
  col = ifelse(abs(ranked_DMPs) > beta_cut, "red", "black"),
  pch=16, bty='l', cex=0.3, 
  ylab = "Beta-value fold change", 
  main = "Distribution of beta-value FC of individually significant CpGs used for DMRcate (TCGA, n=322)"
);
abline(h=c(-beta_cut, beta_cut), lty=2);
text(6000, 0.2, "175 CpGs", col=2);
text(1500, -0.2, "216 CpGs", col=2);

## Heat map:
heat_annot <- data.frame(
  row.names = my_samples$patients,
  SVM_BRCA1 = my_samples$SVM_BRCA1,
  Stage = my_samples$Stage,
  HER2 = my_samples$HER2,
  PR = my_samples$PR, 
  ER = my_samples$ER
);
row_annot <- data.frame(
  row.names = annot.450k$Name,
  Context = annot.450k$Relation_to_Island,
  Promoter = annot.450k$isPromoter,
  Enhancer = annot.450k$isEnhancer,
  DNase = ifelse(annot.450k$DHS=="TRUE", "Yes", "No")
);
row_annot$Context <- gsub("S_", "", row_annot$Context);
row_annot$Context <- gsub("N_", "", row_annot$Context);
colnames(row_annot) <- gsub(".", " ", colnames(row_annot), fixed=TRUE);
ann_colors <- list(
  SVM_BRCA1 = c(`BRCA1-like`="mediumorchid3", `non-BRCA1-like`="darkolivegreen3"),
  Stage = c(`Stage I-II`="lavender", `Stage III-IV`="darkred"),
  PAM50 = c(Basal="purple",  Her2="aquamarine1", LumA="burlywood1", LumB="coral1", Normal="darkgreen"),
  HER2 = c(Positive="black", Negative="lightgray"),
  PR = c(Positive="black", Negative="lightgray"),
  ER = c(Positive="black", Negative="lightgray"),
  
  Promoter = c(Yes="black", No="lightgray"),
  Enhancer = c(Yes="black", No="lightgray"),
  DNase = c(Yes="black", No="lightgray"),
  Context = c(Island="black", Shore="mediumorchid3", Shelf="orange", OpenSea="darkgray")
);
myDMPs <- as.character(dmps$ID[dmps$indfdr < 0.05 & abs(dmps$betafc) >= beta_cut]);
mat <- betas_TCGA450k[rownames(betas_TCGA450k) %in% myDMPs, 
                      colnames(betas_TCGA450k) %in% my_samples$patients];

pheatmap(
  mat,
  cutree_cols = 2,
  show_rownames = FALSE, #CpGs
  show_colnames = FALSE, #samples
  annotation_col = heat_annot,
  annotation_row = row_annot,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow", "blue"))(2048),
  fontsize = 9, #labels,titles
  border_color = NA
);
