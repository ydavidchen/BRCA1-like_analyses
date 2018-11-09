# Differentially methylated CpGs in TCGA hormone receptor-positive BRCA1-like tumors
# Script author: David Chen
# Script maintainer: David Chen
# Notes:

rm(list=ls());
library(DMRcate);
library(matrixStats);
library(pheatmap);
library(reshape2);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("plot_themes.R");

HEAT_COLORS <- colorRampPalette(c("yellow","black","blue"))(2048); 
PATH_450K <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/betas.TCGAC450k.RData"; 
PATH_DMRcate_RES <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/030418_DMRcate_Non-TNBCs.RData";
BETA_CUT <- 0.125;
FDR_CUT <- 0.05; 

## Load data & results:
load(PATH_450K);
load(PATH_DMRcate_RES);

## 450K annotation:
library(IlluminaHumanMethylation450kanno.ilmn12.hg19);
annot.450k <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19));
annot.450k$isEnhancer <- ifelse(annot.450k$Enhancer=="TRUE" | annot.450k$Phantom != "", "Yes", "No");
annot.450k$isPromoter <- ifelse(grepl("TSS", annot.450k$UCSC_RefGene_Group), "Yes", "No");

#----------------------------------------------Heat map of individually significant CpGs----------------------------------------------
sum(dmps$betafc > BETA_CUT); sum(dmps$betafc < -BETA_CUT);

## Rank-ordered mean beta FC distribution:
ranked_DMPs <- sort(dmps$betafc[dmps$indfdr < FDR_CUT]);
plot(
  ranked_DMPs, 
  col = ifelse(abs(ranked_DMPs) > BETA_CUT, "red", "black"),
  pch=16, bty='l', cex=0.3, 
  ylab = "Beta-value fold change", 
  main = "Distribution of beta-value FC of individually significant CpGs used for DMRcate (TCGA, n=322)"
);
abline(h=c(-BETA_CUT, BETA_CUT), lty=2);
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
  Stage = c(`Stage I-II`="lightgray", `Stage III-IV`="black"),
  HER2 = c(Positive="black", Negative="lightgray"),
  PR = c(Positive="black", Negative="lightgray"),
  ER = c(Positive="black", Negative="lightgray"),
  
  Promoter = c(Yes="black", No="lightgray"),
  Enhancer = c(Yes="black", No="lightgray"),
  DNase = c(Yes="black", No="lightgray"),
  Context = c(Island="black", Shore="gray30", Shelf="gray70", OpenSea="gray90")
);

myDMPs <- as.character(dmps$ID[dmps$indfdr < FDR_CUT & abs(dmps$betafc) >= BETA_CUT]);
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
  color = HEAT_COLORS,
  fontsize = 10, #labels,titles
  border_color = NA
);

#----------------------------------Response to BCR Reviewer #1, point d----------------------------------
table(my_samples$HER2, seNA="ifany")
her2neg_tumors <- my_samples$patients[my_samples$HER2=="Negative"];
pheatmap(
  mat[ , colnames(mat) %in% her2neg_tumors],
  cutree_cols = 2,
  show_rownames = FALSE, #CpGs
  show_colnames = FALSE, #samples
  annotation_col = heat_annot,
  annotation_row = row_annot,
  annotation_colors = ann_colors,
  color = HEAT_COLORS,
  fontsize = 9, #labels,titles
  border_color = NA,
  main = "ER+/PR+, HER2- Tumors Only"
);

#----------------------------------Response to BCR Reviewer #1, point e----------------------------------
## step 1. Extract clusters
ph <- pheatmap(
  mat,
  cutree_cols = 2,
  silent = TRUE
);
sampleClust <- as.data.frame(cutree(ph$tree_col, k=2));
colnames(sampleClust) <- "Cluster";
table(sampleClust$Cluster); 
sampInBLClust <- rownames(sampleClust)[sampleClust$Cluster == 1];

## step 2. Compute intersample Var(beta) for each cluster
intSampVarBL <- rowVars(mat[ , colnames(mat) %in% sampInBLClust]);
intSampVarNB <- rowVars(mat[ , ! colnames(mat) %in% sampInBLClust]);

## step 3a. Compute average
avgVarBL <- mean(intSampVarBL);
avgVarNBL <- mean(intSampVarNB);

## step 3b. Plot rank-ordered inter-sample variance distributions 
par(mar=c(5,5,4,2));
plot(
  sort(intSampVarBL, decreasing=TRUE), 
  xlab = "CpG index", ylab="Var(beta-value)", 
  cex = 0.75,
  cex.lab = 1.5,
  ylim = c(0,0.15), 
  bty = "l", 
  col = "mediumorchid3",
  pch = 16
);
points(sort(intSampVarNB, decreasing=TRUE), col="darkolivegreen3", cex=0.75, pch=16);
abline(h=avgVarBL, col="mediumorchid3", lwd=2, lty=2);
abline(h=avgVarNBL, col="darkolivegreen3", lwd=2, lty=2);
legend(
  150, 0.15, pch=16, bty="n", cex=1.5,
  legend = c(paste("BRCA1-like Cluster, mean =", round(avgVarBL,3)), paste("non-BRCA1-like Cluster, mean =", round(avgVarNBL,3) )), 
  col = c("mediumorchid3","darkolivegreen3")
);

# intSampVar <- data.frame(
#   CpG = rownames(mat),
#   BRCA1_like = rowVars(mat[ , colnames(mat) %in% sampInBLClust]),
#   non_BRCA1_like = rowVars(mat[ , ! colnames(mat) %in% sampInBLClust])
# );
# intSampVar <- melt(intSampVar, variable.name="SVM_BRCA1", value.name="Variance");
# intSampVar$SVM_BRCA1 <- gsub("_", "-", intSampVar$SVM_BRCA1); 
# 
# ggplot(intSampVar, aes(reorder(CpG, Variance), Variance)) +
#   geom_point(aes(color=SVM_BRCA1)) +
#   facet_wrap(~ SVM_BRCA1) +
#   myWaterfallTheme
