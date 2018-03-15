###############################################################################################################
# Test **enriched** genomic contexts, stratified by direction of change of differential Non-TNBC CpGs
# Script author: David Chen
# Notes:
###############################################################################################################

rm(list=ls())

library(ggplot2)
library(gridExtra)
library(matrixStats)
library(reshape2)
library(doParallel); registerDoParallel(detectCores() - 1)

## Load data & results:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/BRCA1ness_TCGA_all_types/030418_DMRcate_Non-TNBCs.RData");

## 450K annotation:
library(IlluminaHumanMethylation450kanno.ilmn12.hg19);
annot.450k <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19));
annot.450k$isEnhancer <- annot.450k$Enhancer=="TRUE" | annot.450k$Phantom != "";
annot.450k$isPromoter <- grepl("TSS", annot.450k$UCSC_RefGene_Group);

#------------------------------------------------Data/results loading------------------------------------------------
DMPs <- dmps;
colnames(DMPs) <- gsub("ID", "Name", colnames(DMPs)); #for merging
DMPs <- merge(DMPs, annot.450k, by="Name");

## Define differential set:
pThresh <- 0.05;
esThresh <- 0.125;
DMPs$isConsideredDifferential[DMPs$betafc > esThresh & DMPs$indfdr < pThresh] <- "HYPER";
DMPs$isConsideredDifferential[DMPs$betafc < -esThresh & DMPs$indfdr < pThresh] <- "HYPO";
DMPs$isConsideredDifferential[is.na(DMPs$isConsideredDifferential)] <- FALSE;

## Add 0/1 indicator columns for transcriptional context:
DMPs$Promoter <- ifelse(DMPs$isPromoter, 1, 0);
DMPs$Enhancer_for_enrichment <- ifelse(DMPs$isEnhancer, 1, 0); #450K only
DMPs$Gene.body <- ifelse(grepl("Body",DMPs$UCSC_RefGene_Group), 1, 0);
DMPs$DNase <- ifelse(DMPs$DHS=="TRUE", 1, 0);

## Add 0/1 indicator columns for CGI (each category is mutually exclusive of others):
DMPs$Island <- ifelse(DMPs$Relation_to_Island=="Island", 1, 0);
DMPs$OpenSea <- ifelse(DMPs$Relation_to_Island=="OpenSea", 1, 0);
DMPs$Shore <- ifelse(DMPs$Relation_to_Island %in% c("N_Shore","S_Shore"), 1, 0);
DMPs$Shelf <- ifelse(DMPs$Relation_to_Island %in% c("N_Shelf", "S_Shelf"), 1, 0);

#------------------------------------------------Stratified Forest Plots------------------------------------------------
table(DMPs$isConsideredDifferential, useNA='ifany') #quick check

plotList <- list();
for(directionInCF in c("202 HYPOmethylated","148 HYPERmethylated") ){ #"350 differentially methylated"
  ## Enriched subset & associated color:
  ## DIFFERENT FOR EACH ITERATION!
  if(directionInCF == "148 HYPERmethylated"){
    DMPs.sig <- subset(DMPs, isConsideredDifferential=="HYPER");
    myColor <- "darkred";
  } else if(directionInCF == "202 HYPOmethylated") {
    DMPs.sig <- subset(DMPs, isConsideredDifferential=="HYPO");
    myColor <- "navy";
  } else if(directionInCF == "350 differentially methylated") {
    DMPs.sig <- subset(DMPs, isConsideredDifferential != "FALSE") #both directions
    myColor <- "black";
  }
  ## Contingency table for Promoter CpGs:
  promoterFisher <- rbind(
    Differential = table(DMPs.sig$Promoter),
    Input = table(DMPs$Promoter)
  );
  promoterFisher <- promoterFisher[ , c(2,1)]; #Set `1` (yes) as first column
  
  ## Contingency table for Enhancer CpGs:
  enhancerFisher <- rbind(
    Differential = table(DMPs.sig$Enhancer_for_enrichment),
    Input = table(DMPs$Enhancer_for_enrichment)
  );
  enhancerFisher <- enhancerFisher[ , c(2,1)]; #Set `1` (yes) as first column
  
  ## Contingency table for DNase I Hypersensitivity Sites:
  dhsFisher <- rbind(
    Differential = table(DMPs.sig$DNase),
    Input = table(DMPs$DNase)
  );
  dhsFisher <- dhsFisher[ , c(2,1)]; #Set `1` (yes) as first column
  
  ## Contingency table for Island CpGs:
  cgiFisher <- rbind(
    Differential = table(DMPs.sig$Island),
    Input = table(DMPs$Island)
  );
  cgiFisher <- cgiFisher[ , c(2,1)]; #Set `1` (yes) as first column
  
  ## Contingency table for Open-Sea CpGs:
  OpenSeaFisher <- rbind(
    Differential = table(DMPs.sig$OpenSea),
    Input = table(DMPs$OpenSea)
  ); 
  OpenSeaFisher <- OpenSeaFisher[ , c(2,1)]; #Set `1` (yes) as first column

  ## Fisher tests:
  x <- fisher.test(promoterFisher);
  y <- fisher.test(enhancerFisher);
  z <- fisher.test(dhsFisher);
  w <- fisher.test(cgiFisher);
  u <- fisher.test(OpenSeaFisher);
  
  ## Summary of test results (order of items in each column MUST be consistent!)
  summFisher <- data.frame(
    row.names = c("Promoter", "Enhancer", "DNase", "Island", "OpenSea"),
    OR = c(x$estimate, y$estimate, z$estimate, w$estimate, u$estimate),
    CI.lower = c(x$conf.int[[1]], y$conf.int[[1]], z$conf.int[[1]], w$conf.int[[1]], u$conf.int[[1]]),
    CI.upper = c(x$conf.int[[2]], y$conf.int[[2]], z$conf.int[[2]], w$conf.int[[2]], u$conf.int[[2]]),
    P = c(x$p.value, y$p.value, z$p.value, w$p.value, u$p.value)
  );

  ## Forest ggplot:
  plt.summFisher <- summFisher; #copy
  plt.summFisher$Category <- factor( rownames(plt.summFisher), levels=rev(c("Promoter","Enhancer","DNase","Island","OpenSea")) );
  plt.summFisher$P <- signif(plt.summFisher$P, 3);
  plotList[[directionInCF]] <- ggplot(plt.summFisher, aes(x=Category, y=OR, ymin=CI.lower, ymax=CI.upper)) +
    geom_pointrange(size=1, color=myColor) +
    geom_hline(yintercept=1.00,size=0.2,linetype="dashed",color=myColor) +
    coord_flip() + #order: left to right becomes bottom to up
    scale_y_continuous(limits=c(-0.25,4.5)) + #consistent scale for all
    theme_classic() +
    theme(axis.line=element_line(color=myColor), axis.ticks=element_line(color=myColor),
          axis.text=element_text(size=20,color=myColor), title=element_text(size=20,face="bold",color=myColor),
          axis.title.x=element_text(size=20,face="bold",color=myColor), axis.title.y=element_blank()) +
    labs(x="", y="Odds Ratio (OR)", title=paste(directionInCF,"CpGs in BRCA1-like Tumors")) +
    annotate("text", 5.15, plt.summFisher$OR[1], label=paste("P =", plt.summFisher$P[1]),size=7,color=myColor) +
    annotate("text", 4.15, plt.summFisher$OR[2], label=paste("P =", plt.summFisher$P[2]),size=7,color=myColor) +
    annotate("text", 3.15, plt.summFisher$OR[3], label=paste("P =", plt.summFisher$P[3]),size=7,color=myColor) +
    annotate("text", 2.15, plt.summFisher$OR[4], label=paste("P =", plt.summFisher$P[4]),size=7,color=myColor) +
    annotate("text", 1.15, plt.summFisher$OR[5], label=paste("P =", plt.summFisher$P[5]),size=7,color=myColor)
}
grid.arrange(grobs=plotList, ncol=2)
