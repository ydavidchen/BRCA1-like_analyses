## Helper Functions for SVM BRCA1-like Downstream Statistical Analyses
## Script author: David Chen
## Script maintainer: David Chen
## Note: Each script in the `src` subdirectory is like a "main", but to be executed interactively.

library(doParallel); 
registerDoParallel(detectCores() - 1);

#-------------------------------------Data-loading methods-------------------------------------
loadReceptorPositiveTumors <- function(path="~/repos/BRCA1ness_by_SVM/annotations_and_backups/030418_TCGA_study_population.txt", 
                                       receptorPosOnly=TRUE) {
  #'@description Load clinical annotation for TCGA breast tumors that are positive for ER, PR, and/or HER2
  #'@param path Path to clinical annotation saved in TXT format
  #'@param receptorPosOnly Logical, defaults to TRUE. Should data be subsetted to ER+, PR+, and/or HER2+ only?
  my_samples <- read.csv(path, header=TRUE, sep="\t", stringsAsFactors=FALSE);
  if(receptorPosOnly){ 
    my_samples <- subset(my_samples, TNBC=="Non-TNBC");
  }
  return(my_samples); 
}

loadHRDmetrics <- function(path="~/Dropbox (Christensen Lab)/Pan-cancer-analyses/Marquard AM pan-cancer HRD index/40364_2015_33_MOESM8_ESM.txt") {
  #'@description Load existing TCGA HR deficiency metrics published in Marquard et al. 2015
  #'@param path Path to authors' data file saved in TXT format
  MarquardHRD <- read.table(path, header=TRUE);
  MarquardHRD <- subset(MarquardHRD, Tumor %in% my_samples$patients);
  colnames(MarquardHRD)[colnames(MarquardHRD)=="Tumor"] <- "patients";
  return(MarquardHRD);
}

loadMutationMeasures <- function(path="~/Dropbox (Christensen Lab)/Pan-cancer-analyses/Kandoth2013Nature/Supplementary_Table_3a.xls") {
  #'@description Load mutational count & rates data published in Kandoth et al. 2013 Nature
  #'@param path Path to authors' data file saved in Excel format
  require(gdata);
  mutMeasure <- read.xls(path, stringsAsFactors=FALSE);
  mutMeasure <- subset(mutMeasure, Cancer.Type == "brca");
  mutMeasure$patients <- substr(mutMeasure$TCGA.ID, 1, 12);
  stopifnot( anyDuplicated(mutMeasure$patients) == 0 );
  return(mutMeasure);
}

loadMutationalSignatures <- function(path="~/Dropbox (Christensen Lab)/Pan-cancer-analyses/TCGA_pan_cancer_signature_Rosenthal.txt") {
  #'@description Load TCGA-BRCA mutational signatures published in Rosenthal et al. 2016
  #'@param path Path to authors' data file saved in TXT format
  rosenthalTCGA <- read.table(path, header=TRUE);
  rosenthalTCGA$id <- gsub(".", "-", rosenthalTCGA$id, fixed=TRUE);
  rosenthalTCGA <- subset(rosenthalTCGA, cancer=="BRCA");
  colnames(rosenthalTCGA)[1] <- "patients"; #for merging
  return(rosenthalTCGA);
}

loadTCGAGeneLeveExpr <- function(path="~/Dropbox (Christensen Lab)/Breast_Cancer_Data_Sets/All_TCGA_breast_RNAseq/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt") {
  #'@description Load TCGA-BRCA gene expression
  #'@param path Path to TCGA breast cancer gene expression count data in TXT format from FireBrowse
  require(data.table);
  require(splitstackshape);
  rnaMedianExpr <- fread(path);
  class(rnaMedianExpr) <- "data.frame";
  rnaMedianExpr <- rnaMedianExpr[-1, ];
  rownames(rnaMedianExpr) <- rnaMedianExpr$`Hybridization REF`; 
  rnaMedianExpr <- rnaMedianExpr[ , grepl("01A", colnames(rnaMedianExpr)) | grepl("01B", colnames(rnaMedianExpr))];
  rnaMedianExpr$`Hybridization REF` <- NULL;
  colnames(rnaMedianExpr) <- substr(colnames(rnaMedianExpr), 1, 12); 
  
  rnaMapping <- data.frame(geneName=rownames(rnaMedianExpr)); 
  rnaMapping <- cSplit(rnaMapping, "geneName", sep="|", drop=F); 
  class(rnaMapping) <- "data.frame";
  rownames(rnaMapping) <- rnaMapping$geneName;
  rnaMapping$geneName <- NULL;
  colnames(rnaMapping) <- c("geneName", "entrezID");
  
  rnaMedianExpr <- merge(rnaMedianExpr, rnaMapping[ , "geneName", drop=FALSE], by="row.names");
  rnaMedianExpr <- subset(rnaMedianExpr, geneName != "?");
  if( sum(duplicated(rnaMedianExpr$geneName) != 0 )) {
    rnaMedianExpr <- subset(rnaMedianExpr, ! duplicated(geneName));
  }
  rownames(rnaMedianExpr) <- rnaMedianExpr$geneName; 
  rnaMedianExpr$Row.names <- rnaMedianExpr$geneName <- NULL; 
  rnaMedianExpr <- as.matrix(rnaMedianExpr); 
  if(mode(rnaMedianExpr) != "numeric") { 
    mode(rnaMedianExpr) <- "numeric"; 
  };
  return(rnaMedianExpr);
}

loadTCGAmiRNA <- function(path="~/Dropbox (Christensen Lab)/Breast_Cancer_Data_Sets/All_TCGA_breast_RNAseq/bcgsc.ca_BRCA_IlluminaGA_miRNASeq.miRNAExp.tsv") {
  #'@description Load TCGA primary breast tumor miRNA data & convert to Z scores
  #'@param path Path to TCGA miRNA expression data downloaded from SynapseTCGAlive
  require(data.table); 
  miRExpr <- fread(path, stringsAsFactors=FALSE)
  class(miRExpr) <- "data.frame";
  rownames(miRExpr) <- miRExpr$V1; 
  miRExpr$V1 <- NULL;
  miRExpr <- as.matrix(miRExpr);
  stopifnot( mode(miRExpr) == "numeric" );
  miRExpr <- miRExpr[ , substr(colnames(miRExpr), 14, 15)=="01"];
  colnames(miRExpr) <- substr(colnames(miRExpr), 1, 12);
  miRExpr <- scale(miRExpr, scale=TRUE, center=TRUE);
  return(miRExpr);
}

loadTCGASurvivalMeta <- function(path="~/Dropbox (Christensen Lab)/Pan-cancer-analyses/Liu2018_Cell_TCGA_surv/Liu2018_Cell_TCGA-CDR.csv") {
  #'@description Load TCGA pan-cancer survival data from Liu et al. 2018 Cell
  #'@param path Path to authors' data file saved in CSV format
  metaSurv <- read.csv(path, stringsAsFactors=FALSE);
  metaSurv[metaSurv=="#N/A"] <- NA;
  metaSurv[metaSurv=="[Not Applicable]"] <- NA;
  
  metaSurv$OS <- as.logical(as.numeric(metaSurv$OS)); 
  metaSurv$OS.time <- as.numeric(metaSurv$OS.time) / 365.25; 
  
  metaSurv$PFI <- as.logical(as.numeric(metaSurv$PFI)); 
  metaSurv$PFI.time <- as.numeric(metaSurv$PFI.time) / 365.25; 
  
  metaSurv$DFI <- as.logical(as.numeric(metaSurv$DFI));
  metaSurv$DFI.time <- as.numeric(metaSurv$DFI.time) / 365.25;
  
  metaSurv$patients <- toupper(metaSurv$bcr_patient_barcode);
  return(metaSurv); 
}

loadTCGA4CoxOS <- function(path.surv="~/repos/BRCA1ness_by_SVM/annotations_and_backups/TCGA-BRCA_clinical.csv", 
                           ADMIN_CENSOR=NULL, HER2negOnly=TRUE, receptorPosOnly=TRUE) {
  #'@description Load TCGA breast tumor annotation for survival analysis
  #'@param path.surv Path to TCGA data with overall survival information
  #'@param ADMIN_CENSOR Time for administrative censoring in months
  data <- read.csv(path.surv, stringsAsFactors=FALSE);
  colnames(data)[1] <- "patients";
  
  my_samples <- loadReceptorPositiveTumors(receptorPosOnly=receptorPosOnly); 
  if(HER2negOnly) {
    my_samples <- subset(my_samples, HER2=="Negative");
  }
  my_samples <- subset(my_samples, ! (is.na(Stage) | is.na(Age)) );
  my_samples <- merge(
    my_samples[ , c("patients","SVM_BRCA1","Age","Stage","ER","PR","HER2","TNBC","PAM50","PAM50lite")],
    data[ , c("patients","days_to_last_follow_up","days_to_death","vital_status")]
  );
  
  ## For patients who died, fill in days to last follow-up:
  condition <- is.na(my_samples$days_to_last_follow_up) & my_samples$vital_status=="dead";
  my_samples$days_to_last_follow_up[condition] <- my_samples$days_to_death[condition];
  
  ## Define main variable of interest:
  my_samples$group[my_samples$SVM_BRCA1=="BRCA1-like"] <- 1;
  my_samples$group[my_samples$SVM_BRCA1=="non-BRCA1-like"] <- 0;
  my_samples$group <- as.factor(my_samples$group);
  
  ## Define main outcome of interest:
  my_samples$event <- (my_samples$vital_status=="dead");
  
  ## Important: Convert days_to_follow-up to months:
  my_samples$OS_MONTHS <- my_samples$days_to_last_follow_up / 30.44;
  my_samples$days_to_last_follow_up <- my_samples$days_to_death <- NULL;
  
  ## Impose administrative censoring:
  if(! is.null(ADMIN_CENSOR)) {
    my_samples$event[my_samples$OS_MONTHS >= ADMIN_CENSOR] <- FALSE; #not dead by definition
    my_samples$OS_MONTHS[my_samples$OS_MONTHS >= ADMIN_CENSOR] <- ADMIN_CENSOR;
  }
  
  return(my_samples);
}


loadMETABRICtumors <- function(path="~/repos/BRCA1ness_by_SVM/annotations_and_backups/030718_METABRIC_study_population.txt",
                               receptorPosOnly=TRUE) {
  #'@description Load METABRIC breast tumor clinical annotation
  #'@param path Path to METABRIC data file saved in TXT format downloaded via cBioPortal
  sample_clinical <- read.csv(path, sep="\t",stringsAsFactors=FALSE);
  if(receptorPosOnly) {
    sample_clinical <- subset(sample_clinical, TNBC=="Non-TNBC"); 
  }
  return(sample_clinical); 
}

loadMETABRIC4CoxOS <- function(ADMIN_CENSOR=NULL, receptorPosOnly=TRUE, HER2negOnly=TRUE) {
  #'@description Load METABRIC breast tumor annotation for survival analysis
  #'@param ADMIN_CENSOR Time for administrative censoring in months
  #'@param receptorPosOnly Logical, defaults to TRUE. Should data be subsetted to ER+, PR+, and/or HER2+ only?
  #'@param HER2negOnly Logical, defaults to TRUE. Should HER2-positive/unknown tumors be excluded?
  sample_clinical <- loadMETABRICtumors(receptorPosOnly=receptorPosOnly);
  if(HER2negOnly) {
    sample_clinical <- subset(sample_clinical, (ER_STATUS=="+" | PR_STATUS=="+") & HER2_STATUS=="-");
    table(sample_clinical$ER_STATUS, sample_clinical$PR_STATUS, useNA="ifany");
  }
  
  ## Re-code main variable of interest:
  sample_clinical$group[sample_clinical$SVM_BRCA1=="BRCA1-like"] <- 1;
  sample_clinical$group[sample_clinical$SVM_BRCA1=="non-BRCA1-like"] <- 0;
  sample_clinical$group <- as.factor(sample_clinical$group);
  
  ## Define main outcome of interest:
  sample_clinical$event <- (sample_clinical$OS_STATUS=="DECEASED");
  
  ## Impose administrative censoring:
  if(! is.null(ADMIN_CENSOR)) {
    sample_clinical$event[sample_clinical$OS_MONTHS >= ADMIN_CENSOR] <- FALSE; #not dead by definition
    sample_clinical$OS_MONTHS[sample_clinical$OS_MONTHS >= ADMIN_CENSOR] <- ADMIN_CENSOR;    
  }
  
  return(sample_clinical); 
}

loadMETABRICarrayExpr <- function(path="~/Dropbox (Christensen Lab)/Breast_Cancer_Data_Sets/METABRIC_data_set_cBio/data_expression.txt") {
  #'@description Load METABRIC breast tumor normalized gene expression TXT file
  #'@param path Path to METABRIC gene expression data saved in TXT format downloaded via cBioPortal
  require(data.table);
  arrayExpr <- fread(path);
  class(arrayExpr) <- "data.frame";
  anyDuplicated(arrayExpr$Hugo_Symbol)
  rownames(arrayExpr) <- arrayExpr$Hugo_Symbol;
  arrayExpr$Hugo_Symbol <- arrayExpr$Entrez_Gene_Id <- NULL;
  arrayExpr <- as.matrix(arrayExpr);
  stopifnot( mode(arrayExpr) == "numeric" )
  return(arrayExpr);
}

loadCombinedTCGAandMETABRIC <- function(receptorPosOnly, HER2negOnly, ADMIN_CENSOR=NULL) {
  #'@description Combine TCGA and METABRIC data set for survival analysis. 
  #'@details See loadTCGA4CoxOS and loadMETABRIC4CoxOS for detail
  sampsTCGA <- loadTCGA4CoxOS(ADMIN_CENSOR=ADMIN_CENSOR, receptorPosOnly=TRUE, HER2negOnly=TRUE);
  sampsTCGA <- sampsTCGA[ , c("patients","SVM_BRCA1","group","Age","Stage","ER","PR","HER2","event","OS_MONTHS")]; 
  
  sampsMETABRIC <- loadMETABRIC4CoxOS(ADMIN_CENSOR=ADMIN_CENSOR, receptorPosOnly=TRUE, HER2negOnly=TRUE);
  sampsMETABRIC <- sampsMETABRIC[ , c("SAMPLE_ID","SVM_BRCA1","group","AGE_AT_DIAGNOSIS","Stage","ER_STATUS","PR_STATUS","HER2_STATUS","event","OS_MONTHS")];
  sampsMETABRIC[sampsMETABRIC=="+"] <- "Positive";
  sampsMETABRIC[sampsMETABRIC=="-"] <- "Negative";
  colnames(sampsMETABRIC) <- colnames(sampsTCGA);
  
  sampsTCGA$dataset <- "TCGA";
  sampsMETABRIC$dataset <- "METABRIC";
  sampsCombined <- rbind(sampsTCGA, sampsMETABRIC);
  
  return(sampsCombined); 
}

#-------------------------------------Analytical Methods / Pipelines-------------------------------------
checkAgeStrata <- function(data, thresh, runFisher=FALSE, ...) {
  #'@description Check binary age strata in combined TCGA+METABRIC.
  #'@param data R data.frame with a numeric column named "Age"
  #'@param thresh Numeric threshold for age stratification, in years
  #'@param runFisher Logical, defaults to FALSE. Should Fisher's exact test be run?
  #'@param ... Additional parameters passed to fisher.test
  s <- table(data$SVM_BRCA1, data$Age > thresh);
  colnames(s)[colnames(s)=="FALSE"] <- paste("Younger than", thresh, "y");
  colnames(s)[colnames(s)=="TRUE"] <- paste("Older than", thresh, "y");
  if(runFisher) {
    print( fisher.test(s, ...) ); 
  }
  ## Marginal sums:
  s <- cbind(s, Total=rowSums(s));
  s <- rbind(s, Total=colSums(s));
  return(s);
}

calcDescrpStatsByStrata <- function(data, ageThresh=50) {
  #'@description Print descriptive statistics by age strata, stage, BRCAness, ER/PR/HER2
  #'@param data R data.frame with a numeric column named "Age", "Stage", "SVM_BRCA1", "ER", "PR", "HER2"
  #'@param ageThresh Numeric threshold for age stratification, in years
  print(paste("Age greater than",ageThresh))
  print(table(data$Age > ageThresh))
  
  print("Stage:")
  print(table(data$Stage, useNA="ifany"))
  
  print("ER:")
  print(table(data$ER, useNA="ifany"))
  
  print("PR:")
  print(table(data$PR, useNA="ifany"))
  
  print("HER2:")
  print(table(data$HER2, useNA="ifany"))
  
  print("SVM BRCA1-like status:")
  print(table(data$SVM_BRCA1, useNA="ifany"))
}

computeMultipleKaplanMeier <- function(data, pTheme, myColors=c("darkolivegreen3","mediumorchid"), showTab=FALSE) {
  #'@description Kaplan-Meier curves using a list of survival models
  #'@param data R data.frame consisting of 0/1 for group and time in months for OS, PFI, DFI
  #'@param pTheme ggplot theme. Specified in the module plot_themes.R
  #'@param myColors Colors for each strata compared
  #'@param showTab Logical, defaults to FALSE. Should risk table be shown?
  require(survival);
  require(survminer);
  sModels <- list();
  sModels[["OS"]] <- survfit(Surv(time=OS.time, event=OS) ~ group, data);
  sModels[["PFS"]] <- survfit(Surv(time=PFI.time, event=PFI) ~ group, data);
  sModels[["DFS"]] <- survfit(Surv(time=DFI.time, event=DFI) ~ group, data);
  
  plotList <- list(); 
  for(n in names(sModels)) {
    m <- sModels[[n]];
    pTitle <- n;
    plotList[[n]] <- ggsurvplot(
      m,
      data = data,
      ## Risk Table:
      risk.table.col = "strata",
      risk.table = showTab,
      tables.theme = theme_cleantable(),
      ## Plot:
      pval = TRUE,
      pval.size = 7,
      conf.int = FALSE,
      conf.int.style = "step",
      censor = FALSE,
      size = 1.5,
      alpha = 0.7,
      ## Labels:
      legend = c(0.75, 0.15), #rel. position
      xlab = "Time (years)", 
      title = pTitle,
      ## Global elements:
      ggtheme = pTheme,
      palette = myColors
    );
  }
  return(plotList);
}
