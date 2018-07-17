[![Build Status](https://travis-ci.org/ydavidchen/BRCA1-like_analyses.svg?branch=master)](https://travis-ci.org/ydavidchen/BRCA1-like_analyses)

# BRCA1-like Analyses

Downstream bioinformatics and biostatistical investigations of BRCA1-like hormone receptor-positive breast tumors identified by a Support Vector Machine (SVM) copy number classifier

## Overview

The molecular features of homologous recombination (HR)-deficient, BRCA1-like breast tumors as a group is not well characterized. The gene *BRCA1* plays a role in HR DNA repair. *BRCA1*-mutated tumors, as well as some tumors without a *bona fide BRCA1* mutation, can demonstrate HR deficiency and characteristic copy number profiles. Collectively, these tumors are referred to as "BRCA1-like". If identified, these tumors might be selectively and efficiently targeted by certain cancer therapies.

Here, we compared the clinical and molecular profiles of BRCA1-like ER/PR/HER2-positive breast tumors with their non-BRCA1-like counterpart (equivalent to "controls"). Our work provides insights into alternative classification approaches and therapeutic-targeting strategies for this heterogeneous disease subgroup.

## File Structure

All analysis scripts are in the `src` sub-directory.

* `helper_functions.R`: Helper functions for data loading and some analytical approaches
* `plot_themes.R`: List of themes for `ggplot2` that could be conveniently loaded, used, and modified based on specific data-visualization need
* All other scripts (numbered 01 to 04) are for downstream statistical analysis of The Cancer Genome Atlas (TCGA) or Curtis et al./Pereia et al. breast tumor clinical and molecular profiles. Each script resembles a *main* method executed interactively
* Algorithm for copy number genomic mapping and binary classification are available from co-authors upon request.

Please responsibly use materials associated with this repository and provide relevant citations of the manuscript wherever applicable.

Copyright &copy; 2018 ydavidchen and Christensen Lab
