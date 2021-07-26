# PASS-ON-codes

1. System requirements:
  > - R version 4.0.3 & Rstudio Version 1.3.1093 has to be installed.

2. R Dependencies and version:
 > - here 1.0.1
 > - openxlsx 4.2.3 
 > - org.Hs.eg.db 3.12.0
 > - annotate 1.68.0
 > - BiocParallel 1.24.1
 > - DESeq2 1.30.1
 > - fgsea 1.16.0
 > - GSVA 1.38.2
 > - glmnet 4.1-1
 > - pROC 1.17.0.1
 > - ROCR 1.0-11
 > - survival 3.2-7
 > - survminer 0.4.9
 > - ggplot2 3.3.3
 > - ggplot2 3.3.3
 > - ggrepel 0.9.1
 > - plotmo 3.6.0
 > - ggpubr 0.4.0
 > - RColorBrewer 1.1-2
    
3. Installation Guide:
  > Packages install from R CRAN (https://cran.r-project.org/): 

  > - install.packages(c("here", "openxlsx", "glmnet", "pROC","ROCR","survival" "survminer", "ggplot2", "plotmo", "ggpubr", "RColorBrewer"))
  
  > Packages install from Bioconductor (https://www.bioconductor.org/):

  > - if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  
  > - BiocManager::install(c("org.Hs.eg.db","annotate","BiocParallel","DESeq2", "fgsea" "GSVA"))
  
4. Other Public Signatures:
   > IMPRES: https://github.com/noamaus/IMPRES-codes
   
   > CIBERSORT: https://cibersort.stanford.edu/index.php
   
   
  
