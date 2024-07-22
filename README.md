# Overview

This repository contains code to generate figures for the paper by M. Baretti et al. "Entinostat in combination with nivolumab in metastatic pancreatic ductal adenocarcinoma: a phase 2 clinical trial" accepted for publication in Nature Communications. 

## Input data
Most of the input data are available in the `input_data` folder. Restrictions apply to the availability of patient data and so the 
`table_patients.csv` file is not publicly available. However, the patient data files which included limited deidentified clinical data associated with the manuscript will be available upon reasonable request to the corresponding authors and will be accessable through dbGaP (phs003615). 

## Rmarkdown
Rmarkdown files are stored in the `Rmd` folder. The main script is `runRmd.r` which installs all required packages and renders Rmd files to generate output tables in the `output_data` folder and html/docs reports in `Rmd`. 

## Output data
The `output_data` folder contains html/docs reports rendered from Rmarkdown files, as well as output tables with the analysis results.

# System requirements

The code can be run on any computer with installed the R programming environment. Depending on your system, installation of all required packages can take up to an hour.

# How to run

Open the R environment and run the following line:

```
source('runRmd.r')
```
The code was tested on Windows 10 Pro with 32.0 GB RAM and a 2.81 GHz processor. If all required packages are installed, the runtime of the script is about 10 minutes. The script was tested in R version 4.3.2 and the following packages:
```
## R version 4.3.2 (2023-10-31 ucrt)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19045)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=English_United States.utf8 
## [2] LC_CTYPE=English_United States.utf8   
## [3] LC_MONETARY=English_United States.utf8
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.utf8    
## 
## time zone: America/New_York
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] pheatmap_1.0.12             scales_1.3.0               
##  [3] Matrix_1.6-3                data.table_1.14.8          
##  [5] BiocManager_1.30.22         survival_3.5-7             
##  [7] binom_1.1-1.1               flextable_0.9.4            
##  [9] survminer_0.4.9             msigdbr_7.5.1              
## [11] fgsea_1.26.0                EnhancedVolcano_1.18.0     
## [13] ComplexHeatmap_2.16.0       DESeq2_1.40.2              
## [15] SummarizedExperiment_1.30.2 Biobase_2.60.0             
## [17] MatrixGenerics_1.12.3       matrixStats_1.1.0          
## [19] GenomicRanges_1.52.1        GenomeInfoDb_1.36.4        
## [21] IRanges_2.34.1              S4Vectors_0.38.2           
## [23] BiocGenerics_0.46.0         kableExtra_1.3.4           
## [25] pals_1.8                    gridExtra_2.3              
## [27] DT_0.30                     ggrepel_0.9.4              
## [29] ggpubr_0.6.0                ggplot2_3.4.4              
## [31] dplyr_1.1.2                 knitr_1.45                 
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.3.2           later_1.3.1             bitops_1.0-7           
##   [4] tibble_3.2.1            lifecycle_1.0.4         rstatix_0.7.2          
##   [7] mixsqp_0.3-54           doParallel_1.0.17       lattice_0.21-9         
##  [10] crosstalk_1.2.1         backports_1.4.1         magrittr_2.0.3         
##  [13] sass_0.4.8              rmarkdown_2.25          jquerylib_0.1.4        
##  [16] yaml_2.3.7              httpuv_1.6.12           zip_2.3.0              
##  [19] askpass_1.2.0           cowplot_1.1.2           mapproj_1.2.11         
##  [22] RColorBrewer_1.1-3      maps_3.4.1.1            abind_1.4-5            
##  [25] zlibbioc_1.46.0         rvest_1.0.3             purrr_1.0.2            
##  [28] RCurl_1.98-1.13         gdtools_0.3.4           circlize_0.4.15        
##  [31] GenomeInfoDbData_1.2.10 KMsurv_0.1-5            irlba_2.3.5.1          
##  [34] crul_1.4.0              svglite_2.1.2           commonmark_1.9.0       
##  [37] codetools_0.2-19        DelayedArray_0.26.7     ggtext_0.1.2           
##  [40] xml2_1.3.5              tidyselect_1.2.0        shape_1.4.6            
##  [43] httpcode_0.3.0          farver_2.1.1            webshot_0.5.5          
##  [46] jsonlite_1.8.7          GetoptLong_1.0.5        ellipsis_0.3.2         
##  [49] iterators_1.0.14        systemfonts_1.0.5       foreach_1.5.2          
##  [52] tools_4.3.2             ragg_1.2.6              snow_0.4-4             
##  [55] Rcpp_1.0.11             glue_1.6.2              xfun_0.39              
##  [58] withr_2.5.2             fastmap_1.1.1           fansi_1.0.4            
##  [61] openssl_2.1.1           digest_0.6.33           truncnorm_1.0-9        
##  [64] R6_2.5.1                mime_0.12               textshaping_0.3.7      
##  [67] colorspace_2.1-0        dichromat_2.0-0.1       markdown_1.12          
##  [70] utf8_1.2.3              tidyr_1.3.0             generics_0.1.3         
##  [73] fontLiberation_0.1.0    httr_1.4.7              htmlwidgets_1.6.4      
##  [76] S4Arrays_1.0.6          pkgconfig_2.0.3         gtable_0.3.4           
##  [79] XVector_0.40.0          survMisc_0.5.6          htmltools_0.5.7        
##  [82] fontBitstreamVera_0.1.1 carData_3.0-5           clue_0.3-65            
##  [85] png_0.1-8               ashr_2.2-63             km.ci_0.5-6            
##  [88] rstudioapi_0.15.0       reshape2_1.4.4          rjson_0.2.21           
##  [91] uuid_1.1-1              curl_5.1.0              zoo_1.8-12             
##  [94] cachem_1.0.8            GlobalOptions_0.1.2     stringr_1.5.1          
##  [97] parallel_4.3.2          pillar_1.9.0            vctrs_0.6.3            
## [100] promises_1.2.1          car_3.1-2               xtable_1.8-4           
## [103] cluster_2.1.4           evaluate_0.23           invgamma_1.1           
## [106] cli_3.6.1               locfit_1.5-9.8          compiler_4.3.2         
## [109] rlang_1.1.1             crayon_1.5.2            SQUAREM_2021.1         
## [112] ggsignif_0.6.4          labeling_0.4.3          plyr_1.8.9             
## [115] stringi_1.8.1           viridisLite_0.4.2       BiocParallel_1.34.2    
## [118] babelgene_22.9          munsell_0.5.0           fontquiver_0.2.1       
## [121] gfonts_0.2.0            shiny_1.8.0             highr_0.10             
## [124] gridtext_0.1.5          broom_1.0.5             bslib_0.6.1            
## [127] fastmatch_1.1-4         officer_0.6.3

````
