---
title: "RNAseq analysis"
author: "Joseph A. Tandurella, Ludmila Danilova"
date: "2024-01-15"
output:
  html_document:
    toc: true
    toc_float: true
---

# R libraries and versions

Loading of R packages, setting seed for reproducability, and session info to provide version numbers.


```r
knitr::opts_chunk$set(echo = TRUE, message = FALSE,
                      warning = FALSE)
rm(list = ls())

library(tximport)
library('readxl')
library('DESeq2')
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following object is masked from 'package:flextable':
## 
##     width
```

```
## The following object is masked from 'package:gridExtra':
## 
##     combine
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which.max, which.min
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:Matrix':
## 
##     expand, unname
```

```
## The following objects are masked from 'package:data.table':
## 
##     first, second
```

```
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
```

```
## The following object is masked from 'package:utils':
## 
##     findMatches
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:data.table':
## 
##     shift
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## The following object is masked from 'package:grDevices':
## 
##     windows
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: MatrixGenerics
```

```
## Loading required package: matrixStats
```

```
## Warning: package 'matrixStats' was built under R version 4.3.2
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following object is masked from 'package:dplyr':
## 
##     count
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

```r
library('ComplexHeatmap')
```

```
## Loading required package: grid
```

```
## ========================================
## ComplexHeatmap version 2.16.0
## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
## Github page: https://github.com/jokergoo/ComplexHeatmap
## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
## 
## If you use it in published research, please cite either one:
## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
##     genomic data. Bioinformatics 2016.
## 
## 
## The new InteractiveComplexHeatmap package can directly export static 
## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
## 
## This message can be suppressed by:
##   suppressPackageStartupMessages(library(ComplexHeatmap))
## ========================================
## ! pheatmap() has been masked by ComplexHeatmap::pheatmap(). Most of the arguments
##    in the original pheatmap() are identically supported in the new function. You 
##    can still use the original function by explicitly calling pheatmap::pheatmap().
```

```
## 
## Attaching package: 'ComplexHeatmap'
```

```
## The following object is masked from 'package:pheatmap':
## 
##     pheatmap
```

```r
library('EnhancedVolcano')
library('DT')
library('fgsea')
library(ggplot2)
library(ggpubr)
library(tidyverse)
```

```
## Warning: package 'purrr' was built under R version 4.3.2
```

```
## Warning: package 'lubridate' was built under R version 4.3.2
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ lubridate 1.9.3     ✔ tibble    3.2.1
## ✔ purrr     1.0.2     ✔ tidyr     1.3.0
## ✔ readr     2.1.4
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ lubridate::%within%()    masks IRanges::%within%()
## ✖ data.table::between()    masks dplyr::between()
## ✖ readr::col_factor()      masks scales::col_factor()
## ✖ IRanges::collapse()      masks dplyr::collapse()
## ✖ Biobase::combine()       masks BiocGenerics::combine(), gridExtra::combine(), dplyr::combine()
## ✖ purrr::compose()         masks flextable::compose()
## ✖ matrixStats::count()     masks dplyr::count()
## ✖ IRanges::desc()          masks dplyr::desc()
## ✖ purrr::discard()         masks scales::discard()
## ✖ tidyr::expand()          masks S4Vectors::expand(), Matrix::expand()
## ✖ dplyr::filter()          masks stats::filter()
## ✖ S4Vectors::first()       masks data.table::first(), dplyr::first()
## ✖ kableExtra::group_rows() masks dplyr::group_rows()
## ✖ lubridate::hour()        masks data.table::hour()
## ✖ lubridate::isoweek()     masks data.table::isoweek()
## ✖ dplyr::lag()             masks stats::lag()
## ✖ data.table::last()       masks dplyr::last()
## ✖ lubridate::mday()        masks data.table::mday()
## ✖ lubridate::minute()      masks data.table::minute()
## ✖ lubridate::month()       masks data.table::month()
## ✖ tidyr::pack()            masks Matrix::pack()
## ✖ BiocGenerics::Position() masks ggplot2::Position(), base::Position()
## ✖ lubridate::quarter()     masks data.table::quarter()
## ✖ purrr::reduce()          masks GenomicRanges::reduce(), IRanges::reduce()
## ✖ S4Vectors::rename()      masks dplyr::rename()
## ✖ lubridate::second()      masks S4Vectors::second(), data.table::second()
## ✖ lubridate::second<-()    masks S4Vectors::second<-()
## ✖ IRanges::slice()         masks dplyr::slice()
## ✖ purrr::transpose()       masks data.table::transpose()
## ✖ tidyr::unpack()          masks Matrix::unpack()
## ✖ lubridate::wday()        masks data.table::wday()
## ✖ lubridate::week()        masks data.table::week()
## ✖ lubridate::yday()        masks data.table::yday()
## ✖ lubridate::year()        masks data.table::year()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

```r
library(dplyr)
library(tidyr)
library(data.table)
library(msigdbr)
library(fgsea)
library(Matrix)
library(ashr)
```

```
## Warning: package 'ashr' was built under R version 4.3.2
```

```r
set.seed(1234)

# set root folder
rootDir = "../"

# load shared parameters and functions
source(paste0(rootDir,"parameters_and_functions.r"))
```


# Load expression data


```r
exprsDat = read.csv(paste0(rootDir,"input_data/table_RNAseq.csv"), row.names = 1)

dim(exprsDat)
```

```
## [1] 26467    45
```


# Sample Annotations and counts


```r
#=========================
# read in patient annotation and mapping
patTab = read.csv(file = paste0(rootDir,"input_data/table_patients.csv"), row.names = 1)
# switch to R/NR annotation
patTab = patTab %>% mutate(RECIST = factor(RECIST,levels = c("Responder", "NonResponder"), labels = c("R","NR")))

datatable(patTab)
```

```
## Error in path.expand(path): invalid 'path' argument
```

```r
# sample annotation
sampAnnot = data.frame(
  patientID = paste(sapply(strsplit(colnames(exprsDat), '_'), geti,1), sapply(strsplit(colnames(exprsDat), '_'), geti,2), sep = '_'),
  timepoint = sapply(strsplit(colnames(exprsDat), '_'), geti,3), 
  check.names = F)

rownames(sampAnnot) = colnames(exprsDat)

table(sampAnnot$timepoint)
```

```
## 
## Baseline     C1D1     C2D1 
##       23       18        4
```

```r
# find pairs
mat = printPairs(sampAnnot)
```

```
## The number of paired samples between Baseline and C1D1: 16 
## The number of paired samples between Baseline and C2D1: 4 
## The number of paired samples between C2D1 and C1D1: 4
```

# Data normalization and exploratory data analysis


## Evaluation of sample quality based upon distribution of reads

We evaluate sample quality from the distribution of reads as visualized in a boxplot of log counts. 


```r
boxplot(log2(exprsDat+1), ylab='log2(counts+1)', 
        xlab='samples', las=2, cex.axis = 0.25, 
        border = ifelse(apply(log2(exprsDat+1),2,median)==0,
                        'red','black'), main = "Distribution of log-transformed counts")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

## Subset to paired samples

Due to this analysis being paired, samples with matching timepoints were retained across all comparisons. Samples without matching timepoints were removed. 


```r
##BaselinevC1D1
pairedsamples <- mat %>% filter(Baseline == 1 & C1D1 == 1) %>% rownames()
# select paired samples and exclude a time point from comparison
BaselinevC1D1sampAnnotFilter <- sampAnnot %>% filter(patientID %in% pairedsamples & timepoint != "C2D1")
# subset expression data
BaselinevC1D1exprsDatFilter <- exprsDat[,rownames(BaselinevC1D1sampAnnotFilter)]

##BaselinevC2D1
pairedsamples <- mat %>% filter(Baseline == 1 & C2D1 == 1) %>% rownames()
# select paired samples and exclude a time point from comparison
BaselinevC2D1sampAnnotFilter <-  sampAnnot %>% filter(patientID %in% pairedsamples & timepoint != "C1D1")
# subset expression data
BaselinevC2D1exprsDatFilter <- exprsDat[,rownames(BaselinevC2D1sampAnnotFilter)]

##C1D1vC2D1
pairedsamples <- mat %>% filter(C1D1 == 1 & C2D1 == 1) %>% rownames()
# select paired samples and exclude a time point from comparison
C1D1vC2D1sampAnnotFilter <- sampAnnot %>% filter(patientID %in% pairedsamples & timepoint != "Baseline")
# subset expression data
C1D1vC2D1exprsDatFilter <- exprsDat[,rownames(C1D1vC2D1sampAnnotFilter)]

##AllFinalSamples
sampAnnotPaired <- sampAnnot[rownames(sampAnnot) %in% c(rownames(BaselinevC1D1sampAnnotFilter),rownames(BaselinevC2D1sampAnnotFilter),rownames(C1D1vC2D1sampAnnotFilter)),]
exprsDatPaired <- exprsDat[,rownames(sampAnnotPaired)]
```


# Paired differential expression analysis by Time (Treatment)

For the retained samples, we run differential expression analysis across time with [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). Estimated fold changes are shrunk with `ashr` using `lfcShrink` to account for the variation in the samples in this dataset. 


```r
#DE
  #=======================
##Baseline V Entino
# get samples from baseline and C1D1
  samp <- BaselinevC1D1sampAnnotFilter[!is.na(BaselinevC1D1sampAnnotFilter$timepoint),]
  sampIds <- rownames(samp)
# create a DESeq object  
BaselinevC1D1dds <- DESeqDataSetFromMatrix(
    countData =                             round(BaselinevC1D1exprsDatFilter)[,colnames(BaselinevC1D1exprsDatFilter) %in% sampIds], 
    colData = samp, 
    design=~patientID+timepoint)

  BaselinevC1D1dds <- DESeq(BaselinevC1D1dds)
  
 # get the result table 
  BaselinevC1D1res <- results(BaselinevC1D1dds, contrast = c("timepoint","C1D1", "Baseline"))
  BaselinevC1D1res <- lfcShrink(BaselinevC1D1dds, contrast = c("timepoint","C1D1", "Baseline"), res=BaselinevC1D1res, type="ashr")
  BaselinevC1D1res <- BaselinevC1D1res[!is.na(BaselinevC1D1res$padj),]

  #=======================
##Baseline V C2D1
  samp <- BaselinevC2D1sampAnnotFilter[!is.na(BaselinevC2D1sampAnnotFilter$timepoint),]
  sampIds <- rownames(samp) 
  
  ##DE by Time
  BaselinevC2D1dds <- DESeqDataSetFromMatrix(countData =
                                  round(BaselinevC2D1exprsDatFilter)[,colnames(BaselinevC2D1exprsDatFilter) %in% sampIds], 
                              colData = samp, 
                              design=~patientID+timepoint)
  
  BaselinevC2D1dds <- DESeq(BaselinevC2D1dds)
  BaselinevC2D1res <- results(BaselinevC2D1dds, contrast = c("timepoint","C2D1", "Baseline"))
  BaselinevC2D1res <- lfcShrink(BaselinevC2D1dds, contrast = c("timepoint","C2D1", "Baseline"), res=BaselinevC2D1res, type="ashr")
  BaselinevC2D1res <- BaselinevC2D1res[!is.na(BaselinevC2D1res$padj),]

  #=======================
##Entino v C2D1
# DE 
  samp <- C1D1vC2D1sampAnnotFilter[!is.na(C1D1vC2D1sampAnnotFilter$timepoint),]
  sampIds <- rownames(samp)
  
  ##DE by Time
  C1D1vC2D1dds <- DESeqDataSetFromMatrix(countData =
                                  round(C1D1vC2D1exprsDatFilter)[,colnames(C1D1vC2D1exprsDatFilter) %in% sampIds], 
                              colData = samp, 
                              design=~patientID+timepoint)
  C1D1vC2D1dds <- DESeq(C1D1vC2D1dds)
  
  C1D1vC2D1res <- results(C1D1vC2D1dds, contrast = c("timepoint","C2D1", "C1D1"))
  C1D1vC2D1res <- lfcShrink(C1D1vC2D1dds, contrast = c("timepoint","C2D1", "C1D1"), res=C1D1vC2D1res, type="ashr")
  C1D1vC2D1res <- C1D1vC2D1res[!is.na(C1D1vC2D1res$padj),]
  
  
  #Extract Stats for Pathway Analysis
  BaselinevC1D1stats <- BaselinevC1D1res$log2FoldChange
  names(BaselinevC1D1stats) <- row.names(BaselinevC1D1res)

  BaselinevC2D1stats <- BaselinevC2D1res$log2FoldChange
  names(BaselinevC2D1stats) <- row.names(BaselinevC2D1res)
  
  C1D1vC2D1stats <- C1D1vC2D1res$log2FoldChange
  names(C1D1vC2D1stats) <- row.names(C1D1vC2D1res)
```

For each comparison, a datatable of results and volcanoplot will be provided. Plotted on the y-axis is the adjusted p-value, on the x-axis is the Log2 Fold Change.

## Baseline vs C1D1 DE 

Positive Log2 Fold Change = enriched in C1D1 samples.
Negative Log2 Fold Change = enriched in Baseline samples.

Saved are going to be the genes with an absolute Log2FoldChange > 0.5 and a Padj < 0.05


```r
##BaselinevC1D1
datatable(data.frame(BaselinevC1D1res))
```

```
## Error in path.expand(path): invalid 'path' argument
```

```r
#volcano plot
EnhancedVolcano(BaselinevC1D1res, lab=row.names(BaselinevC1D1res), 
                x='log2FoldChange', y='padj', FCcutoff = 0.5,
                pCutoff = 0.05, title = "Paired DE Results by Time -- BaselinevC1D1")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

```r
# significantly DE genes
BaselinevC1D1res_TopDE <- BaselinevC1D1res[BaselinevC1D1res$padj < 0.05 &
                              abs(BaselinevC1D1res$log2FoldChange) > 0.5,]
#order by logFC
BaselinevC1D1res_TopDE = BaselinevC1D1res_TopDE[order(BaselinevC1D1res_TopDE$log2FoldChange, decreasing = T),]

# the total number if significantly DE genes
cat("The number of significantly DE genes:\n",
nrow(BaselinevC1D1res_TopDE))
```

```
## The number of significantly DE genes:
##  1
```

```r
# count how many up- and down-regulated genes
cat("The number of down- and up-regulated genes:\n",
table(BaselinevC1D1res_TopDE$log2FoldChange >0))
```

```
## The number of down- and up-regulated genes:
##  1
```

```r
write.csv(BaselinevC1D1res_TopDE, paste0(rootDir,"output_data/Baseline_vs_C1D1_TopDE.csv"))
```

## Baseline vs C2D1 DE 

Positive Log2 Fold Change = enriched in C2D1 samples.
Negative Log2 Fold Change = enriched in Baseline samples.

Saved are going to be the genes with an absolute Log2FoldChange > 0.5 and a Padj < 0.05


```r
##BaselinevC2D1
datatable(data.frame(BaselinevC2D1res))
```

```
## Error in path.expand(path): invalid 'path' argument
```

```r
# volcano
EnhancedVolcano(BaselinevC2D1res, lab=row.names(BaselinevC2D1res), 
                x='log2FoldChange', y='padj', FCcutoff = 0.5,
                pCutoff = 0.05, title = "Paired DE Results by Time -- BaselinevC2D1")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

```r
# significantly DE genes
BaselinevC2D1res_TopDE <- BaselinevC2D1res[BaselinevC2D1res$padj < 0.05 &
                              abs(BaselinevC2D1res$log2FoldChange) > 0.5,]
#order by logFC
BaselinevC2D1res_TopDE = BaselinevC2D1res_TopDE[order(BaselinevC2D1res_TopDE$log2FoldChange, decreasing = T),]

# the total number if significantly DE genes
cat("The number of significantly DE genes:\n",
nrow(BaselinevC2D1res_TopDE))
```

```
## The number of significantly DE genes:
##  383
```

```r
# count how many up- and down-regulated genes
cat("The number of down- and up-regulated genes:\n",
table(BaselinevC2D1res_TopDE$log2FoldChange >0))
```

```
## The number of down- and up-regulated genes:
##  121 262
```

```r
write.csv(BaselinevC2D1res_TopDE, paste0(rootDir,"output_data/Baseline_vs_C2D1_TopDE.csv"))
```

## C1D1 vs C2D1 DE 

Positive Log2 Fold Change = enriched in C2D1 samples.
Negative Log2 Fold Change = enriched in C1D1 samples.

Saved are going to be the genes with an absolute Log2FoldChange > 0.5 and a Padj < 0.05


```r
##C1D1vC2D1
datatable(data.frame(C1D1vC2D1res))
```

```
## Error in path.expand(path): invalid 'path' argument
```

```r
# volcano plot
EnhancedVolcano(C1D1vC2D1res, lab=row.names(C1D1vC2D1res), 
                x='log2FoldChange', y='padj', FCcutoff = 0.5,
                pCutoffCol = 'padj', pCutoff = 0.05, title = "Paired DE Results by Time -- C1D1vC2D1")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

```r
# significantly DE genes
C1D1vC2D1res_TopDE <- C1D1vC2D1res[C1D1vC2D1res$padj < 0.05 &
                              abs(C1D1vC2D1res$log2FoldChange) > 0.5,]
#order by logFC
C1D1vC2D1res_TopDE = C1D1vC2D1res_TopDE[order(C1D1vC2D1res_TopDE$log2FoldChange, decreasing = T),]

# the total number if significantly DE genes
cat("The number of significantly DE genes:\n",
nrow(C1D1vC2D1res_TopDE))
```

```
## The number of significantly DE genes:
##  346
```

```r
# count how many up- and down-regulated genes
cat("The number of down- and up-regulated genes:\n",
table(C1D1vC2D1res_TopDE$log2FoldChange >0))
```

```
## The number of down- and up-regulated genes:
##  46 300
```

```r
write.csv(C1D1vC2D1res_TopDE, paste0(rootDir,"output_data/C1D1_vs_C2D1_TopDE.csv"))
```

# Hallmark pathway analysis

Gene set statistics are run with [fgsea](https://www.biorxiv.org/content/10.1101/060012v3) with Hallmark pathways from [MSigDb](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3106198/).

## Baseline vs C1D1 

Positive Normalized Enrichment Score = enriched in C1D1 samples.
Negative Normalized Enrichment Score = enriched in Baseline samples.



```r
# Get Pathways
hallmark_df = msigdbr(species = "Homo sapiens", category = "H")
hallmark_list = hallmark_df %>% split(x = .$gene_symbol, f = .$gs_name)

# run fgsea
  BaselinevC1D1gsResults <- as.data.frame(fgsea(pathways=hallmark_list, stats=BaselinevC1D1stats))
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error writing to connection
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error writing to connection
```

```r
#BaselinevC1D1
datatable(BaselinevC1D1gsResults)
```

```
## Error in eval(expr, envir, enclos): object 'BaselinevC1D1gsResults' not found
```

```r
fwrite(BaselinevC1D1gsResults, file=paste0(rootDir,"output_data/BaselinevC1D1_HALLMARK.csv"))
```

```
## Error in eval(expr, envir, enclos): object 'BaselinevC1D1gsResults' not found
```

## Baseline vs C2D1 

Positive Normalized Enrichment Score = enriched in C2D1 samples.
Negative Normalized Enrichment Score = enriched in Baseline samples.


```r
  #Complete Pathway Analysis
  BaselinevC2D1gsResults <- as.data.frame(fgsea(pathways=hallmark_list, stats=BaselinevC2D1stats))
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error writing to connection
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error writing to connection
```

```r
#BaselinevC2D1
datatable(BaselinevC2D1gsResults)
```

```
## Error in eval(expr, envir, enclos): object 'BaselinevC2D1gsResults' not found
```

```r
fwrite(BaselinevC2D1gsResults, file=paste0(rootDir,"output_data/BaselinevC2D1_HALLMARK.csv"))
```

```
## Error in eval(expr, envir, enclos): object 'BaselinevC2D1gsResults' not found
```

## C1D1 vs C2D1 

Positive Normalized Enrichment Score = enriched in C2D1 samples.
Negative Normalized Enrichment Score = enriched in C1D1 samples.


```r
C1D1vC2D1gsResults <- as.data.frame(fgsea(pathways=hallmark_list, stats=C1D1vC2D1stats))
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': wrong args for environment subassignment
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error writing to connection
```

```r
#C1D1vC2D1
datatable(C1D1vC2D1gsResults)
```

```
## Error in eval(expr, envir, enclos): object 'C1D1vC2D1gsResults' not found
```

```r
fwrite(C1D1vC2D1gsResults, file=paste0(rootDir,"output_data/C1D1vC2D1_HALLMARK.csv"))
```

```
## Error in eval(expr, envir, enclos): object 'C1D1vC2D1gsResults' not found
```

## HALLMARK Pathway Dotplot
Create a dotplot with GSEA results for different comparisons. Plot significant pathways only (adjusted p-values < 0.05). Color is the normalized enrichment score. Red/blue means higher/lower expression on the second time point.



```r
# merge pathwyay analysis results for plotting
  pathRes = rbind(cbind(comparison = rep("Baseline_vs_C1D1",nrow(BaselinevC1D1gsResults)),BaselinevC1D1gsResults),
  cbind(comparison = rep("Baseline_vs_C2D1",nrow(BaselinevC2D1gsResults)),BaselinevC2D1gsResults),
  cbind(comparison = rep("C1D1_vs_C2D1",nrow(C1D1vC2D1gsResults)),C1D1vC2D1gsResults))
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'nrow': object 'BaselinevC1D1gsResults' not found
```

```r
 # make NES for not sugnificant pathways as NA
  pathRes_ns = pathRes %>% mutate(NES, NES = ifelse(padj >0.05, NA, NES))
```

```
## Error in eval(expr, envir, enclos): object 'pathRes' not found
```

```r
  # plot: dot plot
  p = ggplot(data = pathRes_ns %>% filter(!is.na(NES)), aes(x = comparison, y = pathway, color = NES)) + 
    geom_point(size = 4) +
  #  scale_color_gradient(low = "red", high = "blue") +
    scale_colour_gradient2(low = "blue", high = "red") +
    scale_x_discrete(labels = c("Baseline vs C1D1", "Baseline vs C2D1","C1D1 vs C2D1")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
```

```
## Error in eval(expr, envir, enclos): object 'pathRes_ns' not found
```

```r
  print(p)
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'print': object 'p' not found
```


```r
sessionInfo()
```

```
## R version 4.3.1 (2023-06-16 ucrt)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19045)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=English_United States.utf8 
## [2] LC_CTYPE=C                            
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
##  [1] ashr_2.2-63                 msigdbr_7.5.1              
##  [3] lubridate_1.9.3             forcats_1.0.0              
##  [5] stringr_1.5.1               purrr_1.0.2                
##  [7] readr_2.1.4                 tidyr_1.3.0                
##  [9] tibble_3.2.1                tidyverse_2.0.0            
## [11] fgsea_1.26.0                EnhancedVolcano_1.18.0     
## [13] ComplexHeatmap_2.16.0       DESeq2_1.40.2              
## [15] SummarizedExperiment_1.30.2 Biobase_2.60.0             
## [17] MatrixGenerics_1.12.3       matrixStats_1.1.0          
## [19] GenomicRanges_1.52.1        GenomeInfoDb_1.36.4        
## [21] IRanges_2.34.1              S4Vectors_0.38.2           
## [23] BiocGenerics_0.46.0         tximport_1.28.0            
## [25] survival_3.5-5              binom_1.1-1.1              
## [27] survminer_0.4.9             flextable_0.9.4            
## [29] Matrix_1.6-3                kableExtra_1.3.4           
## [31] pals_1.8                    pheatmap_1.0.12            
## [33] xlsx_0.6.5                  readxl_1.4.3               
## [35] scales_1.2.1                gridExtra_2.3              
## [37] DT_0.30                     data.table_1.14.8          
## [39] ggrepel_0.9.4               ggpubr_0.6.0               
## [41] ggplot2_3.4.4               dplyr_1.1.2                
## [43] knitr_1.45                 
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.3.1           later_1.3.1             bitops_1.0-7           
##   [4] cellranger_1.1.0        lifecycle_1.0.4         mixsqp_0.3-48          
##   [7] rstatix_0.7.2           doParallel_1.0.17       lattice_0.21-8         
##  [10] crosstalk_1.2.0         backports_1.4.1         magrittr_2.0.3         
##  [13] sass_0.4.7              rmarkdown_2.25          jquerylib_0.1.4        
##  [16] yaml_2.3.7              httpuv_1.6.12           zip_2.3.0              
##  [19] askpass_1.2.0           cowplot_1.1.1           mapproj_1.2.11         
##  [22] RColorBrewer_1.1-3      maps_3.4.1.1            abind_1.4-5            
##  [25] zlibbioc_1.46.0         rvest_1.0.3             RCurl_1.98-1.13        
##  [28] xlsxjars_0.6.1          gdtools_0.3.4           circlize_0.4.15        
##  [31] GenomeInfoDbData_1.2.10 KMsurv_0.1-5            irlba_2.3.5.1          
##  [34] crul_1.4.0              svglite_2.1.2           commonmark_1.9.0       
##  [37] codetools_0.2-19        DelayedArray_0.26.7     ggtext_0.1.2           
##  [40] xml2_1.3.5              shape_1.4.6             tidyselect_1.2.0       
##  [43] httpcode_0.3.0          farver_2.1.1            webshot_0.5.5          
##  [46] jsonlite_1.8.7          GetoptLong_1.0.5        ellipsis_0.3.2         
##  [49] iterators_1.0.14        systemfonts_1.0.5       foreach_1.5.2          
##  [52] tools_4.3.1             ragg_1.2.6              snow_0.4-4             
##  [55] Rcpp_1.0.11             glue_1.6.2              xfun_0.39              
##  [58] withr_2.5.2             fastmap_1.1.1           fansi_1.0.4            
##  [61] openssl_2.1.1           truncnorm_1.0-9         digest_0.6.33          
##  [64] timechange_0.2.0        R6_2.5.1                mime_0.12              
##  [67] textshaping_0.3.7       colorspace_2.1-0        dichromat_2.0-0.1      
##  [70] markdown_1.11           utf8_1.2.3              generics_0.1.3         
##  [73] fontLiberation_0.1.0    httr_1.4.7              htmlwidgets_1.6.2      
##  [76] S4Arrays_1.0.6          pkgconfig_2.0.3         rJava_1.0-6            
##  [79] gtable_0.3.4            XVector_0.40.0          survMisc_0.5.6         
##  [82] htmltools_0.5.5         fontBitstreamVera_0.1.1 carData_3.0-5          
##  [85] clue_0.3-65             png_0.1-8               km.ci_0.5-6            
##  [88] rstudioapi_0.15.0       tzdb_0.4.0              reshape2_1.4.4         
##  [91] rjson_0.2.21            uuid_1.1-1              curl_5.1.0             
##  [94] cachem_1.0.8            zoo_1.8-12              GlobalOptions_0.1.2    
##  [97] parallel_4.3.1          pillar_1.9.0            vctrs_0.6.3            
## [100] promises_1.2.1          car_3.1-2               cluster_2.1.4          
## [103] xtable_1.8-4            evaluate_0.23           invgamma_1.1           
## [106] cli_3.6.1               locfit_1.5-9.8          compiler_4.3.1         
## [109] rlang_1.1.1             crayon_1.5.2            SQUAREM_2021.1         
## [112] ggsignif_0.6.4          labeling_0.4.3          plyr_1.8.9             
## [115] stringi_1.8.1           viridisLite_0.4.2       BiocParallel_1.34.2    
## [118] babelgene_22.9          munsell_0.5.0           fontquiver_0.2.1       
## [121] hms_1.1.3               gfonts_0.2.0            shiny_1.7.5.1          
## [124] highr_0.10              gridtext_0.1.5          broom_1.0.5            
## [127] bslib_0.5.1             fastmatch_1.1-4         officer_0.6.3
```
