# script to run all Rdm files to generate figure for the paper by 
# M. Baretti et al. "Immunomodulation of the tumor microenvironment 
# of pancreatic ductal adenocarcinoma with histone deacetylase inhibition: 
#results of a phase 2 clinical trial of entinostat in combination with nivolumab".

# a list of packages required to run Rmd files
packageList = c("knitr","dplyr","ggplot2", "ggpubr", "ggrepel","DT",
                "gridExtra","pals","kableExtra",'DESeq2','ComplexHeatmap',
                'EnhancedVolcano', 'fgsea', 'msigdbr','survminer','ashr',
                'flextable','binom', 'survival')
# create a list of packages to be installed
pkgsToInstall = c()
for(p in packageList) if(!require(p,character.only = TRUE, quietly = TRUE)) pkgsToInstall = c(pkgsToInstall,p)

# if there are packages to install, do that 
if (!is.null(pkgsToInstall) & !require("BiocManager", quietly = TRUE))
{
    install.packages("BiocManager")
        BiocManager::install(pkgs = pkgsToInstall)
}
for(p in pkgsToInstall) library(p, character.only = TRUE, quietly = TRUE)

# list all rmd files
files = list.files(path = "Rmd/", pattern = ".Rmd", full.names = T)
# render all rmds
for(i in files) rmarkdown::render(i, output_dir = outputFolder)
