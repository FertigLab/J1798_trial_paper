# Code for paper

Code to publish with a paper by M. Baretti et al. "Immunomodulation of the tumor microenvironment of pancreatic ductal adenocarcinoma with histone deacetylase inhibition: results of a phase 2 clinical trial of entinostat in combination with nivolumab".

# Overview

This repository contains input data and code to generate figures that are used in the paper referenced above. All input data are stored in the `input_data` folder. The Rmarkdown files are stored in the `Rmd` folder. The main script is `runRmd.r` that installs all required packages and renders Rmd files to generate output tables in the `output_data` folder and html/docs reports in `Rmd`. 

# System requirements

The code can be run on any computer with installed the R programming enviroment. Depending on your system, installation of all required packages can take up to an hour.

# How to run

Open the R enviroment and run the following line:

```
source('runRmd.r')

```
The code was tested on Windows 10 Pro with 32.0 GB RAM and the 2.81 GHz processor. If all required packages are installed, the runtime of the script is about 15 minutes. 