---
title: "package bugower tutorial"
author: "Liu Cao"
date: "2017-05-16"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{bugpower tutorial}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---

##Introduction
The goal of `bugpower` is to provide a simple and convenient set of R functions to help perform power analysis for different types of microbiome datasets, including metagenomic taxa profiles, DNA pathways profiles in terms of metagenomics shot gun data as well as normalized RNA pathway profiles in terms of metatranscriptomics shot gun data. Given a study design, the package will mainly output a power calculation report for the minimum difference of relative abundance a study can detect, the power curve and false positive rate curve.

The following brief tutorial will outline the main functions employed in `bugpower` and provide some simple examples of how to properly execute the functions. Generally, with a single main function 'write_report', the user can get a .pdf power calculation report.

##Generating the report
With 'write_report', users can specify their own study design. The following code specifies: 
*the output directory of the report with 'outputdir'
*file name of the .Rnw file and .pdf file with 'filename'
*the author name in the report with 'author'
*the number of samples features, covariate included in the analysis with 'N_sample', 'N_feature' and 'N_covarite' respectively,  
*the type covariate of interest in the study with 'contin',
*the percentage of continuous covariates in the simulation with conti_prop,
*the percentage of positive features in the simulation with pos_prop,
*and the smoothness parameter that controls controls the number of piont in the power plot (small increment means large number of point).

```
library(bugpower) #Load package
write_report(outputdir = "/Users/cello/Desktop/",
             filename= "test",
             author="Liu Cao",
             datatype="species",
             N_sample=100,
             N_feature=300,
             N_covariate=8,
             N_repeat=1,
             contin="discrete",
             conti_prop=0.5,
             pos_prop=0.2,
             increment=0.005
             )
```



##More unitilities provided in the package


##References
