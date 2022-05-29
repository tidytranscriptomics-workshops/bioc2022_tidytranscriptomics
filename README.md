<!-- badges: start -->
[![DOI](https://zenodo.org/badge/379767139.svg)](https://zenodo.org/badge/latestdoi/379767139)
[![.github/workflows/basic_checks.yaml](https://github.com/tidytranscriptomics-workshops/iscb2021_tidytranscriptomics/workflows/.github/workflows/basic_checks.yaml/badge.svg)](https://github.com/tidytranscriptomics-workshops/iscb2021_tidytranscriptomics/actions) 	
<!-- badges: end -->

# Introduction to Tidy Transcriptomics
<p float="left">
<img height="100" width="300" alt="iscbacademy" src="man/figures/ISCBacademy.png"/>
<img height="100" alt="tidybulk" src="https://github.com/Bioconductor/BiocStickers/blob/master/tidybulk/tidybulk.png?raw=true"/>
</p>

## Instructor names and contact information

* Maria Doyle <Maria.Doyle at petermac.org>  
* Stefano Mangiola <mangiola.s at wehi.edu.au>

## Syllabus

Material [web page](https://tidytranscriptomics-workshops.github.io/bioc2022_tidytranscriptomics/articles/tidytranscriptomics_case_study.html).

More details on the workshop are below.

## Workshop package installation

For the bioc2022 workshop, an RStudio in the cloud will be provided with everything installed, all that participants will need is a web browser.

If you want to install the packages and material post-workshop, the instructions are below. The workshop is designed for R `4.1` and Bioconductor 3.14.

```
#install.packages('remotes')

# Need to set this to prevent installation erroring due to even tiny warnings, similar to here: https://github.com/r-lib/remotes/issues/403#issuecomment-748181946
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

# Install same versions used in the workshop
remotes::install_github("stemangiola/tidySingleCellExperiment@v1.7.2")

# Install workshop package

remotes::install_github("tidytranscriptomics-workshops/bioc2022_tidytranscriptomics", build_vignettes = TRUE)

# To view vignettes
library(bioc2022tidytranscriptomics)
browseVignettes("bioc2022tidytranscriptomics")
```

To run the code, you could then copy and paste the code from the workshop vignette or [R markdown file](https://raw.githubusercontent.com/tidytranscriptomics-workshops/bioc2022_tidytranscriptomics/master/vignettes/tidytranscriptomics.Rmd) into a new R Markdown file on your computer.

## Workshop Description

This tutorial will present how to perform analysis of single-cell RNA sequencing data following the tidy data paradigm. The tidy data paradigm provides a standard way to organise data values within a dataset, where each variable is a column, each observation is a row, and data is manipulated using an easy-to-understand vocabulary. Most importantly, the data structure remains consistent across manipulation and analysis functions.

This can be achieved with the integration of packages present in the R CRAN and Bioconductor ecosystem, including [tidySingleCellExperiment](https://stemangiola.github.io/tidySingleCellExperiment/) and [tidyverse](https://www.tidyverse.org/). These packages are part of the tidytranscriptomics suite that introduces a tidy approach to RNA sequencing data representation and analysis. For more information see the [tidy transcriptomics blog](https://stemangiola.github.io/tidytranscriptomics/).

### Pre-requisites

* Basic familiarity with single cell transcriptomic analyses
* Basic familiarity with tidyverse


### Workshop Participation

The workshop format is a 2 hour session consisting of hands-on demos, exercises and Q&A.

### _R_ / _Bioconductor_ packages used

* tidySingleCellExperiment
* org.Hs.eg.db
* ggrepel
* GGally
* plotly



### Workshop goals and objectives

In exploring and analysing RNA sequencing data, there are a number of key concepts, such as filtering, scaling, dimensionality reduction, hypothesis testing, clustering and visualisation, that need to be understood. These concepts can be intuitively explained to new users, however, (i) the use of a heterogeneous vocabulary and jargon by methodologies/algorithms/packages, (ii) the complexity of data wrangling, and (iii) the coding burden, impede effective learning of the statistics and biology underlying an informed RNA sequencing analysis.

The tidytranscriptomics approach to RNA sequencing data analysis abstracts out the coding-related complexity and provides tools that use an intuitive and jargon-free vocabulary, enabling focus on the statistical and biological challenges.

#### Learning goals

* To approach data representation and analysis though a tidy data paradigm, integrating tidyverse, tidySingleCellExperiment and tidyHeatmap.

#### What you will learn

* Basic tidy operations possible with tidySingleCellExperiment
* The differences between SingleCellExperiment representation and tidy representation
* How to interface SingleCellExperiment with tidy manipulation and visualisation
* A real-world case study that will showcase the power of tidy single-cell methods compared with base/ad-hoc methods

#### What you will not learn

* The molecular technology of single-cell sequencing
* The fundamentals of single-cell data analysis
* The fundamentals of tidy data analysis
