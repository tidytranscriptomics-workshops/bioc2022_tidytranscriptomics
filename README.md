<!-- badges: start -->
[![DOI](https://zenodo.org/badge/379767139.svg)](https://zenodo.org/badge/latestdoi/379767139)
[![.github/workflows/basic_checks.yaml](https://github.com/tidytranscriptomics-workshops/bioc2022_tidytranscriptomics/workflows/.github/workflows/basic_checks.yaml/badge.svg)](https://github.com/tidytranscriptomics-workshops/bioc2022_tidytranscriptomics/actions) 	
<!-- badges: end -->

# Introduction to Tidy Transcriptomics
<p float="left">
<img style="height:100px;" alt="BioC2022" src="https://bioc2022.bioconductor.org/img/carousel/BioC2022.png"/>
<img style="height:100px;" alt="tidybulk" src="https://github.com/Bioconductor/BiocStickers/blob/master/tidybulk/tidybulk.png?raw=true"/>
</p>

### Pre-requisites

* Basic familiarity with single-cell transcriptomic analyses
* Basic familiarity with tidyverse

## Workshop Description

This tutorial will present how to perform analysis of single-cell RNA sequencing data following the tidy data paradigm. The tidy data paradigm provides a standard way to organise data values within a dataset, where each variable is a column, each observation is a row, and data is manipulated using an easy-to-understand vocabulary. Most importantly, the data structure remains consistent across manipulation and analysis functions.

This can be achieved with the integration of packages present in the R CRAN and Bioconductor ecosystem, including [tidySingleCellExperiment](https://stemangiola.github.io/tidySingleCellExperiment/) and [tidyverse](https://www.tidyverse.org/). These packages are part of the tidytranscriptomics suite that introduces a tidy approach to RNA sequencing data representation and analysis. For more information see the [tidy transcriptomics blog](https://stemangiola.github.io/tidytranscriptomics/).

## Workshop goals

* To approach single-cell data representation and analysis though a tidy data paradigm, integrating tidyverse with tidySingleCellExperiment.
* Compare SingleCellExperiment and tidy representation  
* Apply tidy functions to SingleCellExperiment objects  
* Reproduce a real-world case study that showcases the power of tidy single-cell methods

### What you will learn

* Basic tidy operations possible with tidySingleCellExperiment
* The differences between SingleCellExperiment representation and tidy representation
* How to interface SingleCellExperiment with tidy manipulation and visualisation
* A real-world case study that will showcase the power of tidy single-cell methods compared with base/ad-hoc methods

### What you will not learn

* The molecular technology of single-cell sequencing
* The fundamentals of single-cell data analysis
* The fundamentals of tidy data analysis

### Workshop Participation

The workshop format is a 1.5 hour session consisting of hands-on demos, exercises and Q&A.

## Syllabus

Material [web page](https://tidytranscriptomics-workshops.github.io/bioc2022_tidytranscriptomics/articles/tidytranscriptomics_case_study.html).

More details on the workshop are below.

## Workshop package installation 

For the BioC2022 workshop, an RStudio in the cloud will be provided with everything installed, all that participants will need is a web browser. 

If you want to install the packages and material post-workshop, they can be installed using one of the two ways below. The workshop is designed for R `4.2` and Bioconductor 3.16.

### Via Docker image

If you're familiar with [Docker](https://docs.docker.com/get-docker/) you could use the Docker image which has all the software pre-configured to the correct versions.

```
docker run -e PASSWORD=abc -p 8787:8787 ghcr.io/tidytranscriptomics-workshops/bioc2022_tidytranscriptomics
```

Once running, navigate to <http://localhost:8787/> and then login with
`Username:rstudio` and `Password:abc`.

You should see the Rmarkdown file with all the workshop code which you can run.

### Via GitHub

Alternatively, you could install the workshop using the commands below in R `4.2`.

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

## Instructor names and contact information

* Maria Doyle <Maria.Doyle at petermac.org>  
* Stefano Mangiola <mangiola.s at wehi.edu.au>

