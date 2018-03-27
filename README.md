# phylotaR
[![Build Status](https://travis-ci.org/DomBennett/phylotaR.svg?branch=master)](https://travis-ci.org/DomBennett/phylotaR) [![Coverage Status](https://coveralls.io/repos/github/DomBennett/phylotaR/badge.svg?branch=master)](https://coveralls.io/github/DomBennett/phylotaR?branch=master) [![](https://badges.ropensci.org/187_status.svg)](https://github.com/ropensci/onboarding/issues/187)

> Please note: this package is in its development and is not fully tested.

R implementation of the PhyLoTa sequence cluster pipeline (http://phylota.net/).

## Install
Currently only devlopment package is available:

```r
devtools::install_github(repo='DomBennett/phylotaR', build_vignettes=TRUE)
```

**Full functionality depends on a local copy of BLAST+**. For details on downloading and compiling BLAST+ on your machine please visit the [NCBI website](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

## Pipeline

`phylotaR` runs the PhyLoTa pipeline in four automated stages: identify and retrieve taxonomic information on all descendent nodes of the taxonomic group of interest (`taxise`), download sequence data for every identified node (`download`), identify orthologous clusters using BLAST (`cluster`), and identify sister clusters for sets of clusters identified in the previous stage (`cluster^2`) After these stages are complete, `phylotaR` provides tools for exploring, identifying and exporting suitable clusters for subsequent analysis.

![phylotaR pipeline: internet dependencies are indicated with green circles, external software dependcies are indiciated with orange circles.](https://github.com/DomBennett/phylotaR/raw/master/other/stages.png)

## Running

At a minimum all a user need do is provide the taxonomic ID of their chosen taxonomic group of interest. For example, if you were interested in primates, you can visit the [NCBI taxonomy home page](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/) and search primates to look up their ID. After identifying the ID, the `phylotaR` pipeline can be run with the following script.

```r
library(phylotaR)
wd <- '[FILEPATH TO WORKING DIRECTORY]'
ncbi_dr <- '[FILEPATH TO COMPILED BLAST+ TOOLS]'
txid <- 9443  # primates ID
setUp(wd=wd, txid=txid, ncbi_dr=ncbi_dr)
run(wd=wd)
```

The pipeline can be stopped and restarted at any point without loss of data. For more details on this script, how to change parameters, check the log and details of the pipeline, please check out the package vignette.

```r
library(phylotaR)
vignette("phylotaR")
```

## License

GPL-2

## Version

Development version 0.1

## Authors

Dom Bennett (maintainer, R package dev), Hannes Hettling (workhouse code dev), Rutger Vos, Alexander Zizka and Alexandre Antonelli

## Reference

Sanderson, M. J., Boss, D., Chen, D., Cranston, K. A., & Wehe, A. (2008). The PhyLoTA Browser: Processing GenBank for molecular phylogenetics research. *Systematic Biology*, **57**(3), 335–346. [DOI:10.1080/10635150802158688](https://doi.org/10.1080/10635150802158688)
