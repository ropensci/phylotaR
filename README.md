# Automated Retrieval of Orthologous DNA Sequences from GenBank <img src="https://raw.githubusercontent.com/ropensci/phylotaR/master/logo.png" height="300" align="right"/>
[![Build Status](https://travis-ci.org/ropensci/phylotaR.svg?branch=master)](https://travis-ci.org/ropensci/phylotaR) [![Coverage Status](https://coveralls.io/repos/github/ropensci/phylotaR/badge.svg?branch=master)](https://coveralls.io/github/ropensci/phylotaR?branch=master) [![](https://badges.ropensci.org/187_status.svg)](https://github.com/ropensci/onboarding/issues/187)

R implementation of the [PhyLoTa sequence cluster pipeline](http://phylota.net/). For more information see the accompanying website. Tested and demonstrated on Unix and Windows. **Find out more by visiting the [phylotaR website](https://ropensci.github.io/phylotaR/).**

## Install
Currently only the devlopment package is available:

```r
devtools::install_github(repo='ropensci/phylotaR',
                         build_vignettes=TRUE)
```

**Full functionality depends on a local copy of BLAST+ (>= 2.0.0)**. For details on downloading and compiling BLAST+ on your machine please visit the [NCBI website](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

## Pipeline

`phylotaR` runs the PhyLoTa pipeline in four automated stages: identify and retrieve taxonomic information on all descendent nodes of the taxonomic group of interest (`taxise`), download sequence data for every identified node (`download`), identify orthologous clusters using BLAST (`cluster`), and identify sister clusters for sets of clusters identified in the previous stage (`cluster^2`) After these stages are complete, `phylotaR` provides tools for exploring, identifying and exporting suitable clusters for subsequent analysis.

![phylotaR pipeline](https://raw.githubusercontent.com/ropensci/phylotaR/master/other/stages.png)

For more information on the pipeline and how it works see the publication, [phylotaR: An Automated Pipeline for Retrieving Orthologous DNA Sequences from GenBank in R](https://doi.org/10.3390/life8020020).

## Running

At a minimum all a user need do is provide the taxonomic ID of their chosen taxonomic group of interest. For example, if you were interested in primates, you can visit the [NCBI taxonomy home page](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/) and search primates to look up their ID. After identifying the ID, the `phylotaR` pipeline can be run with the following script.

```r
library(phylotaR)
wd <- '[FILEPATH TO WORKING DIRECTORY]'
ncbi_dr <- '[FILEPATH TO COMPILED BLAST+ TOOLS]'
txid <- 9443  # primates ID
setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr)
run(wd = wd)
```

The pipeline can be stopped and restarted at any point without loss of data. For more details on this script, how to change parameters, check the log and details of the pipeline, please check out the package vignette.

```r
library(phylotaR)
vignette("phylotaR")
```

## Timings

How long does it take for a phylotaR pipeline to complete? Below is a table listing the runtimes in minutes for different demonstration, taxonomic groups. 

Taxon|N. taxa|N. sequences|N. clusters|Taxise (mins.)|Download (mins.)|Cluster (mins.)|Cluster2 (mins.)|Total (mins.)|
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
Anisoptera|1175|11432|796|1.6|23|48|0.017|72|
Acipenseridae|51|2407|333|0.1|6.9|6.4|0.017|13|
Tinamiformes|25|251|98|0.067|2.4|0.18|0.017|2.7|
Aotus|13|1499|193|0.067|3.2|0.6|0|3.9|
Bromeliaceae|1171|9833|724|1.2|28|37|0.033|66|
Cycadidae|353|8331|540|0.32|19|18|0.033|37|
Eutardigrada|261|960|211|0.3|11|1.8|0.05|14|
Kazachstania|40|623|101|0.1|20|3|0.05|23|
Platyrrhini|212|12731|3112|0.35|51|6.9|1.2|60|

To run these same demonstrations see [demos/demo_run.R](https://github.com/ropensci/phylotaR/blob/master/demos/demo_run.R).

## License

MIT

## Version

Version 1.

## Authors

Dom Bennett (maintainer, R package dev), Hannes Hettling (workhouse code dev), Rutger Vos, Alexander Zizka and Alexandre Antonelli

## Reference

Bennett, D., Hettling, H., Silvestro, D., Zizka, A., Bacon, C., Faurby, S., … Antonelli, A. (2018). phylotaR: An Automated Pipeline for Retrieving Orthologous DNA Sequences from GenBank in R. *Life*, **8**(2), 20. [DOI:10.3390/life8020020](https://doi.org/10.3390/life8020020)

Sanderson, M. J., Boss, D., Chen, D., Cranston, K. A., & Wehe, A. (2008). The PhyLoTA Browser: Processing GenBank for molecular phylogenetics research. *Systematic Biology*, **57**(3), 335–346. [DOI:10.1080/10635150802158688](https://doi.org/10.1080/10635150802158688)

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
