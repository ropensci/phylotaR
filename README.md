# phylotaR
[![Build Status](https://travis-ci.org/AntonelliLab/phylotaR.svg?branch=master)](https://travis-ci.org/AntonelliLab/phylotaR) [![Coverage Status](https://coveralls.io/repos/github/AntonelliLab/phylotaR/badge.svg?branch=master)](https://coveralls.io/github/AntonelliLab/phylotaR?branch=master) [![](https://badges.ropensci.org/187_status.svg)](https://github.com/ropensci/onboarding/issues/187)

R implementation of the PhyLoTa sequence cluster pipeline (http://phylota.net/).

## Install
Currently only the devlopment package is available:

```r
devtools::install_github(repo='AntonelliLab/phylotaR', build_vignettes=TRUE)
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

Taxon|Taxise|Download|Cluster|Cluster2|Total|
|:--|--:|--:|--:|--:|--:|
Anisoptera|0.8|12|21|0.017|33|
Acipenseridae|0.067|NA|NA|NA|NA|
Tinamiformes|0.083|1.3|0.22|0|1.6|
Aotus|0.05|NA|NA|NA|NA|
Bromeliaceae|0.67|15|14|0.017|30|
Cycadidae|0.18|11|5.3|0.017|17|
Eutardigrada|0.28|5|0.82|0|6|
Kazachstania|0.067|9.3|1.1|0.033|10|
Platyrrhini|0.2|30|1.6|0.28|32|

To run these same demonstrations see [´demos/demo_run.R´](demos/demos_run.R). 

## License

MIT

## Version

Development version 0.1

## Authors

Dom Bennett (maintainer, R package dev), Hannes Hettling (workhouse code dev), Rutger Vos, Alexander Zizka and Alexandre Antonelli

## Reference

Bennett, D., Hettling, H., Silvestro, D., Zizka, A., Bacon, C., Faurby, S., … Antonelli, A. (2018). phylotaR: An Automated Pipeline for Retrieving Orthologous DNA Sequences from GenBank in R. Life, 8(2), 20. [DOI:10.3390/life8020020](https://doi.org/10.3390/life8020020)

Sanderson, M. J., Boss, D., Chen, D., Cranston, K. A., & Wehe, A. (2008). The PhyLoTA Browser: Processing GenBank for molecular phylogenetics research. *Systematic Biology*, **57**(3), 335–346. [DOI:10.1080/10635150802158688](https://doi.org/10.1080/10635150802158688)
