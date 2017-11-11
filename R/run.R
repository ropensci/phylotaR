source("blast.R")
source("ci_gi.R")
source("cl.R")
source("ncbi-remote.R")
source("nodes.R")

library(foreach)
library(doMC)
library(igraph)
library(CHNOSZ)
library(rentrez)
library(data.table)

set.seed(111)

options(error=recover)

## Adjustable parameters:
## Maximum number of sequences per species
MODEL.THRESHOLD <<- 3000
## Maximum number of sequences to blast in a single run; if taxon has more subtree sequences
## than that, its children will get clustered
MAX.BLAST.SEQS <<- 10000
## Maximum characters in one sequence
MAX.SEQUENCE.LENGTH <<- 25000
## directory for sequence cache; will be created if does not exist
SEQS.CACHE.DIR <<- "./sequences/"
## Download file ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz,
## unzip and specify directory where file 'nodes.dmp' is located
taxdir <- '/Users/hettling/taxdump'
## Number of processing units
CORES <<- 4

registerDoMC(CORES)

## Do analysis for Bromeliaceae family
taxid <- 4613
## nodes.create(taxid, taxdir=taxdir, file.name='dbfiles-bromeliaceae-nodes.tsv')

clusters.ci_gi.seqs.create(15123, 'dbfiles-bromeliaceae-nodes.tsv',
                           files=list(clusters='dbfiles-bromeliaceae-clusters.tsv',
                                      ci_gi='dbfiles-bromeliaceae-ci_gi.tsv',
                                      seqs='dbfiles-bromeliaceae-seqs.tsv'))
