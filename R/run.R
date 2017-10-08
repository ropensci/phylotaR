## script to make the phylota tables
source('nodes.R')
source('clusters.R')
source('ci_gi.R')
source('cl.R')
library(foreach)
library(doMC)

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

## Do analysis for Felidae family
taxid <- 9681
nodes.create(taxid, taxdir=taxdir, file.name='dbfiles-felidae-nodes.tsv')

clusters.ci_gi.seqs.create(338152, 'dbfiles-felidae-nodes.tsv', files=list(clusters='dbfiles-felidae-clusters.tsv',
                                                            ci_gi='dbfiles-felidae-ci_gi.tsv',
                                                            seqs='dbfiles-felidae-seqs.tsv'))
