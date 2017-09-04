## script to make the phylota tables

source('nodes.R')
source('clusters.R')
source('ci_gi.R')

options(error=recover)

## parse command args for root taxon
#args <- commandArgs(trailingOnly=TRUE)
##root.taxon <- args[1]

dbname <- 'phylota.sqlite'
taxdir <- '/Users/hettling/taxdump'##'/Users/hettling/ftp.ncbi.nih.gov/pub/taxonomy'

## make accession2taxid table
#accession2taxid.create(dbname, taxdir, overwrite=T)

library(foreach)
library(doMC)
registerDoMC(3)

ncbi.names <<- getnames(taxdir)

## root <- "Mammalia"
root <- "Carnivora"
root.taxa <- ncbi.names$id[which(ncbi.names$name == root)]
cat("Taxid(s) to analyse : ", root.taxa)

##remote.nodes.create(root.taxa=root.taxa, file.name='nodes-panthera.tsv')
nodes.create(taxdir, root.taxa=root.taxa, file.name='nodes-mammalia.tsv', model.threshold=50000)

##.db(dbname)
##cl <- clusters.ci_gi.seqs.create(root.taxa, 'nodes-panthera.tsv',
##                                 files=list(clusters='dbfiles-panthera-clusters.tsv',
##                                     ci_gi='dbfiles-panthera-ci_gi.tsv',
##                                    seqs='dbfiles-panthera-seqs.tsv'), max.seqs=130000)
