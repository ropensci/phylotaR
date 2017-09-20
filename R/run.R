## script to make the phylota tables

source('nodes.R')
source('clusters.R')
source('ci_gi.R')

options(error=recover)

## parse command args for root taxon
#args <- commandArgs(trailingOnly=TRUE)
##root.taxon <- args[1]

dbname <- 'phylota.sqlite'
taxdir <- '/home/hettling/taxdump'##'/Users/hettling/ftp.ncbi.nih.gov/pub/taxonomy'

## make accession2taxid table
#accession2taxid.create(dbname, taxdir, overwrite=T)

library(foreach)
library(doMC)
registerDoMC(3)

##ncbi.names <<- getnames(taxdir)

## root <- "Mammalia"
##root <- "Felidae"
##root.taxa <- ncbi.names$id[which(ncbi.names$name == root)]
##cat("Taxid(s) to analyse : ", root.taxa)

##remote.nodes.create(root.taxa=root.taxa, file.name='nodes-panthera.tsv')
#root.taxa <- 9681
#nodes.create(taxdir, root.taxa=root.taxa, file.name='nodes-felidae.tsv', model.threshold=10000)

source('cl.R')
nodes <- read.table('nodes-felidae.tsv', header=T)
##cl <- cluster(9681, nodes)
cl <- cluster(37028, nodes)



##.db(dbname)
##cl <- clusters.ci_gi.seqs.create(root.taxa, 'nodes-felinae.tsv',
##                                 files=list(clusters='dbfiles-felinae-clusters.tsv',
##                                     ci_gi='dbfiles-felinae-ci_gi.tsv',
##                                    seqs='dbfiles-felinae-seqs.tsv'), max.seqs=10000)
