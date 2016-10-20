## script to make the phylota tables

source('nodes.R')
source('clusters.R')
source('ci_gi.R')

options(error=recover)

## parse command args for root taxon
#args <- commandArgs(trailingOnly=TRUE)
##root.taxon <- args[1]

dbname <- 'phylota-02'
taxdir <- '/home/hettling/ftp.ncbi.nih.gov/pub/taxonomy'

## make accession2taxid table
#accession2taxid.create(dbname, taxdir, overwrite=T)

library(foreach)
library(doMC)
registerDoMC(4)

## make nodes table
root.taxon <- 3313 ##1445964  ##1437180 #1445964 ## ##1437180##1445964##3394#3297 #3683 #3684#3683# 
##nodes.create(dbname, taxdir, root.taxa=root.taxon, overwrite=T, append=F)
##cl <- clusters.create(root.taxon, dbname)

.db(dbname)

dfs <- clusters.ci_gi.seqs.create(root.taxon, 'dbfiles')
##clusters <- .make.clusters(root.taxon)
##cldf <- .make.cluster.entries(clusters)
##cigidf <- ci_gi.create(clusters)
