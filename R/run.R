## script to make the phylota tables

source('accession2taxid.R')
source('nodes.R')
source('clusters.R')

## parse command args for root taxon
#args <- commandArgs(trailingOnly=TRUE)
##root.taxon <- args[1]

dbname <- 'phylota-02'
taxdir <- '/home/hettling/ftp.ncbi.nih.gov/pub/taxonomy'

## make accession2taxid table
#accession2taxid.create(dbname, taxdir, overwrite=T)

## make nodes table
root.taxon= 1445964 #3683
#nodes.create(dbname, taxdir, root.taxa=root.taxon, overwrite=F, append=T)
##cl <- clusters.create(root.taxon, dbname)
cl <- .make.clusters(root.taxon)
