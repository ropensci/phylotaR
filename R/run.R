## script to make the phylota tables

source('accession2taxid.R')
source('nodes.R')

## parse command args for root taxon
args <- commandArgs(trailingOnly=TRUE)
root.taxon <- args[1]

dbname <- 'phylota-01'
taxdir <- '~/ftp.ncbi.nih.gov/pub/taxonomy'

## make accession2taxid table
accession2taxid.create(dbname, taxdir, overwrite=F)

## make nodes table
nodes.create(dbname, taxdir, root.taxa=root, overwrite=T)
