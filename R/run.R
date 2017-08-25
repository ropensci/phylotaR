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
registerDoMC(1)

## make nodes table
##root.taxon <- 167485##4055##167489 ##3297 ##1445963##1437180 ##3313 ##1445964  ##1437180 #1445964 ## ##1437180##1445964##3394#3297 #3683 #3684#3683#

##root.names = scan('woody-families.txt', what=character())
##ncbi.names <<- getnames(taxdir)

##root.taxa = unlist(sapply(root.names, function(x)ncbi.names$id[which(ncbi.names$name==x)]))
root.taxa <- 9688 ##9682
##remote.nodes.create(root.taxa=root.taxa, file.name='nodes-panthera.tsv')
##nodes.create(taxdir, root.taxa=root.taxa, file.name='nodes-panthera.tsv', model.threshold=130000)

##.db(dbname)
##cl <- clusters.ci_gi.seqs.create(root.taxa, 'nodes-panthera.tsv',
                                 files=list(clusters='dbfiles-panthera-clusters.tsv',
                                     ci_gi='dbfiles-panthera-ci_gi.tsv',
                                    seqs='dbfiles-panthera-seqs.tsv'), max.seqs=130000)


##dfs <- clusters.ci_gi.seqs.create(root.taxon, 'dbfiles')
##clusters <- .make.clusters(root.taxon)
##cldf <- .make.cluster.entries(clusters)
##cigidf <- ci_gi.create(clusters)
