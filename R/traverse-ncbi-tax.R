## INPUT: Root taxon ID
## This script takes a root taxon id as an argument and
## creates a table with NCBI nodes for which we will later do clustering
## with all-vs-all blast. These nodes are the ones with manageable
## number of sequences (< 20000). If a node is found, it is written to the table,
## but not it's descendants because we don't need to blast the seqs for the descendants
## again.
## OUTPUT: table (<root taxon if>.tsv), containing the NCBI nodes and the information:
## taxid, name, nof.seqs, rank, cluster type (subtree or node)

source('ncbi-tax.R')

## parse command args for root taxon
args <- commandArgs(trailingOnly=TRUE)
root.taxon <- args[1]

cat("Obtaining node table for taxon ", root.taxon, "\n")

## get table
tab <- get.cluster.nodes(root.taxon)
cat("Retreived table with ", nrow(tab), " rows\n")

## write to file
file.name <- paste0(root.taxon, '.tsv')
write.table(tab, file=file.name, sep="\t", quote=F, row.names=F)
cat("Table written to ", file.name, "\n")

