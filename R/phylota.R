require('taxize')
require('ape')
require('plyr')
require('CHNOSZ')
require('data.table')
require('igraph')
source('ncbi-tax.R')
source('ncbi-gen.R')
source('blast.R')
source('cluster.R')

ROOT.TAXON <- 8461

## Phylota parameters
LENGTH.CUTOFF <<- 25000

## Query and retreive sequences
cat("Querying GIs for taxon ID ", ROOT.TAXON, "\n")
gis <- gis.for.taxid(ROOT.TAXON)
cat("Found ", length(gis), " hits for query\n")

cat("Retrieving sequences\n")

seqs <- seqs.for.gis( gis )
lengths <- sapply(seqs, nchar)

## omit sequences that are too long
seqs <- seqs[which(lengths < LENGTH.CUTOFF)]
cat("Using ", length(seqs), " sequences in BLAST search\n")

## make database and BLAST
dbfile <- make.blast.db(seqs)
blast.results <- blast.all.vs.all(dbfile)
cat("Number of BLAST results ", nrow(blast.results), "\n")
filtered.blast.results <- filter.blast.results(blast.results)

## get sequence clusters
clusters <- cluster.seqs(filtered.blast.results)

##br2 = blast.results.chelidae[which(blast.results.chelidae$subject.id %in% gis),]
##br.filtered = br2[which(br2$query.id %in% gis),]
