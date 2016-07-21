require('taxize')
require('rentrez')
require('ape')
require('plyr')
require('igraph')
source('ncbi-tax.R')

ROOT.TAXON <- 124867

## Phylota parameters
LENGTH.CUTOFF <- 25000
MIN.COVERAGE <- 0.51

#workdir <- paste0('taxid-', ROOT.TAXON)
#dir.create(workdir)
#setwd(workdir)

##sink(paste0(ROOT.TAXON, ".txt"))

expand.taxa <- function(taxa, rank='species') {
    db <- 'ncbi'

    queue <- taxa
    result <- data.frame(ncbi.id=integer(0), name=character(0), rank=character(0))
    while ( length(queue) > 0 ) {
        taxon <- queue[1]
        ##cat("Taxon : ", taxon, "\n")
        
        ## remove first element from queue
        queue <- tail(queue, length(queue)-1)        

        ch <- children(taxon, db=db)
        table <- ch[[1]]       
        ##        result <- rbind(result, subset(table, childtaxa_rank==rank))
        result <- rbind(result, table)
        queue <- c(queue, table$childtaxa_name)        
    }
    return (result)
}

## Query and retreive sequences
cat("Querying GIs for taxon ID ", ROOT.TAXON, "\n")
ids <- acc.genbank.seqids.for.taxid(ROOT.TAXON)
gis <- as.vector(sapply(ids, '[[', "gi"))

cat("Found ", length(gis), " hits for query\n")

cat("Retrieving sequences\n")
print(system.time({all.recs <- lapply(gis, function(gi)entrez_fetch(db="nuccore",id=gi, rettype='fasta'))}))

all.recs.cleaned <- gsub("^>.*?\\n|\\n", "", all.recs)
seqs <- all.recs.cleaned
names(seqs) <- gis
lengths <- sapply(seqs, nchar)

## omit sequences that are too long
seqs <- seqs[which(lengths < LENGTH.CUTOFF)]
cat("Using ", length(seqs), " sequences in BLAST search\n")

## write all sequences to BLAST database file
dbfile <- "blastdb.fa"
file.create(dbfile)
for (gi in names(seqs)) {
    write(paste0("> ", gi, "\n", seqs[[as.character(gi)]],"\n"), file=dbfile, append=T)
}

## Now make database and blast all vs all
cat("Making database and performing all-vs-all BLAST\n")
print(system.time(
{
    system('makeblastdb -in blastdb.fa -dbtype nucl')
    system('blastn -query blastdb.fa -db blastdb.fa -outfmt 6 -out blastout.txt')
}))

blast.results <- read.table('blastout.txt')
colnames(blast.results) <- c('query.id', 'subject.id', 'identity', 'alignment.length',
                             'mismatches', 'gap.opens', 'q.start', 'q.end', 's.start',
                             's.end', 'evalue', 'bit.score')


cat("Number of BLAST results ", nrow(blast.results), "\n")

## collapse HSPs such that we end up with unique query-subject pairs
cat("Filtering BLAST results\n")
print(system.time(
{
    result.subset <- ddply(blast.results, c("query.id", "subject.id"), function(x)colSums(x['alignment.length']))
    result.subset['query.length'] <- sapply(result.subset$query.id, function(id){nchar(seqs[[as.character(id)]])})
    result.subset['subject.length'] <- sapply(result.subset$subject.id, function(id){nchar(seqs[[as.character(id)]])})
    coverages <- apply(result.subset, 1, function(x)x['alignment.length'] / max(x['query.length'], x['subject.length']))
    num.discarded.hits <- sum(coverages < MIN.COVERAGE)

    ## keep only the gis for which there is enough overlapping hits
    result.subset <- result.subset[which(coverages >= MIN.COVERAGE),]
}))

cat("Discarding ", num.discarded.hits, " BLAST hits due to insufficient coverage \n")

## make graph object and produce single-linkage clusters
cat("Clustering results\n")
print(system.time(
{
    g <- graph.data.frame(result.subset[,c("query.id", "subject.id")], directed=F)
    clusters <- clusters(g)
}))

cat("Number of sequence clusters : ", length(clusters$member), "\n")

## filter for phylogenetically informative clusters
informative.clusters <- clusters$membership[which(clusters$membership %in% which (clusters$csize > 2))]
cat("Number of informative sequence clusters : ", length(informative.clusters), "\n")
cluster.seqs <- names(informative.clusters)
cat("Number of sequences in cluster : ", length(cluster.seqs), "\n")


## get species etc. for root taxon
taxa.table <- expand.taxa(ROOT.TAXON)
nof.species <- sum(taxa.table$childtaxa_rank == 'species')
cat("Number of species for root taxon : ", nof.species, "\n")

save.image(paste0(ROOT.TAXON, ".Rdata"))
##sink()
