## script to populate the phylota 'clusters' table
## Dependencies: tables 'accession2taxid', 'nodes'

library("RSQLite")
library('igraph')
source('blast.R')
source('db.R')



## Schema in phylota database:
## CREATE TABLE "clusters_194" (
##  "ti_root" int(10)  DEFAULT NULL,
##  "ci" int(10)  DEFAULT NULL,
##  "cl_type" text  DEFAULT NULL,
##  "n_gi" int(10)  DEFAULT NULL,
##  "n_ti" int(10)  DEFAULT NULL,
##  "PI" tinyint(1) DEFAULT NULL,
##  "MinLength" int(10)  DEFAULT NULL,
##  "MaxLength" int(10)  DEFAULT NULL,
##  "MaxAlignDens" float DEFAULT NULL,
##  "ci_anc" int(10)  DEFAULT NULL,
##  "seed_gi" bigint(20)  DEFAULT NULL,
##  "Q" float DEFAULT NULL,
##  "TC" float DEFAULT NULL,
##  "clustalw_tree" longtext,
##  "muscle_tree" longtext,
##  "strict_tree" longtext,
##  "clustalw_res" float DEFAULT NULL,
##  "muscle_res" float DEFAULT NULL,
##  "strict_res" float DEFAULT NULL,
##  "ortho" tinyint(4) DEFAULT NULL,
##  "n_gen" int(10)  DEFAULT NULL,
##  "n_child" int(10)  DEFAULT NULL,
##  "muscle_tree_mid" longtext,
##  "muscle_tree_mid_chron" longtext,
##  "muscle_tree_aln_length" int(10)  DEFAULT NULL,
##  "PD_data" float DEFAULT NULL,
##  "PD_time" float DEFAULT NULL,
##  "brlen_time_max" float DEFAULT NULL,
##  "mp_score" int(10)  DEFAULT NULL,
##  "uid" int(11) NOT NULL ,
##  "uid_anc" int(10)  DEFAULT NULL,
##  "root_cluster" int(10)  DEFAULT NULL,
##  "lca" int(10)  DEFAULT NULL,
##  PRIMARY KEY ("uid")
##);

clusters.create <- function(root.taxon, max.seqs=25000) {

    all.clusters <- list()
    queue <- root.taxon
    while(length(queue) > 0) {
        currentid <- head(queue, 1)
        queue <- tail(queue, length(queue)-1)

        if (.num.seqs.for.taxid(currentid) > max.seqs) {
            queue <- c(queue, .children(currentid))
        }
        else {
            cl <- .make.clusters(currentid)
            all.clusters <- rbind(all.clusters, cl)
        }
    }
    return (all.clusters)
}

.make.clusters <- function(taxon) {
    cat("Retrieving sequences\n")
    seqs <- .seqs.for.taxid( taxon )
    ##gis <- names(seqs)
    dbfile <- make.blast.db(seqs)
    blast.results <- blast.all.vs.all(dbfile)
    cat("Number of BLAST results ", nrow(blast.results), "\n")
    filtered.blast.results <- filter.blast.results(blast.results)
    ## get sequence clusters
    clusters <- cluster.blast.results(filtered.blast.results)

}

## input: List of clusters
.make.cluster.entries <- function(clusters, taxid, parent, seqs) {
    ## take the column names in the database as names here
    result <- data.frame()
    for (i in 1:length(clusters)) {
        cl <- clusters[[i]]
        entry <- vector()
        entry['ti_root'] <- taxid
        entry['ci'] <- i-1
        entry['cl_type'] <- ifelse(length(.children(taxid)) > 0, 'subtree', 'node')
        entry['n_gi'] <- length(cl)
        ## all taxon ids for gis in cluster
        tis <- unique(sapply(cl, .taxid.for.gi))
        entry['n_ti'] <- length(tis)
        ## get sequences and determine their lengths
        cl.seqs <- lapply(cl, function(gi)seqs[[as.character(gi)]])
        l <- sapply(cl.seqs, nchar)
        entry['MinLength'] <- min(l)
        entry['MaxLength'] <- max(l)
        entry['ci_anc'] <- parent
        entry['seed_gi'] <- names(clusters)[i]
        ## get number of genera
        entry['n_gen'] <- length(unique(sapply(tis, .genus.for.taxid)))

    }

}

.make.child.clusters <- function(cluster, taxid) {

}

cluster.blast.results <- function(blast.results, informative=T) {
    g <- graph.data.frame(filtered.blast.results[,c("query.id", "subject.id")], directed=F)
    clusters <- clusters(g)

    ## filter for phylogenetically informative clusters
    if (informative) {
        clusters <- clusters$membership[which(clusters$membership %in% which (clusters$csize > 2))]
    }
    ## we will return a list, one entry with sequence IDs for each cluster
    cluster.list <- lapply(unique(clusters), function(x)sort(names(clusters)[which(clusters==x)]))

    ## Get the seed gi, we will chose it to be the sequence in the cluster that has
    ## the most hits with the ther members in the cluster; i.e. the most connected
    ## node in the graph
    degrees <- degree(g)
    seed.gis <- sapply(cluster.list, function(c){
        idx <- order(degrees[c], decreasing=T)[1] ## index of most connected component
        c[idx]
    })
    ## set seed gis as names to cluster list
    names(cluster.list) <- as.character(seed.gis)
    ##names(cluster.list) <- paste('cluster', 1:length(cluster.list), sep='')
    return (cluster.list)
}

.children <- function(taxid) {
    db <- .db()
    query <- paste('select ti from nodes where ti_anc =', taxid)
    l <- dbGetQuery(db, query)
    return(l[[1]])
}

.gis.for.taxid <- function(taxid) {
    search <- entrez_search(db='nucleotide', term=paste0('txid', taxid, '[Organism:exp]', '1:25000[SLEN]'), use_history=T, retmax=1e9-1)
    return(search$ids)
}

.num.seqs.for.taxid <- function(taxid) {
    search <- entrez_search(db='nucleotide', term=paste0('txid', taxid, '[Organism:exp]', '1:25000[SLEN]'), use_history=T, retmax=1)
    return(search$count)
}

.genera.for.taxids <- function(taxids) {
    db <- .db()
    query <- paste('select ti_genus from nodes where ti in (', paste(taxids, collapse=','), ')')
    l <- dbGetQuery(db, query)
    return(l[[1]])
}

.seqs.for.taxid <- function(taxid) {
    search <- entrez_search(db='nucleotide', term=paste0('txid', taxid, '[Organism:exp]', '1:25000[SLEN]'), use_history=T, retmax=1e9-1)
    gis <- search$ids;
    cat("Going to retrieve ", length(gis), " sequences\n")
    allseqs <- numeric()
    for (seq.idx in seq(0, length(gis), 10000)) {
        cat("Retreiving seqs ", seq.idx+1, " to ", ifelse(seq.idx+10000<length(gis), seq.idx+10000, length(gis)), "\n");
        res <- entrez_fetch(db="nuccore", rettype="fasta", web_history=search$web_history, retmax=10000, retstart=seq.idx)
        ## split fasta string on '>'
        seqs <- unlist(strsplit(res, "(?<=[^>])(?=>)", perl=T))
        ## clear defline
        seqs <- gsub("^>.*?\\n", "", seqs)
        ## clear newlines
        seqs <- gsub("\n", "", seqs, fixed=TRUE)
        allseqs <- c(allseqs, seqs)
    }
    names(allseqs) <- gis
    return(allseqs)
}

#.taxid.for.gi <- function(gi) {
#    cat("Retreiving taxid for gi ", gi, "\n")
#    search <- entrez_search(db='nucleotide', term=gi, use_history=T)
#    res <- entrez_fetch(db="nuccore", rettype="xml", web_history=search$web_history, complexity=4)
#    ti <- as.numeric(gsub(".*taxon.*?id\\s([0-9]+).*", '\\1', res))
#    return(ti)
#}

.taxid.for.gi <- function(gi) {
    cat("Retreiving taxid for gi ", gi, "\n")
    db <- .db()
    query <- paste('select taxid from accession2taxid where gi=', gi)
    l <- dbGetQuery(db, query)
    return(l[[1]])
}

.genus.for.taxid <- function(taxid) {
    cat("Retreiving genus for taxid ", taxid, "\n")
    db <- .db()
    query <- paste('select ti_genus from nodes where ti=', taxid)
    l <- dbGetQuery(db, query)
    return(l[[1]])
}
