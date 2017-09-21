

clusters.ci_gi.seqs.create <- function(root.taxon, nodesfile,
                                       files=list(clusters='clusters.tsv', ci_gi='ci_gi.tsv', seqs='seqs.tsv'),
                                       informative=FALSE) {
    ## load nodes table
    nodes <- read.table(nodesfile, header=T)
    ## load clusters that are already calculated, if any
    clusters.file <- files[['clusters']]
    ci_gi.file <- files[['ci_gi']]
    seqs.file <- files[['seqs']]

    ## get all nodes that will be processed:
    ## all nodes for which all children contain < MAX.BLAST.SEQS sequences
    taxa.to.process <- vector()
    queue <- root.taxon
    while(length(queue) > 0) {
        current.taxon <- head(queue, 1)
        queue <- tail(queue, length(queue)-1)

        ## get number of sequences to determine if it is manageable to calculate clusters
        ## TODO: We should get the counts from the nodes table!!!
        num.seqs <- .rec.num.seqs(current.taxon, nodes, max.len=MAX.SEQUENCE.LENGTH, max.seqs.per.spec=MODEL.THRESHOLD)
        cat("Number of sequences for taxon", current.taxon, ": ", num.seqs, "\n")

        ## if sequence count is smaller than MAX.BLAST.SEQS, add it to the nodes to process, otherwise get children
        if (num.seqs == 0) {
            cat("No sequences for taxid", current.taxon, "\n")
        }
        else if (num.seqs <= MAX.BLAST.SEQS) {
            cat("Will process taxon", current.taxon, "\n")
            taxa.to.process <- c(taxa.to.process, current.taxon)
        }
        else {
            cat("Too many seqs to blast for taxid", current.taxon, ", retrieveing children\n")
            queue <- c(queue, .children(current.taxon, nodes))
        }
    }

    ##foreach (i=seq_along(taxa.to.process), .verbose=T) %dopar% {
    for (i in seq_along(taxa.to.process)) {
        taxid <- taxa.to.process[i]
        cat("Processing taxid ", taxid, " # ", i, " / ", length(taxa.to.process), "\n")

        ## Get the sequences
        seqs <- .rec.retreive.seqs(taxid, nodes, max.len=MAX.SEQUENCE.LENGTH, max.seqs.per.spec=MODEL.THRESHOLD)
        seqdf <- .make.seq.entries(seqs)
        clusters <- .make.clusters(taxid, nodes, seqs, informative=informative)
        cldf <- .make.cluster.entries(clusters)
        cigidf <- .make.ci_gi.entries(clusters)

        ## Write all data to file
        cat("Taxid", taxid, ": writing", nrow(cldf), "clusters,", nrow(seqdf), "sequences,", nrow(cigidf), "to file\n")
        write.table(cldf, file=clusters.file, append=file.exists(clusters.file), quote=F, sep="\t", row.names=F, col.names=!file.exists(clusters.file))
        write.table(seqdf, file=seqs.file, append=file.exists(seqs.file), quote=F, sep="\t", row.names=F, col.names=!file.exists(seqs.file))
        write.table(cigidf, file=ci_gi.file, append=file.exists(ci_gi.file), quote=F, sep="\t", row.names=F, col.names=!file.exists(ci_gi.file))
        cat("Finished processing taxid ", taxid, " # ", i, " / ", length(taxa.to.process), "\n")
    }
}


cluster <- function(taxon, nodes, seqs=NULL, blast.results=NULL, informative=FALSE) {

    ## Retrieve sequences, if not already done
    current.seqs <- list()

    ## get species rank
    rank <- as.character(nodes$rank[match(taxon, nodes$ti)])
    cat ("Processing taxid", taxon, "of rank", rank, "\n")

    if (is.null(seqs)) {
        seqs <- .rec.retreive.seqs(taxon, nodes, max.len=25000, max.seqs.per.spec=10000)
        current.seqs <- seqs
    }
    else {
        ## get GIs for taxon

        ## For rank 'species' it is a special case: sequences should not include the 'subtree'
        ## sequences, only the 'direct' sequences. Otherwise the sequences of subspecies would
        ## be included while they shouldn't.
        gis <- .gis.for.taxid(taxon, direct=(rank=='species'))
        current.seqs <- seqs[gis]
    }

    ## Get blast results if not already calculated
    current.blast.results <- list()

    if (is.null(blast.results)) {
        current.blast.results <- .get.blast.results(taxon, seqs)
    }
    else {
        ##  reduce blast results such that include only the gis for the current taxon
        cat("Using BLAST results from parent cluster\n")
        current.blast.results <- blast.results[which(blast.results$subject.id %in% gis),]
        current.blast.results <- current.blast.results[which(current.blast.results$query.id %in% gis),]

    }

    ## Sequence clusters are stored in a list of lists, named by gi
    ## Get top-level clusters
    cat("Clustering BLAST results for taxid", taxon, "\n")
    clusters <- cluster.blast.results(current.blast.results, informative=informative)
    cat("Finished clustering BLAST results for taxid", taxon, "\n")

    ## make dataframe with fields as in PhyLoTa database
    clusters <- .add.cluster.info(clusters, nodes, taxon, seqs)
    cat("Generated ", length(clusters), "clusters\n")

    ## iterate over taxon's children to retrieve clusters
    for (ch in .children(taxon, nodes)) {
        cat("Processing child taxon of ", taxon, " ", ch, "\n")

        ## calculate clusters for child taxon
        child.clusters <- cluster(ch, nodes, seqs, current.blast.results)
        ## two numbers can only be calculated if we have cluster info on multiple taxonomic levels:
        ##  n_child, the number of child clusters, and ci_anc, the parent cluster.
        ##  Since we are dealing with single-likeage clusters, we can identify the parent cluster of a
        ##  cluster if at least one gi of both clusters is the same.
        child.clusters <- lapply(child.clusters, function(cc) {
            ## Get the indices of the parent clusters that contain a gi of the child clusters
            ## If a gi is in two clusters (e.g. parent and grandparent), take the lower one, by taking
            ## the higher index
            idx <- max(which(sapply(clusters, function(c) any(cc$gis %in% c$gis))))
            ## set parent cluster to child cluster
            cc$ci_anc <- clusters[[idx]]$ci
            ## update child count of parent clusters
            clusters[[idx]]$n_child <- clusters[[idx]]$n_child + 1
                    cc
        })
        ## append child clusters to result cluster list
        clusters <- c(clusters, child.clusters)
    }

    return(clusters)
}

.rec.retreive.seqs <- function(taxid, nodes, max.len=25000, max.seqs.per.spec=100000) {
    ## if rank is species or below, retreive the sequences, and retreive the
    ## sequences of children. Make sure to retreive the 'direct link' sequenes
    cat("Attempting to retrieve sequences for taxid", taxid, "\n")
    seqs <- list()
    node.ranks <- c('species', 'subspecies', 'varietas', 'forma')

    ## get subtree counts. If that is smaller than max.seqs.per.spec, we can take that number
    subtree.count <- .num.seqs.for.taxid(taxid, direct=FALSE, max.len=max.len)
    if (subtree.count <= max.seqs.per.spec) {
        cat(subtree.count, "seqs for taxon", taxid, ", less than maximum of ", max.seqs.per.spec, " ")
        cat("sequences. Retreiving sequences for whole subtree\n")
        seqs <- .seqs.for.taxid(taxid, direct=FALSE, max.len=max.len, max.seqs=max.seqs.per.spec)
        return (seqs)
    }

    current.rank <- .rank.for.taxid(taxid)
    cat("Rank : ", current.rank, "\n")
    if (current.rank %in% node.ranks) {
        ## retreive sequences
        current.seqs <- .seqs.for.taxid(taxid, direct=TRUE, max.len=max.len, max.seqs=max.seqs.per.spec)
        seqs <- c(seqs, current.seqs)
    }
    for (ch in .children(taxid, nodes)) {
        seqs <- c(seqs, .rec.retreive.seqs(ch, nodes, max.len, max.seqs.per.spec))
    }
    return(seqs)
}

.rec.num.seqs <- function(taxid, nodes, max.len=25000, max.seqs.per.spec=100000) {
    cat("Retreiving sequence counts for taxid", taxid, "\n")
    count <- 0

    ## get subtree counts. If that is smaller than max.seqs.per.spec, we can take that number
    subtree.count <- .num.seqs.for.taxid(taxid, direct=FALSE, max.len=max.len)
    if (subtree.count <= max.seqs.per.spec) {
        return (subtree.count)
    }

    ## if it is more sequences, we have to get the children's counts
    node.ranks <- c('species', 'subspecies', 'varietas', 'forma')
    current.rank <- .rank.for.taxid(taxid)
    if (current.rank %in% node.ranks) {
        ## retreive sequence counts
        c <- .num.seqs.for.taxid(taxid, direct=TRUE, max.len=max.len)
        ## if exceeded maximum amount per spec (or lower), return the count that will be retreived
        if (c > max.seqs.per.spec) {
            c <- max.seqs.per.spec
        }
        count <- count + c
    }
    for (ch in .children(taxid, nodes)) {
        count <- count + .rec.num.seqs(ch, nodes, max.len, max.seqs.per.spec)
    }
    return(count)
}

.get.blast.results <- function(taxon, seqs) {
    cat("Performing all vs all BLAST for", length(seqs), "sequences\n")
    ## make unique name for BLAST files
    ## TODO: Change this to tempdir later
    dir = './blast'
    if (! file.exists(dir)) {
        dir.create(dir)
    }
    dbfile <- paste0('taxon-', taxon, '-db.fa')
    outfile <- paste0('taxon-', taxon, '-blastout.txt')
    make.blast.db(seqs, dbfile=dbfile, dir=dir)
    blast.results <- blast.all.vs.all(dbfile, outfile=outfile, dir=dir)
    cat("Number of BLAST results ", nrow(blast.results), "\n")
    cat("Filtering BLAST results\n")
    filtered.blast.results <- filter.blast.results(blast.results, seqs)

    return(filtered.blast.results)
}
