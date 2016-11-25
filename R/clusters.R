## script to populate the phylota 'clusters' table
## Dependencies: tables 'accession2taxid', 'nodes'

library('RSQLite')
library('igraph')
source('blast.R')
source('db.R')
source('ncbi-remote.R')
source('query-local.R')

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

clusters.ci_gi.seqs.create <- function(root.taxon, nodesfile,
                                       files=list(clusters='clusters.tsv', ci_gi='ci_gi.tsv', seqs='seqs.tsv'),
                                       max.seqs=20000, informative=FALSE) {
    
    ## load nodes table    
    nodes <- read.table(nodesfile, header=T)
    
    ## load clusters that are already calculated, if any
    clusters.file <- files[['clusters']]
    ci_gi.file <- files[['ci_gi']]
    seqs.file <- files[['seqs']]
    written.clusters=NULL
    if (file.exists(clusters.file)) {
        written.clusters <- read.table(clusters.file, header=T)
    }

    ## get all nodes that will be processed:
    ## all nodes for which all children contain < max.seqs sequences
    taxa.to.process <- vector()
    queue <- root.taxon
    while(length(queue) > 0) {

        currentid <- head(queue, 1)
        queue <- tail(queue, length(queue)-1)

        ## check if cluster was calculated before
        if (! is.null(written.clusters)) {
            if (currentid %in% written.clusters$ti_root) {
                cat("Clusters for taxid", currentid, "already calculated -skipping-\n")
                next
            }
        }

        ## get number of sequences to determine if it is manageable to calculate clusters        
        num.nonmodel.seqs <- nodes[match(currentid, nodes$ti),'n_gi_sub_nonmodel']
        
        cat("Number of non-model seqs for taxid ", currentid, ":", num.nonmodel.seqs,  "\n")

        if ( num.nonmodel.seqs > max.seqs) {
            cat("Too many seqs; retreiving children \n")
            queue <- c(queue, .children(currentid, nodes))
        }
        else if (num.nonmodel.seqs == 0) {
            cat("No sequences for taxid ", currentid, "\n")
        }
        else {
            cat("Will process ", currentid, "\n")
            taxa.to.process <- c(taxa.to.process, currentid)
        }
    }
    
##    foreach (i=seq_along(taxa.to.process), .verbose=T) %dopar% {
    for (i in seq_along(taxa.to.process)) {
        taxid <- taxa.to.process[i]
        cat("Processing taxid ", taxid, " # ", i, " / ", length(taxa.to.process), "\n")
        ## generate data frames with entries for seqs, clusters, and ci_gi
        
        leaves <- .leaves.with.seqs(taxid, nodes)                
        cat("Getting sequences for", length(leaves), "terminal child taxa of taxid ", taxid, " :\n", paste(leaves, collapse="\n"), "\n")
        seqs <- unlist(lapply(leaves, .seqs.for.taxid), recursive=FALSE)
        
        seqdf <- .make.seq.entries(seqs)
        clusters <- .make.clusters(taxid, nodes, seqs, informative=informative)
        cldf <- .make.cluster.entries(clusters)
        cigidf <- .make.ci_gi.entries(clusters)

        cat("Taxid", taxid, ": writing", nrow(cldf), "clusters,", nrow(seqdf), "sequences,", nrow(cigidf), "to file\n")
        write.table(cldf, file=clusters.file, append=file.exists(clusters.file), quote=F, sep="\t", row.names=F, col.names=!file.exists(clusters.file))
        write.table(seqdf, file=seqs.file, append=file.exists(seqs.file), quote=F, sep="\t", row.names=F, col.names=!file.exists(seqs.file))
        write.table(cigidf, file=ci_gi.file, append=file.exists(ci_gi.file), quote=F, sep="\t", row.names=F, col.names=!file.exists(ci_gi.file))
        cat("Finished processing taxid ", taxid, " # ", i, " / ", length(taxa.to.process), "\n")
    }
}

.leaves.with.seqs <- function(taxid, nodes, max.seqs=10000) {
    queue <- taxid
    result <- numeric()
    while (length(queue) > 0) {
        currentid <- head(queue, 1)
        queue <- tail(queue, length(queue)-1)
        children <- .children(currentid, nodes)
        if (length(children) > 0) {
            child.nodes <- nodes[match(children, nodes$ti),]        
            ## Get nodes that have n_gi_node seqs but less than max.seqs
            leaves <- child.nodes[child.nodes$n_gi_node > 0 & child.nodes$n_gi_node < max.seqs ,'ti']        
            result <- c(result, leaves)
            queue <- c(queue, children)
        }
    }
    return(result)
}

.make.clusters <- function(taxon, nodes, seqs=NULL, blast.results=NULL, recursive=TRUE, informative=FALSE) {

    cat("Making clusters for taxid ", taxon, "\n")
    ## retreive sequences if not done so before
    gis <- vector()
    filtered.blast.results <- data.frame()
    if (is.null(seqs)) {
        cat("Retrieving sequences for taxon ", taxon, "\n")
        ## get the leaves that have sequences, but less sequences than model organisms
        leaves <- .leaves.with.seqs(taxon, nodes)        
        seqs <- unlist(lapply(leaves, .seqs.for.taxid), recursive=FALSE)
        gis <- names(seqs)
    }
    else {
        gis <- .gis.for.taxid(taxon)
    }
    ## gis will also contain the ones of model orgenisms, therefore only  take the ones that are given in our sequences
    seqs <- seqs[gis[which(gis %in% names(seqs))]]
    gis <- names(seqs)

    ## Handle empty and singleton clusters
    if (length(gis)==0) {
        return(list())
    }
    if (length(gis)==1) {
        singleton <- list(list(gis=gis, seed_gi=gis))
        singleton <- .add.cluster.info(singleton, nodes, taxon, seqs)
        return(singleton)
    }

    ## if we have not yet blasted the sequences, do so
    if (is.null(blast.results)) {
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
    }
    else {
        ##  reduce blast results such that include only the gis for the current taxon
        cat("Using BLAST results from parent cluster\n")
        filtered.blast.results <- blast.results[which(blast.results$subject.id %in% gis),]
        filtered.blast.results <- filtered.blast.results[which(filtered.blast.results$query.id %in% gis),]
    }
    ## sequence clusters are stored in a list of lists, named by gi
    ## Get top-level clusters
    cat("Clustering BLAST results for taxid", taxon, "\n")
    clusters <- cluster.blast.results(filtered.blast.results, informative=informative)
    cat("Finished clustering BLAST results for taxid", taxon, "\n")
    ## exit if there are no clusters
    if (length(clusters) < 1) {
        return(list())
    }

    ## make dataframe with fields as in PhyLoTa database
    clusters <- .add.cluster.info(clusters, nodes, taxon, seqs)
    cat("Generated ", length(clusters), "clusters\n")

    ## retrieve clusters for child taxa
    if (recursive) {
        for (ch in .children(taxon, nodes)) {
            cat("Processing child ", ch, "\n")
            
            ## child could be a model organism. In this case we don't have sequences and
            ## a separate blast search has to be conducted.
            if ( nodes[match(ch, nodes$ti),'n_gi_node'] > 10000 ) {
                clusters <- .process.model.organism(ch, clusters, nodes, seqs)
            } 
            else {
                ## calculate clusters for child taxon
                child.clusters <- .make.clusters(ch, nodes, seqs, filtered.blast.results)
                ## two parameters can only be calculated if we have cluster info on multiple taxonomic levels:
                ##  n_child, the number of child clusters, ci_anc, the parent cluster.
                ##  Since we are dealing with single-likeage clusters, we can identify the parent cluster of a
                ##  cluster if at least one gi of both clusters is the same.
                child.clusters <- lapply(child.clusters, function(cc) {
                    ## indices of the parent clusters that contain a gi of the child clusters
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
        }
    }
    cat("Finished making ", length(clusters), " clusters for taxid ", taxon, "\n")
    return(clusters)
}

.process.model.organism <- function(taxid, clusters, nodes, seqs) {


    cluster.df <- .make.cluster.entries(clusters)
    ci_gi.df <- .make.ci_gi.entries(clusters)
    ## Retreive sequences for taxid. If there are too many, pick random ones
    
    ## Make blast db for all clusters of parent node
    ## pick the first parent for which we have cluster data
    parent <- taxid
    while( (! is.na(parent)) && (! parent %in% ci_gi.df$ti) ) {
        parent <- nodes[match(parent, nodes$ti), 'ti_anc']
    }
    ci_gi.subset <- ci_gi.df[which(ci_gi.df$ti==parent),]

    ## iterate over all clusters for parent node
#    for (clustid in unique(ci_gi.subset$clustid)) {
        
#    }
    
}

## returns a data frame with columns matchine the fields in PhyLoTA's cluster table
.make.cluster.entries <- function(clusters) {
    ## remove gi and id fields which won't be part of cluster table
    l <- lapply(clusters, function(cl)cl[which(!names(cl)%in%c('gis', 'tis', 'unique_id'))])
    ## create data frame
    df <- do.call(rbind, lapply(l, as.data.frame))
    return(df)
}

.make.seq.entries <- function(seqs) {
    do.call(rbind, lapply(seqs, as.data.frame))
}

## input: List of clusters
.add.cluster.info <- function(clusters, nodes, taxid, seqs) {

    ## we will take the column names in the database as names in the list
    for (i in 1:length(clusters)) {
        cl <- clusters[[i]]
        cl.seqs <- lapply(cl$gis, function(gi)seqs[[as.character(gi)]])
        cl$ti_root <- taxid
        cl$ci <- i-1
        cl$cl_type <- ifelse(length(.children(taxid, nodes)) > 0, 'subtree', 'node')
        cl$n_gi <- length(cl$gis)
        ## all taxon ids for gis in cluster
        cl$tis <- sapply(cl.seqs, '[[', 'ti')
        cl$n_ti <- length(unique(cl$tis))
        ## sequence lengths
        l <- sapply(cl.seqs, '[[', 'length')
        cl$MinLength <- min(l)
        cl$MaxLength <- max(l)
        ## get number of genera
        cl$n_gen <- length(unique(sapply(cl$tis, .genus.for.taxid, nodes)))
        ## n_child and ci_anc will be set (or incremented) later, when multiple hierarchies are calculated
        cl$n_child <- 0
        cl$ci_anc <- NA
        ## make a unique cluster id consisting of seed gi, taxon id, cluster id, cluster type
        unique.id <- paste0(cl$seed_gi, '-', taxid, '-', cl$ci, '-', cl$cl_type)
        cl$unique_id <- unique.id
        clusters[[i]] <- cl
    }
    return(clusters)
}

cluster.blast.results <- function(blast.results, informative=T) {
    g <- graph.data.frame(blast.results[,c("query.id", "subject.id")], directed=F)
    clusters <- clusters(g)

    ## filter for phylogenetically informative clusters
    if (informative) {
        clusters <- clusters$membership[which(clusters$membership %in% which (clusters$csize > 2))]
    } else {
        clusters <- clusters$membership
    }
    ## we will return a list, one entry with sequence IDs for each cluster
    cluster.list <- lapply(unique(clusters), function(x) {
        list(gis=sort(names(clusters)[which(clusters==x)]))
    })

    ## Get the seed gi, we will chose it to be the sequence in the cluster that has
    ## the most hits with the other members in the cluster; i.e. the most connected
    ## node in the graph
    degrees <- degree(g)
    ## get seed gis and as field to clusters
    cluster.list <- lapply(cluster.list, function(cl){
        idx <- order(degrees[cl$gis], decreasing=T)[1] ## index of most connected component
        cl$seed_gi <- cl$gis[idx]
        cl
    })
    return (cluster.list)
}

.children <- function(taxid, nodes) {
    nodes[nodes[,'ti_anc']==taxid,'ti']
}

.genus.for.taxid <- function(taxid, nodes) {
    nodes[nodes[,'ti_anc']==taxid,'ti_genus']
}

##.genera.for.taxids <- function(taxids) {
##    db <- .db()
##    query <- paste('select ti_genus from nodes where ti in (', paste(taxids, collapse=','), ')')
##    l <- dbGetQuery(db, query)
##    return(l[[1]])
##}

