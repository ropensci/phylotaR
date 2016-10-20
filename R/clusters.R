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

clusters.ci_gi.seqs.create <- function(root.taxon, file.stem=NULL, max.seqs=25000) {
    db <- .db(dbname)

    cluster.entries <- data.frame()
    seq.entries <- data.frame()
    ci_gi.entries <- data.frame()

    ## get all nodes that will be processed:
    ## all nodes for which all children contain < max.seqs sequences
    taxa.to.process <- vector()

    queue <- root.taxon
    while(length(queue) > 0) {
        
        currentid <- head(queue, 1)        
        queue <- tail(queue, length(queue)-1)
        cat("Getting seq counts for taxid ", currentid, "\n")
        if (.local.num.seqs.for.taxid(currentid) > max.seqs) {
            cat("Too many seqs; retreiving children \n")
            queue <- c(queue, .children(currentid))
        }
        else {
            cat("Will process ", currentid, "\n")
            taxa.to.process <- c(taxa.to.process, currentid)
            ## generate data frames with entries for seqs, clusters, and ci_gi
            #seqs <- .seqs.for.taxid(currentid)
            #seqdf <- .make.seq.entries(seqs)
            #clusters <- .make.clusters(currentid, seqs)
            #cldf <- .make.cluster.entries(clusters)
            #cigidf <- .make.ci_gi.entries(clusters)

            ## either append data frames to result which will be returned or
            ## append results to file imedeately
            #if (is.null(file.stem)){
            #    cluster.entries <- rbind(cluster.entries, cldf)
            #    seq.entries <- rbind(seq.entries, seqdf)
            #    ci_gi.entries <- rbind(ci_gi.entries, cigidf)
            #}
            #else {                
            #    clusters.file <- paste0(file.stem, '-clusters.tsv')
            #    write.table(cldf, file=clusters.file, append=file.exists(clusters.file), quote=F, sep="\t", row.names=F)
            #    seqs.file <- paste0(file.stem, '-seqs.tsv')
            #    write.table(seqdf, file=seqs.file, append=file.exists(seqs.file), quote=F, sep="\t", row.names=F)
            #    ci_gi.file <- paste0(file.stem, '-ci_gi.tsv')
            #    write.table(cigidf, file=ci_gi.file, append=file.exists(ci_gi.file), quote=F, sep="\t", row.names=F)
            #}
        }
##        cat("Processed taxid ", currentid, "\n")
    }    
    return (list(clusters=cluster.entries, seqs=seq.entries, ci_gi=ci_gi.entries))
}

.make.clusters <- function(taxon, seqs=NULL, parent=NULL, blast.results=NULL, recursive=TRUE) {

    cat("Making clusters for taxid ", taxon, "\n")
    ## retreive sequences if not done so before
    gis <- vector()
    filtered.blast.results <- data.frame()
    if (is.null(seqs)) {
        cat("Retrieving sequences for taxon ", taxon, "\n")
        seqs <- .seqs.for.taxid(taxon)
        gis <- names(seqs)
    }
    else {
        gis <- .gis.for.taxid(taxon)
    }
    if (length(gis)<1) {
        return(list())
    }
    seqs <- seqs[gis]
    ## if we have not yet blasted the sequences, do so
    if (is.null(blast.results)) {
        cat("Performing all vs all BLAST\n")
        dbfile <- make.blast.db(seqs)
        blast.results <- blast.all.vs.all(dbfile)
        cat("Number of BLAST results ", nrow(blast.results), "\n")
        cat("Filtering BLAST results\n")
        filtered.blast.results <- filter.blast.results(blast.results, seqs)
    }
    else {
        ##  reduce blast results such that include only the gis for the current taxon
        cat("Taking BLAST results from parent cluster\n")
        filtered.blast.results <- blast.results[which(blast.results$subject.id %in% gis),]
        filtered.blast.results <- filtered.blast.results[which(filtered.blast.results$query.id %in% gis),]
    }        
    ## sequence clusters are stored in a list of lists, named by gi
    ## Get top-level clusters
    cat("Clustering BLAST results\n")
    clusters <- cluster.blast.results(filtered.blast.results)

    ## exit if there are no clusters
    if (length(clusters) < 1) {
        return(list())
    }

    ## make dataframe with fields as in PhyLoTa database
    clusters <- .add.cluster.info(clusters, taxon, seqs)
    cat("Generated ", length(clusters), "clusters\n")
    
    ## retrieve clusters for child taxa
    if (recursive) {
        for (ch in .children(taxon)) {
            cat("Processing child ", ch, "\n")
            ## calculate clusters for child taxon
            child.clusters <- .make.clusters(ch, seqs, taxon, filtered.blast.results)
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
    cat("Finished making ", length(clusters), " clusters for taxid ", taxon, "\n")
    return(clusters)
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
.add.cluster.info <- function(clusters, taxid, seqs) {

    ## we will take the column names in the database as names in the list
    for (i in 1:length(clusters)) {
        cl <- clusters[[i]]        
        cl.seqs <- lapply(cl$gis, function(gi)seqs[[as.character(gi)]])       
        cl$ti_root <- taxid
        cl$ci <- i-1
        cl$cl_type <- ifelse(length(.children(taxid)) > 0, 'subtree', 'node')
        cl$n_gi <- length(cl$gis)
        ## all taxon ids for gis in cluster                
        cl$tis <- sapply(cl.seqs, '[[', 'ti')
        cl$n_ti <- length(unique(cl$tis))
        ## sequence lengths
        l <- sapply(cl.seqs, '[[', 'length')
        cl$MinLength <- min(l)
        cl$MaxLength <- max(l)
        ## get number of genera
        cl$n_gen <- length(unique(sapply(cl$tis, .genus.for.taxid)))
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

.children <- function(taxid) {
    db <- .db()
    query <- paste('select ti from nodes where ti_anc =', taxid)
    l <- dbGetQuery(db, query)
    return(l[[1]])
}

.genus.for.taxid <- function(taxid) {
    ##cat("Retreiving genus for taxid ", taxid, "\n")
    db <- .db()
    query <- paste('select ti_genus from nodes where ti=', taxid)
    l <- dbGetQuery(db, query)
    return(l[[1]])
}

##.genera.for.taxids <- function(taxids) {
##    db <- .db()
##    query <- paste('select ti_genus from nodes where ti in (', paste(taxids, collapse=','), ')')
##    l <- dbGetQuery(db, query)
##    return(l[[1]])
##}

