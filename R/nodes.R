## script to populate the phylota 'nodes' table
## Dependency: Table 'accession2taxid'

library("CHNOSZ")
library('rentrez')
source('db.R')
source('ncbi-remote.R')

## Schema in phylota database:
## CREATE TABLE "nodes_194" (
##  "ti" int(10)  NOT NULL,
##  "ti_anc" int(10)  DEFAULT NULL,
##  "terminal_flag" tinyint(1) DEFAULT NULL,
##  "rank_flag" tinyint(1) DEFAULT NULL,
##  "model" tinyint(1) DEFAULT NULL,
##  "taxon_name" varchar(128) DEFAULT NULL,
##  "common_name" varchar(128) DEFAULT NULL,
##  "rank" varchar(64) DEFAULT NULL,
##  "n_gi_node" int(10)  DEFAULT NULL,
##  "n_gi_sub_nonmodel" int(10)  DEFAULT NULL,
##  "n_gi_sub_model" int(10)  DEFAULT NULL,
##  "n_clust_node" int(10)  DEFAULT NULL,
##  "n_clust_sub" int(10)  DEFAULT NULL,
##  "n_PIclust_sub" int(10)  DEFAULT NULL,
##  "n_sp_desc" int(10)  DEFAULT NULL,
##  "n_sp_model" int(10)  DEFAULT NULL,
##  "n_leaf_desc" int(10)  DEFAULT NULL,
##  "n_otu_desc" int(10)  DEFAULT NULL,
##  "ti_genus" int(10)  DEFAULT NULL,
##  "n_genera" int(10)  DEFAULT NULL,
##  PRIMARY KEY ("ti")
##);

remote.nodes.create <- function(root.taxa = c(33090, 4751, 33208), file.name='nodes.tsv') {
    require('taxize')
    nodes <- do.call(rbind, lapply(root.taxa, ncbi_get_taxon_summary))
    nodes <- transform(nodes, uid=as.numeric(uid))

    for (tax in nodes$uid) {
        n <- .add.stats.rem(tax, nodes, recursive=TRUE)
    }
}

## The table is created from the 'nodes' table in the NCBI taxonomy
## Therefore, a directory with the path where the NCBI taxonomy dump
## is located must be provided
nodes.create <- function(taxdir, root.taxa = c(33090, 4751, 33208), file.name='nodes.tsv', model.threshold=10000) {

    ## Get data from NCBI taxonomy dump
    if (! exists('ncbi.nodes'))
        ncbi.nodes <<- getnodes(taxdir)
    if (! exists('ncbi.names'))
        ncbi.names <<- getnames(taxdir)

    set <- get.manageable.node.set(root.taxa, ncbi.nodes, max.descendants=10000, timeout=10,
                                   nodesfile=file.name)
    ids.to.process <- set$manageable.nodes
    ids.not.processed <- set$rejected.nodes

    ## Nodes that are alreay in the table can be skipped.
    ## Try to load nodes table and collect ids of nodes
    ## that are already processed
    #ids.skip <- vector()
    #if (file.exists(file.name)) {
    #    df <- read.table(file.name, header=T)
    #    ids.skip <- df$ti
    #    ## remove from taxa to process
    #    ids.to.process <- ids.to.process[! ids.to.process %in% ids.skip]
    #    ids.not.processed <- ids.not.processed[! ids.not.processed %in% ids.skip]
    #}

    cat("Adding ", length(ids.to.process), " nodes recursively, in parallel \n")
    for (i in seq_along(ids.to.process)) {
##    foreach (i=seq_along(ids.to.process)) %dopar% {
        tax <- ids.to.process[i]
        cat("Recursively processing taxid ", tax, " # ", i, " / ", length(ids.to.process), "\n")
        start <- data.frame()
        if (file.exists(file.name)) {
            start <- read.table(file.name, header=T)
            start <- start[,c('ti','ti_anc','rank','n_gi_node','n_gi_sub_nonmodel','n_gi_sub_model',
                                              'n_sp_desc', 'n_sp_model', 'n_leaf_desc', 'n_otu_desc', 'ti_genus',
                                              'n_genera')]
        }
        else {
            start <- .init.nodes.df()
        }

        n <- .add.stats(tax, start, recursive=TRUE, model.threshold)
        ids.to.write <- n$ti[which(!n$ti %in% start$ti)]
        stopifnot(tax %in% ids.to.write)
        ## only return rows which were not written to file before
        subset <-n[match(ids.to.write, n$ti), ]
        .write.row(subset, ncbi.names, file.name)
        cat("Finished processing taxid ", tax, " # ", i, " / ", length(ids.to.process), "\n")
    }

    ## Add the top nodes non-recursively. This has to be in reversed order to make sure a parent is not
    ## inserted before it's children, since we need the info from the children. Therefore, this must not
    ## happen in parallel!
    cat("Adding ", length(ids.not.processed), " top-level nodes non-recursively, sequentially \n")
    for (i in seq_along(ids.not.processed)) {
        tax <- rev(ids.not.processed)[i]
        ##foreach(i=seq_along(ids.not.processed)) %dopar% {
        ## reload nodes that have been written before
        nodes.written <- read.table(file.name, header=T)
        nodes.written <- nodes.written[,c('ti','ti_anc','rank','n_gi_node','n_gi_sub_nonmodel','n_gi_sub_model',
                                          'n_sp_desc', 'n_sp_model', 'n_leaf_desc', 'n_otu_desc', 'ti_genus',
                                          'n_genera')]

        cat("Processing taxid ", tax, " # ", i, " / ", length(ids.not.processed), "\n")
        n <- .add.stats(tax, nodes.written, recursive=FALSE, model.threshold)
        .write.row(n[match(tax, n$ti),], ncbi.names, file.name)
        cat("Finished processing taxid ", tax, " # ", i, " / ", length(ids.not.processed), "\n")
    }
}


## Given a root taxon, returns a set of nodes that only have up to a maximum
## number of descendant nodes
get.manageable.node.set <- function(root.taxa, ncbi.nodes, max.descendants=10000, timeout=10, nodesfile="") {
    queue <- root.taxa
    manageable.nodes <- vector()
    num.descendants <- vector()
    rejected.nodes <- vector()

    nodes.written <- data.frame()
    if (file.exists(nodesfile)) {
        nodes.written <- read.table(nodesfile, header=T)
    }

    total.count <- 0
    while(length(queue) > 0) {
        id <- head(queue, 1)

        ## remove current node from queue
        queue <- tail(queue, length(queue)-1)

        if (nrow(nodes.written) > 0) {
            if (id %in% nodes.written$ti) {
                cat("Taxon", id, "already in file", nodesfile, " -skipping- \n")
                next
            }
        }

        ## If there are more descendants than max.descendants or retrieveing
        ## the number takes longer than <timeout> secounds, discard the node
        ## and add its children to the queue
        n.desc <- .eval.time.limit(num.descendants(id, ncbi.nodes), cpu=timeout)
        if (is.null(n.desc) || n.desc > max.descendants) {
            queue <- c(queue, children(id, ncbi.nodes))
            rejected.nodes <- c(rejected.nodes, id)
            cat("Taxon ", id, " has too many descendants or timeout reached counting descendants. Processing child taxa.\n")
        }
        else {
            manageable.nodes <- c(manageable.nodes, id)
            num.descendants <- c(num.descendants, n.desc)
            cat("Taxon ", id, " has maneagable number of descendants: ", n.desc, "\n")
            total.count <- total.count + n.desc
            cat("Current number of nodes to be processed: ", total.count, "\n")
        }
    }

    return(list(manageable.nodes=manageable.nodes,
                rejected.nodes=rejected.nodes,
                num.descendants=num.descendants))
}

.write.row <- function(nodes, ncbi.names, file.name) {
    ## Prepare for writing to file
    ## add other column values
    taxon_name<- as.character(ncbi.names[match(nodes$ti, ncbi.names$id),"name"])
    commons <- ncbi.names[grep ('common', ncbi.names$type),]
    common_name <- as.character(commons$name[match(nodes$ti, commons$id)])

    ## Create dataframe with correct column names that will be
    ## inserted into the database. Some fields (such as cluster info) will need additional
    ## information and will be filled later
    nodes.df <- data.frame(ti=as.integer(nodes$ti),
                           ti_anc=as.integer(nodes$ti_anc),
                           terminal_flag=NA,
                           rank_flag=as.integer(! nodes$rank=='no rank'),
                           model=NA,
                           taxon_name=taxon_name,
                           common_name=common_name,
                           rank=as.character(nodes$rank),
                           n_gi_node=as.integer(nodes$n_gi_node),
                           n_gi_sub_nonmodel=as.integer(nodes$n_gi_sub_nonmodel),
                           n_gi_sub_model=as.integer(nodes$n_gi_sub_model),
                           n_clust_node=NA,
                           n_clust_sub=NA,
                           n_PIclust_sub=NA,
                           n_sp_desc=as.integer(nodes$n_sp_desc),
                           n_sp_model=as.integer(nodes$n_sp_model),
                           n_leaf_desc=as.integer(nodes$n_leaf_desc),
                           n_otu_desc=as.integer(nodes$n_otu_desc),
                           ti_genus=as.integer(nodes$ti_genus),
                           n_genera=as.integer(nodes$n_genera)
                           )
    write.table(nodes.df, file=file.name, sep="\t", row.names=F, append=file.exists(file.name), col.names=!file.exists(file.name))
    cat("Wrote ", nrow(nodes.df), " entries to file\n")
}

## TODO: make functionality to insert in DB
##    repeat {
##        ret <- try(dbWriteTable(conn=db, name='nodes', value=nodes.df, row.names=F, overwrite=overwrite, append=append))
##        if(!is(ret, "try-error")) break
##    }
##    dbDisconnect(db)
##    rm(db)
##    return(ret)

.init.nodes.df <- function() {
    nodes.df <- data.frame(ti=numeric(),
                           ti_anc=numeric(),
                           rank=character(),
                           n_gi_node=numeric(),
                           n_gi_sub_nonmodel=numeric(),
                           n_gi_sub_model=numeric(),
                           n_sp_desc=numeric(),
                           n_sp_model=numeric(),
                           n_leaf_desc=numeric(),
                           n_otu_desc=numeric(),
                           ti_genus=numeric(),
                           n_genera=numeric()
                           )
    return(nodes.df)
}

## Recursive function to add number of sequences for each taxon to table as produced by getnodes().
## We set sequence counts for 'node' and 'subtree', if a node is of rank species or lower. For subtree,
## this includes all desendants. We also distinguish between model (>=XXX seqs per node) and nonmodel organisms.
.add.stats <- function(taxid, nodes, recursive=FALSE, model.threshold) {

    if (taxid %in% nodes$ti) {
        cat("Stats for node", taxid, "already there. Skipping \n")
        return(nodes)
    }

    cat("Going to add species counts for taxid ", taxid, "\n")
    ## create node vector storing all relevant info
    nodeinfo <- ncbi.nodes[match(taxid, ncbi.nodes$id),]

    stats <- data.frame(ti=taxid,
                        ti_anc=nodeinfo$parent,
                        rank=nodeinfo$rank,
                        n_gi_node=0,
                        n_gi_sub_model=0,
                        n_gi_sub_nonmodel=0,
                        n_leaf_desc=0, ## also counts itsself, for a leaf it is 1
                        n_sp_desc=0, ## also counts itsself, for a species it is 1
                        n_sp_model=0,
                        ti_genus=.get.genus(taxid, ncbi.nodes),
                        n_genera=0,
                        n_otu_desc=1)

    ## Ranks for which we directly get the species count;
    ## for ranks above this, we will get the sum of their children.
    ## For taxa with these ranks, there will also be a 'n_gi_node'
    node.ranks <- c('species', 'subspecies', 'varietas', 'forma')

    ## Recursion anchor
    rank <- getrank(taxid, nodes=ncbi.nodes)

    if (rank %in% node.ranks) {
        curr.num.seqs <- .num.seqs.for.taxid(taxid)
        n_gi_node <- curr.num.seqs
        stats['n_gi_node'] <- curr.num.seqs
        if (curr.num.seqs >= model.threshold) {
            n_gi_sub_model <- curr.num.seqs
            stats['n_gi_sub_model'] <- curr.num.seqs
        } else {
            n_gi_sub_nonmodel <- curr.num.seqs
            stats['n_gi_sub_nonmodel'] <- curr.num.seqs
        }
        if (rank == 'species') {
            stats['n_sp_desc'] <- 1
            if (curr.num.seqs >= model.threshold) {
                stats['n_sp_model'] <- 1
            }
        }
    }

    ch <- children(taxid, ncbi.nodes)
    if (length(ch)==0) {
        stats['n_leaf_desc'] <- stats['n_leaf_desc'] + 1
    }

    ## add counts from children
    for (child in ch) {

        if (recursive == TRUE) {
            ## call function recursively for children
            nodes <- .add.stats(child, nodes, recursive=TRUE, model.threshold)
        }

        ## increase number of genera, if children are genera
        if (getrank(child, nodes=ncbi.nodes) == 'genus') {
            stats['n_genera'] <- stats['n_genera'] + 1
        }

        ## compile columns for which counts are added from the child nodes
        cols <- c('n_leaf_desc', 'n_otu_desc', 'n_sp_desc', 'n_sp_model', 'n_genera')

        ## only add subtree sequence counts if we have a higher-level node
        if (! rank %in% node.ranks) {
            cols <- c(cols, 'n_gi_sub_model', 'n_gi_sub_nonmodel')
        }

        ## add the counts
        stats[cols] <- stats[cols] + nodes[which(nodes$ti==child), cols]
    }
    ## add current node to nodes to be returned

    nodes <- rbind(nodes, stats[names(nodes)])

    cat("Added stats for taxid ", taxid, "\n")
    return(nodes)
}

.add.stats.rem <- function(taxid, nodes, recursive=FALSE) {

    stats <- c(n_gi_node=0,
               n_gi_sub_model=0,
               n_gi_sub_nonmodel=0,
               n_leaf_desc=0, ## also counts itsself, for a leaf it is 1
               n_sp_desc=0, ## also counts itsself, for a species it is 1
               n_sp_model=0,
               #Hannes
               ti_genus=-1,## .get.genus(taxid),
               n_genera=0,
               n_otu_desc=1)

    ## Ranks for which we directly get the species count;
    ## for ranks above this, we will get the sum of their children.
    ## For taxa with these ranks, there will also be a 'n_gi_node'
    node.ranks <- c('species', 'subspecies', 'varietas', 'forma')

    ## Recursion anchor
    rank <- nodes[match(taxid, nodes$uid),'rank']
    if (rank %in% node.ranks) {
        curr.num.seqs <- .num.seqs.for.taxid(taxid)
        n_gi_node <- curr.num.seqs
        stats['n_gi_node'] <- curr.num.seqs
        if (curr.num.seqs >= 10000) {
            n_gi_sub_model <- curr.num.seqs
            stats['n_gi_sub_model'] <- curr.num.seqs
        } else {
            n_gi_sub_nonmodel <- curr.num.seqs
            stats['n_gi_sub_nonmodel'] <- curr.num.seqs
        }
        if (rank == 'species') {
            stats['n_sp_desc'] <- 1
            if (curr.num.seqs >= 10000) {
                stats['n_sp_model'] <- 1
            }
        }
    }

    ch <- taxize::children(taxid, db='ncbi')[[1]]
    if (nrow(ch)==0) {
        stats['n_leaf_desc'] <- stats['n_leaf_desc'] + 1
    }

    ## add counts from children
    for (child in ch$childtaxa_id) {

        ## extend nodes for child
        nodes <- rbind(nodes, ncbi_get_taxon_summary(child))

        if (recursive == TRUE) {
            ## call function recursively for children
            nodes <- .add.stats.rem(child, nodes, recursive=TRUE)
        }
        ## add subtree sequence counts for higher level nodes
        if (! rank %in% node.ranks) {
            stats['n_gi_sub_model'] <- stats['n_gi_sub_model'] + nodes[which(nodes$uid==child),'n_gi_sub_model']
            stats['n_gi_sub_nonmodel'] <- stats['n_gi_sub_nonmodel'] + nodes[which(nodes$uid==child),'n_gi_sub_nonmodel']
        }
        if (nodes[match(taxid, nodes$uid),'rank'] == 'genus') {
            stats['n_genera'] <- stats['n_genera'] + 1
        }
        stats['n_leaf_desc'] <- stats['n_leaf_desc'] + nodes[which(nodes$uid==child),'n_leaf_desc']
        stats['n_otu_desc'] <- stats['n_otu_desc'] + nodes[which(nodes$uid==child),'n_otu_desc']
        stats['n_sp_desc'] <- stats['n_sp_desc'] + nodes[which(nodes$uid==child),'n_sp_desc']
        stats['n_sp_model'] <- stats['n_sp_model'] + nodes[which(nodes$uid==child),'n_sp_model']
        stats['n_genera'] <- stats['n_genera'] + nodes[which(nodes$uid==child),'n_genera']
    }

    ## add stats to data frame
    for (n in names(stats)) {
        nodes[which(nodes$uid==taxid),n] <- stats[n]
    }

    cat("Added stats for taxid ", taxid, "\n")
    return(nodes)
}

.get.genus <- function(taxid, ncbi.nodes) {

    ## if rank is already higher than genus, do not search for it,
    lower.ranks <- c('genus', 'subgenus', 'species group',
                     'species subgroup', 'species',
                     'subspecies','varietas', 'forma')
    higher.ranks <- c('superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum',
                      'subphylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder',
                      'order', 'suborder', 'infraorder', 'parvorder', 'superfamily', 'family',
                      'subfamily', 'tribe', 'subtribe')
    current.rank <- getrank(taxid, nodes=ncbi.nodes)
    if (! current.rank  %in% lower.ranks) {
        return(NA)
    }

    node <- taxid
    while (current.rank != 'genus') {
        node <- CHNOSZ::parent(node, nodes=ncbi.nodes)
        current.rank <- getrank(node, nodes=ncbi.nodes)

        ## some hybrids between genera can have s subfamily as a parent...
        if (current.rank %in% higher.ranks) {
            break
        }

    }
    return(node)
}

descendants <- function(id, nodes) {
    queue <- id
    result <- numeric()
    while (length(queue)>0) {
        currentid <- head(queue, 1)
        queue <- tail(queue, length(queue)-1)
        currentchildren <- children(currentid, nodes)
        result <- c(result, currentchildren)
        queue <- c(queue, currentchildren)
    }
    return(result)
}

num.descendants <- function(id, nodes) {
    queue <- id
    result <- 0
    while (length(queue)>0) {
        result <- result + length(queue)
        newqueue <- foreach (i=seq_along(queue), .combine=c) %dopar% {
            children(queue[i], nodes)
        }
        queue <- newqueue
    }
    return(result-1)
}

children <- function(id, nodes) {
    return(nodes$id[which(nodes$parent==id)])
}

.eval.time.limit <- function(expr, cpu = Inf, elapsed = Inf) {
    y <- try({setTimeLimit(cpu, elapsed); expr}, silent = TRUE)
    setTimeLimit()
    if(inherits(y, "try-error")) NULL else y
}
