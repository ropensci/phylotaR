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

remote.nodes.create <- function(root.taxa = c(33090, 4751, 33208), file.name='nodes.tsv', cores=3) {
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
nodes.create <- function(dbloc, taxdir, root.taxa = c(33090, 4751, 33208), file.name='nodes.tsv', cores=3) {
    db <- .db(dbloc)

    ## Get data from NCBI taxonomy dump
    ncbi.nodes <- getnodes(taxdir)
    ncbi.names <- getnames(taxdir)

    ## Nodes that are alreay in the table can be skipped.
    ## Try to load nodes table and collect ids of nodes
    ## that are already processed
    ids.skip <- vector()
    if (file.exists(file.name)) {
        df <- read.table(file.name, header=T)
        ids.skip <- df$ti
    }
    
    ## Do calculation in parallel; pick roughly a number of taxon ids
    ## equal to twice the number of CPUs available, to make
    ## sure we won't have unused CPUs
    ids.to.process <- root.taxa
    ids.not.processed <- vector()
       
    ## CAUTION: loop below could get stuck if there are not enough echildren and many cores!
    while (length(ids.to.process) < (cores*2) && length(ids.to.process > 0)) {
        ## skip if node entry was already generated: 
        currentid <- head(ids.to.process, 1)
        cat("CurrentID : ", currentid, "\n")        

        ## do not go further down if we are at genus or lower, because we would't get the field ti_genus of the children!
        if (getrank(currentid, NULL, nodes=ncbi.nodes) %in% c('genus', 'species', 'subspecies', 'varietas', 'forma'))
            break

        ## remove from queue
        ids.to.process <- tail(ids.to.process, length(ids.to.process)-1)

        if (currentid %in% ids.skip) {
            cat("Node with ID ", currentid, " already processed. Skipping \n");
            next
        }

        for (ch in children(currentid, ncbi.nodes)) {
            ids.to.process <- c(ids.to.process, ch)
        }
        ids.not.processed <- c(ids.not.processed, currentid)
    }

    ##    for (tax in ids.to.process) {
    ##nodes.written <- data.frame()
    cat("Adding ", length(ids.to.process), " nodes recursively, in parallel \n")
    ##for (i in seq_along(ids.to.process)) {
    nodes.written <- foreach (i=seq_along(ids.to.process), .combine=rbind) %dopar% {
        tax <- ids.to.process[i]
        cat("Recursively processing taxid ", tax, " # ", i, " / ", length(ids.to.process), "\n")
        n <- .add.stats(tax, ncbi.nodes, recursive=TRUE)
        ## only return rows for which we collected stats, all others have NA in num.seqs.node
        current.nodes <- n[which(! is.na(n$num.seqs.node)),]
        .write.row(current.nodes, ncbi.names, file.name)
        cat("Finished processing taxid ", tax, " # ", i, " / ", length(ids.to.process), "\n")
        ##nodes.written <- rbind(nodes.written, current.nodes)
        current.nodes
    }
    ## add the top nodes non-recursively
    cat("Adding ", length(ids.not.processed), " top-level nodes non-recursively, in parallel \n")
    ##for (i in seq_along(ids.not.processed)) {
    foreach(i=seq_along(ids.not.processed)) %dopar% {
        tax <- ids.to.process[i]
        cat("Processing taxid ", tax, " # ", i, " / ", length(ids.not.processed), "\n")
        n <- .add.stats(tax, nodes.written, recursive=FALSE)
        current.nodes <- n[match(tax, n$id),]
        .write.row(current.nodes, ncbi.names, file.name)
        cat("Finished processing taxid ", tax, " # ", i, " / ", length(ids.not.processed), "\n")
    }
}

.write.row <- function(nodes, ncbi.names, file.name) {
    ## Prepare for writing to file
    ## add other column values
    taxon.name <- as.character(ncbi.names[match(nodes$id, ncbi.names$id),"name"])
    commons <- ncbi.names[grep ('common', ncbi.names$type),]
    common.name <- as.character(commons$name[match(nodes$id, commons$id)])

    ## Create dataframe with correct column names that will be
    ## inserted into the database. Some fields (such as cluster info) will need additional
    ## information and will be filled later
    nodes.df <- data.frame(ti=as.integer(nodes$id),
                           ti_anc=as.integer(nodes$parent),
                           terminal_flag=NA,
                           rank_flag=as.integer(! nodes$rank=='no rank'),
                           model=NA,
                           taxon_name=taxon.name,
                           common_name=common.name,
                           rank=as.character(nodes$rank),
                           n_gi_node=as.integer(nodes$num.seqs.node),
                           n_gi_sub_nonmodel=as.integer(nodes$num.seqs.subtree.nonmodel),
                           n_gi_sub_model=as.integer(nodes$num.seqs.subtree.model),
                           n_clust_node=NA,
                           n_clust_sub=NA,
                           n_PIclust_sub=NA,
                           n_sp_desc=as.integer(nodes$num.spec.desc),
                           n_sp_model=as.integer(nodes$num.spec.model),
                           n_leaf_desc=as.integer(nodes$num.leaf.desc),
                           n_otu_desc=as.integer(nodes$num.otu.desc),
                           ti_genus=as.integer(nodes$ti.genus),
                           n_genera=as.integer(nodes$num.genera)
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

## Recursive function to add number of sequences for each taxon to table as produced by getnodes().
## We set sequence counts for 'node' and 'subtree', if a node is of rank species or lower. For subtree,
## this includes all desendants. We also distinguish between model (>=20000 seqs per node) and nonmodel organisms.
.add.stats <- function(taxid, nodes, recursive=FALSE) {

    ## ti_genus might have been passed down from parent node
    ti.genus <- nodes[which(nodes$id==taxid), 'ti.genus']
    stats <- c(num.seqs.node=0,
               num.seqs.subtree.model=0,
               num.seqs.subtree.nonmodel=0,
               num.leaf.desc=0, ## also counts itsself, for a leaf it is 1
               num.spec.desc=0, ## also counts itsself, for a species it is 1
               num.spec.model=0,
               ti.genus= ifelse(is.null(ti.genus), NA, ti.genus),
               num.genera=0,
               num.otu.desc=1)
    stats <- sapply(stats, as.numeric)
    ## Ranks for which we directly get the species count;
    ## for ranks above this, we will get the sum of their children.
    ## For taxa with these ranks, there will also be a 'n_gi_node'
    node.ranks <- c('species', 'subspecies', 'varietas', 'forma')

    ## Recursion anchor
    rank <- getrank(taxid, NULL, nodes=nodes)
    if (rank %in% node.ranks) {
        curr.num.seqs <- .num.seqs.for.taxid(taxid)
        num.seqs.node <- curr.num.seqs
        stats['num.seqs.node'] <- curr.num.seqs
        if (curr.num.seqs >= 20000) {
            num.seqs.subtree.model <- curr.num.seqs
            stats['num.seqs.subtree.model'] <- curr.num.seqs
        } else {
            num.seqs.subtree.nonmodel <- curr.num.seqs
            stats['num.seqs.subtree.nonmodel'] <- curr.num.seqs
        }
        if (rank == 'species') {
            stats['num.spec.desc'] <- 1
            if (curr.num.seqs >= 20000) {
                stats['num.spec.model'] <- 1
            }
        }
    }

    if (rank == 'genus') {
        stats['ti.genus'] <- taxid
    }

    ch <- children(taxid, nodes)
    if (length(ch)==0) {
        stats['num.leaf.desc'] <- stats['num.leaf.desc'] + 1
    }

    ## add counts from children
    for (child in ch) {
        ## ti.genus has to be passed down the tree and not upwards...
        if (! is.na(stats['ti.genus'])) {
             nodes[which(nodes$id==child),'ti.genus'] <- stats['ti.genus']
        }
        if (recursive == TRUE) {
            ## call function recursively for children
            nodes <- .add.stats(child, nodes, recursive=TRUE)
        }
        ## add subtree sequence counts for higher level nodes
        if (! rank %in% node.ranks) {
            stats['num.seqs.subtree.model'] <- stats['num.seqs.subtree.model'] + nodes[which(nodes$id==child),'num.seqs.subtree.model']
            stats['num.seqs.subtree.nonmodel'] <- stats['num.seqs.subtree.nonmodel'] + nodes[which(nodes$id==child),'num.seqs.subtree.nonmodel']
        }
        if (getrank(child, NULL, nodes) == 'genus') {
            stats['num.genera'] <- stats['num.genera'] + 1
        }
        stats['num.leaf.desc'] <- stats['num.leaf.desc'] + nodes[which(nodes$id==child),'num.leaf.desc']
        stats['num.otu.desc'] <- stats['num.otu.desc'] + nodes[which(nodes$id==child),'num.otu.desc']
        stats['num.spec.desc'] <- stats['num.spec.desc'] + nodes[which(nodes$id==child),'num.spec.desc']
        stats['num.spec.model'] <- stats['num.spec.model'] + nodes[which(nodes$id==child),'num.spec.model']
        stats['num.genera'] <- stats['num.genera'] + nodes[which(nodes$id==child),'num.genera']
    }
    ## add stats to data frame
    for (n in names(stats)) {
        nodes[which(nodes$id==taxid),n] <- stats[n]
    }

    cat("Added species counts for taxid ", taxid, "\n")
    return(nodes)
}

.add.stats.rem <- function(taxid, nodes, recursive=FALSE) {

    ## ti_genus might have been passed down from parent node

    ti.genus <- nodes[which(nodes$id==taxid), 'ti.genus']
    
    stats <- c(num.seqs.node=0,
               num.seqs.subtree.model=0,
               num.seqs.subtree.nonmodel=0,
               num.leaf.desc=0, ## also counts itsself, for a leaf it is 1
               num.spec.desc=0, ## also counts itsself, for a species it is 1
               num.spec.model=0,
               ti.genus= ifelse(is.null(ti.genus), NA, ti.genus),
               num.genera=0,
               num.otu.desc=1)

    ## Ranks for which we directly get the species count;
    ## for ranks above this, we will get the sum of their children.
    ## For taxa with these ranks, there will also be a 'n_gi_node'
    node.ranks <- c('species', 'subspecies', 'varietas', 'forma')

    ## Recursion anchor
    rank <- nodes[match(taxid, nodes$uid),'rank']
    if (rank %in% node.ranks) {
        curr.num.seqs <- .num.seqs.for.taxid(taxid)
        num.seqs.node <- curr.num.seqs
        stats['num.seqs.node'] <- curr.num.seqs
        if (curr.num.seqs >= 20000) {
            num.seqs.subtree.model <- curr.num.seqs
            stats['num.seqs.subtree.model'] <- curr.num.seqs
        } else {
            num.seqs.subtree.nonmodel <- curr.num.seqs
            stats['num.seqs.subtree.nonmodel'] <- curr.num.seqs
        }
        if (rank == 'species') {
            stats['num.spec.desc'] <- 1
            if (curr.num.seqs >= 20000) {
                stats['num.spec.model'] <- 1
            }
        }
    }

    if (rank == 'genus') {
        stats['ti.genus'] <- taxid
    }

    ch <- taxize::children(taxid, db='ncbi')[[1]]
    if (nrow(ch)==0) {
        stats['num.leaf.desc'] <- stats['num.leaf.desc'] + 1
    }
    
    ## add counts from children
    for (child in ch$childtaxa_id) {

        ## extend nodes for child
        nodes <- rbind(nodes, ncbi_get_taxon_summary(child))
        ## ti.genus has to be passed down the tree and not upwards...
        if (! is.na(stats['ti.genus'])) {
             nodes[which(nodes$uid==child),'ti.genus'] <- stats['ti.genus']
        }
        if (recursive == TRUE) {
            ## call function recursively for children
            nodes <- .add.stats.rem(child, nodes, recursive=TRUE)
        }
        ## add subtree sequence counts for higher level nodes
        if (! rank %in% node.ranks) {
            stats['num.seqs.subtree.model'] <- stats['num.seqs.subtree.model'] + nodes[which(nodes$uid==child),'num.seqs.subtree.model']
            stats['num.seqs.subtree.nonmodel'] <- stats['num.seqs.subtree.nonmodel'] + nodes[which(nodes$uid==child),'num.seqs.subtree.nonmodel']
        }
        if (nodes[match(taxid, nodes$uid),'rank'] == 'genus') {
            stats['num.genera'] <- stats['num.genera'] + 1
        }
        stats['num.leaf.desc'] <- stats['num.leaf.desc'] + nodes[which(nodes$uid==child),'num.leaf.desc']
        stats['num.otu.desc'] <- stats['num.otu.desc'] + nodes[which(nodes$uid==child),'num.otu.desc']
        stats['num.spec.desc'] <- stats['num.spec.desc'] + nodes[which(nodes$uid==child),'num.spec.desc']
        stats['num.spec.model'] <- stats['num.spec.model'] + nodes[which(nodes$uid==child),'num.spec.model']
        stats['num.genera'] <- stats['num.genera'] + nodes[which(nodes$uid==child),'num.genera']
    }
    ## add stats to data frame
    for (n in names(stats)) {
        nodes[which(nodes$uid==taxid),n] <- stats[n]
    }

    cat("Added species counts for taxid ", taxid, "\n")
    return(nodes)
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

children <- function(id, nodes) {
    return(nodes$id[which(nodes$parent==id)])
}

