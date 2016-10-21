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

## The table is created from the 'nodes' table in the NCBI taxonomy
## Therefore, a directory with the path where the NCBI taxonomy dump
## is located must be provided
nodes.create <- function(dbloc, taxdir, root.taxa = c(33090, 4751, 33208), file.name='nodes.tsv', cores=3) {
    db <- .db(dbloc)

    ## Get data from NCBI taxonomy dump
    ncbi.nodes <- getnodes(taxdir)
    ncbi.names <- getnames(taxdir)

    ## do calculation in parallel; pick roughly a number of taxon ids
    ## equal to twice the number of CPUs available, to make
    ## sure we won't have unused CPUs
    ids.to.process <- root.taxa
    ids.not.processed <- vector()

    ## CAUTION: loop below could get stuck if there are not enough echildren and many cores!
    while (length(ids.to.process) < (cores*2)) {
        currentid <- head(ids.to.process, 1)

        ## do not go further down if we are at genus or lower, because we would't get the ti_genus of the children!
        if (getrank(currentid, NULL, nodes=ncbi.nodes) %in% c('genus', 'species', 'subspecies', 'varietas', 'forma'))
            break

        cat("CurrentID : ", currentid, "\n")
        ids.to.process <- tail(ids.to.process, length(ids.to.process)-1)

        for (ch in children(currentid, ncbi.nodes)) {
            ids.to.process <- c(ids.to.process, ch)
        }
        ids.not.processed <- c(ids.not.processed, currentid)
    }

    ##    for (tax in ids.to.process) {
    foreach (i=1:length(ids.to.process)) %dopar% {
        tax <- ids.to.process[i]
        cat("Processing taxid ", tax, " # ", i, " / ", length(ids.to.process), "\n")
        n <- .rec.add.stats(tax, ncbi.nodes)
        ## only return rows for which we collected stats, all others have NA in num.seqs.node
        current.nodes <- n[which(! is.na(n$num.seqs.node)),]

        ## Prepare result for writing to file


        ## add other column values
        taxon.name <- as.character(ncbi.names[match(current.nodes$id, ncbi.names$id),"name"])
        commons <- ncbi.names[grep ('common', ncbi.names$type),]
        common.name <- as.character(commons$name[match(current.nodes$id, commons$id)])
        ## rank_flag is 0 if rank higher than genus
        ##rank.flag=ifelse(current.nodes$rank %in% c( 'genus', 'subgenus', 'species group',
        ##                                        'species subgroup', 'species', 'subspecies',
        ##                                        'varietas', 'forma'), 1, 0)

        ## Create dataframe with correct column names that will be
        ## inserted into the database. Some fields will need additional
        ## information and will be filled later
        nodes.df <- data.frame(ti=as.integer(current.nodes$id),
                               ti_anc=as.integer(current.nodes$parent),
                               terminal_flag=NA,
                               rank_flag=as.integer(! current.nodes$rank=='no rank'),
                               model=NA,
                               taxon_name=taxon.name,
                               common_name=common.name,
                               rank=as.character(current.nodes$rank),
                               n_gi_node=as.integer(current.nodes$num.seqs.node),
                               n_gi_sub_nonmodel=as.integer(current.nodes$num.seqs.subtree.nonmodel),
                               n_gi_sub_model=as.integer(current.nodes$num.seqs.subtree.model),
                               n_clust_node=NA,
                               n_clust_sub=NA,
                               n_PIclust_sub=NA,
                               n_sp_desc=as.integer(current.nodes$num.spec.desc),
                               n_sp_model=as.integer(current.nodes$num.spec.model),
                               n_leaf_desc=as.integer(current.nodes$num.leaf.desc),
                               n_otu_desc=as.integer(current.nodes$num.otu.desc),
                               ti_genus=as.integer(current.nodes$ti.genus),
                               n_genera=as.integer(current.nodes$num.genera)
                               )
        write.table(nodes.df, file=file.name, sep="\t", row.names=F, append=file.exists(file.name), col.names=!file.exists(file.name))
        cat("Wrote ", nrow(nodes.df), " entries to file\n")
        cat("Finished processing taxid ", tax, " # ", i, " / ", length(ids.to.process), "\n")
    }
}

## TODO: make functionality to insert in DB
##    repeat {
##        ret <- try(dbWriteTable(conn=db, name='nodes', value=nodes.df, row.names=F, overwrite=overwrite, append=append))
##        if(!is(ret, "try-error")) break
##    }
##    dbDisconnect(db)
##    rm(db)
##    return(ret)


add.stats <- function(taxids, nodes) {
    result <- data.frame()

    ## do calculation in parallel; pick roughly a number of nodes
    ## equal to twice the number of CPUs available, to make
    ## sure we won't have unused CPUs
    ids.to.process <- taxids
    ids.not.processed <- vector()
    cores <- 3

    ## CAUTION: loop below could get stuck if there are not enough echildren and many cores!
    while (length(ids.to.process) < (cores*2)) {
        currentid <- head(ids.to.process, 1)
        if (length(currentid)<1) break;
        cat("CurrentID : ", currentid, "\n")
        ids.to.process <- tail(ids.to.process, length(ids.to.process)-1)

        for (ch in children(currentid, nodes)) {
            ids.to.process <- c(ids.to.process, ch)
        }
        ids.not.processed <- c(ids.not.processed, currentid)

    }

    ##    for (tax in ids.to.process) {
    result <- foreach (i=1:length(ids.to.process), .combine=rbind) %dopar% {
        tax <- ids.to.process[i]
        cat("Processing taxid ", tax, " # ", i, " / ", length(ids.to.process), "\n")
        n <- .rec.add.stats(tax, nodes)
        cat("Finished processing taxid ", tax, " # ", i, " / ", length(ids.to.process), "\n")
        ## only return rows for which we collected stats, all others have NA in num.seqs.node
        n[which(! is.na(n$num.seqs.node)),]
    }
    recover()
    return(result)
}

## Recursive function to add number of sequences for each taxon to table as produced by getnodes().
## We set sequence counts for 'node' and 'subtree', if a node is of rank species or lower. For subtree,
## this includes all desendants. We also distinguish between model (>=20000 seqs per node) and nonmodel organisms.
.rec.add.stats <- function(taxid, nodes) {

    stats <- c(num.seqs.node=0,
               num.seqs.subtree.model=0,
               num.seqs.subtree.nonmodel=0,
               num.leaf.desc=0, ## also counts itsself, for a leaf it is 1
               num.spec.desc=0, ## also counts itsself, for a species it is 1
               num.spec.model=0,
               ti.genus=nodes[which(nodes$id==taxid), 'ti.genus'],
               num.genera=0,
               num.otu.desc=1)

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

    ## call function recursively for children
    for (child in ch) {

        ## ti.genus has to be passed down the tree and not upwards...
        if (! is.null(stats['ti.genus'])) {
            nodes[which(nodes$id==child),'ti.genus'] <- stats['ti.genus']
        }
        nodes <- .rec.add.stats(child, nodes)
        ## add sequence counts for higher level nodes
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

