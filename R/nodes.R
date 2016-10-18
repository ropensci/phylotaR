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
nodes.create <- function(dbloc, taxdir, root.taxa = c(33090, 4751, 33208), overwrite = FALSE, append = TRUE) {
    db <- .db(dbloc)

    ## Get data from NCBI taxonomy dump
    ncbi.nodes <- getnodes(taxdir)
    ncbi.names <- getnames(taxdir)

    ## add sequence counts etc
    nodes <- add.stats(root.taxa, ncbi.nodes)

    ## add other column values
    taxon.name <- as.character(ncbi.names[match(nodes$id, ncbi.names$id),"name"])
    commons <- ncbi.names[grep ('common', ncbi.names$type),]
    common.name <- as.character(commons$name[match(nodes$id, commons$id)])
    ## rank_flag is 0 if rank higher than genus
    ##rank.flag=ifelse(nodes$rank %in% c( 'genus', 'subgenus', 'species group',
    ##                                        'species subgroup', 'species', 'subspecies',
    ##                                        'varietas', 'forma'), 1, 0)

    ## Create dataframe with correct column names that will be
    ## inserted into the database. Some fields will need additional
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
    ret <- FALSE
    write.table(nodes.df, file='nodes.tsv', sep="\t", row.names=F)
    repeat {
    	   ret <- try(dbWriteTable(conn=db, name='nodes', value=nodes.df, row.names=F, overwrite=overwrite, append=append))
	   if(!is(ret, "try-error")) break
    }
    dbDisconnect(db)
    rm(db)
    return(ret)
}

add.stats <- function(taxids, nodes) {
    result <- data.frame()

    for (tax in taxids) {
        n <- .rec.add.stats(tax, nodes)
        ## only add rows for which we collected stats
        result <- rbind(result, n[which(! is.na(n$num.seqs.node)),])
    }
    return(result)
}

## Recursive function to add number of sequences for each taxon to table as produced by getnodes().
## We set sequences for 'node' and 'subtree', a node is of rank species or lower, subtree includes all
## desendants. We also distinguish between model (>=20000 seqs per node) and nonmodel organisms.
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

#num.seqs.for.taxid <- function(taxid, nodes) {
#    db <- .db()
#    ids <- c(taxid, descendants(taxid, nodes))
#    idstr <- paste0(ids, collapse=',')
#    str <- paste('select count(*) from accession2taxid where taxid in (', idstr, ')')
#    l <- dbGetQuery(db, str)
#    return(l[[1]])
#}

.num.seqs.for.taxid <- function(taxid) {
    search <- entrez_search(db='nucleotide', term=paste0('txid', taxid, '[Organism:exp]', '1:25000[SLEN]'), use_history=T, retmax=1)
    return(search$count)
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

