## script to populate the phylota 'nodes' table

library("RSQLite")
require('CHNOSZ')

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
nodes.create <- function(dbloc, taxdir, overwrite = FALSE) {
    db <- .db(dbloc)

    ## Get data from NCBI taxonomy dump
    ncbi.nodes <- getnodes(taxdir)
    ncbi.names <- getnames(taxdir)

    ## Prepare column values
    ## terminal.flag=lapply(ncbi.nodes$id, function(n)length(children(n, ncbi.nodes)==0))
    taxon.name <- ncbi.names[match(ncbi.nodes$id, ncbi.names$id),"name"]
    commons <- ncbi.names[grep ('common', ncbi.names$type),]
    common.name <- commons$name[match(ncbi.nodes$id, commons$id)]
    ## rank_flag is 0 if rank higher than genus
    rank.flag=ifelse(ncbi.nodes$rank %in% c( 'genus', 'subgenus', 'species group',
                                            'species subgroup', 'species', 'subspecies',
                                            'varietas', 'forma'), 1, 0)
    n.gi <- get.n_gi(ncbi.nodes)
    
    ## Create dataframe with correct column names that will be
    ## inserted into the database. Some fields will need additional
    ## information and will be filled later
    nodes.df <- data.frame(ti=ncbi.nodes$id,
                           ti_anc=ncbi.nodes$parent,
                           terminal_flag=NA,
                           rank_flag=rank.flag,
                           model=NA,
                           taxon_name=taxon.name,
                           common_name=common.name,
                           rank=ncbi.nodes$rank,
                           n_gi_node=n.gi,
                           n_gi_sub_nonmodel=NA,
                           n_gi_sub_model=NA,
                           n_clust_node=NA,
                           n_clust_sub=NA,
                           n_PIclust_sub=NA,
                           n_sp_desc=NA,
                           n_sp_model=NA,
                           n_leaf_desc=NA,
                           n_otu_desc=NA,
                           ti_genus=NA,
                           n_genera=NA
                           )
    ## correct data types for integers
    nodes.df <- transform(nodes.df, ti=as.integer(ti), ti_anc=as.integer(ti_anc), rank_flag=as.integer(rank_flag))
    ret <- dbWriteTable(conn=db, name='nodes', value=nodes.df, row.names=F, overwrite=overwrite)
    dbDisconnect(db)
    rm(db)
    return(ret)
}

## Gets the fields n_gi_node,
## Need table accession2taxid for doing so
get.n_gi <- function(nodes) {
    ## get data for kingdoms viridiplantae, fungi, metazoa
    ##kingdoms <- c(33090, 4751, 33208)

    kingdoms <- c(9688, 3655)

    for (k in kingdoms) {        
        nodes <- add.num.seqs(k, nodes)        
    }
    
    ## get n_gi_node
    n.gi.node <- rep(NA, nrow(nodes))    
    n <- which(nodes$rank %in% c('species', 'subspecies','varietas', 'forma'))    
    n.gi.node[n] <- nodes$seqs[n]


    
    return(nodes$seqs)
}

## Recursive function to add number of sequences for each taxon to table as produced by getnodes().
## We set sequences for 'node' and 'subtree', a node is of rank species or lower, subtree includes all
## desendants. We also distinguish between model (>=20000 seqs per node) and nonmodel organisms.
add.num.seqs <- function(taxid, nodes) {
    num.seqs.node <- 0
    num.seqs.subtree <- 0
    num.seqs.subtree.model <- 0
    num.seqs.subtree.nonmodel <- 0    
    num.leaf.desc <- 0 
    num.spec.desc <- 0
    ## number of model species
    num.spec.model <- 0
    ti.genus <- nodes[which(nodes$id==taxid), 'ti.genus']
    n.genera <- 0
    num.otu.desc <- 1

    stats <- c(num.seqs.node=0,
               num.seqs.subtree.model=0,
               num.seqs.subtree.nonmodel=0,
               num.leaf.desc=0, ## also counts itsself, for a leaf it is 1
               num.spec.desc=0, ## also counts itsself, for a species it is 1
               num.spec.model=0,
               ti.genus=nodes[which(nodes$id==taxid), 'ti.genus'],
               n.genra=0,
               num.otu.desc=0)
    
    ## Ranks for which we directly get the species count;
    ## for ranks above this, we will get the sum of their children.
    ## For taxa with these ranks, there will also be a 'n_gi_node'
    node.ranks <- c('species', 'subspecies', 'varietas', 'forma')

    ## Recursion anchor
    rank <- getrank(taxid, NULL, nodes=nodes)
    if (rank %in% node.ranks) {
        curr.num.seqs <- num.seqs.for.taxid(taxid, nodes)
        num.seqs.node <- curr.num.seqs
        stats['num.seqs.node'] <- curr.num.seqs
        num.seqs.subtree <- curr.num.seqs
        stats['num.seqs.subtree'] <- curr.num.seqs
        if (curr.num.seqs >= 20000) {
            num.seqs.subtree.model <- curr.num.seqs
            stats['num.seqs.subtree.model'] <- curr.num.seqs
        } else {
            num.seqs.subtree.nonmodel <- curr.num.seqs
            stats['num.seqs.subtree.nonmodel'] <- curr.num.seqs
        }
        if (rank == 'species') {        
            num.spec.desc <- 1
            stats['num.spec.desc']
            if (curr.num.seqs >= 20000) {
                num.spec.model <- 1
                stats['num.spec.model'] <- 1
            }
        }        
    }

    if (rank == 'genus') {
        ti.genus <- taxid
        stats['ti.genus'] <- taxid
    }
    
    ch <- children(taxid, nodes)
    if (length(ch)==0) {
        num.leaf.desc <- num.leaf.desc + 1
        stats['num.leaf.desc'] <- stats['num.leaf.desc'] + 1
    }
    ## call recursively for children
    for (child in ch) {
        ## ti.genus has to be passed down the tree and not upwards...
        if (! is.null(ti.genus)) {
            nodes[which(nodes$id==child),'ti.genus'] <- stats['ti.genus']
        }        
        nodes <- add.num.seqs(child, nodes)
        if (! rank %in% node.ranks) {
            num.seqs.subtree <- num.seqs.subtree + nodes[which(nodes$id==child),'num.seqs.subtree']
            #stats['num.seqs.subtree'] <- stats['num.seqs.subtree'] + nodes[which(nodes$id==child),'num.seqs.subtree']
            num.seqs.subtree.model <- num.seqs.subtree.model + nodes[which(nodes$id==child),'num.seqs.subtree.model']
            stats['num.seqs.subtree.model'] <- stats['num.seqs.subtree.model'] + nodes[which(nodes$id==child),'num.seqs.subtree.model']
            num.seqs.subtree.nonmodel <- num.seqs.subtree.nonmodel + nodes[which(nodes$id==child),'num.seqs.subtree.nonmodel']
            stats['num.seqs.subtree.nonmodel'] <- stats['num.seqs.subtree.nonmodel'] + nodes[which(nodes$id==child),'num.seqs.subtree.nonmodel']
        }
        if (getrank(child, NULL, nodes) == 'genus') {
            n.genera <- n.genera + 1
            stats['n.genera'] <- stats['n.genera'] + 1
        }
        num.leaf.desc <- num.leaf.desc + nodes[which(nodes$id==child),'num.leaf.desc']
        stats['num.leaf.desc'] <- stats['num.leaf.desc'] + nodes[which(nodes$id==child),'num.leaf.desc']
        num.otu.desc <- num.otu.desc + nodes[which(nodes$id==child),'num.otu.desc']
        stats['num.otu.desc'] <- stats['num.otu.desc'] + nodes[which(nodes$id==child),'num.otu.desc']
        num.spec.desc <- num.spec.desc + nodes[which(nodes$id==child),'num.spec.desc']
        stats['num.spec.desc'] <- stats['num.spec.desc'] + nodes[which(nodes$id==child),'num.spec.desc']
        num.spec.model <- num.spec.model + nodes[which(nodes$id==child),'num.spec.model']
        stats['num.spec.model'] <- stats['num.spec.model'] + nodes[which(nodes$id==child),'num.spec.model']
        n.genera <- n.genera + nodes[which(nodes$id==child),'n.genera']
        stats['n.genera'] <- stats['n.genera'] + nodes[which(nodes$id==child),'n.genera']
    }       
    nodes[which(nodes$id==taxid), names(stats)] <- stats
    ## Set the sequence numbers as columns to the nodes table
#    nodes[which(nodes$id==taxid),'num.seqs.node'] <- num.seqs.node
#    nodes[which(nodes$id==taxid),'num.seqs.subtree'] <- num.seqs.subtree
#    nodes[which(nodes$id==taxid),'num.seqs.subtree.model'] <- num.seqs.subtree.model
#    nodes[which(nodes$id==taxid),'num.seqs.subtree.nonmodel'] <- num.seqs.subtree.nonmodel    

#    nodes[which(nodes$id==taxid),'num.leaf.desc'] <- num.leaf.desc
#    nodes[which(nodes$id==taxid),'num.otu.desc'] <- num.otu.desc
#    nodes[which(nodes$id==taxid),'num.spec.desc'] <- num.spec.desc
#    nodes[which(nodes$id==taxid),'num.spec.model'] <- num.spec.model
#    nodes[which(nodes$id==taxid),'n.genera'] <- n.genera    
#    nodes[which(nodes$id==taxid),'ti.genus'] <- ifelse(is.null(ti.genus, NA, ti.genus))
    
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
library('rentrez')
num.seqs.for.taxid <- function(taxid, nodes) {
    search <- entrez_search(db='nucleotide', term=paste0('txid', taxid, '[Organism:exp]'), retmax=999999999)
    return (search$count)
}

## function to get a singleton global database object
.db <- function(dbloc = NULL) {
    if (exists('db')){
        return(db)
    }
    else {
        if (is.null(dbloc)) {
            stop('Need location of database')
        }
        db <<- dbConnect(RSQLite::SQLite(), dbloc)
    }
    return(db)        
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


