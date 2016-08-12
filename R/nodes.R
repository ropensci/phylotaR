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

## add column with number of seqs for subtree or terminals to a nodes table
add.num.seqs <- function(root.taxon, nodes) {
    return(.num.seqs.for.descendants2(root.taxon, nodes))
}

.num.seqs.for.descendants2 <- function(taxid, nodes) {

    num.seqs.node <- 0
    num.seqs.subtree <- 0
    
    ##anchor
    if ( getrank(taxid, NULL, nodes=nodes) %in% c('species', 'subspecies') ) {
        curr.num.seqs <- num.seqs.for.taxid(taxid, nodes)
        num.seqs.node <- curr.num.seqs
        num.seqs.subtree <- num.seqs.subtree + num.seqs.node
    }

    if ( getrank (taxid, NULL, nodes=nodes) == 'subspecies') {
        num.seqs.subtree <- 0
    }
    
    for (child in children(taxid, nodes)) {
        nodes <- .num.seqs.for.descendants2(child, nodes)

        ##if (! getrank(taxid, NULL, nodes=nodes) %in% c('species', 'subspecies') ) {
        num.seqs.subtree <- num.seqs.subtree + nodes[which(nodes$id==child),'num.seqs.subtree']
        ##}
    }    

    cat("Setting for taxid ", taxid, "  ", num.seqs.subtree, "\n")
    nodes[which(nodes$id==taxid),'num.seqs.node'] <- num.seqs.node
    nodes[which(nodes$id==taxid),'num.seqs.subtree'] <- num.seqs.subtree
    return(nodes)
}

## recursive function to add number of sequences for each taxon to table as produced by getnodes()
.num.seqs.for.descendants <- function(taxid, nodes) {
    total <- nrow(nodes)
    num.seqs <- 0
    num.seqs.node <- 0
    num.seqs.subtree <- 0
    num.seqs.subtree.model <- 0
    num.seqs.subtree.nonmodel <- 0
    cat("Taxid : ", taxid, "\n")
    ## ANCHOR OF RECURSION ... set sequence counts for the below levels
    if (getrank(taxid, NULL, nodes=nodes) %in% c('species', 'subspecies','varietas', 'forma')){
        curr.num.seqs <- num.seqs.for.taxid(taxid, nodes)
        num.seqs <- num.seqs + curr.num.seqs

        num.seqs.node <- curr.num.seqs
        num.seqs.subtree <- num.seqs.subtree + curr.num.seqs
        cat("Num seqs subtree : ", num.seqs.subtree , "\n")
        if (curr.num.seqs > 20000) {
            num.seqs.subtree.model <- num.seqs.subtree.model + curr.num.seqs
        } else {
            num.seqs.subtree.nonmodel <- num.seqs.subtree.nonmodel + curr.num.seqs
        }
    }

    for (child in children(taxid, nodes)) {
        ## 1 is its own parent, therefore omit
        if (child==1) {
            next
        }
        l <- .num.seqs.for.descendants(child, nodes)
        nodes <- l$nodes
                
        ## If children are (possibly) terminal nodes, do not add seq count again, because this was done above!
        if (! getrank(taxid, NULL, nodes=nodes) %in% c('species', 'subspecies','varietas', 'forma')){
            num.seqs <- num.seqs + l$num.seqs        
            num.seqs.subtree <- num.seqs.subtree + nodes[which(nodes$id==child),'num.seqs.subtree']
            cat("Now num seqs subtree : ", num.seqs.subtree, "\n")
            num.seqs.subtree.model <- num.seqs.subtree.model + nodes[which(nodes$id==child),'num.seqs.subtree.model']    
            num.seqs.subtree.nonmodel <- num.seqs.subtree.nonmodel + nodes[which(nodes$id==child),'num.seqs.subtree.nonmodel']                
        }        


##        cat("Child : ", child, " model : ", num.seqs.subtree.model, " nonmodel: ", num.seqs.subtree.nonmodel, "\n")

    }
    nodes[which(nodes$id==taxid),'num.seqs'] <- num.seqs
    nodes[which(nodes$id==taxid),'num.seqs.node'] <- num.seqs.node
    nodes[which(nodes$id==taxid),'num.seqs.subtree'] <- num.seqs.subtree + num.seqs.node
    nodes[which(nodes$id==taxid),'num.seqs.subtree.model'] <- num.seqs.subtree.model
    nodes[which(nodes$id==taxid),'num.seqs.subtree.nonmodel'] <- num.seqs.subtree.nonmodel
    ret <- list(nodes=nodes, num.seqs=num.seqs)
    return(ret)
}

num.seqs.for.taxid <- function(taxid, nodes) {
    db <- .db()
    ids <- c(taxid, descendants(taxid, nodes))
    idstr <- paste0(ids, collapse=',')
    str <- paste('select count(*) from accession2taxid where taxid in (', idstr, ')')
    l <- dbGetQuery(db, str)
    return(l[[1]])
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


