library("RSQLite")
require('CHNOSZ')
require('rentrez')
require('data.table')

taxdir <<- '~/ncbi-taxonomy/ftp.ncbi.nih.gov/pub/taxonomy'
nodes <<- getnodes(taxdir)
names <<- getnames(taxdir)

##SQLITE3commands
##.mode tabs
##.import nucl_gb.accession2taxid accession2taxid
##create index taxid on accession2taxid(taxid);

con <<- dbConnect(RSQLite::SQLite(),'~/ncbi-taxonomy/ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/mytest')


ncbi.taxonomic.ranks <<- c('superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum',
					   'subphylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder',
					   'order', 'suborder', 'infraorder', 'parvorder', 'superfamily', 'family',
					   'subfamily', 'tribe', 'subtribe', 'genus', 'subgenus', 'species group',
					   'species subgroup', 'species', 'subspecies','varietas', 'forma')

###  DEPRECATED, fiels will be discontinued after Sep 2016
## Will get all GIss for a (vector of) taxon ID('s) from
## the NCBI file gi_taxid_nucl.dmp, the dir it is located
## can be specified by argument 'dir'
## Returns list named by taxids containing the gis
## Will get *all* GIS (gb, gss, wgs, est, etc..)
gis.all.for.taxid <- function( taxid, dir='.' ) {
    gis <- list()
    
    file.name <- 'gi_taxid_nucl.dmp'
    path <- file.path(dir, file.name )

    for ( ti in taxid ) {
        ## cannot read file at once, need to filter for taxid    
        cmd <- paste0("awk '$2 == ", ti, "' ", path)
        tab <- fread(cmd)
        colnames(tab) = c('gi', 'ti')        
        curr.gis <- tab$gi
        gis[ti] <- curr.gis
    }    
    return( gis )
}

## For given taxon id(s), returns a list with accessions and GIs for
## each taxon id. Will only deliver accessions and GIs that are in GenBank
## Returns list of lists for each taxid, list containing accessions and
## GI codes
## When dir=NULL (the default), the NCBI server is queried.
## If the 'local' flag is set and a directory containing the NCBI taxonomy flatfiles
## (available at) is given
## as argument 'dir', the data is parsed locally
acc.genbank.seqids.for.taxid <- function( taxids, dir=NULL, local=FALSE ) {    
    if ( local ) {
        if ( is.null(dir) ) {
            stop( "Need 'dir' argument if 'local==TRUE'" )
        }
        else {
            return( .local.acc.genbank.seqids.for.taxid( taxids, dir ) )
        }
    }
    else {
        return( .remote.acc.genbank.seqids.for.taxid( taxids ) )
    }
}

.local.acc.genbank.seqids.for.taxid <- function( taxids, dir ) {
    ret <- list()
    
    file.name <- 'nucl_gb.accession2taxid'
    subdir <- 'accession2taxid'
    path <- file.path(dir, subdir, file.name)
    
    for ( ti in taxids ) {
        cat("Collecting Genbank accessions for ti", ti, "from local data\n")

        ## need to expand taxids to all lower levels in order to find
        ## seq Ids in the file
        taxa.table <- get.children(ti, depth=1e6, dir)
        
        ## keep only ids for species and lower
        taxranks <- ncbi.taxonomic.ranks[which(ncbi.taxonomic.ranks=='species'):length(ncbi.taxonomic.ranks)]
        expanded.ids <- taxa.table[which(taxa.table[,'rank'] %in% taxranks),'id']
        
        ## cannot read file at once, need to filter for taxid,
        ## get all rows from, file which contain the taxid in the third column
        cmd <- paste0("awk '$3 ~ /", paste0('^', expanded.ids, '$', collapse='|'), "/' ", path)
        tab <- try(fread(cmd))
        if (class(tab) == 'try-error') {
            cat("Could not retrieve sequence IDs for taxon", ti, ", probably there aren't any. Skipping.\n")
        }
        else {
            cat("found", nrow(tab), "accessions\n")
            colnames(tab) <- c('accession', 'accession.version', 'taxid', 'gi') 
            l <- list()
            l[['accession.version']] <- tab$accession.version
            l[['gi']] <- tab$gi
            ret[[as.character(ti)]] <- l
        }
    }
    return (ret)    
}

.remote.acc.genbank.seqids.for.taxid <- function( taxids ) {
    ret <- list()

    for ( ti in taxids ) {        
        cat("Collecting Genbank accessions for ti", ti, "from ncbi server\n")
        search <- entrez_search(db='nucleotide', term=paste0('txid', ti, '[Organism:exp]'), retmax=999999999)
        gis <- search$ids
        cat("found", length(gis), "accessions\n")
        l <- list()
        l[['gi']] <- gis        
        ## query for accessions
        accessions <- sapply( gis, function(gi) {
            s <- entrez_summary(db='nucleotide', id=gi)
            s$acc
        })
        l[['accession.version']] <- accessions
        ret[[as.character(ti)]] <- l        
    }
    return(ret)    
}


get.children <- function(taxa, depth=1e6, dir=NULL, local=FALSE) {
    if (local && is.null(dir)) {
        stop( "Need 'dir' argument if 'local==TRUE'" )
    }
    ret <- do.call(rbind, lapply(taxa, .get.children, depth, dir))
    if (nrow(ret)) {
        rownames(ret) <- 1:nrow(ret)
    }
    return(ret)
}

.get.children <- function(taxon, depth, dir) {    
    queue <- c(taxon)
    result <- data.frame(id=integer(0), rank=character(0))
    while ( length(queue) > 0  & depth > 0) {
#        cat("Depth : ", depth, "\n")
        depth <- depth - 1
        curr.taxon <- queue[1]
        
        ## remove first element from queue
        queue <- tail(queue, length(queue)-1)                
        
        ## get children
        ch <- data.frame()
        if (is.null(dir)) {
            ch <- .remote.children(curr.taxon)
        } else {
            ch <- .local.children(curr.taxon, dir)
        }
        if (nrow(ch)){
            result <- rbind(result, ch)
            ## id 1 is it's own parent, therefore have to remove from queue
            queue <- c(queue, ch[which(ch[,'id']!=1), 'id'])        
        }
    }
    return (result)
}

.remote.children <- function(tax) {
    ret <- data.frame(id=numeric(), rank=character())
    db <- 'ncbi'   
    ch <- taxize::children(tax, db)
    ## Check for unsuccesful query
    if (! is.na (ch)) {
        df <- ch[[1]]
        if (nrow(df)) {
            ret <- data.frame(id=df[,'childtaxa_id'], rank=df[,'childtaxa_rank'])
        }
    }
    return(ret)
}

.local.children <- function(tax, dir) {
    ## buffer nodes object if possible
    if (! exists('nodes')) {
        nodes <<- CHNOSZ::getnodes(dir)
    }
    ch <- nodes[which(nodes['parent']==tax),c('id', 'rank')]
    return(ch[,c('id', 'rank')])
}

nodes.max.seq <- function(root.taxon, max.seqs=20000, dir) {
    ## get descendants of root.taxon which have < max.seq sequences
    ## if a node has > max.seqs sequences, further traverse down the tree
    result <- data.frame(id=numeric(), name=character(), rank=character(), numseq=numeric())

    queue <- c(root.taxon)    
    while( length(queue) > 0 ) {
        cat("Queue : ", paste(queue, sep=','), "\n")
        curr.taxon <- queue[1]        
        cat("Current taxon : ", curr.taxon, "\n")
        queue <- tail(queue, length(queue)-1)                      
        curr.nof.seqs <- num.seqs.for.taxid(curr.taxon)
        cat("Number of seqs : ", curr.nof.seqs, "\n")
        if (curr.nof.seqs > max.seqs) {
            cat("Number of seqs larger than max.seqs, getting children\n")
            taxa.table <- get.children(curr.taxon, depth=1, dir)
            cat("Number of children : ", nrow(taxa.table), "\n")
            queue <- c(queue, taxa.table$id)
        }
        else {
            rank <- getrank(curr.taxon, dir, nodes)
            name <- sciname(curr.taxon, dir, names)
            cat("Found manageable node for taxon ", name, "\n")
            result <- rbind(result, c(curr.taxon, name, rank, curr.nof.seqs ))
            
        }        
    }
    return (result)        
}

num.seqs.for.taxid <- function(taxid) {
    search <- entrez_search(db='nucleotide', term=paste0('txid', taxid, '[Organism:exp]'), retmax=999999999)
    return (search$count)
}

num.seqs.for.taxid2 <- function(taxid, taxdir = NULL, nodes = NULL) {
    length(gis.for.taxid(taxid, taxdir, nodes))
}

num.seqs.for.taxid3 <- function(taxid) {
    ids <- c(taxid, allchildren(taxid, taxdir, nodes))
    idstr <- paste0(ids, collapse=',')
    str <- paste('select count(*) from accession2taxid where taxid in (', idstr, ')')
    l <- dbGetQuery(con, str)
    return(l[[1]])
}

gis.for.taxid <- function(id, taxdir = NULL, nodes = NULL) {
    ch <- c(id, allchildren(id, taxdir, nodes))
    ##ranks <- sapply(children, getrank, taxdir, nodes)
    ## filter for 
    ##valid.ranks <- ncbi.taxonomic.ranks[which(ncbi.taxonomic.ranks=='species'):length(ncbi.taxonomic.ranks)]
    file.name <- 'nucl_gb.accession2taxid'
    subdir <- 'accession2taxid'
    path <- file.path(taxdir, subdir, file.name)

    cmd <- paste0("awk '$3 ~ /", paste0('^', ch, '$', collapse='|'), "/' ", path)
    tab <- try(fread(cmd))
    if (class(tab) == 'try-error') {
        cat("Could not retrieve sequence IDs for taxon", id, ", probably there aren't any. Skipping.\n")
    }
    else {
            cat("found", nrow(tab), "accessions\n")
#            colnames(tab) <- c('accession', 'accession.version', 'taxid', 'gi')             
    }
    return (tab$gi)
}
    


children <- function(id, taxdir = NULL, nodes = NULL) {
    if (is.null(nodes)) 
        nodes <- getnodes(taxdir)
    return(nodes$id[which(nodes$parent==id)])
}

allchildren <- function(id, taxdir = NULL, nodes = NULL) {
    if (is.null(nodes)) 
        nodes <- getnodes(taxdir)
    queue <- id
    result <- numeric()
    while (length(queue)>0) {
        currentid <- head(queue, 1)
        queue <- tail(queue, length(queue)-1)
        currentchildren <- children(currentid, taxdir, nodes)
        result <- c(result, currentchildren)
        queue <- c(queue, currentchildren)      
    }
    return(result)
}

