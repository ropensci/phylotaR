require('CHNOSZ')
require('data.table')

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
acc.genbank.seqids.for.taxid <- function( taxids, local=FALSE, dir=NULL ) {    
    if ( local ) {
        if ( is.null(dir) ) {
            stop( "Need 'dir' argument if 'local=TRUE'" )
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
        cat("Collecting Genbank accessions for ti", ti, "... ")
        ## cannot read file at once, need to filter for taxid            
        cmd <- paste0("awk '$3 == ", ti, "' ", path)
        tab <- fread(cmd)
        cat("found", nrow(tab), "accessions\n")
        colnames(tab) = c('accession', 'accession.version', 'taxid', 'gi') 
        l <- list()
        l[['accession.version']] <- tab$accession.version
        l[['gi']] <- tab$gi
        ret[[as.character(ti)]] <- l        
    }
    return ( ret )    
}

.remote.acc.genbank.seqids.for.taxid <- function( taxids ) {
    ret <- list()

    for ( ti in taxids ) {
        search <- entrez_search(db='nucleotide', term=paste0('txid', ti, '[Organism:exp]'), retmax=999999999)
        gis <- search$ids
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
    return( ret )    
}
