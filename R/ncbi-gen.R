## require('rentrez')
require('RCurl')

## get sequences for given GIs
## returns a list named by GIs containing the sequences
seqs.for.gis <- function( gis, local=FALSE, dir=NULL, omit.defline=TRUE ) {

    seqs <- list()
    if ( local ) {
        if ( is.null(dir) ) {
            stop( "Need 'dir' argument if 'local=TRUE'" )
        }
        else {
            seqs <- .local.seqs.for.gis( gis, dir, omit.defline=omit.defline )
        }
    }
    else {
        seqs <- .remote.seqs.for.gis( gis, omit.defline ) 
    }    
    
    if ( omit.defline ) {
        seqs <- gsub("^>.*?\\n|\\n", "", seqs)        
    }

    return(seqs)
}

## retrieve seqs for GIs from ncbi server
.remote.seqs.for.gis <- function( gis, omit.defline=TRUE ) {
    ## get FASTA records
    ncbi.url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%d&rettype=fasta&retmode=text"
    counter <- 0
    num.gis <- length(gis)
    seqs <- lapply ( gis, function(gi) {                        
        counter <<- counter + 1        
        msg <- "Retrieving sequence for gi %d ( # %d of %d, %.2f %% finished )"
        cat(sprintf(msg, gi, counter, num.gis, counter*100/length(gis)), "\n")
        url <- sprintf(ncbi.url, gi)
        getURL(url)       
    } )
    names(seqs) <- gis

    return (seqs)
}

## get GIs using a local copy of Genbank
.local.seqs.for.gis <- function( gis, dir, omit.defline=TRUE ) {
    
}
