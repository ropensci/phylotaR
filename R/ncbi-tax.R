library("RSQLite")
require('CHNOSZ')
require('rentrez')
require('data.table')

taxdir <<- '~/ftp.ncbi.nih.gov/pub/taxonomy'
if (!exists('ncbi.nodes'))
    ncbi.nodes <<- getnodes(taxdir)
if (!exists('ncbi.names'))
    ncbi.names <<- getnames(taxdir)

##SQLITE3commands
##.mode tabs
##.import nucl_gb.accession2taxid accession2taxid
##create index taxid on accession2taxid(taxid);

con <<- dbConnect(RSQLite::SQLite(),'~/ncbi-taxonomy/ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/mytest')

get.cluster.nodes <- function(root.taxon, max.seqs=20000) {
    ## get descendants of root.taxon which have < max.seq sequences
    ## if a node has > max.seqs sequences, further traverse down the tree
    result <- as.data.frame(matrix(ncol=4, nrow=0))
    queue <- c(root.taxon)    
    while( length(queue) > 0 ) {
        cat("Queue : ", paste(queue, sep=','), "\n")
        curr.taxon <- queue[1]        
        cat("Current taxon : ", curr.taxon, "\n")
        queue <- tail(queue, length(queue)-1)                      
        curr.nof.seqs <- num.seqs.for.taxid(curr.taxon)
        cat("Number of seqs : ", curr.nof.seqs, "\n")        
        ch <- children(curr.taxon)
        rank <- getrank(curr.taxon, taxdir, ncbi.nodes)
        if (curr.nof.seqs > max.seqs && length(ch) > 0  && ! rank %in% c('species', 'subspecies','varietas', 'forma')) {
            cat("Number of seqs larger than max.seqs, getting children\n")

            ##cat("Number of children : ", length(ch), "\n")
            queue <- c(queue, ch)
        }
        else {
            name <- sciname(curr.taxon, taxdir, ncbi.names)
            cat("Found manageable node for taxon ", name, "\n")
            result[nrow(result)+1,] <- c(curr.taxon, name, rank, curr.nof.seqs)
        }        
    }
    colnames(result) <- c("taxid", "name", "rank", "seqs")
    return (result)        
}

num.seqs.for.taxid <- function(taxid) {
    ## taxid could be higher level, therefore include all descendants
    ids <- c(taxid, descendants(taxid))
    idstr <- paste0(ids, collapse=',')
    str <- paste('select count(*) from accession2taxid where taxid in (', idstr, ')')
    l <- dbGetQuery(con, str)
    return(l[[1]])
}

gis.for.taxid <- function(taxid) {
    ## taxid could be higher level, therefore include all descendants
    ch <- c(id, descendants(taxid))
    idstr <- paste0(ch, collapse=',')
    str <- paste('select gi from accession2taxid where taxid in (', idstr, ')')
    l <- dbGetQuery(con, str)
    return(l[[1]])        
}
    
children <- function(id, taxdir) {
    return(ncbi.nodes$id[which(ncbi.nodes$parent==id)])
}

descendants <- function(id) {
    queue <- id
    result <- numeric()
    while (length(queue)>0) {
        currentid <- head(queue, 1)
        queue <- tail(queue, length(queue)-1)
        currentchildren <- children(currentid)
        result <- c(result, currentchildren)
        queue <- c(queue, currentchildren)      
    }
    return(result)
}

