library("RSQLite")

## function to get a singleton global database object
.db <- function(dbloc = NULL) {
    if (exists('db') && dbIsValid(db)){
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


## DB queries
.children <- function(taxid) {
    db <- .db()
    query <- paste('select ti from nodes where ti_anc =', taxid)
    l <- dbGetQuery(db, query)
    return(l[[1]])
}

.genus.for.taxid <- function(taxid) {
    ##cat("Retreiving genus for taxid ", taxid, "\n")
    db <- .db()
    query <- paste('select ti_genus from nodes where ti=', taxid)
    l <- dbGetQuery(db, query)
    return(l[[1]])
}
