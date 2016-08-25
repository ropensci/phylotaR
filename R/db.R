library("RSQLite")

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
