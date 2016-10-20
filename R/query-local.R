## Query the 'nodes' table to get number of sequences for taxid from loacl databases

.local.num.seqs.for.taxid <- function(taxid) {
    db <- .db()
    query <- paste('select n_gi_sub_nonmodel, n_gi_sub_model from nodes where ti =', taxid)
    l <- dbGetQuery(db, query)
    return(sum(l))
}
