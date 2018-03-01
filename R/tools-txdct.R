#' @name getIDFrmTxdct
#' @title Get taxonomic IDs by rank
#' @description
#' @param txdct Taxonomic dictionary
#' @param id Taxon IDs
#' @param ret Return ID or Name?
#' @param rank Rank of output
getIDFrmTxdct <- function(txdct, id, ret=c('TaxId', 'ScientificName'),
                          rank=c("superkingdom", "kingdom", "phylum",
                                 "subphylum", "class", "superorder",
                                 "order", "suborder", "infraorder",
                                 "parvorder", "family", "genus",
                                 "species", "subspecies")) {
  calc <- function(id) {
    rcrd <- txdct[[id]]
    rnks <- sapply(rcrd[['LineageEx']],
                   function(x) x[['Rank']])
    ids <- sapply(rcrd[['LineageEx']],
                  function(x) x[[ret]])
    mtch <- which(rank == rnks)
    if(length(mtch) == 0) {
      # if no match, use lowest available rank
      mtch <- length(rnks)
    }
    ids[[mtch[[1]]]]
  }
  rank <- match.arg(rank)
  ret <- match.arg(ret)
  sapply(as.character(id), calc)
}
