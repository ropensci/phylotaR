#' @name prtDwnldTxRcrds
#' @title Download Taxonomic Records
#' @description Downloads one batch of taxonomic
#' records.
#' @return list of list
#' @param txids Vector of taxonomic IDs
#' @param ps Parameter list
prtDwnldTxRcrds <- function(txids, ps) {
  rcrds <- vector('list', length=length(txids))
  args <- list(db="taxonomy", id=txids, rettype='xml')
  rw_rcrds <- safeSrch(func=rentrez::entrez_fetch,
                       fnm='fetch', args=args, ps=ps)
  rw_rcrds <- XML::xmlToList(rw_rcrds)
  for(i in seq_along(rw_rcrds)) {
    rcrd <- rw_rcrds[[i]]
    rcrd[['Lineage']] <- strsplit(rcrd[['Lineage']],
                                  split='; ')[[1]]
    itslf <- list('TaxId'=rcrd[['TaxId']],
                  'ScientificName'=rcrd[['ScientificName']],
                  'Rank'=rcrd[['Rank']])
    rcrd[['LineageEx']] <- c(rcrd[['LineageEx']],
                             list('Taxon'=itslf))
    rcrds[[i]] <- rcrd
  }
  rcrds
}

#' @name btchDwnldTxRcrds
#' @title Download taxonomic records
#' @description Takes a vector of txids and returns
#' a list of taxonomic records.
#' @return List of taxonomic records
#' @param txids Vector of taxonomic IDs
#' @param ps Parameter list
btchDwnldTxRcrds <- function(txids, ps) {
  all_rcrds <- vector('list', length=length(txids))
  names(all_rcrds) <- txids
  btch <- 500
  for(i in seq(0, length(txids)-1, btch)) {
    lower <- i+1
    upper <- ifelse(i+btch<length(txids), i+btch, length(txids))
    crrnt_ids <- txids[lower:upper]
    info(lvl=2, ps=ps, "[", lower, "-", upper, "]");
    rcrds <- prtDwnldTxRcrds(txids=crrnt_ids, ps=ps)
    all_rcrds[lower:upper] <- rcrds
  }
  all_rcrds
}