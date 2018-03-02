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