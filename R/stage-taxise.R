
#' @name getTxids
#' @title Searches for descendent taxonomic IDs
#' @description Searches NCBI taxonomy for all descendent
#' taxonomic nodes.
#' @return Vector of txids
#' @param ps Parameter list
getTxids <- function(ps) {
  # TODO: handle multiple txids
  # TODO: use safe search
  trm <- paste0('txid', ps[['txid']],'[Subtree]')
  srch_rs <- rentrez::entrez_search(db='taxonomy',
                                    term=trm, retmax=10000)
  srch_rs[['ids']]
}

#' @name dwnldTxRcrds
#' @title Download taxonomic records
#' @description Takes a vector of txids and returns
#' a list of taxonomic records.
#' @return List of taxonomic records
#' @param txids Vector of taxonomic IDs
#' @param ps Parameter list
dwnldTxRcrds <- function(txids, ps) {
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

#' @name genTxDct
#' @title Generate taxonomic dictionary
#' @description Takes a vector of txids and a list
#' of taxonomic records and returns a taxonomic dictionary.
#' @return TxDct
#' @param txids Vector of taxonomic IDs
#' @param rcrds List of taxonomic records
genTxDct <- function(txids, rcrds) {
  # TODO: allow paraphyly
  prids <- vapply(txids, function(x) rcrds[[x]][['ParentTaxId']], '')
  names(prids) <- NULL
  root <- txids[!prids %in% txids]
  prinds <- match(prids, txids)
  prinds[is.na(prinds)] <- which(is.na(prinds))
  prinds <- as.integer(prinds)
  txtr <- genTxTr(prinds=prinds, txids=txids, root=root)
  new('TxDct', txids=txids, rcrds=list2env(rcrds), txtr=txtr,
      prnt=root)
}