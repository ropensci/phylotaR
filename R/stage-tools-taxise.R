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
    rw_rcrd <- rw_rcrds[[i]]
    rw_lng <- rw_rcrd[['LineageEx']]
    lng_ids <- vapply(rw_lng, function(x) x[['TaxId']], '')
    names(lng_ids) <- NULL
    lng_ids <- c(lng_ids, rw_rcrd[['TaxId']])
    lng_rnks <- vapply(rw_lng, function(x) x[['Rank']], '')
    names(lng_rnks) <- NULL
    lng_rnks <- c(lng_rnks, rw_rcrd[['Rank']])
    lng <- list('ids'=lng_ids, 'rnks'=lng_rnks)
    cmnms <- rw_rcrd[['OtherNames']]
    if('GenbankCommonName' %in% names(cmnms)) {
      cmnm <- cmnms[['GenbankCommonName']]
    } else if ('CommonName' %in% names(cmnms)) {
      cmnm <- cmnms[['CommonName']]
    } else {
      cmnm <- ''
    }
    rcrd <- new('TxRcrd', id=rw_rcrd[['TaxId']],
                scnm=rw_rcrd[['ScientificName']],
                cmnm=cmnm, rnk=rw_rcrd[['Rank']],
                lng=lng, prnt=rw_rcrd[['ParentTaxId']])
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
  btch <- 100
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