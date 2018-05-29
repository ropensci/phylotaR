#' @name taxise_run
#' @title Run taxise stage
#' @description TODO
#' @param wd Working directory
#' @details Object will be cached.
#' @export
taxise_run <- function(wd) {
  # TODO: allow a user to have their own taxids and/or tax tree
  ps <- ldPrmtrs(wd)
  msg <- paste0('Starting stage TAXISE: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
  info(lvl = 1, ps = ps, 'Searching taxonomic IDs ...')
  txids <- txids_get(ps = ps)
  info(lvl = 1, ps = ps, 'Downloading taxonomic records ...')
  rcrds <- batcher(ids = txids, func = tax_download, ps = ps,
                   lvl = 2)
  info(lvl = 1, ps = ps, 'Generating taxonomic dictionary ...')
  txdct <- taxdict_gen(rcrds = rcrds, txids = txids)
  svObj(wd = wd, obj = txdct, nm = 'txdct')
  msg <- paste0('Completed stage TAXISE: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
}

#' @name txids_get
#' @title Searches for descendent taxonomic IDs
#' @description Searches NCBI taxonomy for all descendent
#' taxonomic nodes.
#' @return Vector of txids
#' @param ps Parameter list
#' @param retmax integer, maximum number of IDs to return per query
txids_get <- function(ps, retmax = 1E4) {
  # TODO: handle multiple txids
  trm <- paste0('txid', ps[['txid']],'[Subtree]')
  args <- list(db = 'taxonomy', term = trm, retmax = retmax)
  srch_rs <- search_and_cache(func = rentrez::entrez_search,
                      args = args, fnm = 'search', ps = ps)
  txcnt <- srch_rs[['count']]
  txids <- srch_rs[['ids']]
  if (txcnt < retmax) {
    return(txids)
  }
  ret_strts <- seq(retmax, txcnt, retmax)
  for (ret_strt in ret_strts) {
    args <- list(db = 'taxonomy', term = trm, retmax = retmax,
                 retstart = ret_strt)
    srch_rs <- search_and_cache(func = rentrez::entrez_search,
                                args = args, fnm = 'search', ps = ps)
    txids <- c(txids, srch_rs[['ids']])
  }
  txids
}

#' @name taxdict_gen
#' @title Generate taxonomic dictionary
#' @description Takes a vector of txids and a list
#' of taxonomic records and returns a taxonomic dictionary.
#' @return TaxDct
#' @param txids Vector of taxonomic IDs
#' @param rcrds List of taxonomic records
taxdict_gen <- function(txids, rcrds) {
  # TODO: allow paraphyly
  # identify pre-node IDs
  # based upon: https://github.com/DomBennett/treeman/wiki/trmn-format
  # drop singletons
  # create index to recover original IDs, indx
  prids <- vapply(txids, function(x) rcrds[[x]]@prnt, '')
  names(prids) <- NULL
  root_bool <- !prids %in% txids
  root <- txids[root_bool]
  prids[root_bool] <- root
  root_trid <- as.character(which(root_bool))
  trids <- seq_along(txids)
  n_ptnds <- table(prids)
  to_drop <- NULL
  while(any(n_ptnds == 1)) {
    # get singletons
    sngltns <- names(n_ptnds)[n_ptnds == 1]
    sngltns_indx <- match(sngltns, txids)
    sngltns_trids <- match(sngltns, prids)
    # replace singleton prid in tree IDs
    trids[sngltns_indx] <- sngltns_trids
    # replace singleton prid in prids
    prids[sngltns_trids] <- prids[sngltns_indx]
    to_drop <- c(to_drop, sngltns)
    n_ptnds <- table(prids)
  }
  indx <- trids
  to_drop <- !txids %in% to_drop
  prids <- prids[to_drop]
  trids <- as.character(trids[to_drop])
  prinds <- match(prids, txids[to_drop])
  prinds <- as.integer(prinds)
  # create tax tree
  txtr <- taxonomic_tree_generate(prinds=prinds, trids=trids, root=root_trid)
  # create tax dict
  new('TaxDct', txids=txids, rcrds=list2env(rcrds), txtr=txtr,
      prnt=root, indx=indx)
}
