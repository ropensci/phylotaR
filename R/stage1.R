#' @name taxise_run
#' @title Run taxise stage
#' @description Run the first stage of phylotaR, taxise. This looks up
#' all descendant taxonomic nodes for a given taxonomic ID. It then
#' looks up relevant taxonomic information and generates a taxonomic
#' dictionary for user interaction after phylotaR has completed.
#' @param wd Working directory
#' @details Objects will be cached.
#' @return NULL
#' @family run-public
#' @export
#' @example examples/taxise_run.R
taxise_run <- function(wd) {
  # TODO: allow a user to have their own taxids and/or tax tree
  ps <- parameters_load(wd)
  msg <- paste0('Starting stage TAXISE: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
  info(lvl = 1, ps = ps, 'Searching taxonomic IDs ...')
  txids <- txids_get(ps = ps)
  info(lvl = 1, ps = ps, 'Downloading taxonomic records ...')
  recs <- batcher(ids = txids, func = tax_download, ps = ps, lvl = 2)
  info(lvl = 1, ps = ps, 'Generating taxonomic dictionary ...')
  txdct <- taxdict_gen(recs = recs, txids = txids, ps = ps)
  obj_save(wd = wd, obj = txdct, nm = 'txdct')
  msg <- paste0('Completed stage TAXISE: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
}

#' @name txids_get
#' @title Searches for descendant taxonomic IDs
#' @description Searches NCBI taxonomy for all descendant taxonomic nodes.
#' @return Vector of txids
#' @template ps
#' @param retmax integer, maximum number of IDs to return per query
#' @return vector of ids
#' @family run-private
txids_get <- function(ps, retmax = 1E4) {
  # TODO: handle multiple txids
  trm <- paste0('txid', ps[['txid']],'[Subtree]')
  args <- list(db = 'taxonomy', term = trm, retmax = retmax)
  srch_rs <- search_and_cache(func = rentrez::entrez_search,
                              args = args, fnm = 'search', ps = ps)
  txcnt <- srch_rs[['count']]
  txids <- srch_rs[['ids']]
  if (txcnt <= retmax) {
    return(txids)
  }
  ret_strts <- seq(from = retmax, to = txcnt, by = retmax)
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
#' @param txids Vector of taxonomic IDs
#' @param recs List of taxonomic records
#' @template ps
#' @family run-private
#' @return TaxDict
taxdict_gen <- function(txids, recs, ps) {
  # TODO: allow paraphyly
  # identify pre-node IDs
  # based upon: https://github.com/DomBennett/treeman/wiki/trmn-format
  # create index to recover original IDs, indx
  prids <- vapply(txids, function(x) recs[[x]]@prnt, '')
  names(prids) <- NULL
  root_bool <- !prids %in% txids
  root <- txids[root_bool]
  prids[root_bool] <- root
  prinds <- match(prids, txids)
  prinds <- as.integer(prinds)
  # create tax tree
  txtr <- taxtree_gen(prinds = prinds, ids = txids, root = root,
                      ps = ps)
  # create tax dict
  new('TaxDict', txids = txids, recs = list2env(recs), txtr = txtr,
      prnt = root)
}
