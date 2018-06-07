#' @name tax_download
#' @title Download taxonomic records
#' @description Downloads one batch of taxonomic
#' records.
#' @return list of list
#' @param ids Vector of taxonomic IDs
#' @param ps Parameter list
#' @noRd
tax_download <- function(ids, ps) {
  rcrds <- vector('list', length = length(ids))
  names(rcrds) <- ids
  args <- list(db = "taxonomy", id = ids, rettype = 'xml')
  rw_rcrds <- search_and_cache(func = rentrez::entrez_fetch,
                               fnm = 'fetch', args = args, ps = ps)
  rw_rcrds <- XML::xmlToList(rw_rcrds)
  for (i in seq_along(rw_rcrds)) {
    rw_rcrd <- rw_rcrds[[i]]
    rw_lng <- rw_rcrd[['LineageEx']]
    lng_ids <- vapply(rw_lng, function(x) x[['TaxId']], '')
    names(lng_ids) <- NULL
    lng_ids <- c(lng_ids, rw_rcrd[['TaxId']])
    lng_rnks <- vapply(rw_lng, function(x) x[['Rank']], '')
    names(lng_rnks) <- NULL
    lng_rnks <- c(lng_rnks, rw_rcrd[['Rank']])
    lng <- list('ids' = lng_ids, 'rnks' = lng_rnks)
    cmnms <- rw_rcrd[['OtherNames']]
    if ('GenbankCommonName' %in% names(cmnms)) {
      cmnm <- cmnms[['GenbankCommonName']]
    } else if ('CommonName' %in% names(cmnms)) {
      cmnm <- cmnms[['CommonName']]
    } else {
      cmnm <- ''
    }
    rcrd <- new('TaxRec', id = rw_rcrd[['TaxId']],
                scnm = rw_rcrd[['ScientificName']],
                cmnm = cmnm, rnk = rw_rcrd[['Rank']],
                lng = lng, prnt = rw_rcrd[['ParentTaxId']])
    rcrds[[rcrd@id]] <- rcrd
  }
  rcrds
}

