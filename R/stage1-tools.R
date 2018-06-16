#' @name tax_download
#' @title Download taxonomic records
#' @description Downloads one batch of taxonomic
#' records.
#' @return list of list
#' @param ids Vector of taxonomic IDs
#' @template ps
#' @family run-private
tax_download <- function(ids, ps) {
  recs <- vector('list', length = length(ids))
  names(recs) <- ids
  args <- list(db = "taxonomy", id = ids, rettype = 'xml')
  raw_recs <- search_and_cache(func = rentrez::entrez_fetch,
                               fnm = 'fetch', args = args, ps = ps)
  raw_recs <- XML::xmlToList(raw_recs)
  for (i in seq_along(raw_recs)) {
    raw_rec <- raw_recs[[i]]
    rw_lng <- raw_rec[['LineageEx']]
    lng_ids <- vapply(rw_lng, function(x) x[['TaxId']], '')
    names(lng_ids) <- NULL
    lng_ids <- c(lng_ids, raw_rec[['TaxId']])
    lng_rnks <- vapply(rw_lng, function(x) x[['Rank']], '')
    names(lng_rnks) <- NULL
    lng_rnks <- c(lng_rnks, raw_rec[['Rank']])
    lng <- list('ids' = lng_ids, 'rnks' = lng_rnks)
    cmnms <- raw_rec[['OtherNames']]
    if ('GenbankCommonName' %in% names(cmnms)) {
      cmnm <- cmnms[['GenbankCommonName']]
    } else if ('CommonName' %in% names(cmnms)) {
      cmnm <- cmnms[['CommonName']]
    } else {
      cmnm <- ''
    }
    rec <- new('TaxRec', id = raw_rec[['TaxId']],
                scnm = raw_rec[['ScientificName']],
                cmnm = cmnm, rnk = raw_rec[['Rank']],
                lng = lng, prnt = raw_rec[['ParentTaxId']])
    recs[[rec@id]] <- rec
  }
  recs
}
