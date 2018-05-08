# GENERATE TEST DATA

# FUNCTIONS
# TODO: use internal phylotaR functions
dwnldTxRcrds <- function(txids) {
  rcrds <- vector('list', length = length(txids))
  rw_rcrds <- rentrez::entrez_fetch(db = "taxonomy", id = txids,
                                    rettype = 'xml')
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
    rcrd <- new('RecordTaxon', id = rw_rcrd[['TaxId']],
                scnm = rw_rcrd[['ScientificName']],
                cmnm = cmnm, rnk = rw_rcrd[['Rank']],
                lng = lng, prnt = rw_rcrd[['ParentTaxId']])
    rcrds[[i]] <- rcrd
  }
  rcrds
}

# '6962' '1445963' '4613'
demo_txids <- c('7900', '71245', '9504',
                '42242', '9479', '8802')
test_data <- list()
for (demo_txid in demo_txids) {
  cat('... ', demo_txid, '\n', sep = '')
  term <- paste0('txid', demo_txid, '[Subtree]')
  search_results <- rentrez::entrez_search(db = 'taxonomy',
                                           term = term,
                                           retmax = 10000)
  txids <- search_results[['ids']]
  records <- dwnldTxRcrds(txids)
  names(records) <- sapply(records, function(x) x@id)
  test_data[[demo_txid]] <- list('txids' = names(records), 'records' = records)
}
saveRDS(object = test_data, file = file.path('tests', 'testthat', 'data',
                                             'taxonomy', 'txids_records.RData'))
