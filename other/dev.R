# VARS
id <- 'AJ428944'
# AJ428944, AF434826
mn_lngth <- 250
mx_lngth <- 2500

# RETRIEVE
raw_rcrd <- rentrez::entrez_fetch(db="nucleotide",
                                 rettype='gbwithparts',
                                 retmode='text', id=id)
rcrd <- XML::xmlToList(raw_rcrd)

# KEY INFO
id <- rcrd[['GBSeq']][['GBSeq_primary-accession']]
ml_typ <- rcrd[['GBSeq']][['GBSeq_moltype']]

# FEATURE
fltr_tbl <- rcrd[['GBSeq']][['GBSeq_feature-table']]
fltr <- fltr_tbl[[2]]
ids <- nms <- strts <- stps <- NULL
for(i in seq_along(fltr_tbl)) {
  fltr <- fltr_tbl[[i]]
  print(i)
  if(fltr[['GBFeature_key']] == 'gene') {
    lctn <- fltr[['GBFeature_location']]
    if(!grepl('^[0-9\\.]+$', lctn)) {
      next
    }
    strtstp <- strsplit(x=lctn, split='\\.\\.')[[1]]
    strtstp <- as.numeric(strtstp)
    sqlngth <- strtstp[2] - strtstp[1]
    if(sqlngth > mn_lngth & sqlngth < mx_lngth) {
      nm <- fltr[['GBFeature_quals']][['GBQualifier']][['GBQualifier_value']]
      nms <- c(nms, nm)
      strts <- c(strts, strtstp[1])
      stps <- c(stps, strtstp[2])
      ids <- c(ids, paste0(id, '/', lctn))
    }
  }
}
