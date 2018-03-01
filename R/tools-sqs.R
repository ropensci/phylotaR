#' @name rwRcrd2SqRcrd
#' @title Convert raw Entrez gbwithparts record to SqRcrds
#' @description Parses returned sequences features with Entrez,
#' returns a SqRcrd for each raw record.
#' @param rw_rcrds Raw records returned from Entrez fetch
#' @param gis GIs of each fectched record
#' @param ps Parameters list
rwRcrd2SqRcrd <- function(rw_rcrds, gis, ps) {
  res <- NULL
  rcrds <- XML::xmlToList(rw_rcrds)
  for(i in seq_along(rcrds)) {
    rcrd <- rcrds[[i]]
    wftrs <-  FALSE
    # key info
    accssn <- rcrd[['GBSeq_primary-accession']]
    vrsn <- rcrd[["GBSeq_accession-version"]]
    ml_typ <- rcrd[['GBSeq_moltype']]
    sq <- rcrd[['GBSeq_sequence']]
    if(is.null(sq)) {
      # master records have no sequences
      # e.g. https://www.ncbi.nlm.nih.gov/nuccore/1283191328
      next
    }
    dfln <- rcrd[['GBSeq_definition']]
    orgnsm <- rcrd[['GBSeq_organism']]
    ftr_tbl <- rcrd[['GBSeq_feature-table']]
    kywrds <- vapply(rcrd[['GBSeq_keywords']], '[', '')
    kywrds <- paste0(kywrds, collapse=' | ')
    # extract features
    for(ftr in ftr_tbl) {
      if(ftr[['GBFeature_key']] == 'gene') {
        lctn <- ftr[['GBFeature_location']]
        if(!grepl('^[0-9\\.]+$', lctn)) {
          next
        }
        strtstp <- strsplit(x=lctn, split='\\.\\.')[[1]]
        strtstp <- as.numeric(strtstp)
        sqlngth <- strtstp[2] - strtstp[1]
        if(sqlngth > ps[['mnsql']] & sqlngth < ps[['mxsql']]) {
          ftr_nm <- ftr[['GBFeature_quals']][[
            'GBQualifier']][['GBQualifier_value']]
          tmp <- substr(x=sq, start=strtstp[1], stop=strtstp[2])
          sq_rcrd <- genSqRcrd(accssn=accssn, lctn=lctn,
                               orgnsm=orgnsm, txid='',
                               vrsn=vrsn, sq=tmp, dfln=dfln,
                               gi=gis[[i]], ml_typ=ml_typ,
                               nm=ftr_nm, rcrd_typ='feature')
          res <- c(res, sq_rcrd)
          wftrs <- TRUE
        }
      }
    }
    # if no features, add entire sequence
    if(!wftrs) {
      sqlngth <- as.numeric(rcrd[['GBSeq_length']])
      if(sqlngth < ps[['mxsql']]) {
        sq_rcrd <- genSqRcrd(accssn=accssn, txid='',
                             nm=kywrds, gi=gis[[i]],
                             orgnsm=orgnsm, sq=sq,
                             dfln=dfln, ml_typ=ml_typ,
                             rcrd_typ='whole', vrsn=vrsn)
        res <- c(res, sq_rcrd)
      }
    }
  }
  res
}

#' @name genSqRcrd
#' @title Generate SqRcrd
#' @description Creates an S4 SqRcrd
#' @param accssn Accession ID
#' @param GI GI
#' @param nm Sequence name
#' @param txid Taxonomic ID of source organism
#' @param sq Sequence
#' @param dfln Defintion line
#' @param orgnsm Source organism name
#' @param ml_typ Molecule type
#' @param rcrd_typ Sequence record type
#' @param vrsn Accession version
#' @param lctn Location numbers for features, e.g. '1..200'
genSqRcrd <- function(accssn, gi, nm, txid, sq, dfln,
                      orgnsm, ml_typ, rcrd_typ,
                      vrsn, lctn=NULL) {
  if(is.null(lctn)) {
    id <- vrsn
  } else {
    id <- paste0(vrsn, '/', lctn)
  }
  url <- paste0('https://www.ncbi.nlm.nih.gov/nuccore/', id)
  nncltds <- nchar(sq)
  nambgs <- gregexpr(pattern='[^atcgATCG]', text=sq)[[1]]
  nambgs <- sum(nambgs != -1)
  new('SqRcrd', accssn=accssn, gi=gi, nm=nm, sq=charToRaw(sq),
      dfln=dfln, ml_typ=ml_typ, rcrd_typ=rcrd_typ, vrsn=vrsn, id=id,
      url=url, nncltds=nncltds, nambgs=nambgs, orgnsm=orgnsm,
      txid=txid)
}

#' @name genSqRcrdBx
#' @title Generate SqRcrdBx
#' @description Creates an S4 SqRcrdBx from list of SqRcrds
#' @param sqs List of SqRcrds
genSqRcrdBx <- function(sqs) {
  nambgs <- vapply(sqs, function(x) x@nambgs, 1)
  nncltds <- vapply(sqs, function(x) x@nncltds, 1)
  ids <- vapply(sqs, function(x) x@id, '')
  txids <- vapply(sqs, function(x) x@txid, '')
  new('SqRcrdBx', ids=ids, nncltds=nncltds, nambgs=nambgs,
      txids=txids, sqs=sqs)
}
