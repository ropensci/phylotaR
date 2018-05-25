#' @name rwRcrd2SqRcrd
#' @title Convert raw Entrez gbwithparts record to SqRcrds
#' @description Parses returned sequences features with Entrez,
#' returns a SqRcrd for each raw record.
#' @param rw_rcrds Raw records returned from Entrez fetch
#' @param gis GIs of each fectched record
#' @param ps Parameters list
rwRcrd2SqRcrd <- function(rw_rcrds, gis, ps) {
  # gis <- '728952670'
  # rw_rcrds <- rentrez::entrez_fetch(db="nucleotide",
  #                                   rettype='gbwithparts',
  #                                   retmode='xml',
  #                                   id=gis)
  res <- NULL
  rcrds <- try(XML::xmlToList(rw_rcrds), silent=TRUE)
  if(inherits(rcrds, 'try-error')) {
    msg <- paste0('XML parsing error with one or more of the following GIs:\n[',
                  paste0(gis, collapse=','), ']\nSkipping...\n',
                  '(This can occur with GenBank data. ',
                  'Trying restarting pipeline to see if warning reoccurs. ',
                  'If it does, contact maintainer.)')
    warn(ps=ps, msg)
    return(NULL)
  }
  for(i in seq_along(rcrds)) {
    # for debugging.... save last GI
    svObj(wd=ps[['wd']], obj=gis[[i]], nm='last_gi')
    rcrd <- rcrds[[i]]
    wftrs <-  FALSE
    # key info
    accssn <- rcrd[['GBSeq_primary-accession']]
    vrsn <- rcrd[["GBSeq_accession-version"]]
    if (is.null(vrsn)) {
      next
    }
    ml_typ <- rcrd[['GBSeq_moltype']]
    if(is.null(ml_typ)) {
      # ml_typ not always recorded, e.g. NR_040059
      ml_typ <- ''
    }
    sq <- rcrd[['GBSeq_sequence']]
    create_date <- rcrd[["GBSeq_create-date"]]
    create_date <- as.Date(create_date,
                           format="%d-%b-%Y")
    age <- as.integer(difftime(ps[['date']],
                               create_date, units = 'days'))
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
      # TODO: perhaps make sure locations don't overlap
      if(grepl('^source$', ftr[['GBFeature_key']],
               ignore.case=TRUE)) {
        next
      }
      lctn <- ftr[['GBFeature_location']]
      lctn <- gsub('[^0-9\\.]', '', lctn)
      if(nchar(lctn) == 0) {
        next
      }
      strtstp <- strsplit(x=lctn, split='\\.\\.')[[1]]
      if(length(strtstp) != 2) {
        next
      }
      strtstp <- as.numeric(strtstp)
      sqlngth <- strtstp[2] - strtstp[1]
      if(all(strtstp == c(1, nchar(sq)))) {
        # feature should not span entire sequence
        next
      }
      if(sqlngth > ps[['mnsql']] & sqlngth < ps[['mxsql']]) {
        ftr_nm <- ftr[['GBFeature_quals']][[
          'GBQualifier']][['GBQualifier_value']]
        if(is.null(ftr_nm)) {
          ftr_nm <- ftr[['GBFeature_key']]
        }
        tmp <- substr(x=sq, start=strtstp[1], stop=strtstp[2])
        sq_rcrd <- genSqRcrd(accssn=accssn, lctn=lctn, age=age,
                             orgnsm=orgnsm, txid='',
                             vrsn=vrsn, sq=tmp, dfln=dfln,
                             gi=gis[[i]], ml_typ=ml_typ,
                             nm=ftr_nm, rcrd_typ='feature')
        res <- c(res, sq_rcrd)
        wftrs <- TRUE
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
                             rcrd_typ='whole', vrsn=vrsn,
                             age=age)
        res <- c(res, sq_rcrd)
      }
    }
  }
  ids <- vapply(res, function(x) x@id, '')
  nondups <- !duplicated(ids)
  res[nondups]
}

#' @name genSqRcrd
#' @title Generate SqRcrd
#' @description Creates an S4 SqRcrd
#' @param accssn Accession ID
#' @param gi GI
#' @param nm Sequence name
#' @param txid Taxonomic ID of source organism
#' @param sq Sequence
#' @param dfln Defintion line
#' @param orgnsm Source organism name
#' @param ml_typ Molecule type
#' @param rcrd_typ Sequence record type
#' @param vrsn Accession version
#' @param age Number of days since upload
#' @param lctn Location numbers for features, e.g. '1..200'
genSqRcrd <- function(accssn, gi, nm, txid, sq, dfln,
                      orgnsm, ml_typ, rcrd_typ,
                      vrsn, age, lctn=NULL) {
  if(is.null(lctn)) {
    id <- vrsn
  } else {
    id <- paste0(vrsn, '/', lctn)
  }
  url <- paste0('https://www.ncbi.nlm.nih.gov/nuccore/', id)
  nncltds <- nchar(sq)
  unambsq <- gsub(pattern='[^atcgATCG]', replacement='',
                  x=sq)[[1]]
  nambgs <- nncltds - nchar(unambsq)
  pambgs <- nambgs/nncltds
  gcr <- gregexpr(pattern='[^cgCG]', text=unambsq)[[1]]
  gcr <- length(gcr)/nchar(unambsq)
  new('SqRcrd', accssn=accssn, gi=gi, nm=nm, sq=charToRaw(sq),
      dfln=dfln, ml_typ=ml_typ, rcrd_typ=rcrd_typ, vrsn=vrsn, id=id,
      url=url, nncltds=nncltds, nambgs=nambgs, orgnsm=orgnsm,
      txid=txid, pambgs=pambgs, gcr=gcr, age=age)
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
