#' @name seqrec_convert
#' @title Convert raw Entrez gbwithparts record to SeqRecs
#' @description Parses returned sequences features with Entrez,
#' returns a SeqRec for each raw record.
#' @param raw_recs Raw records returned from Entrez fetch
#' @template ps
#' @family run-private
#' @return SeqRecs
seqrec_convert <- function(raw_recs, ps) {
  # gis <- '558614019'
  # raw_recs <- rentrez::entrez_fetch(db = "nucleotide",
  #                                   rettype = 'gbwithparts',
  #                                   retmode = 'xml',
  #                                   id = gis)
  res <- NULL
  recs <- try(XML::xmlToList(raw_recs), silent = TRUE)
  if (inherits(recs, 'try-error')) {
    msg <- paste0('XML parsing error:',
                  'This can occur with GenBank data. ',
                  'Try pipeline again to see if warning reoccurs. ',
                  'If it does, contact maintainer.')
    warn(ps = ps, msg)
    return(NULL)
  }
  for (i in seq_along(recs)) {
    rec <- recs[[i]]
    # for debugging.... save last record
    obj_save(wd = ps[['wd']], obj = rec, nm = 'last_seqrec')
    wftrs <-  FALSE
    # key info
    accssn <- rec[['GBSeq_primary-accession']]
    vrsn <- rec[["GBSeq_accession-version"]]
    ml_typ <- rec[['GBSeq_moltype']]
    if (is.null(vrsn)) {
      next
    }
    if(is.null(ml_typ)) {
      # ml_typ not always recorded, e.g. NR_040059
      ml_typ <- ''
    }
    sq <- rec[['GBSeq_sequence']]
    create_date <- rec[["GBSeq_create-date"]]
    create_date <- as.Date(create_date,
                           format="%d-%b-%Y")
    age <- as.integer(difftime(ps[['date']],
                               create_date, units = 'days'))
    if(is.null(sq)) {
      # master records have no sequences
      # e.g. https://www.ncbi.nlm.nih.gov/nuccore/1283191328
      next
    }
    dfln <- rec[['GBSeq_definition']]
    orgnsm <- rec[['GBSeq_organism']]
    ftr_tbl <- rec[['GBSeq_feature-table']]
    kywrds <- vapply(rec[['GBSeq_keywords']], '[', '')
    kywrds <- paste0(kywrds, collapse = ' | ')
    # extract features
    for (ftr in ftr_tbl) {
      # TODO: perhaps make sure locations don't overlap
      if (grepl('^source$', ftr[['GBFeature_key']],
               ignore.case = TRUE)) {
        next
      }
      lctn <- ftr[['GBFeature_location']]
      lctn <- gsub('[^0-9\\.]', '', lctn)
      if (nchar(lctn) == 0) {
        next
      }
      strtstp <- strsplit(x = lctn, split = '\\.\\.')[[1]]
      if (length(strtstp) != 2) {
        next
      }
      strtstp <- as.numeric(strtstp)
      sqlngth <- strtstp[2] - strtstp[1]
      if (all(strtstp == c(1, nchar(sq)))) {
        # feature should not span entire sequence
        next
      }
      if (sqlngth > ps[['mnsql']] & sqlngth < ps[['mxsql']]) {
        ftr_nm <- ftr[['GBFeature_quals']][[
          'GBQualifier']][['GBQualifier_value']]
        if (is.null(ftr_nm)) {
          ftr_nm <- ftr[['GBFeature_key']]
        }
        tmp <- substr(x = sq, start = strtstp[1], stop = strtstp[2])
        seqrec <- seqrec_gen(accssn = accssn, lctn = lctn, age = age,
                             orgnsm = orgnsm, txid = '', vrsn = vrsn, sq = tmp,
                             dfln = dfln, ml_typ = ml_typ, nm = ftr_nm,
                             rec_typ = 'feature')
        res <- c(res, seqrec)
        wftrs <- TRUE
      }
    }
    # if no features, add entire sequence
    if (!wftrs) {
      sqlngth <- as.numeric(rec[['GBSeq_length']])
      if (sqlngth < ps[['mxsql']]) {
        seqrec <- seqrec_gen(accssn = accssn, txid = '', nm = kywrds,
                             orgnsm = orgnsm, sq = sq, dfln = dfln,
                             ml_typ = ml_typ, rec_typ = 'whole',
                             vrsn = vrsn, age = age)
        res <- c(res, seqrec)
      }
    }
  }
  # filter out ambiguous sequences and dups
  nonambgs <- vapply(res, function(x) x@pambgs < 0.1, logical(1))
  res <- res[nonambgs]
  ids <- vapply(res, function(x) x@id, '')
  nondups <- !duplicated(ids)
  res[nondups]
}

#' @name seqrec_gen
#' @title Generate sequence record
#' @description Creates an S4 SeqRec
#' @param accssn Accession ID
#' @param nm Sequence name
#' @param txid Taxonomic ID of source organism
#' @param sq Sequence
#' @param dfln Definition line
#' @param orgnsm Source organism name
#' @param ml_typ Molecule type
#' @param rec_typ Sequence record type
#' @param vrsn Accession version
#' @param age Number of days since upload
#' @param lctn Location numbers for features, e.g. '1..200'
#' @return SeqRec
#' @family run-private
seqrec_gen <- function(accssn, nm, txid, sq, dfln, orgnsm, ml_typ,
                       rec_typ, vrsn, age, lctn=NULL) {
  if (is.null(lctn)) {
    id <- vrsn
  } else {
    id <- paste0(vrsn, '/', lctn)
  }
  url <- paste0('https://www.ncbi.nlm.nih.gov/nuccore/', id)
  nncltds <- nchar(sq)
  unambsq <- gsub(pattern = '[^atcgATCG]', replacement = '',
                  x = sq)[[1]]
  nambgs <- nncltds - nchar(unambsq)
  pambgs <- nambgs/nncltds
  gcr <- gregexpr(pattern = '[^cgCG]', text = unambsq)[[1]]
  gcr <- length(gcr)/nchar(unambsq)
  new('SeqRec', accssn = accssn, nm = nm, sq = charToRaw(sq),
      dfln = dfln, ml_typ = ml_typ, rec_typ = rec_typ, vrsn = vrsn, id = id,
      url = url, nncltds = nncltds, nambgs = nambgs, orgnsm = orgnsm,
      txid = txid, pambgs = pambgs, gcr = gcr, age = age)
}

#' @name seqarc_gen
#' @title Generate sequence archive
#' @description Creates an S4 SeqArc from list of SeqRecs
#' @param seqrecs List of SeqRecs
#' @return SeqArc
#' @family run-private
seqarc_gen <- function(seqrecs) {
  nambgs <- vapply(seqrecs, function(x) x@nambgs, 1)
  nncltds <- vapply(seqrecs, function(x) x@nncltds, 1)
  ids <- vapply(seqrecs, function(x) x@id, '')
  txids <- vapply(seqrecs, function(x) x@txid, '')
  new('SeqArc', ids = ids, nncltds = nncltds, nambgs = nambgs, txids = txids,
      sqs = seqrecs)
}
