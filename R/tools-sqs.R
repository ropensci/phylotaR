#' @name gb_extract
#' @title Extract elements from a raw GenBank record
#' @description Returns a list of elements from a GenBank record such as
#' 'organism', 'sequence' and features.
#' @details Uses restez extract functions. See restez package for more details.
#' @param record raw GenBank text record
#' @family run-private
#' @return list of GenBank elements
gb_extract <- function(record) {
  accession <- restez::gb_extract(record, what = 'accession')
  vrsn <- restez::gb_extract(record, what = 'version')
  locus <- restez::gb_extract(record, what = 'locus')
  organism <- restez::gb_extract(record, what = 'organism')
  seq <- restez::gb_extract(record, what = 'sequence')
  def <- restez::gb_extract(record, what = 'definition')
  keywords <- restez::gb_extract(record, what = 'keywords')
  features <- restez::gb_extract(record, what = 'features')
  date <- locus[['date']]
  moltype <- paste0(locus[['mol']], '-', locus[['type']])
  list('accession' = accession, 'version' = vrsn, 'definition' = def,
       'sequence' = seq, 'moltype' = moltype, 'features' = features,
       'date' = date, 'organism' = organism, 'keywords' = keywords)
}

#' @name rawseqrec_breakdown
#' @title Breakdown a sequence record into its features
#' @description Takes GenBank record's elements and returns a SeqRec. For
#' sequences with lots of features, the sequence is broken down into these
#' features provided they are of the right size. Sequences are either returned
#' as features or whole sequence records, never both.
#' @param record_parts list of record elements from a GenBank record
#' @template ps
#' @family run-private
#' @return list of SeqRecs
rawseqrec_breakdown <- function(record_parts, ps) {
  # get details
  seq <- record_parts[['sequence']]
  accssn <- record_parts[['accession']]
  orgnsm <- record_parts[['organism']]
  create_date <- as.Date(record_parts[['date']], format = "%d-%b-%Y")
  age <- as.integer(difftime(ps[['date']], create_date, units = 'days'))
  vrsn <- record_parts[['version']]
  dfln <- record_parts[['definition']]
  moltype <- record_parts[['moltype']]
  kywrds <- record_parts[['keywords']]
  # filter features
  features <- record_parts[['features']]
  types <- vapply(X = features, FUN = '[[', FUN.VALUE = character(1), 'type')
  pull <- types != 'source'
  while (any(pull)) {
    locations <- vapply(X = features, FUN = '[[', FUN.VALUE = character(1),
                        'location')
    pull <- grepl(pattern = '^[0-9]+\\.\\.[0-9]+$', x = locations) & pull
    features <- features[pull]
    locations <- locations[pull]
    # remove dups
    pull <- !duplicated(locations)
    features <- features[pull]
    locations <- locations[pull]
    # filter by length
    locations <- strsplit(x = locations, split = '\\.\\.')
    locations <- lapply(X = locations, FUN = as.integer)
    lengths <- vapply(X = locations, FUN = function(x) x[[2]] - x[[1]],
                      FUN.VALUE = integer(1))
    pull <- lengths > ps[['mnsql']] & lengths < ps[['mxsql']]
    features <- features[pull]
    locations <- locations[pull]
    pull <- FALSE
  }
  # extract sequences
  if (any(pull)) {
    seqrecs <- vector(mode = 'list', length = length(features))
    for (i in seq_along(features)) {
      strt <- locations[[i]][[1]]
      stp <- locations[[i]][[2]]
      ftr_nm <- features[[i]][[3]]
      lctn <- paste0(strt, '..', stp)
      tmp <- substr(x = seq, start = strt, stop = stp)
      seqrecs[[i]] <- seqrec_gen(accssn = accssn, lctn = lctn, age = age,
                                 orgnsm = orgnsm, txid = '', vrsn = vrsn,
                                 sq = tmp, dfln = dfln, ml_typ = moltype,
                                 nm = ftr_nm, rec_typ = 'feature')
    }
  } else {
    seqrecs <- seqrec_gen(accssn = accssn, txid = '', nm = kywrds, age = age,
                          orgnsm = orgnsm, sq = seq, dfln = dfln, vrsn = vrsn,
                          ml_typ = moltype, rec_typ = 'whole')
    seqrecs <- list(seqrecs)
  }
  seqrecs
}

#' @name seqrec_convert
#' @title Convert raw Entrez gb text record to SeqRecs
#' @description Parses returned sequences features with Entrez, returns one or
#' more SeqRec objects for each raw record.
#' @param raw_recs Raw text records returned from Entrez fetch
#' @template ps
#' @family run-private
#' @return SeqRecs
seqrec_convert <- function(raw_recs, ps) {
  # gis <- c('558614019', '558051', '1283191328', 'NR_040059')
  # raw_recs <- rentrez::entrez_fetch(db = "nucleotide",
  #                                   rettype = 'gbwithparts',
  #                                   retmode = 'text',
  #                                   id = gis)
  #  write(x = raw_recs, file = 'raw_recs.txt')
  recs <- strsplit(x = raw_recs, split = '//\n\n')[[1]]
  recs <- paste0(recs, '//') # record not complete without //
  records_parts <- lapply(X = recs, FUN = gb_extract)
  pull <- vapply(X = records_parts, FUN = function(x) x[['sequence']] != '',
                 logical(1))
  records_parts <- records_parts[pull]
  seqrecs <- unlist(lapply(X = records_parts, FUN = rawseqrec_breakdown,
                           ps = ps))
  seqrecs
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
