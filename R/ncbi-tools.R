#' @name safeSrch
#' @title Safely run rentrez function
#' @description Safely run a rentrez function.
#' If the query fails, the function will retry.
#' @param func rentrez function
#' @param args rentrez function arguments, list
#' @param fnm rentrez function name
#' @param verbose Verbose? T/F
#' @param mx_retry Maximum number of attempts
#' @export
safeSrch <- function(func, args,
                     fnm, verbose=FALSE,
                     mx_retry=1000) {
  # TODO: wait times?
  res <- NULL
  for(i in 1:mx_retry) {
    query <- try(do.call(func, args),
                 silent=TRUE)
    if(!is(query, "try-error") ) {
      res <- query
      break
    } else {
      .cp(v=verbose,
          "Retry [", i, "] for [",
          fnm, ']')
    }
  }
  res
}

#' @name nSqs
#' @title Count number of sequences for txid
#' @description Return the number of sequences
#' associated with a taxonomic ID on NCBI GenBank.
#' @param txid Taxonomic ID
#' @param direct exp or noexp? T/F
#' @param mx_len Maximum number of sequences, numeric.
#' @param verbose Verbose? T/F
#' @export
# @Hannes: eqv to .num.seqs.for.txid
nSqs <- function(txid, direct=FALSE,
                 mx_len=25000, verbose=FALSE) {
  org_term <- ifelse(direct, '[Organism:noexp]',
                     '[Organism:exp]' )
  args <- list(db='nucleotide',
               term=paste0('txid', txid,
                           org_term, '1:',
                           mx_len, '[SLEN]'),
               use_history=TRUE, retmax=1)
  res <- safeSrch(func=rentrez::entrez_search,
                  args=args, fnm='entrez_search',
                  verbose=verbose)
  res[['count']]
}

#' @title dwnldFrmNCBI
#' @description Given a taxon ID, queries NCBI
#' for sequences that match the taxid.
#' @param txid NCBI taxon identifier
#' @param direct Whether to download all 'directly' linked
#' sequences ([Organism:noexp] in Genbank search) or all 'subtree'
#' sequences ([Organism:exp])
#' @param mx_lngth Maximum length of sequence to be included in GI list
#' @param mx_sqs Take a random subset of mx_sqs if the number of
#' seqs exceeds it.
#' Only makes sense when direct=TRUE
#' @return list of lists containing sequence objects
#' @export
dwnldFrmNCBI <- function(txid, direct=FALSE,
                         mx_lngth=25000, mx_sqs=100000,
                         verbose=FALSE) {
  # test w/ golden moles 9389
  org_trm <- ifelse(direct, '[Organism:noexp]',
                     '[Organism:exp]')
  allseqs <- numeric()
  args <- list(db='nucleotide',
               term=paste0('txid', txid, org_trm, '1:',
                           mx_lngth, '[SLEN]'),
               use_history=TRUE,
               ## XXX Without 1e9-1 one could save time...
               # but are the hits random???
               # @Hannes: from exp. they're not random,
               #  easily tested though.
               retmax=1e9-1)
  srch <- safeSrch(func=rentrez::entrez_search,
                   args=args, fnm='search',
                   verbose=FALSE)
  gis <- srch[['ids']];
  if(length(gis) < 1) {
    return(list())
  }
  # If the maximum amount of number is exceeded,
  # randomly pick mx_sqs gis
  if(direct && length(gis) > mx_sqs) {
    .cp(v=verbose, "Choosing [", mx_sqs,
        "] random sequences from [",
        length(gis), "] available")
    gis <- gis[sample(1:mx_sqs, replace=FALSE)]
  }
  .cp(v=verbose, "Retrieving [", length(gis),
      "] sequences for taxid [", txid, "]")
  # Fetch sequences in batches
  btch <- 500
  for(i in seq(0, length(gis)-1, btch)) {
    ## Get FASTA strings for IDs in the specified segment
    .cp(v=verbose, "Retreiving seqs [",
        i+1, "-", ifelse(i+btch<length(gis),
               i+btch, length(gis)),
        "] for taxid [", txid, "]");
    seqargs <- list(db="nuccore",
                    rettype="fasta",
                    web_history=srch[['web_history']],
                    retmax=btch,
                    retstart=i)
    res <- safeSrch(func=rentrez::entrez_fetch,
                    fnm='fetch', args=seqargs,
                    verbose=verbose)
    # split fasta string on '>'
    seqstrs <- unlist(strsplit(res, "(?<=[^>])(?=\n>)",
                               perl=TRUE))
    # clear defline
    seqstrs <- gsub("^\n?>.*?\\n", "", seqstrs)
    # clear newlines
    seqstrs <- gsub("\n", "", seqstrs, fixed=TRUE)
    # Get summary object to retreive taxid, accession,
    # and other parameters for seqs
    summaries <- safeSrch(func=rentrez::entrez_summary,
                          fnm='summary',
                          args=list(db='nucleotide',
                                    web_history=srch[['web_history']],
                                    retmax=btch,
                                    retstart=i),
                          verbose=verbose)
    # cnvrt to list when only one result
    if(length(seqstrs)==1) {
      summaries <- list(summaries)
    }
    seqs <- lapply(1:length(seqstrs), function(i) {
      sm <- summaries[[i]]
      se <- seqstrs[[i]]
      seq <- list(gi=sm$uid,
                  ti=sm$taxid,
                  acc=sm$caption,
                  acc_vers=sm$accessionversion,
                  length=sm$slen,
                  ##TODO: Division not in summary object, how to get?
                  division=NA,
                  acc_date=sm$createdate,
                  ##TODO: GBrel not in summary object, how to get?
                  gbrel=NA,
                  def=sm$title,
                  seq=se)
      seq
    })
    .cp(v=verbose, "Done retreiving [", length(seqs),
        "(", i+1, "-",
        ifelse(i+btch<length(gis),
               i+btch, length(gis)),
        ")] seqs for taxid [", txid, "]")
    
    allseqs <- c(allseqs, seqs)
  }
  .cp(v=verbose, "Finished retreiving [",
      length(allseqs), "] sequences for taxid [",
      txid, "]")
  names(allseqs) <- gis
  allseqs
}
