#' @name safeSrch
#' @title Safely run rentrez function
#' @description Safely run a rentrez function.
#' If the query fails, the function will retry.
#' All query results are cached. To remove cached data
#' use hard reset.
#' @param func rentrez function
#' @param args rentrez function arguments, list
#' @param fnm rentrez function name
#' @param ps Parameters
#' @export
safeSrch <- function(func, args, fnm, ps) {
  res <- ldNcbiCch(fnm=fnm, args=args, wd=ps[['wd']])
  if(!is.null(res)) {
    return(res)
  }
  for(wt_tm in ps[['wt_tms']]) {
    query <- try(do.call(func, args),
                 silent=TRUE)
    if(!is(query, "try-error")) {
      res <- query
      break
    } else {
      # ctrl+c
      if(grepl('Operation was aborted by an application callback',
            query[[1]])) {
        stop(query[[1]])
      }
      info(lvl=1, ps=ps, "Retrying in [", wt_tm, "s] for [",
           fnm, ']')
      Sys.sleep(wt_tm)
    }
  }
  svNcbiCch(fnm=fnm, args=args, wd=ps[['wd']], obj=res)
  res
}

#' @name nSqs
#' @title Count number of sequences for txid
#' @description Return the number of sequences
#' associated with a taxonomic ID on NCBI GenBank.
#' @param txid Taxonomic ID
#' @param direct exp or noexp? T/F
#' @export
# @Hannes: eqv to .num.seqs.for.txid
nSqs <- function(txid, ps, direct=FALSE) {
  org_term <- ifelse(direct, '[Organism:noexp]',
                     '[Organism:exp]' )
  term <- paste0('txid', txid, org_term, '1:',
                 ps[['mxsql']], '[SLEN]')
  args <- list(db='nucleotide', retmax=0, term=term)
  res <- safeSrch(func=rentrez::entrez_search,
                  args=args, fnm='search',
                  ps=ps)
  res[['count']]
}


#' @title Return random set of sequence IDs
#' @description For a given txid return a random set of 
#' sequences associated.
#' @details For model organisms downloading all IDs can take long
#' or even cause an xml parsing error. For any search with more than
#' hrdmx sequences, this function we will run multiple small searches
#' downloading retmax seq IDs at a time with different retstart values
#' to generate a semi-random vector of sequence IDs. For all other
#' searches, all IDs will be retrieved.
#' Note, it makes no sense for mxsqs in parameters to be greater
#' than hrdmx in this function.
#' @param txid NCBI taxon identifier
#' @param direct Whether to download all 'directly' linked
#' sequences ([Organism:noexp] in Genbank search) or all 'subtree'
#' sequences ([Organism:exp])
#' @param ps Parameters
#' @param sqcnt Sequence count as determined with nSqs()
#' @param retmax Maximum number of sequences when querying model
#' organisms. The smaller the more random, the larger the faster.
#' @param hrdmx Absolute maximum number of sequence IDs to download
#' in a single query.
#' @return vector ot IDs
#' @export
getGIs <- function(txid, direct, sqcnt, ps, retmax=1000,
                   hrdmx=100000) {
  org_term <- ifelse(direct, '[Organism:noexp]',
                     '[Organism:exp]' )
  term <- paste0('txid', txid, org_term, '1:',
                 ps[['mxsql']], '[SLEN]')
  if(sqcnt > hrdmx) {
    ids <- NULL
    retstarts <- sample(1:(sqcnt-retmax),
                       round(hrdmx/retmax))
    for(retstart in retstarts) {
      args <- list(db='nucleotide', retmax=retmax,
                   term=term, retstart=retstart)
      srch <- safeSrch(func=rentrez::entrez_search,
                       args=args, fnm='search', ps=ps)
      ids <- c(ids, srch[['ids']])
    }
  } else {
    args <- list(db='nucleotide', retmax=sqcnt, term=term)
    srch <- safeSrch(func=rentrez::entrez_search,
                     args=args, fnm='search',
                     ps=ps)
    ids <- srch[['ids']]
  }
  ids
}

#' @title dwnldFrmNCBI
#' @description Given a taxon ID, queries NCBI
#' for sequences that match the taxid.
#' @param txid NCBI taxon identifier
#' @param direct Whether to download all 'directly' linked
#' sequences ([Organism:noexp] in Genbank search) or all 'subtree'
#' sequences ([Organism:exp])
#' seqs exceeds it.
#' Only makes sense when direct=TRUE
#' @return list of lists containing sequence objects
#' @export
dwnldFrmNCBI <- function(txid, ps, direct=FALSE) {
  # test w/ golden moles 9389
  allseqs <- numeric()
  sqcnt <- nSqs(txid=txid, ps=ps, direct=direct)
  if(sqcnt < 1) {
    return(list())
  }
  gis <- getGIs(txid=txid, direct=direct, sqcnt=sqcnt,
                ps=ps)
  # If the maximum amount of number is exceeded,
  # randomly pick ps[['mxsqs']] gis
  if(direct && length(gis) > ps[['mxsqs']]) {
    info(lvl=1, ps=ps, "Choosing [", ps[['mxsqs']],
        "] random sequences from [",
        length(gis), "] available")
    gis <- gis[sample(1:ps[['mxsqs']], replace=FALSE)]
  }
  info(lvl=1, ps=ps, "Retrieving [", length(gis),
      "] sequences for taxid [", txid, "]")
  # Fetch sequences in batches
  btch <- 500
  for(i in seq(0, length(gis)-1, btch)) {
    # Get FASTA strings for IDs in the specified segment
    lower <- i+1
    upper <- ifelse(i+btch<length(gis), i+btch, length(gis))
    crrnt_ids <- gis[lower:upper]
    info(lvl=2, ps=ps, "Retreiving seqs [",
        lower, "-", upper, "] for taxid [", txid, "]");
    args <- list(db="nuccore", rettype="fasta",
                 id=crrnt_ids)
    res <- safeSrch(func=rentrez::entrez_fetch,
                    fnm='fetch', args=args, ps=ps)
    seqstrs <- unlist(strsplit(res, "(?<=[^>])(?=\n>)",
                               perl=TRUE))
    seqstrs <- gsub("^\n?>.*?\\n", "", seqstrs)
    seqstrs <- gsub("\n", "", seqstrs, fixed=TRUE)
    args=list(db='nucleotide', id=crrnt_ids)
    summaries <- safeSrch(func=rentrez::entrez_summary,
                          fnm='summary', args=args, ps=ps)
    if(length(seqstrs) == 1) {
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
    info(lvl=1, ps=ps, "Retreived [", length(seqs),
        "(", lower, "-", upper,
        ")] seqs for taxid [", txid, "]")
    
    allseqs <- c(allseqs, seqs)
  }
  info(lvl=1, ps=ps, "Finished retreiving [",
      length(allseqs), "] sequences for taxid [",
      txid, "]")
  names(allseqs) <- gis
  allseqs
}
