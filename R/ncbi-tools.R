#' @name nNcbiNds
#' @title Count number of descending taxonomic nodes
#' @description Searches NCBI taxonomy and returns
#' number of descendents taxonomic nodes (species, genera ...)
#' of ID.
#' @param txid Taxonomic ID
#' @param ps Parameters
#' @export
nNcbiNds <- function(txid, ps) {
  trm <- paste0('txid', txid, '[Subtree]')
  args <- list(db='taxonomy', retmax=0, term=trm)
  res <- srchNCch(func=rentrez::entrez_search,
                  args=args, fnm='search',
                  ps=ps)
  res[['count']]
}

#' @name nSqs
#' @title Count number of sequences for txid
#' @description Return the number of sequences
#' associated with a taxonomic ID on NCBI GenBank.
#' @param txid Taxonomic ID
#' @param direct exp or noexp? T/F
#' @export
nSqs <- function(txid, ps, direct=FALSE) {
  org_term <- ifelse(direct, '[Organism:noexp]',
                     '[Organism:exp]' )
  term <- paste0('txid', txid, org_term, ps[['mnsql']], ':',
                 ps[['mxsql']], '[SLEN]')
  args <- list(db='nucleotide', retmax=0, term=term)
  res <- srchNCch(func=rentrez::entrez_search,
                  args=args, fnm='search',
                  ps=ps)
  res[['count']]
}


#' @title Return random set of sequence IDs
#' @description For a given txid return a random set of 
#' sequences associated.
#' @details For model organisms downloading all IDs can a take long time
#' or even cause an xml parsing error. For any search with more than
#' hrdmx sequences, this function we will run multiple small searches
#' downloading retmax seq IDs at a time with different retstart values
#' to generate a semi-random vector of sequence IDs. For all other
#' searches, all IDs will be retrieved.
#' Note, it makes no sense for mdlthrs in parameters to be greater
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
getGIs <- function(txid, direct, sqcnt, ps, retmax=100,
                   hrdmx=100000) {
  org_term <- ifelse(direct, '[Organism:noexp]',
                     '[Organism:exp]' )
  term <- paste0('txid', txid, org_term,
                 ps[['mnsql']], ':', ps[['mxsql']],
                 '[SLEN]')
  if(sqcnt <= hrdmx) {
    args <- list(db='nucleotide', retmax=sqcnt, term=term)
    srch <- srchNCch(func=rentrez::entrez_search,
                     args=args, fnm='search',
                     ps=ps)
    return(srch[['ids']])
  }
  srch_args <- list(db='nucleotide', term=term,
                    use_history=TRUE, retmax=0)
  ids <- NULL
  ret_strts <- sample(seq(1, (sqcnt-retmax), retmax),
                      round(hrdmx/retmax), replace=FALSE)
  rsrch <- TRUE
  for(ret_strt in ret_strts) {
    mxtry <- 3
    for(i in 1:mxtry) {
      # use webobject, search without cache
      if(rsrch) {
        srch <- safeSrch(func=rentrez::entrez_search,
                         args=srch_args, fnm='search',
                         ps=ps)
      }
      ftch_args <- list(db='nucleotide',
                        web_history=srch[['web_history']],
                        retmax=retmax,
                        rettype='GI',
                        retstart=ret_strt)
      id_ftch <- try(do.call(what=rentrez::entrez_fetch,
                             args=ftch_args), silent=TRUE)
      if(inherits(id_ftch, 'try-error')) {
        rsrch <- TRUE
        if(i == mxtry) {
          error(ps=ps, 'Failed to get IDs for model organism [',
                txid, ']')
        }
      }
      rsrch <- FALSE
      break
    }
    ids <- c(ids, strsplit(id_ftch, '\n')[[1]])
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
dwnldFrmNCBI <- function(txid, ps, direct=FALSE, lvl=0) {
  # test w/ golden moles 9389
  allseqs <- numeric()
  sqcnt <- nSqs(txid=txid, ps=ps, direct=direct)
  if(sqcnt < 1) {
    return(list())
  }
  gis <- getGIs(txid=txid, direct=direct, sqcnt=sqcnt,
                ps=ps)
  if(direct && length(gis) > ps[['mdlthrs']]) {
    info(lvl=lvl+3, ps=ps, "More than [", ps[['mdlthrs']], ' sqs] available.',
         ' Choosing at random.')
    gis <- gis[sample(1:length(gis), size=ps[['mdlthrs']],
                      replace=FALSE)]
  }
  info(lvl=lvl+3, ps=ps, "Getting [", length(gis), " sqs] ...")
  # Fetch sequences in batches
  btch <- 500
  for(i in seq(0, length(gis)-1, btch)) {
    # Get FASTA strings for IDs in the specified segment
    lower <- i+1
    upper <- ifelse(i+btch<length(gis), i+btch, length(gis))
    crrnt_ids <- gis[lower:upper]
    info(lvl=lvl+4, ps=ps, "[", lower, "-", upper, "]");
    args <- list(db="nucleotide", rettype="fasta",
                 id=crrnt_ids)
    res <- srchNCch(func=rentrez::entrez_fetch,
                    fnm='fetch', args=args, ps=ps)
    seqstrs <- unlist(strsplit(res, "(?<=[^>])(?=\n>)",
                               perl=TRUE))
    seqstrs <- gsub("^\n?>.*?\\n", "", seqstrs)
    seqstrs <- gsub("\n", "", seqstrs, fixed=TRUE)
    args=list(db='nucleotide', id=crrnt_ids)
    summaries <- srchNCch(func=rentrez::entrez_summary,
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
    allseqs <- c(allseqs, seqs)
  }
  names(allseqs) <- gis
  allseqs
}
