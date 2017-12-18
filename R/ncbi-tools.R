#' @name safeSrch
#' @title Safely run rentrez function
#' @description Safely run a rentrez function.
#' If the query fails, the function will retry.
#' @param func rentrez function
#' @param args rentrez function arguments, list
#' @param fnm rentrez function name
#' @param ps Parameters
#' @export
safeSrch <- function(func, args, fnm, ps) {
  # TODO: wait times?
  res <- NULL
  for(i in 1:ps[['mxretry']]) {
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
      info(lvl=1, ps=ps, "Retry [", i, "] for [", fnm, ']')
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
#' @export
# @Hannes: eqv to .num.seqs.for.txid
nSqs <- function(txid, ps, direct=FALSE) {
  org_term <- ifelse(direct, '[Organism:noexp]',
                     '[Organism:exp]' )
  args <- list(db='nucleotide',
               term=paste0('txid', txid,
                           org_term, '1:',
                           ps[['mxsql']], '[SLEN]'),
               use_history=TRUE, retmax=1)
  res <- safeSrch(func=rentrez::entrez_search,
                  args=args, fnm='entrez_search',
                  ps=ps)
  res[['count']]
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
  org_trm <- ifelse(direct, '[Organism:noexp]',
                     '[Organism:exp]')
  allseqs <- numeric()
  args <- list(db='nucleotide',
               term=paste0('txid', txid, org_trm, '1:',
                           ps[['mxsql']], '[SLEN]'),
               use_history=TRUE,
               ## XXX Without 1e9-1 one could save time...
               # but are the hits random???
               # @Hannes: from exp. they're not random,
               #  easily tested though.
               retmax=1e9-1)
  srch <- safeSrch(func=rentrez::entrez_search,
                   args=args, fnm='search',
                   ps=ps)
  gis <- srch[['ids']];
  if(length(gis) < 1) {
    return(list())
  }
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
    info(lvl=1, ps=ps, "Retreiving seqs [",
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
                    ps=ps)
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
                          ps=ps)
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
    info(lvl=1, ps=ps, "Done retreiving [", length(seqs),
        "(", i+1, "-",
        ifelse(i+btch<length(gis),
               i+btch, length(gis)),
        ")] seqs for taxid [", txid, "]")
    
    allseqs <- c(allseqs, seqs)
  }
  info(lvl=1, ps=ps, "Finished retreiving [",
      length(allseqs), "] sequences for taxid [",
      txid, "]")
  names(allseqs) <- gis
  allseqs
}
