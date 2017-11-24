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
