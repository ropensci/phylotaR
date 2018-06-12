#' @name searchterm_gen
#' @title Construct GenBank Search Term
#' @description Construct search term for searching GenBank's
#' nucleotide database. Limits the maximum size of sequences, avoids
#' whole genome shotguns, predicted, unverified and synthetic
#' sequences.
#' @param txid Taxonomic ID
#' @template ps
#' @param direct Node-level only or subtree as well? Default FALSE.
#' @family run-private
#' @return character, search term
searchterm_gen <- function(txid, ps, direct = FALSE) {
  org_trm <- ifelse(direct, '[Organism:noexp]', '[Organism:exp]' )
  avd1 <- ' NOT predicted[TI] NOT "whole genome shotgun"[TI]'
  avd2 <- ' NOT unverified[TI] NOT "synthetic construct"[Organism]'
  paste0('(txid', txid, org_trm, ' AND ', ps[['mnsql']],
         ':', ps[['mxsql']], '[SLEN])', avd1, avd2)
}


#' @name txnds_count
#' @title Count number of descending taxonomic nodes
#' @description Searches NCBI taxonomy and returns number of descendants
#' taxonomic nodes (species, genera ...) of ID.
#' @param txid Taxonomic ID
#' @template ps
#' @family run-private
#' @return integer
txnds_count <- function(txid, ps) {
  trm <- paste0('txid', txid, '[Subtree]')
  args <- list(db = 'taxonomy', retmax = 0, term = trm)
  res <- search_and_cache(func = rentrez::entrez_search,
                          args = args, fnm = 'search',
                          ps = ps)
  res[['count']]
}

#' @name sqs_count
#' @title Count number of sequences for txid
#' @description Return the number of sequences associated with a
#' taxonomic ID on NCBI GenBank.
#' @param txid Taxonomic ID
#' @template ps
#' @param direct Node-level only or subtree as well? Default FALSE.
#' @family run-private
#' @return integer
sqs_count <- function(txid, ps, direct=FALSE) {
  trm <- searchterm_gen(txid = txid, ps = ps, direct = direct)
  args <- list(db = 'nucleotide', retmax = 0, term = trm)
  res <- search_and_cache(func = rentrez::entrez_search,
                          args = args, fnm = 'search',
                          ps = ps)
  res[['count']]
}

#' @title Return random set of sequence IDs
#' @description For a given txid return a random set of 
#' sequences associated.
#' @details For model organisms downloading all IDs can a take long
#' time or even cause an xml parsing error. For any search with more
#' than hrdmx sequences, this function we will run multiple small
#' searches downloading retmax seq IDs at a time with different
#' retstart values to generate a semi-random vector of sequence IDs.
#' For all other searches, all IDs will be retrieved. Note, it makes
#' no sense for mdlthrs in parameters to be greater than hrdmx in this
#' function.
#' @param txid NCBI taxon identifier
#' @param direct Node-level only or subtree as well? Default FALSE.
#' @template ps
#' @param retmax Maximum number of sequences when querying model
#' organisms. The smaller the more random, the larger the faster.
#' @param hrdmx Absolute maximum number of sequence IDs to download
#' in a single query.
#' @return vector of IDs
#' @family run-private
sids_get <- function(txid, direct, ps, retmax=100, hrdmx=100000) {
  if (sids_check(wd = ps[['wd']], txid = txid)) {
    return(sids_load(wd = ps[['wd']], txid = txid))
  }
  # search
  trm <- searchterm_gen(txid = txid, ps = ps, direct = direct)
  args <- list(db = 'nucleotide', retmax = 0, term = trm,
               use_history = TRUE)
  srch <- safely_connect(func = rentrez::entrez_search, args = args,
                         fnm = 'search', ps = ps)
  nsqs <- srch[['count']]
  if (nsqs == 0) {
    return(NULL)
  }
  ids <- NULL
  if (nsqs > hrdmx) {
    ret_strts <- sample(seq(1, (nsqs - retmax), retmax),
                        round(hrdmx/retmax), replace = FALSE)
  } else {
    if (nsqs <= retmax) {
      ret_strts <- 1
    } else {
      ret_strts <- seq(1, (nsqs - retmax), retmax)
    }
  }
  rsrch <- FALSE
  for (ret_strt in ret_strts) {
    mxtry <- 3
    for (i in 1:mxtry) {
      # use webobject, search without cache
      if (rsrch) {
        srch <- safely_connect(func = rentrez::entrez_search, args = args,
                               fnm = 'search', ps = ps)
      }
      ftch_args <- list(db = 'nucleotide', web_history = srch[['web_history']],
                        retmax = retmax, rettype = 'acc', retstart = ret_strt)
      id_ftch <- try(do.call(what = rentrez::entrez_fetch, args = ftch_args),
                     silent = TRUE)
      if (inherits(id_ftch, 'try-error')) {
        rsrch <- TRUE
        if (i == mxtry) {
          error(ps = ps, 'Failed to get IDs for [', txid, ']')
        }
      }
      rsrch <- FALSE
      break
    }
    ids <- c(ids, strsplit(id_ftch, '\n')[[1]])
  }
  sids_save(wd = ps[['wd']], txid = txid, sids = ids)
  ids
}
