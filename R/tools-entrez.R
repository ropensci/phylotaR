#' @name mkSrchTrm
#' @title Construct GenBank Search Term
#' @description Construct search term for searching
#' GenBank's nucleotide database.
#' Limits the maximum size of sequences, avoids
#' whole genome shotguns, predicted, unverified and
#' synthetic sequences.
#' @param txid Taxonomic ID
#' @param ps Parameter list
#' @param drct Node-level only or subtree as well? Default FALSE.
searchterm_gen <- function(txid, ps, drct=FALSE) {
  org_trm <- ifelse(drct, '[Organism:noexp]',
                     '[Organism:exp]' )
  avd1 <- ' NOT predicted[TI] NOT "whole genome shotgun"[TI]'
  avd2 <- ' NOT unverified[TI] NOT "synthetic construct"[Organism]'
  paste0('(txid', txid, org_trm, ' AND ', ps[['mnsql']],
         ':100000[SLEN])', avd1, avd2)
}


#' @name nNds
#' @title Count number of descending taxonomic nodes
#' @description Searches NCBI taxonomy and returns
#' number of descendents taxonomic nodes (species, genera ...)
#' of ID.
#' @param txid Taxonomic ID
#' @param ps Parameters
txnds_count <- function(txid, ps) {
  trm <- paste0('txid', txid, '[Subtree]')
  args <- list(db='taxonomy', retmax=0, term=trm)
  res <- search_and_cache(func=rentrez::entrez_search,
                          args=args, fnm='search',
                          ps=ps)
  res[['count']]
}

#' @name nSqs
#' @title Count number of sequences for txid
#' @description Return the number of sequences
#' associated with a taxonomic ID on NCBI GenBank.
#' @param txid Taxonomic ID
#' @param ps parameters
#' @param drct Node-level only or subtree as well? Default FALSE.
seqs_count <- function(txid, ps, drct=FALSE) {
  trm <- searchterm_gen(txid=txid, ps=ps, drct=drct)
  args <- list(db='nucleotide', retmax=0, term=trm)
  res <- search_and_cache(func=rentrez::entrez_search,
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
#' @param drct Node-level only or subtree as well? Default FALSE.
#' @param ps Parameters
#' @param sqcnt Sequence count as determined with nSqs()
#' @param retmax Maximum number of sequences when querying model
#' organisms. The smaller the more random, the larger the faster.
#' @param hrdmx Absolute maximum number of sequence IDs to download
#' in a single query.
#' @return vector ot IDs
sids_get <- function(txid, drct, ps, retmax=100,
                     hrdmx=100000) {
  if(chckGIs(wd=ps[['wd']], txid=txid)) {
    return(ldGIs(wd=ps[['wd']], txid=txid))
  }
  # search
  trm <- searchterm_gen(txid=txid, ps=ps, drct=drct)
  args <- list(db='nucleotide', retmax=0, term=trm,
               use_history=TRUE)
  srch <- safely_connect(func=rentrez::entrez_search,
                         args=args, fnm='search',
                         ps=ps)
  nseqs <- srch[['count']]
  ids <- NULL
  if (nseqs > hrdmx) {
    ret_strts <- sample(seq(1, (nseqs-retmax), retmax),
                        round(hrdmx/retmax), replace=FALSE)
  } else {
    ret_strts <- seq(1, (nseqs-retmax), retmax)
  }
  rsrch <- FALSE
  for(ret_strt in ret_strts) {
    mxtry <- 3
    for(i in 1:mxtry) {
      # use webobject, search without cache
      if(rsrch) {
        srch <- safely_connect(func=rentrez::entrez_search,
                               args=args, fnm='search',
                               ps=ps)
      }
      ftch_args <- list(db='nucleotide',
                        web_history=srch[['web_history']],
                        retmax=retmax,
                        rettype='acc',
                        retstart=ret_strt)
      id_ftch <- try(do.call(what=rentrez::entrez_fetch,
                             args=ftch_args), silent=TRUE)
      if(inherits(id_ftch, 'try-error')) {
        rsrch <- TRUE
        if(i == mxtry) {
          error(ps=ps, 'Failed to get IDs for [', txid, ']')
        }
      }
      rsrch <- FALSE
      break
    }
    ids <- c(ids, strsplit(id_ftch, '\n')[[1]])
  }
  svGIs(wd=ps[['wd']], txid=txid, gis=ids)
  ids
}
