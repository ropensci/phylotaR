#' @name blstSeeds
#' @title BLAST seed sequences
#' @description Runs all-v-all blast for seed sequences.
#' @param sqs All seed sequences to be BLASTed
#' @param ps Parameters
blstSeeds <- function(sqs, ps) {
  info(lvl=2, ps=ps, "BLASTing [", length(sqs@ids),
       " sqs]")
  dbfl <- 'seeds-db.fa'
  outfl <- 'seeds-db-blastout.txt'
  file.path(ps[['wd']], 'blast', dbfl)
  mkBlstDB(sqs, dbfl=dbfl, ps=ps)
  blst_rs <- blstN(dbfl=dbfl, outfl=outfl, ps=ps)
  blst_rs
}

#' @name jnClstrs
#' @title Join clusters for merging
#' @description Uses seed sequence BLAST results and IDs
#' to join clusters identifed as sisters into single clusters.
#' Resulting object is of joined clusters, merging is required
#' to reformat the clusters for subsequent analysis.
#' @param blst_rs Seed sequence BLAST results
#' @param seed_ids Seed sequence IDs
#' @param all_clstrs List of all clusters
#' @param ps Parameters
jnClstrs <- function(blst_rs, seed_ids, all_clstrs, ps) {
  join <- function(x) {
    pull <- seed_ids %in% x[['sids']]
    jnd_clstr <- all_clstrs[pull]
    nms <- slotNames(jnd_clstr[[1]])
    clstr <- lapply(nms, function(nm)
      unlist(lapply(jnd_clstr, function(cl) slot(cl, nm))))
    names(clstr) <- nms
    clstr[['typ']] <- 'merged'
    clstr[['seed']] <- x[['seed']]
    # ensure no dups seqs in joined cluster
    pull <- !duplicated(clstr[['sids']])
    clstr[['sids']] <- clstr[['sids']][pull]
    clstr[['txids']] <- clstr[['txids']][pull]
    clstr
  }
  pull <- blst_rs[['query.id']] != blst_rs[['subject.id']] &
    blst_rs[['qcovs']] > ps[['mncvrg']]
  blst_rs <- blst_rs[pull, ]
  clstr_lst <- clstrBlstRs(blst_rs=blst_rs)
  info(lvl=2, ps=ps, "Identified [", length(clstr_lst),
       "] clusters")
  lapply(clstr_lst, join)
}

#' @name mrgClstrs
#' @title Merge joined clusters
#' @description Takes a list of joined clusters and computes
#' each data slot to create a single merged cluster.
#' @param jnd_clstrs List of joined clusters
mrgClstrs <- function(jnd_clstrs) {
  mrg_clstrs <- vector('list', length=length(jnd_clstrs))
  for(i in seq_along(jnd_clstrs)) {
    # TODO: look up parent
    cl <- jnd_clstrs[[i]]
    clstr <- new('Clstr', sids=cl[['sids']],
                 txids=cl[['txids']], typ='merged',
                 seed=cl[['seed']])
    mrg_clstrs[[i]] <- clstr
  }
  mrg_clstrs
}

#' @name rnmbrClstrs
#' @title Renumbers cluster IDs
#' @description Returns a list of clusters with
#' ID determined by the number of sequences
#' in each cluster.
#' @param clstrs List of clusters
rnmbrClstrs <- function(clstrs) {
  nsqs <- vapply(clstrs, function(x) length(x@sids),
                 numeric(1))
  ord <- order(nsqs, decreasing=TRUE)
  clstrs <- clstrs[ord]
  for(i in seq_along(clstrs)) {
    clstrs[[i]]@id <- as.integer(i - 1)
  }
  clstrs
}
