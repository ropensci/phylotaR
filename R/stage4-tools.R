#' @name seeds_blast
#' @title BLAST seed sequences
#' @description Runs all-v-all blast for seed sequences.
#' @param sqs All seed sequences to be BLASTed
#' @template ps
#' @family run-private
#' @return blast res data.frame
seeds_blast <- function(sqs, ps) {
  info(lvl = 2, ps = ps, "BLASTing [", length(sqs@ids), " sqs]")
  dbfl <- 'seeds-db.fa'
  outfl <- 'seeds-db-blastout.txt'
  file.path(ps[['wd']], 'blast', dbfl)
  blastdb_gen(sqs, dbfl = dbfl, ps = ps)
  blast_res <- blastn_run(dbfl = dbfl, outfl = outfl, ps = ps)
  blast_res
}

#' @name clstrs_join
#' @title Join clusters for merging
#' @description Uses seed sequence BLAST results and IDs to join clusters
#' identified as sisters into single clusters. Resulting object is of joined
#' clusters, merging is required to reformat the clusters for subsequent
#' analysis.
#' @param blast_res Seed sequence BLAST results
#' @param seed_ids Seed sequence IDs
#' @param all_clstrs List of all clusters
#' @template ps
#' @family run-private
#' @return list of joined clusters
clstrs_join <- function(blast_res, seed_ids, all_clstrs, ps) {
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
  pull <- blast_res[['query.id']] != blast_res[['subject.id']] &
    blast_res[['qcovs']] > ps[['mncvrg']]
  blast_res <- blast_res[pull, ]
  clstr_list <- blast_clstr(blast_res = blast_res)
  info(lvl = 2, ps = ps, "Identified [", length(clstr_list), "] clusters")
  lapply(clstr_list, join)
}

#' @name clstrs_merge
#' @title Merge joined clusters
#' @description Takes a list of joined clusters and computes each
#' data slot to create a single merged cluster. txdct is required for
#' parent look-up.
#' @param jnd_clstrs List of joined clusters
#' @param txdct Taxonomic dictionary
#' @return list of ClstrRecs
#' @family run-private
clstrs_merge <- function(jnd_clstrs, txdct) {
  mrg_clstrs <- vector('list', length = length(jnd_clstrs))
  for (i in seq_along(jnd_clstrs)) {
    cl <- jnd_clstrs[[i]]
    prnt <- parent_get(id = cl[['txids']], txdct = txdct)
    nsqs <- length(cl[['sids']])
    ntx <- length(unique(cl[['txids']]))
    clstrrec <- new('ClstrRec', sids = cl[['sids']],
                      txids = cl[['txids']], nsqs = nsqs,
                      ntx = ntx, typ = 'merged', prnt = prnt,
                      seed = cl[['seed']])
    mrg_clstrs[[i]] <- clstrrec
  }
  mrg_clstrs
}

#' @name clstrs_renumber
#' @title Renumber cluster IDs
#' @description Returns a ClstrArc with ID determined by the number
#' of sequences in each cluster.
#' @param clstrrecs List of clusters
#' @return ClstrArc
#' @family run-private
clstrs_renumber <- function(clstrrecs) {
  nsqs <- vapply(clstrrecs, function(x) length(x@sids), numeric(1))
  ord <- order(nsqs, decreasing = TRUE)
  clstrrecs <- clstrrecs[ord]
  for (i in seq_along(clstrrecs)) {
    clstrrecs[[i]]@id <- as.integer(i - 1)
  }
  clstrarc_gen(clstrrecs)
}
