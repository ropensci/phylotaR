#' @name blstSeeds
#' @title BLAST seed sequences
#' @description Runs all-v-all blast for seed sequences.
#' @param sqs All seed sequences to be BLASTed
#' @param ps Parameters
seeds_blast <- function(sqs, ps) {
  info(lvl=2, ps=ps, "BLASTing [", length(sqs@ids),
       " sqs]")
  dbfl <- 'seeds-db.fa'
  outfl <- 'seeds-db-blastout.txt'
  file.path(ps[['wd']], 'blast', dbfl)
  blastdb_gen(sqs, dbfl=dbfl, ps=ps)
  blst_res <- blastn_run(dbfl=dbfl, outfl=outfl, ps=ps)
  blst_res
}

#' @name jnclusters
#' @title Join clusters for merging
#' @description Uses seed sequence BLAST results and IDs
#' to join clusters identifed as sisters into single clusters.
#' Resulting object is of joined clusters, merging is required
#' to reformat the clusters for subsequent analysis.
#' @param blst_rs Seed sequence BLAST results
#' @param seed_ids Seed sequence IDs
#' @param all_clusters List of all clusters
#' @param ps Parameters
clusters_join <- function(blast_res, seed_ids, all_clusters, ps) {
  join <- function(x) {
    pull <- seed_ids %in% x[['sids']]
    jnd_cluster <- all_clusters[pull]
    nms <- slotNames(jnd_cluster[[1]])
    cluster <- lapply(nms, function(nm)
      unlist(lapply(jnd_cluster, function(cl) slot(cl, nm))))
    names(cluster) <- nms
    cluster[['typ']] <- 'merged'
    cluster[['seed']] <- x[['seed']]
    # ensure no dups seqs in joined cluster
    pull <- !duplicated(cluster[['sids']])
    cluster[['sids']] <- cluster[['sids']][pull]
    cluster[['txids']] <- cluster[['txids']][pull]
    cluster
  }
  pull <- blast_res[['query.id']] != blast_res[['subject.id']] &
    blast_res[['qcovs']] > ps[['mncvrg']]
  blast_res <- blast_res[pull, ]
  cluster_list <- blast_cluster(blast_res = blast_res)
  info(lvl = 2, ps = ps, "Identified [", length(cluster_list),
       "] clusters")
  lapply(cluster_list, join)
}

#' @name mrgclusters
#' @title Merge joined clusters
#' @description Takes a list of joined clusters and computes
#' each data slot to create a single merged cluster.
#' txdct is required for parent look-up.
#' @param jnd_clusters List of joined cluster records
#' @param txdct Taxonomic dictionary
clusters_merge <- function(jnd_clusters, txdct) {
  mrg_clusters <- vector('list', length = length(jnd_clusters))
  for (i in seq_along(jnd_clusters)) {
    cl <- jnd_clusters[[i]]
    prnt <- parent_get(id = cl[['txids']], txdct = txdct)
    nsqs <- length(cl[['sids']])
    ntx <- length(unique(cl[['txids']]))
    cl_rcrd <- new('ClusterRec', sids = cl[['sids']],
                   txids = cl[['txids']], nsqs = nsqs,
                   ntx = ntx, typ = 'merged',
                   prnt = prnt, seed = cl[['seed']])
    mrg_clusters[[i]] <- cl_rcrd
  }
  mrg_clusters
}

#' @name rnmbrclusters
#' @title Renumbers cluster IDs
#' @description Returns a ClRcrdBx with
#' ID determined by the number of sequences
#' in each cluster.
#' @param cluster_rcrds List of clusters
clusters_renumber <- function(cluster_rcrds) {
  nsqs <- vapply(cluster_rcrds, function(x) length(x@sids),
                 numeric(1))
  ord <- order(nsqs, decreasing=TRUE)
  cluster_rcrds <- cluster_rcrds[ord]
  for (i in seq_along(cluster_rcrds)) {
    cluster_rcrds[[i]]@id <- as.integer(i - 1)
  }
  clusterarc_gen(cluster_rcrds)
}
