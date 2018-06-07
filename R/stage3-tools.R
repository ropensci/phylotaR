#' @name clusterAll
#' @title Hierarchically cluster all sequences of a txid
#' @description Identifies all direct and subtree clusters
#' for a taxonomic ID.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param txdct PhyLoTa table
#' @param ps Parameters
#' @param lvl Log level
cluster_all <- function(txid, sqs, txdct, ps, lvl=0) {
  dds <- descendants_get(id = txid, txdct = txdct, direct = TRUE)
  all_clusters <- cluster_subtree(txid = txid, sqs = sqs, txdct = txdct,
                                ps = ps, dds = dds, lvl = lvl + 1)
  for (dd in dds) {
    info(lvl = lvl + 2, ps = ps, "Processing [id ", txid,
         "] child [id ", dd, "]")
    dd_clusters <- cluster_all(txid = dd, txdct = txdct,
                             sqs = sqs, ps = ps, lvl = lvl + 1)
    all_clusters <- clusterarc_join(all_clusters, dd_clusters)
  }
  all_clusters
}

#' @name clusterSbtr
#' @title Cluster all sequences descending from a txid
#' @description Identifies clusters from sequences associated
#' with a txid and all its descendants. Clusters returned by
#' this function will thus be of cl_type 'subtree'.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param txdct Taxonomic dictionary
#' @param dds Vector of direct descendants
#' @param ps Parameters
#' @param lvl Log level
cluster_subtree <- function(txid, sqs, txdct, dds, ps, lvl) {
  all_clusters <- clusterarc_gen(list())
  rnk <- rank_get(txid = txid, txdct = txdct)
  info(lvl = lvl + 1, ps = ps, "Generating subtree clusters for [id ",
       txid, "(", rnk, ")]")
  if (length(dds) > 0) {
    drct_clusters <- cluster_direct(txid, ps = ps, txdct = txdct,
                                  sqs = sqs, lvl = lvl)
    all_clusters <- clusterarc_join(all_clusters, drct_clusters)
  }
  txids <- descendants_get(id = txid, txdct = txdct, direct = FALSE)
  all_sq_txids <- sqs@txids
  sids <- sqs@ids[which(all_sq_txids %in% as.character(txids))]
  if (length(sids) < 3) {
    info(lvl = lvl + 3, ps = ps, "[", length(sids), " sqs]",
         " -- too few sequences, cannot make clusters")
    return(all_clusters)
  }
  sqs_prt <- sqs[sids]
  sbtr_clusters <- cluster_sqs(txid = txid, sqs = sqs_prt,
                             typ = 'subtree', ps = ps, lvl = lvl)
  clusterarc_join(all_clusters, sbtr_clusters)
}

#' @name clusterDrct
#' @title Cluster sequences directly associated with txid
#' @description In GenBank certain sequences may only be associated
#' with a higher level taxon (e.g. genus, family ...). This function
#' generates clusters from these sequences, alone. This function
#' identifies such sequences in the sequence object and generates
#' a list of clusters of cl_type 'Node'.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param txdct PhyLoTa table
#' @param ps Parameters
#' @param lvl Log level
cluster_direct <- function(txid, sqs, txdct, ps, lvl) {
  all_clusters <- clusterarc_gen(list())
  rnk <- rank_get(txid = txid, txdct = txdct)
  info(lvl = lvl + 1, ps = ps, "Generating direct clusters for [id ",
       txid, "(", rnk, ")]")
  all_sq_txids <- sqs@txids
  sids <- sqs@ids[all_sq_txids %in% as.character(txid)]
  info(lvl = lvl + 2, ps = ps, "[", length(sids), " sqs]")
  if (length(sids) < 3) {
    info(lvl = lvl + 3, ps = ps, "Too few sequences, cannot make clusters")
    return(all_clusters)
  }
  sqs_prt <- sqs[sids]
  cluster_sqs(txid = txid, sqs = sqs_prt, typ = 'direct',
              ps = ps, lvl = lvl)
}

#' @name clusterSqs
#' @title Identify clusters from sequences
#' @description Given a sequence object, this function will generate
#' a list of cluster objects using BLAST
#' @param txid Taxonomic ID
#' @param sqs Sequence object of sequences to be BLASTed
#' @param ps Parameters
#' @param typ Direct or Subtree?
#' @param lvl Log level
cluster_sqs <- function(txid, sqs, ps, lvl,
                         typ=c('direct', 'subtree', 'paraphyly')) {
  typ <- match.arg(typ)
  info(lvl = lvl + 1, ps = ps, "BLASTing [", length(sqs@ids),
       " sqs] ....")
  blast_res <- blast_sqs(txid = txid, typ = typ, sqs = sqs, ps = ps, lvl = lvl)
  if (is.null(blast_res)) {
    return(NULL)
  }
  cluster_list <- blast_cluster(blast_res = blast_res)
  cl_rcrds <- clusterrec_gen(cluster_list = cluster_list, txid = txid, sqs = sqs, typ = typ)
  info(lvl = lvl + 1, ps = ps, "Identified [", length(cl_rcrds@ids), "] clusters")
  cl_rcrds
}

#' @name blstSqs
#' @title BLAST All vs All
#' @description Return blast results from
#' BLASTing all vs all for given sequences
#' @param txid Taxonomic node ID, numeric
#' @param typ Cluster type, 'direct' or 'subtree'
#' @param sqs Sequences
#' @param lvl Log level
#' @param ps Parameters
blast_sqs <- function(txid, typ, sqs, ps, lvl) {
  blast_res <- ldBlstCch(sqs@ids, wd = ps[['wd']])
  if (is.null(blast_res)) {
    dbfl <- paste0('taxon-', txid, '-typ-', typ,
                   '-db.fa')
    outfl <- paste0('taxon-', txid, '-typ-', typ,
                    '-blastout.txt')
    blastdb_gen(sqs = sqs, dbfl = dbfl, ps = ps)
    blast_res <- blastn_run(dbfl = dbfl, outfl = outfl, ps = ps)
    if (is.null(blast_res)) {
      blast_res <- NA
    }
    svBlstCch(sqs@ids, wd = ps[['wd']], obj = blast_res)
  }
  # TODO: Not so elegant
  if (any(is.na(blast_res))) {
    return(NULL)
  }
  blast_filter(blast_res = blast_res, ps = ps)
}

#' @name clusterBlstRs
#' @title Cluster BLAST Results
#' @description Find single-linkage clusters from
#' BLAST results. Identifies seed sequence.
#' @return List of list
#' @param blst_rs BLAST results
blast_cluster <- function(blast_res) {
  g <- igraph::graph.data.frame(blast_res[ ,c("query.id", "subject.id")],
                                directed = FALSE)
  clusters <- igraph::clusters(g)
  clusters <- clusters[['membership']]
  cluster_list <- lapply(unique(clusters), function(x) {
    list('sids' = sort(names(clusters)[which(clusters == x)]))
  })
  degrees <- igraph::degree(g)
  cluster_list <- lapply(cluster_list, function(cl){
    idx <- order(degrees[cl[['sids']]], decreasing = TRUE)[1]
    # index of most connected component
    cl[['seed']] <- cl[['sids']][idx]
    cl
  })
  cluster_list
}

#' @name genClRcrds
#' @title Generate list of clusters
#' @description Takes a list of lists of cluster descriptions,
#' returns a ClRcrdBx.
#' @param cluster_lst List of list of cluster descriptions
#' @param txid Taxonomic node ID
#' @param sqs Sequnece records
#' @param typ Subtree of direct?
clusterrec_gen <- function(cluster_list, txid, sqs, typ) {
  cl_rcrds <- vector('list', length = length(cluster_list))
  for (i in seq_along(cluster_list)) {
    cl <- cluster_list[[i]]
    cl_sqs <- sqs[cl[['sids']]]
    nsqs <- length(cl[['sids']])
    ntx <- length(unique(cl_sqs@txids))
    cl_rcrd <- new('ClusterRec', sids = cl[['sids']],
                   txids = cl_sqs@txids, nsqs = nsqs,
                   ntx = ntx, typ = typ,
                   prnt = as.character(txid),
                   seed = cl[['seed']])
    cl_rcrds[[i]] <- cl_rcrd
  }
  clusterarc_gen(cl_rcrds)
}

#' @name archive_cluster_gen
#' @title Generate ArchiveCluster container class
#' @description Takes a list of RecordCluster classes, returns an ArchiveCluster.
#' @param cluster_rcrds list of RecordCluster classes
#' @return ArchiveCluster
#' @noRd
clusterarc_gen <- function(cluster_rcrds) {
  ids <- as.character(seq_along(cluster_rcrds) - 1)
  names(cluster_rcrds) <- ids
  new('ClusterArc', ids = ids, cls = cluster_rcrds)
}

#' @name archive_cluster_join
#' @title Join two ArchiveCluster classes
#' @description Take two ArchiveCluster classes and join them into
#' a single ArchiveCluster
#' @param ac_1 ArchiveCluster
#' @param ac_2 ArchiveCluster
#' @return ArchiveCluster
#' @noRd
clusterarc_join <- function(ac_1, ac_2) {
  clusterarc_gen(c(ac_1@cls, ac_2@cls))
}