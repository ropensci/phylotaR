#' @name clstr_all
#' @title Hierarchically cluster all sequences of a txid
#' @description Identifies all direct and subtree clusters for a taxonomic ID.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param txdct Taxonomic dictionary
#' @template ps
#' @template lvl
#' @return ClstrArc
#' @family run-private
clstr_all <- function(txid, sqs, txdct, ps, lvl=0) {
  dds <- descendants_get(id = txid, txdct = txdct, direct = TRUE)
  all_clstrs <- clstr_subtree(txid = txid, sqs = sqs, txdct = txdct, ps = ps,
                              dds = dds, lvl = lvl + 1)
  for (dd in dds) {
    info(lvl = lvl + 2, ps = ps, "Processing [id ", txid, "] child [id ", dd,
         "]")
    dd_clstrs <- clstr_all(txid = dd, txdct = txdct, sqs = sqs, ps = ps,
                           lvl = lvl + 1)
    all_clstrs <- clstrarc_join(all_clstrs, dd_clstrs)
  }
  all_clstrs
}

#' @name clstr_subtree
#' @family run-private
#' @title Cluster all sequences descending from a txid
#' @description Identifies clusters from sequences associated with a
#' txid and all its descendants. Clusters returned by this function
#' will thus be of cl_type 'subtree'.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param txdct Taxonomic dictionary
#' @param dds Vector of direct descendants
#' @template ps
#' @template lvl
#' @return ClstrArc
clstr_subtree <- function(txid, sqs, txdct, dds, ps, lvl) {
  all_clstrs <- clstrarc_gen(list())
  rnk <- rank_get(txid = txid, txdct = txdct)
  info(lvl = lvl + 1, ps = ps, "Generating subtree clusters for [id ", txid,
       " (", rnk, ")]")
  if (length(dds) > 0) {
    drct_clstrs <- clstr_direct(txid, ps = ps, txdct = txdct,
                                sqs = sqs, lvl = lvl)
    all_clstrs <- clstrarc_join(all_clstrs, drct_clstrs)
  }
  txids <- descendants_get(id = txid, txdct = txdct, direct = FALSE)
  all_sq_txids <- sqs@txids
  sids <- sqs@ids[which(all_sq_txids %in% as.character(txids))]
  if (length(sids) < 3) {
    # info(lvl = lvl + 3, ps = ps, "[", length(sids), " sqs]",
    #      " -- too few sequences, cannot make clusters")
    return(all_clstrs)
  }
  sqs_prt <- sqs[sids]
  sbtr_clstrs <- clstr_sqs(txid = txid, sqs = sqs_prt, typ = 'subtree', ps = ps,
                           lvl = lvl)
  clstrarc_join(all_clstrs, sbtr_clstrs)
}

#' @name clstr_direct
#' @family run-private
#' @title Cluster sequences directly associated with txid
#' @description In GenBank certain sequences may only be associated
#' with a higher level taxon (e.g. genus, family ...). This function
#' generates clusters from these sequences, alone. This function
#' identifies such sequences in the sequence object and generates
#' a list of clusters of cl_type 'direct'.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param txdct Taxonomic dictionary
#' @template ps
#' @template lvl
#' @return ClstrArc
clstr_direct <- function(txid, sqs, txdct, ps, lvl) {
  all_clstrs <- clstrarc_gen(list())
  rnk <- rank_get(txid = txid, txdct = txdct)
  info(lvl = lvl + 1, ps = ps, "Generating direct clusters for [id ", txid,
       "(", rnk, ")]")
  all_sq_txids <- sqs@txids
  sids <- sqs@ids[all_sq_txids %in% as.character(txid)]
  info(lvl = lvl + 2, ps = ps, "[", length(sids), " sqs]")
  if (length(sids) < 3) {
    info(lvl = lvl + 3, ps = ps, "Too few sequences, cannot make clusters")
    return(all_clstrs)
  }
  sqs_prt <- sqs[sids]
  clstr_sqs(txid = txid, sqs = sqs_prt, typ = 'direct', ps = ps, lvl = lvl)
}

#' @name clstr_sqs
#' @title Identify clusters from sequences
#' @description Given a sequence object, this function will generate
#' a list of cluster objects using BLAST
#' @param txid Taxonomic ID
#' @param sqs Sequence object of sequences to be BLASTed
#' @template ps
#' @param typ Direct, subtree or paraphyly?
#' @template lvl
#' @family run-private
clstr_sqs <- function(txid, sqs, ps, lvl,
                      typ=c('direct', 'subtree', 'paraphyly')) {
  typ <- match.arg(typ)
  info(lvl = lvl + 1, ps = ps, "BLASTing [", length(sqs@ids), " sqs] ....")
  blast_res <- blast_sqs(txid = txid, typ = typ, sqs = sqs, ps = ps, lvl = lvl)
  if (is.null(blast_res)) {
    return(NULL)
  }
  clstr_list <- blast_clstr(blast_res = blast_res)
  clstrrecs <- clstrrec_gen(clstr_list = clstr_list, txid = txid, sqs = sqs,
                            typ = typ)
  info(lvl = lvl + 1, ps = ps, "Identified [", length(clstrrecs@ids),
       "] clusters")
  clstrrecs
}

# TODO: update this function to make it generic
#' @name blast_sqs
#' @title BLAST All vs All
#' @description Return BLAST results from BLASTing all vs all for
#' given sequences. Returns NULL if no BLAST results generated.
#' @param txid Taxonomic node ID, numeric
#' @param typ Cluster type, 'direct' or 'subtree'
#' @param sqs Sequences
#' @template lvl
#' @template ps
#' @family run-private
#' @return blast_res data.frame or NULL
blast_sqs <- function(txid, typ, sqs, ps, lvl) {
  blast_res <- blastcache_load(sqs@ids, wd = ps[['wd']])
  if (is.null(blast_res)) {
    dbfl <- paste0('taxon-', txid, '-typ-', typ, '-db.fa')
    outfl <- paste0('taxon-', txid, '-typ-', typ, '-blastout.txt')
    blastdb_gen(sqs = sqs, dbfl = dbfl, ps = ps)
    blast_res <- blastn_run(dbfl = dbfl, outfl = outfl, ps = ps)
    if (is.null(blast_res)) {
      blast_res <- NA
    }
    blastcache_save(sqs@ids, wd = ps[['wd']], obj = blast_res)
  }
  # TODO: Not so elegant
  if (any(is.na(blast_res))) {
    return(NULL)
  }
  blast_filter(blast_res = blast_res, ps = ps)
}

#' @name blast_clstr
#' @title Cluster BLAST Results
#' @description Find single-linkage clusters from BLAST results.
#' Identifies seed sequence.
#' @return List of list
#' @param blast_res BLAST results
#' @return list of cluster descriptions
#' @family run-private
blast_clstr <- function(blast_res) {
  g <- igraph::graph.data.frame(blast_res[ ,c("query.id", "subject.id")],
                                directed = FALSE)
  clstrs <- igraph::clusters(g)
  clstrs <- clstrs[['membership']]
  clstr_list <- lapply(unique(clstrs), function(x) {
    list('sids' = sort(names(clstrs)[which(clstrs == x)]))
  })
  degrees <- igraph::degree(g)
  clstr_list <- lapply(clstr_list, function(cl){
    idx <- order(degrees[cl[['sids']]], decreasing = TRUE)[1]
    # index of most connected component
    cl[['seed']] <- cl[['sids']][idx]
    cl
  })
  clstr_list
}

#' @name clstrrec_gen
#' @title Generate list of clusters
#' @description Takes a list of lists of cluster descriptions
#' and generates ClstrRecs.
#' @param clstr_list List of list of cluster descriptions
#' @param txid Taxonomic node ID
#' @param sqs Sequence records
#' @param typ Cluster type
#' @family run-private
#' @return list of ClstrRecs
clstrrec_gen <- function(clstr_list, txid, sqs, typ) {
  clstrrecs <- vector('list', length = length(clstr_list))
  for (i in seq_along(clstr_list)) {
    clstr <- clstr_list[[i]]
    clstr_sqs <- sqs[clstr[['sids']]]
    nsqs <- length(clstr[['sids']])
    ntx <- length(unique(clstr_sqs@txids))
    clstrrec <- new('ClstrRec', sids = clstr[['sids']], txids = clstr_sqs@txids,
                    nsqs = nsqs, ntx = ntx, typ = typ,
                    prnt = as.character(txid), seed = clstr[['seed']])
    clstrrecs[[i]] <- clstrrec
  }
  clstrarc_gen(clstrrecs)
}

#' @name clstrarc_gen
#' @title Generate cluster archive container class
#' @description Takes a list of ClstrRecs, returns a ClstrArc.
#' @param clstrrecs list of ClstrRecs
#' @return ClstrArc
#' @family run-private
clstrarc_gen <- function(clstrrecs) {
  ids <- as.character(seq_along(clstrrecs) - 1)
  names(clstrrecs) <- ids
  new('ClstrArc', ids = ids, clstrs = clstrrecs)
}

#' @name clstrarc_join
#' @title Join two cluster archive
#' @description Take two ClstrArc classes and join them into a single
#' ClstrArc.
#' @param clstrarc_1 ClstrArc
#' @param clstrarc_2 ClstrArc
#' @return ClstrArc
#' @family run-private
clstrarc_join <- function(clstrarc_1, clstrarc_2) {
  clstrarc_gen(c(clstrarc_1@clstrs, clstrarc_2@clstrs))
}
