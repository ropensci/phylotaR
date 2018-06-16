#' @name drop_sqs
#' @title Drop sequences in a cluster
#' @description Drop all sequences in a cluster except those identified by user.
#' @param phylota Phylota object
#' @param cid Cluster ID
#' @param sid Sequence ID(s) to be kept
#' @return phylota
#' @export
#' @family tools-public
#' @example examples/drop_sqs.R
drop_sqs <- function(phylota, cid, sid) {
  indx <- which(phylota@clstrs@ids == cid)
  clstr <- phylota@clstrs@clstrs[[indx]]
  pull <- clstr@sids %in% sid
  clstr@sids <- clstr@sids[pull]
  clstr@nsqs <- length(clstr@sids)
  clstr@txids <- clstr@txids[pull]
  clstr@ntx <- length(unique(clstr@txids))
  phylota@clstrs@clstrs[[indx]] <- clstr
  update_phylota(phylota)
}

#' @name drop_clstrs
#' @title Drop cluster records from phylota object
#' @description Drops all clusters except those
#' identified by user.
#' @param phylota Phylota object
#' @param cid Cluster ID(s) to be kept
#' @return phylota
#' @export
#' @family tools-public
#' @example examples/drop_clstrs.R
drop_clstrs <- function(phylota, cid) {
  phylota@clstrs <- phylota@clstrs[cid]
  phylota@cids <- cid
  update_phylota(phylota)
}

#' @name drop_by_rank
#' @title Reduce clusters to specific rank
#' @description Identifies higher level taxa for each sequence in clusters for
#' given rank. Selects representative sequences for each unique taxon using the
#' choose_by functions. By default, the function will choose the top ten
#' sequences by first sorting by those with fewest number of ambiguous
#' sequences, then by youngest, then by sequence length.
#' @param phylota Phylota object
#' @param rnk Taxonomic rank
#' @param keep_higher Keep higher taxonomic ranks?
#' @param n Number of sequences per taxon
#' @param choose_by Vector of selection functions
#' @param greatest Greatest of lowest for each choose_by function
#' @return phylota
#' @export
#' @family tools-public
#' @example examples/drop_by_rank.R
drop_by_rank <- function(phylota, rnk = 'species', keep_higher = FALSE, n = 10,
                         choose_by = c('pambgs', 'age', 'nncltds'),
                         greatest = c(FALSE, FALSE, TRUE)) {
  slct <- function(txid) {
    pssbls <- sids[txid == txids]
    for (i in seq_along(choose_by)) {
      vals <- get_sq_slot(phylota = phylota, sid = pssbls,
                          slt_nm = choose_by[[i]])
      names(vals) <- pssbls
      mx_n <- ifelse(length(vals) > n, n, length(vals))
      vals <- sort(x = vals, decreasing = greatest[i])[1:mx_n]
      pssbls <- names(vals)
    }
    pssbls
  }
  pull <- !choose_by %in% list_seqrec_slots()
  if (any(pull)) {
    stop(paste0('[', choose_by[pull], '] not in SeqRec'))
  }
  for (cid in phylota@cids) {
    txids <- get_txids(phylota = phylota, cid = cid, rnk = rnk,
                       keep_higher = keep_higher)
    sids <- phylota@clstrs[[cid]]@sids
    pull <- txids != ''
    sids <- sids[pull]
    txids <- txids[pull]
    unqids <- unique(txids)
    keep <- unlist(lapply(unqids, slct))
    phylota <- drop_sqs(phylota = phylota, cid = cid, sid = keep)
  }
  phylota
}
