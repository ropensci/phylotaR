#' @name hierarchic_download
#' @title Hierarchically get sequences for a txid
#' @description Looks up and downloads sequences for a taxonomic ID.
#' @param txid Taxonomic node ID, numeric
#' @param txdct Taxonomic dictionary
#' @param ps parameters
#' @param lvl Log level
#' @noRd
#' @return Vector of SeqRecs
hierarchic_download <- function(txid, txdct, ps, lvl=0) {
  # get subtree counts if that is smaller than ps[['mdlthrs']]
  # or if there are no direct descendants
  dds <- descendants_get(id = txid, txdct = txdct, direct = TRUE)
  subtree_count <- sqs_count(txid, drct = FALSE, ps = ps)
  if (subtree_count <= ps[['mdlthrs']] | length(dds) == 0) {
    info(lvl = 2 + lvl, ps = ps, "+ whole subtree ...")
    sqs <- seqrec_get(txid = txid, drct = FALSE, ps = ps, lvl = lvl)
    return(sqs)
  }
  # 1st direct sqs from focal taxon, then from DDs
  info(lvl = 2 + lvl, ps = ps, "+ direct ...")
  sqs <- seqrec_get(txid = txid, drct = TRUE, ps = ps, lvl = lvl)
  info(lvl = 2 + lvl, ps = ps, "+ by child ...")
  for (dd in dds) {
    lvl <- lvl + 1
    info(lvl = 2 + lvl, ps = ps, "Working on child [id ", dd,"]")
    sqs <- c(sqs, hierarchic_download(txid = dd, txdct = txdct,
                              ps = ps, lvl = lvl))
  }
  sqs
}

#' @name seqrec_augment
#' @title Augment sequence records list
#' @description Add taxids to records and convert to archive.
#' @param sqs List of SqRcrds
#' @param txdct Taxonomic Dictionary
#' @return SeqArc
#' @noRd
seqrec_augment <- function(sqs, txdct) {
  # txids are not downloaded as part of sequence, added here
  txdct_nms <- vapply(txdct@rcrds, function(x) x@scnm, '')
  txdct_ids <- vapply(txdct@rcrds, function(x) x@id, '')
  sqs_nms <- vapply(sqs, function(x) x@orgnsm, '')
  sqs_ids <- txdct_ids[match(sqs_nms, txdct_nms)]
  for (i in seq_along(sqs)) {
    sqs[[i]]@txid <- as.character(sqs_ids[[i]])
  }
  seqarc_gen(sqs)
}

#' @title seqrec_get
#' @description Downloads sequences from GenBank in batches.
#' @param txid NCBI taxonomic ID
#' @param drct Node-level only or subtree as well? Default FALSE.
#' @param ps parameters
#' @param lvl Log level
#' @return Vector of sequence records
#' @noRd
seqrec_get <- function(txid, ps, drct=FALSE, lvl=0) {
  # test w/ golden moles 9389
  downloader <- function(ids, ps) {
    ftch_args <- list(db = "nucleotide",
                      rettype = 'gbwithparts',
                      retmode = 'xml', id = ids)
    search_and_cache(func = rentrez::entrez_fetch,
                     args = ftch_args, fnm = 'fetch',
                     ps = ps)
  }
  # get accessions
  sids <- sids_get(txid = txid, drct = drct, ps = ps)
  if (length(sids) < 1) {
    return(list())
  }
  if (length(sids) > ps[['mdlthrs']]) {
    info(lvl = lvl + 3, ps = ps, "More than [", ps[['mdlthrs']],
         ' sqs] available. Choosing at random.')
    sids <- sids[sample(1:length(sids), size = ps[['mdlthrs']],
                        replace = FALSE)]
  }
  info(lvl = lvl + 3, ps = ps, "Getting [", length(sids), " sqs] ...")
  # return whole GB record
  raw_recs <- batcher(sids, func = downloader, ps = ps, lvl = lvl + 4)
  #  split up by feature, return SeqRecs
  seqrecs <- numeric()
  for (i in seq_along(raw_recs)) {
    seqrecs <- c(seqrec_convert(raw_recs = raw_recs[[i]], ps = ps),
                 seqrecs)
  }
  seqrecs
}
