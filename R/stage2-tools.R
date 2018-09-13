#' @name hierarchic_download
#' @title Hierarchically get sequences for a txid
#' @description Looks up and downloads sequences for a taxonomic ID.
#' @param txid Taxonomic node ID, numeric
#' @param txdct Taxonomic dictionary
#' @template ps
#' @template lvl
#' @family run-private
#' @return Vector of SeqRecs
hierarchic_download <- function(txid, txdct, ps, lvl=0) {
  # get subtree counts if that is smaller than ps[['mdlthrs']]
  # or if there are no direct descendants
  dds <- descendants_get(id = txid, txdct = txdct, direct = TRUE)
  subtree_count <- sqs_count(txid, direct = FALSE, ps = ps)
  if (subtree_count <= ps[['mdlthrs']] | length(dds) == 0) {
    info(lvl = 2 + lvl, ps = ps, "+ whole subtree ...")
    sqs <- seqrec_get(txid = txid, direct = FALSE, ps = ps, lvl = lvl)
    return(sqs)
  }
  # 1st direct sqs from focal taxon, then from DDs
  info(lvl = 2 + lvl, ps = ps, "+ direct ...")
  sqs <- seqrec_get(txid = txid, direct = TRUE, ps = ps, lvl = lvl)
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
#' @param sqs List of SeqRecs
#' @param txdct Taxonomic Dictionary
#' @return SeqArc
#' @family run-private
seqrec_augment <- function(sqs, txdct) {
  # txids are not downloaded as part of sequence, added here
  # hopefully, restez integration will make this function redundant
  txdct_nms <- vapply(txdct@recs, function(x) x@scnm, '')
  txdct_ids <- vapply(txdct@recs, function(x) x@id, '')
  sqs_nms <- vapply(sqs, function(x) x@orgnsm, '')
  sqs_ids <- txdct_ids[match(sqs_nms, txdct_nms)]
  # tax name may be missing from taxdict
  pull <- !is.na(sqs_ids)
  sqs <- sqs[pull]
  sqs_ids <- sqs_ids[pull]
  # Look up with entrez if missing -- needs developing
  # if (any(pull)) {
  #   sq_smmrs <- rentrez::entrez_fetch(db = 'nucleotide', id = sqs_ids[pull])
  #   sqs_ids <- vapply(X = sq_smmrs, FUN = '[[', FUN.VALUE = character(1),
  #                     i = 'taxid')
  # }
  for (i in seq_along(sqs)) {
    sqs[[i]]@txid <- as.character(sqs_ids[[i]])
  }
  seqarc_gen(sqs)
}

#' @title seqrec_get
#' @description Downloads sequences from GenBank in batches.
#' @details If a restez database is available and the number of sequences to
#' retrieve is less than 'btchsz', the function will look the sequences up
#' from the database rather than download.
#' @param txid NCBI taxonomic ID
#' @param direct Node-level only or subtree as well? Default FALSE.
#' @template ps
#' @template lvl
#' @return Vector of sequence records
#' @family run-private
seqrec_get <- function(txid, ps, direct=FALSE, lvl=0) {
  # test w/ golden moles 9389
  downloader <- function(ids, ps) {
    ftch_args <- list(db = "nucleotide", rettype = 'gb', retmode = 'text',
                      id = ids)
    search_and_cache(func = rentrez::entrez_fetch, args = ftch_args,
                     fnm = 'fetch', ps = ps)
  }
  # get accessions
  sids <- sids_get(txid = txid, direct = direct, ps = ps)
  if (length(sids) < 1) {
    return(list())
  }
  if (length(sids) > ps[['mdlthrs']]) {
    info(lvl = lvl + 3, ps = ps, "More than [", ps[['mdlthrs']],
         ' sqs] available. Choosing at random.')
    sids <- sids[sample(seq_along(sids), size = ps[['mdlthrs']],
                        replace = FALSE)]
  }
  if (restez::restez_ready()) {
    info(lvl = lvl + 3, ps = ps, "Getting [", length(sids),
         " sqs] from restez database...")
    unversioned <- sub(pattern = '\\.[0-9]+', replacement = '', x = sids)
    raw_recs <- restez::gb_record_get(id = unversioned)
    sids <- sids[!unversioned %in% names(raw_recs)]
  } else {
    raw_recs <- NULL
  }
  if (length(sids) > 0) {
    info(lvl = lvl + 3, ps = ps, "Downloading [", length(sids), " sqs] ...")
    raw_recs <- c(batcher(sids, func = downloader, ps = ps, lvl = lvl + 4),
                  raw_recs)
  }
  # split up by feature, return SeqRecs
  seqrecs <- numeric()
  for (i in seq_along(raw_recs)) {
    seqrecs <- c(seqrec_convert(raw_recs = raw_recs[[i]], ps = ps), seqrecs)
  }
  seqrecs
}
